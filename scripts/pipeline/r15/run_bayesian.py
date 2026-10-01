"""Hierarchical marginal beta-binomial strategy profiles and pressure associations.

Distinct outcomes are fitted in a batch. No cross-outcome joint-event inference.
"""
from pathlib import Path
import argparse,json,time,sys,platform
import numpy as np
import pandas as pd
from scipy.special import expit
import pymc as pm
import arviz as az

HERE=Path(__file__).resolve().parent; IN=HERE/'inputs'
META=json.loads((IN/'metadata.json').read_text())
TRAITS=META['traits']; DENOM=np.array(META['denominators'],int)

def save(d,path):d.to_csv(path,sep='\t',index=False)
def summarize(v):
 return {'median':float(np.median(v)),'mean':float(np.mean(v)),
  'low95':float(np.quantile(v,.025)),'high95':float(np.quantile(v,.975)),
  'p_positive':float(np.mean(v>0)),**{f'p_gt_{k}pp':float(np.mean(v>k/100)) for k in [1,5,10]},
  **{f'p_lt_minus_{k}pp':float(np.mean(v< -k/100)) for k in [1,5,10]},
  **{f'p_abs_lt_{k}pp':float(np.mean(abs(v)<k/100)) for k in [1,5,10]}}

def prepare(args):
 if args.model in ['global','simulation']:
  file='global_shared_families.tsv' if args.cohort=='shared' else 'global_strict.tsv'
  d=pd.read_csv(IN/file,sep='\t');sources=sorted(d.dataset_id.unique());families=sorted(d.family_path.unique())
  X=d[['z_'+c for c in META['global_covariates']]].to_numpy(float)
  s=pd.Categorical(d.dataset_id,categories=sources).codes;g=pd.Categorical(d.family_path,categories=families).codes
  k=d[TRAITS].to_numpy(int);n=np.broadcast_to(DENOM,k.shape).copy();truth=None
  covs=META['global_covariates']
  if args.model=='simulation':
   rng=np.random.default_rng(90927);sources=['sim_source_'+str(i) for i in range(7)];families=['sim_family_'+str(i) for i in range(12)]
   # Every family occurs in every source with two independent observations.
   pairs=np.array([(i,j) for i in range(7) for j in range(12) for _ in range(2)])
   s=pairs[:,0];g=pairs[:,1];X=rng.normal(size=(len(pairs),3))
   alpha=np.linspace(-2,1,len(TRAITS));sf=rng.normal(0,.4,(7,len(TRAITS)));sf-=sf.mean(axis=0)
   ff=rng.normal(0,.4,(12,len(TRAITS)));beta=rng.normal(0,.15,(3,len(TRAITS)));phi=np.ones(len(TRAITS))*18
   mu=expit(alpha+sf[s]+ff[g]+X@beta);n=np.broadcast_to(DENOM,mu.shape).copy();k=rng.binomial(n,rng.beta(mu*phi,(1-mu)*phi))
   d=pd.DataFrame({'dataset_id':np.array(sources)[s],'family_path':np.array(families)[g]})
   truth={'alpha':alpha.tolist(),'source_effect':sf.tolist(),'beta':beta.tolist(),'phi':phi.tolist()}
  return dict(d=d,X=X,k=k,n=n,s=s,g=g,sources=sources,groups=families,covs=covs,truth=truth,kind='global')
 cohort='strict' if args.cohort=='strict' else 'primary'
 d=pd.read_csv(IN/f'qaidam_{cohort}.tsv',sep='\t');cfg=META['environment_scales'][cohort]
 covs=cfg['environment_columns']+cfg['control_columns'];X=d[['z_'+c for c in covs]].to_numpy(float)
 groups=sorted(d.site_id.unique());g=pd.Categorical(d.site_id,categories=groups).codes
 k=d[TRAITS].to_numpy(int);n=d.n_MAGs.to_numpy(int)[:,None]*DENOM[None,:]
 return dict(d=d,X=X,k=k,n=n,g=g,groups=groups,covs=covs,env_cfg=cfg,kind='environment',truth=None)

def build(data,scale,pressure_unpooled=False):
 coords={'trait':TRAITS,'covariate':data['covs'],'group':data['groups'],'obs':np.arange(len(data['k']))}
 if data['kind']=='global':coords['source']=data['sources']
 with pm.Model(coords=coords) as model:
  alpha=pm.Normal('alpha',0,1.5,dims='trait')
  group_sigma=pm.HalfNormal('group_sigma',1. if data['kind']=='global' else .5,dims='trait')
  z=pm.Normal('group_z',0,1,dims=('group','trait'))
  group_effect=pm.Deterministic('group_effect',z*group_sigma,dims=('group','trait'))
  phi=pm.LogNormal('phi',np.log(15),1,dims='trait')
  if data['kind']=='global':
   beta=pm.Normal('beta',0,.5,dims=('covariate','trait'))
   source_sigma=pm.HalfNormal('source_sigma',.5*scale,dims='trait')
   sz=pm.Normal('source_z',0,1,dims=('source','trait'))
   source_effect=pm.Deterministic('source_effect',(sz-sz.mean(axis=0))*source_sigma,dims=('source','trait'))
   eta=alpha+source_effect[data['s']]+group_effect[data['g']]+pm.math.dot(data['X'],beta)
   pm.Deterministic('reference_probability',pm.math.sigmoid(alpha+source_effect),dims=('source','trait'))
  else:
   # A shrinkage scale for each pressure/level, without forcing signs to agree.
   model.add_coord('pressure',data['covs'][:6]);model.add_coord('control',data['covs'][6:])
   pressure_sigma=(pm.Deterministic('pressure_sigma',pm.math.ones(6)*.5*scale,dims='pressure')
    if pressure_unpooled else pm.HalfNormal('pressure_sigma',.5*scale,dims='pressure'))
   bz=pm.Normal('pressure_z',0,1,dims=('pressure','trait'))
   bp=pm.Deterministic('pressure_beta',bz*pressure_sigma[:,None],dims=('pressure','trait'))
   bc=pm.Normal('control_beta',0,.5,dims=('control','trait'))
   beta=pm.Deterministic('beta',pm.math.concatenate([bp,bc],axis=0),dims=('covariate','trait'))
   eta=alpha+group_effect[data['g']]+pm.math.dot(data['X'],beta)
  prob=pm.math.sigmoid(eta)
  pm.BetaBinomial('observed_counts',n=data['n'],alpha=prob*phi,beta=(1-prob)*phi,
   observed=data['k'],dims=('obs','trait'))
 return model

def predictive_summary(reps,observed,n,groups,labels,path):
 # reps: posterior/prior draws x observations x traits
 rows=[];obs_prop=observed/n;rep_prop=reps/n[None,:,:]
 for group in ['ALL']+list(labels):
  mask=np.ones(len(observed),bool) if group=='ALL' else np.asarray(groups)==group
  for j,t in enumerate(TRAITS):
   metrics={'mean':(obs_prop[mask,j].mean(),rep_prop[:,mask,j].mean(axis=1)),
    'sd':(obs_prop[mask,j].std(),rep_prop[:,mask,j].std(axis=1)),
    'zero_fraction':((observed[mask,j]==0).mean(),(reps[:,mask,j]==0).mean(axis=1)),
    'full_fraction':((observed[mask,j]==n[mask,j]).mean(),(reps[:,mask,j]==n[mask,j]).mean(axis=1))}
   for metric,(value,v) in metrics.items():
    lo,med,hi=np.quantile(v,[.025,.5,.975]);rows.append(dict(group=group,trait=t,metric=metric,observed=value,
     predictive_low95=lo,predictive_median=med,predictive_high95=hi,observed_in_95=bool(lo<=value<=hi)))
 save(pd.DataFrame(rows),path)

def finish(idata,data,out,args):
 post=idata['posterior'].to_dataset()
 core=['alpha','phi','beta','group_sigma']+(['source_effect','source_sigma'] if data['kind']=='global' else ['pressure_beta']+([] if getattr(args,'pressure_unpooled',False) else ['pressure_sigma']))
 summary=az.summary(idata,var_names=core,kind='all',round_to="none",ci_prob=.95)
 save(summary.rename_axis('parameter').reset_index(),out/'parameter_summary.tsv')
 all_summary=az.summary(idata,var_names=core+['group_effect'],kind='diagnostics',round_to="none")
 save(all_summary.rename_axis('parameter').reset_index(),out/'all_diagnostics.tsv')
 stats=idata['sample_stats'].to_dataset();divergence_name=next((c for c in stats if 'diverg' in c),None)
 rhat=next(c for c in summary if 'r_hat' in c or c=='rhat')
 essb=next(c for c in summary if 'ess_bulk' in c);esst=next(c for c in summary if 'ess_tail' in c)
 diag={'rhat_max_core':float(summary[rhat].max()),'ess_bulk_min_core':float(summary[essb].min()),
  'ess_tail_min_core':float(summary[esst].min()),'divergences':int(stats[divergence_name].sum()) if divergence_name else None,
  'rhat_max_all':float(all_summary[rhat].max()),'ess_bulk_min_all':float(all_summary[essb].min()),
  'chains':int(post.sizes['chain']),'draws_per_chain':int(post.sizes['draw']),
  'stat_variables':list(stats.data_vars)}
 diag['core_gate_pass']=bool(diag['rhat_max_core']<1.01 and diag['ess_bulk_min_core']>400 and diag['ess_tail_min_core']>400 and diag['divergences']==0)
 (out/'diagnostics.json').write_text(json.dumps(diag,indent=2))
 np.savez_compressed(out/'posterior_arrays.npz',**{v:post[v].values for v in core+['group_effect']})
 flat={v:post[v].values.reshape((-1,)+post[v].values.shape[2:]) for v in core+['group_effect']}
 if data['kind']=='global':
  ref=expit(flat['alpha'][:,None,:]+flat['source_effect']);base=expit(flat['alpha']);rows=[];pairs=[];ranges=[]
  for i,source in enumerate(data['sources']):
   for j,trait in enumerate(TRAITS):
    r=summarize(ref[:,i,j]-base[:,j]);ps=summarize(ref[:,i,j])
    rows.append(dict(source=source,trait=trait,n_source=int((data['d'].dataset_id==source).sum()),**r,
     reference_probability_median=ps['median'],reference_probability_low95=ps['low95'],reference_probability_high95=ps['high95']))
  for i in range(len(data['sources'])):
   for j in range(i+1,len(data['sources'])):
    for k,trait in enumerate(TRAITS):pairs.append(dict(source_a=data['sources'][i],source_b=data['sources'][j],trait=trait,**summarize(ref[:,i,k]-ref[:,j,k])))
  for j,trait in enumerate(TRAITS):
   v=ref[:,:,j].max(axis=1)-ref[:,:,j].min(axis=1)
   ranges.append(dict(trait=trait,**summarize(v)))
  save(pd.DataFrame(rows),out/'source_profile_posteriors.tsv');save(pd.DataFrame(pairs),out/'source_pairwise_posteriors.tsv');save(pd.DataFrame(ranges),out/'source_range_posteriors.tsv')
  if data['truth']:
   rows=[]
   for i,s in enumerate(data['sources']):
    for j,t in enumerate(TRAITS):
     truth=data['truth']['source_effect'][i][j];v=flat['source_effect'][:,i,j];lo,med,hi=np.quantile(v,[.025,.5,.975])
     rows.append(dict(source=s,trait=t,truth=truth,median=med,low95=lo,high95=hi,covered=bool(lo<=truth<=hi)))
   save(pd.DataFrame(rows),out/'simulation_recovery.tsv')
 else:
  rows=[]
  for i,cov in enumerate(data['covs'][:6]):
   step=1/data['env_cfg']['sds'][cov]
   for j,trait in enumerate(TRAITS):
    b=flat['beta'][:,i,j];v=expit(flat['alpha'][:,j]+b*step/2)-expit(flat['alpha'][:,j]-b*step/2)
    rows.append(dict(exposure=cov,trait=trait,**summarize(v),logit_beta_median=float(np.median(b))))
  save(pd.DataFrame(rows),out/'environment_posteriors.tsv')
 # Exact numpy beta-binomial PPC from posterior parameters. No fitted-value plug-in.
 rng=np.random.default_rng(args.seed+71);ids=rng.choice(len(flat['alpha']),size=min(400,len(flat['alpha'])),replace=False);reps=[]
 for ix in ids:
  eta=flat['alpha'][ix]+flat['group_effect'][ix,data['g']]+data['X']@flat['beta'][ix]
  if data['kind']=='global':eta+=flat['source_effect'][ix,data['s']]
  p=expit(eta);phi=flat['phi'][ix];latent=rng.beta(p*phi,(1-p)*phi);reps.append(rng.binomial(data['n'],latent))
 reps=np.asarray(reps)
 groups=data['d'].dataset_id.to_numpy() if data['kind']=='global' else data['d'].site_id.to_numpy()
 labels=data['sources'] if data['kind']=='global' else data['groups']
 predictive_summary(reps,data['k'],data['n'],groups,labels,out/'posterior_predictive_checks.tsv')
 # Only a thin trace is needed for visual chain inspection; full posterior is saved.
 np.savez_compressed(out/'trace_core.npz',alpha=post.alpha.values,beta=post.beta.values,
  scale=post['source_sigma' if data['kind']=='global' else 'pressure_sigma'].values)
 print('DIAGNOSTICS',json.dumps(diag),flush=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--model',choices=['global','environment','simulation'],required=True)
 ap.add_argument('--cohort',default='primary');ap.add_argument('--scale',type=float,default=1);ap.add_argument('--tag',required=True)
 ap.add_argument('--draws',type=int,default=1000);ap.add_argument('--tune',type=int,default=1000)
 ap.add_argument('--seed',type=int,default=260927);ap.add_argument('--target-accept',type=float,default=.95)
 ap.add_argument('--pressure-unpooled',action='store_true')
 args=ap.parse_args();out=HERE/'results'/args.tag
 if out.exists() and (out/'posterior.nc').exists():raise SystemExit('Refusing to overwrite posterior; use a new tag')
 out.mkdir(parents=True,exist_ok=True);start=time.time();data=prepare(args)
 assert (data['k']>=0).all() and (data['k']<=data['n']).all()
 (out/'run_config.json').write_text(json.dumps(dict(vars(args),n_observations=len(data['k']),kind=data['kind'],traits=TRAITS,
  sources=data.get('sources',[]),groups=data['groups'],covariates=data['covs'],python=platform.python_version(),pymc=pm.__version__,arviz=az.__version__),indent=2))
 if data['truth']:(out/'simulation_truth.json').write_text(json.dumps(data['truth'],indent=2))
 save(data['d'],out/'analysis_rows.tsv')
 model=build(data,args.scale,args.pressure_unpooled)
 print('PRIOR',args.tag,len(data['k']),flush=True)
 with model:
  prior=pm.sample_prior_predictive(draws=200,var_names=['observed_counts'],random_seed=args.seed)
 values=prior['prior_predictive']['observed_counts'].values.reshape((-1,)+data['k'].shape)
 groups=data['d'].dataset_id.to_numpy() if data['kind']=='global' else data['d'].site_id.to_numpy()
 labels=data['sources'] if data['kind']=='global' else data['groups']
 predictive_summary(values,data['k'],data['n'],groups,labels,out/'prior_predictive_checks.tsv')
 print('SAMPLING',args.tag,flush=True)
 keep=['alpha','beta','phi','group_sigma','group_effect']+(['source_effect','source_sigma','reference_probability'] if data['kind']=='global' else ['pressure_beta','pressure_sigma','control_beta'])
 with model:
  idata=pm.sample(draws=args.draws,tune=args.tune,chains=4,cores=4,nuts_sampler='nutpie',target_accept=args.target_accept,
   random_seed=args.seed,progressbar=False,var_names=keep,idata_kwargs={'log_likelihood':False},compute_convergence_checks=True)
 print('SAVING',args.tag,flush=True)
 idata.to_netcdf(out/'posterior.nc',engine='h5netcdf')
 finish(idata,data,out,args)
 (out/'elapsed_seconds.txt').write_text(str(time.time()-start)+'\n')
 print('DONE',args.tag,'seconds',round(time.time()-start,1),flush=True)

if __name__=='__main__':main()
