"""Count-distribution revision justified by beta-binomial PPC; see amendment."""
from pathlib import Path
import argparse,json,time
import numpy as np
import pandas as pd
from scipy.special import expit,logit
import pymc as pm
import arviz as az
from run_bayesian import prepare,TRAITS,DENOM,summarize,predictive_summary,save

H=Path(__file__).resolve().parent
def simulate(data,seed=90928):
 rng=np.random.default_rng(seed);truth=data['truth'];a=np.array(truth['alpha'])
 sf=np.array(truth['source_effect']);beta=np.array(truth['beta'])
 ff=rng.normal(0,.4,(len(data['groups']),len(TRAITS)))
 eta=sf[data['s']]+ff[data['g']]+data['X']@beta
 truth['q']=[]
 for j,n in enumerate(DENOM):
  q=rng.dirichlet(np.ones(n+1)*2);c=logit(np.cumsum(q)[:-1])+a[j]
  q=np.diff(np.r_[0,expit(c),1]);truth['q'].append(q.tolist())
  cdf=expit(c[None,:]-eta[:,j,None]);u=rng.random(len(eta))[:,None]
  data['k'][:,j]=(u>cdf).sum(axis=1)

def build(data,scale,baseline_alpha=1.):
 coords={'trait':TRAITS,'source':data['sources'],'group':data['groups'],
  'covariate':data['covs'],'obs':np.arange(len(data['k']))}
 with pm.Model(coords=coords) as m:
  beta=pm.Normal('beta',0,.5,dims=('covariate','trait'))
  ss=pm.HalfNormal('source_sigma',.5*scale,dims='trait')
  sz=pm.Normal('source_z',0,1,dims=('source','trait'))
  sf=pm.Deterministic('source_effect',(sz-sz.mean(axis=0))*ss,dims=('source','trait'))
  gs=pm.HalfNormal('group_sigma',1,dims='trait')
  gz=pm.Normal('group_z',0,1,dims=('group','trait'))
  gf=pm.Deterministic('group_effect',gz*gs,dims=('group','trait'))
  eta=sf[data['s']]+gf[data['g']]+pm.math.dot(data['X'],beta)
  for j,n in enumerate(DENOM):
   m.add_coord(f'category_{j}',np.arange(n+1))
   q=pm.Dirichlet(f'q_{j}',a=np.ones(n+1)*baseline_alpha,dims=f'category_{j}')
   cumulative=pm.math.cumsum(q)[:-1]
   if baseline_alpha<1: cumulative=pm.math.clip(cumulative,1e-12,1-1e-12)
   c=pm.math.log(cumulative)-pm.math.log(1-cumulative)
   pm.OrderedLogistic(f'counts_{j}',eta=eta[:,j],cutpoints=c,
    observed=data['k'][:,j],compute_p=False,dims='obs')
 return m

def finish(dt,data,out,args):
 post=dt['posterior'].to_dataset();core=['beta','source_sigma','source_effect','group_sigma']+[f'q_{j}' for j in range(len(TRAITS))]
 summary=az.summary(dt,var_names=core,round_to='none',ci_prob=.95)
 save(summary.rename_axis('parameter').reset_index(),out/'parameter_summary.tsv')
 allsum=az.summary(dt,var_names=core+['group_effect'],kind='diagnostics',round_to='none')
 save(allsum.rename_axis('parameter').reset_index(),out/'all_diagnostics.tsv')
 diag=dict(rhat_max_core=float(summary.r_hat.max()),ess_bulk_min_core=float(summary.ess_bulk.min()),
  ess_tail_min_core=float(summary.ess_tail.min()),rhat_max_all=float(allsum.r_hat.max()),
  ess_bulk_min_all=float(allsum.ess_bulk.min()),divergences=int(dt['sample_stats']['diverging'].sum()),
  chains=int(post.sizes['chain']),draws_per_chain=int(post.sizes['draw']))
 diag['core_gate_pass']=bool(diag['rhat_max_core']<1.01 and diag['ess_bulk_min_core']>400 and diag['ess_tail_min_core']>400 and diag['divergences']==0)
 (out/'diagnostics.json').write_text(json.dumps(diag,indent=2))
 variables=core+['group_effect'];flat={v:post[v].values.reshape((-1,)+post[v].values.shape[2:]) for v in variables}
 np.savez_compressed(out/'posterior_arrays.npz',**{v:post[v].values for v in variables})
 ns=len(flat['beta']);ref=np.empty((ns,len(data['sources']),len(TRAITS)));base=np.empty((ns,len(TRAITS)));cuts=[]
 for j,n in enumerate(DENOM):
  c=logit(np.clip(np.cumsum(flat[f'q_{j}'],axis=1)[:,:-1],1e-12,1-1e-12));cuts.append(c)
  ref[:,:,j]=expit(flat['source_effect'][:,:,j,None]-c[:,None,:]).sum(axis=2)/n
  base[:,j]=expit(-c).sum(axis=1)/n
 rows=[];pairs=[];ranges=[]
 for i,s in enumerate(data['sources']):
  for j,t in enumerate(TRAITS):
   r=summarize(ref[:,i,j]-base[:,j]);q=summarize(ref[:,i,j]);rows.append(dict(source=s,trait=t,
    n_source=int((data['d'].dataset_id==s).sum()),**r,reference_probability_median=q['median'],
    reference_probability_low95=q['low95'],reference_probability_high95=q['high95']))
 for i,s in enumerate(data['sources']):
  for k in range(i+1,len(data['sources'])):
   for j,t in enumerate(TRAITS):pairs.append(dict(source_a=s,source_b=data['sources'][k],trait=t,**summarize(ref[:,i,j]-ref[:,k,j])))
 for j,t in enumerate(TRAITS):ranges.append(dict(trait=t,**summarize(ref[:,:,j].max(axis=1)-ref[:,:,j].min(axis=1))))
 save(pd.DataFrame(rows),out/'source_profile_posteriors.tsv');save(pd.DataFrame(pairs),out/'source_pairwise_posteriors.tsv');save(pd.DataFrame(ranges),out/'source_range_posteriors.tsv')
 if data['truth']:
  rows=[]
  for i,s in enumerate(data['sources']):
   for j,t in enumerate(TRAITS):
    v=flat['source_effect'][:,i,j];lo,med,hi=np.quantile(v,[.025,.5,.975]);true=data['truth']['source_effect'][i][j]
    rows.append(dict(source=s,trait=t,truth=true,median=med,low95=lo,high95=hi,covered=bool(lo<=true<=hi)))
  save(pd.DataFrame(rows),out/'simulation_recovery.tsv')
 rng=np.random.default_rng(args.seed+71);ids=rng.choice(ns,min(400,ns),replace=False);reps=[]
 for ix in ids:
  eta=flat['source_effect'][ix,data['s']]+flat['group_effect'][ix,data['g']]+data['X']@flat['beta'][ix]
  rep=np.empty_like(data['k'])
  for j,n in enumerate(DENOM):
   cdf=expit(cuts[j][ix,None,:]-eta[:,j,None]);rep[:,j]=(rng.random((len(rep),1))>cdf).sum(axis=1)
  reps.append(rep)
 reps=np.asarray(reps)
 predictive_summary(reps,data['k'],data['n'],data['d'].dataset_id.to_numpy(),data['sources'],out/'posterior_predictive_checks.tsv')
 catrows=[]
 for source in ['ALL']+data['sources']:
  mask=np.ones(len(data['k']),bool) if source=='ALL' else data['d'].dataset_id.to_numpy()==source
  for j,n in enumerate(DENOM):
   for value in range(n+1):
    observed=(data['k'][mask,j]==value).mean();v=(reps[:,mask,j]==value).mean(axis=1)
    lo,med,hi=np.quantile(v,[.025,.5,.975]);catrows.append(dict(group=source,trait=TRAITS[j],category=value,
     observed=observed,predictive_low95=lo,predictive_median=med,predictive_high95=hi,observed_in_95=bool(lo<=observed<=hi)))
 save(pd.DataFrame(catrows),out/'posterior_predictive_categories.tsv')
 np.savez_compressed(out/'trace_core.npz',beta=post.beta.values,scale=post.source_sigma.values)
 print('DIAGNOSTICS',json.dumps(diag),flush=True)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--model',default='global',choices=['global','simulation']);ap.add_argument('--cohort',default='primary')
 ap.add_argument('--scale',type=float,default=1);ap.add_argument('--tag',required=True);ap.add_argument('--seed',type=int,default=260927)
 ap.add_argument('--baseline-alpha',type=float,default=1.)
 ap.add_argument('--draws',type=int,default=1500);ap.add_argument('--tune',type=int,default=1500);ap.add_argument('--target-accept',type=float,default=.95)
 args=ap.parse_args();out=H/'results'/args.tag
 if (out/'posterior.nc').exists():raise SystemExit('Refusing posterior overwrite')
 out.mkdir(parents=True,exist_ok=True);start=time.time();data=prepare(args)
 if args.model=='simulation':simulate(data)
 (out/'run_config.json').write_text(json.dumps(dict(vars(args),likelihood='ordered_logistic',n_observations=len(data['k']),
  sources=data['sources'],groups=data['groups'],traits=TRAITS,covariates=data['covs']),indent=2))
 if data['truth']:(out/'simulation_truth.json').write_text(json.dumps(data['truth'],indent=2))
 save(data['d'],out/'analysis_rows.tsv');m=build(data,args.scale,args.baseline_alpha)
 print('PRIOR',args.tag,len(data['k']),flush=True)
 with m:prior=pm.sample_prior_predictive(draws=200,var_names=[f'counts_{j}' for j in range(len(TRAITS))],random_seed=args.seed)
 reps=np.stack([prior['prior_predictive'][f'counts_{j}'].values.reshape((-1,len(data['k']))) for j in range(len(TRAITS))],axis=-1)
 predictive_summary(reps,data['k'],data['n'],data['d'].dataset_id.to_numpy(),data['sources'],out/'prior_predictive_checks.tsv')
 keep=['beta','source_sigma','source_effect','group_sigma','group_effect']+[f'q_{j}' for j in range(len(TRAITS))]
 print('SAMPLING',args.tag,flush=True)
 with m:dt=pm.sample(draws=args.draws,tune=args.tune,chains=4,cores=4,nuts_sampler='nutpie',target_accept=args.target_accept,
  random_seed=args.seed,progressbar=False,var_names=keep,idata_kwargs={'log_likelihood':False})
 dt.to_netcdf(out/'posterior.nc',engine='h5netcdf');finish(dt,data,out,args)
 (out/'elapsed_seconds.txt').write_text(str(time.time()-start)+'\n');print('DONE',args.tag,round(time.time()-start,1),flush=True)

if __name__=='__main__':main()
