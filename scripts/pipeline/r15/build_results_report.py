"""Render the final interpretation from completed, named model outputs."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from summarize_results import GLOBAL,ENV,md_table,TRAITS
H=Path(__file__).resolve().parent;R=H/'results'
def read(tag,name):return pd.read_csv(R/tag/name,sep='\t')
def interval(r,sign=1):
 lo,hi=(r.low95,r.high95) if sign==1 else (-r.high95,-r.low95)
 return f'{sign*r["median"]*100:.2f} [{lo*100:.2f}, {hi*100:.2f}]'

def main():
 diag=pd.read_csv(R/'diagnostics_summary.tsv',sep='\t');pp=pd.read_csv(R/'predictive_check_summary.tsv',sep='\t')
 selected=diag[diag.model.isin(GLOBAL+ENV)]
 assert selected.core_gate_pass.all(), 'Unresolved sampling diagnostics; do not publish final results'
 pairs=[]
 for tag in GLOBAL:
  a=read(tag,'source_pairwise_posteriors.tsv')
  r=a[a.trait.eq('osmotic_desiccation_salt')&a.source_a.eq('qaidam_basin_mmag_1773')&a.source_b.eq('stordalen_mire_2019_hybrid_mags')].iloc[0]
  pairs.append({'模型':tag,'Stordalen−Qaidam差值pp [95%CrI]':interval(r,-1),'P差值>5pp':round(r.p_lt_minus_5pp,4),'P差值>10pp':round(r.p_lt_minus_10pp,4)})
 rows=[]
 for tag in ENV:
  a=read(tag,'environment_posteriors.tsv')
  rows.append({'模型':tag,'60项中位数范围pp':f'{a["median"].min()*100:.2f} 至 {a["median"].max()*100:.2f}',
   '最大方向后验概率':round(np.maximum(a.p_positive,1-a.p_positive).max(),3),
   '95%区间不含0的项数':int(((a.low95>0)|(a.high95<0)).sum()),'最小P绝对效应<5pp':round(a.p_abs_lt_5pp.min(),3)})
 # Actual fit checks, not declarations of adequacy based on convergence alone.
 dtable=selected[['model','rhat_max_core','ess_bulk_min_core','ess_tail_min_core','divergences']].copy().round(3)
 ptable=pp[pp.model.isin(GLOBAL+ENV+['global_main'])][['model','overall_outside95','overall_checks','all_groups_outside95','all_group_checks','category_all_outside95','category_all_checks','category_group_outside95','category_group_checks']].fillna('—')
 recovery=[]
 for tag in ['simulation','ordinal_simulation']:
  a=read(tag,'simulation_recovery.tsv');recovery.append(f'{tag}：{int(a.covered.sum())}/{len(a)}个真实来源效应被边际95%区间覆盖')
 prior=read('ordinal_baseline_diffuse_long','prior_predictive_checks.tsv');prior=prior[prior.group.eq('ALL')]
 baseline_prior_out=int((~prior.observed_in_95).sum())
 text=f'''# 贝叶斯环境策略分析与未知基因连接

2026-09-27；同一PaperSpine任务 `nee-mars-20260927`。本轮是新问题下的贝叶斯实算及有界文献对照，不重复旧KO修正、蛋白AI预测或DPFQ病例对照，不覆盖原稿。

## 当前可以使用的结果

**可以保留“不同来源存在谱系与质量条件下仍可见的功能配置差异”这一现象。其中渗透/干燥标记配置的来源差异最清楚；但还不能将它解释成盐度、干旱造成的适应机制，也没有证据把DPFQ008接成其必要基因。**

这里的“条件下”指所列协变量和科随机截距的模型条件，不是环境/研究混杂被完全消除。科随机项不是完整系统发育协方差；控制基因组大小也改变了估计对象，不能解释为环境的总效应。工作模型的局部预测失配仍需保留。

Qaidam的实测水分、EC和TOC在本轮多变量模型中的关联总体较弱；这与全域来源差异是两个不同问题，不能用前者否定后者，也不能拿后者声称这些实测压力已经得到机制验证。贝叶斯分析并未把旧p值换个名字，而是估计效应分布、检查收缩和模型适配。

## 1. 分析对象和实际模型

- 全域：879个严格90/5细菌ANI物种代表、7来源、272科路径。来源MAG数为Alaska216、Atacama卤石7、盐壳3、Australia150、Mauna Loa31、Qaidam339、Stordalen133。仅跨至少2来源的77科子集含560代表；小来源不能与大来源同等精度解读。
- 全部固定10终点：8个可变core65面板及12类氨基酸/8类辅因子的路线支持计数。光能面板在严格细菌队列全0，不能估计，不能称无光能功能。
- 本轮“功能轮廓”由上述10个面板的支持比例组成；它不是同样KO数量下的基因身份差异，也没有估计跨面板协同或功能权衡。不能把这里的计数差值与07/08旧配置距离互换。
- 全域采用有序logistic计数似然、来源部分收缩、科随机截距及完整度、污染、log2基因组大小控制。初始Beta-binomial无法充分描述部分欠离散计数，作为诊断保留。修订理由在 `ordinal_amendment_zh.md`，不是伪称事前注册。
- Qaidam主队列1484MAG/54样本/16地点；严格330MAG/39样本/15地点。按样本聚合支持计数，分母=MAG数×面板上界；拟合Beta-binomial工作似然、地点随机截距，压力按地点内/间分解，调整pH、深度、质量、大小与3个科组成PC。不是每个MAG一份独立环境重复，也不是丰度加权群落功能。
- 4链；全域主模型每链1500后验抽样，基线稀疏先验敏感性延长至3000；Qaidam每链2000。主/敏感性具体参数、随机种子、输入哈希、完整后验和软件锁定文件均保留。

## 2. 来源轮廓：敏感性分析中保留的差异

来源各自具有自己的效应，不要求七种环境同向。主有序模型中，渗透/干燥标记支持相对模型参考的中位差为：Alaska +4.61pp、Stordalen +7.82pp、Qaidam −9.82pp、Mauna Loa −9.31pp。其余来源及全部10终点、70轮廓和210两两对比完整保留，不以这四个数代替全结果。

全域缺少完整一致的采样单元层级，未另加入全球样本/地点随机项。879个基因组不是879份独立环境重复；后验区间不涵盖未建模的采样设计、空间结构或注释批次不确定性。这一部分应作为所收集基因组的条件来源轮廓，而非环境处理效应。

最直观的一个探索性对比是Stordalen与Qaidam的渗透标记支持。表内不是同一物种在两环境中的实验响应；共有科分析是改变分析人群的敏感性，也不是独立复现。

{md_table(pd.DataFrame(pairs))}

全部概率为有限MCMC抽样估计，显示1.0000不表示数学上确定。5/10pp是报告用的效应幅度，并非实验验证的生物阈值。来源最大−最小的区间天然倾向正值，不能据其不含0宣布“异质性显著”。

较干环境的所选渗透标记较少，不等于其耐旱能力较差：面板仅覆盖已定义的基因，谱系、替代实现、研究设计和检测仍可能共同贡献。当前也不能把这解释成某个未知基因弥补了缺口。此处最有价值的是待解释的配置差异，不是预先把环境排成极端程度等级。

## 3. Qaidam：明确区分压力关联与来源差异

效应是参考协变量及地点效应为0时，压力倍增所对应的预期支持比例变化；within和between分别报告。每个模型60项完整保留。

{md_table(pd.DataFrame(rows))}

主模型跨10终点共享每种压力的收缩尺度。先验尺度半倍/双倍检查与固定独立Normal(0,0.5)系数的“不共享收缩”检查必须一起阅读；即便P绝对效应<5pp较高，也只是所用工作模型、参考条件和面板下的条件概率，不能改写成生态等效或压力没有作用。

不共享收缩时出现较明显但不确定的方向线索，例如地点间TOC倍增对应AA支持+1.59pp（95%CrI −0.46至3.71），地点内水分倍增对应冷/蛋白质量面板+1.85pp（−0.67至4.35）。它们是全60项中的说明性例子，均保留负向或零附近的可能，不单独升级成新发现。主层级模型更靠近零，说明“效应很小”的程度依赖收缩假设。

不能用“60项95%区间都跨0”再做一次非黑即白的否决。例如不共享收缩时，地点间EC倍增与痕量气体面板支持下降的后验概率为97.0%，中位差−1.04pp；但主层级模型对应为71.6%、−0.069pp。这是有方向线索但对收缩设定敏感的结果，适合保留为不确定性明确的候选关联，尚不是稳定主结论。

压力支持域有限：主队列每种压力只有10/16地点有地点内变化，严格队列8/15；地点内TOC跨度达到一倍增的地点主队列5个、严格仅2个。水分对应5/5个、EC8/7个。倍增是标准化的模型比较，部分地点没有覆盖完整倍增，不能作逐地点普适预测。设计矩阵含截距15列，主/严格均满秩，条件数约6.07/4.88；满秩不等于无混杂或生物上充分识别。逐地点范围及全设计列相关矩阵见support_audit。

## 4. 采样、分布检查与仍有的限制

{md_table(dtable)}

以上为未舍入数值判定后的展示；门槛Rhat<1.01、bulk/tail ESS>400且无发散。原短跑environment_main的tail ESS约394未通过，已用长跑替代；稀疏基线先验初次1500抽样中4个极小类别概率Rhat/ESS不足，保留ordinal_baseline_diffuse并延长至3000，不隐去失败记录。极小/稀少功能和参数的后验仍可能依赖先验。

**先验预测没有被当作“通过”。** 主有序模型总体40项中25项观测在95%先验预测区间之外，5个终点的均值低于先验范围；Qaidam主模型为14/40项，3个低支持终点的均值低于先验范围。Dirichlet(1)类别先验诱导的支持均值偏向0.5，Normal截距先验也低估部分面板的稀少程度。数据更新后的PPC改善不抹去这种张力。根据独立复核，另将基线Dirichlet每类别浓度降至0.25并实际重拟合（上述ordinal_baseline_diffuse_long），其总体先验PPC越界{baseline_prior_out}/40；来源对比的变化与其他敏感性并列报告。实际检查仍未覆盖科效应、注释检测模型等全部先验选择，不能全面声称“先验稳健”。完整先验预测表保留。

下表列“观测值落在95%后验预测区间外的诊断条目数”；不是多重检验p值或自动通过/失败比率。总体包括均值、SD、零比例、满比例共40项；来源/地点分层包含总体。类别检查适用于全域有序计数模型。

{md_table(ptable)}

主有序模型总体40摘要及93类别检查均覆盖观测；但来源层面仍有局部失配，尤其Stordalen渗透计数的SD、痕量气体零比例，以及Mauna Loa辅因子SD。完整类别分布也有偏离，见 `posterior_predictive_categories.tsv`，不能称模型完全拟合或验证了比例优势假设。Qaidam主休眠SD、严格队列若干SD偏离，结果仍是工作模型的条件估计。PPC是拟合诊断，不是外部泛化检验。

代数检查：7000组随机类别概率的参考均值公式与直接类别加权均值最大差4.44e−16，切点均有序。参数恢复：{'；'.join(recovery)}。每模型仅一个模拟数据集，不能称完成全面校准；模拟覆盖不等于真实生态结论正确。各终点均为边际模型，没有残差跨终点关联，不能计算“所有功能共同响应”的联合后验，或声称完全消除了筛选带来的多重性。

## 5. 文献与AI新基因如何进入主线

详见 [文献对照](literature_gene_bridge_zh.md)。最贴近的路径来自FESNov（未知家族质量、保守性、邻域/结构与验证）、FUGAsseM（多证据集成与留出评价）以及南极供能策略研究（生态策略与候选酶功能测量）。NEE生态位宽度研究则示范了如何用可量化生态概念连接基因组策略。以上是方法借鉴，不是我们研究的外部复现。

建议保留的研究假说是：**环境相关功能配置可能由不同谱系的不同遗传实现支撑，其中未注释家族可能是条件相关的辅助成分或替代实现。** 此假说尚未验证，不应改写为本轮中心结论。

DPFQ008当前只可作为候选例子：既有88代表中qcr三位置齐全50个，其中25个未检出候选；严格子集13个中5个未检出。六个原固定面板全部在 [覆盖审计](candidate_coverage_zh.md) 中报告。没有“标记阳性必有DPFQ008”的观测对应；但未检出不等于真实缺失，不能据此排除条件特异作用。qcr邻域参与过选基因，不能再独立验证其功能，且没有证据把它连接到本轮最清楚的渗透面板差异。

可用的证据链应是：独立定义环境与功能配置 → 在可比谱系中检验未知家族关联 → 用留出地点/家族或外部数据检验可迁移性 → 由AI序列/结构提出具体功能假说 → 根据条件相关功能或必要性问题选择测量。必要性需与条件明确的基因扰动、适合度损失和互补恢复对应；不必强求所有生态相关基因都是必需基因。

纯计算仍可继续检验生态关联和预测可迁移性，不需要先承诺湿实验才能推进。只有将结论升级为具体生化功能或条件必需性时，才需要与该命题相称的功能证据。论文主问题可以是环境相关遗传实现，不能因为尚未证明必需性就把所有未知基因研究否定。

## 6. 对NEE主线的实际判断

本轮不是“没有可用结果”：来源相关配置差异可作为现象层结果，完整输入和模型可复用；文献支持用生态策略组织未知基因研究。另一方面，仅凭来源差异、AI结构及既有qcr邻域，目前仍不足以宣称新的生态适应机制或NEE中心创新已经成立。

下一段应围绕已经定义的生态问题做未知家族的独立筛选/功能解释可行性审计，先识别可比较谱系和检测灵敏度；不继续更换统计方法去争取一个更好看的概率。也不把DPFQ008强接到渗透面板上。火星/月球的联系限于环境约束下的生命检测靶标或候选生物过程；ISRU生产能力还需底物、供水/供氧、通量和产率证据。

## 文件与复核

- `complete_posterior_tables_zh.md`：全70来源轮廓及主/严格/不共享收缩的全60压力效应。
- `results/all_source_posteriors.tsv`、`all_environment_posteriors.tsv`：完整精度及敏感性结果；各模型目录保留全部两两比较、后验和PPC。
- `figures/source_profiles.*`、`qaidam_pressure_associations.*`、`ordinal_count_checks.*`与链轨迹：PNG/PDF/SVG分析图。
- 独立方法/模型修订及已完成结果的有界复核见本目录review文件；其时间范围须按原文区分，早期复核不自动覆盖后来生成的全域后验。本轮不是全稿审稿或投稿就绪认证。
'''
 (H/'bayesian_results_zh.md').write_text(text)
 print('Wrote report',len(text),'characters')

if __name__=='__main__':main()
