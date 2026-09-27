## Abstract


## Introduction
#### Background
Spatial transcriptomics has transformed tissue biology by simultaneously capturing genome-wide gene expression and its native spatial coordinates, thereby preserving the positional context that bulk and dissociated single-cell assays discard. (REF) Because biological processes are intrinsically spatially organized, this technology has enabled discoveries that are inaccessible to non-spatial assays, including spatially resolved differential expression, the delineation of tissue niches and domains, and the mapping of ligand–receptor-mediated cell–cell communication within intact tissue architecture. (REF) As both resolution and throughput continue to increase, spatial transcriptomics is becoming a routine assay in biomedical research.

Yet biological interpretation rarely rests on individual genes. Phenotypic states and disease processes are executed by coordinated gene programs—pathways and regulatory modules whose constituent genes act in concert rather than in isolation. This is particularly relevant for spatial data, because pathway-level readouts are more robust to the technical noise, dropouts, and cell-type admixture that characterize spot-level measurements. (REF) Moreover, within a single pathway, individual genes may exert opposite directional effects on the biological process of interest; aggregating gene-level evidence into a signed pathway score is often the only way to obtain a coherent picture of the underlying activity. 
#### Problem
Existing pathway approaches have an inference limitation. 现有 pathway analysis 往往ignore spatial dependence 没有充分处理 spatial autocorrelation 以及由此带来的：dependence in pathway-level evidence
#### Method
We developed SpaPath：a linear mixed model framework to quantify pathway-level expression differences across tissue regions. 每个 pathway 在所有 regions 中联合拟合，不再把均值差异和方差差异作为两个主要生物学输出。它们保留为模型参数和诊断量，并共同决定最终的 probability score。
#### Benchmark
Simulation: 并与GSDensity/STPathway/naïve analysis比较。
#### Application
在 spatial brain transcriptomics 数据上识别：cell-type- and region-specific pathway activity
#### Conclusion


## Results 
### Method overview
SpaPath uses a Gaussian linear mixed model (LMM) for normalized expression data and a negative binomial generalized linear mixed model (GLMM) for raw UMI counts. 两者都包含随机效应。Gaussian 版本用正态误差描述观测波动, SCT/normalized read counts 已经校正 library size，因此默认 additive offset 为 0；不再次加入 log library size；NB 版本通过负二项分布描述计数波动，通常不再另加一个 Gaussian 观测误差项 $\epsilon_{ig}$.

Gaussian / Normalized reads：
设 $i$ 表示 spot, $g$ 表示基因, $r(i)$ 表示 spot 所属的区域。对于一条 pathway，代码拟合：

$$
Y_{ig}=X_i^T\beta+\mu_{r_i}+b_{g,r_i}+w_g(s_i)+\epsilon_{ig},
\quad
\epsilon_{ig}\sim N(0,\sigma_\epsilon^2), 
\quad
b_{gr}\mid\tau_r^2\sim N(0,\tau_r^2).
$$

- $X_i$: spot-level covariates. (batch effect)
- $\beta$: 跨 genes 共享的 fixed effects
- $\mu_r$：区域 $r$ 的通路平均水平
- $b_{gr}$：基因 $g$ 在区域 $r$ 中，相对于通路平均水平的偏差
- $w_g(s_i)$: 个基因在组织位置 $s_i$ 附近的连续空间变化
- $\epsilon_{ig}$: 其余未解释的观测波动

NB / Raw reads：
每个基因有自己的偏离，但这些偏离共享一个分布，其离散程度由 $\tau_r^2$ 控制。
Small $\tau_r^2$：gene expression比较接近通路共同水平。Large $\tau_r^2$：更分散。
数据较少、噪声较大的基因，其 $b_{gr}$ 会更多地被shrinkage向 0。

$$
Y_{ig}\sim NB(m_{ig},n_g),
$$

$$
\log m_{ig}=\log L_i+X_i^T\beta+\mu_{r_i}+b_{g,r_i}+w_g(s_i)
$$

- $L_i$: 第 $i$ 个cell(cell-level-resolution)或spot(sopt-level-resolution)的 library size，也就是这个样本总共测到了多少 reads/UMIs. 
- $\log L_i$: offset. 因为测得越深，一个基因被观察到的 count 自然越大，account for differences in library size across cells.

$$
E(Y_{ig})=m_{ig},\mathrm{Var}(Y_{ig})=m_{ig}+\alpha_g m_{ig}^{2}
$$

在相同的平均表达下, $\alpha_g$ 越大，计数越分散。这与 INLA 使用的 size 参数 $n_g$ 对应，关系为 $\alpha_g=1/n_g$.

$\alpha_g^{rel}$ 由全 transcriptome 的 DESeq2 dispersion trend/ maximum a posteriori (MAP) shrinkage estimation 预估并进行 geometric-mean normalization。假设全转录组中用于归一化的有效基因，其 dispersion 几何均值为：

$$
G=\exp\left(\frac{1}{J}\sum_{g=1}^{J}\log\widehat\alpha_g\right), \quad\alpha_g^{rel}=\frac{\widehat\alpha_g}{G}
$$

这里的 dispersion trend 可以理解为：
1.先根据每个基因的数据，初步估计它自己的 dispersion。
2.利用大量基因，拟合平均表达水平与 dispersion 之间的总体趋势。
3.结合这个趋势，稳定各基因的估计。信息不足的基因通常更受趋势约束；信息充分的基因保留更多自己的数据特征。
所以，它利用其他基因的信息，减少单个基因估计过于不稳定的问题。收缩的参照与平均表达水平有关，并不是把所有基因都拉向同一个 dispersion。
整个数据集只做一次 DESeq2 预处理，保存每个基因的 $\alpha_g^{rel}$。分析不同通路时，取出对应基因的相对值重复使用。

Gene-specific overdispersion was modeled using a reduced-parameter negative-binomial formulation,

$$
\alpha_g=\alpha_0\alpha_g^{\mathrm{rel}},
$$

where $\alpha_g^{\mathrm{rel}}$ denotes the pre-estimated relative dispersion of gene $g$, and $\alpha_0$ is a pathway-level overall dispersion scale shared across genes. Under the parameterization $\alpha_g=1/n_g$, we define

$$
S_g=\frac{1}{\alpha_g^{\mathrm{rel}}}
$$

and write

$$
\alpha_0=\exp(-\theta_{\mathrm{NB}}),
$$

which gives

$$
n_g=S_g\exp(\theta_{\mathrm{NB}}).
$$

Thus, gene-to-gene differences in dispersion are fixed by the pre-estimated relative factors $S_g$, while the pathway model estimates only the single global calibration parameter $\theta_{\mathrm{NB}}$.

Spatial hierarchy:

$$
w_g(\cdot)\overset{iid}{\sim}GP(0,\sigma_s^2C_{\kappa_0}),
\qquad g=1,\ldots,K.
$$

- 每个 gene 有独立 spatial realization。
- 所有 genes 共用 mesh、固定的 $\kappa_0$ 和一个 spatial variance $\sigma_s^2$。
- spatial field 跨越完整组织，不在 region boundary 重新开始。
- 使用 INLA replicated SPDE；不使用 shared $w(s)$，不估计 gene-specific spatial variance。

因此自由 hyperparameter 数量为 $R+2$: $R$ regional variances、a spatial variance，and a Gaussian residual variance or a NB calibration parameter.

For region $r$ and $s$, 定义随机 pathway states：

$$
Z_r\sim N(\mu_r,\tau_r^2),
\qquad
Z_s\sim N(\mu_s,\tau_s^2),
$$

以及：

$$
\pi_{rs}=P( Z_r > Z_s\mid\mu_r,\tau_r^2,\mu_s,\tau_s^2).
$$



### Benchmarking on simulated spatial data

- GSVA etc single cell pathway method in ST data
- STpathway [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503753/)
  > employing the Brownian distance correlation (dCor) test in a novel set-up to evaluate the dependency of two random vectors with different dimensionality.
- GSDensity [paper](https://www.nature.com/articles/s41467-023-44206-x)
  > Here, we present GSDensity, a graph-modeling approach that allows users to obtain pathway-centric interpretation and dissection of single-cell and spatial transcriptomics (ST) data without performing clustering. Using pathway gene sets, we show that GSDensity can accurately detect biologically distinct cells and reveal novel cell-pathway associations ignored by existing methods. 
  
We generated two simulated datasets with different random seeds to benchmark SpaPath against three representative pathway-analysis workflows: GSVA, the gene-set scoring approach popularized in single-cell RNA sequencing, tested here with a spatial block bootstrap region test; STpathway, a distance-correlation (dCor) test for region-associated pathways; and GSDensity, evaluated by its published region enrichment test on dichotomized spot activity. Each dataset was built on the real spot coordinates of DLPFC section 151673 and comprised 50 non-overlapping 100-gene sets spanning five spatial activity patterns, with the number K of genes responding to the activity map prespecified for each gene set (Methods).

Across the ten groups defined by K and effect distribution, the existing methods did not provide discriminative region-level evidence. 
GSDensity's region test reflected strong simulated pathway activity: as activity increased, the difference in high-activity spot proportions between target and rest increased whereas the region p-value decreased. The pass rate of GSDensity's coherence pre-filter increased with the prespecified pathway activity; among the 100 prespecified gene sets, 35 passed this filter, which checks whether a set's genes cluster in the embedding space and whether the activity scores are reliable. Following the original publication, we report its region test only for those passing sets. The pass rate of GSDensity's coherence pre-filter increased with the prespecified pathway activity; however, of the 100 prespecified gene sets, only 35 passed this filter, which checks whether a set's genes cluster in the embedding space and whether its PAL activity scores are reliable. Following the original publication, we report its region test only for those passing sets. Among these sets, the median region p in the broad, weak-effect setting was 3.6 × 10⁻⁵ at K = 31 to 60 and 0.011 at K = 61 to 100, but not in the narrow, strong-effect setting, where the median p was 0.247 at K = 31 to 60 and 0.082 at K = 61 to 100. Additionally, in the high-K groups of the narrow, strong-effect setting, fewer than 1% of spots were labeled highly active in either target or rest, and the region p-value was therefore not meaningful.
STpathway's region dCor decreased with increasing K, indicating that it can reflect strong regional differences. However, its group medians spanned only 0.94 to 0.98, with an overall range of 0.87 to 0.99, offering limited discrimination of pathway activity. GSVA's region score difference increased with K, and its spatial block bootstrap test gave BH-adjusted median q-values below 0.05 in all groups with K of 11 or more, but it failed to discriminate activity strength among weakly responsive pathways. Furthermore, under the domain and interface patterns, the two pattern types with the largest target regions, the test returned extremely small p-values even when the prespecified pathway activity was low or absent: the K = 0 domain pathways attained p = 0.0005 and remained significant after BH correction in the datasets from both random seeds. This indicates that under these spatial patterns the test yields strong significance regardless of the prespecified activity level and cannot separate pattern-driven structure from true pathway responses.

SpaPath's posterior evidence instead increased monotonically as the prespecified pathway activity increased. For the Gaussian model, the median pai rose from 0.37 in the K = 0 to 4 group to 0.81 in the K = 61 to 100 group in the broad, weak-effect setting, and from 0.40 to 0.75 in the narrow, strong-effect setting. Posterior evidence reached about 0.98 already in the K = 11 to 30 group and was close to 1 in the highest-K groups. Importantly, pai and evidence also recovered the relative activity strength among weak pathways, where the existing methods could not: between the K = 0 to 4 and K = 5 to 10 groups, the median pai of the Gaussian model rose from 0.37 to 0.48 in the broad, weak-effect setting and from 0.40 to 0.51 in the narrow, strong-effect setting, and the median evidence rose from 0.15 to 0.39 and from 0.24 to 0.60, respectively. Thus, pathways with little or no prespecified activity received weak posterior support, yet the ordering of their activity levels was preserved. \[NB 的结果虽然也递增，但是在设定的弱活性通路下 evidence 约 0.6，原因待查\] This pattern was consistent across all five spatial activity patterns. 

[Table 1. Median region-level statistics of the four methods across the 100 simulated gene sets.](benchmark/analysis/simulation_data/result/simulation_summary.xlsx)

[Supplementary Table 1. Per-pathway simulation results for the simulation dataset.](benchmark/analysis/simulation_data/result/simulation_per_pathway_result.xlsx)
### Layer-resolved pathway activity in the human DLPFC

To evaluate the performance of these methods on real data, we analyzed a human dorsolateral prefrontal cortex (DLPFC) Visium dataset because its six cortical layers and white matter were manually annotated and the molecular functions of each layer are well documented. [Maynard et al., 2021] We applied SpaPath, GSVA, STpathway, and GSDensity to section 151673, which contains 3,639 spots. We tested a prespecified panel of 30 MSigDB gene sets: 9 for synaptic transmission, including glutamatergic, GABAergic, and NMDA/AMPA receptor programs; 3 for ion-channel activity; 1 for oxidative phosphorylation; 2 for myelination; 2 for lipid metabolism; and 3 for blood-brain barrier and endothelial programs, together with 10 reference sets for basal RNA and protein maintenance and DNA repair. 
[Supplementary Table 2. The 30 prespecified MSigDB gene sets for the DLPFC analysis.](benchmark/analysis/DLPFC/data/geneset_30.xlsx)

We compared the layer-versus-rest results of the four methods. With the spatial block bootstrap, the GSVA region test detected significant elevation of the synaptic transmission gene sets in layers 2 to 5 but not in layer 6 or white matter. The sets elevated in white matter were instead axon ensheathment, myelin maintenance, sphingolipid metabolism, and endothelial programs, with the largest region score difference, Δ = 0.40, observed for axon ensheathment. Almost no gene set was elevated in layer 1. When the spatial block bootstrap was replaced by a permutation bootstrap, the number of significant contrasts rose from 78 to 94 of the 210 layer-by-set tests, and the extra significant results included housekeeping gene sets in layers with no significant set under the spatial test. In layer 4 alone, the number of significant gene sets rose from 12 to 20, including the basal sets for translation, rRNA processing, mRNA splicing, and nucleotide excision repair. 
STpathway's dCor ranking identified regulation of synaptic plasticity as the most distinctive gene set in white matter, with dCor = 0.828, consistent with activity-dependent myelination, but its top-ranked sets in layers 1 and 6 were dominated by translation, rRNA processing, and mRNA splicing, showing limited layer specificity. 
Of the 30 prespecified gene sets, 19 passed GSDensity's coherence pre-filter. Among them, axon ensheathment and myelin maintenance were enriched in white matter, with p = 1.8 × 10⁻¹⁸⁶ and 3.0 × 10⁻²²², respectively. However, even among these passing sets, many contrasts had essentially no highly active spots in either region, and the region p-value was therefore not meaningful; regulation of synaptic plasticity, translation, and rRNA processing each showed this in six of the seven regions. In addition, 109 of the 210 contingency-table cells had expected counts below 5, violating the condition for its chi-square approximation.

### Whole Mouse brain
Data: [Single-cell spatial transcriptomic atlas of the whole mouse brain. (Han et al.)](https://www.cell.com/neuron/abstract/S0896-6273(25)00133-3)

Using single-nucleus RNA sequencing (snRNA-seq) and Stereo-seq, genome-wide coverage. 
> Coronal sections (10-μm thick) of the left hemisphere were prepared at 100-μm intervals, producing 123 sections (after quality control) with expression profiles for 29,655 genes, covering 95.5% of annotated protein-coding and non-coding genes. We defined the boundary of each cell based on the image of single-stranded DNA (ssDNA). A total of 4,229,623 cells were characterized in one mouse brain, with an average of 1,267 molecular identifiers (MIDs) and 668 detected genes per cell. Moreover, we manually segmented the brain areas for each section based on cytoarchitectural pattern and brain region delineation in ABA common coordinate framework (CCFv3)27.

The processed data can be downloaded in [Whole Mouse brain data](https://www.braindatacenter.cn/datacenter/web/#/dataSet/details?id=1888062851995131906). Raw data available upon request. 

The author has performed quality control steps, including:

- 去除 7 张捕获异常切片、5 张低质量切片，最终保留 195 张
- 用 SAW 处理原始序列，过滤低质量 UMI、进行基因组比对和 UMI 去重
- 使用检测基因数 100、线粒体比例 15% 两个阈值进行筛选，去除 20,019 个细胞，保留 4,209,604 个高质量细胞
- Bin100: 过滤低基因数、高线粒体比例的 bin，然后进行 SCTransform

For the normalized data using SCTransform, we model the transformed expression values using a Gaussian distribution. For the rae data, use a NB distribution.


### Tumor microenvironment
Data: 
1. [A single-cell and spatially resolved atlas of human breast cancers. (Wu et al.)](https://www.nature.com/articles/s41588-021-00911-1)
ST 部分包含 6 个原发乳腺癌样本：4 个 TNBC、2 个 ER+，使用 Visium，spot 直径约 55 μm。有独立的病理区域标注, 例如，可以把标注映射成肿瘤区、基质区、淋巴细胞富集区，再研究这些区域的通路差异，但是sample少。[data](https://zenodo.org/records/4739739)
3. 另有 sample 多的：[Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer with potential clinical implications. (Wang et al.)](https://www.nature.com/articles/s41467-024-54145-w)

   > In this retrospective study, 94 early-stage TNBC patients treated at the Institut Jules Bordet (Brussels, Belgium) with standard-of-care therapies between 2000 and 2016 were included. TNBC patients were selected based on their negative status for estrogen (ER) and progesterone (PR) receptors and the absence of HER2 amplification.
   > Information regarding germline mutations in BRCA1 and BRCA2 genes was extracted from the patients’ medical records. All patients underwent initial surgery, followed by adjuvant chemotherapy and/or radiotherapy.
   > Each TNBC sample was profiled for ST as triplicates consisting of three consecutive sections (16 μm-thick). As such, a 3D-stack of spatial gene expression could be obtained for each sample after image alignment.
   
   可能有点复杂....






## Discussion
#### Summary
SpaPath addresses a fundamental dependence problem in spatial pathway analysis.
#### 和现有方法区别
#### Biology insight 
#### Limitations
- computational cost (INLA对hyperparameter sensitive)
- sensitivity to spatial covariance assumptions
- pathway gene dependence (假设i.i.d 但实际是稀疏正负相关)
- pathway overlap
- uncertainty in cell type labels
- future work

## Methods
### Simulation framework

We simulate raw spatial transcriptomic counts at the 3,639 tissue spots of the spatialLIBD DLPFC 151673 section, using its observed spot coordinates. The primary dataset contains 10,000 genes: 50 non-overlapping test pathways of 100 genes each, 1,000 background spatially variable genes (SVGs), and 4,000 background genes without a spatial effect. Background SVGs receive independently generated spatial activity maps unrelated to the test-pathway maps.

We designed 50 pathways spanning five response-size groups, five spatial pattern classes, and two effect-size distributions. For each responsive gene, $E_g>0$ denotes the absolute target-versus-rest expected log2 fold change, adjusted for spot depth:

$$
E_g \sim
\begin{cases}
\mathrm{Lognormal}(\text{s}=0.3,\text{ scale}=0.4), & \text{broad and weak},\\
\mathrm{Lognormal}(\text{s}=0.4,\text{ scale}=0.9), & \text{narrow and strong}.
\end{cases}
$$

For each pathway $p$, $K_p$ is the number of the responsive members among the 100 genes, sampled from five ranges: 0–4, 5–10, 11–30, 31–60, or 61–100; the remaining $100-K_p$ members are not responsive. Response direction is drawn before count generation, with $d_g=+1$ for target-upregulated genes and $d_g=-1$ for target-downregulated genes. Genes in weaker response group are more likely to have their direction reversed, with reversal probabilities of 30–50%, 15–30%, and 10–15% in the weak, medium, anad strong groups, respectively. $E_g$ denotes the effect size of an individual responsive gene.

Each pathway is assigned an activity map $A_i^{(p)}$, normalized to $[0,1]$, from one of five classes: discrete domains, smooth gradients, radial or niche-centered patterns, localized hotspots, or interface-specific patterns. Domain activity is piecewise constant across $D_p$ spatial regions; gradients vary smoothly across the tissue; radial activity decays with distance from a niche center; hotspot activity is localized; and interface activity follows
$A_i^{(p)} = \mathbb{1}(h_i \leq r)\exp(-h_i/\lambda)$, where $h_i$ is the distance to the relevant boundary. Responsive genes receive their pathway's activity map through $H_{ig}=A_i^{(p)}$.

Raw read counts $Y_{ig}$ are simulated from a negative bionomial distribution:

$$
Y_{ig} \sim \mathrm{NB}(\mu_{ig}, \theta_g), \qquad
\mathrm{Var}(Y_{ig}) = \mu_{ig} + \frac{\mu_{ig}^2}{\theta_g},
$$

$$
\log \mu_{ig} = \log c_i + \log m_g + d_g b_g H_{ig}.
$$

Here $c_i$ is the spot-specific depth factor and $(m_g, \theta_g)$ is a jointly sampled empirical baseline-expression and dispersion pair. 

SpaPath NB receives raw counts and depth factors. GSDensity receives a Seurat object built from raw counts and uses the author function's internal NormalizeData() step. SpaPath Gaussian and STpathway receive the SCTransform Pearson residuals.
### Model detail










