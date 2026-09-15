## Abstract


## Introduction
#### Background
spatial transcriptomics 的意义: spatial transcriptomics preserves molecular + spatial information, biological processes are spatially organized. 空间转录组已经可以研究空间 DEG, cell-cell communication, etc... However, biological interpretation often relies on coordinated gene programs and pathways rather than individual genes. 很多pathway中不同gene可能也有正向，负向，只考虑gene-level的不能提供全面视角。
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
PathwayLGM uses a Gaussian linear mixed model (LMM) for normalized expression data and a negative binomial generalized linear mixed model (GLMM) for raw UMI counts. 两者都包含随机效应。Gaussian 版本用正态误差描述观测波动, SCT/normalized read counts 已经校正 library size，因此默认 additive offset 为 0；不再次加入 log library size；NB 版本通过负二项分布描述计数波动，通常不再另加一个 Gaussian 观测误差项 $\epsilon_{ig}$.

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

### Benchmark
simulation result

- GSVA etc single cell pathway method in ST data
- STpathway [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503753/)
  > employing the Brownian distance correlation (dCor) test in a novel set-up to evaluate the dependency of two random vectors with different dimensionality.

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
#### summary
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


### Simulation design

We simulated expression data for 10,000 genes, organized into 20 non-overlapping test pathways of 100 genes each; the remaining 8,000 genes formed a background pool containing both spatially variable genes (SVGs) and non-SVGs, with the spatial patterns of background SVGs generated independently of the test pathways. To assess the sensitivity of each method to the size of the gene universe, we additionally simulated datasets comprising 5,000 and 20,000 genes.

The 20 test pathways fell into three categories: (i) ten spatially coherent pathways, in which a specified proportion of responsive genes shared a common pathway-level spatial activity map and a common response direction; (ii) six directionally heterogeneous pathways, in which responsive genes shared the same spatial activity map but were assigned opposing response directions, with the proportion of genes in each direction set to 50:50, 60:40, or 75:25 (two pathways per ratio); and (iii) four null pathways that carried no pathway-specific spatial signal and served as negative controls. Coherent and heterogeneous pathways were constructed in matched sets: pathways under comparison shared the same spatial pattern, the same proportion of responsive genes, and the same effect-size distribution, and differed only in the composition of response directions, thereby isolating directional heterogeneity as the sole variable of interest. Within each heterogeneous pathway, genes assigned to the two opposing directions were additionally matched in baseline expression and effect magnitude. As a result, at the 50:50 ratio the contributions of the two gene sets to the pathway-level mean canceled exactly, whereas at the 60:40 and 75:25 ratios they canceled only partially, producing graded levels of directional imbalance.

For each spatially responsive pathway $p$—that is, each coherent or heterogeneous pathway—we defined a ground-truth spatial activity map $A_i^p \in [0,1]$, denoting the relative activity of pathway $p$ at spatial location $i$. We considered four canonical spatial patterns: regional, in which activity is elevated within a defined region; gradient, in which activity varies continuously along a spatial axis; boundary, in which activity peaks at the interface between two regions and decays with distance; and multi-hotspot, in which activity is elevated in multiple disconnected regions. Each responsive gene $g$ in pathway $p$ inherited the shared spatial map, $H_{ig} = A_i^p$, whereas background SVGs were assigned independently generated patterns $H_{ig}$ unrelated to the test pathways. The response direction of gene $g$ was encoded as $d_g \in \{-1, +1\}$.

#### Raw-count simulations

Raw counts $Y_{ig}$ for gene $g$ at spatial location $i$ were drawn from a negative binomial distribution,

$$Y_{ig} \sim \mathrm{NB}(\mu_{ig}, \theta_g), \qquad \mathrm{Var}(Y_{ig}) = \mu_{ig} + \frac{\mu_{ig}^2}{\theta_g},$$

with mean structure

$$\mu_{ig} = s_i\, m_g\, \bigl(1 + d_g\, b_g\, H_{ig}\bigr),$$

where $m_g > 0$ is the gene-specific baseline expression level, $s_i$ the location-specific sequencing-depth factor (fixed at 1 throughout), $\theta_g$ the gene-specific overdispersion parameter, and $b_g \in (0,1)$ the spatial effect magnitude. Baseline expression levels $m_g$ and dispersions $\theta_g$ were sampled from empirical gene-level distributions estimated from reference spatial transcriptomic data, and spatial effect magnitudes for responsive genes were drawn from $b_g \sim \mathrm{Beta}(\alpha_b, \beta_b)$.

#### Gaussian simulations

To emulate expression measurements on a transformed, approximately Gaussian scale, we generated

$$X_{ig} \sim \mathcal{N}\left(a_g + d_g\, \delta_g\, H_{ig},\; \sigma_g^2\right),$$

where $a_g$ is the baseline expression level, $\sigma_g^2$ the residual variance, and $\delta_g > 0$ the absolute spatial effect magnitude; $a_g$ and $\sigma_g^2$ were likewise sampled from empirical reference distributions. To place effect strength on a common scale across genes, we defined the standardized spatial effect $r_g = \delta_g / \sigma_g$ and drew $r_g \sim \mathrm{Beta}(\alpha_r, \beta_r)$.



### Model detail










