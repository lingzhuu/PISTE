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
SpaPath uses a Gaussian linear mixed model (LMM) for normalized expression data, 正态误差描述观测波动.
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

To evaluate the performance of these methods on real data, we analyzed a human dorsolateral prefrontal cortex (DLPFC) Visium dataset because its six cortical layers and white matter were manually annotated and the molecular functions of each layer are well documented. [Maynard et al., 2021] We applied SpaPath, GSVA, STpathway, and GSDensity to section 151673, which contains 3,639 spots. We tested a prespecified panel of 30 MSigDB gene sets: 9 for synaptic transmission (including glutamatergic, GABAergic, and NMDA/AMPA receptor programs), 3 for voltage-gated ion-channel activity, 1 for oxidative phosphorylation, 3 for myelination, 4 for lipid metabolism, and 2 for dendrite development and morphogenesis, together with 8 reference sets for basal RNA/protein maintenance and DNA repair. Sets were selected for appropriate size (15–229 detected genes in the transcriptome of the selected tissue section) and minimal redundancy (mean pairwise Jaccard index = 0.010; the only substantial overlap was between the two intentionally related dendrite sets, J = 0.60).

[Supplementary Table 2. The 30 prespecified MSigDB gene sets for the DLPFC analysis.](benchmark/analysis/DLPFC/data/geneset_30.xlsx)

We compared the layer-versus-rest results of the four methods. Using the WM-versus-rest comparison as an example, most methods recovered the strong expected signals: lipid metabolism and myelin-related pathways were inferred as up-regulated in WM, consistent with myelin being 70–80% lipid by dry weight (Baumann and Pham-Dinh, 2001; Simons and Nave, 2015), whereas synaptic transmission pathways were inferred as down-regulated, consistent with the signaling-dominated energy budget of gray matter (Harris, Jolivet and Attwell, 2012). The methods diverged on weaker signals. STpathway's region dCor ranked lipid and myelin pathways as the most regionally divergent but spanned only 0.83 to 0.99 across all 30 pathways, offering limited discrimination between strong and weak signals. GSDensity reported only 20 of the 30 gene sets after its coherence filter, and this reported number was unstable with respect to panel composition: when we restricted the analysis to a 20-set subset, only 12 sets were reported. The reason is that the coherence screen is applied to the panel as a whole and discards sets whose high-activity spots are spatially diffuse, so whether a set is retained depends on which other sets are in the panel. Moreover, the discarded sets were not random. They included voltage-gated sodium channel activity, myelin maintenance, cholesterol biosynthesis, sphingolipid metabolism, and very-long-chain fatty acid metabolism, several of which all other methods identified as strongly WM-up-regulated. Even among the reported sets, because the GSDensity test compares proportions of high-activity spots, pathways with sparse high-activity spots were judged non-significant even when regional differences were strong, like glutamatergic synaptic transmission (p = 0.27). For GSVA, we tested WM-versus-rest differences with a spatial block bootstrap (1,999 resamples over 85 blocks of side 0.10 in normalized coordinates, preserving within-block scores and labels), with Benjamini–Hochberg correction across the 30 tests; plain permutation tests with the same correction gave identical conclusions. These tests confirmed the myelin and lipid up-regulation (q = 0.0021) and additionally supported two basal DNA-repair sets (base excision repair, q = 0.019; mismatch repair, q = 0.0021) but gave no support for the weaker synaptic, ion-channel, and oxidative phosphorylation signals.

In contrast, SpaPath recovered all strong signals in the correct direction while remaining sensitive to weaker, biologically expected ones (Fig. X). When the 22 brain-related gene sets were positioned by their estimated effect proportion (π) and posterior evidence, they resolved into two distinct, biologically coherent groups (Fig. X). All seven myelin- and lipid-related sets were assigned to white matter: three with evidence > 0.99 (myelin sheath, sphingolipid metabolism, and monocarboxylic acid transport), three with evidence > 0.95 (axon ensheathment, myelin maintenance, and very-long-chain fatty acid metabolism), and one, cholesterol biosynthesis, with evidence ≥ 0.90 (0.912), consistent with reduced cholesterol synthesis in the adult brain (Dietschy & Turley, 2004). Collectively, these sets recapitulate the oligodendrocyte lipid-synthesis program that underlies myelination. Conversely, all synaptic transmission (n = 9), voltage-gated channel (n = 3), oxidative phosphorylation, and dendrite development/morphogenesis (n = 2) sets were assigned to gray matter (cortical layers L1–L6) with evidence ≤ 0.05, reflecting three defining properties of this compartment: dense synaptic connectivity and dendritic arborization, neuronal excitability, and high energetic demand. Notably, oxidative phosphorylation received strong support for lower WM activity (evidence = 2.5 × 10⁻⁵) despite a modest effect size, consistent with a recent human brain map of mitochondrial respiratory capacity showing that gray matter contains over 50% more mitochondria than white matter (Mosharov et al., 2025), whereas other benchmark methods could not distinguish this signal from noise.

We additionally noted that several pathways representing basal cellular processes also showed a tendency toward higher inferred activity in WM, such as mRNA splicing, mismatch repair, nucleotide excision repair, and base excision repair (evidence > 0.9). We hypothesized that these inferred pathway differences arise from regional differences in cell composition. To test this hypothesis, we estimated cell composition by RCTD deconvolution (Cable et al., 2022) using the LIBD spatialDLPFC snRNA-seq reference (Huuki-Myers et al., 2024), which showed a higher oligodendrocyte proportion in WM than in cortical layers (Fig. Xa). To quantify the contribution expected from cell composition, we derived cell-type-specific pathway activities from the same snRNA-seq reference by GSVA (Fig. Xb) and weighted these reference scores by the cell-type proportions. Across all 30 pathways, the resulting WM-versus-rest score differences correlated with the π estimates (Pearson r = 0.872, R² = 0.760; Fig. Xc).

[Figure 1. SpaPath recovers regional pathway activity differences in DLPFC](benchmark/analysis/DLPFC/figure/celltype.png)

These results indicate that a substantial fraction of the inferred WM–gray matter differences in basal cellular pathways can be explained by regional differences in cell composition. These considerations motivated us to distinguish composition-driven pathway differences from those that persist after accounting for cell composition. To this end, we extended SpaPath with a composition covariate. For each spot, we used its deconvolved cell-type proportions to compute the expected expression of each gene and included this quantity as a covariate with a gene-specific coefficient (Methods). The regional effects therefore capture pathway differences beyond those predicted by cell composition.

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
[Integrative spatial omics reveals distinct tumor-promoting multicellular niches and immunosuppressive mechanisms in Black American and White American patients with TNBC. (Zhu et al.)](https://www.nature.com/articles/s41467-025-61034-3)  [data](https://zenodo.org/records/15247064) 
> This dataset contains the 10X Visium Spatial Transcriptomic data of TNBC patients. There are two cohorts. Baylor Scott and White (BSW) cohort: `10x.visium.tar.gz`, containing 10 patients with TNBC from Baylor Scott and White affiliated Hospital. Each sample is made of Space Ranger processed spot-separated gene expression data (processed to HDF5 AnnData file). There are also H&E images, and spot coordinate files available. 
As well in our paper, we integrated a published racial TNBC cohort for deriving some of initial results in the paper. This refers to the Bassiouni et al (Cancer Research) paper in Carpten's group. `GSM_giotto_processed.tar.gz` refers to this dataset, which we deposit here. The data were normalized by Giotto using standard procedure.

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


**Spatial hierarchical model.** We model spatial variation in gene expression with a replicated stochastic partial differential equation (SPDE) framework, fitted via INLA. Each gene $g = 1, \ldots, K$ has an independent spatial field

$$w_g(\cdot) \sim_{iid} \mathrm{GP}\left(0,\; \sigma_s^2 C_{\kappa_0}\right),$$

where all genes share one triangulation mesh, a fixed range parameter $\kappa_0$, and a single global spatial variance $\sigma_s^2$. Each field spans the entire tissue and is not re-initialized at region boundaries. We specify gene-specific spatial fields rather than a shared field $w(s)$ across genes, and assume a single global spatial variance $\sigma_s^2$, so that spatial structure is gene-specific in pattern but common in scale.

Each region $r$ carries a latent pathway state

$$Z_r \sim \mathcal{N}\left(\mu_r,\; \tau_r^2\right),$$

capturing region-level activity complementary to the gene-specific fields. 
Notably, all variance components are shared across genes, so the number of hyperparameters in the INLA formulation grows only with the number of regions $R$ and is independent of the number of genes $K$.






