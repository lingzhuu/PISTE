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

The central question addressed by SpaPath is whether a pathway is more active in one tissue region than in another. Answering it is complicated by two features of spatial transcriptomics data: the genes of a pathway may respond with different magnitudes or even in opposite directions, and neighboring spots tend to have similar expression, so that spots cannot be treated as independent observations (Fig. X). SpaPath takes as input normalized expression, spot coordinates, region labels and a gene set, optionally together with deconvolved cell-type proportions and a single-cell or single-nucleus RNA-seq reference. It then fits a single statistical model, a Gaussian linear mixed model, jointly to all genes of the pathway across all regions (Methods). In this model, each region is described by the average expression level of the pathway's genes and by how widely the individual genes spread around that average. In addition, each gene has its own smooth spatial pattern across the whole tissue, which captures gradual trends in expression so that they are not mistaken for differences between regions. Because these spatial patterns share a common scale across genes, the number of parameters to estimate grows only with the number of regions and not with the size of the pathway. The model is fitted in a Bayesian framework with integrated nested Laplace approximation (INLA), a fast and deterministic alternative to sampling-based (Markov chain Monte Carlo) inference (REF).

From the fitted model, SpaPath summarizes each region by the expression of a randomly chosen gene of the pathway and defines the directional score $\pi_{rs}$ as the probability that this expression in region $r$ exceeds that in region $s$. A pathway whose genes shift consistently in one direction therefore receives a score near 0 or 1, whereas one whose genes move in opposite directions receives a score near 0.5. Because the model parameters are estimated with uncertainty, SpaPath reports two summaries: the π score, which is the expected value of this probability, and the posterior evidence, which is the probability that the pathway is more active in region $r$ than in region $s$. To compare one region with the rest of the tissue, the scores against all other regions are averaged, weighted by the number of spots in each region.

Finally, because regional differences can arise either from differences in cell-type composition or from changes within the same cell types, SpaPath offers two modes. The observed mode reports the total regional signal, whereas the composition-adjusted mode additionally accounts for the expression expected from each spot's cell-type composition and thus reports the signal beyond composition (Methods). Since both modes share all other parameters, their results can be compared directly.

### Benchmarking on simulated spatial data

- GSVA etc single cell pathway method in ST data
- STpathway [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503753/)
  > employing the Brownian distance correlation (dCor) test in a novel set-up to evaluate the dependency of two random vectors with different dimensionality.
- GSDensity [paper](https://www.nature.com/articles/s41467-023-44206-x)
  > Here, we present GSDensity, a graph-modeling approach that allows users to obtain pathway-centric interpretation and dissection of single-cell and spatial transcriptomics (ST) data without performing clustering. Using pathway gene sets, we show that GSDensity can accurately detect biologically distinct cells and reveal novel cell-pathway associations ignored by existing methods. 
  
We generated two simulated datasets with different random seeds to benchmark SpaPath against three representative pathway-analysis workflows: GSVA, the gene-set scoring approach popularized in single-cell RNA sequencing, tested here with a spatial block bootstrap region test; STpathway, a distance-correlation (dCor) test for region-associated pathways; and GSDensity, evaluated by its published region enrichment test on dichotomized spot activity. Each dataset was built on the real spot coordinates of DLPFC section 151673 and comprised 50 non-overlapping 100-gene sets spanning five spatial activity patterns, with the number K of genes responding to the activity map prespecified for each gene set (Methods).

Across the ten groups defined by K and effect distribution, the existing methods did not provide discriminative region-level evidence. 
GSDensity's region test reflected strong simulated pathway activity: as activity increased, the difference in high-activity spot proportions between target and rest increased whereas the region p-value decreased. The pass rate of GSDensity's coherence pre-filter increased with the prespecified pathway activity; however, of the 100 prespecified gene sets, only 35 passed this filter, which checks whether a set's genes cluster in the embedding space and whether its PAL activity scores are reliable. Following the original publication, we report its region test only for those passing sets. Among these sets, the median region p in the broad, weak-effect setting was 3.6 × 10⁻⁵ at K = 31 to 60 and 0.011 at K = 61 to 100, but not in the narrow, strong-effect setting, where the median p was 0.247 at K = 31 to 60 and 0.082 at K = 61 to 100. Additionally, in the high-K groups of the narrow, strong-effect setting, fewer than 1% of spots were labeled highly active in either target or rest, and the region p-value was therefore not meaningful.
STpathway's region dCor decreased with increasing K, indicating that it can reflect strong regional differences. However, its group medians spanned only 0.94 to 0.98, with an overall range of 0.87 to 0.99, offering limited discrimination of pathway activity. GSVA's region score difference increased with K, and its spatial block bootstrap test gave BH-adjusted median q-values below 0.05 in all groups with K of 11 or more, but it failed to discriminate activity strength among weakly responsive pathways. Furthermore, under the domain and interface patterns, the two pattern types with the largest target regions, the test returned extremely small p-values even when the prespecified pathway activity was low or absent: the K = 0 domain pathways attained p = 0.0005 and remained significant after BH correction in the datasets from both random seeds. This indicates that under these spatial patterns the test yields strong significance regardless of the prespecified activity level and cannot separate pattern-driven structure from true pathway responses.

SpaPath's posterior evidence instead increased monotonically as the prespecified pathway activity increased. The median π score rose from 0.37 in the K = 0 to 4 group to 0.81 in the K = 61 to 100 group in the broad, weak-effect setting, and from 0.40 to 0.75 in the narrow, strong-effect setting. Posterior evidence reached about 0.98 already in the K = 11 to 30 group and was close to 1 in the highest-K groups. Importantly, the π score and evidence also recovered the relative activity strength among weak pathways, where the existing methods could not: between the K = 0 to 4 and K = 5 to 10 groups, the median π score rose from 0.37 to 0.48 in the broad, weak-effect setting and from 0.40 to 0.51 in the narrow, strong-effect setting, and the median evidence rose from 0.15 to 0.39 and from 0.24 to 0.60, respectively. Thus, pathways with little or no prespecified activity received weak posterior support, yet the ordering of their activity levels was preserved. This pattern was consistent across all five spatial activity patterns. 

[Table 1. Median region-level statistics of the four methods across the 100 simulated gene sets.](benchmark/analysis/simulation_data/result/simulation_summary.xlsx)

[Supplementary Table 1. Per-pathway simulation results for the simulation dataset.](benchmark/analysis/simulation_data/result/simulation_per_pathway_result.xlsx)
### Layer-resolved pathway activity in the human DLPFC

To evaluate the performance of these methods on real data, we analyzed a human dorsolateral prefrontal cortex (DLPFC) Visium dataset because its six cortical layers and white matter were manually annotated and the molecular functions of each layer are well documented. [Maynard et al., 2021] We applied SpaPath, GSVA, STpathway, and GSDensity to section 151673, which contains 3,639 spots. We tested a prespecified panel of 20 MSigDB gene sets: 5 for synaptic transmission (including glutamatergic and NMDA/AMPA receptor programs), 2 for ion-channel activity, 1 for oxidative phosphorylation, 3 for myelination, and 4 for lipid metabolism, together with 5 reference sets for basal RNA/protein maintenance and DNA repair. Sets were selected for appropriate size (16–210 detected genes in the transcriptome of the selected tissue section) and minimal redundancy (mean pairwise Jaccard index = 0.008; maximum = 0.19).

[Supplementary Table 2. The 20 prespecified MSigDB gene sets for the DLPFC analysis.](benchmark/analysis/DLPFC/data/geneset_20.xlsx)

We compared the layer-versus-rest results of the four methods. Using the WM-versus-rest comparison as an example, most methods recovered the strong expected signals: lipid metabolism and myelin-related pathways were inferred as up-regulated in WM, consistent with myelin being 70–80% lipid by dry weight (Baumann and Pham-Dinh, 2001; Simons and Nave, 2015), whereas synaptic transmission pathways were inferred as down-regulated, consistent with the signaling-dominated energy budget of gray matter (Harris, Jolivet and Attwell, 2012). The methods diverged on weaker signals. STpathway's region dCor ranked lipid and myelin pathways as the most regionally divergent but spanned only 0.83 to 0.99 across all pathways, offering limited discrimination between strong and weak signals. GSDensity reported 11 of the 20 gene sets after its coherence filter, and because its test compares proportions of high-activity spots, pathways with sparse high-activity spots were judged non-significant even when regional differences were strong, like glutamatergic synaptic transmission (p = 0.27). GSVA's spatial block bootstrap tests confirmed the myelin and lipid upregulation (q = 0.0017) but gave no support for weaker signals. 

In contrast, SpaPath recovered all strong signals in the correct direction while remaining sensitive to weaker, biologically expected ones (Fig. X). Here the target region is WM and the rest consists of the gray-matter cortical layers L1–L6, so the π score measures how often a WM pathway state exceeds the spot-weighted gray-matter state, and the posterior evidence is the probability that it does (π > 0.5): evidence near 1 indicates higher activity in WM and evidence near 0 higher activity in gray matter. When the 15 brain-related gene sets were positioned by their WM-versus-rest π score and posterior evidence, they resolved into two distinct, biologically coherent groups (Fig. X). All seven myelin- and lipid-related sets were inferred to be more active in WM: three with evidence > 0.99 (myelin sheath, sphingolipid metabolism, and monocarboxylic acid transport), three with evidence > 0.95 (axon ensheathment, myelin maintenance, and very-long-chain fatty acid metabolism), and one, cholesterol biosynthesis, with evidence ≥ 0.90 (0.912), consistent with reduced cholesterol synthesis in the adult brain (Dietschy & Turley, 2004). Collectively, these sets recapitulate the oligodendrocyte lipid-synthesis program that underlies myelination. Conversely, all synaptic transmission (n = 5), voltage-gated channel (n = 2), and oxidative phosphorylation sets were inferred to be more active in gray matter (cortical layers L1–L6), reflecting three defining properties of this compartment: dense synaptic connectivity, neuronal excitability, and high energetic demand, respectively. Several of these sets had very low evidence (e.g., synaptic vesicle exocytosis, 3.2 × 10⁻⁷; glutamatergic synaptic transmission, 4.4 × 10⁻⁷). Notably, oxidative phosphorylation received strong support for lower WM activity (evidence = 2.5 × 10⁻⁵) despite a modest effect size, consistent with gray matter containing over 50% more mitochondria than white matter (Mosharov et al., 2025), whereas other benchmark methods could not distinguish this signal from noise.

We additionally noted that several pathways representing basal cellular processes also showed a tendency toward higher inferred activity in WM, such as mRNA splicing, nucleotide excision repair, base excision repair, and mismatch repair (evidence > 0.9). We hypothesized that these inferred pathway differences arise from regional differences in cell composition. To test this hypothesis, we examined cell composition by cell2location deconvolution, which estimated a higher oligodendrocyte proportion in WM than in cortical layers(Fig. Xa). To quantify the contribution expected from cell composition, we derived cell-type-specific pathway activities from the human MTG snRNA-seq reference by GSVA (Fig. Xb) and weighted these reference scores by the cell-type proportions. The resulting WM-versus-rest score differences were positive for all five selected pathways and correlated with the π scores (Pearson r = 0.652, R² = 0.425; Fig. Xc).

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

For the normalized data using SCTransform, we model the transformed expression values using a Gaussian distribution.


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
### SpaPath method

#### Design principles

SpaPath asks whether a pathway is more active in one tissue region than in another. Because the genes of a pathway can respond with different magnitudes and even in opposite directions, we do not test a single aggregated score; instead, we model all genes of the pathway jointly and describe each region by the distribution of its gene-level states. Spatial transcriptomics data are, moreover, spatially autocorrelated, so treating spots as independent overstates the evidence for regional differences. SpaPath therefore combines regional pathway effects with gene-specific spatial random fields in one Gaussian linear mixed model, fits it with integrated nested Laplace approximation (INLA), which approximates Bayesian posterior distributions deterministically and much faster than Markov chain Monte Carlo sampling (REF), and converts the posterior into a directional score with an accompanying posterior evidence. Optionally, the model can also adjust for differences in cell composition between spots.

#### Data input

SpaPath requires normalized expression of the pathway's genes, spot coordinates, region labels and a triangulation mesh of the tissue; spot-level covariates such as batch are optional. Normalization and quality control are performed beforehand, and in this study we used SCTransform Pearson residuals. The composition-adjusted mode additionally requires deconvolved cell-type proportions for each spot and a single-cell or single-nucleus RNA-seq reference of cell-type mean expression.

#### Model specification

For spot $i=1,\ldots,N$, gene $g=1,\ldots,K$ and region $r(i)\in\{1,\ldots,R\}$, SpaPath fits

$$
Y_{ig}=X_i^T\beta+\mu_{r(i)}+b_{g,r(i)}+w_g(s_i)+\epsilon_{ig},
\qquad
b_{gr}\mid\tau_r^2\sim N(0,\tau_r^2),
\qquad
\epsilon_{ig}\sim N(0,\sigma_\epsilon^2),
$$

where $\beta$ are covariate effects shared by all genes and $\mu_r$ is the mean pathway level in region $r$. The gene effects $b_{gr}$ are conditionally independent across genes, and each region has its own variance $\tau_r^2$, which measures how consistently the pathway's genes behave in that region. Each gene has an independent Matérn field $w_g(\cdot)\overset{iid}{\sim}\mathrm{GP}(0,\sigma_s^2C_{\kappa_0})$, represented through a replicated stochastic partial differential equation (SPDE) on a shared mesh, which makes the spatial fields computationally tractable (REF). All fields share the scale $\kappa_0=\sqrt8/\rho_0$, where $\rho_0$ is a user-specified practical range, and a single variance $\sigma_s^2$, so spatial structure is gene-specific in pattern but common in scale; furthermore, the fields span the whole tissue and do not restart at region boundaries. The model thus has $R+2$ hyperparameters regardless of the number of genes.

#### Adjustment for cell composition

In the composition-adjusted mode, let $\omega_{ik}$ be the proportion of cell type $k$ in spot $i$ and $\nu_{kg}$ the depth-normalized mean expression of gene $g$ in cell type $k$ in the reference. The log expression expected from composition alone, $o_{ig}=\log(\sum_k\omega_{ik}\nu_{kg}+c_{\mathrm{ref}})$ with a pseudo-count $c_{\mathrm{ref}}$ on the scale of $\nu_{kg}$, is centered within each gene, $\tilde o_{ig}=o_{ig}-N^{-1}\sum_{i'}o_{i'g}$, and the model gains the term $\lambda_g\tilde o_{ig}$. The coefficient $\lambda_g$ is estimated separately for each gene, since the reference and the spatial data differ in scale and platform. Because centering leaves each gene's overall level in $b_{gr}$, all other parameters keep their meaning, and the two modes can be compared directly: the observed mode reports the total regional signal, whereas the adjusted mode reports the part not explained by composition. This adjustment treats $\omega_{ik}$ and $\nu_{kg}$ as known and assumes that a cell type has the same expression profile in every region.

#### Priors and posterior inference

Fixed effects, including the regional means, have $N(0,10^2)$ priors, and $\lambda_g\overset{iid}{\sim}N(0,10^2)$ with a fixed variance, so the composition term adds no hyperparameters. The regional log precisions follow $\log(1/\tau_r^2)\sim N(0,3^2)$, the residual precision has a penalized-complexity prior with $P(\sigma_\epsilon>2)=0.05$, and $\log\sigma_s\sim N(0,1)$. INLA integrates over the hyperparameters with the central composite design (CCD); because the number of design points grows quickly with the number of hyperparameters, SpaPath warns when it exceeds a user-set limit (10 by default). For each design point $k$ with normalized weight $q_k$, INLA gives a Gaussian posterior $N(m_k,\Sigma_k)$ for the regional means together with the variances $\tau_{r,k}^2$, and the joint posterior of the regional means is approximated by the mixture of these Gaussians.

#### Directional pathway score and posterior evidence

Each region is summarized by the state of a randomly chosen pathway gene, $Z_r\sim N(\mu_r,\tau_r^2)$, and the directional score $\pi_{rs}=P(Z_r>Z_s)=\Phi\{(\mu_r-\mu_s)/(\tau_r^2+\tau_s^2)^{1/2}\}$ is the probability that the pathway state in region $r$ exceeds that in region $s$. We report its posterior mean, the π score, which has the closed form

$$
\bar\pi_{rs}=\sum_k q_k\,\Phi\!\left(\frac{m_{k,r}-m_{k,s}}{\sqrt{\tau_{r,k}^2+\tau_{s,k}^2+\Sigma_{k,rr}+\Sigma_{k,ss}-2\Sigma_{k,rs}}}\right),
$$

together with the posterior evidence $P(\pi_{rs}>0.5\mid Y)$, which equals $P(\mu_r>\mu_s\mid Y)$. To compare a region with the rest of the tissue, we average the pairwise scores with spot-count weights, $\pi_{r,\mathrm{rest}}=\sum_{s\neq r}a_{rs}\pi_{rs}$ with $a_{rs}\propto n_s$. With two regions the evidence is computed analytically; with more regions it is computed by deterministic quasi-Monte Carlo integration over the mixture, doubling the number of Sobol points until two consecutive refinements change the evidence by less than 0.002, so no posterior sampling is involved. Regional mean contrasts and variance components are reported as diagnostics only.

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

For each pathway $p$, $K_p$ is the number of the responsive members among the 100 genes, sampled from five ranges: 0–4, 5–10, 11–30, 31–60, or 61–100; the remaining $100-K_p$ members are not responsive. Response direction is drawn before count generation, with $d_g=+1$ for target-upregulated genes and $d_g=-1$ for target-downregulated genes. Genes in weaker response group are more likely to have their direction reversed, with reversal probabilities of 30–50%, 15–30%, and 10–15% in the weak, medium, and strong groups, respectively. $E_g$ denotes the effect size of an individual responsive gene.

Each pathway is assigned an activity map $A_i^{(p)}$, normalized to $[0,1]$, from one of five classes: discrete domains, smooth gradients, radial or niche-centered patterns, localized hotspots, or interface-specific patterns. Domain activity is piecewise constant across $D_p$ spatial regions; gradients vary smoothly across the tissue; radial activity decays with distance from a niche center; hotspot activity is localized; and interface activity follows
$A_i^{(p)} = \mathbb{1}(h_i \leq r)\exp(-h_i/\lambda)$, where $h_i$ is the distance to the relevant boundary. Responsive genes receive their pathway's activity map through $H_{ig}=A_i^{(p)}$.

Raw read counts $Y_{ig}$ are simulated from a negative binomial distribution:

$$
Y_{ig} \sim \mathrm{NB}(\mu_{ig}, \theta_g), \qquad
\mathrm{Var}(Y_{ig}) = \mu_{ig} + \frac{\mu_{ig}^2}{\theta_g},
$$

$$
\log \mu_{ig} = \log c_i + \log m_g + d_g b_g H_{ig}.
$$

Here $c_i$ is the spot-specific depth factor and $(m_g, \theta_g)$ is a jointly sampled empirical baseline-expression and dispersion pair. 

GSDensity receives a Seurat object built from raw counts and uses the author function's internal NormalizeData() step. SpaPath and STpathway receive the SCTransform Pearson residuals.
