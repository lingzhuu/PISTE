## Abstract


## Introduction
#### Background
Spatial transcriptomics has transformed tissue biology by simultaneously capturing genome-wide gene expression and its native spatial coordinates, thereby preserving the positional context that bulk and dissociated single-cell assays discard. (REF) Because biological processes are intrinsically spatially organized, this technology has enabled discoveries that are inaccessible to non-spatial assays, including spatially resolved differential expression, the delineation of tissue niches and domains, and the mapping of ligand–receptor-mediated cell–cell communication within intact tissue architecture. (REF) As both resolution and throughput continue to increase, spatial transcriptomics is becoming a routine assay in biomedical research.

Yet biological interpretation rarely rests on individual genes. Phenotypic states and disease processes are executed by coordinated gene programs—pathways and regulatory modules whose constituent genes act in concert rather than in isolation. This is particularly relevant for spatial data, because pathway-level readouts are more robust to the technical noise, dropouts, and cell-type admixture that characterize spot-level measurements. (REF) Moreover, within a single pathway, individual genes may exert opposite directional effects on the biological process of interest; aggregating gene-level evidence into a signed pathway score is often the only way to obtain a coherent picture of the underlying activity. 
#### Problem
Existing pathway approaches have two limitations. First, spatial transcriptomics data are spatially autocorrelated, because neighbouring spots share cell types and microenvironments, yet existing pathway analyses largely treat spots as independent observations. Ignoring this dependence underestimates the uncertainty of regional differences in pathway activity, which overstates statistical significance and inflates false-positive rates. (REF) Second, existing approaches have limited power to detect pathways with weak regional differences, and such pathways are therefore often missed even when they are biologically relevant. (REF)
#### Method
We developed SpaPath：a linear mixed model framework to quantify pathway-level expression differences across tissue regions. 每个 pathway 在所有 regions 中联合拟合，不再把均值差异和方差差异作为两个主要生物学输出。它们保留为模型参数和诊断量，并共同决定最终的 probability score。
#### Benchmark
Simulation: 并与GSDensity/STPathway/naïve analysis比较。
#### Application
在 spatial brain transcriptomics 数据上识别：cell-type- and region-specific pathway activity
#### Conclusion


## Results 
### Method overview

SpaPath is designed to evaluate whether a pathway is more active in one tissue region than in another. Because the genes of a pathway can respond with different magnitudes and even in opposite directions, SpaPath models all genes of the pathway jointly and describes each region by the distribution of its gene-level expression states. Spatial transcriptomics data are also spatially autocorrelated, so SpaPath assigns each gene a smooth spatial random field that captures the dependence among neighbouring spots. SpaPath describes the expression of each gene in each spot through a regional pathway level shared by all genes in that region, a gene-specific deviation from this level whose spread may differ between regions, and the gene-specific spatial field. These components form a Bayesian hierarchical model, specifically a latent Gaussian model (Methods). The model is fitted with integrated nested Laplace approximation (INLA) (REF), and its posterior is converted into a directional score with an accompanying posterior evidence. As input, SpaPath requires spot coordinates, region labels and a gene set, together with normalized expression in the observed mode, or raw counts, deconvolved cell-type proportions and a single-cell or single-nucleus RNA-seq reference in the composition-adjusted mode.

SpaPath compares two regions $r$ and $s$ through the directional score $\pi_{rs}$, defined as the probability that a gene of the pathway is expressed at a higher level in region $r$ than in region $s$, given the regional pathway levels and the variability among genes. Its posterior mean, which we term the π score, lies between 0 and 1. Values near 1 indicate that the pathway is consistently more active in region $r$, values near 0 indicate the opposite, and values near 0.5 indicate no consistent direction, as happens when the genes of a pathway change in opposite directions. SpaPath also reports the posterior evidence $P(\pi_{rs}>0.5)$, the posterior probability that the pathway is more active in region $r$ than in region $s$. A region can be compared with the rest of the tissue in the same way.

SpaPath can be run in two modes. The observed mode fits a Gaussian model to normalized expression and quantifies the total difference in pathway activity between regions, including differences that arise because the regions contain different cell types. The composition-adjusted mode fits a negative binomial model to raw counts. It includes the expression of each gene predicted from the cell-type proportions of each spot as an offset, so that the regional difference it reports is the part that cannot be explained by cell composition. We recommend running both modes, because together they distinguish pathways whose regional differences are driven by cell composition, which are supported only in the observed mode, from pathways that change within cell types, which remain supported after adjustment.

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

To evaluate the performance of these methods on real data, we analyzed a human dorsolateral prefrontal cortex (DLPFC) Visium dataset because its six cortical layers and white matter were manually annotated and the molecular functions of each layer are well documented. [Maynard et al., 2021] We applied SpaPath, GSVA, STpathway, and GSDensity to section 151673, which contains 3,639 spots. We tested a prespecified panel of 20 MSigDB gene sets: 6 for synaptic transmission (including glutamatergic, GABAergic, and NMDA/AMPA receptor programs), 1 for voltage-gated potassium channel activity, 1 for oxidative phosphorylation, 3 for myelination, and 4 for lipid metabolism, together with 5 reference sets for basal RNA/protein maintenance and DNA repair. Sets were selected for appropriate size (16–210 detected genes in the transcriptome of the selected tissue section) and minimal redundancy (mean pairwise Jaccard index = 0.008; maximum = 0.19).

[Supplementary Table 2. The 20 prespecified MSigDB gene sets for the DLPFC analysis.](benchmark/analysis/DLPFC/data/geneset_20.xlsx)

We compared the layer-versus-rest results of the four methods. Using the WM-versus-rest comparison as an example, most methods recovered the strong expected signals: lipid metabolism and myelin-related pathways were inferred as up-regulated in WM, consistent with myelin being 70–80% lipid by dry weight (Baumann and Pham-Dinh, 2001; Simons and Nave, 2015), whereas synaptic transmission pathways were inferred as down-regulated, consistent with the signaling-dominated energy budget of gray matter (Harris, Jolivet and Attwell, 2012). The methods diverged on weaker signals. STpathway's region dCor ranked lipid and myelin pathways as the most regionally divergent but spanned only [TBD] to [TBD] across all 20 pathways, offering limited discrimination between strong and weak signals. GSDensity reported only [TBD] of the 20 gene sets after its coherence filter, and this reported number was unstable with respect to panel composition: [TBD: comparison panel and number of sets reported]. The reason is that the coherence screen is applied to the panel as a whole and discards sets whose high-activity spots are spatially diffuse, so whether a set is retained depends on which other sets are in the panel. Moreover, the discarded sets were not random. They included [TBD: discarded sets], several of which all other methods identified as strongly WM-up-regulated. Even among the reported sets, because the GSDensity test compares proportions of high-activity spots, pathways with sparse high-activity spots were judged non-significant even when regional differences were strong, like glutamatergic synaptic transmission (p = 0.27). For GSVA, we tested WM-versus-rest differences with a spatial block bootstrap (1,999 resamples over 85 blocks of side 0.10 in normalized coordinates, preserving within-block scores and labels), with Benjamini–Hochberg correction across the 20 tests; plain permutation tests with the same correction gave identical conclusions. These tests confirmed the myelin and lipid up-regulation (q = [TBD]) [TBD: basal sets also supported, if any] but gave no support for the weaker synaptic, ion-channel, and oxidative phosphorylation signals.

In contrast, SpaPath recovered all strong signals in the correct direction while remaining sensitive to weaker, biologically expected ones (Fig. X). Here the target region is WM and the rest consists of the gray-matter cortical layers L1–L6, so the π score measures how often a WM pathway state exceeds the spot-weighted gray-matter state, and the posterior evidence is the probability that it does (π > 0.5): evidence near 1 indicates higher activity in WM and evidence near 0 higher activity in gray matter. When the 15 brain-related gene sets were positioned by their WM-versus-rest π score and posterior evidence, they resolved into two distinct, biologically coherent groups (Fig. X). All seven myelin- and lipid-related sets were inferred to be more active in WM: three with evidence > 0.99 (myelin sheath, sphingolipid metabolism, and monocarboxylic acid transport), three with evidence > 0.95 (axon ensheathment, myelin maintenance, and very-long-chain fatty acid metabolism), and one, cholesterol biosynthesis, with evidence ≥ 0.90 (0.912). Collectively, these sets recapitulate the oligodendrocyte lipid-synthesis program that underlies myelination. Conversely, all synaptic transmission (n = 6), voltage-gated potassium channel, and oxidative phosphorylation sets were inferred to be more active in gray matter (cortical layers L1–L6) with evidence ≤ 0.05, reflecting three defining properties of this compartment: dense synaptic connectivity, neuronal excitability, and high energetic demand. Notably, oxidative phosphorylation received strong support for lower WM activity (evidence = 2.5 × 10⁻⁵) despite a modest effect size, consistent with a recent human brain map of mitochondrial respiratory capacity showing that gray matter contains over 50% more mitochondria than white matter (Mosharov et al., 2025), whereas other benchmark methods could not distinguish this signal from noise.

We additionally noted that several pathways representing basal cellular processes also showed a tendency toward higher inferred activity in WM, such as mRNA splicing and nucleotide excision repair (evidence > 0.9). We hypothesized that these inferred pathway differences arise from regional differences in cell composition. To test this hypothesis, we estimated cell composition by RCTD deconvolution (Cable et al., 2022) using the LIBD spatialDLPFC snRNA-seq reference (Huuki-Myers et al., 2024), which showed a higher oligodendrocyte proportion in WM than in cortical layers (Fig. Xa). To quantify the contribution expected from cell composition, we derived cell-type-specific pathway activities from the same snRNA-seq reference by GSVA (Fig. Xb) and weighted these reference scores by the cell-type proportions. Across all 20 pathways, the resulting WM-versus-rest score differences correlated with the π scores (Pearson r = [TBD], R² = [TBD]; Fig. Xc).

[Figure 1. SpaPath recovers regional pathway activity differences in DLPFC](benchmark/analysis/DLPFC/figure/celltype.png)

These results indicate that a substantial fraction of the inferred WM–gray matter differences in basal cellular pathways can be explained by regional differences in cell composition. These considerations motivated us to distinguish composition-driven pathway differences from those that persist after accounting for cell composition. To this end, we extended SpaPath with a composition covariate. For each spot, we used its deconvolved cell-type proportions to compute the expected expression of each gene and included this quantity as a covariate with a gene-specific coefficient (Methods). The regional effects therefore capture pathway differences beyond those predicted by cell composition.

### Composition-adjusted pathway activity in the DLPFC

Cell composition acts on the expected counts of a spot, which are a mixture of cell-type expression profiles, and this relationship is lost after normalization such as Pearson residuals. The composition-adjusted mode therefore models the raw counts with a negative binomial likelihood and a sequencing-depth offset, and adds a composition term $\lambda_g\log m_{ig}$, where $m_{ig}=\sum_k\omega_{ik}\nu_{kg}$ is the expression of gene $g$ in spot $i$ expected from its cell-type proportions $\omega_{ik}$ and the reference expression $\nu_{kg}$ (Methods). The coefficient $\lambda_g\sim N(1,0.25^2)$ is estimated for each gene rather than fixed at 1, and each gene has its own intercept. The regional levels therefore describe the pathway differences that remain after accounting for cell composition.

We applied this mode to the WM-versus-rest comparison for the same 20 gene sets. The myelin and lipid sets remained more active in WM, with evidence between 0.86 and 0.997, indicating that these differences are not explained by the higher oligodendrocyte proportion of WM alone. The exception was cholesterol biosynthesis, whose π score fell from 0.61 to 0.42. Oxidative phosphorylation remained lower in WM with stronger support (evidence 2.2 × 10⁻⁹). Synaptic vesicle exocytosis, the glutamate release cycle and NMDA receptor activation remained lower in WM (evidence ≤ 0.023), whereas glutamatergic transmission, GABAergic transmission and voltage-gated potassium channel activity moved close to no direction. Among the reference sets, all except mRNA splicing no longer showed higher activity in WM, supporting the hypothesis that their observed-mode differences arise from cell composition.

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

#### Data input

SpaPath requires the expression of the pathway's genes in each spot, spot coordinates and region labels, and spot-level covariates such as batch are optional. In the observed mode the expression input must be normalized beforehand, and we recommend SCTransform (REF) for this purpose. The composition-adjusted mode instead takes raw counts and additionally requires deconvolved cell-type proportions for each spot and a single-cell or single-nucleus RNA-seq reference of cell-type mean expression.

#### Gaussian model for the observed mode

In the observed mode, SpaPath fits a Gaussian model to the normalized expression $Y_{ig}$. For spot $i=1,\ldots,N$, gene $g=1,\ldots,K$ and region $r(i)\in\{1,\ldots,R\}$,

$$
Y_{ig}=X_i^T\beta+\mu_{r(i)}+b_{g,r(i)}+w_g(s_i)+\epsilon_{ig},
\qquad
b_{gr}\mid\tau_r^2\sim N(0,\tau_r^2),
\qquad
\epsilon_{ig}\sim N(0,\sigma_\epsilon^2),
$$

where $\beta$ are covariate effects shared by all genes and $\mu_r$ is the mean pathway level in region $r$. The gene effects $b_{gr}$ are conditionally independent across genes, and each region has its own variance $\tau_r^2$, which measures how consistently the pathway's genes behave in that region. Each gene has an independent Matérn field $w_g(\cdot)\overset{iid}{\sim}\mathrm{GP}(0,\sigma_s^2C_{\kappa_0})$, represented through a replicated stochastic partial differential equation (SPDE) on a shared mesh (REF). All fields share the scale $\kappa_0=\sqrt8/\rho_0$, where $\rho_0$ is a user-specified practical range at which the spatial correlation drops to about 0.1; fixing it rather than estimating it keeps the number of hyperparameters small. The fields also share a single variance $\sigma_s^2$, so spatial structure is gene-specific in pattern but common in scale. Furthermore, the fields span the whole tissue and do not restart at region boundaries. The model thus has $R+2$ hyperparameters regardless of the number of genes.

#### Negative binomial model for the composition-adjusted mode

In the composition-adjusted mode, SpaPath models the raw counts $C_{ig}$. Let $\omega_{ik}$ be the proportion of cell type $k$ in spot $i$ and $\nu_{kg}$ the depth-normalized mean expression of gene $g$ in cell type $k$ in the reference. The log expression expected from composition alone and its within-gene centered version are

$$
o_{ig}=\log\Big(\sum_k\omega_{ik}\nu_{kg}+\varepsilon\Big),
\qquad
\tilde o_{ig}=o_{ig}-\frac{1}{N}\sum_{i'=1}^{N}o_{i'g},
$$

where $\varepsilon$ is a small constant that keeps the logarithm finite for genes not detected in the reference. SpaPath then fits

$$
C_{ig}\sim\mathrm{NB}(\eta_{ig},\phi),
\qquad
\log\eta_{ig}=\log d_i+\tilde o_{ig}+X_i^T\beta+a_g+\mu_{r(i)}+b_{g,r(i)}+w_g(s_i),
$$

where $\eta_{ig}$ is the expected count, $d_i$ is the total count of spot $i$, $\phi$ is an overdispersion parameter shared by all genes and $a_g$ is a gene intercept, with all other terms as in the Gaussian model. Because composition determines expression on the count scale, it enters the log link exactly, and both $\log d_i$ and $\tilde o_{ig}$ are offsets with coefficients fixed at 1. A gene-specific difference in scale between the reference and the spatial data multiplies the expected count by a constant and is absorbed by $a_g$, which also removes each gene's baseline level, so that $b_{gr}$ describes deviations from the regional level rather than overall abundance. The regional means $\mu_r$ thus describe the regional difference that composition does not explain. This adjustment treats $\omega_{ik}$ and $\nu_{kg}$ as known and assumes that a cell type has the same expression profile in every region. Like the Gaussian model, it has $R+2$ hyperparameters.

#### Priors and posterior inference

Fixed effects, including the regional means and the gene intercepts, have $N(0,10^2)$ priors. The regional log precisions follow $\log(1/\tau_r^2)\sim N(0,3^2)$ and $\log\sigma_s\sim N(0,1)$ in both models, and the residual precision of the Gaussian model has a penalized-complexity prior with $P(\sigma_\epsilon>2)=0.05$. INLA integrates over the hyperparameters with the central composite design (CCD) (REF); because the number of design points grows quickly with the number of hyperparameters, SpaPath warns when it exceeds a user-set limit (10 by default). For each design point $k$ with normalized weight $q_k$, INLA gives a Gaussian posterior $N(m_k,\Sigma_k)$ for the regional means together with the variances $\tau_{r,k}^2$, and the joint posterior of the regional means is approximated by the mixture of these Gaussians.

#### Directional pathway score and posterior evidence

The expression levels of the pathway's genes in region $r$, on the scale of the normalized expression in the Gaussian model and of the log expected count in the negative binomial model, follow $Z_r\sim N(\mu_r,\tau_r^2)$, centered on the regional mean and spread by the gene-to-gene variance, and the directional score $\pi_{rs}=P(Z_r>Z_s)=\Phi\{(\mu_r-\mu_s)/(\tau_r^2+\tau_s^2)^{1/2}\}$ is the probability that a gene of the pathway is expressed at a higher level in region $r$ than a gene of the pathway in region $s$. We report its posterior mean, the π score, which has the closed form

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
