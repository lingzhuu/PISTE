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

SpaPath quantifies how the activity of a gene set (pathway) differs between tissue regions in spatial transcriptomics data. It takes a normalized expression matrix (for example, SCTransform Pearson residuals), spot coordinates, region labels and a gene set, and fits one Gaussian linear mixed model jointly to all genes of the pathway and all regions.

Let $i=1,\ldots,N$ index spots, $g=1,\ldots,K$ the genes of the pathway, and $r(i)\in\{1,\ldots,R\}$ the region of spot $i$. For one pathway, SpaPath fits

$$
Y_{ig}=X_i^T\beta+\mu_{r(i)}+b_{g,r(i)}+w_g(s_i)+\epsilon_{ig},
\qquad
b_{gr}\mid\tau_r^2\sim N(0,\tau_r^2),
\qquad
\epsilon_{ig}\sim N(0,\sigma_\epsilon^2).
$$

Here $X_i$ are spot-level covariates such as batch, and their coefficients $\beta$ are shared by all genes. The regional mean $\mu_r$ gives the typical pathway level in region $r$, while $b_{gr}$ describes how far gene $g$ departs from that level, so that $\tau_r^2$ measures how widely the pathway's genes spread within region $r$. Beyond these regional terms, $w_g(s_i)$ is a smooth spatial field for gene $g$ evaluated at the location $s_i$ of spot $i$, and $\epsilon_{ig}$ is the remaining measurement noise. The spatial fields are independent across genes but share one covariance function,

$$
w_g(\cdot)\overset{iid}{\sim}\mathrm{GP}(0,\sigma_s^2C_{\kappa_0}),
\qquad g=1,\ldots,K,
$$

so that each gene has its own spatial pattern on a common scale. Because these fields span the whole tissue and do not restart at region boundaries, they absorb smooth spatial trends rather than regional differences. As a result, the model has only $R+2$ free hyperparameters, namely $\tau_1^2,\ldots,\tau_R^2$, $\sigma_s^2$ and $\sigma_\epsilon^2$, and this number does not depend on the number of genes.

To turn the fitted model into a statement about pathway activity, we summarize each region by the state of a randomly chosen pathway gene once covariate and spatial effects are removed, $Z_r\sim N(\mu_r,\tau_r^2)$. For two regions $r$ and $s$, the directional score

$$
\pi_{rs}=P(Z_r>Z_s\mid\mu_r,\tau_r^2,\mu_s,\tau_s^2)
=\Phi\!\left(\frac{\mu_r-\mu_s}{\sqrt{\tau_r^2+\tau_s^2}}\right)
$$

is then the probability that the pathway state in region $r$ exceeds that in region $s$. Since $\pi_{rs}$ depends on unknown parameters, SpaPath summarizes its posterior in two complementary ways. The π score $\bar\pi_{rs}=E(\pi_{rs}\mid Y)$ describes the expected direction and size of the difference, whereas the posterior evidence $P(\pi_{rs}>0.5\mid Y)$ gives the posterior probability that the pathway is more active in region $r$ than in region $s$. When one region is compared with the rest of the tissue, the pairwise scores are averaged over the other regions with weights proportional to their numbers of spots, $\pi_{r,\mathrm{rest}}=\sum_{s\neq r}a_{rs}\pi_{rs}$ with $a_{rs}\propto n_s$, and the π score and evidence of this average are reported in the same way. Regional mean differences and variance components, in contrast, are reported only as diagnostics.

Regional differences in pathway activity can arise either because regions contain different cell types or because the same cell types behave differently across regions. SpaPath therefore offers two analysis modes, and the user chooses between them. By default, the observed mode fits the model above and reports the total regional signal, including the part driven by cell composition. Alternatively, the composition-adjusted mode adds a term $\lambda_g\tilde o_{ig}$, in which $\tilde o_{ig}$ is the expression of gene $g$ expected in spot $i$ from its deconvolved cell-type composition (Methods), so that the regional effects describe pathway differences beyond those predicted by cell composition. Because the two modes share all other parameters, their π scores and evidence can be compared directly.

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
### Model detail

SpaPath models normalized expression with Gaussian errors, and normalization and quality control are therefore performed beforehand; in this study we used SCTransform Pearson residuals. For one pathway, the input consists of the $N\times K$ matrix $Y$ of its genes, the spot coordinates $s_i$, the region labels $r(i)$, optional spot-level covariates $X_i$, and a triangulation mesh of the tissue.

Within each region $r$, the gene effects $b_{gr}$ are conditionally independent across genes given $\tau_r^2$, and every region has its own variance $\tau_r^2$. Together with $\mu_r$, this defines the regional pathway state $Z_r\sim N(\mu_r,\tau_r^2)$, in which $\mu_r$ is the typical level of the pathway's genes in region $r$ and $\tau_r^2$ reflects how consistently they behave. Consequently, a pathway whose genes all shift in the same direction yields a large $|\mu_r-\mu_s|$ relative to $\tau_r^2+\tau_s^2$ and hence a π score near 0 or 1, whereas a pathway whose genes move in opposite directions yields a score near 0.5. In the implementation, $\mu_r$ is estimated relative to a baseline region and the regional means are recovered as linear combinations; this choice does not affect the π scores or evidence.

Spatial variation in gene expression is modeled with a replicated stochastic partial differential equation (SPDE) framework fitted with INLA. Each gene $g=1,\ldots,K$ has an independent Matérn field $w_g(\cdot)\overset{iid}{\sim}\mathrm{GP}(0,\sigma_s^2C_{\kappa_0})$ with smoothness $\nu=1$. All genes share one triangulation mesh, a fixed scale parameter $\kappa_0$ and a single spatial variance $\sigma_s^2$, so spatial structure is gene-specific in pattern but common in scale. Moreover, each field spans the entire tissue and is not re-initialized at region boundaries. The scale is fixed as $\kappa_0=\sqrt8/\rho_0$, where $\rho_0$ is a user-specified practical range, that is, the distance at which the spatial correlation falls to about 0.1.

In the composition-adjusted mode, the model additionally accounts for differences in cell composition between spots. Let $\omega_{ik}$ be the proportion of cell type $k$ in spot $i$, estimated by deconvolution (for example, cell2location), and let $\nu_{kg}$ be the mean expression of gene $g$ in cell type $k$ in a depth-normalized single-cell reference. The expression of gene $g$ expected in spot $i$ from its composition alone is then

$$
o_{ig}=\log\Big(\sum_k\omega_{ik}\nu_{kg}+c_{\mathrm{ref}}\Big),
$$

where the pseudo-count $c_{\mathrm{ref}}$ is chosen on the scale of $\nu_{kg}$, for example $c_{\mathrm{ref}}=1$ for counts per million. We center this quantity within each gene, $\tilde o_{ig}=o_{ig}-N^{-1}\sum_{i'}o_{i'g}$, and add a single term to the model:

$$
Y_{ig}=X_i^T\beta+\mu_{r(i)}+b_{g,r(i)}+\lambda_g\tilde o_{ig}+w_g(s_i)+\epsilon_{ig}.
$$

The coefficient $\lambda_g$ measures how closely gene $g$ follows the expression predicted from cell composition, and it is estimated separately for each gene because the reference and the spatial data differ in scale and platform. Centering removes each gene's overall level from $o_{ig}$, so that $\lambda_g$ reflects only how expression varies across spots with composition, while the gene's overall level remains in $b_{gr}$ exactly as in the observed mode. All other parameters therefore keep their meaning, and $\mu_r$ becomes the regional pathway level at the tissue-average composition. It follows that the π scores and evidence of the two modes can be compared directly: the observed mode gives the total regional signal, and the composition-adjusted mode gives the part not explained by cell composition. This adjustment nevertheless treats $\omega_{ik}$ and $\nu_{kg}$ as known and assumes that a cell type has the same expression profile in every region.

All fixed effects, that is, the covariate coefficients $\beta$ and the regional means (parameterized relative to the baseline region), have $N(0,10^2)$ priors. The regional log precisions have $\log(1/\tau_r^2)\sim N(0,3^2)$, a proper approximation to a log-flat prior on the scale. The residual precision has a penalized-complexity prior with $P(\sigma_\epsilon>2)=0.05$, and the spatial standard deviation has $\log\sigma_s\sim N(0,1)$. In the composition-adjusted mode, $\lambda_g\overset{iid}{\sim}N(0,10^2)$ with a fixed prior variance, so this mode adds $K$ latent coefficients but no hyperparameters.

The model is fitted with INLA, which integrates over the hyperparameters using the central composite design (CCD). Both modes have $R+2$ free hyperparameters, and because the number of CCD configurations, and with it the run time, grows quickly with this number, SpaPath warns when it exceeds a user-set limit (10 by default); in that case merging regions is advisable. For each CCD configuration $k$ with normalized weight $q_k$, INLA provides a conditional Gaussian posterior for the regional means, $(\mu_1,\ldots,\mu_R)\mid\theta_k,Y\sim N(m_k,\Sigma_k)$, together with the regional variances $\tau_{r,k}^2$. The joint posterior of the regional means is accordingly approximated by the mixture of these Gaussians.

Given configuration $k$, the difference $Z_r-Z_s$ is Gaussian, so the π score has a closed form,

$$
\bar\pi_{rs}=\sum_k q_k\,\Phi\!\left(\frac{m_{k,r}-m_{k,s}}{\sqrt{\tau_{r,k}^2+\tau_{s,k}^2+\Sigma_{k,rr}+\Sigma_{k,ss}-2\Sigma_{k,rs}}}\right).
$$

Furthermore, because $\pi_{rs}>0.5$ exactly when $\mu_r>\mu_s$, the pairwise evidence equals $P(\mu_r>\mu_s\mid Y)$, which is computed analytically from the mixture when there are two regions. With three or more regions, however, the region-versus-rest evidence $P(\pi_{r,\mathrm{rest}}>0.5\mid Y)$ has no closed form, and the pairwise and region-versus-rest evidence are therefore computed together by deterministic quasi-Monte Carlo integration over the mixture using antithetic, unrandomized Sobol points. The number of points is doubled until two consecutive refinements change the evidence by less than a tolerance (0.002 by default) and reproduce the analytic π score. Since no posterior sampling is involved, the results are fully reproducible.

Finally, SpaPath reports posterior summaries of $\mu_r$, $\tau_r^2$, $\sigma_s^2$, $\sigma_\epsilon^2$ and, in the composition-adjusted mode, $\lambda_g$, together with regional mean contrasts $\mu_r-\sum_{s\neq r}a_{rs}\mu_s$, their posterior probability of being positive and their local false sign rate. These quantities serve as diagnostics and are not used to rank pathways.
