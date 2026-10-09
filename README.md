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

SpaPath can be run in two modes. The observed mode fits a Gaussian model to normalized expression and quantifies the total difference in pathway activity between regions, including differences that arise because the regions contain different cell types. The composition-adjusted mode fits a negative binomial model to raw counts. It includes the expression of each gene predicted from the cell-type proportions of each spot as a covariate with a gene-specific coefficient, so that the regional difference it reports is the part that cannot be explained by cell composition. We recommend running both modes, because together they distinguish pathways whose regional differences are driven by cell composition, which are supported only in the observed mode, from pathways that change within cell types, which remain supported after adjustment.

### SpaPath detects weak pathway differences

- GSVA etc single cell pathway method in ST data
- STpathway [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11503753/)
  > employing the Brownian distance correlation (dCor) test in a novel set-up to evaluate the dependency of two random vectors with different dimensionality.
- GSDensity [paper](https://www.nature.com/articles/s41467-023-44206-x)
  > Here, we present GSDensity, a graph-modeling approach that allows users to obtain pathway-centric interpretation and dissection of single-cell and spatial transcriptomics (ST) data without performing clustering. Using pathway gene sets, we show that GSDensity can accurately detect biologically distinct cells and reveal novel cell-pathway associations ignored by existing methods. 
  
We generated seven simulated datasets with different random seeds to benchmark SpaPath against three representative pathway-analysis workflows: GSVA, the gene-set scoring approach popularized in single-cell RNA sequencing, tested here with a spatial block bootstrap region test; STpathway, a distance-correlation (dCor) test for region-associated pathways; and GSDensity, a region enrichment test on dichotomized spot activity. Each dataset was built on the real spot coordinates of DLPFC section 151673 and comprised 40 non-overlapping 100-gene sets spanning five spatial activity patterns and two effect-size distributions, with the number K of genes responding to the activity map prespecified for each gene set (Methods). Every method compared the target region of each gene set's activity map with the rest of the tissue, giving 280 gene sets in total and 35 in each of the eight groups defined by K and effect distribution.

The existing methods detected strong pathway differences but could not separate them from weaker ones. GSVA's target-minus-rest score difference increased with K, and its spatial block bootstrap test gave median q-values below 0.05 in all groups with K of 11 or more. However, the test saturated. Of the 210 gene sets with K of 11 or more, 183 reached the smallest p-value the bootstrap can attain and received the same q-value, so gene sets with 11 and with 100 responsive genes could not be ranked against each other. STpathway's region dCor decreased with K in the narrow, strong-effect setting, from a median of 0.98 at K = 0 to 10 to 0.94 and 0.95 at K of 31 or more, but stayed at 0.98 in all groups of the broad, weak-effect setting. Across all 280 gene sets it spanned only 0.87 to 0.99, offering limited discrimination of pathway activity. GSDensity's coherence filter retained 137 of the 280 gene sets and none with K of 10 or less, so weak pathways were not evaluated. Among the retained sets, its region test gave q < 0.05 for only 62, and its significance did not strengthen with K. In the broad, weak-effect setting the median q-value was 0.003 at K = 31 to 60 but 0.11 at K = 61 to 100. The test compares the proportions of spots labelled as highly active, and GSDensity sets this label by assuming that pathway activity has two modes across spots and cutting at the low point between them (Methods). When the activity has a single mode with a long upper tail, the cut falls in the tail and almost no spot is labelled, whatever the regional difference. For example, for a gradient gene set with 94 of its 100 genes responding in the first simulated dataset, the mean GSDensity activity was 1.5 times higher in the target than in the rest, yet the cut lay above the 99.9th percentile of the activity, so only 2 of 728 target spots and none of the 2,911 rest spots were labelled highly active and the test was not significant (q = 0.14). In the narrow, strong-effect setting with K of 31 or more, fewer than 1% of spots were labelled highly active in both the target and the rest for 37 of the 67 retained sets.

SpaPath's π score and posterior evidence instead increased steadily with K in both settings. The median π score rose from 0.48 at K = 0 to 10 to 0.85 at K = 61 to 100 in the broad, weak-effect setting, and from 0.56 to 0.79 in the narrow, strong-effect setting. The median evidence was 0.40 and 0.75 in the two K = 0 to 10 groups, exceeded 0.98 in the K = 11 to 30 groups and was close to 1 in the higher-K groups. Unlike GSVA's q-value, the π score continued to rank pathways that already had strong support. Among the 194 gene sets that reached GSVA's smallest p-value, the π score ranged from 0.34 to 0.95 and its group median rose with K. Of the 75 retained gene sets that GSDensity's region test did not detect, 74 received posterior evidence above 0.99 from SpaPath. The increase of the median π score with K held for each of the five spatial activity patterns.

[Table 1. Median region-level statistics of the four methods across the 280 simulated gene sets.](benchmark/analysis/simulation_data/result/simulation_summary.xlsx)

[Supplementary Table 1. Per-pathway simulation results for each of the seven simulated datasets.](benchmark/analysis/simulation_data/result/simulation_per_pathway_result.xlsx)
### Resolving layer-specific pathways in human DLPFC

To evaluate the performance of these methods on real data, we analyzed a human dorsolateral prefrontal cortex (DLPFC) Visium dataset because its six cortical layers and white matter were manually annotated and the molecular functions of each layer are well documented. [Maynard et al., 2021] We applied SpaPath, GSVA, STpathway, and GSDensity to section 151673, which contains 3,639 spots. We tested a prespecified panel of 20 MSigDB gene sets (Methods; Supplementary Table 2).

[Supplementary Table 2. The 20 prespecified MSigDB gene sets for the DLPFC analysis.](benchmark/analysis/DLPFC/data/geneset_20.xlsx)

We compared the layer-versus-rest results of the four methods. Using the WM-versus-rest comparison as an example, most methods recovered the strong expected signals: lipid metabolism and myelin-related pathways were inferred as up-regulated in WM, consistent with myelin being 70–80% lipid by dry weight (Baumann and Pham-Dinh, 2001; Simons and Nave, 2015), whereas synaptic transmission pathways were inferred as down-regulated, consistent with the signaling-dominated energy budget of gray matter (Harris, Jolivet and Attwell, 2012). The methods diverged on weaker signals. STpathway's region dCor ranked lipid and myelin pathways as the most regionally divergent but spanned only 0.88 to 0.99 across all 20 pathways, offering limited discrimination between strong and weak signals. GSDensity reported only 12 of the 20 gene sets after its coherence filter. Moreover, the discarded sets were not random. They included myelin maintenance, sphingolipid metabolism and very-long-chain fatty acid metabolism, which GSVA and SpaPath both identified as more active in WM. Among the reported sets, GSDensity compares the proportions of spots labelled as highly active, so it missed differences that did not change these proportions, as for myelin sheath, where almost no spot was labelled in either region (p = 0.31), and glutamatergic synaptic transmission, where similar proportions were labelled in both (31% versus 29%, p = 0.27). For GSVA, we tested WM-versus-rest score differences with a two-sided spatial block bootstrap (Methods). The test recovered the expected directions, with the six myelin and lipid sets higher and the synaptic, potassium channel and oxidative phosphorylation sets lower in WM, but it could not distinguish strong from weak signals, because 15 of the 20 sets reached the smallest p-value the bootstrap can attain and received the same q = 0.00067.

In contrast, SpaPath recovered all strong signals in the correct direction while remaining sensitive to weaker, biologically expected ones (Fig. X). Here the target region is WM and the rest consists of the gray-matter cortical layers L1–L6, so the π score measures how often a WM pathway state exceeds the spot-weighted gray-matter state, and the posterior evidence is the probability that it does (π > 0.5): evidence near 1 indicates higher activity in WM and evidence near 0 higher activity in gray matter. When the 15 brain-related gene sets were positioned by their WM-versus-rest π score and posterior evidence, they resolved into two distinct, biologically coherent groups (Fig. X). All seven myelin- and lipid-related sets were inferred to be more active in WM: three with evidence > 0.99 (myelin sheath, sphingolipid metabolism, and monocarboxylic acid transport), three with evidence > 0.95 (axon ensheathment, myelin maintenance, and very-long-chain fatty acid metabolism), and one, cholesterol biosynthesis, with evidence ≥ 0.90 (0.912). Collectively, these sets recapitulate the oligodendrocyte lipid-synthesis program that underlies myelination. Conversely, all synaptic transmission (n = 6), voltage-gated potassium channel, and oxidative phosphorylation sets were inferred to be more active in gray matter (cortical layers L1–L6) with evidence ≤ 0.05, reflecting three defining properties of this compartment: dense synaptic connectivity, neuronal excitability, and high energetic demand. Notably, oxidative phosphorylation received strong support for lower WM activity (evidence = 2.5 × 10⁻⁵) despite a modest effect size, consistent with a recent human brain map of mitochondrial respiratory capacity showing that gray matter contains over 50% more mitochondria than white matter (Mosharov et al., 2025), whereas other benchmark methods could not distinguish this signal from noise.

We also observed that several pathways representing basic cellular processes, including mRNA splicing and nucleotide excision repair, tended toward higher inferred activity in WM (evidence > 0.9). We hypothesized that these differences arise from regional variation in cellular composition. WM is dominated by oligodendrocytes, which produce and maintain large amounts of myelin membrane and protein (Pfeiffer et al., 1993; Nave and Werner, 2014), rely on extensive QKI-regulated alternative splicing of myelin transcripts (Wu et al., 2002), and are prone to oxidative DNA damage owing to their high iron content and metabolic rate (Bradl and Lassmann, 2010). Elevated splicing and DNA repair activity in WM may therefore reflect oligodendrocyte enrichment rather than region-intrinsic regulation. To test this, we estimated cell-type composition by RCTD deconvolution (Cable et al., 2022) using the LIBD spatialDLPFC snRNA-seq reference (Huuki-Myers et al., 2024), which confirmed a higher oligodendrocyte proportion in WM than in cortical layers (Fig. Xa). To quantify the expected contribution of composition, we derived cell-type-specific pathway activities from the same reference using GSVA (Fig. Xb) and weighted these scores by the estimated cell-type proportions. The resulting composition-predicted differences correlated with π scores across all 20 pathways (Pearson r = 0.856, R² = 0.732; Fig. Xc), yet composition did not account for all signals. Most notably, oxidative phosphorylation showed opposite directions: the snRNA-seq reference predicted higher activity in WM, whereas SpaPath inferred lower activity from the spatial data (π < 0.5). Because nuclear transcriptomes underrepresent the cytoplasmic and mitochondrial transcripts enriched in the gray-matter neuropil, the composition-based prediction misses this signal, whereas whole-tissue spatial profiling retains it, consistent with the higher mitochondrial density of gray matter (Mosharov et al., 2025).

[Figure 1. SpaPath recovers regional pathway activity differences in DLPFC](benchmark/analysis/DLPFC/figure/celltype.png)

These results indicate that a substantial fraction of the inferred WM–gray matter differences in basal cellular pathways can be explained by regional differences in cell composition. These considerations motivated us to distinguish composition-driven pathway differences from those that persist after accounting for cell composition. To this end, we extended SpaPath with a composition covariate. For each spot, we used its deconvolved cell-type proportions to compute the expected expression of each gene and included this quantity as a covariate with a gene-specific coefficient (Methods). The regional effects therefore capture pathway differences beyond those predicted by cell composition.

### Disentangling cell composition from pathway activity

Cell composition acts on the expected counts of a spot, which are a mixture of cell-type expression profiles, and this relationship is lost after normalization such as Pearson residuals. The composition-adjusted mode therefore models the raw counts with a negative binomial likelihood and a sequencing-depth offset, and adds a composition term $\lambda_g\log m_{ig}$, where $m_{ig}=\sum_k\omega_{ik}\nu_{kg}$ is the expression of gene $g$ in spot $i$ expected from its cell-type proportions $\omega_{ik}$ and the reference expression $\nu_{kg}$ (Methods). The coefficient $\lambda_g\sim N(1,0.25^2)$ is estimated for each gene rather than fixed at 1, and each gene has its own intercept. The regional levels therefore describe the pathway differences that remain after accounting for cell composition.

We applied this mode to the WM-versus-rest comparison for the same 20 gene sets (Fig. 2). The myelin and lipid sets remained more active in WM, with evidence between 0.86 and 0.997, indicating that these differences are not explained by the higher oligodendrocyte proportion of WM alone. The exception was cholesterol biosynthesis, whose π score decreased to 0.42 after adjustment. Oxidative phosphorylation remained lower in WM (evidence 2.2 × 10⁻⁹). Synaptic vesicle exocytosis, the glutamate release cycle and NMDA receptor activation remained lower in WM (evidence ≤ 0.023), whereas glutamatergic transmission, GABAergic transmission and voltage-gated potassium channel activity moved close to no direction. Among the five basic cellular process sets, all except mRNA splicing no longer showed higher activity in WM, supporting the hypothesis that their observed-mode differences arise from cell composition.

[Figure 2. WM-versus-rest π scores and posterior evidence of the 20 DLPFC gene sets in the observed and composition-adjusted modes](benchmark/analysis/DLPFC/figure/adjust_vs_observe.png)

### Mapping pathway activity in 3D mouse brain

Brain regions are three-dimensional structures, whereas a single coronal section captures only a cross-section of each region and only the regions present in that section. A three-dimensional analysis describes each region over its full extent and compares pathway activity among all regions within the brain as a whole, providing a more comprehensive view of pathway activity across the brain. The SPDE representation of the Matérn field used by SpaPath (REF) is defined in any spatial dimension, so SpaPath extends to three dimensions by replacing the triangular mesh of a single section with a tetrahedral mesh of the whole tissue, while the rest of the model remains unchanged.

We applied SpaPath to the Stereo-seq whole-brain atlas of Han et al. (REF), in which consecutive coronal sections of an adult mouse brain were registered to the Allen Mouse Brain Common Coordinate Framework (CCFv3) (REF), with every spot assigned to an anatomical region. We analyzed the left hemisphere across 122 sections spanning 0.25–12.65 mm along the anterior–posterior axis, with adjacent sections a median of 100 µm apart (Methods). By analyzing all sections jointly, we reconstructed a three-dimensional atlas of regional pathway activity across the brain.

[TODO: results]

### SpaPath reveals [TODO] in breast cancer
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

In the observed mode, SpaPath fits a Gaussian model to the normalized expression $Y_{ig}$. For spot $i=1,\ldots,N$, gene $g=1,\ldots,G$ and region $r(i)\in\{1,\ldots,R\}$,

$$
Y_{ig}=X_i^T\beta+\mu_{r(i)}+b_{g,r(i)}+w_g(s_i)+\epsilon_{ig},
\qquad
b_{gr}\mid\tau_r^2\sim N(0,\tau_r^2),
\qquad
\epsilon_{ig}\sim N(0,\sigma_\epsilon^2),
$$

where $\beta$ are covariate effects shared by all genes and $\mu_r$ is the mean pathway level in region $r$. The gene effects $b_{gr}$ are conditionally independent across genes, and each region has its own variance $\tau_r^2$, which measures how consistently the pathway's genes behave in that region. Each gene has an independent Matérn field $w_g(\cdot)\overset{iid}{\sim}\mathrm{GP}(0,\sigma_s^2C_{\kappa_0})$, represented through a replicated stochastic partial differential equation (SPDE) on a shared mesh (REF). All fields share the scale $\kappa_0=\sqrt8/\rho_0$, where $\rho_0$ is a user-specified practical range at which the spatial correlation drops to about 0.1; fixing it rather than estimating it keeps the number of hyperparameters small. The fields also share a single variance $\sigma_s^2$, so spatial structure is gene-specific in pattern but common in scale. Furthermore, the fields span the whole tissue and do not restart at region boundaries. The model thus has $R+2$ hyperparameters regardless of the number of genes.

#### Negative binomial model for the composition-adjusted mode

In the composition-adjusted mode, SpaPath models the raw counts $C_{ig}$. Let $\omega_{ik}$ be the proportion of cell type $k$ in spot $i$ and $\nu_{kg}$ the depth-normalized mean expression of gene $g$ in cell type $k$ in the reference. The expression of gene $g$ in spot $i$ expected from composition alone is

$$
m_{ig}=\sum_k\omega_{ik}\nu_{kg},
$$

and SpaPath fits

$$
C_{ig}\sim\mathrm{NB}(\eta_{ig},\phi_g),
\qquad
\log\eta_{ig}=\log d_i+\lambda_g\log m_{ig}+X_i^T\beta+a_g+\mu_{r(i)}+b_{g,r(i)}+w_g(s_i),
$$

where $\eta_{ig}$ is the expected count, $d_i$ is the total count of spot $i$, $\lambda_g$ is a composition coefficient and $a_g$ is a gene intercept, with all other terms as in the Gaussian model. Because composition determines expression on the count scale, it enters through the log link, and $\log d_i$ is an offset. The coefficient $\lambda_g$ is estimated for each gene with a prior centered at 1, so that the model starts from the exact mixture but allows the relation between the reference and the spatial data to differ across genes. A gene-specific difference in scale between the reference and the spatial data multiplies the expected count by a constant and is absorbed by $a_g$, which also removes each gene's baseline level, so that $b_{gr}$ describes deviations from the regional level rather than overall abundance. Genes with $m_{ig}=0$ in any spot are removed before fitting. The size parameter of gene $g$ is $\phi_g=\phi S_g$, where $S_g$ is the inverse of its DESeq2 dispersion estimated on the whole transcriptome (REF), scaled to a geometric mean of 1, and $\phi$ is shared by all genes. The regional means $\mu_r$ thus describe the regional difference that composition does not explain. This adjustment treats $\omega_{ik}$ and $\nu_{kg}$ as known and assumes that a cell type has the same expression profile in every region. Like the Gaussian model, it has $R+2$ hyperparameters.

#### Priors and posterior inference

Fixed effects, including the regional means, have $N(0,10^2)$ priors. Each gene intercept $a_g$ is the sum of a common intercept and a gene-specific deviation with a $N(0,10^2)$ prior, and the prior of the common intercept is set so that the mean of the $a_g$, the baseline pathway level, also has a $N(0,10^2)$ prior. The composition coefficients have $\lambda_g\sim N(1,0.25^2)$. The regional log precisions follow $\log(1/\tau_r^2)\sim N(0,3^2)$ and $\log\sigma_s\sim N(0,1)$ in both models, the residual precision of the Gaussian model has a penalized-complexity prior with $P(\sigma_\epsilon>2)=0.05$, and the shared size parameter of the negative binomial model has $\log\phi\sim N(\log 10,1.5^2)$. INLA integrates over the hyperparameters with the central composite design (CCD) (REF); because the number of design points grows quickly with the number of hyperparameters, SpaPath warns when it exceeds a user-set limit (10 by default). For each design point $k$ with normalized weight $q_k$, INLA gives a Gaussian posterior $N(m_k,\Sigma_k)$ for the regional means together with the variances $\tau_{r,k}^2$, and the joint posterior of the regional means is approximated by the mixture of these Gaussians.

#### Directional pathway score and posterior evidence

The expression levels of the pathway's genes in region $r$, on the scale of the normalized expression in the Gaussian model and of the log expected count in the negative binomial model, follow $Z_r\sim N(\mu_r,\tau_r^2)$, centered on the regional mean and spread by the gene-to-gene variance, and the directional score $\pi_{rs}=P(Z_r>Z_s)=\Phi\{(\mu_r-\mu_s)/(\tau_r^2+\tau_s^2)^{1/2}\}$ is the probability that a gene of the pathway is expressed at a higher level in region $r$ than a gene of the pathway in region $s$. We report its posterior mean, the π score, which has the closed form

$$
\bar\pi_{rs}=\sum_k q_k\,\Phi\!\left(\frac{m_{k,r}-m_{k,s}}{\sqrt{\tau_{r,k}^2+\tau_{s,k}^2+\Sigma_{k,rr}+\Sigma_{k,ss}-2\Sigma_{k,rs}}}\right),
$$

together with the posterior evidence $P(\pi_{rs}>0.5\mid Y)$, which equals $P(\mu_r>\mu_s\mid Y)$. To compare a region with the rest of the tissue, we average the pairwise scores with spot-count weights, $\pi_{r,\mathrm{rest}}=\sum_{s\neq r}a_{rs}\pi_{rs}$ with $a_{rs}\propto n_s$. With two regions the evidence is computed analytically; with more regions it is computed by deterministic quasi-Monte Carlo integration over the mixture, doubling the number of Sobol points until two consecutive refinements change the evidence by less than 0.002, so no posterior sampling is involved. Regional mean contrasts and variance components are reported as diagnostics only.

### Simulation framework

#### Tissue layout and genes

We simulated raw spatial transcriptomic counts at the 3,639 tissue-covered spots of the dorsolateral prefrontal cortex (DLPFC) section 151673 from the spatialLIBD resource (REF), using the observed spot coordinates rescaled to $[0,1]$ with a common scale for both axes, so that all simulated spatial patterns are defined on a realistic tissue layout. Each simulated dataset contains 10,000 genes: 40 non-overlapping test pathways of 100 genes each, 1,000 background spatially variable genes (SVGs) and 5,000 background genes with no spatial effect. Seven independent datasets were generated with random seeds 2020 to 2026.

#### Pathway design

The 40 test pathways cover the full factorial of four response-size groups, five spatial-pattern classes and two effect-size distributions. For pathway $p$, the number of responsive genes $K_p$ among its 100 members was drawn uniformly from one of four ranges, 0–10, 11–30, 31–60 or 61–100, and the remaining $100-K_p$ members were left unresponsive. Pathways in the 0–10 group serve as null and weak-signal controls.

#### Spatial activity patterns

Each pathway was assigned a spatial activity map $A_i^{(p)}$, scaled to $[0,1]$, from one of five classes chosen to represent common tissue architectures. Discrete domains represent compartmentalized structures such as cortical layers, with piecewise-constant activity across three bands separated by randomly shifted, wavy boundaries. Smooth gradients represent continuous processes such as metabolic zonation and were drawn from a Gaussian random field with a Matérn covariance of smoothness 3/2 and range 0.25. Radial or niche-centered patterns represent niches centered on vessels or ducts, with activity decaying exponentially with distance from a randomly chosen spot at a scale drawn uniformly from 0.12 to 0.24. Localized hotspots represent focal structures such as tertiary lymphoid structures and consist of two Gaussian hotspots of widths 0.06 and 0.09 centered on randomly chosen spots. Interface-specific patterns represent tissue boundaries such as the tumor–stroma interface or an invasive front, with

$$
A_i^{(p)} = \mathbb{1}(h_i \le 0.14)\exp(-h_i/0.06),
$$

where $h_i$ is the distance of spot $i$ to a sinusoidal boundary. The target region of each pathway is the most active band for domain maps, the 20% most active spots for gradient and radial maps, the spots with $A_i^{(p)}\ge 0.5$ for hotspot maps and the spots with $A_i^{(p)}>0$ for interface maps. Responsive genes inherit their pathway's activity map, $H_{ig}=A_i^{(p)}$, whereas unresponsive pathway genes and background genes without a spatial effect have $H_{ig}=0$. Each background SVG received its own radial activity map, generated independently of the test-pathway maps, to mimic the spatially structured background variation of real tissue.

#### Effect sizes and response direction

For each responsive gene $g$, the absolute target-versus-rest log2 fold change of the depth-corrected expected expression, $E_g>0$, was drawn from one of two log-normal distributions,

$$
\log E_g \sim \begin{cases}
N(\log 0.4,\ 0.3^2), & \text{broad and weak},\\
N(\log 0.9,\ 0.4^2), & \text{narrow and strong}.
\end{cases}
$$

Response direction was assigned before count generation, with $d_g=+1$ for target-upregulated and $d_g=-1$ for target-downregulated genes. To mimic the heterogeneous composition of real pathways, $\lfloor 0.4\sqrt{K_p}\rfloor$ randomly chosen responsive genes of each pathway were downregulated and the others upregulated, so that pathways with fewer responsive genes contain a larger share of genes responding in the opposite direction. Background SVGs were assigned a random direction and an effect from the broad and weak distribution.

#### Count generation

Raw read counts were then generated from a negative binomial model,

$$
Y_{ig} \sim \mathrm{NB}(\mu_{ig}, \theta_g), \qquad
\mathrm{Var}(Y_{ig}) = \mu_{ig} + \frac{\mu_{ig}^2}{\theta_g},
$$

$$
\log \mu_{ig} = \log c_i + \log m_g + d_g b_g H_{ig},
$$

where $b_g\ge 0$ was calibrated numerically so that the target-versus-rest log2 fold change of the depth-corrected mean $m_g\exp(d_g b_g H_{ig})$ equals $d_g E_g$, and $b_g=0$ for genes without a spatial effect. The baseline mean $m_g$ and dispersion $\theta_g$ were estimated jointly from the observed counts of section 151673 by an intercept-only DESeq2 fit (REF), and each simulated gene received one intact empirical pair, sampled without replacement from the genes detected in at least 100 spots, so that the simulated data preserve the mean–dispersion relationship of Visium data. The spot depth factors $c_i$ are the DESeq2 size factors of the same section, randomly permuted across spots and rescaled to a geometric mean of 1. The simulated counts were then normalized with SCTransform (REF) for all downstream analyses, because STpathway and GSDensity both take SCTransform output as input.

#### Evaluation on simulated data

Each simulated dataset was analyzed with the four methods as described in Compared methods and benchmarking, with the target region of each pathway defined by its activity map. All methods used random seed 42. Each statistic was summarized by its median over the 35 gene sets of each of the eight groups defined by effect-size distribution and response-size group, without pooling tests across gene sets.

### Compared methods and benchmarking

#### Benchmarking design

We benchmarked SpaPath against three representative pathway-analysis methods, GSVA (REF), STpathway (REF) and GSDensity (REF). [TODO: software versions and code sources.] Because the methods return statistics that are not directly comparable, each method's output was mapped onto a common question, whether a pathway differs between a target region and the remaining tissue.

The benchmarking protocol was identical for simulated and real data, and only the dataset and the definition of the target region differed. On simulated data, the target region of each pathway was the one defined by its activity map, and each method's region-level statistic was compared with the prespecified pathway activity across the eight strata defined by response size and effect-size distribution, and separately across the five spatial-pattern classes.

#### DLPFC data

On the DLPFC data, all methods were applied to the 3,611 morphologically annotated spots of section 151673, with the 513 white matter (WM) spots as the target region and the remaining 3,098 spots as the rest. All methods used the Space Ranger filtered matrix and kept every spot with a layer annotation. No additional spot-level filtering was applied, because it can preferentially remove spots from biologically distinct regions and the genuine signal they carry (REF: Totty et al., 2025). WM serves as a partial biological reference because myelin-related and oligodendrocyte-associated pathways are expected a priori to be WM-enriched (REF). The DLPFC panel consisted of 20 MSigDB gene sets spanning synaptic transmission (n = 6), voltage-gated potassium channel activity (n = 1), oxidative phosphorylation (n = 1), myelination (n = 3), lipid metabolism (n = 4), and basal RNA/protein maintenance and DNA repair (n = 5). Sets were selected for size (16–210 detected genes in this section) and minimal redundancy (mean pairwise Jaccard index 0.008, maximum 0.19; Supplementary Table 2).

#### Normalization

All methods started from the same SCTransform fit of the raw counts (REF). GSDensity received the SCTransform-corrected counts, which its functions log-normalize before the embedding. GSVA, STpathway and SpaPath, run in the observed mode, received the Pearson residuals.

#### GSVA

GSVA scored each pathway in each spot from the SCTransform Pearson residuals, estimating the cumulative distribution of each gene across spots directly from the empirical distribution without kernel smoothing. Because the scores of neighbouring spots are spatially autocorrelated, we tested target-versus-rest differences in mean score with a two-sided spatial block bootstrap. Spot coordinates were scaled to the unit square and divided into 85 square blocks of side 0.10, and the blocks were resampled with replacement 1,999 times. P-values were adjusted by the Benjamini–Hochberg procedure across the gene sets tested in each dataset.

#### STpathway

STpathway adapts the Brownian distance covariance test (Székely and Rizzo, 2009) to pathway analysis. In its predefined-region mode, used here, the Pearson residuals of a set's genes are split into a target-region matrix and a rest matrix, gene-by-gene Euclidean distance matrices are computed within each region, and the dependence between the two is summarized by the distance correlation (dCor) and assessed with 199 permutations. A dCor close to zero indicates that the genes are organized differently in the two regions, so a lower dCor indicates a larger regional difference. The statistic carries no direction.

#### GSDensity

The GSDensity pathway activity of each spot was computed by random walk with restart on a 300-nearest-neighbour graph, the coherence of each set was tested against 100 size-matched random gene sets, and sets with Benjamini–Hochberg q < 0.05 within each dataset were retained. For these sets, spots were labelled as highly active or not with the GSDensity labelling function, which estimates the density of the pathway activity across spots with a Gaussian kernel at the critical bandwidth that yields two modes and labels spots above the antimode as highly active (locmodes in the multimode R package, REF). The activity of each set was divided by its mean before labelling, which leaves the labels unchanged and avoids numerical tolerance problems at the original scale. The proportions of highly active spots in the target region and the remaining spots were compared by a chi-squared test with continuity correction. P-values were adjusted by the Benjamini–Hochberg procedure across the retained gene sets in each dataset.

#### SpaPath

For SpaPath, we report the π score and posterior evidence of the target region versus the rest of the tissue, as defined in Directional pathway score and posterior evidence.

### Mouse whole-brain data

We used mouse 1 of the adult mouse whole-brain Stereo-seq atlas (REF), whose 123 coronal sections were registered to CCFv3, with each spot assigned to a CCFv3 structure and hemisphere. We kept the spots of the left hemisphere, because the right hemisphere was truncated unevenly across sections, whereas the left-hemisphere tissue area of every section was 0.87 to 1.14 times the CCFv3 hemisphere area at the same position. Spots without a region label, which lie outside the registered brain, were removed. One section with an inconsistent registered position was excluded. The remaining 122 sections are 40 to 270 µm apart (median 100 µm). Spots were aggregated into 100-µm bins within each section. [TODO: normalization, region comparisons and gene sets.]

The spatial field of each gene was defined on a tetrahedral mesh of the left hemisphere, built from the CCFv3 left-hemisphere mask together with all observed bins so that every bin lies inside the mesh. [TODO: mesh resolution, range of the spatial field and handling of section effects.]
