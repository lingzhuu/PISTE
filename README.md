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
PathwayLGM uses a Gaussian linear mixed model (LMM) for normalized expression data and a negative binomial generalized linear mixed model (GLMM) for raw UMI counts. 两者都包含随机效应。Gaussian 版本用正态误差描述观测波动, SCT/normalized read counts 已经校正 library size，因此默认 additive offset 为 0；不再次加入 log library size；NB 版本通过负二项分布描述计数波动，通常不再另加一个 Gaussian 观测误差项 \(\epsilon_{ig}\).

Gaussian / Normalized reads：
设 $i$ 表示 spot，$g$ 表示基因，$r(i)$ 表示 spot 所属的区域。对于一条 pathway，代码拟合：

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

在相同的平均表达下，$\alpha_g$ 越大，计数越分散。这与 INLA 使用的 size 参数 $n_g$ 对应，关系为 $\alpha_g=1/n_g$.
NB dispersion 使用降维参数化：

$$
\alpha_g=\alpha_0\alpha_g^{rel},
\qquad
n_g=S_g\exp(\theta_{NB}),
\qquad
S_g=1/\alpha_g^{rel}.
$$

$\alpha_g^{rel}$ 由全 transcriptome 的 DESeq2 dispersion trend/ maximum a posteriori (MAP) shrinkage estimation 预估并进行 geometric-mean normalization。假设全转录组中用于归一化的有效基因，其 dispersion 几何均值为：

$$
G=\exp\left(\frac{1}{J}\sum_{g=1}^{J}\log\widehat\alpha_g\right), \quad\alpha_g^{rel}=\frac{\widehat\alpha_g}{G}
$$

这里的 dispersion trend 可以理解为：
1.先根据每个基因的数据，初步估计它自己的 dispersion。
2.利用大量基因，拟合平均表达水平与 dispersion 之间的总体趋势。
3.结合这个趋势，稳定各基因的估计。信息不足的基因通常更受趋势约束；信息充分的基因保留更多自己的数据特征。
所以，它利用其他基因的信息，减少单个基因估计过于不稳定的问题。收缩的参照与平均表达水平有关，并不是把所有基因都拉向同一个 dispersion。
整个数据集只做一次 DESeq2 预处理，保存每个基因的 \(\alpha_g^{rel}\)。分析不同通路时，取出对应基因的相对值重复使用。

进入通路模型后，只估计 \(\alpha_0\), 每条通路可以估计自己的 $\alpha_0$ 和 \(\theta_{NB}\)在该通路内共同使用。Therefore每个 pathway 的 INLA 模型只估计一个全局 NB calibration parameter。
INLA 使用\(n_g=\frac{1}{\alpha_g}\), 代入前面的分解：

$$
n_g
=\frac{1}{\alpha_0\alpha_g^{rel}}
=\underbrace{\frac{1}{\alpha_g^{rel}}}_{S_g}
\underbrace{\frac{1}{\alpha_0}}_{\exp(\theta_{NB})}
$$

$$
\boxed{\alpha_0=\exp(-\theta_{NB})}
$$

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

## Discussion
#### summary
SpaPath addresses a fundamental dependence problem in spatial pathway analysis.
#### 和现有方法区别
#### Biology insight 
#### Limitations
- computational cost (INLA对hyperparameter sensitive)
- sensitivity to spatial covariance assumptions
- pathway gene dependence (i.i.d 但实际是稀疏正负相关)
- pathway overlap
- uncertainty in cell type labels
- future work

## Methods

### Data prepration
QC, SCTransform








