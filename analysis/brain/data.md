## Whole Mouse brain
#### Data Description 
[Single-cell spatial transcriptomic atlas of the whole mouse brain. (Han et al.)](https://www.cell.com/neuron/abstract/S0896-6273(25)00133-3)

Using single-nucleus RNA sequencing (snRNA-seq) and **Stereo-seq**, genome-wide coverage. 
> Coronal sections (10-μm thick) of the left hemisphere were prepared at 100-μm intervals, producing 123 sections (after quality control) with expression profiles for 29,655 genes, covering 95.5% of annotated protein-coding and non-coding genes. We defined the boundary of each cell based on the image of single-stranded DNA (ssDNA). A total of 4,229,623 cells were characterized in one mouse brain, with an average of 1,267 molecular identifiers (MIDs) and 668 detected genes per cell. Moreover, we manually segmented the brain areas for each section based on cytoarchitectural pattern and brain region delineation in ABA common coordinate framework (CCFv3)27.

The processed data can be downloaded in [Whole Mouse brain data](https://www.braindatacenter.cn/datacenter/web/#/dataSet/details?id=1888062851995131906). Raw data available upon request. 

The author has performed quality control steps, including:

1. 去除 7 张捕获异常切片、5 张低质量切片，最终保留 195 张
2. 用 SAW 处理原始序列，过滤低质量 UMI、进行基因组比对和 UMI 去重
3. 使用检测基因数 100、线粒体比例 15% 两个阈值进行筛选，去除 20,019 个细胞，保留 4,209,604 个高质量细胞
4. Bin100: 过滤低基因数、高线粒体比例的 bin，然后进行 SCTransform

#### Data organization
主要数据共分2类，分别为：
- 连续整脑系列是冠状的：以 100 µm 间隔、沿冠状方向切 10 µm 厚切面，质控后 195 张可用（mouse #1 的 123 张 + mouse #2 的 72 张）。
- 发育系列每个时间点只有 1 张矢状切面，一共 7 张跨 E12.5–P77。

命名规则：

- total_gene_P14_mouse_2D-20221120.txt.gz
发育系列的切片，出生后第 14 天的小鼠单张矢状切面。矢状切面在一个平面内同时切到了大部分主要脑区（嗅球、皮层、海马、丘脑、小脑、脑干都在里面）。未标准化。
- P14_bin100_cluster.h5ad
上述切片bin100后+聚类结果
- total_gene_T239_mouse_f001_2D_mouse1-20230119.txt.gz
单张嗅球切片，gene / x / y / umi_count / cell_label / gene_area / rx / ry
  - cell_label：0 = 胞外背景，1–3,587 = 分割出的细胞。只有 8.4 M UMI（26%）落在细胞内，其余 23.6 M 是背景/游离 RNA。
  - gene_area：人工脑区注释。0 = 未注释。
  - rx, ry：配准后坐标，用于把切面对齐到统一解剖框架。
- T239_T248_SCT_merge_bayes_anno.rds
名字写着叫SCT，其实是raw counts (不懂何意味..)。对应T239–T248 十张连续嗅球切面。

大量胞外游离RNA (Fig.1c)来源是 Stereo-seq 的固有特性：pepsin 透化时释放的 RNA 在被 DNB 阵列捕获前会侧向扩散，加上全芯片的环境 RNA 背景。This prevents the library-size exposure $L_i$ from being inflated by ambient RNA outside the segmented cells.
