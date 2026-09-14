## Breast cancer
#### Data Description 
[A single-cell and spatially resolved atlas of human breast cancers. (Wu et al.)](https://www.nature.com/articles/s41588-021-00911-1)

[data](https://zenodo.org/records/4739739)

ST 部分包含 6 个原发乳腺癌样本：4 个 TNBC、2 个 ER+，使用 Visium，spot 直径约 55 μm。有独立的病理区域标注, 例如，可以把标注映射成肿瘤区、基质区、淋巴细胞富集区，再研究这些区域的通路差异，但是sample少。

#### Data organization
filtered_count_matrices.tar.gz
└─ filtered_count_matrices/CID4535_filtered_count_matrix/
   ├─ matrix.mtx.gz       表达计数矩阵：gene × spot
   ├─ features.tsv.gz     每一行对应的 gene symbol
   └─ barcodes.tsv.gz     每一列对应的 spot barcode

metadata.tar.gz
└─ metadata/CID4535_metadata.csv
                          spot 信息与病理注释

spatial.tar.gz
└─ spatial/CID4535_spatial/
   ├─ tissue_positions_list.csv   barcode 与空间坐标
   ├─ scalefactors_json.json      图像坐标缩放系数
   ├─ tissue_hires_image.png      H&E 图像
   └─ tissue_lowres_image.png     低分辨率 H&E 图像
