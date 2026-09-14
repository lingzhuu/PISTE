## Breast cancer
#### Data Description 
[A single-cell and spatially resolved atlas of human breast cancers. (Wu et al.)](https://www.nature.com/articles/s41588-021-00911-1)

[data](https://zenodo.org/records/4739739)

Processed count matrices, brightfield HE-images (plain and annotated) and meta-data (containing clinical information and spot pathological details) for 6 primary breast cancers profiled using the Visium assay (10X Genomics).
The content of the files are:
`raw_count_matrices.tar.gz` - spaceranger processed raw count matrices.
`spatial.tar.gz` - spaceranger processed spatial files (images, scalefactors, aligned fiducials, position lists)
`filtered_count_matrices.tar.gz` - filtered count matrices.
`metadata.tar.gz` - metadata for tissues and spots of filtered count matrices, including clinical subtype and pathological annotation of each spot.
`images.pdf` - pdf detailing the H&E and annotation images.

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
