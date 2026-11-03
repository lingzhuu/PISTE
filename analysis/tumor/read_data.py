# %% 1. Import libraries and specify the data directory and sample
import tarfile
import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from PIL import Image

DATA_DIR = "data"
SAMPLE = "CID4535"

# %% 2. Read the author-provided filtered matrix, gene names, and spot barcodes
# The .gz files inside this archive are plain text and can be read directly.
# Some filtered values differ from the corresponding raw entries; additional processing has not been verified.
with tarfile.open(f"{DATA_DIR}/filtered_count_matrices.tar.gz") as archive:
    folder = f"filtered_count_matrices/{SAMPLE}_filtered_count_matrix"
    counts = mmread(archive.extractfile(folder + "/matrix.mtx.gz"), spmatrix=True).T.tocsr()
    genes = pd.read_csv(
        archive.extractfile(folder + "/features.tsv.gz"),
        sep="\t", header=None, names=["gene_symbol"],
    )
    barcodes = pd.read_csv(
        archive.extractfile(folder + "/barcodes.tsv.gz"), header=None,
    )[0].tolist()

# Rows of counts represent spots; columns represent genes.
var = genes.set_index("gene_symbol", drop=False)

# %% 3. Read pathology annotations and align them to the matrix barcode order
with tarfile.open(f"{DATA_DIR}/metadata.tar.gz") as archive:
    metadata = pd.read_csv(
        archive.extractfile(f"metadata/{SAMPLE}_metadata.csv"),
        index_col=0, keep_default_na=False,
    )

obs = metadata.loc[barcodes].copy()
obs["Classification"] = pd.Categorical(
    obs["Classification"].replace({"": np.nan, "NA": np.nan})
)

# %% 4. Read spatial coordinates, image scale factors, and H&E images
with tarfile.open(f"{DATA_DIR}/spatial.tar.gz") as archive:
    folder = f"spatial/{SAMPLE}_spatial"
    positions = pd.read_csv(
        archive.extractfile(folder + "/tissue_positions_list.csv"),
        header=None, index_col="barcode",
        names=["barcode", "in_tissue", "array_row", "array_col",
               "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )
    scalefactors = json.load(archive.extractfile(folder + "/scalefactors_json.json"))
    images = {}
    for resolution in ["hires", "lowres"]:
        with Image.open(archive.extractfile(folder + f"/tissue_{resolution}_image.png")) as im:
            images[resolution] = np.array(im.convert("RGB"))

positions = positions.loc[barcodes]
positions["in_tissue"] = positions["in_tissue"].astype(bool)
obs = obs.join(positions)

# %% 5. Create AnnData; coordinate order is x = pixel column, y = pixel row
adata = ad.AnnData(X=counts, obs=obs, var=var)
adata.obsm["spatial"] = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
adata.uns["spatial"] = {SAMPLE: {"images": images, "scalefactors": scalefactors}}

print("AnnData Construction:")
print(adata)
print("Spot counts by pathology classification:")
print(adata.obs["Classification"].value_counts(dropna=False))

# Save data to .h5ad
# adata.write_h5ad(f"{DATA_DIR}/{SAMPLE}_filtered.h5ad")