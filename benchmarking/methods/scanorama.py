import scanpy as sc
import gc
import scanorama
import numpy as np

adata = sc.read('/lustre/scratch126/cellgen/team292/ha10/data/Heart_Atlas/adata_Heart_Reichart_HV_train.h5ad')

donors = list(adata.obs['donor_id'].unique())

alldata = {}
for batch in donors:
    alldata[batch] = adata[adata.obs['donor_id'] == batch,]

adatas = list(alldata.values())

scanorama.integrate_scanpy(adatas, dimred = 25, approx=False)

# Get all the integrated matrices.
scanorama_int = [ad.obsm['X_scanorama'] for ad in adatas]

# make into one matrix.
all_s = np.concatenate(scanorama_int)
print(all_s.shape)

# add to the AnnData object, create a new object first
adata_sc = adata.copy()
adata_sc.obsm["X_scanorama"] = all_s

np.savetxt("data3.csv", all_s, delimiter = ",")

adata_sc.write('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/adata_scanorama_notfinal_oct24_d25.h5ad', compression='gzip')
