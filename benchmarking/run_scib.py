import numpy as np
import scanpy as sc
import gc
import os
import sys

from scib_metrics.benchmark import Benchmarker


adata = sc.read('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/objs_heart/_adata_allresults6.h5ad')

labelkey=sys.argv[1]

bm = Benchmarker(
    adata,
    batch_key="donor_id",
    label_key=labelkey,
    embedding_obsm_keys=['X_scanorama', 'X_scmerge2_pca', 'X_harmony', 'X_combat_pca', 'X_pca', 'X_scVI', 'X_fastmnn'],
    n_jobs=8,
)
bm.benchmark()


os.mkdir(f"/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_{labelkey}_oct24")
os.mkdir(f"/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_{labelkey}_oct24/min_max_scale_f_")
bm.plot_results_table(save_dir=f"/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_{labelkey}_oct24")
df = bm.get_results().transpose()
df.to_csv(f"/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_{labelkey}_oct24/scib_results.csv")

bm.plot_results_table(min_max_scale=False, save_dir=f"/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_{labelkey}_oct24/min_max_scale_f_")
df = bm.get_results(min_max_scale=False).transpose()
df.to_csv(f"/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_{labelkey}_oct24/min_max_scale_f_/scib_results.csv")

