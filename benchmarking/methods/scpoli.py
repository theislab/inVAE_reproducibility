import scanpy as sc
from scarches.models.scpoli import scPoli

adata = sc.read('/lustre/scratch126/cellgen/team292/ha10/data/Heart_Atlas/adata_Heart_Reichart_HV_train.h5ad')

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 20,
    "reduce_lr": True,
    "lr_patience": 13,
    "lr_factor": 0.1,
}

scpoli_model = scPoli(
    adata=adata,
    condition_keys='donor_id',
    cell_type_keys='disease',
    embedding_dims=5,
    recon_loss='nb',
)
scpoli_model.train(
    n_epochs=50,
    pretraining_epochs=40,
    early_stopping_kwargs=early_stopping_kwargs,
    eta=5,
)

adata.obsm['X_scpoli'] = scpoli_model.get_latent(
    adata,
    mean=True
)

adata.write('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/adata_scpoli_notfinal.h5ad', compression='gzip')

