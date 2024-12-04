library(schard)
library(batchelor)

sce1 = schard::h5ad2sce('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/objs_heart/_adata_allresults2.h5ad')

sce1

f.out1 <- fastMNN(sce1, batch=sce1$donor_id, assay.type='X', d=25)

saveRDS(f.out1, '/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/objs_heart/_sce_batchelor_oct24.rds')

# sce1 <- readRDS('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/_sce_batchelor.rds')

write.csv(reducedDim(f.out1, 'corrected'), file='/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/objs_heart/fastmnn_dims_oct24.csv')
