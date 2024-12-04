setwd("/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/actions")

library(scMerge)

sce1 = schard::h5ad2sce('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/objs_heart/_adata_allresults2.h5ad')

data(segList, package = 'scMerge')

scMerge2_res <- scMerge2(exprsMat = assay(sce1),
                         batch = sce1$donor_id,
                         k_celltype = 20,
                         condition = sce1$disease_short,
                         verbose = TRUE)

assay(sce1, "scMerge2") <- scMerge2_res$newY

write.csv(scMerge2_res$newY, file='scmerge2_newmatrix_sep24_noctl.csv')
saveRDS(scMerge2_res, 'obj_scmerge2_sep24_noctl.rds')

