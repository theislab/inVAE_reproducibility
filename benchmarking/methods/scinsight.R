setwd("/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/actions")

library(scINSIGHT)

sce1 = schard::h5ad2sce('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/_adata_allresults2.h5ad')
test = readRDS('sim_count.rds')

donors_all = sce1$donor_id
donors = unique(donors_all)

donor_mats <- list()
for (i in 1:length(donors)){
  donor_mats[[i]] <- as.matrix(assay(sce1)[,donors_all == donors[i]])
}

donors_disease <- data.frame(sce1$donor_id, sce1$disease)
donors_disease <- donors_disease[!duplicated(donors_disease),]

names(donor_mats) = donors
condition = donors_disease$sce1.disease
names(condition) = donors

print(paste("Length of donor_mats:", length(donor_mats)))
print(paste("Length of condition:", length(condition)))

scobj = create_scINSIGHT(norm.data = donor_mats, condition = condition)

scobj = run_scINSIGHT(scobj, out.dir = './scinsight_out/', num.cores = 8)

saveRDS(scobj, 'scinsight_obj.rds')
