library(tibble)
library(RColorBrewer)
library(stringr)
library(plyr)
library(ggplot2)
library(dplyr)

res1 = read.csv('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_disease_aug24/min_max_scale_f_/scib_results.csv', check.names = FALSE, row.names = 1)
res1 = res1[, !colnames(res1) %in% c('X_fastmnn', 'X_scVI', 'X_scanorama', 'X_combat_pca', 'X_harmony')] # remove fastmnn, scvi and scanorama
## second run with recently tested methods
res2 = read.csv('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_disease_sep24/min_max_scale_f_/scib_results.csv', check.names = FALSE, row.names = 1)
res2 = res2[, !colnames(res2) %in% c('X_scVI', 'X_pca', 'X_scmerge_pca')] # leave scdisinfect only
# third run with previously problematic methods
res3 = read.csv('/lustre/scratch127/cellgen/cellgeni/tickets/tic-3135/results_heart/scib_results_disease_oct24/min_max_scale_f_/scib_results.csv', check.names = FALSE, row.names = 1)

# when checked, it was seen that only kbet is different with very minimal difference (0.002)

res1 = t(res1)
metrics = res1['Metric Type',]
res1 = res1[1:(nrow(res1)-1),]
res2 = t(res2)
res2 = res2[1:(nrow(res2)-1),]
res3 = t(res3)
res3 = res3[1:(nrow(res3)-1),]

# concat both dataframes
res = rbind(res1, res2, res3)
rownames(res)[rownames(res) == 'res2'] = 'X_disinfact'

# add metrics row to new df (recalculate scores after removing pcr)
res0 = rbind(res, metrics)
rownames(res0)[rownames(res0) == 'metrics'] = 'Metric type'


for (i in 1:(nrow(res0)-1)){
  res0[i, colnames(res0) == 'Batch correction'] <- mean(as.numeric(res0[i, (res0['Metric type',] == 'Batch correction') & (! colnames(res0) %in% c("KBET", "iLISI", "PCR comparison"))]))
  res0[i, colnames(res0) == 'Total'] <- as.numeric(res0[i, 'Batch correction']) * 0.4 + as.numeric(res0[i, 'Bio conservation']) * 0.6
}

res0 <- res0[1:(nrow(res0)-1),]
res0 <- as.data.frame(res0)
res0 <- mutate_all(res0, function(x) as.numeric(as.character(x)))
res0 <- tibble::rownames_to_column(res0, "Method")

res0$Method <- plyr::mapvalues(res0$Method, 
                                from = c("X_harmony", "X_scVI", "X_scanorama", "X_FinVAE", "X_scpoli", "X_combat_pca", "X_fastmnn", "X_disinfact", "X_scmerge2_pca", "X_pca"), 
                                to = c("Harmony", "scVI", "Scanorama", "inVAE", "scPoli", "Combat (PCA)", "fastMNN", "scDisInFact", "scMerge (PCA)", "PCA"))

res0 <- res0[,c("Method", "Batch correction", "Bio conservation", "Total")]
res0 <- res0[order(res0$Total,decreasing=TRUE),]

metrics <- tidyr::pivot_longer(res0, cols =-c("Method"), names_to = "Testing", values_to = "Score")
metrics$Score <- round(metrics$Score, digits=3)

metrics2 <- metrics
metrics2$Testing <- str_replace_all(metrics2$Testing, " ", "\n")

library(ggplot2)

n_methods <- length(unique(metrics$Method))

# default plot
ggplot(metrics, aes(x=reorder(Method, Score), y=Score, label=Score, fill = Testing)) + geom_bar(stat='identity') + 
  geom_text(color='white', hjust = 1.2, position = position_dodge(width = .75)) +
  facet_wrap(~Testing) + coord_flip() + xlab('Score') + theme_bw() + 
  theme(legend.position = "none", text=element_text(size=15))

# more clear plot with removed grid 
ggplot(metrics2, aes(x=reorder(Method, Score), y=Score, label=Score, fill = Testing)) + geom_bar(stat='identity') + 
  facet_wrap(~~factor(Testing, c("Total", "Batch\ncorrection", "Bio\nconservation")), strip.position = "bottom", ) + coord_flip() + theme_minimal() +
  theme(legend.position = "none", text=element_text(size=15), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.y = element_text(hjust = 0), strip.text.x.bottom = element_text(angle = 90, hjust=1))

# create a column for alternating color
metrics2 <- metrics %>%
  group_by(Method) %>%
  mutate(total_score = Score[Testing == "Total"]) %>%
  ungroup() %>%
  arrange(total_score) %>%
  mutate(Method = factor(Method, levels = unique(Method)),
         row_color = if_else(as.numeric(Method) %% 2 == 0, "even", "odd"))

metrics2$Testing <- str_replace_all(metrics2$Testing, " ", "\n")


# Then, use this modified data frame in your ggplot code
ggplot(metrics2, aes(x = Method, y = Score, fill = Testing)) + 
  geom_tile(aes(y = 0.5, height = Inf, fill = row_color), alpha = 0.2) +
  geom_col(position = position_dodge(width = 0.9)) +
  facet_wrap(~factor(Testing, levels = c("Total", "Batch\ncorrection", "Bio\nconservation")), 
             strip.position = "bottom", scales = "free_x") + 
  coord_flip() + 
  theme_minimal() +
  scale_fill_manual(values = c("even" = "gray70", "odd" = "white", 
                               "Total" = "skyblue", "Batch\ncorrection" = "violet", "Bio\nconservation" = "orange"), 
                    guide = "none") +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(
    legend.position = "none", 
    text = element_text(size = 15), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0), 
    strip.text.x.bottom = element_text(angle = 90, hjust = 1),
    panel.spacing = unit(1, "lines")
  )

ggsave("invae_heart_plot_oct24_disease.pdf", width = 5, height = 10)

write.csv(res, 'scib_scores_heart_disease.csv')
write.csv(res0, 'scib_scores_heart_disease_summary.csv')


metrics