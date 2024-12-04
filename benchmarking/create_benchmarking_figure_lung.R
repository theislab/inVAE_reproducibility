library(tibble)
library(RColorBrewer)
library(stringr)
library(plyr)
library(ggplot2)
library(dplyr)

res_name = strsplit(res,'_')[[1]]
res_name = paste(res_name[3:length(res_name)],collapse = '_')

## first run with most methods
res0 = read.csv(paste0('results_lung/run1/scib_results_newcelltype_stage/min_max_scale_f_/scib_results.csv'), check.names = FALSE, row.names = 1)
res0 = t(res0)
metrics = res0['Metric Type',]
res0 = res0[1:(nrow(res0)-1),]

res1 = read.csv(paste0('results_lung/run2/scib_results_newcelltype_stage/min_max_scale_f_/scib_results.csv'), check.names = FALSE, row.names = 1)
res1 = t(res1)
res1 = res1[1:(nrow(res1)-1),]

res0 = rbind(res0, res1)
rownames(res0)[rownames(res0) == 'res1'] = 'X_FinVAE'

# add metrics row to new df (recalculate scores after removing pcr)
res0 = rbind(res0, metrics)
rownames(res0)[rownames(res0) == 'metrics'] = 'Metric Type'

for (i in 1:(nrow(res0)-1)){
  res0[i, colnames(res0) == 'Batch correction'] <- mean(as.numeric(res0[i, (res0['Metric Type',] == 'Batch correction') & (! colnames(res0) %in% c("KBET", "iLISI", "PCR comparison"))]))
  res0[i, colnames(res0) == 'Total'] <- as.numeric(res0[i, 'Batch correction']) * 0.4 + as.numeric(res0[i, 'Bio conservation']) * 0.6
}

res0 <- res0[1:(nrow(res0)-1),]
res0 <- as.data.frame(res0)
res0 <- mutate_all(res0, function(x) as.numeric(as.character(x)))
res0 <- tibble::rownames_to_column(res0, "Method")

res0$Method <- plyr::mapvalues(res0$Method, 
                                from = c("X_scVI", "X_FinVAE", "X_scpoli", "X_combat_pca", "X_pca"), 
                                to = c("scVI", "inVAE", "scPoli", "Combat (PCA)", "PCA"))

res0 <- res0[,c("Method", "Batch correction", "Bio conservation", "Total")]
res0 <- res0[order(res0$Total,decreasing=TRUE),]

metrics <- tidyr::pivot_longer(res0, cols =-c("Method"), names_to = "Testing", values_to = "Score")
metrics$Score <- round(metrics$Score, digits=3)

metrics2 <- metrics
metrics2$Testing <- str_replace_all(metrics2$Testing, " ", "\n")

library(ggplot2)

n_methods <- length(unique(metrics$Method))
  
# create a column for alternating color
metrics2 <- metrics %>%
  group_by(Method) %>%
  mutate(total_score = Score[Testing == "Total"]) %>%
  ungroup() %>%
  arrange(total_score) %>%
  mutate(Method = factor(Method, levels = unique(Method)),
          row_color = if_else(as.numeric(Method) %% 2 == 0, "even", "odd"))

metrics2$Testing <- str_replace_all(metrics2$Testing, " ", "\n")

# Then, use this modified data frame in ggplot code
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

ggsave(paste0('invae_lung_plot_',res_name,'_oct24.pdf'), width = 5, height = 10)

write.csv(res, paste0('scib_scores_lung_',res_name,'.csv'))
write.csv(res0, paste0('scib_scores_lung_',res_name,'_summary.csv'))
