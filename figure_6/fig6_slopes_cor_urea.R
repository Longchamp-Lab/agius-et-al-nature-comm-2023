library(tidyverse)
library(reshape2)
library(openxlsx)

urea <- read.csv("urea.csv")
counts <- read.csv("counts_logCPM_grouped.csv")

## delete useless column
counts <- counts[,-c(2:5)]

## melt colums and create new group column
counts.melt <- melt(counts,id.vars = c("gene"),measure.vars = c("avg.Ctrl","avg.Ctrl.IGF1","avg.HC","avg.LP","avg.LPHC","avg.LPHC.IGF1"))

urea.mean <- urea %>%
  group_by(group) %>%
  summarise_at(vars(AUC.urea), list(mean.AUC.urea = mean))

## create group column
counts.melt$group[counts.melt$variable == "avg.Ctrl"] = "Ctrl"
counts.melt$group[counts.melt$variable == "avg.HC"] = "HC"
counts.melt$group[counts.melt$variable == "avg.LP"] = "LP"
counts.melt$group[counts.melt$variable == "avg.LPHC"] = "LPHC"
counts.melt$group[counts.melt$variable == "avg.LPHC.IGF1"] = "LPHC.IGF1"
counts.melt$group[counts.melt$variable == "avg.Ctrl.IGF1"] = "Ctrl.IGF1"

## add AUC urea depending on the group
counts.melt$urea[counts.melt$group == "Ctrl"] <- urea.mean$mean.AUC.urea[urea.mean$group == "Ctrl"]
counts.melt$urea[counts.melt$group == "HC"] <- urea.mean$mean.AUC.urea[urea.mean$group == "HC"]
counts.melt$urea[counts.melt$group == "LP"] <- urea.mean$mean.AUC.urea[urea.mean$group == "LP"]
counts.melt$urea[counts.melt$group == "LPHC"] <- urea.mean$mean.AUC.urea[urea.mean$group == "LPHC"]
counts.melt$urea[counts.melt$group == "Ctrl.IGF1"] <- urea.mean$mean.AUC.urea[urea.mean$group == "Ctrl-IGF1"]
counts.melt$urea[counts.melt$group == "LPHC.IGF1"] <- urea.mean$mean.AUC.urea[urea.mean$group == "LPHC-IGF1"]

## export
#write.csv(DEG_all, "DEG.all.csv")

## Calculate pearson r2 for logFC and FDR vs protein intake
counts.melt$slope = 0
counts.melt$r2 = 0
counts.melt$up.down = ""

counts.melt_list <- counts.melt %>% group_split(gene)

for(i in 1:length(counts.melt_list)){
  counts.melt_list[[i]]$r2 <- abs(cor(counts.melt_list[[i]]$value, counts.melt_list[[i]]$urea, use = "all.obs")^2)
  
  res <- lm(counts.melt_list[[i]]$value~counts.melt_list[[i]]$urea)
  counts.melt_list[[i]]$slope <- as.numeric(abs(res$coefficients[2]))
  
  coef <- as.numeric(res$coefficients[2])
  if ( coef >= 0) {
    counts.melt_list[[i]]$up.down <- "up"
  } else {
    counts.melt_list[[i]]$up.down <- "down"
  }
}

counts.final <- combine(counts.melt_list)

write.xlsx(counts.final, "counts.final_slope_r2.xlsx")

counts.final.sum <- counts.final %>%
  group_by(gene, up.down) %>%
  summarise_at(vars(slope, r2), mean)

write.xlsx(counts.final.sum, "gene_slope_r2.xlsx")

################################################################################
## Scatter plots up and down
################################################################################
## Graph ggpaired for cluster representation
library(ggpubr)
library("ggsci")
library(ragg)

relevant <- read.xlsx("gene_slope_r2.xlsx" , sheet = "relevant")

relevant_up <- relevant %>% filter(up.down == "down")
relevant_down <- relevant %>% filter(up.down == "up")

## filter relevant genes from counts
geneList.up <- counts %>% filter(gene %in% relevant_up$gene)
geneList.down <- counts %>% filter(gene %in% relevant_down$gene)

## delete useless column
geneList.up <- geneList.up[,-c(2:5)]
geneList.down <- geneList.down[,-c(2:5)]

## melt colums and create new group column
geneList.up <- melt(geneList.up,id.vars = c("gene"),measure.vars = c("avg.Ctrl","avg.Ctrl.IGF1","avg.HC","avg.LP","avg.LPHC","avg.LPHC.IGF1"))
geneList.down <- melt(geneList.down,id.vars = c("gene"),measure.vars = c("avg.Ctrl","avg.Ctrl.IGF1","avg.HC","avg.LP","avg.LPHC","avg.LPHC.IGF1"))

## create group colunm
geneList.up$group[geneList.up$variable == "avg.Ctrl"] = "Ctrl"
geneList.up$group[geneList.up$variable == "avg.HC"] = "HC"
geneList.up$group[geneList.up$variable == "avg.LP"] = "LP"
geneList.up$group[geneList.up$variable == "avg.LPHC"] = "LPHC"
geneList.up$group[geneList.up$variable == "avg.LPHC.IGF1"] = "LPHC.IGF1"
geneList.up$group[geneList.up$variable == "avg.Ctrl.IGF1"] = "Ctrl.IGF1"

geneList.down$group[geneList.down$variable == "avg.Ctrl"] = "Ctrl"
geneList.down$group[geneList.down$variable == "avg.HC"] = "HC"
geneList.down$group[geneList.down$variable == "avg.LP"] = "LP"
geneList.down$group[geneList.down$variable == "avg.LPHC"] = "LPHC"
geneList.down$group[geneList.down$variable == "avg.LPHC.IGF1"] = "LPHC.IGF1"
geneList.down$group[geneList.down$variable == "avg.Ctrl.IGF1"] = "Ctrl.IGF1"

## with GGplot
agg_png("genes_protein_dilution_plot_Up.png", width = 6, height = 8, res = 300, units = "cm", background = "transparent")
geneList.up %>% 
  mutate(group=factor(group)) %>% 
  mutate(group=fct_relevel(group,c("Ctrl.IGF1","LPHC.IGF1","Ctrl","HC","LP","LPHC"))) %>%
  dplyr::arrange(gene, group) %>% 
  dplyr::group_by(gene) %>% 
  ggplot2::ggplot(., aes(y = value, x = group), show.legend = FALSE) +
  geom_path(
    aes(group = gene), 
    size = 0.1, alpha = 0.6, 
    position = position_jitter(width = 0.1, seed = 3922)
  ) + 
  geom_point(aes(color= group), alpha = 0.6, position = position_jitter(width = 0.1, seed = 3922), show.legend = FALSE) +  
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot") + 
  #scale_y_continuous(breaks = seq(200, 650, 50)) + 
  #coord_cartesian(ylim = c(200, 650)) + 
  labs(y = "Average mRNA expression (logCPM)") +
  scale_color_npg() +
  theme_light() + 
  theme(
    axis.title = element_text(face = "bold"), 
    legend.title = element_text(face = "bold"),
    axis.line = element_line(colour = 'black', size = .25),
    axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = 'italic', size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent')
    #aspect.ratio = 3/1.6
  )
invisible(dev.off())

## with GGplot
agg_png("genes_protein_dilution_plot_Down.png", width = 6, height = 8, res = 300, units = "cm", background = "transparent")
geneList.down %>% 
  mutate(group=factor(group)) %>% 
  mutate(group=fct_relevel(group,c("Ctrl.IGF1","LPHC.IGF1","Ctrl","HC","LP","LPHC"))) %>%
  dplyr::arrange(gene, group) %>% 
  dplyr::group_by(gene) %>% 
  ggplot2::ggplot(., aes(y = value, x = group), show.legend = FALSE) +
  geom_path(
    aes(group = gene), 
    size = 0.1, alpha = 0.6, 
    position = position_jitter(width = 0.1, seed = 3922)
  ) + 
  geom_point(aes(color= group), alpha = 0.6, position = position_jitter(width = 0.1, seed = 3922), show.legend = FALSE) +  
  stat_summary(geom = "pointrange", fun.data = "mean_cl_boot") + 
  #scale_y_continuous(breaks = seq(200, 650, 50)) + 
  #coord_cartesian(ylim = c(200, 650)) + 
  labs(y = "") +
  scale_color_npg() +
  theme_light() + 
  theme(
    axis.title = element_text(face = "bold"), 
    legend.title = element_text(face = "bold"),
    axis.line = element_line(colour = 'black', size = .25),
    axis.title.x=element_blank(),
    #axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = 'italic', size = 11),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    legend.background = element_rect(fill='transparent')
    #aspect.ratio = 3/1.6
  )
invisible(dev.off())