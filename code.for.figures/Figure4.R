## =========================
## Figure_4
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68","RaceID2" = "#3A9D54",
                 "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3","scCAD_1%"  = "#1C4E73", "scCAD_5%" = "#3E7FBF", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 4a ----
## =========================
data$SampleID_Cellname <- paste(data$SampleID, data$Cellname, sep = "_")

benchmark_global <- data %>%
  group_by(Method, Scenario, SampleID, Celltype, Cellname, PrecRare, PrecRare_num, SampleID_Cellname) %>%
  dplyr::summarise(F1        = median(F1, na.rm = TRUE),
                   MCC       = median(MCC, na.rm = TRUE),
                   precision = median(precision, na.rm = TRUE),
                   recall    = median(recall, na.rm = TRUE),
                   .groups = "drop")

dataset_difficulty <- benchmark_global %>%
  group_by(Scenario, SampleID_Cellname, SampleID, Cellname) %>%
  dplyr::summarise(
    scenario_median_F1 = median(F1, na.rm = TRUE),
    scenario_best_F1   = max(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SampleID_Cellname, SampleID) %>%
  dplyr::summarise(median_scenario_F1 = median(scenario_median_F1, na.rm = TRUE),
                   median_best_F1     = median(scenario_best_F1, na.rm = TRUE),
                   n_scenarios        = n(),
                   .groups = "drop") %>%
  
  arrange(median_scenario_F1)

dataset_levels <- sort(unique(benchmark_global$SampleID))
dataset_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(dataset_levels)),dataset_levels)

similarity_meta <- data %>%
  group_by(SampleID_Cellname, SampleID, Cellname) %>%
  dplyr::summarise(
    Rare_cell_similarity = median(Rare_cell_similarity, na.rm = TRUE),
    .groups = "drop")

difficulty_similarity <- dataset_difficulty %>%
  left_join(similarity_meta, by = "SampleID_Cellname") %>%
  filter(!is.na(median_scenario_F1),
         !is.na(Rare_cell_similarity),
         is.finite(median_scenario_F1),
         is.finite(Rare_cell_similarity))

cor_spearman <- cor.test(difficulty_similarity$Rare_cell_similarity,difficulty_similarity$median_scenario_F1,method = "spearman")

rho  <- round(cor_spearman$estimate, 3)
pval <- signif(cor_spearman$p.value, 3)

p1 <- ggplot(difficulty_similarity, aes(x = Rare_cell_similarity,y = median_scenario_F1)) +
  geom_point(aes(fill = SampleID.x),shape = 21,  color = "grey30", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "#2C3E50", fill = "#BFC9D4", linewidth = 0.8, alpha = 0.25) +
  #geom_text_repel(aes(label = SampleID), size = 3, max.overlaps = 20) +
  annotate("text",
           x = Inf,
           y = Inf,
           hjust = 1.1,
           vjust = 1.5,
           label = paste0("Spearman rho = ", rho, "\np = ", pval),
           size = 4) +
  scale_fill_manual(values = dataset_colors) +
  labs(x = "Rare-cell similarity", y = "Median scenario-level F1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"))
p1    


## =========================
## Figure 4b ----
## =========================
cell_class_map <- tribble(
  ~Cellname, ~Cell_class,
  
  ## -------------------------
  ## Epithelial
  ## -------------------------
  "Alveolar Epithelial Type 2", "Epithelial",
  "basal cell", "Epithelial",
  "basal cell of prostate epithelium", "Epithelial",
  "bladder urothelial cell", "Epithelial",
  "club cell", "Epithelial",
  "colon epithelial cell", "Epithelial",
  "epithelial cell", "Epithelial",
  "epithelial cell of proximal tubule", "Epithelial",
  "intestine goblet cell", "Epithelial",
  "kidney collecting duct intercalated cell", "Epithelial",
  "kidney collecting duct principal cell", "Epithelial",
  "kidney connecting tubule epithelial cell", "Epithelial",
  "kidney distal convoluted tubule epithelial cell", "Epithelial",
  "kidney loop of Henle thick ascending limb epithelial cell", "Epithelial",
  "kidney loop of Henle thin ascending limb epithelial cell", "Epithelial",
  "luminal cell of prostate epithelium", "Epithelial",
  "luminal epithelial cell of mammary gland", "Epithelial",
  "mammary gland epithelial cell", "Epithelial",
  "secretory cell", "Epithelial",
  
  ## -------------------------
  ## Endothelial
  ## -------------------------
  "blood vessel endothelial cell", "Endothelial",
  "Capillary", "Endothelial",
  "endothelial cell", "Endothelial",
  "endothelial cell of coronary artery", "Endothelial",
  
  ## -------------------------
  ## Stromal / mesenchymal / muscle-associated
  ## -------------------------
  "bronchial smooth muscle cell", "Stromal",
  "fibroblast", "Stromal",
  "fibroblast of cardiac tissue", "Stromal",
  "fibroblast of connective tissue of prostate", "Stromal",
  "granulosa cell", "Stromal",
  "kidney interstitial cell", "Stromal",
  "myofibroblast cell", "Stromal",
  "pericyte", "Stromal",
  "regular ventricular cardiac myocyte", "Stromal",
  "smooth muscle cell", "Stromal",
  "stromal cell", "Stromal",
  "theca cell", "Stromal",
  "vascular associated smooth muscle cell", "Stromal",
  
  ## -------------------------
  ## Lymphoid immune
  ## -------------------------
  "B", "Lymphoid immune",
  "B cell", "Lymphoid immune",
  "CD4-positive, alpha-beta memory T cell", "Lymphoid immune",
  "CD4-positive, alpha-beta T cell", "Lymphoid immune",
  "CD8 T", "Lymphoid immune",
  "CD8-positive, alpha-beta memory T cell", "Lymphoid immune",
  "CD8+ Naive T", "Lymphoid immune",
  "gamma-delta T cell", "Lymphoid immune",
  "IgA plasma cell", "Lymphoid immune",
  "leukocyte", "Lymphoid immune",
  "mature NK T cell", "Lymphoid immune",
  "memory B cell", "Lymphoid immune",
  "Memory CD4 T", "Lymphoid immune",
  "Naive CD4 T", "Lymphoid immune",
  "Natural Killer", "Lymphoid immune",
  "natural killer cell", "Lymphoid immune",
  "NK", "Lymphoid immune",
  "plasma cell", "Lymphoid immune",
  "T cell", "Lymphoid immune",
  "native cell", "Lymphoid immune",
  
  ## -------------------------
  ## Myeloid immune
  ## -------------------------
  "CD14+ Mono", "Myeloid immune",
  "classical monocyte", "Myeloid immune",
  "Classical Monocyte", "Myeloid immune",
  "erythrocyte", "Myeloid immune",
  "FCGR3A+ Mono", "Myeloid immune",
  "intermediate monocyte", "Myeloid immune",
  "macrophage", "Myeloid immune",
  "Macrophage", "Myeloid immune",
  "monocyte", "Myeloid immune",
  "neutrophil", "Myeloid immune",
  "non-classical monocyte", "Myeloid immune",
  "phagocyte", "Myeloid immune",
  "tissue-resident macrophage", "Myeloid immune"
)



difficulty_similarity_class <- difficulty_similarity %>%
  left_join(cell_class_map, by = "Cellname")



cellclass_color <- c("Epithelial"      = "#4C78A8",
                     "Endothelial"    = "#72B7B2",
                     "Stromal"        = "#54A24B",
                     "Lymphoid immune"= "#E45756",
                     "Myeloid immune" = "#F58518")

cellclass_order_f1 <- difficulty_similarity_class %>%
  group_by(Cell_class) %>%
  dplyr::summarise(
    median_F1 = median(median_scenario_F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_F1) %>% 
  pull(Cell_class)

difficulty_similarity_class <- difficulty_similarity_class %>%
  mutate(Cell_class = factor(Cell_class, levels = cellclass_order_f1))

p2 <- ggplot(difficulty_similarity_class,
             aes(x = Cell_class, y = median_scenario_F1, fill = Cell_class)) +
  geom_boxplot(width = 0.65,alpha = 0.65,linewidth = 0.35) +
  geom_point(aes(color = Cell_class),size = 1.8,alpha = 0.75,
             position = position_jitter(width = 0.15, height = 0)) +
  scale_fill_manual(values = cellclass_color) +
  scale_color_manual(values = cellclass_color) +
  ylim(c(0,1))+
  labs(x = "", y = "Median scenario-level F1") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12, angle = 35, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
p2


## =========================
## Figure S5a ----
## =========================
p3 <- ggplot(difficulty_similarity_class,
             aes(x = Rare_cell_similarity, y = median_scenario_F1)) +
  geom_point(aes(fill = Cell_class),alpha = 0.75, shape = 21,color = "grey30",size = 2.8) +
  geom_smooth(method = "lm",se = TRUE,color = "#2C3E50",fill = "#BFC9D4",linewidth = 0.8,alpha = 0.25) +
  annotate("text",x = Inf,
           y = Inf,
           hjust = 1.1,
           vjust = 1.5,
           label = paste0("Spearman rho = ", rho, "\np = ", pval),
           size = 4) +
  scale_fill_manual(values = cellclass_color) +
  labs(x = "Rare-cell similarity",y = "Median scenario-level F1",fill = "Cell class") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 11, color = "black"))
p3


## =========================
## Figure S5b ----
## =========================
p4 <- ggplot(difficulty_similarity_class,
             aes(x = Cell_class, y = Rare_cell_similarity, fill = Cell_class)) +
  geom_boxplot(width = 0.65,alpha = 0.65,linewidth = 0.35) +
  geom_point(aes(color = Cell_class),size = 1.8,alpha = 0.75,
             position = position_jitter(width = 0.15, height = 0)) +
  scale_fill_manual(values = cellclass_color) +
  scale_color_manual(values = cellclass_color) +
  labs(x = "",y = "Rare-cell similarity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12, angle = 35, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
p4


## =========================
## Figure 4c ----
## =========================
plot_data <- data %>%
  filter(SampleID == "PBMC3k") %>%  
  group_by(Cellname, Method, PrecRare) %>%
  dplyr::summarise(median_F1 = median(F1, na.rm = TRUE),
                   Rare_cell_similarity = unique(Rare_cell_similarity),
                   .groups = "drop")

plot_data$Cellname <- reorder(plot_data$Cellname, plot_data$Rare_cell_similarity)

p5 <- ggplot(plot_data, aes(x = Cellname, y = median_F1, fill = Rare_cell_similarity)) +
  geom_col(width = 0.7) +
  scale_fill_distiller(palette = "Spectral", direction = -1, name = "Rare cell similarity") +
  facet_grid(rows = vars(PrecRare), cols = vars(Method), scales = "free_x") + 
  theme_classic() +
  labs(x = "", y = "Median F1 score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        strip.text.x = element_text(color = "black", size = 10),           
        strip.text.y = element_text(color = "black", angle = 0, size = 10), 
        strip.background = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(1.0, "lines"),
        plot.title = element_text(color = "black", hjust = 0.5))

p5


