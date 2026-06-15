## =========================
## Figure_2
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68","RaceID2" = "#3A9D54",
                 "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3","scCAD_1%"  = "#1C4E73", "scCAD_5%" = "#3E7FBF", "SCISSORS"  = "#8EC9E6")


## =========================
## Figure 2 ----
## =========================
benchmark_global <- data %>%
  group_by(Method, Scenario, SampleID, Celltype, Cellname, PrecRare, PrecRare_num) %>%
  dplyr::summarise(F1        = median(F1, na.rm = TRUE),
                   precision = median(precision, na.rm = TRUE),
                   recall    = median(recall, na.rm = TRUE),
                   .groups = "drop")

relative_F1 <- benchmark_global %>%
  group_by(Scenario) %>%
  mutate(scenario_median_F1 = median(F1, na.rm = TRUE),
         scenario_mad_F1    = mad(F1, constant = 1, na.rm = TRUE),
         relative_delta_F1  = F1 - scenario_median_F1,
         robust_zF1 = case_when(scenario_mad_F1 > 0 ~ 0.6745 * (F1 - scenario_median_F1) / scenario_mad_F1,
                                scenario_mad_F1 == 0 & F1 == scenario_median_F1 ~ 0,
                                TRUE ~ NA_real_),
         above_median = relative_delta_F1 > 0) %>%
  ungroup()

method_order_rel <- relative_F1 %>%
  group_by(Method) %>%
  dplyr::summarise(median_relative_delta_F1 = median(relative_delta_F1, na.rm = TRUE),
                   .groups = "drop") %>%
  arrange(desc(median_relative_delta_F1)) %>%
  pull(Method)

relative_F1 <- relative_F1 %>%
  mutate(Method = factor(Method, levels = method_order_rel))


scenario_meta <- relative_F1 %>%
  distinct(Scenario, SampleID, Cellname, PrecRare, PrecRare_num) %>%
  arrange(PrecRare_num, SampleID, Cellname)

scenario_order <- scenario_meta$Scenario

delta_mat <- relative_F1 %>%
  dplyr::select(Method, Scenario, relative_delta_F1) %>%
  pivot_wider(names_from = Scenario, values_from = relative_delta_F1) %>%
  tibble::column_to_rownames("Method") %>%
  as.matrix()

delta_mat <- delta_mat[intersect(method_order_rel, rownames(delta_mat)),scenario_order]

col_group <- scenario_meta$PrecRare
dataset_group <- scenario_meta$SampleID

prop_levels <- unique(col_group)
dataset_levels <- unique(dataset_group)

ha <- HeatmapAnnotation(
  Dataset = dataset_group,
  Proportion = col_group,
  col = list(Dataset = structure(colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(dataset_levels)),
                                 names = dataset_levels),
             Proportion = structure(RColorBrewer::brewer.pal(length(prop_levels), "Set3"),
                                    names = prop_levels)),
  annotation_name_gp = grid::gpar(fontsize = 9),
  simple_anno_size = unit(0.35, "cm"),
  gap = unit(1, "mm"))

lim <- quantile(abs(delta_mat), 0.95, na.rm = TRUE)

Heatmap(delta_mat,
        name = "ΔF1",
        #rect_gp = grid::gpar(col = "white", lwd = 0.05),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_split = col_group,
        column_gap = unit(3, "mm"),
        bottom_annotation = ha,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = 10),
        row_names_max_width = unit(10, "cm"),
        col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")))


## =========================
## Figure S2a ----
## =========================
wilcox_relative <- relative_F1 %>% 
  group_by(Method) %>% 
  group_modify(~ {
    x <- .x$relative_delta_F1 
    x <- x[!is.na(x)] 
    
    if (length(x) < 3) {
      return(tibble( n_scenario = length(x), 
                     median_relative_delta_F1 = median(x, na.rm = TRUE), 
                     p_value = NA_real_ )) 
    } 
    if (all(x == 0)) {
      p_val <- 1 
    }
    else {
      p_val <- wilcox.test(x, mu = 0, alternative = "two.sided", exact = FALSE )$p.value 
    } 
    tibble( n_scenario = length(x), 
            median_relative_delta_F1 = median(x, na.rm = TRUE), 
            p_value = p_val ) 
  }) %>% 
  ungroup() %>% 
  mutate( p_adj_BH = p.adjust(p_value, method = "BH"), 
          significance = case_when( 
            p_adj_BH < 0.001 ~ "***", 
            p_adj_BH < 0.01 ~ "**", 
            p_adj_BH < 0.05 ~ "*", 
            TRUE ~ "ns" ) ) %>% 
  arrange(desc(median_relative_delta_F1)) 

wilcox_relative

relative_F1_plot <- relative_F1 %>%
  mutate(Method = factor(Method, levels = method_order_rel))

wilcox_anno <- wilcox_relative %>%
  mutate(Method = factor(Method, levels = method_order_rel),
         label = significance) %>%
  left_join(relative_F1_plot %>%
              group_by(Method) %>%
              summarise(y_pos = max(relative_delta_F1, na.rm = TRUE),.groups = "drop"),
            by = "Method") %>%
  mutate(y_pos = y_pos + 0.03)

p1 <- ggplot(relative_F1_plot, aes(x = Method, y = relative_delta_F1, fill = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(width = 0.7, alpha = 0.8, outlier.size = 0.4) +
  geom_point(aes(color = Method),alpha = 0.18,size = 0.35,position = position_jitter(width = 0.15)) +
  geom_text(data = wilcox_anno,
            aes(x = Method, y = y_pos, label = label),
            inherit.aes = FALSE,
            size = 5) +
  scale_fill_manual(values = Methodcolor) +
  scale_color_manual(values = Methodcolor) +
  labs(y = expression(Delta * "F1 relative to scenario median"),x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
p1  



## =========================
## Figure S2b ----
## =========================
above_median_df <- relative_F1 %>%
  group_by(Method) %>%
  dplyr::summarise(above_median_rate = mean(above_median, na.rm = TRUE),
                   n_above = sum(above_median, na.rm = TRUE),
                   n_scenario = sum(!is.na(relative_delta_F1)),
                   .groups = "drop") %>%
  arrange(desc(above_median_rate)) %>%
  mutate(Method = factor(Method, levels = Method))

p2 <- ggplot(above_median_df, aes(x = Method, y = above_median_rate, fill = Method)) +
  geom_col(color = "black",alpha = 0.9, linewidth = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = paste0(round(above_median_rate * 100, 1), "%")),vjust = -0.3,size = 4) +
  scale_fill_manual(values = Methodcolor) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(y = "Proportion of scenarios above median F1",x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"))
p2    


