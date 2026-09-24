## =========================
## Figure_2
## =========================
library(grid)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

load("/data/benchmarking.result.11method.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68",
                 "RaceID2" = "#3A9D54", "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3", "scCAD"  = "#1C4E73", "SCISSORS"  = "#8EC9E6")

## =========================
## Figure 2a ----
## =========================
relative_F1 <- data %>%
  group_by(Scenario, Repeat) %>%
  mutate(
    dataset_median_F1 = median(F1, na.rm = TRUE),
    relative_delta_F1 = F1 - dataset_median_F1
  ) %>%
  ungroup()


benchmark_global <- relative_F1 %>%
  group_by(Method,Scenario,SampleID,Celltype,Cellname,PrecRare,PrecRare_num) %>%
  dplyr::summarise(
    relative_delta_F1 = median(relative_delta_F1, na.rm = TRUE),
    F1 = median(F1, na.rm = TRUE),
    precision = median(precision, na.rm = TRUE),
    recall = median(recall, na.rm = TRUE),
    .groups = "drop"
  )


method_order_rel <- benchmark_global %>%
  group_by(Method) %>%
  dplyr::summarise(
    median_relative_delta_F1 =
      median(relative_delta_F1, na.rm = TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(median_relative_delta_F1)) %>%
  pull(Method)


benchmark_global <- benchmark_global %>%
  mutate(Method = factor(Method, levels = method_order_rel))


scenario_meta <- benchmark_global %>%
  distinct(Scenario, SampleID, Cellname, PrecRare, PrecRare_num) %>%
  arrange(PrecRare_num, SampleID, Cellname)

scenario_order <- scenario_meta$Scenario

delta_mat <- benchmark_global %>%
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

median_delta <- apply(delta_mat, 1, median, na.rm = TRUE)

method_order_delta <- names(sort(median_delta,decreasing = TRUE))
delta_mat    <- delta_mat[method_order_delta,]
median_delta <- median_delta[method_order_delta]


relative_F1_scenario <- relative_F1 %>%
  group_by(Scenario,SampleID,Celltype,Cellname,PrecRare,PrecRare_num,Method) %>%
  summarise(deltaF1 = median(relative_delta_F1, na.rm = TRUE),
            .groups = "drop")
    
method_order <- relative_F1_scenario %>%
  group_by(Method) %>%
  summarise(
    median_delta = median(deltaF1, na.rm = TRUE),
    mean_delta   = mean(deltaF1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(median_delta, mean_delta) %>% 
  pull(Method)

df <- relative_F1_scenario

df$Method <- factor(df$Method, levels = method_order)

pairwise_adjacent <- data.frame()

for(i in 1:(length(method_order)-1)){
  
  m1 <- method_order[i]
  m2 <- method_order[i+1]
  
  tmp <- df %>%
    filter(Method %in% c(m1, m2)) %>%
    dplyr::select(Scenario, Method, deltaF1) %>%
    pivot_wider(
      names_from = Method,
      values_from = deltaF1
    )
  
  wt <- wilcox.test(
    tmp[[m1]],
    tmp[[m2]],
    paired = TRUE
  )
  
  pairwise_adjacent <- rbind(
    pairwise_adjacent,
    data.frame(
      group1 = m1,
      group2 = m2,
      p = wt$p.value
    )
  )
}

pairwise_adjacent$label <- case_when(
  pairwise_adjacent$p < 0.001 ~ "***",
  pairwise_adjacent$p < 0.01  ~ "**",
  pairwise_adjacent$p < 0.05  ~ "*",
  TRUE                            ~ "ns"
)

pairwise_adjacent <- pairwise_adjacent %>%
  mutate(
    xmin = group1,
    xmax = group2,
    y.position = seq(
      max(df$deltaF1, na.rm = TRUE) + 0.02,
      by = 0.03,
      length.out = n()
    )
  )

method_p_label <- rep("-", length(method_order_delta))
names(method_p_label) <- method_order_delta

for(i in 2:length(method_order_delta)){
  
  m1 <- method_order_delta[i-1]
  m2 <- method_order_delta[i]
  
  tmp <- pairwise_adjacent %>%
    filter(
      (group1 == m1 & group2 == m2) |
        (group1 == m2 & group2 == m1)
    )
  
  if(nrow(tmp) > 0){
    method_p_label[m2] <- tmp$label
  }
}


right_anno <- rowAnnotation(
  
  "Median ΔF1" = anno_points(
    median_delta,
    pch = 18,
    size = unit(4, "mm"),
    gp = gpar(col = "#8B0000"),
    axis = TRUE
  ),
  
  "Adjacent\ncomparison" = anno_text(
    method_p_label,
    gp = gpar(fontsize = 10),
    just = "left"
  ),
  
  width = unit(3.5, "cm")
)

pdf("/results/Fig2a.heatmap.deltaF1.sig.pdf", width = 12, height = 8)
Heatmap(delta_mat,
        name = "ssr-F1",
        #rect_gp = grid::gpar(col = "white", lwd = 0.05),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_split = col_group,
        column_gap = unit(3, "mm"),
        bottom_annotation = ha,
        right_annotation = right_anno,
        row_names_side = "left",
        row_names_gp = grid::gpar(fontsize = 10),
        row_names_max_width = unit(10, "cm"),
        col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")))
dev.off()


## =========================
## Figure 2b ----
## =========================
overall_rank <- relative_F1_scenario %>%
  group_by(Method) %>%
  summarise(
    median_delta = median(deltaF1, na.rm = TRUE),
    mean_delta   = mean(deltaF1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_delta), desc(mean_delta)) %>%
  mutate(OverallRank = row_number())


rank_by_prop <- relative_F1_scenario %>%
  group_by(PrecRare, Method) %>%
  summarise(
    median_delta = median(deltaF1, na.rm = TRUE),
    mean_delta   = mean(deltaF1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(PrecRare) %>%
  arrange(
    desc(median_delta),
    desc(mean_delta),
    .by_group = TRUE
  ) %>%
  mutate(Rank = row_number()) %>%
  ungroup()


rank_overall <- overall_rank %>%
  transmute(
    Method,
    Rank = OverallRank,
    Group = "Overall"
  )

rank_prop <- rank_by_prop %>%
  transmute(
    Method,
    Rank,
    Group = PrecRare
  )

rank_all <- bind_rows(
  rank_overall,
  rank_prop
)

rank_matrix <- rank_all %>%
  dplyr::select(Method, Group, Rank) %>%
  pivot_wider(
    names_from = Group,
    values_from = Rank
  ) %>%
  tibble::column_to_rownames("Method") %>%
  as.matrix()

rank_matrix <- rank_matrix[, c("Overall","0.125%","0.25%","0.5%","1%","2%","3%","4%","5%")]
cor_mat <- cor(rank_matrix, method = "spearman")

cor_plot <- cor_mat
cor_plot[upper.tri(cor_plot)] <- NA

pdf("/results/Fig2b.heatmap.prop.corr.pdf", width = 7, height = 5)
Heatmap(cor_plot,
        name = "Spearman correlation coefficient (\u03c1)",
        na_col = "white",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "white"),
        row_names_side = "left",
        column_names_side = "bottom",
        column_names_rot = 45,
        col = colorRamp2(
          c(0, 1),
          c("white", "#B2182B")
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (!is.na(cor_plot[i, j])) {
            grid.text(
              sprintf("%.2f", cor_plot[i, j]),
              x, y,
              gp = gpar(fontsize = 9)
            )
          }
        }
)
dev.off()

