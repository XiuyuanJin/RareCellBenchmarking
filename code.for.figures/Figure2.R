## =========================
## Figure_2
## =========================

setwd("/home/jinxiuyuan/Proj_scCellFishing/paper/section2/figure")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.11method.RData")

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
  summarise(
    deltaF1 = median(relative_delta_F1, na.rm = TRUE),
    .groups = "drop"
  )

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
    select(Scenario, Method, deltaF1) %>%
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

pairwise_adjacent$p.adj <- p.adjust(pairwise_adjacent$p, method = "BH")

pairwise_adjacent$label <- case_when(
  pairwise_adjacent$p.adj < 0.001 ~ "***",
  pairwise_adjacent$p.adj < 0.01  ~ "**",
  pairwise_adjacent$p.adj < 0.05  ~ "*",
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

pdf("Fig2a.heatmap.deltaF1.sig.pdf", width = 12, height = 8)
Heatmap(delta_mat,
        name = "median ΔF1",
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
## Figure S1a-h ----
## =========================
plot_deltaF1_box <- function(prop){
  
  df <- relative_F1_scenario %>% filter(PrecRare == prop)
  
  method_order <- df %>%
    group_by(Method) %>%
    dplyr::summarise(
      median_delta = median(deltaF1, na.rm = TRUE),
      mean_delta   = mean(deltaF1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_delta,mean_delta) %>%
    pull(Method)
  
  df$Method <- factor(df$Method, levels = method_order)
  
  p <- ggplot(df, aes(Method, deltaF1, fill = Method)) +
    geom_boxplot(width = 0.7, outlier.size = 0.3) +
    stat_summary(fun = median,geom = "point",shape = 18,size = 3,colour = "red") +
    geom_hline(yintercept = 0,linetype = 2,colour = "grey50") +
    coord_flip() +
    scale_fill_manual(values = Methodcolor) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          axis.text = element_text(size = 12, color = "black"),
          axis.title.x = element_text(size = 14, color = "black"),
          axis.title.y = element_blank()) +
    labs(title = paste("Rare-cell proportion =", prop),
         y = expression(Delta*F1))
  
  return(p)
}


p1 <- plot_deltaF1_box("0.125%")
p1
ggsave(p1, filename = "FigS1a.boxplot.0.125.pdf",width= 5, height= 5, units='in')

p2 <- plot_deltaF1_box("0.25%")
p2
ggsave(p2, filename = "FigS1b.boxplot.0.25.pdf",width= 5, height= 5, units='in')

p3 <- plot_deltaF1_box("0.5%")
p3
ggsave(p3, filename = "FigS1c.boxplot.0.5.pdf",width= 5, height= 5, units='in')

p4 <- plot_deltaF1_box("1%")
p4
ggsave(p4, filename = "FigS1d.boxplot.1.pdf",width= 5, height= 5, units='in')

p5 <- plot_deltaF1_box("2%")
p5
ggsave(p5, filename = "FigS1e.boxplot.2.pdf",width= 5, height= 5, units='in')

p6 <- plot_deltaF1_box("3%")
p6
ggsave(p6, filename = "FigS1f.boxplot.3.pdf",width= 5, height= 5, units='in')

p7 <- plot_deltaF1_box("4%")
p7
ggsave(p7, filename = "FigS1g.boxplot.4.pdf",width= 5, height= 5, units='in')

p8 <- plot_deltaF1_box("5%")
p8
ggsave(p8, filename = "FigS1h.boxplot.5.pdf",width= 5, height= 5, units='in')


## =========================
## Figure S1i-j ----
## =========================
relative_F1_scenario <- relative_F1_scenario %>%
  mutate(Species = ifelse(grepl("^Mouse", SampleID), "Mouse", "Human"))

plot_deltaF1_species <- function(species){
  
  df <- relative_F1_scenario %>% filter(Species == species)
  
  method_order <- df %>%
    group_by(Method) %>%
    dplyr::summarise(
      median_delta = median(deltaF1, na.rm = TRUE),
      mean_delta   = mean(deltaF1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_delta,mean_delta) %>%
    pull(Method)
  
  df$Method <- factor(df$Method, levels = method_order)
  
  ggplot(df, aes(Method, deltaF1, fill = Method)) +
    geom_boxplot(width = 0.7, outlier.size = 0.3) +
    stat_summary(fun = median,geom = "point",shape = 18,size = 3,colour = "red") +
    geom_hline(yintercept = 0,linetype = 2,colour = "grey50") +
    coord_flip() +
    scale_fill_manual(values = Methodcolor) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none",
          axis.text = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_blank()) +
    labs(title = species, y = expression(Delta*F1))
}


p9 <- plot_deltaF1_species("Human")
p9
ggsave(p9, filename = "FigS1i.boxplot.human.pdf",width= 5, height= 5, units='in')

p10 <- plot_deltaF1_species("Mouse")
p10
ggsave(p10, filename = "FigS1j.boxplot.mouse.pdf",width= 5, height= 5, units='in')


## =========================
## Figure S1k ----
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


## Human
rank_human <- relative_F1_scenario %>%
  filter(!startsWith(SampleID, "Mouse")) %>%
  group_by(Method) %>%
  summarise(
    median_delta = median(deltaF1, na.rm = TRUE),
    mean_delta   = mean(deltaF1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_delta), desc(mean_delta)) %>%
  mutate(Rank = row_number()) %>%
  mutate(Group = "Human")


## Mouse
rank_mouse <- relative_F1_scenario %>%
  filter(startsWith(SampleID, "Mouse")) %>%
  group_by(Method) %>%
  summarise(
    median_delta = median(deltaF1, na.rm = TRUE),
    mean_delta   = mean(deltaF1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_delta), desc(mean_delta)) %>%
  mutate(Rank = row_number()) %>%
  mutate(Group = "Mouse")

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
  rank_human %>% select(Method, Rank, Group),
  rank_mouse %>% select(Method, Rank, Group),
  rank_prop
)

rank_matrix <- rank_all %>%
  select(Method, Group, Rank) %>%
  pivot_wider(
    names_from = Group,
    values_from = Rank
  ) %>%
  tibble::column_to_rownames("Method") %>%
  as.matrix()

rank_matrix <- rank_matrix[, c("Overall","Human","Mouse","0.125%","0.25%","0.5%","1%","2%","3%","4%","5%")]
cor_mat <- cor(rank_matrix, method = "spearman")

cor_plot <- cor_mat
cor_plot[upper.tri(cor_plot)] <- NA

pdf("FigS1k.heatmap.corr.pdf", width = 7, height = 6)
Heatmap(
  cor_plot,
  name = "Spearman rho",
  na_col = "white",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "white"),     
  row_names_side = "left",
  column_names_side = "bottom",
  
  col = colorRamp2(
    c(0, 0.5, 1),
    c("#FFFFFF", "#FDB863", "#B2182B")
  ),
  cell_fun = function(j,i,x,y,width,height,fill){
    if(!is.na(cor_plot[i,j])){
      grid.text(
        sprintf("%.2f", cor_plot[i,j]),
        x, y,
        gp = gpar(fontsize = 9)
      )
    }
  }
)
dev.off()

