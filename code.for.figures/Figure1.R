## =========================
## Figure_1
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scmamp)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.RData")

Methodcolor <- c("aKNNO" = "#9E2A2B", "CellSIUS" = "#E8768A","CIARA" = "#C96E2D", "EDGE" = "#F39C45", "GiniClust3" = "#EFCB68","RaceID2" = "#3A9D54",
                 "RaceID3" = "#A8D8A2", "RareQ" = "#7B43A6", "SCA" = "#C8A6E3","scCAD_1%"  = "#1C4E73", "scCAD_5%" = "#3E7FBF", "SCISSORS"  = "#8EC9E6")


## =========================
## Friedman test ----
## =========================
benchmark_global <- data %>%
  group_by(Method, Scenario) %>%
  dplyr::summarise(F1      = median(F1, na.rm = TRUE),
                   MCC     = median(MCC, na.rm = TRUE),
                   .groups = "drop")

benchmark_test <- benchmark_global %>%
  mutate(F1 = ifelse(is.na(F1), 0, F1)) %>%
  group_by(Scenario) %>%
  filter(n_distinct(Method)==12) %>%
  ungroup()

friedman.test(F1 ~ Method | Scenario, data=benchmark_test)


## =========================
## Figure 1b ----
## =========================
rank_df <- benchmark_test %>%
  group_by(Scenario) %>%
  mutate(Rank = rank(-F1, ties.method = "average")) %>%
  ungroup()

method_order_df <- rank_df %>%
  group_by(Method) %>%
  dplyr::summarise(
    median_rank = median(Rank, na.rm = TRUE),
    mean_rank   = mean(Rank, na.rm = TRUE),
    q1_rank     = quantile(Rank, 0.25, na.rm = TRUE),
    q3_rank     = quantile(Rank, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    median_rank,
    q3_rank,
    mean_rank,
    as.character(Method)
  )

method_order <- method_order_df %>%
  pull(Method)


rank_df$Method <- factor(rank_df$Method, levels = method_order)

p1 <- ggplot(rank_df, aes(Method, Rank, fill = Method)) +
  geom_boxplot(width = 0.7, alpha = 0.9, outlier.size = 0.6) +
  scale_fill_manual(values = Methodcolor) +
  scale_y_reverse(breaks = 1:12) +
  labs(y = "Rank", x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y =element_text(size = 14, color = "black"))
p1


## =========================
## Figure 1d ----
## =========================
all_methods <- unique(benchmark_test$Method)
top1 <- benchmark_test %>%
  group_by(Scenario) %>%
  mutate(rk = rank(-F1, ties.method = "min")) %>%
  filter(rk == 1) %>%
  ungroup() %>%
  dplyr::count(Method) %>%
  complete(Method = all_methods, fill = list(n = 0))

p2 <- ggplot(top1, aes(x = reorder(Method, -n), y = n, fill = Method)) +
  geom_col(color = "black", alpha = 0.9, linewidth = 0.2) +
  scale_fill_manual(values = Methodcolor) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  labs(y = "Top-rank occurrences", x= "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y =element_text(size = 14, color = "black"))
p2


## =========================
## Figure S1a ----
## =========================
perf_matrix <- benchmark_test %>%
  dplyr::select(Scenario, Method, F1) %>%
  tidyr::pivot_wider(names_from = Method, values_from = F1)

perf_matrix <- as.data.frame(perf_matrix)
rownames(perf_matrix) <- perf_matrix$Scenario
perf_matrix$Scenario <- NULL

plotCD(results.matrix = perf_matrix, alpha = 0.05)


## =========================
## Figure 1c ----
## =========================
perf_wide <- benchmark_test %>%
  dplyr::select(Scenario, Method, F1) %>%
  tidyr::pivot_wider(names_from = Method, values_from = F1) %>%
  as.data.frame()
rownames(perf_wide) <- perf_wide$Scenario
perf_wide$Scenario  <- NULL

perf_wide <- perf_wide[, method_order]
methods   <- colnames(perf_wide)

wilcox_pairwise <- combn(methods, 2, simplify = FALSE) %>%
  purrr::map_dfr(function(x) {
    
    m1 <- x[1]
    m2 <- x[2]
    
    tmp <- perf_wide[, c(m1, m2)]
    tmp <- tmp[complete.cases(tmp), ]
    
    test_res <- wilcox.test(
      tmp[[m1]],
      tmp[[m2]],
      paired = TRUE,
      exact = FALSE
    )
    
    data.frame(
      Method1 = m1,
      Method2 = m2,
      n_pair = nrow(tmp),
      median_Method1 = median(tmp[[m1]], na.rm = TRUE),
      median_Method2 = median(tmp[[m2]], na.rm = TRUE),
      median_diff = median(tmp[[m1]] - tmp[[m2]], na.rm = TRUE),
      statistic = as.numeric(test_res$statistic),
      p_value = test_res$p.value
    )
  })


wilcox_pairwise <- wilcox_pairwise %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
         significance = case_when(
           p_adj_BH < 0.001 ~ "***",
           p_adj_BH < 0.01  ~ "**",
           p_adj_BH < 0.05  ~ "*",
           TRUE ~ "ns")) %>%
  arrange(p_adj_BH)


p_mat <- matrix(
  NA,
  nrow = length(methods),
  ncol = length(methods),
  dimnames = list(methods, methods)
)

diff_mat <- matrix(
  NA,
  nrow = length(methods),
  ncol = length(methods),
  dimnames = list(methods, methods)
)

for (i in seq_len(nrow(wilcox_pairwise))) {
  
  m1 <- wilcox_pairwise$Method1[i]
  m2 <- wilcox_pairwise$Method2[i]
  
  p_mat[m1, m2] <- wilcox_pairwise$p_adj_BH[i]
  p_mat[m2, m1] <- wilcox_pairwise$p_adj_BH[i]
  
  diff_mat[m1, m2] <- wilcox_pairwise$median_diff[i]
  diff_mat[m2, m1] <- -wilcox_pairwise$median_diff[i]
}

diag(p_mat) <- 1
diag(diff_mat) <- 0

wilcox_plot_df <- wilcox_pairwise %>%
  mutate(Method1 = factor(Method1, levels = method_order),
         Method2 = factor(Method2, levels = method_order),
         neg_log10_padj = -log10(p_adj_BH))

p3 <- ggplot(wilcox_plot_df, aes(x = Method1, y = Method2, fill = neg_log10_padj)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = significance), size = 4) +
  scale_fill_gradient(low = "white", high = "firebrick", name = "-log10(adj.P.Val)") +
  labs(x = "", y = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
p3    


## =========================
## Figure S1c ----
## =========================
benchmark_test <- benchmark_global %>%
  mutate(MCC = ifelse(is.na(MCC), 0, MCC)) %>%
  group_by(Scenario) %>%
  filter(n_distinct(Method)==12) %>%
  ungroup()

friedman.test(MCC ~ Method | Scenario, data=benchmark_test)


rank_df <- benchmark_test %>%
  group_by(Scenario) %>%
  mutate(Rank = rank(-MCC, ties.method = "average")) %>%
  ungroup()

method_order_df <- rank_df %>%
  group_by(Method) %>%
  dplyr::summarise(
    median_rank = median(Rank, na.rm = TRUE),
    mean_rank   = mean(Rank, na.rm = TRUE),
    q1_rank     = quantile(Rank, 0.25, na.rm = TRUE),
    q3_rank     = quantile(Rank, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    median_rank,
    q3_rank,
    mean_rank,
    as.character(Method)
  )

method_order <- method_order_df %>%
  pull(Method)


rank_df$Method <- factor(rank_df$Method, levels = method_order)
p4 <- ggplot(rank_df, aes(Method, Rank, fill = Method)) +
  geom_boxplot(width = 0.7, alpha = 0.9, outlier.size = 0.6) +
  scale_fill_manual(values = Methodcolor) +
  scale_y_reverse(breaks = 1:12) +
  labs(y = "Rank", x = "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y =element_text(size = 14, color = "black"))
p4


## =========================
## Figure S1d ----
## =========================
all_methods <- unique(benchmark_test$Method)
top1 <- benchmark_test %>%
  group_by(Scenario) %>%
  mutate(rk = rank(-MCC, ties.method = "min")) %>%
  filter(rk == 1) %>%
  ungroup() %>%
  dplyr::count(Method)%>%
  complete(Method = all_methods, fill = list(n = 0))

p5 <- ggplot(top1, aes(x = reorder(Method, -n), y = n, fill = Method)) +
  geom_col(color = "black", alpha = 0.9, linewidth = 0.2) +
  scale_fill_manual(values = Methodcolor) +
  geom_text(aes(label = n), vjust = -0.3, size = 4) +
  labs(y = "Top-rank occurrences", x= "") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y =element_text(size = 14, color = "black"))
p5


## =========================
## Figure S1b ----
## =========================
perf_matrix <- benchmark_test %>%
  dplyr::select(Scenario, Method, MCC) %>%
  tidyr::pivot_wider(names_from = Method, values_from = MCC)

perf_matrix <- as.data.frame(perf_matrix)
rownames(perf_matrix) <- perf_matrix$Scenario
perf_matrix$Scenario <- NULL

plotCD(results.matrix = perf_matrix, alpha = 0.05)

