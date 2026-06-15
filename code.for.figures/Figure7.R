## =========================
## Figure_7
## =========================

library(dplyr)
library(ggplot2)

load("/home/jinxiuyuan/Proj_scCellFishing/output/benchmarking.result.RData")
load("/home/jinxiuyuan/Proj_scCellFishing/output/HLCA.result.RData")


## =========================
## Figure 7a ----
## =========================
df_overall <- data %>%
  group_by(Method) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    index = "Overall performance",
    rank = rank(-value, ties.method = "min")
  )


df_F1 <- data %>%
  mutate(
    index = cut(
      PrecRare_num,
      breaks = c(0, 0.0025, 0.02, 0.05),
      labels = c(
        "Extremely rare\n(0.125–0.25%)",
        "Moderately rare\n(0.5–2%)",
        "Slightly rare\n(3–5%)"
      ),
      include.lowest = TRUE
    )
  ) %>%
  group_by(Method, index) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(index) %>%
  mutate(
    rank = rank(-value, ties.method = "min")
  ) %>%
  ungroup()


df_similarity <- data %>%
  group_by(Method) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    index = "High-similarity\ntargets",
    rank = rank(-value, ties.method = "min")
  )


df_ann <- anno_data %>%
  group_by(method, ann_level) %>%
  dplyr::summarise(
    value = median(F1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(ann_level) %>%
  mutate(
    rank = rank(-value, ties.method = "min")
  ) %>%
  ungroup() %>%
  rename(
    Method = method,
    index = ann_level
  ) %>%
  mutate(
    index = recode(
      index,
      "ann_level_2" = "Coarse annotation",
      "ann_level_3" = "Intermediate annotation",
      "ann_level_4" = "Fine annotation"
    )
  )


df_time <- data %>%
  group_by(Method) %>%
  dplyr::summarise(
    value = median(Time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    index = "Runtime efficiency",
    rank = rank(value, ties.method = "min")
  )


df_all <- bind_rows(
  df_overall,
  df_F1,
  df_similarity,
  df_ann,
  df_time
)

context_order <- c(
  "Overall performance",
  "Extremely rare\n(0.125–0.25%)",
  "Moderately rare\n(0.5–2%)",
  "Slightly rare\n(3–5%)",
  "High-similarity\ntargets",
  "Coarse annotation",
  "Intermediate annotation",
  "Fine annotation",
  "Runtime efficiency"
)

df_all <- df_all %>%
  mutate(index = factor(index, levels = context_order))


method_order <- df_all %>%
  group_by(Method) %>%
  dplyr::summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_rank) %>%
  pull(Method)

df_all <- df_all %>%
  mutate(Method = factor(Method, levels = rev(method_order)))


p1 <- ggplot(df_all, aes(x = index, y = Method, fill = rank)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = rank), size = 3.2, color = "black") +
  scale_fill_distiller(
    palette = "Spectral",
    direction = 1,
    limits = c(1, 12),
    breaks = c(1, 3, 6, 9, 12),
    name = "Rank") +
  labs(x = "", y = "") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 11),
        panel.grid = element_blank(),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 9, color = "black"))
    

p1


