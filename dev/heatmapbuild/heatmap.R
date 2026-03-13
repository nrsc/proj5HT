library(dplyr)
library(tidyr)
library(ggplot2)

genes <- c("HTR2A", "HTR2B", "HTR2C")
species_keep <- c("M.nemestrina", "M.fascicularis")

df <- brain@meta.data %>%
  dplyr::filter(Species %in% species_keep)

count_df <- df %>%
  count(cluster_Corr_label, Species, name = "N")


ggplot(count_df, aes(
  x = Species,
  y = cluster_Corr_label,
  fill = N
)) +
  geom_tile(color = "grey80") +
  geom_text(aes(label = N), size = 3) +
  scale_fill_viridis_c(option = "C", direction = 1) +
  labs(
    title = "Cell counts by cluster and macaque species",
    x = "Species",
    y = "Cluster (Corr)",
    fill = "N cells"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


library(dplyr)
library(tidyr)
library(ggplot2)

obj <- subset(
  brain,
  subset = subclass_Corr == "L5_ET" & Species %in% species_keep
)

expr <- FetchData(
  obj,
  vars = c(genes, "Species", "cluster_Corr_label")
)

expr_long <- expr %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "expr"
  )

mean_expr <- expr_long %>%
  group_by(cluster_Corr_label, Species, gene) %>%
  summarise(
    mean_expr = mean(expr, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(mean_expr, aes(
  x = Species,
  y = cluster_Corr_label,
  fill = mean_expr
)) +
  geom_tile(color = "grey80") +
  facet_wrap(~gene, nrow = 1) +
  scale_fill_viridis_c(option = "C") +
  labs(
    title = "Mean HTR2 expression by cluster and macaque species",
    fill = "Mean expression"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


obj <- subset(
  brain,
  subset = subclass_Corr == "L5_ET" & Species %in% species_keep
)

expr <- FetchData(
  obj,
  vars = c(genes, "Species", "cluster_Corr_label")
)

expr_long <- expr %>%
  pivot_longer(
    cols = all_of(genes),
    names_to = "gene",
    values_to = "expr"
  )
expr_long <- expr_long %>%
  mutate(expr_pos = expr > 0)
mean_expr <- expr_long %>%
  group_by(cluster_Corr_label, Species, gene) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop")

delta_mean <- mean_expr %>%
  pivot_wider(
    names_from = Species,
    values_from = mean_expr
  ) %>%
  mutate(
    delta = `M.nemestrina` - `M.fascicularis`
  )
ggplot(delta_mean, aes(
  x = gene,
  y = cluster_Corr_label,
  fill = delta
)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0
  ) +
  labs(
    title = "Δ Mean HTR2 expression (Nemestrina − Fascicularis)",
    fill = "Δ mean expr"
  ) +
  theme_bw()

N_df <- expr_long %>%
  group_by(cluster_Corr_label, Species, gene) %>%
  summarise(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = Species, values_from = N)

delta_frac <- delta_frac %>%
  left_join(N_df, by = c("cluster_Corr_label", "gene"))

ggplot(delta_frac, aes(
  x = gene,
  y = cluster_Corr_label,
  fill = delta
)) +
  geom_tile(color = "grey80") +
  # geom_text(
  #   aes(label = paste0("n=", `M.nemestrina`, "/", `M.fascicularis`)),
  #   size = 3
  # ) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    labels = scales::percent_format()
  ) +
  labs(
    title = "Δ HTR2 expression with cell counts (Nemestrina / Fascicularis)",
    fill = "Δ % cells +"
  ) +
  theme_bw()

