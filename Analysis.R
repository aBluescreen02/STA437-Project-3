library(tidyverse)
library(ggcorrplot)
library(patchwork)

df <- read.csv('sdss_data.csv')

cat("Objects:", nrow(df), ", Variables:", ncol(df), "\n")
table(df$class)

# Photometric features for multivariate analysis
photo_cols <- c('u', 'g', 'r', 'i', 'z')
X <- as.matrix(df[, photo_cols])
y <- df$class

class_colors <- c(GALAXY = "#378ADD", STAR = "#1D9E75", QSO = "#D85A30")

band_labels <- c(
  u = "u band (ultraviolet, ~354 nm)",
  g = "g band (green, ~477 nm)",
  r = "r band (red, ~623 nm)",
  i = "i band (near-infrared, ~763 nm)",
  z = "z band (infrared, ~913 nm)"
)

n_labels <- c(
  GALAXY = "Galaxy (n=4,998)",
  STAR   = "Star (n=4,152)",
  QSO    = "QSO (n=850)"
)

df_long <- df %>%
  select(all_of(c(photo_cols, "class"))) %>%
  filter(if_all(all_of(photo_cols), ~ . < 35)) %>%   # drop extreme outliers
  pivot_longer(cols = all_of(photo_cols),
               names_to  = "band",
               values_to = "magnitude") %>%
  mutate(
    band  = factor(band,  levels = photo_cols,          labels = band_labels),
    class = factor(class, levels = c("GALAXY", "STAR", "QSO"), labels = n_labels)
  )

ggplot(df_long, aes(x = magnitude, fill = class, color = class)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.45, linewidth = 0.3) +
  scale_fill_manual(values  = setNames(class_colors, n_labels)) +
  scale_color_manual(values = setNames(class_colors, n_labels)) +
  facet_wrap(~ band, ncol = 2, scales = "free") +
  labs(
    x     = "Magnitude (AB system)",
    y     = "Count",
    fill  = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "top",
    legend.key.size  = unit(0.5, "cm"),
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey92"),
    axis.text        = element_text(size = 8, colour = "grey50"),
    axis.title       = element_text(size = 9, colour = "grey40"),
    legend.text      = element_text(size = 9),
    plot.margin      = margin(10, 10, 10, 10)
  )

# Helper: correlation matrix + ggcorrplot
make_corrplot <- function(data, title, n) {
  R <- cor(data[, photo_cols])
  ggcorrplot(R,
             method   = "square",
             type     = "full",
             lab      = TRUE,
             lab_size = 3,
             colors   = c("#2166AC", "#F7F7F7", "#B2182B"),
             title    = paste0(title, " (n = ", scales::comma(n), ")"),
             ggtheme  = theme_minimal(base_size = 10)) +
    theme(plot.title       = element_text(size = 10, hjust = 0.5),
          legend.position  = "none",
          axis.text        = element_text(size = 9))
}

p_all  <- make_corrplot(df,                          "All objects", nrow(df))
p_gal  <- make_corrplot(filter(df, class=="GALAXY"), "Galaxy",      sum(df$class=="GALAXY"))
p_star <- make_corrplot(filter(df, class=="STAR"),   "Star",        sum(df$class=="STAR"))
p_qso  <- make_corrplot(filter(df, class=="QSO"),    "QSO",         sum(df$class=="QSO"))

# Shared legend from the full-data plot (add it back just for extraction)
p_legend <- make_corrplot(df, "", nrow(df)) +
  theme(legend.position = "right")

layout <- "
AB
CD
"
p_combined <- (p_all + plot_spacer()) / (p_gal + p_star + p_qso) +
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom")

# Simpler layout: all four panels in a 2x2 grid
(p_all | plot_spacer()) / (p_gal | p_star | p_qso) +
  plot_annotation(
    title   = "Photometric band correlations — overall vs. by class",
    caption = "Simpson's paradox check: compare overall correlations to within-class correlations",
    theme   = theme(plot.title   = element_text(size = 12, face = "plain"),
                    plot.caption = element_text(size = 8, colour = "grey50"))
  )

df_adj <- df %>%
  mutate(
    `u-g` = u - g,
    `g-r` = g - r,
    `r-i` = r - i,
    `i-z` = i - z
  )

ci_cols <- c('u-g', 'g-r', 'r-i', 'i-z')

# ── Histograms ────────────────────────────────────────────────────────────────
df_adj_hist <- df_adj %>%
  select(all_of(c(ci_cols, "class"))) %>%
  pivot_longer(cols = all_of(ci_cols), names_to = "index", values_to = "value") %>%
  mutate(index = factor(index, levels = ci_cols))

p_hist <- ggplot(df_adj_hist, aes(x = value, fill = class, color = class)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.42, linewidth = 0.3) +
  scale_fill_manual(values  = class_colors) +
  scale_color_manual(values = class_colors) +
  facet_wrap(~ index, ncol = 2, scales = "free") +
  labs(x = "Color index (mag)", y = "Count", fill = NULL, color = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "top",
        strip.text       = element_text(size = 9),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(size = 8, colour = "grey50"))

# ── Correlation heatmaps ──────────────────────────────────────────────────────
make_corrplot <- function(data, title) {
  R <- cor(data[, ci_cols])
  ggcorrplot(R,
             method   = "square",
             type     = "full",
             lab      = TRUE,
             lab_size = 3,
             colors   = c("#2166AC", "#F7F7F7", "#B2182B"),
             title    = title,
             ggtheme  = theme_minimal(base_size = 9)) +
    theme(plot.title      = element_text(size = 9, hjust = 0.5),
          legend.position = "none",
          axis.text       = element_text(size = 8))
}

p_all  <- make_corrplot(df_adj,                            "All objects")
p_gal  <- make_corrplot(filter(df_adj, class == "GALAXY"), "Galaxy")
p_star <- make_corrplot(filter(df_adj, class == "STAR"),   "Star")
p_qso  <- make_corrplot(filter(df_adj, class == "QSO"),    "QSO")

p_corr <- (p_all | plot_spacer()) / (p_gal | p_star | p_qso) +
  plot_annotation(title = "Color index correlations — overall vs. by class",
                  theme = theme(plot.title = element_text(size = 11)))

# ── Combine ───────────────────────────────────────────────────────────────────
p_hist / p_corr + plot_layout(heights = c(1, 1.2))

