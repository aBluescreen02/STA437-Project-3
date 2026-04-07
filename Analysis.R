library(MASS)
library(ggcorrplot)
library(patchwork)
library(tidyverse)
library(robustbase)
library(MVN)
library(psych)
library(cowplot)
library(nnet)

# ── Data loading ──────────────────────────────────────────────────────────────
df <- read.csv('sdss_data.csv')

cat("Objects:", nrow(df), ", Variables:", ncol(df), "\n")
table(df$class)

photo_cols   <- c('u', 'g', 'r', 'i', 'z')
ci_cols      <- c('u-g', 'g-r', 'r-i', 'i-z')
class_order  <- c("GALAXY", "STAR", "QSO")
class_colors <- c(GALAXY = "#378ADD", STAR = "#1D9E75", QSO = "#D85A30")

n_labels <- c(
  GALAXY = "Galaxy (n=4,998)",
  STAR   = "Star (n=4,152)",
  QSO    = "QSO (n=850)"
)

X_bands <- as.matrix(df[, photo_cols])
y        <- df$class

# ── Shared helper: four-quadrant D-D plot ─────────────────────────────────────
make_dd_plot <- function(dd_df, cutoff, p, feature_label) {
  dd_df <- dd_df %>%
    mutate(
      outlier_type = case_when(
        robust > cutoff & classical > cutoff  ~ "True outlier",
        robust > cutoff & classical <= cutoff ~ "Masked outlier",
        robust <= cutoff & classical > cutoff ~ "Swamping",
        TRUE                                  ~ "Regular"
      ),
      outlier_type = factor(outlier_type,
                            levels = c("Regular", "Swamping",
                                       "Masked outlier", "True outlier"))
    )
  
  type_colors <- c(
    "Regular"        = "grey70",
    "Swamping"       = "#7FBCD2",
    "Masked outlier" = "#9B72CF",
    "True outlier"   = "#D85A30"
  )
  type_shapes <- c(
    "Regular"        = 16,
    "Swamping"       = 16,
    "Masked outlier" = 18,
    "True outlier"   = 18
  )
  
  counts        <- dd_df %>% count(outlier_type)
  legend_labels <- setNames(
    paste0(counts$outlier_type, " (n=", counts$n, ")"),
    counts$outlier_type
  )
  
  x_max <- max(quantile(dd_df$classical, 0.999) * 1.05, cutoff * 2.5)
  y_max <- max(quantile(dd_df$robust,    0.999) * 1.05, cutoff * 2.5)
  
  ggplot(dd_df, aes(x = classical, y = robust,
                    color = outlier_type, shape = outlier_type)) +
    annotate("rect", xmin = cutoff, xmax = Inf,  ymin = cutoff, ymax = Inf,
             fill = "#D85A30", alpha = 0.04) +
    annotate("rect", xmin = -Inf,  xmax = cutoff, ymin = cutoff, ymax = Inf,
             fill = "#9B72CF", alpha = 0.04) +
    annotate("rect", xmin = cutoff, xmax = Inf,  ymin = -Inf,  ymax = cutoff,
             fill = "#7FBCD2", alpha = 0.04) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = cutoff,
               linetype = "dotted", linewidth = 0.4, color = "grey40") +
    geom_hline(yintercept = cutoff,
               linetype = "dotted", linewidth = 0.4, color = "grey40") +
    geom_point(data = filter(dd_df, outlier_type == "Regular"),
               alpha = 0.25, size = 0.9) +
    geom_point(data = filter(dd_df, outlier_type != "Regular"),
               alpha = 0.75, size = 1.6) +
    annotate("text", x = x_max * 0.97, y = y_max * 0.97,
             label = "True outlier", hjust = 1, vjust = 1, size = 2.8,
             color = "#D85A30", fontface = "bold") +
    annotate("text", x = cutoff * 0.95, y = y_max * 0.97,
             label = "Masked outlier\n(classical misses)",
             hjust = 1, vjust = 1, size = 2.8,
             color = "#9B72CF", fontface = "bold") +
    annotate("text", x = x_max * 0.97, y = cutoff * 0.95,
             label = "Swamping\n(false alarm)", hjust = 1, vjust = 1, size = 2.8,
             color = "#7FBCD2", fontface = "bold") +
    annotate("text", x = cutoff * 0.95, y = cutoff * 0.95,
             label = "Regular", hjust = 1, vjust = 1, size = 2.8,
             color = "grey55", fontface = "plain") +
    scale_color_manual(values = type_colors, labels = legend_labels) +
    scale_shape_manual(values = type_shapes, labels = legend_labels) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    labs(
      x       = "Classical Mahalanobis distance",
      y       = "Robust (MCD) Mahalanobis distance",
      color   = NULL, shape = NULL,
      title   = feature_label
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "top",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7.5, colour = "grey50"),
      plot.title       = element_text(size = 11, hjust = 0.5)
    )
}

# ── Pooled D-D plot: points coloured by class ─────────────────────────────────
make_dd_plot_pooled <- function(dd_df, cutoff, p, feature_label) {
  x_max <- max(quantile(dd_df$classical, 0.999) * 1.05, cutoff * 2.5)
  y_max <- max(quantile(dd_df$robust,    0.999) * 1.05, cutoff * 2.5)
  
  ggplot(dd_df, aes(x = classical, y = robust, color = class)) +
    annotate("rect", xmin = cutoff, xmax = Inf,  ymin = cutoff, ymax = Inf,
             fill = "#D85A30", alpha = 0.04) +
    annotate("rect", xmin = -Inf,  xmax = cutoff, ymin = cutoff, ymax = Inf,
             fill = "#9B72CF", alpha = 0.04) +
    annotate("rect", xmin = cutoff, xmax = Inf,  ymin = -Inf,  ymax = cutoff,
             fill = "#7FBCD2", alpha = 0.04) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = cutoff,
               linetype = "dotted", linewidth = 0.4, color = "grey40") +
    geom_hline(yintercept = cutoff,
               linetype = "dotted", linewidth = 0.4, color = "grey40") +
    geom_point(alpha = 0.30, size = 0.9) +
    annotate("text", x = x_max * 0.97, y = y_max * 0.97,
             label = "True outlier", hjust = 1, vjust = 1, size = 2.8,
             color = "grey30", fontface = "bold") +
    annotate("text", x = cutoff * 0.95, y = y_max * 0.97,
             label = "Masked outlier\n(classical misses)",
             hjust = 1, vjust = 1, size = 2.8,
             color = "grey30", fontface = "bold") +
    annotate("text", x = x_max * 0.97, y = cutoff * 0.95,
             label = "Swamping\n(false alarm)", hjust = 1, vjust = 1, size = 2.8,
             color = "grey30", fontface = "bold") +
    annotate("text", x = cutoff * 0.95, y = cutoff * 0.95,
             label = "Regular", hjust = 1, vjust = 1, size = 2.8,
             color = "grey55", fontface = "plain") +
    scale_color_manual(values = class_colors, name = "Class") +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    labs(
      x     = "Classical Mahalanobis distance",
      y     = "Robust (MCD) Mahalanobis distance",
      title = feature_label
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "top",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7.5, colour = "grey50"),
      plot.title       = element_text(size = 11, hjust = 0.5)
    )
}

get_outlier_flags <- function(dd_df, cutoff) {
  dd_df %>%
    mutate(
      outlier_type = case_when(
        robust > cutoff & classical > cutoff  ~ "True outlier",
        robust > cutoff & classical <= cutoff ~ "Masked outlier",
        robust <= cutoff & classical > cutoff ~ "Swamping",
        TRUE                                  ~ "Regular"
      ),
      outlier = outlier_type %in% c("True outlier", "Masked outlier")
    )
}

summarise_outliers <- function(dd_typed, label) {
  dd_typed %>%
    group_by(class, outlier_type) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = outlier_type, values_from = n, values_fill = 0) %>%
    mutate(
      total    = rowSums(across(where(is.numeric))),
      outliers = `True outlier` + `Masked outlier`,
      pct      = round(100 * outliers / total, 1)
    ) %>%
    mutate(feature_set = label)
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — RAW PHOTOMETRIC BANDS (u, g, r, i, z) #####
# ═══════════════════════════════════════════════════════════════════════════════

band_labels <- c(
  u = "u band (~354 nm)", g = "g band (~477 nm)", r = "r band (~623 nm)",
  i = "i band (~763 nm)", z = "z band (~913 nm)"
)

df_bands_long <- df %>%
  select(all_of(c(photo_cols, "class"))) %>%
  pivot_longer(cols = all_of(photo_cols), names_to = "band", values_to = "magnitude") %>%
  mutate(band  = factor(band,  levels = photo_cols, labels = band_labels),
         class = factor(class, levels = class_order, labels = n_labels))

p_bands_hist <- ggplot(df_bands_long, aes(x = magnitude, fill = class, color = class)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.45, linewidth = 0.3) +
  scale_fill_manual(values  = setNames(class_colors, n_labels)) +
  scale_color_manual(values = setNames(class_colors, n_labels)) +
  facet_wrap(~ band, ncol = 1, scales = "free") +
  labs(x = "Magnitude (AB system)", y = "Count", fill = NULL, color = NULL,
       title = "Distributions") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none", strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92"),
        axis.text  = element_text(size = 7, colour = "grey50"),
        axis.title = element_text(size = 8, colour = "grey40"),
        plot.title = element_text(size = 9, face = "plain"),
        plot.margin = margin(4, 4, 4, 4))

make_corrplot_bands <- function(data, title, n = NULL) {
  R     <- cor(data[, photo_cols])
  label <- if (!is.null(n)) paste0(title, "\n(n=", scales::comma(n), ")") else title
  ggcorrplot(R, method = "square", type = "full", lab = TRUE, lab_size = 2.5,
             colors = c("#2166AC", "#F7F7F7", "#B2182B"), title = label,
             ggtheme = theme_minimal(base_size = 9)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          legend.position = "none", axis.text = element_text(size = 7),
          plot.margin = margin(2, 2, 2, 2))
}

p_bands_all  <- make_corrplot_bands(df,                            "All",    nrow(df))
p_bands_gal  <- make_corrplot_bands(filter(df, class == "GALAXY"), "Galaxy", sum(df$class == "GALAXY"))
p_bands_star <- make_corrplot_bands(filter(df, class == "STAR"),   "Star",   sum(df$class == "STAR"))
p_bands_qso  <- make_corrplot_bands(filter(df, class == "QSO"),    "QSO",    sum(df$class == "QSO"))

p_bands_corr <- (p_bands_all | p_bands_gal) / (p_bands_star | p_bands_qso) +
  plot_annotation(title = "B",
                  theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0)))

p_class_legend <- ggplot(
  tibble(class = factor(n_labels, levels = n_labels), x = 1, y = seq_along(n_labels)),
  aes(x = x, y = y, fill = class, color = class)
) +
  geom_point(shape = 21, size = 4, stroke = 0.8) +
  scale_fill_manual(values  = setNames(class_colors, n_labels), labels = n_labels) +
  scale_color_manual(values = setNames(class_colors, n_labels), labels = n_labels) +
  guides(fill  = guide_legend(override.aes = list(size = 4, alpha = 0.6),
                              ncol = 1, title = "Class"), color = "none") +
  labs(fill = NULL) + theme_void() +
  theme(legend.position = "right", legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8), legend.key.size = unit(0.45, "cm"),
        legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0))

p_corr_legend <- ggplot(data.frame(x = 1, y = 1, z = 0), aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation",
                       guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                              barwidth = unit(0.4, "cm"), barheight = unit(2.5, "cm"),
                                              ticks = TRUE, label = TRUE)) +
  theme_void() +
  theme(legend.position = "right", legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8), legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(0,0,0,0))

p_leg_class <- get_legend(p_class_legend)
p_leg_corr  <- get_legend(p_corr_legend)

p_legend_combined <- plot_grid(p_leg_class, p_leg_corr, ncol = 1,
                               rel_heights = c(0.55, 0.45), align = "v")

p_bands_hist_labelled <- ggdraw() +
  draw_plot(p_bands_hist) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 11)

p_bands_combined <- plot_grid(p_legend_combined, p_bands_hist_labelled, p_bands_corr,
                              nrow = 1, rel_widths = c(0.2, 0.42, 0.58),
                              align = "h", axis = "tb")
plot_grid(p_bands_combined, ncol = 1)  # Figure: raw dists

# ── D-D plots: raw bands ──────────────────────────────────────────────────────
p_bands      <- ncol(X_bands)
cutoff_bands <- sqrt(qchisq(0.975, df = p_bands))

mcd_bands_pooled <- covMcd(X_bands, alpha = 0.95)
dd_bands_pooled  <- tibble(
  classical = sqrt(mahalanobis(X_bands, center = colMeans(X_bands), cov = cov(X_bands))),
  robust    = sqrt(mahalanobis(X_bands, center = mcd_bands_pooled$center,
                               cov = mcd_bands_pooled$cov)),
  class     = df$class
)
dd_bands_pooled_typed <- get_outlier_flags(dd_bands_pooled, cutoff_bands)
cat(sprintf("[Bands pooled] Outliers: %d / %d (%.1f%%)\n",
            sum(dd_bands_pooled_typed$outlier), nrow(dd_bands_pooled_typed),
            100 * mean(dd_bands_pooled_typed$outlier)))

dd_bands_typed <- map_dfr(class_order, function(cls) {
  idx  <- which(df$class == cls)
  Xsub <- X_bands[idx, ]
  d_cl <- sqrt(mahalanobis(Xsub, center = colMeans(Xsub), cov = cov(Xsub)))
  mcd  <- covMcd(Xsub, alpha = 0.95)
  d_ro <- sqrt(mahalanobis(Xsub, center = mcd$center, cov = mcd$cov))
  tibble(classical = d_cl, robust = d_ro, class = cls, row_idx = idx)
}) %>% get_outlier_flags(cutoff_bands)

cat(sprintf("[Bands] Within-class outliers: %d / %d (%.1f%%)\n",
            sum(dd_bands_typed$outlier), nrow(dd_bands_typed),
            100 * mean(dd_bands_typed$outlier)))

# Combined D-D figure: pooled (top, spanning full width) + within-class (bottom)
p_dd_bands_pooled <- make_dd_plot_pooled(dd_bands_pooled, cutoff_bands, p_bands,
                                         "Pooled")
p_dd_bands_within <- map(class_order, function(cls) {
  make_dd_plot(filter(dd_bands_typed, class == cls), cutoff_bands, p_bands,
               paste0(cls, " (n=", sum(df$class == cls), ")"))
})

p_dd_bands_pooled_labelled <- ggdraw() +
  draw_plot(p_dd_bands_pooled) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_dd_bands_within_b <- ggdraw() +
  draw_plot(p_dd_bands_within[[1]]) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

blank <- ggplot() + theme_void()

plot_grid(                                                           # Figure: dd (raw bands)
  plot_grid(blank, p_dd_bands_pooled_labelled, blank,
            nrow = 1, rel_widths = c(1, 1, 1)),
  plot_grid(p_dd_bands_within_b, p_dd_bands_within[[2]], p_dd_bands_within[[3]],
            nrow = 1),
  ncol = 1
)


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — ADJACENT COLOR INDICES (u-g, g-r, r-i, i-z) #####
# ═══════════════════════════════════════════════════════════════════════════════

df_ci <- df %>%
  mutate(`u-g` = u - g, `g-r` = g - r, `r-i` = r - i, `i-z` = i - z)
X_ci <- as.matrix(df_ci[, ci_cols])

df_ci_long <- df_ci %>%
  select(all_of(c(ci_cols, "class"))) %>%
  pivot_longer(cols = all_of(ci_cols), names_to = "index", values_to = "value") %>%
  mutate(index = factor(index, levels = ci_cols))

p_ci_hist <- ggplot(df_ci_long, aes(x = value, fill = class, color = class)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.42, linewidth = 0.3) +
  scale_fill_manual(values = class_colors) +
  scale_color_manual(values = class_colors) +
  facet_wrap(~ index, ncol = 1, scales = "free") +
  labs(x = "Color index (mag)", y = "Count", fill = NULL, color = NULL, title = "A") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "none", strip.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        axis.text  = element_text(size = 7, colour = "grey50"),
        axis.title = element_text(size = 8, colour = "grey40"),
        plot.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.margin = margin(4, 4, 4, 4))

make_corrplot_ci <- function(data, title) {
  R <- cor(data[, ci_cols])
  ggcorrplot(R, method = "square", type = "full", lab = TRUE, lab_size = 2.5,
             colors = c("#2166AC", "#F7F7F7", "#B2182B"), title = title,
             ggtheme = theme_minimal(base_size = 9)) +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          legend.position = "none", axis.text = element_text(size = 7),
          plot.margin = margin(2, 2, 2, 2))
}

p_ci_all  <- make_corrplot_ci(df_ci,                            "All\n(n=10,000)")
p_ci_gal  <- make_corrplot_ci(filter(df_ci, class == "GALAXY"), "Galaxy\n(n=4,998)")
p_ci_star <- make_corrplot_ci(filter(df_ci, class == "STAR"),   "Star\n(n=4,152)")
p_ci_qso  <- make_corrplot_ci(filter(df_ci, class == "QSO"),    "QSO\n(n=850)")

p_ci_corr <- (p_ci_all | p_ci_gal) / (p_ci_star | p_ci_qso) +
  plot_annotation(title = "B",
                  theme = theme(plot.title = element_text(size = 11, face = "bold", hjust = 0)))

plot_grid(plot_grid(p_legend_combined, p_ci_hist, p_ci_corr, nrow = 1,
                    rel_widths = c(0.2, 0.42, 0.58), align = "h", axis = "tb"), ncol = 1)  # Figure: adj dists


# ── D-D plots: color indices ──────────────────────────────────────────────────
p_ci      <- ncol(X_ci)
cutoff_ci <- sqrt(qchisq(0.975, df = p_ci))

mcd_ci_pooled <- covMcd(X_ci, alpha = 0.95)
dd_ci_pooled  <- tibble(
  classical = sqrt(mahalanobis(X_ci, center = colMeans(X_ci), cov = cov(X_ci))),
  robust    = sqrt(mahalanobis(X_ci, center = mcd_ci_pooled$center,
                               cov = mcd_ci_pooled$cov)),
  class     = df_ci$class
)
dd_ci_pooled_typed <- get_outlier_flags(dd_ci_pooled, cutoff_ci)
cat(sprintf("[CI pooled] Outliers: %d / %d (%.1f%%)\n",
            sum(dd_ci_pooled_typed$outlier), nrow(dd_ci_pooled_typed),
            100 * mean(dd_ci_pooled_typed$outlier)))

dd_ci_typed <- map_dfr(class_order, function(cls) {
  idx  <- which(df_ci$class == cls)
  Xsub <- X_ci[idx, ]
  d_cl <- sqrt(mahalanobis(Xsub, center = colMeans(Xsub), cov = cov(Xsub)))
  mcd  <- covMcd(Xsub, alpha = 0.95)
  d_ro <- sqrt(mahalanobis(Xsub, center = mcd$center, cov = mcd$cov))
  tibble(classical = d_cl, robust = d_ro, class = cls, row_idx = idx)
}) %>% get_outlier_flags(cutoff_ci)

cat(sprintf("[CI] Within-class outliers: %d / %d (%.1f%%)\n",
            sum(dd_ci_typed$outlier), nrow(dd_ci_typed),
            100 * mean(dd_ci_typed$outlier)))

# Combined D-D figure: pooled (top, centred) + within-class (bottom)
p_dd_ci_pooled <- make_dd_plot_pooled(dd_ci_pooled, cutoff_ci, p_ci, "Pooled")
p_dd_ci_within <- map(class_order, function(cls) {
  make_dd_plot(filter(dd_ci_typed, class == cls), cutoff_ci, p_ci, cls)
})

p_dd_ci_pooled_labelled <- ggdraw() +
  draw_plot(p_dd_ci_pooled) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_dd_ci_within_b <- ggdraw() +
  draw_plot(p_dd_ci_within[[1]]) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

plot_grid(                                                           # Figure: dd (color indices)
  plot_grid(blank, p_dd_ci_pooled_labelled, blank,
            nrow = 1, rel_widths = c(1, 1, 1)),
  plot_grid(p_dd_ci_within_b, p_dd_ci_within[[2]], p_dd_ci_within[[3]],
            nrow = 1),
  ncol = 1
)


# ── Outlier summaries ─────────────────────────────────────────────────────────
cat("\n--- Pooled outlier summary ---\n")
bind_rows(summarise_outliers(dd_bands_pooled_typed, "Raw bands (pooled)"),
          summarise_outliers(dd_ci_pooled_typed,    "Color indices (pooled)")) %>%
  select(feature_set, class,
         any_of(c("True outlier", "Masked outlier", "Swamping", "Regular")),
         outliers, total, pct) %>%
  print(n = Inf)

cat("\n--- Within-class outlier summary ---\n")
bind_rows(summarise_outliers(dd_bands_typed, "Raw bands (within-class)"),
          summarise_outliers(dd_ci_typed,    "Color indices (within-class)")) %>%
  select(feature_set, class,
         any_of(c("True outlier", "Masked outlier", "Swamping", "Regular")),
         outliers, total, pct) %>%
  print(n = Inf)

# ── Carry forward datasets ────────────────────────────────────────────────────
make_clean_df <- function(base_df, flags_df) {
  base_df %>%
    mutate(row_idx = row_number()) %>%
    left_join(flags_df %>%
                select(row_idx, outlier_type, outlier,
                       d_classical = classical, d_robust = robust) %>%
                arrange(row_idx), by = "row_idx") %>%
    select(-row_idx) %>%
    filter(!outlier)
}

# 2SD / 95% MCD — main analysis datasets
df_bands_no_out <- make_clean_df(df,    dd_bands_typed)
df_ci_no_out    <- make_clean_df(df_ci, dd_ci_typed)

# 3SD / 99.7% MCD — appendix sensitivity only (not used in downstream sections)
cutoff_bands_3sd <- sqrt(qchisq(0.9985, df = p_bands))
cutoff_ci_3sd    <- sqrt(qchisq(0.9985, df = p_ci))

make_robust_flags <- function(X_mat, class_vec, alpha_mcd, cutoff) {
  map_dfr(class_order, function(cls) {
    idx  <- which(class_vec == cls)
    Xsub <- X_mat[idx, ]
    d_cl <- sqrt(mahalanobis(Xsub, center = colMeans(Xsub), cov = cov(Xsub)))
    mcd  <- covMcd(Xsub, alpha = alpha_mcd)
    d_ro <- sqrt(mahalanobis(Xsub, center = mcd$center, cov = mcd$cov))
    tibble(classical = d_cl, robust = d_ro, class = cls, row_idx = idx)
  }) %>% get_outlier_flags(cutoff)
}

dd_bands_3sd <- make_robust_flags(X_bands, df$class,    0.997,     cutoff_bands_3sd)
dd_ci_3sd    <- make_robust_flags(X_ci,    df_ci$class, 0.997,     cutoff_ci_3sd)

df_bands_no_out_extreme <- make_clean_df(df,    dd_bands_3sd)  # appendix only
df_ci_no_out_extreme    <- make_clean_df(df_ci, dd_ci_3sd)     # appendix only

# 5SD — main analysis datasets
cutoff_bands_5sd <- sqrt(qchisq(1 - 5.7e-7, df = p_bands))
cutoff_ci_5sd    <- sqrt(qchisq(1 - 5.7e-7, df = p_ci))

dd_bands_5sd <- make_robust_flags(X_bands, df$class,    0.9999994, cutoff_bands_5sd)
dd_ci_5sd    <- make_robust_flags(X_ci,    df_ci$class, 0.9999994, cutoff_ci_5sd)

df_bands_no_out_5sd <- make_clean_df(df,    dd_bands_5sd)
df_ci_no_out_5sd    <- make_clean_df(df_ci, dd_ci_5sd)

cat(sprintf("Full data          — Bands n: %d  |  CI n: %d\n", nrow(df),                  nrow(df_ci)))
cat(sprintf("2SD / 95%% MCD     — Bands n: %d  |  CI n: %d\n", nrow(df_bands_no_out),     nrow(df_ci_no_out)))
cat(sprintf("3SD / 99.7%% (app) — Bands n: %d  |  CI n: %d\n", nrow(df_bands_no_out_extreme), nrow(df_ci_no_out_extreme)))
cat(sprintf("5SD hard flag      — Bands n: %d  |  CI n: %d\n", nrow(df_bands_no_out_5sd), nrow(df_ci_no_out_5sd)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — PCA #####
# ═══════════════════════════════════════════════════════════════════════════════

run_pca <- function(data, cols) prcomp(data[, cols], center = TRUE, scale. = TRUE)

pca_variance <- function(pca_obj) {
  var_exp <- pca_obj$sdev^2
  tibble(pc = seq_along(var_exp), variance = var_exp,
         prop = var_exp / sum(var_exp), cumulative = cumsum(var_exp / sum(var_exp)))
}

make_scree <- function(var_df, title, color, data = NULL, cols = NULL, n_sim = 500) {
  n_pc    <- nrow(var_df)
  max_eig <- max(var_df$variance) * 1.2
  scale_r <- max_eig
  
  # Parallel analysis if data and cols supplied
  pa_line <- NULL
  if (!is.null(data) && !is.null(cols)) {
    X <- scale(as.matrix(data[, cols]))
    n <- nrow(X); p <- ncol(X)
    set.seed(437)
    sim_eigvals <- replicate(n_sim, {
      eigen(cor(matrix(rnorm(n * p), n, p)))$values
    })
    pa_95 <- apply(sim_eigvals, 1, quantile, probs = 0.95)
    pa_df  <- tibble(pc = seq_len(p), pa_95 = pa_95)
    n_pa   <- sum(var_df$variance > pa_95)
    pa_line <- list(df = pa_df, n_pa = n_pa)
  }
  
  n_kaiser <- sum(var_df$variance > 1)
  
  p_out <- ggplot(var_df, aes(x = pc)) +
    geom_col(aes(y = variance), fill = color, alpha = 0.75, width = 0.5) +
    geom_text(aes(y = variance - max_eig * 0.05, label = round(variance, 2)),
              size = 2.6, vjust = 0, color = "grey30") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "#B2182B") +
    annotate("text", x = n_pc + 0.45, y = 1, label = "Kaiser Criterion",
             hjust = 1, vjust = -0.3, size = 2.4, color = "#B2182B") +
    geom_line(aes(y = cumulative * scale_r), color = "grey50", linewidth = 0.6, linetype = "dotted") +
    geom_point(aes(y = cumulative * scale_r), color = "grey50", size = 1.8) +
    geom_text(aes(y = cumulative * scale_r + max_eig * 0.03,
                  label = paste0(round(cumulative * 100, 1), "%")),
              size = 2.4, color = "grey50", vjust = 0) +
    scale_x_continuous(breaks = seq_len(n_pc), labels = paste0("PC", seq_len(n_pc))) +
    scale_y_continuous(name = "Eigenvalue", limits = c(0, max_eig),
                       sec.axis = sec_axis(transform = ~ . / scale_r,
                                           name = "Cumulative variance explained",
                                           labels = scales::percent_format(accuracy = 1))) +
    labs(x = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text = element_text(size = 8, colour = "grey50"),
          axis.title.y.left  = element_text(size = 8, colour = "grey40"),
          axis.title.y.right = element_text(size = 8, colour = "grey50"),
          plot.title = element_text(size = 9, face = "plain"))
  
  if (!is.null(pa_line)) {
    p_out <- p_out +
      geom_line(data  = pa_line$df, aes(x = pc, y = pa_95),
                color = "#7B2D8B", linewidth = 0.6, linetype = "longdash") +
      geom_point(data = pa_line$df, aes(x = pc, y = pa_95),
                 color = "#7B2D8B", size = 1.8, shape = 1) +
      annotate("text", x = 1 - 0.45, y = pa_line$df$pa_95[1],
               label = "Parallel Analysis", hjust = 0, vjust = -0.3,
               size = 2.4, color = "#7B2D8B")
  }
  
  p_out
}

# ── 3a. PCA — raw bands (full + 2SD + 5SD) ───────────────────────────────────
pca_bands_full  <- run_pca(df,                 photo_cols)
pca_bands_clean <- run_pca(df_bands_no_out,    photo_cols)
pca_bands_5sd   <- run_pca(df_bands_no_out_5sd, photo_cols)

var_bands_full  <- pca_variance(pca_bands_full)
var_bands_clean <- pca_variance(pca_bands_clean)
var_bands_5sd   <- pca_variance(pca_bands_5sd)

p_scree_bands_full  <- make_scree(var_bands_full,  paste0("Raw bands — full (n=", nrow(df), ")"),                class_colors["GALAXY"], data = df,                  cols = photo_cols)
p_scree_bands_clean <- make_scree(var_bands_clean, paste0("Raw bands — 2SD/95% (n=", nrow(df_bands_no_out), ")"), class_colors["GALAXY"], data = df_bands_no_out,     cols = photo_cols)
p_scree_bands_5sd   <- make_scree(var_bands_5sd,   paste0("Raw bands — 5SD (n=", nrow(df_bands_no_out_5sd), ")"),  class_colors["GALAXY"], data = df_bands_no_out_5sd, cols = photo_cols)

# ── 3b. PCA — color indices (full + 2SD + 5SD) ───────────────────────────────
pca_ci_full  <- run_pca(df_ci,              ci_cols)
pca_ci_clean <- run_pca(df_ci_no_out,       ci_cols)
pca_ci_5sd   <- run_pca(df_ci_no_out_5sd,   ci_cols)

var_ci_full  <- pca_variance(pca_ci_full)
var_ci_clean <- pca_variance(pca_ci_clean)
var_ci_5sd   <- pca_variance(pca_ci_5sd)

p_scree_ci_full  <- make_scree(var_ci_full,  paste0("Color indices — full (n=", nrow(df_ci), ")"),              class_colors["QSO"], data = df_ci,              cols = ci_cols)
p_scree_ci_clean <- make_scree(var_ci_clean, paste0("Color indices — 2SD/95% (n=", nrow(df_ci_no_out), ")"),    class_colors["QSO"], data = df_ci_no_out,       cols = ci_cols)
p_scree_ci_5sd   <- make_scree(var_ci_5sd,   paste0("Color indices — 5SD (n=", nrow(df_ci_no_out_5sd), ")"),    class_colors["QSO"], data = df_ci_no_out_5sd,   cols = ci_cols)

# ── 3c. Combined scree plots (2x2: full and 5SD only) ────────────────────────
p_scree_grid <- (p_scree_bands_full | p_scree_bands_5sd) /
  (p_scree_ci_full  | p_scree_ci_5sd)


# ── 3e. Score biplots ─────────────────────────────────────────────────────────
make_biplot <- function(pca_obj, class_vec, var_df, title,
                        circle = NULL, pc1_name = NULL, pc2_name = NULL) {
  # circle = list(cx, cy, rx, ry) for an ellipse, or list(cx, cy, r) for a circle
  scores <- as_tibble(pca_obj$x[, 1:2]) %>% mutate(class = class_vec)
  loads  <- as_tibble(pca_obj$rotation[, 1:2]) %>%
    mutate(variable = rownames(pca_obj$rotation))
  score_range <- max(abs(c(scores$PC1, scores$PC2)))
  load_scale  <- score_range * 0.45 / max(abs(c(loads$PC1, loads$PC2)))
  loads <- loads %>% mutate(PC1s = PC1 * load_scale, PC2s = PC2 * load_scale)
  xlab <- if (!is.null(pc1_name))
    paste0("PC1 — ", pc1_name, " (", round(var_df$prop[1] * 100, 1), "%)")
  else paste0("PC1 (", round(var_df$prop[1] * 100, 1), "%)")
  ylab <- if (!is.null(pc2_name))
    paste0("PC2 — ", pc2_name, " (", round(var_df$prop[2] * 100, 1), "%)")
  else paste0("PC2 (", round(var_df$prop[2] * 100, 1), "%)")
  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_point(data = scores %>% group_by(class),
               aes(x = PC1, y = PC2, color = class), alpha = 0.25, size = 0.7) +
    geom_point(data = scores %>% group_by(class) %>%
                 summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop"),
               aes(x = PC1, y = PC2, color = class),
               size = 5, shape = 21, fill = "white", stroke = 2) +
    geom_segment(data = loads, aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                 linewidth = 0.6, color = "grey20") +
    geom_text(data = loads, aes(x = PC1s * 1.12, y = PC2s * 1.12, label = variable),
              size = 2.8, color = "grey20", fontface = "bold") +
    scale_color_manual(values = class_colors, name = "Midpoint") +
    labs(x = xlab, y = ylab, title = title, color = "Midpoint") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right", legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8), panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey93"),
          axis.text = element_text(size = 8, colour = "grey50"),
          plot.title = element_text(size = 9, face = "plain"))
  if (!is.null(circle)) {
    rx <- if (!is.null(circle$rx)) circle$rx else circle$r
    ry <- if (!is.null(circle$ry)) circle$ry else circle$r
    t  <- seq(0, 2 * pi, length.out = 200)
    p  <- p + annotate("path",
                       x = circle$cx + rx * cos(t),
                       y = circle$cy + ry * sin(t),
                       color = "#B2182B", linewidth = 1.2)
  }
  p
}

p_bi_bands_full  <- make_biplot(pca_bands_full,  df$class,               var_bands_full,  paste0("Raw bands — full (n=",   nrow(df), ")"))
p_bi_bands_clean <- make_biplot(pca_bands_clean, df_bands_no_out$class,  var_bands_clean, paste0("Raw bands — 2SD/95% (n=", nrow(df_bands_no_out), ")"))
p_bi_bands_5sd   <- make_biplot(pca_bands_5sd,   df_bands_no_out_5sd$class, var_bands_5sd,
                                paste0("Raw bands — 5SD (n=", nrow(df_bands_no_out_5sd), ")"),
                                pc1_name = "Distance/Brightness", pc2_name = "UV Excess")

p_bi_ci_full  <- make_biplot(pca_ci_full,  df_ci$class,            var_ci_full,  paste0("Color indices — full (n=",   nrow(df_ci), ")"))
p_bi_ci_clean <- make_biplot(pca_ci_clean, df_ci_no_out$class,     var_ci_clean, paste0("Color indices — 2SD/95% (n=", nrow(df_ci_no_out), ")"))
p_bi_ci_5sd   <- make_biplot(pca_ci_5sd,   df_ci_no_out_5sd$class, var_ci_5sd,
                             paste0("Color indices — 5SD (n=", nrow(df_ci_no_out_5sd), ")"),
                             circle = list(cx = 5.5, cy = 0.5, rx = 1.8, ry = 1.2),
                             pc1_name = "Redness", pc2_name = "UV Excess")

# ── 3f. Combined biplots (2x2: full and 5SD only) ────────────────────────────
p_biplot_grid <- (p_bi_bands_full | p_bi_bands_5sd) /
  (p_bi_ci_full  | p_bi_ci_5sd) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ── 3g. PCA loading heatmaps ─────────────────────────────────────────────────
make_pca_loadings_heatmap <- function(pca_obj, title) {
  loads <- as_tibble(pca_obj$rotation[, 1:2]) %>%
    mutate(variable = rownames(pca_obj$rotation)) %>%
    pivot_longer(-variable, names_to = "component", values_to = "loading") %>%
    mutate(variable = factor(variable, levels = rev(unique(variable))))
  
  ggplot(loads, aes(x = component, y = variable, fill = loading)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", loading),
                  color = abs(loading) > 0.5),
              size = 2.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Loading") +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey30"),
                       guide = "none") +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 8) +
    theme(panel.grid        = element_blank(),
          axis.text.x       = element_text(size = 7, colour = "grey40"),
          axis.text.y       = element_text(size = 7, colour = "grey40"),
          legend.position   = "none",
          plot.title        = element_text(size = 8, face = "plain"))
}

p_load_bands_full <- make_pca_loadings_heatmap(pca_bands_full, "Raw bands — full")
p_load_bands_5sd  <- make_pca_loadings_heatmap(pca_bands_5sd,  "Raw bands — 5SD")
p_load_ci_full    <- make_pca_loadings_heatmap(pca_ci_full,    "Color indices — full")
p_load_ci_5sd     <- make_pca_loadings_heatmap(pca_ci_5sd,     "Color indices — 5SD")

p_load_grid <- (p_load_bands_full | p_load_bands_5sd) /
  (p_load_ci_full | p_load_ci_5sd)

# ── 3h. Scree (A) + Loadings (B) + Biplot (C) stacked ───────────────────────
p_scree_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_scree_grid)) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_load_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_load_grid)) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_biplot_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_biplot_grid)) +
  draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

plot_grid(p_scree_labelled, p_load_labelled, p_biplot_labelled,
          nrow = 1, rel_widths = c(1, 0.35, 1))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — ASSUMPTION TESTING FOR FA #####
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n--- KMO: Raw bands — full ---\n");     KMO(cor(df[, photo_cols]))
cat("\n--- KMO: Raw bands — 2SD/95% ---\n");  KMO(cor(df_bands_no_out[, photo_cols]))
cat("\n--- KMO: Raw bands — 5SD ---\n");      KMO(cor(df_bands_no_out_5sd[, photo_cols]))
cat("\n--- KMO: Color indices — full ---\n"); KMO(cor(df_ci[, ci_cols]))
cat("\n--- KMO: Color indices — 2SD/95% ---\n"); KMO(cor(df_ci_no_out[, ci_cols]))
cat("\n--- KMO: Color indices — 5SD ---\n");  KMO(cor(df_ci_no_out_5sd[, ci_cols]))

# ── Bartlett's test of sphericity ─────────────────────────────────────────────
run_bartlett <- function(data, cols, label) {
  R <- cor(scale(as.matrix(data[, cols])))
  n <- nrow(data)
  p <- length(cols)
  # chi-sq statistic: -(n - 1 - (2p+5)/6) * log(det(R))
  chi_sq <- -(n - 1 - (2 * p + 5) / 6) * log(det(R))
  df_val <- p * (p - 1) / 2
  p_val  <- pchisq(chi_sq, df = df_val, lower.tail = FALSE)
  cat(sprintf("%-35s  chi2 = %8.2f  df = %2d  p = %.2e\n",
              label, chi_sq, df_val, p_val))
}

cat("\n--- Bartlett's test of sphericity ---\n")
run_bartlett(df,                  photo_cols, "Raw bands — full")
run_bartlett(df_bands_no_out,     photo_cols, "Raw bands — 2SD/95%")
run_bartlett(df_bands_no_out_5sd, photo_cols, "Raw bands — 5SD")
run_bartlett(df_ci,               ci_cols,    "Color indices — full")
run_bartlett(df_ci_no_out,        ci_cols,    "Color indices — 2SD/95%")
run_bartlett(df_ci_no_out_5sd,    ci_cols,    "Color indices — 5SD")

set.seed(437)
v   <- rnorm(100000)
EGv <- mean(log(cosh(v)))

negentropy_vars <- function(data, cols, label) {
  X <- scale(as.matrix(data[, cols]))
  neg <- sapply(cols, function(col) {
    s <- X[, col]; round((mean(log(cosh(s))) - EGv)^2, 6)
  })
  tibble(dataset = label, variable = cols, negentropy = neg)
}

neg_results <- bind_rows(
  negentropy_vars(df,                  photo_cols, "Raw bands — full"),
  negentropy_vars(df_bands_no_out,     photo_cols, "Raw bands — 2SD/95%"),
  negentropy_vars(df_bands_no_out_5sd, photo_cols, "Raw bands — 5SD"),
  negentropy_vars(df_ci,               ci_cols,    "Color indices — full"),
  negentropy_vars(df_ci_no_out,        ci_cols,    "Color indices — 2SD/95%"),
  negentropy_vars(df_ci_no_out_5sd,    ci_cols,    "Color indices — 5SD")
)
cat("\n--- Negentropy per variable ---\n")
print(neg_results %>% pivot_wider(names_from = variable, values_from = negentropy), n = Inf)

# ── Mardia's test of multivariate normality ───────────────────────────────────
run_mardia <- function(data, cols, label) {
  cat(sprintf("\n--- Mardia's test: %s ---\n", label))
  result <- mvn(data[, cols], mvnTest = "mardia")
  print(result$multivariateNormality)
}

run_mardia(df,                  photo_cols, "Raw bands — full")
run_mardia(df_bands_no_out,     photo_cols, "Raw bands — 2SD/95%")
run_mardia(df_bands_no_out_5sd, photo_cols, "Raw bands — 5SD")
run_mardia(df_ci,               ci_cols,    "Color indices — full")
run_mardia(df_ci_no_out,        ci_cols,    "Color indices — 2SD/95%")
run_mardia(df_ci_no_out_5sd,    ci_cols,    "Color indices — 5SD")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — FACTOR ANALYSIS (ML, 2 factors) #####
# ═══════════════════════════════════════════════════════════════════════════════

make_fa_scree <- function(data, cols, title, color, n_sim = 500) {
  X      <- scale(as.matrix(data[, cols]))
  n      <- nrow(X); p <- ncol(X)
  eigval <- eigen(cor(X))$values
  
  set.seed(437)
  sim_eigvals <- replicate(n_sim, eigen(cor(matrix(rnorm(n * p), n, p)))$values)
  pa_95    <- apply(sim_eigvals, 1, quantile, probs = 0.95)
  pa_df    <- tibble(factor = seq_len(p), pa_95 = pa_95)
  n_kaiser <- sum(eigval > 1)
  n_pa     <- sum(eigval > pa_95)
  
  max_eig  <- max(eigval) * 1.2
  scale_r  <- max_eig
  
  # Build cumulative variance for secondary axis
  var_df <- tibble(
    factor     = seq_len(p),
    observed   = eigval,
    prop       = eigval / sum(eigval),
    cumulative = cumsum(eigval / sum(eigval))
  )
  
  ggplot(var_df, aes(x = factor)) +
    geom_col(aes(y = observed), fill = color, alpha = 0.75, width = 0.5) +
    geom_text(aes(y = observed - max_eig * 0.05, label = round(observed, 2)),
              size = 2.6, vjust = 0, color = "grey30") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "#B2182B") +
    annotate("text", x = p + 0.45, y = 1, label = "Kaiser Criterion",
             hjust = 1, vjust = -0.3, size = 2.4, color = "#B2182B") +
    geom_line(aes(y = cumulative * scale_r),
              color = "grey50", linewidth = 0.6, linetype = "dotted") +
    geom_point(aes(y = cumulative * scale_r), color = "grey50", size = 1.8) +
    geom_text(aes(y = cumulative * scale_r + max_eig * 0.03,
                  label = paste0(round(cumulative * 100, 1), "%")),
              size = 2.4, color = "grey50", vjust = 0) +
    geom_line(data = pa_df, aes(x = factor, y = pa_95),
              color = "#7B2D8B", linewidth = 0.6, linetype = "longdash") +
    geom_point(data = pa_df, aes(x = factor, y = pa_95),
               color = "#7B2D8B", size = 1.8, shape = 1) +
    annotate("text", x = 1 - 0.45, y = pa_95[1],
             label = "Parallel Analysis", hjust = 0, vjust = -0.3,
             size = 2.4, color = "#7B2D8B") +
    scale_x_continuous(breaks = seq_len(p), labels = paste0("F", seq_len(p))) +
    scale_y_continuous(
      name     = "Eigenvalue",
      limits   = c(0, max_eig),
      sec.axis = sec_axis(transform = ~ . / scale_r,
                          name      = "Cumulative variance explained",
                          labels    = scales::percent_format(accuracy = 1))
    ) +
    labs(x = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor    = element_blank(),
          panel.grid.major.x  = element_blank(),
          axis.text            = element_text(size = 8, colour = "grey50"),
          axis.title.y.left   = element_text(size = 8, colour = "grey40"),
          axis.title.y.right  = element_text(size = 8, colour = "grey50"),
          plot.title           = element_text(size = 9, face = "plain"))
}

# ── 5a. FA scree plots (2x3) ─────────────────────────────────────────────────
p_fa_scree_bands_full  <- make_fa_scree(df,                  photo_cols, paste0("Raw bands — full (n=",   nrow(df), ")"),               class_colors["GALAXY"])
p_fa_scree_bands_clean <- make_fa_scree(df_bands_no_out,     photo_cols, paste0("Raw bands — 2SD/95% (n=", nrow(df_bands_no_out), ")"), class_colors["GALAXY"])
p_fa_scree_bands_5sd   <- make_fa_scree(df_bands_no_out_5sd, photo_cols, paste0("Raw bands — 5SD (n=",    nrow(df_bands_no_out_5sd), ")"), class_colors["GALAXY"])
p_fa_scree_ci_full     <- make_fa_scree(df_ci,               ci_cols,    paste0("Color indices — full (n=",   nrow(df_ci), ")"),              class_colors["QSO"])
p_fa_scree_ci_clean    <- make_fa_scree(df_ci_no_out,        ci_cols,    paste0("Color indices — 2SD/95% (n=", nrow(df_ci_no_out), ")"),    class_colors["QSO"])
p_fa_scree_ci_5sd      <- make_fa_scree(df_ci_no_out_5sd,    ci_cols,    paste0("Color indices — 5SD (n=",    nrow(df_ci_no_out_5sd), ")"),  class_colors["QSO"])

(p_fa_scree_bands_full | p_fa_scree_bands_clean | p_fa_scree_bands_5sd) /
  (p_fa_scree_ci_full  | p_fa_scree_ci_clean   | p_fa_scree_ci_5sd) +
  plot_annotation(title = "FA scree plots with parallel analysis",
                  caption = "Solid = observed  |  Dashed = PA 95th percentile (500 sims)  |  Red = Kaiser",
                  theme = theme(plot.title   = element_text(size = 12),
                                plot.caption = element_text(size = 7.5, colour = "grey50")))

# ── Helper: FA loadings ───────────────────────────────────────────────────────
run_fa_loadings <- function(data, cols, n_factors, dataset_label) {
  X <- scale(as.matrix(data[, cols]))
  map2_dfr(c("none","varimax","promax"), c("No Rotation","Varimax","Promax"),
           function(rot, rot_lab) {
             fa_obj <- fa(X, nfactors = n_factors, rotate = rot, fm = "ml", scores = "regression")
             loads  <- as.data.frame(unclass(fa_obj$loadings))
             loads$variable <- rownames(loads)
             loads %>% pivot_longer(-variable, names_to = "factor", values_to = "loading") %>%
               mutate(rotation = rot_lab, dataset = dataset_label,
                      factor = recode(factor, MR1 = "ML1", MR2 = "ML2"))
           })
}

N_FACTORS <- 2

fa_bands_full  <- run_fa_loadings(df,                  photo_cols, N_FACTORS, "Raw bands — full")
fa_bands_clean <- run_fa_loadings(df_bands_no_out,     photo_cols, N_FACTORS, "Raw bands — 2SD/95%")
fa_bands_5sd   <- run_fa_loadings(df_bands_no_out_5sd, photo_cols, N_FACTORS, "Raw bands — 5SD")
fa_ci_full     <- run_fa_loadings(df_ci,               ci_cols,    N_FACTORS, "Color indices — full")
fa_ci_clean    <- run_fa_loadings(df_ci_no_out,        ci_cols,    N_FACTORS, "Color indices — 2SD/95%")
fa_ci_5sd      <- run_fa_loadings(df_ci_no_out_5sd,    ci_cols,    N_FACTORS, "Color indices — 5SD")

make_fa_heatmap <- function(fa_df, title) {
  fa_df <- fa_df %>%
    mutate(rotation = factor(rotation, levels = c("No Rotation", "Varimax", "Promax")),
           variable = factor(variable, levels = rev(unique(variable))))
  ggplot(fa_df, aes(x = factor, y = variable, fill = loading)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", loading), color = abs(loading) > 0.5),
              size = 2.3, fontface = "plain") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1), name = "Loading") +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey30"), guide = "none") +
    facet_wrap(~ rotation, nrow = 1) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 8) +
    theme(strip.text = element_text(size = 7, face = "plain"), panel.grid = element_blank(),
          axis.text.x = element_text(size = 7, colour = "grey40"),
          axis.text.y = element_text(size = 7, colour = "grey40"),
          legend.position = "none",
          plot.title = element_text(size = 8, face = "plain"))
}


# ── 5c. Communalities, uniquenesses, Phi ─────────────────────────────────────
print_fa_diagnostics <- function(data, cols, label) {
  X          <- scale(as.matrix(data[, cols]))
  fa_varimax <- fa(X, nfactors = N_FACTORS, rotate = "varimax", fm = "ml")
  fa_promax  <- fa(X, nfactors = N_FACTORS, rotate = "promax",  fm = "ml")
  result <- tibble(variable    = cols,
                   communality = round(fa_varimax$communality, 3),
                   uniqueness  = round(fa_varimax$uniquenesses, 3))
  cat(sprintf("\n--- Communalities & Uniquenesses: %s ---\n", label))
  print(result, n = Inf)
  if (!is.null(fa_promax$Phi))
    cat(sprintf("Promax factor correlation (Phi): %.3f\n", fa_promax$Phi[1, 2]))
}


# ── 5d. FA score biplots ──────────────────────────────────────────────────────
make_fa_biplot <- function(data, cols, class_vec, n_factors, title,
                           rotation = "varimax", circle = NULL,
                           f1_name = NULL, f2_name = NULL) {
  X      <- scale(as.matrix(data[, cols]))
  fa_obj <- fa(X, nfactors = n_factors, rotate = rotation, fm = "ml", scores = "regression")
  scores <- as_tibble(fa_obj$scores) %>% setNames(paste0("F", seq_len(n_factors))) %>%
    mutate(class = class_vec)
  load_mat <- if (rotation == "promax" && !is.null(fa_obj$Structure))
    as.data.frame(fa_obj$Structure[, 1:2]) else
      as.data.frame(unclass(fa_obj$loadings)[, 1:2])
  loads <- load_mat %>% setNames(c("F1", "F2")) %>% mutate(variable = cols)
  score_range <- max(abs(c(scores$F1, scores$F2)))
  load_scale  <- score_range * 0.45 / max(abs(c(loads$F1, loads$F2)))
  loads <- loads %>% mutate(F1s = F1 * load_scale, F2s = F2 * load_scale)
  var_exp <- fa_obj$Vaccounted
  xlab <- if (!is.null(f1_name))
    paste0("F1 — ", f1_name, " (", round(var_exp["Proportion Var", 1] * 100, 1), "%)")
  else paste0("F1 (", round(var_exp["Proportion Var", 1] * 100, 1), "%)")
  ylab <- if (!is.null(f2_name))
    paste0("F2 — ", f2_name, " (", round(var_exp["Proportion Var", 2] * 100, 1), "%)")
  else paste0("F2 (", round(var_exp["Proportion Var", 2] * 100, 1), "%)")
  p <- ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_point(data = scores %>% group_by(class),
               aes(x = F1, y = F2, color = class), alpha = 0.25, size = 0.7) +
    geom_point(data = scores %>% group_by(class) %>%
                 summarise(F1 = mean(F1), F2 = mean(F2), .groups = "drop"),
               aes(x = F1, y = F2, color = class),
               size = 5, shape = 21, fill = "white", stroke = 2) +
    geom_segment(data = loads, aes(x = 0, y = 0, xend = F1s, yend = F2s),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                 linewidth = 0.6, color = "grey20") +
    geom_text(data = loads, aes(x = F1s * 1.12, y = F2s * 1.12, label = variable),
              size = 2.8, color = "grey20", fontface = "bold") +
    scale_color_manual(values = class_colors, name = "Midpoint") +
    labs(x = xlab, y = ylab, title = title, color = "Midpoint") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right", legend.key.size = unit(0.4, "cm"),
          legend.text = element_text(size = 8), panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey93"),
          axis.text = element_text(size = 8, colour = "grey50"),
          plot.title = element_text(size = 9, face = "plain"))
  if (!is.null(circle)) {
    rx <- if (!is.null(circle$rx)) circle$rx else circle$r
    ry <- if (!is.null(circle$ry)) circle$ry else circle$r
    t  <- seq(0, 2 * pi, length.out = 200)
    p  <- p + annotate("path",
                       x = circle$cx + rx * cos(t),
                       y = circle$cy + ry * sin(t),
                       color = "#B2182B", linewidth = 1.2)
  }
  p
}

# Varimax biplots
p_fa_bi_bands_full_vm  <- make_fa_biplot(df,                  photo_cols, df$class,               N_FACTORS, paste0("Raw bands — full (n=",   nrow(df), ")"),               "varimax")
p_fa_bi_bands_clean_vm <- make_fa_biplot(df_bands_no_out,     photo_cols, df_bands_no_out$class,  N_FACTORS, paste0("Raw bands — 2SD/95% (n=", nrow(df_bands_no_out), ")"), "varimax")
p_fa_bi_bands_5sd_vm   <- make_fa_biplot(df_bands_no_out_5sd, photo_cols, df_bands_no_out_5sd$class, N_FACTORS, paste0("Raw bands — 5SD (n=",    nrow(df_bands_no_out_5sd), ")"), "varimax",
                                         f1_name = "Redness", f2_name = "Distance/Brightness")
p_fa_bi_ci_full_vm     <- make_fa_biplot(df_ci,               ci_cols,    df_ci$class,             N_FACTORS, paste0("Color indices — full (n=",   nrow(df_ci), ")"),           "varimax")
p_fa_bi_ci_clean_vm    <- make_fa_biplot(df_ci_no_out,        ci_cols,    df_ci_no_out$class,      N_FACTORS, paste0("Color indices — 2SD/95% (n=", nrow(df_ci_no_out), ")"),  "varimax")
p_fa_bi_ci_5sd_vm      <- make_fa_biplot(df_ci_no_out_5sd,    ci_cols,    df_ci_no_out_5sd$class,  N_FACTORS, paste0("Color indices — 5SD (n=",    nrow(df_ci_no_out_5sd), ")"), "varimax",
                                         circle = list(cx = -2.5, cy = -1.7, rx = 0.8, ry = 1.1),
                                         f1_name = "Redness", f2_name = "UV Excess")

# Promax biplots
p_fa_bi_bands_full_pm  <- make_fa_biplot(df,                  photo_cols, df$class,               N_FACTORS, paste0("Raw bands — full (n=",   nrow(df), ")"),               "promax")
p_fa_bi_bands_clean_pm <- make_fa_biplot(df_bands_no_out,     photo_cols, df_bands_no_out$class,  N_FACTORS, paste0("Raw bands — 2SD/95% (n=", nrow(df_bands_no_out), ")"), "promax")
p_fa_bi_bands_5sd_pm   <- make_fa_biplot(df_bands_no_out_5sd, photo_cols, df_bands_no_out_5sd$class, N_FACTORS, paste0("Raw bands — 5SD (n=",    nrow(df_bands_no_out_5sd), ")"), "promax")
p_fa_bi_ci_full_pm     <- make_fa_biplot(df_ci,               ci_cols,    df_ci$class,             N_FACTORS, paste0("Color indices — full (n=",   nrow(df_ci), ")"),           "promax")
p_fa_bi_ci_clean_pm    <- make_fa_biplot(df_ci_no_out,        ci_cols,    df_ci_no_out$class,      N_FACTORS, paste0("Color indices — 2SD/95% (n=", nrow(df_ci_no_out), ")"),  "promax")
p_fa_bi_ci_5sd_pm      <- make_fa_biplot(df_ci_no_out_5sd,    ci_cols,    df_ci_no_out_5sd$class,  N_FACTORS, paste0("Color indices — 5SD (n=",    nrow(df_ci_no_out_5sd), ")"), "promax")

# ── 5e. FA scree (A) + varimax loadings (B) + varimax biplots (C) ────────────
p_fa_scree_grid <- (p_fa_scree_bands_full | p_fa_scree_bands_5sd) /
  (p_fa_scree_ci_full  | p_fa_scree_ci_5sd)

p_fa_load_vm_grid <-
  (make_fa_heatmap(filter(fa_bands_full, rotation == "Varimax"), "Raw bands — full") |
     make_fa_heatmap(filter(fa_bands_5sd,  rotation == "Varimax"), "Raw bands — 5SD")) /
  (make_fa_heatmap(filter(fa_ci_full, rotation == "Varimax"), "Color indices — full") |
     make_fa_heatmap(filter(fa_ci_5sd,  rotation == "Varimax"), "Color indices — 5SD"))
p_fa_bi_vm_grid <- (p_fa_bi_bands_full_vm | p_fa_bi_bands_5sd_vm) /
  (p_fa_bi_ci_full_vm  | p_fa_bi_ci_5sd_vm) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

p_fa_scree_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_fa_scree_grid)) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_fa_load_vm_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_fa_load_vm_grid)) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_fa_bi_vm_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_fa_bi_vm_grid)) +
  draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

plot_grid(p_fa_scree_labelled, p_fa_load_vm_labelled, p_fa_bi_vm_labelled,
          nrow = 1, rel_widths = c(1, 0.35, 1))

# ── 5f. FA scree (A) + promax loadings (B) + promax biplots (C) ──────────────
p_fa_load_pm_grid <-
  (make_fa_heatmap(filter(fa_bands_full, rotation == "Varimax"), "Raw bands — full") |
     make_fa_heatmap(filter(fa_bands_5sd,  rotation == "Varimax"), "Raw bands — 5SD")) /
  (make_fa_heatmap(filter(fa_ci_full, rotation == "Varimax"), "Color indices — full") |
     make_fa_heatmap(filter(fa_ci_5sd,  rotation == "Varimax"), "Color indices — 5SD"))

p_fa_bi_pm_grid <- (p_fa_bi_bands_full_pm | p_fa_bi_bands_5sd_pm) /
  (p_fa_bi_ci_full_pm  | p_fa_bi_ci_5sd_pm) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

p_fa_load_pm_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_fa_load_pm_grid)) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_fa_bi_pm_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_fa_bi_pm_grid)) +
  draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

plot_grid(p_fa_scree_labelled, p_fa_load_pm_labelled, p_fa_bi_pm_labelled,
          nrow = 1, rel_widths = c(1, 0.35, 1))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — LDA & LOGISTIC REGRESSION WITH CROSS-VALIDATION #####
# ═══════════════════════════════════════════════════════════════════════════════
set.seed(437)
K_FOLDS <- 10

make_folds_stratified <- function(y, k) {
  folds <- vector("list", k)
  for (cls in unique(y)) {
    idx    <- which(y == cls); idx <- sample(idx)
    splits <- cut(seq_along(idx), breaks = k, labels = FALSE)
    for (i in seq_len(k)) folds[[i]] <- c(folds[[i]], idx[splits == i])
  }
  folds
}

make_folds_unstratified <- function(y, k) {
  n <- length(y); idx <- sample(seq_len(n))
  splits <- cut(seq_len(n), breaks = k, labels = FALSE)
  map(seq_len(k), ~ idx[splits == .x])
}

run_cv <- function(data, cols, y_col, k, method, dataset_label, stratified = TRUE) {
  df_model <- data %>%
    select(all_of(c(cols, y_col))) %>%
    mutate(across(all_of(cols), ~ scale(.)[,1]))
  y       <- df_model[[y_col]]
  folds   <- if (stratified) make_folds_stratified(y, k) else make_folds_unstratified(y, k)
  classes <- sort(unique(y))
  
  fold_results <- map(seq_len(k), function(i) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df_model)), test_idx)
    train <- df_model[train_idx, ]; test <- df_model[test_idx, ]
    formula <- as.formula(paste(y_col, "~",
                                paste(paste0("`", cols, "`"), collapse = " + ")))
    if (method == "lda") {
      fit <- lda(formula, data = train); preds <- predict(fit, newdata = test)$class
    } else {
      fit <- multinom(formula, data = train, trace = FALSE); preds <- predict(fit, newdata = test)
    }
    acc <- mean(preds == test[[y_col]])
    per_class <- map_dfr(classes, function(cls) {
      truth <- test[[y_col]] == cls; pred <- preds == cls
      tibble(class     = cls,
             recall    = ifelse(sum(truth) > 0, sum(truth & pred) / sum(truth), NA_real_),
             precision = ifelse(sum(pred)  > 0, sum(truth & pred) / sum(pred),  NA_real_))
    })
    cm <- as_tibble(table(truth = test[[y_col]], predicted = preds)) %>%
      mutate(fold = i, method = method, dataset = dataset_label)
    list(
      metrics = tibble(fold = i, method = method, dataset = dataset_label, accuracy = acc) %>%
        bind_cols(per_class %>%
                    pivot_wider(names_from = class,
                                values_from = c(recall, precision), names_glue = "{.value}_{class}")),
      cm = cm
    )
  })
  list(metrics = map_dfr(fold_results, "metrics"), cm = map_dfr(fold_results, "cm"))
}

# ── 6a. Run CV — six datasets × two methods × stratified/unstratified ─────────
cv_datasets <- list(
  list(data = df,                  cols = photo_cols, label = "Raw bands — full"),
  list(data = df_bands_no_out,     cols = photo_cols, label = "Raw bands — 2SD/95%"),
  list(data = df_bands_no_out_5sd, cols = photo_cols, label = "Raw bands — 5SD"),
  list(data = df_ci,               cols = ci_cols,    label = "Color indices — full"),
  list(data = df_ci_no_out,        cols = ci_cols,    label = "Color indices — 2SD/95%"),
  list(data = df_ci_no_out_5sd,    cols = ci_cols,    label = "Color indices — 5SD")
)

run_all_cv <- function(stratified) {
  map(cv_datasets, function(d) {
    list(lda      = run_cv(d$data, d$cols, "class", K_FOLDS, "lda",      d$label, stratified),
         logistic = run_cv(d$data, d$cols, "class", K_FOLDS, "logistic", d$label, stratified))
  })
}

cv_raw_strat   <- run_all_cv(stratified = TRUE)
cv_raw_unstrat <- run_all_cv(stratified = FALSE)

extract_results <- function(cv_raw, sampling) {
  metrics <- map_dfr(cv_raw, function(x) bind_rows(x$lda$metrics, x$logistic$metrics)) %>%
    mutate(sampling = sampling)
  cm <- map_dfr(cv_raw, function(x) bind_rows(x$lda$cm, x$logistic$cm)) %>%
    mutate(sampling = sampling)
  list(metrics = metrics, cm = cm)
}

res_strat   <- extract_results(cv_raw_strat,   "Stratified")
res_unstrat <- extract_results(cv_raw_unstrat, "Unstratified")
cv_results  <- bind_rows(res_strat$metrics, res_unstrat$metrics)
cv_cm_all   <- bind_rows(res_strat$cm,      res_unstrat$cm)

# ── 6b. Summarise ─────────────────────────────────────────────────────────────
cv_summary <- cv_results %>%
  group_by(method, dataset, sampling) %>%
  summarise(mean_acc = mean(accuracy), sd_acc = sd(accuracy),
            across(starts_with("recall_"),    ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
            across(starts_with("precision_"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
            .groups = "drop") %>%
  mutate(method   = recode(method, lda = "LDA", logistic = "Logistic"),
         dataset  = factor(dataset, levels = map_chr(cv_datasets, "label")),
         sampling = factor(sampling, levels = c("Stratified", "Unstratified")))


# ── 6c. Accuracy bar chart — full & 5SD only, stratified ─────────────────────
method_colors <- c(LDA = "#378ADD", Logistic = "#D85A30")

acc_label_map <- c(
  "Raw bands — full"        = "Raw bands\nfull",
  "Raw bands — 2SD/95%"     = "Raw bands\n2SD/95%",
  "Raw bands — 5SD"         = "Raw bands\n5SD",
  "Color indices — full"    = "Color indices\nfull",
  "Color indices — 2SD/95%" = "Color indices\n2SD/95%",
  "Color indices — 5SD"     = "Color indices\n5SD"
)

cv_summary_main <- cv_summary %>%
  filter(dataset %in% names(acc_label_map)) %>%
  mutate(dataset  = recode(dataset, !!!acc_label_map),
         dataset  = factor(dataset, levels = unname(acc_label_map)),
         sampling = factor(sampling, levels = c("Unstratified", "Stratified")))

p_acc_main <- ggplot(cv_summary_main, aes(x = dataset, y = mean_acc, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc),
                position = position_dodge(width = 0.7),
                width = 0.25, linewidth = 0.5, color = "grey30") +
  geom_text(aes(label = paste0(round(mean_acc * 100, 1), "%"),
                y = mean_acc + sd_acc + 0.005),
            position = position_dodge(width = 0.7),
            size = 2.4, vjust = 0, color = "grey20") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  facet_wrap(~ sampling, ncol = 1) +
  labs(x = NULL, y = "Mean accuracy (10-fold CV)", fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position    = "top",
        strip.text         = element_text(size = 9, face = "plain"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x  = element_text(size = 7, colour = "grey40"),
        axis.text.y  = element_text(size = 8, colour = "grey50"),
        plot.caption = element_text(size = 7.5, colour = "grey50"))

# ── 6d. Confusion matrices — full & 5SD, stratified, both methods ─────────────
make_cm_plots <- function(method_label, sampling_label = "Stratified",
                          datasets_keep = NULL) {
  cm_data <- cv_cm_all %>%
    filter(method == tolower(method_label), sampling == sampling_label) %>%
    group_by(dataset, truth, predicted) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(dataset, truth) %>%
    mutate(pct = n / sum(n)) %>% ungroup()
  if (!is.null(datasets_keep))
    cm_data <- filter(cm_data, dataset %in% datasets_keep)
  # Apply display labels with line breaks after filtering
  label_map <- c(
    "Raw bands — full"        = "Raw bands\nfull",
    "Raw bands — 5SD"         = "Raw bands\n5SD",
    "Color indices — full"    = "Color indices\nfull",
    "Color indices — 5SD"     = "Color indices\n5SD",
    "Raw bands — 2SD/95%"     = "Raw bands\n2SD/95%",
    "Color indices — 2SD/95%" = "Color indices\n2SD/95%"
  )
  cm_data <- cm_data %>%
    mutate(dataset   = recode(dataset, !!!label_map),
           dataset   = factor(dataset, levels = unname(label_map)[unname(label_map) %in% unique(recode(cm_data$dataset, !!!label_map))]),
           truth     = factor(truth,     levels = class_order),
           predicted = factor(predicted, levels = class_order),
           label     = paste0(round(pct * 100, 1), "%\n(n=", n, ")"))
  ggplot(cm_data, aes(x = predicted, y = truth, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label, color = pct > 0.55), size = 2.5, lineheight = 1.2) +
    scale_fill_gradient2(low = "#F7F7F7", mid = "#7FBCD2", high = "#2166AC",
                         midpoint = 0.5, limits = c(0, 1), name = "Accuracy") +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey30"), guide = "none") +
    facet_wrap(~ dataset, ncol = 2) +
    labs(x = "Predicted", y = "True",
         title = method_label) +
    theme_minimal(base_size = 9) +
    theme(strip.text        = element_text(size = 7),
          panel.grid        = element_blank(),
          axis.text         = element_text(size = 7, colour = "grey40"),
          legend.key.height = unit(0.5, "cm"),
          legend.key.width  = unit(0.3, "cm"),
          legend.text       = element_text(size = 7),
          plot.title        = element_text(size = 9, face = "plain",
                                           hjust = 0.5))
}

cm_datasets <- c("Raw bands — full", "Raw bands — 5SD",
                 "Color indices — full", "Color indices — 5SD")

p_cm_lda  <- make_cm_plots("LDA",      datasets_keep = cm_datasets)
p_cm_log  <- make_cm_plots("Logistic", datasets_keep = cm_datasets)

p_cm_grid <- (p_cm_lda | p_cm_log) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")

# ── 6e. Combined: bars (A) left, confusion matrices (B) right ────────────────
p_acc_labelled <- ggdraw() +
  draw_plot(p_acc_main) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

p_cm_labelled <- ggdraw() +
  draw_plot(patchworkGrob(p_cm_grid)) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
             fontface = "bold", size = 11)

plot_grid(p_acc_labelled, p_cm_labelled, nrow = 1, rel_widths = c(0.5, 1))

# ── 6f. FA score biplot: prediction correctness overlay ──────────────────────
# Fit logistic on full 5SD datasets to get per-observation predictions
make_pred_biplot <- function(data, cols, title, f1_name = NULL, f2_name = NULL) {
  X      <- scale(as.matrix(data[, cols]))
  y      <- data$class
  
  # FA scores (varimax, 2 factors)
  fa_obj <- fa(X, nfactors = 2, rotate = "varimax", fm = "ml", scores = "regression")
  scores <- as_tibble(fa_obj$scores) %>%
    setNames(c("F1", "F2")) %>%
    mutate(class = y)
  
  # Logistic regression predictions
  df_model <- as_tibble(X) %>%
    setNames(cols) %>%
    mutate(class = y)
  formula  <- as.formula(paste("class ~", paste(paste0("`", cols, "`"), collapse = " + ")))
  fit      <- multinom(formula, data = df_model, trace = FALSE)
  preds    <- predict(fit, newdata = df_model)
  
  scores <- scores %>%
    mutate(predicted = preds,
           correct   = factor(ifelse(class == predicted, "Correct", "Incorrect"),
                              levels = c("Correct", "Incorrect")),
           predicted = factor(predicted, levels = class_order))
  
  # Class centroids
  centroids <- scores %>%
    group_by(class) %>%
    summarise(F1 = mean(F1), F2 = mean(F2), .groups = "drop")
  
  var_exp <- fa_obj$Vaccounted
  xlab <- if (!is.null(f1_name))
    paste0("F1 — ", f1_name, " (", round(var_exp["Proportion Var", 1] * 100, 1), "%)")
  else paste0("F1 (", round(var_exp["Proportion Var", 1] * 100, 1), "%)")
  ylab <- if (!is.null(f2_name))
    paste0("F2 — ", f2_name, " (", round(var_exp["Proportion Var", 2] * 100, 1), "%)")
  else paste0("F2 (", round(var_exp["Proportion Var", 2] * 100, 1), "%)")
  
  correct_colors  <- c(Correct = "grey60", Incorrect = "#B2182B")
  correct_alphas  <- c(Correct = 0.20,     Incorrect = 0.85)
  pred_shapes     <- c(GALAXY = 16, STAR = 17, QSO = 15)
  
  ggplot(scores, aes(x = F1, y = F2)) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey85") +
    geom_point(data = filter(scores, correct == "Correct"),
               aes(shape = predicted), color = "grey60",
               alpha = 0.20, size = 0.8) +
    geom_point(data = filter(scores, correct == "Incorrect"),
               aes(color = class, shape = predicted),
               alpha = 0.85, size = 1.6) +
    geom_point(data = centroids,
               aes(x = F1, y = F2, fill = class), color = "white",
               size = 5, shape = 21, stroke = 2,
               inherit.aes = FALSE) +
    scale_color_manual(values = class_colors, name = "True class") +
    scale_fill_manual(values  = class_colors, name = "True class") +
    scale_shape_manual(values = pred_shapes,  name = "Predicted class") +
    labs(x = xlab, y = ylab, title = title) +
    guides(fill  = guide_legend(override.aes = list(size = 3, stroke = 1.5)),
           color = guide_legend(override.aes = list(size = 2.5, alpha = 1)),
           shape = guide_legend(override.aes = list(size = 2.5, alpha = 1,
                                                    color = "#B2182B"))) +
    theme_minimal(base_size = 10) +
    theme(legend.position  = "right",
          legend.key.size  = unit(0.4, "cm"),
          legend.text      = element_text(size = 8),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey93"),
          axis.text        = element_text(size = 8, colour = "grey50"),
          plot.title       = element_text(size = 9, face = "plain"))
}

# ── Helper: extract FA scores + logistic predictions without plotting ─────────
get_pred_scores <- function(data, cols) {
  X      <- scale(as.matrix(data[, cols]))
  y      <- data$class
  fa_obj <- fa(X, nfactors = 2, rotate = "varimax", fm = "ml", scores = "regression")
  scores <- as_tibble(fa_obj$scores) %>%
    setNames(c("F1", "F2")) %>%
    mutate(class = y)
  df_model <- as_tibble(X) %>% setNames(cols) %>% mutate(class = y)
  formula  <- as.formula(paste("class ~", paste(paste0("`", cols, "`"), collapse = " + ")))
  fit      <- multinom(formula, data = df_model, trace = FALSE)
  preds    <- predict(fit, newdata = df_model)
  scores %>% mutate(predicted = preds,
                    correct   = class == predicted)
}

pred_ci_5sd <- get_pred_scores(df_ci_no_out_5sd, ci_cols)

# Accuracy inside the circle (cx=-2, cy=-2, r=1)
circle_mask <- with(pred_ci_5sd, (F1 - (-2))^2 + (F2 - (-2))^2 <= 1^2)
cat(sprintf("\nAccuracy inside circle: %d / %d (%.1f%%)\n",
            sum(pred_ci_5sd$correct[circle_mask]),
            sum(circle_mask),
            100 * mean(pred_ci_5sd$correct[circle_mask])))
cat(sprintf("Class breakdown inside circle:\n"))
print(table(pred_ci_5sd$class[circle_mask], pred_ci_5sd$predicted[circle_mask]))

p_pred_bands_5sd <- make_pred_biplot(df_bands_no_out_5sd, photo_cols,
                                     paste0("Raw bands — 5SD (n=", nrow(df_bands_no_out_5sd), ")"),
                                     f1_name = "Redness", f2_name = "Distance/Brightness")
p_pred_ci_5sd    <- make_pred_biplot(df_ci_no_out_5sd,    ci_cols,
                                     paste0("Color indices — 5SD (n=", nrow(df_ci_no_out_5sd), ")"),
                                     f1_name = "Redness", f2_name = "UV Excess") 
t <- seq(0, 2 * pi, length.out = 200)
circle_acc <- round(100 * mean(pred_ci_5sd$correct[circle_mask]), 1)
p_pred_ci_5sd <- p_pred_ci_5sd +
  annotate("path",
           x = -2 + 1 * cos(t),
           y = -2 + 1 * sin(t),
           color = "#B2182B", linewidth = 1.2) +
  annotate("text", x = -2, y = -3.2,
           label = paste0(circle_acc, "% accuracy"),
           color = "#B2182B", size = 2.8, fontface = "bold", hjust = 0.5)

(p_pred_bands_5sd | p_pred_ci_5sd) +
  plot_layout(guides = "collect") +
  theme(legend.position = "right")