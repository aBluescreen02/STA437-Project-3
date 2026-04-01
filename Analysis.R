library(MASS)
library(ggcorrplot)
library(patchwork)
library(tidyverse)
library(robustbase)

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

# Matrices used throughout
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
    # Quadrant shading
    annotate("rect", xmin = cutoff, xmax = Inf,  ymin = cutoff, ymax = Inf,
             fill = "#D85A30", alpha = 0.04) +
    annotate("rect", xmin = -Inf,  xmax = cutoff, ymin = cutoff, ymax = Inf,
             fill = "#9B72CF", alpha = 0.04) +
    annotate("rect", xmin = cutoff, xmax = Inf,  ymin = -Inf,  ymax = cutoff,
             fill = "#7FBCD2", alpha = 0.04) +
    # Reference lines
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", linewidth = 0.4, color = "grey50") +
    geom_vline(xintercept = cutoff,
               linetype = "dotted", linewidth = 0.4, color = "grey40") +
    geom_hline(yintercept = cutoff,
               linetype = "dotted", linewidth = 0.4, color = "grey40") +
    # Points — regular beneath, outliers on top
    geom_point(data = filter(dd_df, outlier_type == "Regular"),
               alpha = 0.25, size = 0.9) +
    geom_point(data = filter(dd_df, outlier_type != "Regular"),
               alpha = 0.75, size = 1.6) +
    # Quadrant labels
    annotate("text", x = x_max * 0.97, y = y_max * 0.97,
             label = "True outlier",
             hjust = 1, vjust = 1, size = 2.8,
             color = "#D85A30", fontface = "bold") +
    annotate("text", x = cutoff * 0.95, y = y_max * 0.97,
             label = "Masked outlier\n(classical misses)",
             hjust = 1, vjust = 1, size = 2.8,
             color = "#9B72CF", fontface = "bold") +
    annotate("text", x = x_max * 0.97, y = cutoff * 0.95,
             label = "Swamping\n(false alarm)",
             hjust = 1, vjust = 1, size = 2.8,
             color = "#7FBCD2", fontface = "bold") +
    annotate("text", x = cutoff * 0.95, y = cutoff * 0.95,
             label = "Regular",
             hjust = 1, vjust = 1, size = 2.8,
             color = "grey55", fontface = "plain") +
    scale_color_manual(values = type_colors, labels = legend_labels) +
    scale_shape_manual(values = type_shapes, labels = legend_labels) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    labs(
      x       = "Classical Mahalanobis distance",
      y       = "Robust (MCD) Mahalanobis distance",
      color   = NULL, shape = NULL,
      title   = paste0("D-D plot: ", feature_label),
      caption = paste0(feature_label, "  |  n = ", nrow(dd_df),
                       "  |  MCD h = 95%  |  chi-sq cutoff (p = ", p,
                       ", 97.5%): ", round(cutoff, 2))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "top",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7.5, colour = "grey50"),
      plot.title       = element_text(size = 11)
    )
}

# ── Shared helper: extract outlier flags from a dd tibble ─────────────────────
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

# ── Shared helper: outlier counts by class ────────────────────────────────────
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
# SECTION 1 — RAW PHOTOMETRIC BANDS (u, g, r, i, z)
# ═══════════════════════════════════════════════════════════════════════════════

# ── 1a. Histograms ────────────────────────────────────────────────────────────
band_labels <- c(
  u = "u band (ultraviolet, ~354 nm)",
  g = "g band (green, ~477 nm)",
  r = "r band (red, ~623 nm)",
  i = "i band (near-infrared, ~763 nm)",
  z = "z band (infrared, ~913 nm)"
)

df_bands_long <- df %>%
  select(all_of(c(photo_cols, "class"))) %>%
  pivot_longer(cols = all_of(photo_cols),
               names_to  = "band",
               values_to = "magnitude") %>%
  mutate(
    band  = factor(band,  levels = photo_cols, labels = band_labels),
    class = factor(class, levels = class_order, labels = n_labels)
  )

ggplot(df_bands_long, aes(x = magnitude, fill = class, color = class)) +
  geom_histogram(bins = 40, position = "identity",
                 alpha = 0.45, linewidth = 0.3) +
  scale_fill_manual(values  = setNames(class_colors, n_labels)) +
  scale_color_manual(values = setNames(class_colors, n_labels)) +
  facet_wrap(~ band, ncol = 2, scales = "free") +
  labs(x = "Magnitude (AB system)", y = "Count", fill = NULL, color = NULL) +
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

# ── 1b. Correlation heatmaps (raw bands) ──────────────────────────────────────
make_corrplot_bands <- function(data, title, n = NULL) {
  R     <- cor(data[, photo_cols])
  label <- if (!is.null(n)) paste0(title, " (n = ", scales::comma(n), ")") else title
  ggcorrplot(R,
             method   = "square",
             type     = "full",
             lab      = TRUE,
             lab_size = 3,
             colors   = c("#2166AC", "#F7F7F7", "#B2182B"),
             title    = label,
             ggtheme  = theme_minimal(base_size = 10)) +
    theme(plot.title      = element_text(size = 10, hjust = 0.5),
          legend.position = "none",
          axis.text       = element_text(size = 9))
}

p_bands_all  <- make_corrplot_bands(df,                            "All objects", nrow(df))
p_bands_gal  <- make_corrplot_bands(filter(df, class == "GALAXY"), "Galaxy",      sum(df$class == "GALAXY"))
p_bands_star <- make_corrplot_bands(filter(df, class == "STAR"),   "Star",        sum(df$class == "STAR"))
p_bands_qso  <- make_corrplot_bands(filter(df, class == "QSO"),    "QSO",         sum(df$class == "QSO"))

(p_bands_all | plot_spacer()) / (p_bands_gal | p_bands_star | p_bands_qso) +
  plot_annotation(
    title   = "Photometric band correlations — overall vs. by class",
    caption = "Simpson's paradox check: compare overall correlations to within-class correlations",
    theme   = theme(plot.title   = element_text(size = 12, face = "plain"),
                    plot.caption = element_text(size = 8, colour = "grey50"))
  )

# ── 1c. D-D plot: raw bands — pooled (reference) ─────────────────────────────
p_bands      <- ncol(X_bands)
cutoff_bands <- sqrt(qchisq(0.975, df = p_bands))

d_cl_bands_pooled <- sqrt(mahalanobis(X_bands,
                                      center = colMeans(X_bands),
                                      cov    = cov(X_bands)))
mcd_bands_pooled  <- covMcd(X_bands, alpha = 0.95)
d_ro_bands_pooled <- sqrt(mahalanobis(X_bands,
                                      center = mcd_bands_pooled$center,
                                      cov    = mcd_bands_pooled$cov))

dd_bands_pooled <- tibble(
  classical = d_cl_bands_pooled,
  robust    = d_ro_bands_pooled,
  class     = df$class
)

make_dd_plot(dd_bands_pooled, cutoff_bands, p_bands,
             "Raw photometric bands — pooled (reference)")

dd_bands_pooled_typed <- get_outlier_flags(dd_bands_pooled, cutoff_bands)
cat(sprintf("[Bands pooled] Outliers (true + masked): %d / %d (%.1f%%)\n",
            sum(dd_bands_pooled_typed$outlier), nrow(dd_bands_pooled_typed),
            100 * mean(dd_bands_pooled_typed$outlier)))

# ── 1d. D-D plots: raw bands — within-class MCD ───────────────────────────────

dd_bands_typed <- map_dfr(class_order, function(cls) {
  idx  <- which(df$class == cls)
  Xsub <- X_bands[idx, ]
  
  d_cl <- sqrt(mahalanobis(Xsub, center = colMeans(Xsub), cov = cov(Xsub)))
  
  mcd  <- covMcd(Xsub, alpha = 0.95)
  d_ro <- sqrt(mahalanobis(Xsub, center = mcd$center, cov = mcd$cov))
  
  tibble(classical = d_cl, robust = d_ro, class = cls, row_idx = idx)
}) %>%
  get_outlier_flags(cutoff_bands)

# One D-D plot per class, combined with patchwork
dd_plots_bands <- map(class_order, function(cls) {
  make_dd_plot(filter(dd_bands_typed, class == cls),
               cutoff_bands, p_bands,
               paste0("Raw bands — ", cls, " (n=", sum(df$class == cls), ")"))
})

wrap_plots(dd_plots_bands, ncol = 3) +
  plot_annotation(title = "Within-class D-D plots: raw photometric bands",
                  theme = theme(plot.title = element_text(size = 12)))

cat(sprintf("[Bands] Within-class outliers (true + masked): %d / %d (%.1f%%)\n",
            sum(dd_bands_typed$outlier), nrow(dd_bands_typed),
            100 * mean(dd_bands_typed$outlier)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — ADJACENT COLOR INDICES (u-g, g-r, r-i, i-z)
# ═══════════════════════════════════════════════════════════════════════════════

df_ci <- df %>%
  mutate(
    `u-g` = u - g,
    `g-r` = g - r,
    `r-i` = r - i,
    `i-z` = i - z
  )

X_ci <- as.matrix(df_ci[, ci_cols])

# ── 2a. Histograms ────────────────────────────────────────────────────────────
df_ci_long <- df_ci %>%
  select(all_of(c(ci_cols, "class"))) %>%
  pivot_longer(cols = all_of(ci_cols),
               names_to  = "index",
               values_to = "value") %>%
  mutate(index = factor(index, levels = ci_cols))

p_ci_hist <- ggplot(df_ci_long, aes(x = value, fill = class, color = class)) +
  geom_histogram(bins = 50, position = "identity",
                 alpha = 0.42, linewidth = 0.3) +
  scale_fill_manual(values  = class_colors) +
  scale_color_manual(values = class_colors) +
  facet_wrap(~ index, ncol = 2, scales = "free") +
  labs(x = "Color index (mag)", y = "Count", fill = NULL, color = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position  = "top",
        strip.text       = element_text(size = 9),
        panel.grid.minor = element_blank(),
        axis.text        = element_text(size = 8, colour = "grey50"))

# ── 2b. Correlation heatmaps (color indices) ──────────────────────────────────
make_corrplot_ci <- function(data, title) {
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

p_ci_all  <- make_corrplot_ci(df_ci,                            "All objects")
p_ci_gal  <- make_corrplot_ci(filter(df_ci, class == "GALAXY"), "Galaxy")
p_ci_star <- make_corrplot_ci(filter(df_ci, class == "STAR"),   "Star")
p_ci_qso  <- make_corrplot_ci(filter(df_ci, class == "QSO"),    "QSO")

p_ci_corr <- (p_ci_all | plot_spacer()) / (p_ci_gal | p_ci_star | p_ci_qso) +
  plot_annotation(
    title = "Color index correlations — overall vs. by class",
    theme = theme(plot.title = element_text(size = 11))
  )

p_ci_hist / p_ci_corr + plot_layout(heights = c(1, 1.2))

# ── 2c. D-D plot: color indices — pooled (reference) ─────────────────────────
p_ci      <- ncol(X_ci)
cutoff_ci <- sqrt(qchisq(0.975, df = p_ci))

d_cl_ci_pooled <- sqrt(mahalanobis(X_ci,
                                   center = colMeans(X_ci),
                                   cov    = cov(X_ci)))
mcd_ci_pooled  <- covMcd(X_ci, alpha = 0.95)
d_ro_ci_pooled <- sqrt(mahalanobis(X_ci,
                                   center = mcd_ci_pooled$center,
                                   cov    = mcd_ci_pooled$cov))

dd_ci_pooled <- tibble(
  classical = d_cl_ci_pooled,
  robust    = d_ro_ci_pooled,
  class     = df_ci$class
)

make_dd_plot(dd_ci_pooled, cutoff_ci, p_ci,
             "Adjacent color indices — pooled (reference)")

dd_ci_pooled_typed <- get_outlier_flags(dd_ci_pooled, cutoff_ci)
cat(sprintf("[CI pooled] Outliers (true + masked): %d / %d (%.1f%%)\n",
            sum(dd_ci_pooled_typed$outlier), nrow(dd_ci_pooled_typed),
            100 * mean(dd_ci_pooled_typed$outlier)))

# ── 2d. D-D plots: color indices — within-class MCD ──────────────────────────

dd_ci_typed <- map_dfr(class_order, function(cls) {
  idx  <- which(df_ci$class == cls)
  Xsub <- X_ci[idx, ]
  
  d_cl <- sqrt(mahalanobis(Xsub, center = colMeans(Xsub), cov = cov(Xsub)))
  
  mcd  <- covMcd(Xsub, alpha = 0.95)
  d_ro <- sqrt(mahalanobis(Xsub, center = mcd$center, cov = mcd$cov))
  
  tibble(classical = d_cl, robust = d_ro, class = cls, row_idx = idx)
}) %>%
  get_outlier_flags(cutoff_ci)

dd_plots_ci <- map(class_order, function(cls) {
  make_dd_plot(filter(dd_ci_typed, class == cls),
               cutoff_ci, p_ci,
               paste0("Color indices — ", cls, " (n=", sum(df_ci$class == cls), ")"))
})

wrap_plots(dd_plots_ci, ncol = 3) +
  plot_annotation(title = "Within-class D-D plots: adjacent color indices",
                  theme = theme(plot.title = element_text(size = 12)))

cat(sprintf("[CI]    Within-class outliers (true + masked): %d / %d (%.1f%%)\n",
            sum(dd_ci_typed$outlier), nrow(dd_ci_typed),
            100 * mean(dd_ci_typed$outlier)))

# ── Outlier counts by class — pooled ─────────────────────────────────────────
cat("\n--- Pooled outlier summary ---\n")
bind_rows(
  summarise_outliers(dd_bands_pooled_typed, "Raw bands (pooled)"),
  summarise_outliers(dd_ci_pooled_typed,    "Color indices (pooled)")
) %>%
  select(feature_set, class,
         any_of(c("True outlier", "Masked outlier", "Swamping", "Regular")),
         outliers, total, pct) %>%
  print(n = Inf)

# ── Outlier counts by class — within-class ────────────────────────────────────
cat("\n--- Within-class outlier summary ---\n")
bind_rows(
  summarise_outliers(dd_bands_typed, "Raw bands (within-class)"),
  summarise_outliers(dd_ci_typed,    "Color indices (within-class)")
) %>%
  select(feature_set, class,
         any_of(c("True outlier", "Masked outlier", "Swamping", "Regular")),
         outliers, total, pct) %>%
  print(n = Inf)

# ── Carry forward clean datasets for downstream analysis ──────────────────────

# Within-class results are in row order per class; re-join on row_idx
bands_flags <- dd_bands_typed %>%
  select(row_idx, outlier_type, outlier, d_classical = classical, d_robust = robust) %>%
  arrange(row_idx)

ci_flags <- dd_ci_typed %>%
  select(row_idx, outlier_type, outlier, d_classical = classical, d_robust = robust) %>%
  arrange(row_idx)

df_bands_clean  <- df    %>%
  mutate(row_idx = row_number()) %>%
  left_join(bands_flags, by = "row_idx") %>%
  select(-row_idx)

df_ci_clean     <- df_ci %>%
  mutate(row_idx = row_number()) %>%
  left_join(ci_flags, by = "row_idx") %>%
  select(-row_idx)

df_bands_no_out <- filter(df_bands_clean, !outlier)
df_ci_no_out    <- filter(df_ci_clean,    !outlier)

cat(sprintf("Bands — clean n: %d  |  CI — clean n: %d\n",
            nrow(df_bands_no_out), nrow(df_ci_no_out)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — PCA
# ═══════════════════════════════════════════════════════════════════════════════

# ── Helper: run PCA and return tidy variance explained tibble ─────────────────
run_pca <- function(data, cols) {
  X <- data[, cols]
  prcomp(X, center = TRUE, scale. = TRUE)
}

pca_variance <- function(pca_obj) {
  var_exp <- pca_obj$sdev^2
  tibble(
    pc        = seq_along(var_exp),
    variance  = var_exp,
    prop      = var_exp / sum(var_exp),
    cumulative = cumsum(var_exp / sum(var_exp))
  )
}

# ── Helper: scree plot ────────────────────────────────────────────────────────
# Left axis:  eigenvalue (variance)
# Right axis: cumulative % variance explained
# Kaiser criterion: horizontal line at eigenvalue = 1
make_scree <- function(var_df, title, color) {
  n_pc    <- nrow(var_df)
  max_eig <- max(var_df$variance) * 1.2   # headroom above tallest bar
  
  # Scale factor to map cumulative proportion [0,1] onto eigenvalue axis [0, max_eig]
  scale_r <- max_eig
  
  ggplot(var_df, aes(x = pc)) +
    # Eigenvalue bars
    geom_col(aes(y = variance), fill = color, alpha = 0.75, width = 0.5) +
    # Eigenvalue labels on bars
    geom_text(aes(y = variance + max_eig * 0.07,
                  label = round(variance, 2)),
              size = 2.6, vjust = 0, color = "grey30") +
    # Kaiser criterion: eigenvalue = 1
    geom_hline(yintercept = 1, linetype = "dashed",
               linewidth = 0.5, color = "#B2182B") +
    annotate("text", x = n_pc + 0.45, y = 1,
             label = "Kaiser", hjust = 1, vjust = -0.3,
             size = 2.4, color = "#B2182B") +
    # Cumulative % line (scaled to eigenvalue axis)
    geom_line(aes(y = cumulative * scale_r),
              color = "grey50", linewidth = 0.6, linetype = "dotted") +
    geom_point(aes(y = cumulative * scale_r),
               color = "grey50", size = 1.8) +
    # Cumulative % labels on points
    geom_text(aes(y = cumulative * scale_r + max_eig * 0.03,
                  label = paste0(round(cumulative * 100, 1), "%")),
              size = 2.4, color = "grey50", vjust = 0) +
    scale_x_continuous(breaks = seq_len(n_pc),
                       labels = paste0("PC", seq_len(n_pc))) +
    scale_y_continuous(
      name   = "Eigenvalue",
      limits = c(0, max_eig),
      sec.axis = sec_axis(
        transform = ~ . / scale_r,
        name      = "Cumulative variance explained",
        labels    = scales::percent_format(accuracy = 1)
      )
    ) +
    labs(x = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text          = element_text(size = 8, colour = "grey50"),
      axis.title.y.left  = element_text(size = 8, colour = "grey40"),
      axis.title.y.right = element_text(size = 8, colour = "grey50"),
      plot.title         = element_text(size = 9, face = "plain")
    )
}

# ── 3a. PCA — raw photometric bands ──────────────────────────────────────────

# Full dataset
pca_bands_full    <- run_pca(df, photo_cols)
var_bands_full    <- pca_variance(pca_bands_full)

# Outlier-removed dataset (within-class MCD flags)
pca_bands_clean   <- run_pca(df_bands_no_out, photo_cols)
var_bands_clean   <- pca_variance(pca_bands_clean)

p_scree_bands_full  <- make_scree(var_bands_full,
                                  paste0("Raw bands — full (n=", nrow(df), ")"),
                                  class_colors["GALAXY"])
p_scree_bands_clean <- make_scree(var_bands_clean,
                                  paste0("Raw bands — outliers removed (n=",
                                         nrow(df_bands_no_out), ")"),
                                  class_colors["GALAXY"])

# ── 3b. PCA — adjacent color indices ─────────────────────────────────────────

pca_ci_full    <- run_pca(df_ci, ci_cols)
var_ci_full    <- pca_variance(pca_ci_full)

pca_ci_clean   <- run_pca(df_ci_no_out, ci_cols)
var_ci_clean   <- pca_variance(pca_ci_clean)

p_scree_ci_full  <- make_scree(var_ci_full,
                               paste0("Color indices — full (n=", nrow(df_ci), ")"),
                               class_colors["QSO"])
p_scree_ci_clean <- make_scree(var_ci_clean,
                               paste0("Color indices — outliers removed (n=",
                                      nrow(df_ci_no_out), ")"),
                               class_colors["QSO"])

# ── 3c. Combined scree plot (2x2) ─────────────────────────────────────────────
(p_scree_bands_full | p_scree_bands_clean) /
  (p_scree_ci_full  | p_scree_ci_clean) +
  plot_annotation(
    title   = "Scree plots: PCA on raw bands and color indices",
    caption = paste0("Bars = individual variance explained  |  ",
                     "Dashed line = cumulative  |  ",
                     "Data centered and scaled (z-score) before PCA"),
    theme   = theme(plot.title   = element_text(size = 12),
                    plot.caption = element_text(size = 7.5, colour = "grey50"))
  )

# ── 3d. Print variance tables ─────────────────────────────────────────────────
cat("\n--- PCA variance explained: raw bands ---\n")
cat("Full dataset:\n")
print(var_bands_full  %>% mutate(across(prop:cumulative, ~ round(., 3))))
cat("Outliers removed:\n")
print(var_bands_clean %>% mutate(across(prop:cumulative, ~ round(., 3))))

cat("\n--- PCA variance explained: color indices ---\n")
cat("Full dataset:\n")
print(var_ci_full  %>% mutate(across(prop:cumulative, ~ round(., 3))))
cat("Outliers removed:\n")
print(var_ci_clean %>% mutate(across(prop:cumulative, ~ round(., 3))))

# ── 3e. Score biplots — PC1 vs PC2 ───────────────────────────────────────────

# Helper: extract scores + loadings and build biplot
make_biplot <- function(pca_obj, class_vec, var_df, title) {
  
  n_pc  <- ncol(pca_obj$x)
  p_var <- ncol(pca_obj$rotation)
  
  # PC scores
  scores <- as_tibble(pca_obj$x[, 1:2]) %>%
    mutate(class = class_vec)
  
  # Loadings scaled for display (scale to ~20% of score range)
  loads <- as_tibble(pca_obj$rotation[, 1:2]) %>%
    mutate(variable = rownames(pca_obj$rotation))
  
  score_range <- max(abs(c(scores$PC1, scores$PC2)))
  load_scale  <- score_range * 0.45 / max(abs(c(loads$PC1, loads$PC2)))
  
  loads <- loads %>%
    mutate(PC1s = PC1 * load_scale,
           PC2s = PC2 * load_scale)
  
  # Axis labels with variance explained
  xlab <- paste0("PC1 (", round(var_df$prop[1] * 100, 1), "%)")
  ylab <- paste0("PC2 (", round(var_df$prop[2] * 100, 1), "%)")
  
  ggplot() +
    # Zero lines
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
    # Score points — sample 2000 max per class to avoid overplotting
    geom_point(
      data = scores %>%
        group_by(class) %>%
        slice_sample(n = 2000, replace = FALSE) %>%
        ungroup(),
      aes(x = PC1, y = PC2, color = class),
      alpha = 0.25, size = 0.7
    ) +
    # Class centroids
    geom_point(
      data = scores %>%
        group_by(class) %>%
        summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop"),
      aes(x = PC1, y = PC2, color = class),
      size = 3.5, shape = 18
    ) +
    # Loading arrows
    geom_segment(
      data = loads,
      aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      linewidth = 0.6, color = "grey20"
    ) +
    # Loading labels
    geom_text(
      data = loads,
      aes(x = PC1s * 1.12, y = PC2s * 1.12, label = variable),
      size = 2.8, color = "grey20", fontface = "bold"
    ) +
    scale_color_manual(values = class_colors) +
    labs(x = xlab, y = ylab, title = title, color = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position  = "top",
      legend.key.size  = unit(0.4, "cm"),
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93"),
      axis.text        = element_text(size = 8, colour = "grey50"),
      plot.title       = element_text(size = 9, face = "plain")
    )
}

# Raw bands — full
p_bi_bands_full  <- make_biplot(pca_bands_full,
                                df$class,
                                var_bands_full,
                                paste0("Raw bands — full (n=", nrow(df), ")"))

# Raw bands — outliers removed
p_bi_bands_clean <- make_biplot(pca_bands_clean,
                                df_bands_no_out$class,
                                var_bands_clean,
                                paste0("Raw bands — outliers removed (n=",
                                       nrow(df_bands_no_out), ")"))

# Color indices — full
p_bi_ci_full     <- make_biplot(pca_ci_full,
                                df_ci$class,
                                var_ci_full,
                                paste0("Color indices — full (n=", nrow(df_ci), ")"))

# Color indices — outliers removed
p_bi_ci_clean    <- make_biplot(pca_ci_clean,
                                df_ci_no_out$class,
                                var_ci_clean,
                                paste0("Color indices — outliers removed (n=",
                                       nrow(df_ci_no_out), ")"))

# ── 3f. Combined biplot (2x2) — coloured by class ────────────────────────────
(p_bi_bands_full | p_bi_bands_clean) /
  (p_bi_ci_full  | p_bi_ci_clean) +
  plot_layout(guides = "collect") +
  theme(legend.position = "top") +
  plot_annotation(
    title   = "PCA score biplots: PC1 vs PC2 by class",
    caption = paste0("Points = individual observations (max 2,000 per class shown)  |  ",
                     "Diamond = class centroid  |  ",
                     "Arrows = variable loadings (scaled)"),
    theme   = theme(plot.title   = element_text(size = 12),
                    plot.caption = element_text(size = 7.5, colour = "grey50"))
  )

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — ASSUMPTION TESTING: MULTIVARIATE NORMALITY (MARDIA'S TEST)
# ═══════════════════════════════════════════════════════════════════════════════
# Required for FA (assumes MVN of common factors) and to assess whether PCA
# covariance-based inference is appropriate.
# install.packages("MVN") if not yet installed
library(MVN)

# ── Helper: run Mardia's test and return a tidy summary tibble ────────────────
run_mardia <- function(data, cols, label) {
  X   <- as.matrix(data[, cols])
  res <- mvn(X, mvnTest = "mardia", desc = FALSE)$multivariateNormality
  
  tibble(
    dataset   = label,
    test      = res$Test,
    statistic = round(as.numeric(res$`Test Statistic`), 3),
    p_value   = round(as.numeric(res$`p value`), 4),
    result    = res$Result
  )
}

# ── 4a. Mardia's test — all four dataset/feature combinations ─────────────────
mardia_results <- bind_rows(
  run_mardia(df,              photo_cols, "Raw bands — full"),
  run_mardia(df_bands_no_out, photo_cols, "Raw bands — outliers removed"),
  run_mardia(df_ci,           ci_cols,    "Color indices — full"),
  run_mardia(df_ci_no_out,    ci_cols,    "Color indices — outliers removed")
)

cat("\n--- Mardia's MVN test results ---\n")
print(mardia_results, n = Inf)

# ── 4b. Mardia's test — within each class (outliers removed) ─────────────────
cat("\n--- Mardia's test by class: raw bands (outliers removed) ---\n")
mardia_bands_class <- bind_rows(lapply(class_order, function(cls) {
  run_mardia(filter(df_bands_no_out, class == cls), photo_cols,
             paste0("Raw bands — ", cls))
}))
print(mardia_bands_class, n = Inf)

cat("\n--- Mardia's test by class: color indices (outliers removed) ---\n")
mardia_ci_class <- bind_rows(lapply(class_order, function(cls) {
  run_mardia(filter(df_ci_no_out, class == cls), ci_cols,
             paste0("Color indices — ", cls))
}))
print(mardia_ci_class, n = Inf)

# ── 4c. Chi-squared QQ plots — visual MVN check ───────────────────────────────
# Squared Mahalanobis distances should follow chi-sq(p) under MVN.
# Deviations from the diagonal indicate heavy tails or skewness.
make_mvn_qq <- function(data, cols, label) {
  X   <- scale(as.matrix(data[, cols]))
  p   <- ncol(X)
  n   <- nrow(X)
  d2  <- mahalanobis(X, center = rep(0, p), cov = cov(X))
  probs <- (seq_len(n) - 0.5) / n
  
  tibble(
    empirical   = sort(d2),
    theoretical = qchisq(probs, df = p),
    label       = label
  )
}

qq_data <- bind_rows(
  make_mvn_qq(df,              photo_cols, "Raw bands — full"),
  make_mvn_qq(df_bands_no_out, photo_cols, "Raw bands — outliers removed"),
  make_mvn_qq(df_ci,           ci_cols,    "Color indices — full"),
  make_mvn_qq(df_ci_no_out,    ci_cols,    "Color indices — outliers removed")
) %>%
  mutate(label = factor(label, levels = c(
    "Raw bands — full",
    "Raw bands — outliers removed",
    "Color indices — full",
    "Color indices — outliers removed"
  )))

ggplot(qq_data, aes(x = theoretical, y = empirical)) +
  geom_point(alpha = 0.3, size = 0.6, color = "grey40") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", linewidth = 0.5, color = "#B2182B") +
  facet_wrap(~ label, ncol = 2, scales = "free") +
  labs(
    x       = "Theoretical chi-squared quantiles",
    y       = "Squared Mahalanobis distances",
    title   = "Chi-squared QQ plots: multivariate normality check",
    caption = paste0("Points on the red line indicate MVN  |  ",
                     "Upward deviation in tails = heavy-tailed distribution")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text       = element_text(size = 9),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(size = 8, colour = "grey50"),
    plot.caption     = element_text(size = 7.5, colour = "grey50"),
    plot.title       = element_text(size = 11)
  )