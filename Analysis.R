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

# ── 1d. D-D plots: raw bands — within-class MCD ──────────────────────────────
dd_bands_typed <- map_dfr(class_order, function(cls) {
  idx  <- which(df$class == cls)
  Xsub <- X_bands[idx, ]
  
  d_cl <- sqrt(mahalanobis(Xsub, center = colMeans(Xsub), cov = cov(Xsub)))
  
  mcd  <- covMcd(Xsub, alpha = 0.95)
  d_ro <- sqrt(mahalanobis(Xsub, center = mcd$center, cov = mcd$cov))
  
  tibble(classical = d_cl, robust = d_ro, class = cls, row_idx = idx)
}) %>%
  get_outlier_flags(cutoff_bands)

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

run_pca <- function(data, cols) {
  prcomp(data[, cols], center = TRUE, scale. = TRUE)
}

pca_variance <- function(pca_obj) {
  var_exp <- pca_obj$sdev^2
  tibble(
    pc         = seq_along(var_exp),
    variance   = var_exp,
    prop       = var_exp / sum(var_exp),
    cumulative = cumsum(var_exp / sum(var_exp))
  )
}

make_scree <- function(var_df, title, color) {
  n_pc    <- nrow(var_df)
  max_eig <- max(var_df$variance) * 1.2
  scale_r <- max_eig
  
  ggplot(var_df, aes(x = pc)) +
    geom_col(aes(y = variance), fill = color, alpha = 0.75, width = 0.5) +
    geom_text(aes(y = variance - max_eig * 0.05, label = round(variance, 2)),
              size = 2.6, vjust = 0, color = "grey30") +
    geom_hline(yintercept = 1, linetype = "dashed",
               linewidth = 0.5, color = "#B2182B") +
    annotate("text", x = n_pc + 0.45, y = 1,
             label = "Kaiser Criterion", hjust = 1, vjust = -0.3,
             size = 2.4, color = "#B2182B") +
    geom_line(aes(y = cumulative * scale_r),
              color = "grey50", linewidth = 0.6, linetype = "dotted") +
    geom_point(aes(y = cumulative * scale_r), color = "grey50", size = 1.8) +
    geom_text(aes(y = cumulative * scale_r + max_eig * 0.03,
                  label = paste0(round(cumulative * 100, 1), "%")),
              size = 2.4, color = "grey50", vjust = 0) +
    scale_x_continuous(breaks = seq_len(n_pc),
                       labels = paste0("PC", seq_len(n_pc))) +
    scale_y_continuous(
      name     = "Eigenvalue",
      limits   = c(0, max_eig),
      sec.axis = sec_axis(transform = ~ . / scale_r,
                          name      = "Cumulative variance explained",
                          labels    = scales::percent_format(accuracy = 1))
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
pca_bands_full  <- run_pca(df,             photo_cols)
var_bands_full  <- pca_variance(pca_bands_full)
pca_bands_clean <- run_pca(df_bands_no_out, photo_cols)
var_bands_clean <- pca_variance(pca_bands_clean)

p_scree_bands_full  <- make_scree(var_bands_full,
                                  paste0("Raw bands — full (n=", nrow(df), ")"),
                                  class_colors["GALAXY"])
p_scree_bands_clean <- make_scree(var_bands_clean,
                                  paste0("Raw bands — outliers removed (n=",
                                         nrow(df_bands_no_out), ")"),
                                  class_colors["GALAXY"])

# ── 3b. PCA — adjacent color indices ─────────────────────────────────────────
pca_ci_full  <- run_pca(df_ci,        ci_cols)
var_ci_full  <- pca_variance(pca_ci_full)
pca_ci_clean <- run_pca(df_ci_no_out, ci_cols)
var_ci_clean <- pca_variance(pca_ci_clean)

p_scree_ci_full  <- make_scree(var_ci_full,
                               paste0("Color indices — full (n=", nrow(df_ci), ")"),
                               class_colors["QSO"])
p_scree_ci_clean <- make_scree(var_ci_clean,
                               paste0("Color indices — outliers removed (n=",
                                      nrow(df_ci_no_out), ")"),
                               class_colors["QSO"])

# ── 3c. Combined scree plot (2x2) ────────────────────────────────────────────
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

# ── 3d. Print variance tables ────────────────────────────────────────────────
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
make_biplot <- function(pca_obj, class_vec, var_df, title) {
  scores <- as_tibble(pca_obj$x[, 1:2]) %>% mutate(class = class_vec)
  loads  <- as_tibble(pca_obj$rotation[, 1:2]) %>%
    mutate(variable = rownames(pca_obj$rotation))
  
  score_range <- max(abs(c(scores$PC1, scores$PC2)))
  load_scale  <- score_range * 0.45 / max(abs(c(loads$PC1, loads$PC2)))
  loads <- loads %>% mutate(PC1s = PC1 * load_scale, PC2s = PC2 * load_scale)
  
  xlab <- paste0("PC1 (", round(var_df$prop[1] * 100, 1), "%)")
  ylab <- paste0("PC2 (", round(var_df$prop[2] * 100, 1), "%)")
  
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey80") +
    geom_point(
      data = scores %>% group_by(class) %>%
        slice_sample(n = 2000, replace = FALSE) %>% ungroup(),
      aes(x = PC1, y = PC2, color = class), alpha = 0.25, size = 0.7
    ) +
    geom_point(
      data = scores %>% group_by(class) %>%
        summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop"),
      aes(x = PC1, y = PC2, color = class), size = 3.5, shape = 18
    ) +
    geom_segment(data = loads,
                 aes(x = 0, y = 0, xend = PC1s, yend = PC2s),
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                 linewidth = 0.6, color = "grey20") +
    geom_text(data = loads,
              aes(x = PC1s * 1.12, y = PC2s * 1.12, label = variable),
              size = 2.8, color = "grey20", fontface = "bold") +
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

p_bi_bands_full  <- make_biplot(pca_bands_full,    df$class,
                                var_bands_full,
                                paste0("Raw bands — full (n=", nrow(df), ")"))
p_bi_bands_clean <- make_biplot(pca_bands_clean,   df_bands_no_out$class,
                                var_bands_clean,
                                paste0("Raw bands — outliers removed (n=",
                                       nrow(df_bands_no_out), ")"))
p_bi_ci_full     <- make_biplot(pca_ci_full,       df_ci$class,
                                var_ci_full,
                                paste0("Color indices — full (n=", nrow(df_ci), ")"))
p_bi_ci_clean    <- make_biplot(pca_ci_clean,      df_ci_no_out$class,
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
# SECTION 4 — ASSUMPTION TESTING FOR FA
# ═══════════════════════════════════════════════════════════════════════════════
# install.packages(c("MVN", "psych")) if not yet installed
# library(MVN)
library(psych)

# ── 4a. Mardia's test — MVN check ────────────────────────────────────────────
# Mardia's skewness is O(n^2) — subsample to 1000 obs to keep it tractable.
# set.seed(437)
# N_MARDIA <- 1000
# 
# cat("\n--- Raw bands: full dataset ---\n")
# mvn(as.matrix(df[sample(nrow(df), N_MARDIA), photo_cols]),
#     mvnTest = "mardia", desc = FALSE)
# 
# cat("\n--- Raw bands: outliers removed ---\n")
# mvn(as.matrix(df_bands_no_out[sample(nrow(df_bands_no_out), N_MARDIA), photo_cols]),
#     mvnTest = "mardia", desc = FALSE)
# 
# cat("\n--- Color indices: full dataset ---\n")
# mvn(as.matrix(df_ci[sample(nrow(df_ci), N_MARDIA), ci_cols]),
#     mvnTest = "mardia", desc = FALSE)
# 
# cat("\n--- Color indices: outliers removed ---\n")
# mvn(as.matrix(df_ci_no_out[sample(nrow(df_ci_no_out), N_MARDIA), ci_cols]),
#     mvnTest = "mardia", desc = FALSE)

# ── 4b. KMO — factorability check ────────────────────────────────────────────
# KMO measures sampling adequacy for FA. Values >= 0.80 are "meritorious",
# >= 0.90 are "marvelous". Low per-variable MSA flags variables that share
# little common variance with others and may not belong in a factor model.

cat("\n--- KMO: Raw bands — full ---\n")
KMO(cor(df[, photo_cols]))

cat("\n--- KMO: Raw bands — outliers removed ---\n")
KMO(cor(df_bands_no_out[, photo_cols]))

cat("\n--- KMO: Color indices — full ---\n")
KMO(cor(df_ci[, ci_cols]))

cat("\n--- KMO: Color indices — outliers removed ---\n")
KMO(cor(df_ci_no_out[, ci_cols]))

# ── 4c. Negentropy — effect size of non-Gaussianity per variable ─────────────
# Not inflated by large n — complements Mardia's test.
# Approximation: Hyvarinen (1998): J(x) ~ [E{G(x)} - E{G(v)}]^2, G(u) = log(cosh(u))

set.seed(437)
v   <- rnorm(100000)
EGv <- mean(log(cosh(v)))

negentropy_vars <- function(data, cols, label) {
  X <- scale(as.matrix(data[, cols]))
  neg <- sapply(cols, function(col) {
    s   <- X[, col]
    EGs <- mean(log(cosh(s)))
    round((EGs - EGv)^2, 6)
  })
  tibble(dataset = label, variable = cols, negentropy = neg)
}

neg_results <- bind_rows(
  negentropy_vars(df,              photo_cols, "Raw bands — full"),
  negentropy_vars(df_bands_no_out, photo_cols, "Raw bands — outliers removed"),
  negentropy_vars(df_ci,           ci_cols,    "Color indices — full"),
  negentropy_vars(df_ci_no_out,    ci_cols,    "Color indices — outliers removed")
)

cat("\n--- Negentropy per variable (0 = Gaussian, higher = more non-Gaussian) ---\n")
print(neg_results %>% pivot_wider(names_from = variable, values_from = negentropy),
      n = Inf)

neg_class <- bind_rows(
  negentropy_vars(filter(df_bands_no_out, class == "GALAXY"), photo_cols, "Raw bands — GALAXY"),
  negentropy_vars(filter(df_bands_no_out, class == "STAR"),   photo_cols, "Raw bands — STAR"),
  negentropy_vars(filter(df_bands_no_out, class == "QSO"),    photo_cols, "Raw bands — QSO"),
  negentropy_vars(filter(df_ci_no_out,    class == "GALAXY"), ci_cols,    "Color indices — GALAXY"),
  negentropy_vars(filter(df_ci_no_out,    class == "STAR"),   ci_cols,    "Color indices — STAR"),
  negentropy_vars(filter(df_ci_no_out,    class == "QSO"),    ci_cols,    "Color indices — QSO")
)

cat("\n--- Negentropy per variable by class ---\n")
print(neg_class %>% pivot_wider(names_from = variable, values_from = negentropy),
      n = Inf)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — FACTOR ANALYSIS (ML, 2 factors)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Helper: FA scree plot with parallel analysis ──────────────────────────────
make_fa_scree <- function(data, cols, title, color, n_sim = 500) {
  X      <- scale(as.matrix(data[, cols]))
  n      <- nrow(X)
  p      <- ncol(X)
  R      <- cor(X)
  eigval <- eigen(R)$values
  
  # Parallel analysis: simulate random correlation matrices n_sim times
  set.seed(437)
  sim_eigvals <- replicate(n_sim, {
    Xr <- matrix(rnorm(n * p), nrow = n, ncol = p)
    eigen(cor(Xr))$values
  })
  # 95th percentile of simulated eigenvalues per factor (conservative criterion)
  pa_95 <- apply(sim_eigvals, 1, quantile, probs = 0.95)
  # Mean of simulated eigenvalues
  pa_mean <- rowMeans(sim_eigvals)
  
  var_df <- tibble(
    factor   = seq_len(p),
    observed = eigval,
    pa_95    = pa_95,
    pa_mean  = pa_mean
  )
  
  ggplot(var_df, aes(x = factor)) +
    # Parallel analysis 95th percentile line
    geom_line(aes(y = pa_95), color = "grey55", linewidth = 0.6,
              linetype = "longdash") +
    geom_point(aes(y = pa_95), color = "grey55", size = 1.8, shape = 1) +
    # Observed eigenvalue line
    geom_line(aes(y = observed), color = color, linewidth = 0.7) +
    geom_point(aes(y = observed), color = color, size = 2.5) +
    # Kaiser criterion
    geom_hline(yintercept = 1, linetype = "dashed",
               linewidth = 0.5, color = "#B2182B") +
    annotate("text", x = p, y = 1,
             label = "Kaiser Criterion", hjust = 1, vjust = -0.4,
             size = 2.3, color = "#B2182B") +
    scale_x_continuous(breaks = seq_len(p),
                       labels = paste0("F", seq_len(p))) +
    labs(x = NULL, y = "Eigenvalue", title = title) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text          = element_text(size = 8, colour = "grey50"),
      plot.title         = element_text(size = 9, face = "plain")
    )
}

# ── 5a. FA scree plots — all four combinations ───────────────────────────────
p_fa_scree_bands_full  <- make_fa_scree(df,              photo_cols,
                                        paste0("Raw bands — full (n=", nrow(df), ")"),
                                        class_colors["GALAXY"])
p_fa_scree_bands_clean <- make_fa_scree(df_bands_no_out, photo_cols,
                                        paste0("Raw bands — outliers removed (n=", nrow(df_bands_no_out), ")"),
                                        class_colors["GALAXY"])
p_fa_scree_ci_full     <- make_fa_scree(df_ci,           ci_cols,
                                        paste0("Color indices — full (n=", nrow(df_ci), ")"),
                                        class_colors["QSO"])
p_fa_scree_ci_clean    <- make_fa_scree(df_ci_no_out,    ci_cols,
                                        paste0("Color indices — outliers removed (n=", nrow(df_ci_no_out), ")"),
                                        class_colors["QSO"])

(p_fa_scree_bands_full | p_fa_scree_bands_clean) /
  (p_fa_scree_ci_full  | p_fa_scree_ci_clean) +
  plot_annotation(
    title   = "FA scree plots with parallel analysis",
    caption = paste0("Solid line = observed eigenvalues  |  ",
                     "Dashed line = parallel analysis 95th percentile (n=500 simulations)  |  ",
                     "Red dashed = Kaiser criterion (eigenvalue > 1)"),
    theme   = theme(plot.title   = element_text(size = 12),
                    plot.caption = element_text(size = 7.5, colour = "grey50"))
  )

# ── Helper: run FA for all three rotations and extract loadings ───────────────
run_fa_loadings <- function(data, cols, n_factors, dataset_label) {
  X <- scale(as.matrix(data[, cols]))
  
  rotations <- c("none", "varimax", "promax")
  rot_labels <- c("No Rotation", "Varimax", "Promax")
  
  map2_dfr(rotations, rot_labels, function(rot, rot_lab) {
    fa_obj <- fa(X, nfactors = n_factors, rotate = rot, fm = "ml", scores = "regression")
    loads  <- as.data.frame(unclass(fa_obj$loadings))
    loads$variable <- rownames(loads)
    loads %>%
      pivot_longer(-variable, names_to = "factor", values_to = "loading") %>%
      mutate(
        rotation = rot_lab,
        dataset  = dataset_label,
        factor   = recode(factor, MR1 = "ML1", MR2 = "ML2")
      )
  })
}

N_FACTORS <- 2

fa_bands_full  <- run_fa_loadings(df,              photo_cols, N_FACTORS, "Raw bands — full")
fa_bands_clean <- run_fa_loadings(df_bands_no_out, photo_cols, N_FACTORS, "Raw bands — outliers removed")
fa_ci_full     <- run_fa_loadings(df_ci,           ci_cols,    N_FACTORS, "Color indices — full")
fa_ci_clean    <- run_fa_loadings(df_ci_no_out,    ci_cols,    N_FACTORS, "Color indices — outliers removed")

# ── Helper: loading heatmap like the reference image ─────────────────────────
make_fa_heatmap <- function(fa_df, title) {
  fa_df <- fa_df %>%
    mutate(
      rotation = factor(rotation, levels = c("No Rotation", "Varimax", "Promax")),
      variable = factor(variable, levels = rev(unique(variable)))
    )
  
  ggplot(fa_df, aes(x = factor, y = variable, fill = loading)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", loading),
                  color = abs(loading) > 0.5),
              size = 2.8, fontface = "plain") +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "#F7F7F7",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Loading"
    ) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey30"),
                       guide  = "none") +
    facet_wrap(~ rotation, nrow = 1) +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text       = element_text(size = 9, face = "plain"),
      panel.grid       = element_blank(),
      axis.text.x      = element_text(size = 8, colour = "grey40"),
      axis.text.y      = element_text(size = 8, colour = "grey40"),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.text       = element_text(size = 7),
      legend.title      = element_text(size = 8),
      plot.title        = element_text(size = 9, face = "plain")
    )
}

# ── 5b. Loading heatmaps — all four combinations ─────────────────────────────
p_fa_bands_full  <- make_fa_heatmap(fa_bands_full,  "Raw bands — full")
p_fa_bands_clean <- make_fa_heatmap(fa_bands_clean, "Raw bands — outliers removed")
p_fa_ci_full     <- make_fa_heatmap(fa_ci_full,     "Color indices — full")
p_fa_ci_clean    <- make_fa_heatmap(fa_ci_clean,    "Color indices — outliers removed")

p_fa_bands_full
p_fa_bands_clean
p_fa_ci_full
p_fa_ci_clean

# ── 5c. Print uniquenesses and communalities for each combination ─────────────
print_uniquenesses <- function(data, cols, label) {
  X      <- scale(as.matrix(data[, cols]))
  fa_obj <- fa(X, nfactors = N_FACTORS, rotate = "varimax", fm = "ml")
  result <- tibble(
    variable    = cols,
    communality = round(fa_obj$communality, 3),
    uniqueness  = round(fa_obj$uniquenesses, 3)
  )
  cat(sprintf("\n--- Communalities & Uniquenesses: %s ---\n", label))
  print(result, n = Inf)
}

print_uniquenesses(df,              photo_cols, "Raw bands — full")
print_uniquenesses(df_bands_no_out, photo_cols, "Raw bands — outliers removed")
print_uniquenesses(df_ci,           ci_cols,    "Color indices — full")
print_uniquenesses(df_ci_no_out,    ci_cols,    "Color indices — outliers removed")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — LDA & LOGISTIC REGRESSION WITH CROSS-VALIDATION
# ═══════════════════════════════════════════════════════════════════════════════
# install.packages("nnet") if not yet installed — for multinomial logistic regression
library(nnet)

set.seed(437)
K_FOLDS <- 10

# ── Helper: stratified k-fold indices ────────────────────────────────────────
make_folds_stratified <- function(y, k) {
  folds <- vector("list", k)
  for (cls in unique(y)) {
    idx    <- which(y == cls)
    idx    <- sample(idx)
    splits <- cut(seq_along(idx), breaks = k, labels = FALSE)
    for (i in seq_len(k)) folds[[i]] <- c(folds[[i]], idx[splits == i])
  }
  folds
}

# ── Helper: unstratified k-fold indices ──────────────────────────────────────
make_folds_unstratified <- function(y, k) {
  n      <- length(y)
  idx    <- sample(seq_len(n))
  splits <- cut(seq_len(n), breaks = k, labels = FALSE)
  map(seq_len(k), ~ idx[splits == .x])
}

# ── Helper: run CV for LDA or logistic and return per-fold metrics ────────────
run_cv <- function(data, cols, y_col, k, method, dataset_label,
                   stratified = TRUE) {
  df_model <- data %>%
    select(all_of(c(cols, y_col))) %>%
    mutate(across(all_of(cols), ~ scale(.)[,1]))
  
  y       <- df_model[[y_col]]
  folds   <- if (stratified) make_folds_stratified(y, k) else
    make_folds_unstratified(y, k)
  classes <- sort(unique(y))
  
  fold_results <- map(seq_len(k), function(i) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df_model)), test_idx)
    
    train <- df_model[train_idx, ]
    test  <- df_model[test_idx,  ]
    
    formula <- as.formula(paste(y_col, "~",
                                paste(paste0("`", cols, "`"), collapse = " + ")))
    
    if (method == "lda") {
      fit   <- lda(formula, data = train)
      preds <- predict(fit, newdata = test)$class
    } else {
      fit   <- multinom(formula, data = train, trace = FALSE)
      preds <- predict(fit, newdata = test)
    }
    
    # Overall accuracy
    acc <- mean(preds == test[[y_col]])
    
    # Per-class recall
    per_class <- map_dfr(classes, function(cls) {
      truth <- test[[y_col]] == cls
      pred  <- preds == cls
      tibble(
        class     = cls,
        recall    = ifelse(sum(truth) > 0, sum(truth & pred) / sum(truth), NA_real_),
        precision = ifelse(sum(pred)  > 0, sum(truth & pred) / sum(pred),  NA_real_)
      )
    })
    
    # Confusion matrix counts for this fold
    cm <- as_tibble(table(truth = test[[y_col]], predicted = preds)) %>%
      mutate(fold = i, method = method, dataset = dataset_label)
    
    list(
      metrics = tibble(fold = i, method = method, dataset = dataset_label,
                       accuracy = acc) %>%
        bind_cols(per_class %>%
                    pivot_wider(names_from = class,
                                values_from = c(recall, precision),
                                names_glue = "{.value}_{class}")),
      cm = cm
    )
  })
  
  list(
    metrics = map_dfr(fold_results, "metrics"),
    cm      = map_dfr(fold_results, "cm")
  )
}

# ── 6a. Run CV — stratified and unstratified ─────────────────────────────────
cv_datasets <- list(
  list(data = df,              cols = photo_cols, label = "Raw bands — full"),
  list(data = df_bands_no_out, cols = photo_cols, label = "Raw bands — outliers removed"),
  list(data = df_ci,           cols = ci_cols,    label = "Color indices — full"),
  list(data = df_ci_no_out,    cols = ci_cols,    label = "Color indices — outliers removed")
)

run_all_cv <- function(stratified) {
  map(cv_datasets, function(d) {
    list(
      lda      = run_cv(d$data, d$cols, "class", K_FOLDS, "lda",      d$label, stratified),
      logistic = run_cv(d$data, d$cols, "class", K_FOLDS, "logistic", d$label, stratified)
    )
  })
}

cv_raw_strat   <- run_all_cv(stratified = TRUE)
cv_raw_unstrat <- run_all_cv(stratified = FALSE)

extract_results <- function(cv_raw, sampling) {
  metrics <- map_dfr(cv_raw, function(x)
    bind_rows(x$lda$metrics, x$logistic$metrics)) %>%
    mutate(sampling = sampling)
  cm <- map_dfr(cv_raw, function(x)
    bind_rows(x$lda$cm, x$logistic$cm)) %>%
    mutate(sampling = sampling)
  list(metrics = metrics, cm = cm)
}

res_strat   <- extract_results(cv_raw_strat,   "Stratified")
res_unstrat <- extract_results(cv_raw_unstrat, "Unstratified")

cv_results <- bind_rows(res_strat$metrics, res_unstrat$metrics)
cv_cm_all  <- bind_rows(res_strat$cm,      res_unstrat$cm)

# ── 6b. Summarise CV metrics ──────────────────────────────────────────────────
cv_summary <- cv_results %>%
  group_by(method, dataset, sampling) %>%
  summarise(
    mean_acc = mean(accuracy),
    sd_acc   = sd(accuracy),
    across(starts_with("recall_"),    ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(starts_with("precision_"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    .groups = "drop"
  ) %>%
  mutate(
    method   = recode(method, lda = "LDA", logistic = "Logistic"),
    dataset  = factor(dataset, levels = map_chr(cv_datasets, "label")),
    sampling = factor(sampling, levels = c("Stratified", "Unstratified"))
  )

cat("\n--- CV accuracy summary ---\n")
cv_summary %>%
  select(sampling, method, dataset, mean_acc, sd_acc) %>%
  mutate(across(c(mean_acc, sd_acc), ~ round(., 4))) %>%
  arrange(dataset, method, sampling) %>%
  print(n = Inf)

# ── 6c. Bar chart comparing accuracy — stratified vs unstratified ─────────────
method_colors <- c(LDA = "#378ADD", Logistic = "#D85A30")

p_acc <- ggplot(cv_summary,
                aes(x = dataset, y = mean_acc, fill = method)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc),
                position = position_dodge(width = 0.7),
                width = 0.25, linewidth = 0.5, color = "grey30") +
  geom_text(aes(label = paste0(round(mean_acc * 100, 1), "%"),
                y     = mean_acc + sd_acc + 0.005),
            position = position_dodge(width = 0.7),
            size = 2.4, vjust = 0, color = "grey20") +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 16)) +
  facet_wrap(~ sampling, ncol = 1) +
  labs(
    x       = NULL,
    y       = "Mean accuracy (10-fold CV)",
    fill    = NULL,
    title   = "LDA vs. logistic: CV accuracy — stratified vs. unstratified folds",
    caption = paste0("Error bars = ±1 SD across folds  |  ",
                     "Features standardized before fitting  |  ",
                     "Stratified folds preserve class imbalance (50% / 42% / 8%)")
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position    = "top",
    strip.text         = element_text(size = 9, face = "plain"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(size = 8, colour = "grey40"),
    axis.text.y        = element_text(size = 8, colour = "grey50"),
    plot.caption       = element_text(size = 7.5, colour = "grey50"),
    plot.title         = element_text(size = 11)
  )

p_acc

# ── 6d. Confusion matrices — averaged across folds ───────────────────────────
# Normalise by true class (row) so cells show recall per class pair
make_cm_plots <- function(method_label, sampling_label = "Stratified") {
  cm_data <- cv_cm_all %>%
    filter(method == tolower(method_label),
           sampling == sampling_label) %>%
    group_by(dataset, truth, predicted) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(dataset, truth) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      dataset   = factor(dataset, levels = map_chr(cv_datasets, "label")),
      truth     = factor(truth,     levels = class_order),
      predicted = factor(predicted, levels = class_order),
      label     = paste0(round(pct * 100, 1), "%\n(n=", n, ")")
    )
  
  ggplot(cm_data, aes(x = predicted, y = truth, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label,
                  color = pct > 0.55),
              size = 2.5, lineheight = 1.2) +
    scale_fill_gradient2(low = "#F7F7F7", mid = "#7FBCD2", high = "#2166AC",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Recall") +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey30"),
                       guide = "none") +
    facet_wrap(~ dataset, ncol = 2,
               labeller = labeller(dataset = function(x) str_wrap(x, 22))) +
    labs(
      x     = "Predicted class",
      y     = "True class",
      title = paste0("Confusion matrices: ", method_label,
                     " — ", sampling_label,
                     " (10-fold CV, normalised by true class)")
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text    = element_text(size = 8),
      panel.grid    = element_blank(),
      axis.text     = element_text(size = 8, colour = "grey40"),
      legend.key.height = unit(0.6, "cm"),
      plot.title    = element_text(size = 11)
    )
}

# Stratified confusion matrices
make_cm_plots("lda",      "Stratified")
make_cm_plots("logistic", "Stratified")

# Unstratified confusion matrices
make_cm_plots("lda",      "Unstratified")
make_cm_plots("logistic", "Unstratified")

# ── AIC/BIC factor selection ──────────────────────────────────────────────────
fa_ic <- function(data, cols, label, max_factors = 4) {
  X <- scale(as.matrix(data[, cols]))
  n <- nrow(X)
  
  map_dfr(seq_len(max_factors), function(nf) {
    fit <- tryCatch(
      fa(X, nfactors = nf, rotate = "none", fm = "ml"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    
    # psych stores objective (minimised -2*loglik / n), recover loglik
    # BIC and SABIC are available directly in fit
    tibble(
      dataset    = label,
      n_factors  = nf,
      BIC        = fit$BIC,
      SABIC      = fit$SABIC,   # sample-size adjusted BIC — most useful here
      RMSEA      = fit$RMSEA[1],
      TLI        = fit$TLI
    )
  })
}

ic_results <- bind_rows(
  fa_ic(df,              photo_cols, "Raw bands — full"),
  fa_ic(df_bands_no_out, photo_cols, "Raw bands — outliers removed"),
  fa_ic(df_ci,           ci_cols,    "Color indices — full"),
  fa_ic(df_ci_no_out,    ci_cols,    "Color indices — outliers removed")
)

cat("\n--- FA information criteria ---\n")
print(ic_results, n = Inf)