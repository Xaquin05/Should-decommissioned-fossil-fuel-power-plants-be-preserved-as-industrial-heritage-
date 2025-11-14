# ==========================================================
# 1) Load packages & global settings
# 2) Weighted-mean imputation helper (subset-level)
# 3) 3-panel PCA biplots (Spain | Asturias+León+Palencia | As Pontes)
# 4) Robustness: 6-panel PCA biplots (Coruña, Asturias, Teruel + clusters)
# 5) Descriptives (means±SD, medians) for 8 samples
# 6) Eigenvalues table (incl. FactorsExtracted) + PNG
# 7) Scree plots (matching panels)
# 9) Communalities

# ==========================================================
# ---- 1) Packages & settings ----
# ==========================================================

pkgs <- c("haven","dplyr","psych","ggplot2","patchwork","ggrepel",
          "stringi","tidyr","tibble","readr","gridExtra","grid")
to_install <- setdiff(pkgs, installed.packages()[,"Package"])
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, require, character.only = TRUE))

file_sav <- "Dataset_Nature_heritage.sav"
vars     <- c("Reindustrialization","Renaturalization","Preservation")
items    <- vars
dist_thr <- 150
mining_provinces <- c("Asturias","León","Leon","Palencia")

# Visual styles
base_size_small <- 18 * 0.5
axis_title_small<- 16 * 0.5
axis_text_small <- 14 * 0.5
label_size_small<- 5.3 * 0.5
arrow_len_small <- 0.12

base_size_big   <- 16
axis_title_big  <- 15
axis_text_big   <- 13
label_size_big  <- 6.3
arrow_len_big   <- 0.22

# Turn imputation on/off (applies everywhere)
impute_missing <- TRUE

# =============================================================================
# ---- 2) Weighted-mean imputation helper (subset-level) ----
# =============================================================================
impute_weighted_means <- function(df, vars, weight_col = NULL, clamp01 = TRUE, lo = 0, hi = 10) {
  d <- df
  for (v in vars) {
    x <- as.numeric(d[[v]])
    if (!is.null(weight_col) && weight_col %in% names(d)) {
      w <- as.numeric(d[[weight_col]])
      ok <- is.finite(x) & is.finite(w)
      m  <- if (any(ok)) sum(w[ok] * x[ok]) / sum(w[ok]) else mean(x, na.rm = TRUE)
    } else {
      m  <- mean(x, na.rm = TRUE)
    }
    if (!is.finite(m)) m <- 0
    x[is.na(x)] <- m
    if (clamp01) x <- pmax(lo, pmin(hi, x))
    d[[v]] <- x
  }
  d
}

# ---- Data load once; also normalized copies for name-robust filters ----
dat <- haven::read_sav(file_sav)
dat_norm <- dat %>%
  dplyr::mutate(
    Province_norm = tolower(stringi::stri_trans_general(Province, "Latin-ASCII")),
    Cluster_norm  = tolower(stringi::stri_trans_general(PowerPlantCluster, "Latin-ASCII")),
    Dist_to_nearest_plant = suppressWarnings(as.numeric(Dist_to_nearest_plant))
  )

# =============================================================================
# 3) Three-panel PCA biplots (small fonts), with imputation
# =============================================================================
biplot_from_df <- function(df, weight_col = NULL, panel_label = "(a)",
                           base_size = base_size_small, axis_title = axis_title_small,
                           axis_text = axis_text_small, label_size = label_size_small,
                           arrow_len = arrow_len_small) {
  need <- unique(c(vars, weight_col))
  d <- df %>%
    dplyr::select(any_of(need)) %>%
    dplyr::mutate(dplyr::across(all_of(vars), as.numeric))
  
  # drop NA weights if provided
  w <- NULL
  if (!is.null(weight_col) && weight_col %in% names(d)) {
    d <- d %>% dplyr::filter(!is.na(.data[[weight_col]]))
    w <- as.numeric(d[[weight_col]])
  }
  
  # impute within subset; else listwise
  if (isTRUE(impute_missing)) {
    d <- impute_weighted_means(d, vars, weight_col)
  } else {
    d <- d %>% dplyr::filter(dplyr::if_all(all_of(vars), ~ !is.na(.x)))
  }
  
  X <- as.matrix(d[, vars, drop = FALSE])
  
  # Weighted correlation matrix
  R <- stats::cov.wt(X, wt = w, method = "unbiased", cor = TRUE)$cor
  
  # PCA: 2 comps, Varimax (Kaiser)
  pc  <- psych::principal(r = R, nfactors = 2, rotate = "varimax")
  L   <- as.matrix(unclass(pc$loadings))
  pct <- 100 * colSums(L^2) / nrow(L)
  
  # Plot data + keep labels inside
  t <- seq(0, 2*pi, length.out = 400)
  circ <- data.frame(x = cos(t), y = sin(t))
  load_df <- data.frame(Variable = rownames(L), Dim1 = L[,1], Dim2 = L[,2])
  radial_inset <- 0.92
  r <- sqrt(load_df$Dim1^2 + load_df$Dim2^2)
  s <- ifelse(r > radial_inset, radial_inset / pmax(r, 1e-9), 1)
  load_df$TX <- load_df$Dim1 * s
  load_df$TY <- load_df$Dim2 * s
  
  ggplot() +
    geom_path(data = circ, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_segment(data = load_df,
                 aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
                 arrow = arrow(length = unit(arrow_len, "cm"))) +
    ggrepel::geom_text_repel(
      data = load_df, aes(x = TX, y = TY, label = Variable),
      size = label_size, max.overlaps = Inf, box.padding = 0.1,
      point.padding = 0, segment.size = 0, xlim = c(-1, 1), ylim = c(-1, 1), seed = 123
    ) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE, clip = "on") +
    labs(title = panel_label,
         x = sprintf("Dim1 (%.1f%%)", pct[1]),
         y = sprintf("Dim2 (%.1f%%)", pct[2])) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 2, hjust = 0),
      axis.title = element_text(size = axis_title),
      axis.text  = element_text(size = axis_text)
    )
}

# Build 3 panels
p3_a <- biplot_from_df(
  dplyr::filter(dat, Survey == 1),
  weight_col = "WeightSpain", panel_label = "(a)"
)
p3_b <- biplot_from_df(
  dplyr::filter(dat, Survey == 1, Province %in% mining_provinces),
  weight_col = "WeightSpain", panel_label = "(b)"
)
p3_c <- biplot_from_df(
  dplyr::filter(dat, Survey == 2),
  weight_col = if ("WeightAsPontes" %in% names(dat)) "WeightAsPontes" else NULL,
  panel_label = "(c)"
)

final3 <- p3_a | p3_b | p3_c
print(final3)
ggsave("Figure_three_biplots_small_labels_inside_IMPUTED.png",
       final3, width = 10.0 * 0.7, height = 3.8 * 0.7, dpi = 300)


# =============================================================================
# ROBUSTNESS: 3×2 PCA Biplots FULL-BLEED con separación de columnas + márgenes laterales
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(dplyr)
  library(psych)
  library(grid)   # unit()
})

# -------------------- Defaults seguros -------------------------
if (!exists("base_size_big"))  base_size_big  <- 8
if (!exists("axis_title_big")) axis_title_big <- 8
if (!exists("axis_text_big"))  axis_text_big  <- 8

# Tema global sin márgenes (se aplicará en cada panel)
theme_tight <- theme_minimal(base_size = base_size_big) +
  theme(
    plot.margin      = unit(c(0,0,0,0), "pt"),
    panel.spacing    = unit(0, "pt"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.title       = element_text(size = axis_title_big * 0.8,
                                    margin = margin(0,0,0,0)),
    axis.text        = element_text(size = axis_text_big,
                                    margin = margin(0,0,0,0)),
    legend.position  = "none"
  )

# -------------------- Constructor de biplots -------------------
make_biplot_tight <- function(df,
                              weight_col    = "WeightSpain",
                              panel_tag     = "(a)",
                              label_size_mm = 3.2,
                              arrow_len_cm  = 0.18,
                              radial_inset  = 0.70){
  
  need <- unique(c(vars, weight_col))
  d <- df %>%
    dplyr::select(any_of(need)) %>%
    dplyr::mutate(dplyr::across(all_of(vars), as.numeric))
  
  # Pesos
  w <- NULL
  if (!is.null(weight_col) && weight_col %in% names(d)) {
    d <- d %>% dplyr::filter(!is.na(.data[[weight_col]]))
    w <- as.numeric(d[[weight_col]])
  }
  
  # Imputación o filtrado listwise
  if (isTRUE(impute_missing)) {
    d <- impute_weighted_means(d, vars, weight_col)
  } else {
    d <- d %>% dplyr::filter(dplyr::if_all(all_of(vars), ~ !is.na(.x)))
  }
  
  # Correlaciones y PCA
  X <- as.matrix(d[, vars, drop = FALSE])
  R <- stats::cov.wt(X, wt = w, method = "unbiased", cor = TRUE)$cor
  pc  <- psych::principal(r = R, nfactors = 2, rotate = "varimax")
  L   <- as.matrix(unclass(pc$loadings))
  pct <- 100 * colSums(L^2) / nrow(L)
  
  # Círculo y cargas
  t <- seq(0, 2*pi, length.out = 400)
  circ <- data.frame(x = cos(t), y = sin(t))
  load_df <- data.frame(Variable = rownames(L), Dim1 = L[,1], Dim2 = L[,2])
  
  r <- sqrt(load_df$Dim1^2 + load_df$Dim2^2)
  s <- ifelse(r > radial_inset, radial_inset / pmax(r, 1e-9), 1)
  load_df$TX <- load_df$Dim1 * s
  load_df$TY <- load_df$Dim2 * s
  
  tag_x <- -0.98
  tag_y <-  0.98
  
  ggplot() +
    geom_path(data = circ, aes(x, y)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_segment(
      data = load_df,
      aes(x = 0, y = 0, xend = Dim1, yend = Dim2),
      arrow = arrow(length = unit(arrow_len_cm, "cm"))
    ) +
    geom_text(
      data = load_df,
      aes(x = TX, y = TY, label = Variable),
      size = label_size_mm, lineheight = 0.92
    ) +
    annotate("text", x = tag_x, y = tag_y, label = panel_tag,
             hjust = 0, vjust = 1, fontface = "bold", size = 4.2) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1),
                expand = FALSE, clip = "on") +
    labs(
      x = sprintf("Dim1 (%.1f%%)", pct[1]),
      y = sprintf("Dim2 (%.1f%%)", pct[2])
    ) +
    theme_tight
}

# --------------------------- Paneles 3×2 ------------------------
p6_a <- dat_norm %>% filter(Province_norm == "coruna") %>%
  make_biplot_tight(panel_tag = "(a)")
p6_b <- dat_norm %>% filter(Cluster_norm %in% c("a coruna","coruna"),
                            Survey == 1, Dist_to_nearest_plant < dist_thr) %>%
  make_biplot_tight(panel_tag = "(b)")
p6_c <- dat_norm %>% filter(Province_norm %in% c("asturias","leon","palencia")) %>%
  make_biplot_tight(panel_tag = "(c)")
p6_d <- dat_norm %>% filter(Cluster_norm == "asturias",
                            Survey == 1, Dist_to_nearest_plant < dist_thr) %>%
  make_biplot_tight(panel_tag = "(d)")
p6_e <- dat_norm %>% filter(Province_norm == "teruel") %>%
  make_biplot_tight(panel_tag = "(e)")
p6_f <- dat_norm %>% filter(Cluster_norm == "teruel",
                            Survey == 1, Dist_to_nearest_plant < dist_thr) %>%
  make_biplot_tight(panel_tag = "(f)")

# ------------------ Añadir gap entre columnas -------------------
gap_rel <- 0.03  # 3% del ancho como espacio central
row1 <- (p6_a | patchwork::plot_spacer() | p6_b) + plot_layout(widths = c(1, gap_rel, 1))
row2 <- (p6_c | patchwork::plot_spacer() | p6_d) + plot_layout(widths = c(1, gap_rel, 1))
row3 <- (p6_e | patchwork::plot_spacer() | p6_f) + plot_layout(widths = c(1, gap_rel, 1))

final6_patch <- row1 / row2 / row3

# Márgenes: solo 1 cm en izquierda y derecha
final6_patch <- final6_patch + plot_annotation(
  theme = theme(
    plot.margin   = margin(t = 0, r = 10, b = 0, l = 10, unit = "mm"),
    panel.spacing = unit(0, "pt")
  )
)

# ------------------- Full-bleed en el lienzo --------------------
final6_fullbleed <- ggdraw() + draw_plot(final6_patch, x = 0, y = 0, width = 1, height = 1)

print(final6_fullbleed)

# ----------------------------- Export ---------------------------
if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

ggsave("Six_Biplots_3x2_fullbleed_gap_sidemargin.png",
       plot   = final6_fullbleed,
       device = ragg::agg_png,
       width  = 16.5,   # ancho útil Word A4 con márgenes
       height = 24,     # más alto
       units  = "cm",
       dpi    = 600,
       bg     = "white",
       limitsize = FALSE)

ggsave("Six_Biplots_3x2_fullbleed_gap_sidemargin.pdf",
       plot   = final6_fullbleed,
       device = cairo_pdf,
       width  = 16.5,
       height = 24,
       units  = "cm",
       dpi    = 600,
       bg     = "white",
       limitsize = FALSE)


# =============================================================================
# 5) Descriptives (means±SD, medians) for 8 samples — with imputation
# =============================================================================
library(ggplot2)

w_mean <- function(x, w = NULL){ 
  if (is.null(w)) return(mean(x, na.rm = TRUE))
  ok <- is.finite(x) & is.finite(w); if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]; sum(w*x)/sum(w)
}
w_sd <- function(x, w = NULL){
  if (is.null(w)) return(sd(x, na.rm = TRUE))
  ok <- is.finite(x) & is.finite(w); x <- x[ok]; w <- w[ok]
  if (length(x) < 2 || sum(w) == 0) return(NA_real_)
  m <- sum(w*x)/sum(w); sqrt(sum(w*(x-m)^2)/sum(w))
}
w_median <- function(x, w = NULL){
  if (is.null(w)) return(median(x, na.rm = TRUE))
  ok <- is.finite(x) & is.finite(w); x <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)
  o <- order(x); x <- x[o]; w <- w[o]; cw <- cumsum(w); x[which(cw >= sum(w)/2)[1]]
}

summarise_group <- function(df, group_label, weight_col = NULL){
  need <- unique(c(items, weight_col))
  d <- df %>% dplyr::select(any_of(need)) %>% dplyr::mutate(dplyr::across(all_of(items), as.numeric))
  if (!is.null(weight_col) && weight_col %in% names(d)) d <- d %>% dplyr::filter(!is.na(.data[[weight_col]]))
  if (isTRUE(impute_missing)) {
    d <- impute_weighted_means(d, items, weight_col)
  } else {
    d <- d %>% dplyr::filter(dplyr::if_all(all_of(items), ~ !is.na(.x)))
  }
  w <- if (!is.null(weight_col) && weight_col %in% names(d)) as.numeric(d[[weight_col]]) else NULL
  
  tibble::tibble(
    Sample = group_label,
    Item   = factor(items, levels = items),
    n      = nrow(d),
    mean   = sapply(items, \(v) w_mean(d[[v]], w)),
    median = sapply(items, \(v) w_median(d[[v]], w)),
    sd     = sapply(items, \(v) w_sd(d[[v]], w))
  )
}

label_ast_area <- "Asturias/Le\u00F3n/Palencia"
label_ast_cl   <- paste0("Asturias/Le\u00F3n/Palencia < ", dist_thr)
label_cor_area <- "A Coruña (province)"
label_cor_cl   <- paste0("A Coruña < ", dist_thr)
label_ter_area <- "Teruel (province)"
label_ter_cl   <- paste0("Teruel < ", dist_thr)

desc_all <- dplyr::bind_rows(
  summarise_group(dplyr::filter(dat_norm, Survey == 1), "Spain", "WeightSpain"),
  summarise_group(dplyr::filter(dat_norm, Survey == 1, Province_norm %in% c("asturias","leon","palencia")),
                  label_ast_area, "WeightSpain"),
  summarise_group(dplyr::filter(dat_norm, Survey == 2), "As Pontes",
                  if ("WeightAsPontes" %in% names(dat_norm)) "WeightAsPontes" else NULL),
  summarise_group(dplyr::filter(dat_norm, Survey == 1, Cluster_norm %in% c("asturias","leon","palencia"),
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_ast_cl, "WeightSpain"),
  summarise_group(dplyr::filter(dat_norm, Survey == 1, Province_norm == "coruna"),
                  label_cor_area, "WeightSpain"),
  summarise_group(dplyr::filter(dat_norm, Survey == 1, Cluster_norm %in% c("a coruna","coruna"),
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_cor_cl, "WeightSpain"),
  summarise_group(dplyr::filter(dat_norm, Province_norm == "teruel"),
                  label_ter_area, "WeightSpain"),
  summarise_group(dplyr::filter(dat_norm, Survey == 1, Cluster_norm == "teruel",
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_ter_cl, "WeightSpain")
)

# Order panels & build wrapped labels with N
layout_order <- c("Spain", label_ast_area, label_cor_area, label_ter_area,
                  "As Pontes", label_ast_cl, label_cor_cl, label_ter_cl)
desc_all <- desc_all %>% dplyr::mutate(Sample = factor(Sample, levels = layout_order))

letters_vec  <- paste0("(", letters[1:8], ") ")
panel_titles <- paste0(letters_vec, layout_order)
n_by_sample  <- desc_all %>% dplyr::distinct(Sample, n) %>% dplyr::rename(n_fac = n)
label_levels <- paste0(panel_titles, "\n(n=", n_by_sample$n_fac[match(layout_order, n_by_sample$Sample)], ")")

desc_all <- desc_all %>%
  dplyr::left_join(n_by_sample, by = "Sample") %>%
  dplyr::mutate(SampleLab = factor(
    paste0(letters_vec[match(as.character(Sample), layout_order)],
           as.character(Sample), "\n(n=", .data[["n_fac"]], ")"),
    levels = label_levels
  ))

# Facet labeller that auto-wraps long strip titles
wrap_labeller <- labeller(SampleLab = label_wrap_gen(width = 18))

# Plot (tweaked spacing & strip styling to avoid overlaps)
base_size   <- 17
axis_text   <- 15
strip_size  <- 14       # slightly smaller than before
legend_text <- 15

p_desc <- ggplot(desc_all, aes(x = Item, y = mean, color = Item)) +
  geom_hline(yintercept = seq(0, 10, by = 2), color = "grey90", linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, mean - sd), ymax = pmin(10, mean + sd)),
                width = 0.18, linewidth = 1.0, alpha = 0.95) +
  geom_segment(aes(xend = Item, y = 0, yend = mean),
               linewidth = 0.45, alpha = 0.65, color = "grey55") +
  geom_point(size = 3.3) +
  geom_point(aes(y = median), shape = 124, size = 6, color = "black") +
  facet_wrap(~ SampleLab, nrow = 2, ncol = 4, labeller = wrap_labeller) +
  scale_x_discrete(expand = expansion(mult = c(0.10, 0.10))) +
  scale_y_continuous(limits = c(0, 10), breaks = 0:10,
                     expand = expansion(mult = c(0.03, 0.07))) +
  scale_color_manual(values = c(
    Reindustrialization = "#1f77b4",
    Renaturalization    = "#2ca02c",
    Preservation        = "#ff7f0e"
  )) +
  labs(x = NULL, y = "Score (mean ± SD)", color = NULL) +
  theme_minimal(base_size = base_size) +
  theme(
    legend.position    = "bottom",
    legend.text        = element_text(size = legend_text),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(16, "pt"),  # more breathing room
    panel.spacing.y    = unit(20, "pt"),
    strip.placement    = "outside",
    strip.text.x       = element_text(size = strip_size, lineheight = 1.14,
                                      margin = margin(t = 6, b = 12)),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.text.y        = element_text(size = axis_text, margin = margin(r = 4)),
    plot.margin        = margin(20, 24, 20, 30)
  ) +
  coord_cartesian(clip = "off")

print(p_desc)
ggsave("descriptives_all_samples_2x4_shorttitles_noXlabels_IMPUTED.png",
       p_desc, width = 9.0, height = 8.2, dpi = 400)
ggsave("descriptives_all_samples_2x4_shorttitles_noXlabels_IMPUTED.pdf",
       p_desc, width = 9.0, height = 8.2)


# =============================================================================
# 6) Eigenvalues across samples (weighted PCA on correlations) — with imputation
# =============================================================================
weighted_cor <- function(df_items, w = NULL) {
  stats::cov.wt(as.matrix(df_items), wt = w, method = "unbiased", cor = TRUE)$cor
}
prep_group <- function(df, weight_col = NULL) {
  need <- unique(c(items, weight_col))
  d <- df %>% dplyr::select(any_of(need)) %>% dplyr::mutate(dplyr::across(all_of(items), as.numeric))
  if (!is.null(weight_col) && weight_col %in% names(d)) d <- d %>% dplyr::filter(!is.na(.data[[weight_col]]))
  if (isTRUE(impute_missing)) d <- impute_weighted_means(d, items, weight_col) else d <- d %>% dplyr::filter(dplyr::if_all(all_of(items), ~ !is.na(.x)))
  w <- if (!is.null(weight_col) && weight_col %in% names(d)) as.numeric(d[[weight_col]]) else NULL
  list(data = d[, items, drop = FALSE], w = w, n = nrow(d))
}
eigs_for_sample <- function(df, label, weight_col = NULL) {
  pg <- prep_group(df, weight_col)
  if (pg$n == 0) {
    return(tibble(
      Sample = label, n_used = 0, Component = integer(),
      Eigenvalue = numeric(), PropVar = numeric(), CumProp = numeric(),
      FactorsExtracted = integer()
    ))
  }
  R  <- weighted_cor(pg$data, pg$w)
  ev <- sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  prop <- ev / length(ev)
  kaiser_k <- sum(ev > 1)
  tibble(
    Sample = label, n_used = pg$n, Component = seq_along(ev),
    Eigenvalue = ev, PropVar = prop, CumProp = cumsum(prop),
    FactorsExtracted = kaiser_k
  )
}

tbl_list <- list(
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1), "Spain", "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Province_norm %in% c("asturias","leon","palencia")),
                  label_ast_area, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 2), "As Pontes",
                  if ("WeightAsPontes" %in% names(dat_norm)) "WeightAsPontes" else NULL),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1,
                                Cluster_norm %in% c("asturias","leon","palencia"),
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_ast_cl, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Province_norm == "coruna"),
                  label_cor_area, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1,
                                Cluster_norm %in% c("a coruna","coruna"),
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_cor_cl, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Province_norm == "teruel"),
                  label_ter_area, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Cluster_norm == "teruel",
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_ter_cl, "WeightSpain")
)
eigs_all <- dplyr::bind_rows(tbl_list)

layout_order <- c("Spain", label_ast_area, label_cor_area, label_ter_area,
                  "As Pontes", label_ast_cl, label_cor_cl, label_ter_cl)
eigs_all <- eigs_all %>% dplyr::mutate(Sample = factor(Sample, levels = layout_order)) %>%
  dplyr::arrange(Sample, Component)

eigs_fmt <- eigs_all %>%
  dplyr::mutate(
    Component  = paste0("PC", Component),
    Eigenvalue = sprintf("%.3f", Eigenvalue),
    `% Var`    = sprintf("%.1f%%", PropVar * 100),
    `Cum %`    = sprintf("%.1f%%", CumProp * 100)
  ) %>%
  dplyr::select(Sample, n_used, Component, Eigenvalue, `% Var`, `Cum %`, FactorsExtracted)

readr::write_csv(eigs_fmt, "eigenvalues_across_samples_IMPUTED.csv")

tt <- gridExtra::ttheme_minimal(
  core   = list(fg_params = list(fontface = "plain", cex = 0.75), padding = unit(c(6,6), "pt")),
  colhead = list(fg_params = list(fontface = 2, cex = 0.85), padding = unit(c(6,6), "pt"))
)
tbl_grob <- gridExtra::tableGrob(eigs_fmt, rows = NULL, theme = tt)
title_grob <- grid::textGrob("PCA eigenvalues (weighted correlations) by sample — IMPUTED",
                             gp = grid::gpar(fontface = "bold", cex = 0.95))
padding <- unit(6, "pt")
tbl_with_title <- gtable::gtable_add_rows(tbl_grob, heights = grid::grobHeight(title_grob) + padding, pos = 0)
tbl_with_title <- gtable::gtable_add_grob(tbl_with_title, title_grob, t = 1, l = 1, r = ncol(tbl_grob))

png("eigenvalues_across_samples_table_IMPUTED.png", width = 1600, height = 1200, res = 200)
grid::grid.newpage(); grid::grid.draw(tbl_with_title)
dev.off()
print(eigs_fmt)
ggplot2::ggsave("eigenvalues_across_samples_table_IMPUTED.png", plot = ggplotify::as.ggplot(tbl_with_title), width = 8, height = 6, dpi = 200, bg = "white")

# =============================================================================
# 7) Scree plots (matching descriptives) — with imputation
# =============================================================================
eigs_all_scree <- dplyr::bind_rows(
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1), "Spain", "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Province_norm %in% c("asturias","leon","palencia")),
                  label_ast_area, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 2), "As Pontes",
                  if ("WeightAsPontes" %in% names(dat_norm)) "WeightAsPontes" else NULL),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Cluster_norm %in% c("asturias","leon","palencia"),
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_ast_cl, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Province_norm == "coruna"),
                  label_cor_area, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Cluster_norm %in% c("a coruna","coruna"),
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_cor_cl, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Province_norm == "teruel"), label_ter_area, "WeightSpain"),
  eigs_for_sample(dplyr::filter(dat_norm, Survey == 1, Cluster_norm == "teruel",
                                !is.na(Dist_to_nearest_plant), Dist_to_nearest_plant < dist_thr),
                  label_ter_cl, "WeightSpain")
)

letters_vec  <- paste0("(", letters[1:8], ") ")
n_by_sample2 <- eigs_all_scree %>% dplyr::distinct(Sample, n_used)
panel_titles <- paste0(letters_vec, layout_order)
label_levels <- paste0(panel_titles, "\n(n=", n_by_sample2$n_used[match(layout_order, n_by_sample2$Sample)], ")")

plot_df <- eigs_all_scree %>%
  dplyr::mutate(
    Sample = factor(Sample, levels = layout_order),
    SampleLab = factor(
      paste0(letters_vec[match(as.character(Sample), layout_order)],
             as.character(Sample), "\n(n=",
             n_by_sample2$n_used[match(as.character(Sample), n_by_sample2$Sample)], ")"),
      levels = label_levels
    ),
    Component = as.numeric(Component)
  ) %>% dplyr::arrange(Sample)

lvls <- levels(plot_df$SampleLab)
lvls <- sub("Asturias/Le\u00F3n/Palencia( < 150)?",
            "Asturias/Le\u00F3n/\nPalencia\\1", lvls, perl = TRUE)
levels(plot_df$SampleLab) <- lvls

p_scree <- ggplot(plot_df, aes(x = Component, y = Eigenvalue, group = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.7, alpha = 0.7) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.8) +
  facet_wrap(~ SampleLab, nrow = 2, ncol = 4) +
  scale_x_continuous(breaks = 1:length(items)) +
  labs(x = "Component", y = "Eigenvalue") +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.spacing.x  = unit(12, "pt"),
    panel.spacing.y  = unit(14, "pt"),
    strip.placement  = "outside",
    strip.text.x     = element_text(size = 14.5, lineheight = 1.14, margin = margin(t = 6, b = 12)),
    plot.margin      = margin(18, 24, 18, 32)
  ) +
  coord_cartesian(clip = "off")

print(p_scree)
ggsave("scree_plots_matching_descriptives_IMPUTED.png", p_scree, width = 12, height = 7.2, dpi = 350)
ggsave("scree_plots_matching_descriptives_IMPUTED.pdf",  p_scree, width = 12, height = 7.2)

#====================================================================
# 8) Communalities
# Communalities (h^2) for Spain, As Pontes, Asturias/León/Palencia
# PCA: 2 components, varimax, on weighted correlation matrices
# Saves: communalities_plot.png / .pdf and communalities.csv
# ============================================================
library(tidyverse)
library(haven)
library(stringi)
library(psych)

# ---- Settings ----
file_sav <- "Dataset_Nature_heritage.sav"
items    <- c("Reindustrialization","Renaturalization","Preservation")

# ---- Load + normalize province/cluster names ----
dat <- read_sav(file_sav) %>%
  mutate(
    Province_norm = tolower(stri_trans_general(Province, "Latin-ASCII")),
    Cluster_norm  = tolower(stri_trans_general(PowerPlantCluster, "Latin-ASCII"))
  )

# ---- Helper: communalities from weighted correlation PCA (2 comps, varimax) ----
get_communalities <- function(df, label, weight_col = NULL) {
  need <- unique(c(items, weight_col))
  d <- df %>%
    select(any_of(need)) %>%
    mutate(across(all_of(items), as.numeric)) %>%
    filter(if_all(all_of(items), ~ !is.na(.x)))
  if (!is.null(weight_col) && weight_col %in% names(d)) {
    d <- d %>% filter(!is.na(.data[[weight_col]]))
    w <- as.numeric(d[[weight_col]])
  } else {
    w <- NULL
  }
  X <- as.matrix(d[, items, drop = FALSE])
  
  # Weighted correlation matrix
  R <- stats::cov.wt(X, wt = w, method = "unbiased", cor = TRUE)$cor
  
  # PCA: 2 components, varimax
  pc <- psych::principal(r = R, nfactors = 2, rotate = "varimax")
  L  <- as.matrix(unclass(pc$loadings))
  
  # Communalities: sum of squared loadings
  h2 <- rowSums(L^2)
  
  tibble(
    Sample = label,
    Item   = factor(rownames(L), levels = items),
    h2     = as.numeric(h2)
  )
}

# ---- Build the three samples (labels WITHOUT survey tags) ----
comm_spain <- get_communalities(
  df = filter(dat, Survey == 1),
  label = "Spain",
  weight_col = "WeightSpain"
)

comm_aspontes <- get_communalities(
  df = filter(dat, Survey == 2),
  label = "As Pontes",
  weight_col = if ("WeightAsPontes" %in% names(dat)) "WeightAsPontes" else NULL
)

comm_alp <- get_communalities(
  df = filter(dat, Survey == 1, Province_norm %in% c("asturias","leon","palencia")),
  label = "Asturias/León/Palencia",
  weight_col = "WeightSpain"
)

comm_all <- bind_rows(comm_spain, comm_aspontes, comm_alp) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = c("Spain", "Asturias/León/Palencia", "As Pontes")
    )
  )

# ---- Save tidy table ----
readr::write_csv(comm_all %>% arrange(Sample, Item), "communalities.csv")

# ---- Plot: communalities (0–1), faceted by sample; NO x-axis labels or ticks ----
p <- ggplot(comm_all, aes(x = Item, y = h2, fill = Item)) +
  geom_col(width = 0.68, alpha = 0.95) +
  geom_text(aes(label = sprintf("%.2f", h2)), vjust = -0.35, size = 4) +
  facet_wrap(~ Sample, nrow = 1) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2),
                     expand = expansion(mult = c(0.02, 0.10))) +
  scale_fill_manual(values = c(
    Reindustrialization = "#1f77b4",
    Renaturalization    = "#2ca02c",
    Preservation        = "#ff7f0e"
  )) +
  labs(
    title = "",
    subtitle = "",
    x = NULL, y = "h² (0–1)", fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = 13, face = "bold"),
    axis.text.x        = element_blank(),      # remove x-axis labels
    axis.ticks.x       = element_blank(),      # remove x-axis ticks
    plot.title         = element_text(face = "bold")
  )

print(p)

# ---- Save images ----
ggsave("communalities_plot.png", p, width = 12, height = 4.5, dpi = 350)
ggsave("communalities_plot.pdf", p, width = 12, height = 4.5)
