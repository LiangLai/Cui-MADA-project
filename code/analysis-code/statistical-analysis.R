###############################################################################
# statistical-analysis.R
#
# Formal statistical analysis for the 5-body-site HMP microbiome project.
# This script loads the processed phyloseq data and performs:
#   - H1: Alpha diversity comparison across body sites (Kruskal-Wallis,
#         pairwise Wilcoxon with BH correction)
#   - H2: Beta diversity / clustering (Bray-Curtis + Jaccard; PERMANOVA with
#         BETADISPER; PCoA; hierarchical clustering)
#   - H3: Microbial co-occurrence networks per body site (Spearman-based
#         correlation networks with Louvain community detection)
#
# All results (tables, figures, model objects) are written to
# results/tables and results/figures.
#
# Usage:  Run after processingcode.R / processingfile-v2.qmd have produced
#         data/processed-data/ps_filt.rds
###############################################################################

## ---- packages --------
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(phyloseq)
  library(vegan)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(Hmisc)
  library(broom)
  library(RColorBrewer)
  library(tidymodels)
})
tidymodels_prefer()

options(lifecycle_verbosity = "quiet")
theme_set(theme_bw(base_size = 12))

colors_body <- c(
  "Airways"    = "#87CEEB",
  "Gut"        = "#8B4513",
  "Oral"       = "#DC143C",
  "Skin"       = "#FFB6C1",
  "Urogenital" = "#9370DB"
)

fig_dir   <- here("results", "figures")
tab_dir   <- here("results", "tables")
out_dir   <- here("results", "output")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


## ---- loaddata --------
# IMPORTANT: ps_filt.rds from the processing step may contain the full
# filtered HMP dataset (~4500 samples x ~34000 taxa). Running PERMANOVA /
# network analysis on that is impractical (the Bray distance alone has
# ~10M pairs). We therefore apply the stratified subsampling + prevalence
# filter + rarefaction here, mirroring the group's original R script. The
# size parameters are tuned so the full pipeline completes in a few minutes.

# --- Tunable knobs (keep small so the pipeline finishes reliably) -----------
N_PER_SITE       <- 200    # max samples per body site
MIN_TAXA_PREV    <- 10     # keep taxa present in >= this many samples
MIN_TAXA_ABUND   <- 50     # keep taxa with total count >= this
RAREFY_PCTILE    <- 0.05   # rarefy to this percentile of sample depths
N_PERM           <- 99     # permutations for PERMANOVA / BETADISPER
# ---------------------------------------------------------------------------

ps <- readRDS(here("data", "processed-data", "ps_filt.rds"))

# Ensure the simplified body-site column exists
if (!"SampleType" %in% colnames(sample_data(ps))) {
  sample_data(ps)$SampleType <- sample_data(ps)$HMP_BODY_SITE |>
    gsub("Gastrointestinal Tract", "Gut", x = _) |>
    gsub("Urogenital Tract", "Urogenital", x = _)
}

# Keep only the five major body sites
keep_sites <- names(colors_body)
ps <- prune_samples(
  sample_data(ps)$SampleType %in% keep_sites, ps
)
cat(sprintf("[setup] %d samples across %d sites before subsampling\n",
            nsamples(ps), length(unique(sample_data(ps)$SampleType))))

# Stratified subsample: up to N_PER_SITE per site (or all, if fewer)
# NOTE: avoid phyloseq::subset_samples() here — it uses NSE and won't see
# `site` from a closure reliably.
set.seed(123)
site_vec <- as.character(sample_data(ps)$SampleType)
names(site_vec) <- sample_names(ps)
selected <- unlist(lapply(keep_sites, function(s) {
  site_names <- names(site_vec)[site_vec == s]
  if (length(site_names) == 0) return(character(0))
  sample(site_names, min(N_PER_SITE, length(site_names)))
}))
ps <- prune_samples(selected, ps)

# Drop rare taxa (present in too few samples OR too low total abundance)
cat(sprintf("[setup] %d samples, %d taxa before prevalence/abundance filter\n",
            nsamples(ps), ntaxa(ps)))
ps <- filter_taxa(
  ps,
  function(x) sum(x > 0) >= MIN_TAXA_PREV & sum(x) >= MIN_TAXA_ABUND,
  prune = TRUE
)
ps <- prune_samples(sample_sums(ps) > 0, ps)
cat(sprintf("[setup] %d samples, %d taxa after filter\n",
            nsamples(ps), ntaxa(ps)))

# Rarefy to an even depth (RAREFY_PCTILE-th percentile) so that
# alpha/beta diversity are not confounded by sequencing depth.
set.seed(123)
min_depth <- sort(sample_sums(ps))[max(1, floor(nsamples(ps) * RAREFY_PCTILE))]
cat(sprintf("[setup] rarefying to %d reads/sample\n", min_depth))
ps_rare <- rarefy_even_depth(ps, sample.size = min_depth, rngseed = 123,
                             verbose = FALSE)
cat(sprintf("[setup] final dataset: %d samples x %d taxa\n",
            nsamples(ps_rare), ntaxa(ps_rare)))


###############################################################################
## H1: Alpha diversity differs by body site
###############################################################################

## ---- alpha_diversity --------
alpha_div <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Shannon", "Simpson", "Chao1")
)
alpha_div$SampleType <- sample_data(ps_rare)$SampleType

metrics <- c("Observed", "Shannon", "Simpson", "Chao1")

# Overall Kruskal-Wallis tests
kw_table <- purrr::map_dfr(metrics, function(m) {
  kw <- kruskal.test(as.formula(paste(m, "~ SampleType")), data = alpha_div)
  tibble(metric = m,
         chi_sq = unname(kw$statistic),
         df     = unname(kw$parameter),
         p      = kw$p.value)
})
saveRDS(kw_table, file.path(tab_dir, "alpha_kruskal.rds"))

# Pairwise Wilcoxon with Benjamini-Hochberg correction (long format)
pairwise_list <- lapply(metrics, function(m) {
  pw <- pairwise.wilcox.test(alpha_div[[m]], alpha_div$SampleType,
                             p.adjust.method = "BH")
  pw_df <- as.data.frame(as.table(pw$p.value))
  names(pw_df) <- c("site_1", "site_2", "p_adj")
  pw_df$metric <- m
  pw_df |> filter(!is.na(p_adj))
})
pairwise_tbl <- bind_rows(pairwise_list) |>
  select(metric, site_1, site_2, p_adj)
saveRDS(pairwise_tbl, file.path(tab_dir, "alpha_pairwise.rds"))

# Boxplot across all four indices
alpha_long <- alpha_div |>
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(alpha_long,
                  aes(x = SampleType, y = Value, fill = SampleType)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = colors_body) +
  labs(title = "Alpha diversity across five major body sites",
       subtitle = sprintf("Rarefied to %d reads/sample", min_depth),
       x = "Body Site", y = "Diversity value") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "alpha_diversity_5sites.png"),
       p_alpha, width = 10, height = 8, dpi = 300)

###############################################################################
## H2: Samples cluster by body site (beta diversity)
###############################################################################
# Note: the previous version used phyloseq::distance() / ordinate() repeatedly,
# which recomputed the Bray-Curtis distance matrix up to three times and
# sometimes hung on large objects. This version extracts the OTU matrix
# ONCE, uses vegan::vegdist directly, and reuses the same distance matrix for
# PERMANOVA, BETADISPER, PCoA (via cmdscale) and hclust.

## ---- beta_prep --------
cat("[H2] Preparing OTU matrices...\n")

otu_counts <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_counts <- t(otu_counts)   # rows = samples
storage.mode(otu_counts) <- "double"

otu_log <- log1p(otu_counts)                              # log(1+x)
otu_pa  <- (otu_counts > 0) + 0L                          # presence/absence

metadata <- data.frame(sample_data(ps_rare))
metadata <- metadata[rownames(otu_counts), , drop = FALSE]
metadata$SampleType <- factor(metadata$SampleType)

cat(sprintf("[H2] %d samples x %d taxa\n",
            nrow(otu_counts), ncol(otu_counts)))

## ---- beta_permanova_bray --------
cat("[H2] Computing Bray-Curtis distance...\n")
dmat_bray <- vegan::vegdist(otu_log, method = "bray")

cat("[H2] PERMANOVA (Bray)...\n")
set.seed(123)
perm_bray <- vegan::adonis2(dmat_bray ~ SampleType,
                            data = metadata, permutations = N_PERM)

cat("[H2] BETADISPER (Bray)...\n")
set.seed(123)
bdisp_bray <- vegan::betadisper(dmat_bray, metadata$SampleType)
bperm_bray <- vegan::permutest(bdisp_bray, permutations = N_PERM)

## ---- beta_permanova_jaccard --------
cat("[H2] Computing Jaccard distance...\n")
dmat_jaccard <- vegan::vegdist(otu_pa, method = "jaccard", binary = TRUE)

cat("[H2] PERMANOVA (Jaccard)...\n")
set.seed(123)
perm_jaccard <- vegan::adonis2(dmat_jaccard ~ SampleType,
                               data = metadata, permutations = N_PERM)

cat("[H2] BETADISPER (Jaccard)...\n")
set.seed(123)
bdisp_jaccard <- vegan::betadisper(dmat_jaccard, metadata$SampleType)
bperm_jaccard <- vegan::permutest(bdisp_jaccard, permutations = N_PERM)

permanova_tbl <- bind_rows(
  tibble(distance     = "bray",
         R2           = perm_bray$R2[1],
         F            = perm_bray$F[1],
         p_permanova  = perm_bray$`Pr(>F)`[1],
         F_betadisper = bperm_bray$tab$F[1],
         p_betadisper = bperm_bray$tab$`Pr(>F)`[1]),
  tibble(distance     = "jaccard",
         R2           = perm_jaccard$R2[1],
         F            = perm_jaccard$F[1],
         p_permanova  = perm_jaccard$`Pr(>F)`[1],
         F_betadisper = bperm_jaccard$tab$F[1],
         p_betadisper = bperm_jaccard$tab$`Pr(>F)`[1])
)
saveRDS(permanova_tbl, file.path(tab_dir, "permanova.rds"))

## ---- pcoa --------
# Use cmdscale on the already-computed Bray distance; do NOT call ordinate()
# (which would recompute the distance matrix internally).
cat("[H2] PCoA (cmdscale on Bray)...\n")
pcoa_bray <- cmdscale(dmat_bray, k = 2, eig = TRUE)

eig_vals <- pcoa_bray$eig
var1 <- round(100 * eig_vals[1] / sum(abs(eig_vals)), 1)
var2 <- round(100 * eig_vals[2] / sum(abs(eig_vals)), 1)

ord_df <- data.frame(
  Axis.1     = pcoa_bray$points[, 1],
  Axis.2     = pcoa_bray$points[, 2],
  SampleType = metadata$SampleType
)

p_pcoa <- ggplot(ord_df,
                 aes(x = Axis.1, y = Axis.2, color = SampleType)) +
  geom_point(size = 2.5, alpha = 0.75) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_color_manual(values = colors_body) +
  labs(title = "PCoA (Bray-Curtis) - Five body sites",
       subtitle = sprintf(
         "PERMANOVA R2 = %.3f, p = %.3g",
         permanova_tbl$R2[permanova_tbl$distance == "bray"],
         permanova_tbl$p_permanova[permanova_tbl$distance == "bray"]
       ),
       x = sprintf("PC1 (%.1f%% variance)", var1),
       y = sprintf("PC2 (%.1f%% variance)", var2),
       color = "Body Site")

ggsave(file.path(fig_dir, "pcoa_ordination_5sites.png"),
       plot = p_pcoa, width = 9, height = 7, dpi = 300)

## ---- hclust --------
cat("[H2] Hierarchical clustering (Ward D2 on Bray)...\n")
hc <- hclust(dmat_bray, method = "ward.D2")

sample_types  <- as.character(metadata$SampleType[hc$order])
sample_colors <- colors_body[sample_types]

png(file.path(fig_dir, "hierarchical_clustering_5sites.png"),
    width = 14, height = 10, units = "in", res = 300)
layout(matrix(c(1, 2, 3), nrow = 3), heights = c(8, 0.5, 0.5))

par(mar = c(0, 4, 4, 2))
plot(hc, labels = FALSE, xlab = "",
     main = sprintf("Hierarchical clustering - Five body sites (n = %d)",
                    nrow(otu_counts)),
     ylab = "Bray-Curtis distance (Ward's method)")

par(mar = c(0, 4, 0, 2))
barplot(rep(1, length(sample_colors)), col = sample_colors, border = NA,
        space = 0, axes = FALSE, ylab = "")

par(mar = c(1, 0, 0, 0))
plot.new()
legend("center", legend = names(colors_body), fill = colors_body,
       cex = 1.3, title = "Body Site", bty = "n", horiz = TRUE)
dev.off()

cat("[H2] Done.\n")


###############################################################################
## H3: Microbial co-occurrence networks
###############################################################################

## ---- network_function --------
build_cooccurrence_network <- function(ps_subset, site_name,
                                       cor_threshold = 0.5,
                                       p_threshold   = 0.01,
                                       max_taxa      = 500) {

  otu_mat <- as(otu_table(ps_subset), "matrix")
  if (taxa_are_rows(ps_subset)) otu_mat <- t(otu_mat)

  # Prevalence filter: taxa must be present in at least 10% of samples
  prevalence <- colSums(otu_mat > 0) / nrow(otu_mat)
  otu_mat    <- otu_mat[, prevalence >= 0.1, drop = FALSE]

  # Cap the number of taxa by total abundance so rcorr() does not blow up
  # on very dense sites. This keeps runtime bounded regardless of input size.
  if (ncol(otu_mat) > max_taxa) {
    totals <- colSums(otu_mat)
    keep   <- names(sort(totals, decreasing = TRUE))[seq_len(max_taxa)]
    otu_mat <- otu_mat[, keep, drop = FALSE]
  }
  cat(sprintf("[H3] %s: %d samples, %d taxa (after cap)\n",
              site_name, nrow(otu_mat), ncol(otu_mat)))

  if (ncol(otu_mat) < 3) return(NULL)

  cor_result <- Hmisc::rcorr(otu_mat, type = "spearman")
  cor_matrix <- cor_result$r
  p_matrix   <- cor_result$P

  cor_matrix[abs(cor_matrix) < cor_threshold] <- 0
  cor_matrix[p_matrix >= p_threshold]         <- 0
  diag(cor_matrix) <- 0

  g <- igraph::graph_from_adjacency_matrix(
    cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  g <- igraph::delete.vertices(g, which(igraph::degree(g) == 0))

  if (igraph::vcount(g) == 0) return(g)

  # Attach taxonomy
  tax_info     <- as.data.frame(tax_table(ps_subset))
  vertex_names <- igraph::V(g)$name
  for (col in c("Phylum", "Family", "Genus")) {
    if (col %in% colnames(tax_info)) {
      igraph::V(g)$.data <- tax_info[vertex_names, col]
      names(igraph::vertex_attr(g))[length(igraph::vertex_attr_names(g))] <- col
    }
  }

  # Centrality
  igraph::V(g)$degree      <- igraph::degree(g)
  igraph::V(g)$betweenness <- igraph::betweenness(g)
  igraph::V(g)$closeness   <- igraph::closeness(g)
  igraph::V(g)$eigenvector <- igraph::eigen_centrality(g)$vector

  # Edges
  igraph::E(g)$correlation <- igraph::E(g)$weight
  igraph::E(g)$edge_type   <- ifelse(igraph::E(g)$weight > 0,
                                     "Positive", "Negative")

  # Communities (Louvain)
  communities <- igraph::cluster_louvain(g)
  igraph::V(g)$community <- igraph::membership(communities)
  attr(g, "site") <- site_name

  g
}

## ---- build_networks --------
networks <- list()
network_stats <- data.frame()

site_vec_rare <- as.character(sample_data(ps_rare)$SampleType)
names(site_vec_rare) <- sample_names(ps_rare)

for (site in names(colors_body)) {
  site_samples <- names(site_vec_rare)[site_vec_rare == site]
  if (length(site_samples) == 0) next
  ps_site <- prune_samples(site_samples, ps_rare)
  ps_site <- prune_taxa(taxa_sums(ps_site) > 0, ps_site)

  if (nsamples(ps_site) < 20) next

  g <- build_cooccurrence_network(ps_site, site)
  if (is.null(g) || igraph::vcount(g) == 0) next

  networks[[site]] <- g
  network_stats <- rbind(network_stats, data.frame(
    BodySite    = site,
    Nodes       = igraph::vcount(g),
    Edges       = igraph::ecount(g),
    Density     = igraph::edge_density(g),
    AvgDegree   = mean(igraph::degree(g)),
    Clustering  = igraph::transitivity(g, type = "global"),
    Modularity  = igraph::modularity(igraph::cluster_louvain(g)),
    Communities = length(unique(igraph::V(g)$community))
  ))
}

saveRDS(network_stats, file.path(tab_dir, "network_stats.rds"))
saveRDS(networks,      file.path(out_dir, "networks.rds"))

## ---- network_plot --------
stats_long <- network_stats |>
  pivot_longer(cols = c(Nodes, Edges, Density, AvgDegree, Clustering,
                        Modularity, Communities),
               names_to = "Metric", values_to = "Value") |>
  mutate(Metric = factor(Metric,
                         levels = c("Nodes", "Edges", "Density", "AvgDegree",
                                    "Clustering", "Modularity", "Communities")))

p_network_stats <- ggplot(stats_long,
                          aes(x = BodySite, y = Value, fill = BodySite)) +
  geom_bar(stat = "identity", alpha = 0.85) +
  facet_wrap(~Metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = colors_body) +
  labs(title = "Co-occurrence network properties across body sites",
       subtitle = "Spearman |r| > 0.5, p < 0.01",
       x = "Body Site", y = "Value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "network_statistics_comparison.png"),
       p_network_stats, width = 12, height = 8, dpi = 300)


###############################################################################
## H4: Predicting body site with tidymodels (multinomial GLM vs random forest)
###############################################################################
# Uses ps_rare from above. Top-abundance taxa as features, stratified 75/25
# split, 5-fold CV, compare two tuned models, report test-set metrics.

cat("[H4] Preparing feature matrix...\n")
otu_feat <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_feat <- t(otu_feat)
top <- names(sort(colSums(otu_feat), decreasing = TRUE))[seq_len(min(300, ncol(otu_feat)))]
ml_df <- as_tibble(log1p(otu_feat[, top, drop = FALSE]),
                   .name_repair = ~ make.names(.x, unique = TRUE)) |>
  mutate(SampleType = factor(sample_data(ps_rare)$SampleType,
                             levels = names(colors_body)))

set.seed(123)
split <- initial_split(ml_df, prop = 0.75, strata = SampleType)
folds <- vfold_cv(training(split), v = 5, strata = SampleType)

rec <- recipe(SampleType ~ ., data = training(split)) |>
  step_zv(all_predictors()) |> step_normalize(all_numeric_predictors())

glmnet_wf <- workflow() |> add_recipe(rec) |>
  add_model(multinom_reg(penalty = tune(), mixture = tune()) |>
              set_engine("glmnet"))
rf_wf <- workflow() |> add_recipe(rec) |>
  add_model(rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
              set_engine("ranger", importance = "impurity") |>
              set_mode("classification"))

mset <- metric_set(accuracy, roc_auc)

cat("[H4] Tuning elastic-net multinomial GLM...\n")
set.seed(123)
glm_res <- tune_grid(glmnet_wf, folds,
                     grid = grid_regular(penalty(c(-4, 0)), mixture(c(0, 1)),
                                         levels = c(8, 3)),
                     metrics = mset)

cat("[H4] Tuning random forest...\n")
set.seed(123)
rf_res <- tune_grid(rf_wf, folds,
                    grid = grid_regular(mtry(c(5L, 100L)), min_n(c(2L, 20L)),
                                        levels = c(4, 4)),
                    metrics = mset)

best_glm <- select_best(glm_res, metric = "roc_auc")
best_rf  <- select_best(rf_res,  metric = "roc_auc")

cat("[H4] Final fit on test set...\n")
glm_fit <- finalize_workflow(glmnet_wf, best_glm) |> last_fit(split, metrics = mset)
rf_fit  <- finalize_workflow(rf_wf,     best_rf)  |> last_fit(split, metrics = mset)

ml_test_metrics <- bind_rows(
  collect_metrics(glm_fit) |> mutate(model = "multinom_glmnet"),
  collect_metrics(rf_fit)  |> mutate(model = "random_forest")
)
saveRDS(ml_test_metrics, file.path(tab_dir, "ml_test_metrics.rds"))

# Use the better model for the confusion matrix + keep the tuned hyperparams
best_name <- ml_test_metrics |> filter(.metric == "roc_auc") |>
  slice_max(.estimate, n = 1) |> pull(model)
best_fit  <- if (best_name == "random_forest") rf_fit else glm_fit
best_par  <- if (best_name == "random_forest") best_rf else best_glm

saveRDS(bind_rows(
  tibble(model = "multinom_glmnet", param = names(best_glm),
         value = as.character(unlist(best_glm))),
  tibble(model = "random_forest", param = names(best_rf),
         value = as.character(unlist(best_rf)))
) |> filter(param != ".config"),
file.path(tab_dir, "ml_best_params.rds"))

cm <- conf_mat(collect_predictions(best_fit),
               truth = SampleType, estimate = .pred_class)
saveRDS(list(model = best_name, cm = cm),
        file.path(tab_dir, "ml_confusion_matrix.rds"))
ggsave(file.path(fig_dir, "ml_confusion_matrix.png"),
       autoplot(cm, type = "heatmap") +
         labs(title = sprintf("Confusion matrix - %s (test set)", best_name)),
       width = 7, height = 6, dpi = 300)

cat(sprintf("[H4] Best model: %s (test roc_auc = %.3f)\n",
            best_name,
            ml_test_metrics$.estimate[ml_test_metrics$model == best_name &
                                      ml_test_metrics$.metric == "roc_auc"]))

## ---- finish --------
cat("Statistical analysis complete. Outputs written to:\n")
cat("  ", fig_dir, "\n")
cat("  ", tab_dir, "\n")
cat("  ", out_dir, "\n")
