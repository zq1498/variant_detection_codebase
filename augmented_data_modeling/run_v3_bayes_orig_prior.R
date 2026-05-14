suppressPackageStartupMessages({
  library(tidymodels)
  library(rstanarm)
  library(yardstick)
  library(ggplot2)
  library(tidyr)
})

options(mc.cores = parallel::detectCores())
set.seed(2024)

setwd("augmented_data_modeling")

out_dir = "v3_outputs_orig_prior"
dir.create(out_dir, showWarnings = FALSE)

# ============================================================================
# Concept
# ----------------------------------------------------------------------------
# Prior data:    original screen (data/covid_screen_data_cleaned.csv) — the
#                same data used to train the LASSO/RF feature-selection models.
#                Restricted to dAb13's panel (matches the cAb columns that
#                appear in the augmented data).
# Update data:   batch 1, batch 2, or both — the new experimental rounds.
# Eval slice:    20% holdout from each batch, combined. Same eval slice for
#                every scenario so metrics are directly comparable.
#
# Scenarios:
#   O   : original-only (weak prior) — the deployable model before field data
#   O+1 : O posterior as prior, updated on batch-1 train
#   O+2 : O posterior as prior, updated on batch-2 train
#   O+B : O posterior as prior, updated on batch-1+batch-2 train
#   W+B : weak prior, fit on batch-1+batch-2 train (control: how much does
#         the original prior actually contribute vs. having field data alone)
# ============================================================================

# ----- Augmented (field) data -----
aug = read.csv("combined_augment_with_batch.csv") |>
  as_tibble() |>
  mutate(target_vector = factor(target_vector),
         batch         = factor(batch)) |>
  select(-cAb13)  # self-pair excluded throughout the pipeline

aug_cab_cols = grep("^cAb", names(aug), value = TRUE)

# ----- Original data: extract dAb13 panel, drop self-pair, align names -----
orig_raw = read.csv("../data/covid_screen_data_cleaned.csv", check.names = FALSE) |>
  as_tibble()

# Original Replicates is e.g. "R1_WT" — keep only the dAb13_cAbX columns the
# augmented panel uses, then strip the dAb13_ prefix so column names align
# with the augmented frame.
keep_orig = c("Conc.", "target_vector", paste0("dAb13_", aug_cab_cols))
orig = orig_raw |>
  select(all_of(keep_orig)) |>
  rename_with(~ sub("^dAb13_", "", .x), starts_with("dAb13_")) |>
  mutate(target_vector = factor(target_vector, levels = levels(aug$target_vector)))

cat("Original rows:", nrow(orig), "  cols:", ncol(orig), "\n")
cat("Augmented rows:", nrow(aug), "\n")

stopifnot(setequal(intersect(names(orig), names(aug)),
                   c("Conc.", "target_vector", aug_cab_cols)))

# ----- Threshold -----
conc_threshold = 0.9

orig_f = orig |> filter(Conc. > conc_threshold)
aug_f  = aug  |> filter(Conc. > conc_threshold) |> select(-Replicate)

cat("After threshold ", conc_threshold, ":\n",
    "  orig:  ", nrow(orig_f), "\n",
    "  batch1:", sum(aug_f$batch == "1"), "\n",
    "  batch2:", sum(aug_f$batch == "2"), "\n", sep = "")

cab_cols = aug_cab_cols
classes  = levels(aug_f$target_vector)

# ----- Standardize using ORIGINAL stats so the prior scale is meaningful -----
mu_o = colMeans(orig_f[cab_cols])
sd_o = apply(orig_f[cab_cols], 2, sd)
sd_o[sd_o == 0] = 1

scale_x = function(df) {
  z = sweep(as.matrix(df[cab_cols]), 2, mu_o, "-")
  z = sweep(z, 2, sd_o, "/")
  out = as_tibble(z)
  out$target_vector = df$target_vector
  if ("batch" %in% names(df)) out$batch = df$batch
  if ("Conc." %in% names(df)) out$Conc. = df$Conc.
  out
}

orig_z = scale_x(orig_f)
aug_z  = scale_x(aug_f)

# ----- Hold out 20% from each batch as common eval slice -----
b1_z = aug_z |> filter(batch == "1")
b2_z = aug_z |> filter(batch == "2")

set.seed(2024)
b1_split = initial_split(b1_z, prop = 0.80, strata = target_vector)
b2_split = initial_split(b2_z, prop = 0.80, strata = target_vector)

b1_train = training(b1_split); b1_eval = testing(b1_split)
b2_train = training(b2_split); b2_eval = testing(b2_split)

eval_set = bind_rows(b1_eval, b2_eval)
cat("\nEval set:", nrow(eval_set), "rows (batch1:", sum(eval_set$batch == "1"),
    " batch2:", sum(eval_set$batch == "2"), ")\n")

# ----- One-vs-rest fitter -----
ovr_target = function(df, cls) df |> mutate(y = as.integer(target_vector == cls))

fit_binary = function(df, cls,
                      prior_coefs = NULL, prior_sds = NULL,
                      prior_intercept_mean = NULL, prior_intercept_sd = NULL,
                      iter = 2000, chains = 2, seed = 2024) {
  d   = ovr_target(df, cls)
  fml = as.formula(paste("y ~", paste(cab_cols, collapse = " + ")))

  if (is.null(prior_coefs)) {
    pri     = normal(0, 2.5, autoscale = FALSE)
    pri_int = normal(0, 5,   autoscale = FALSE)
  } else {
    pri     = normal(prior_coefs, prior_sds, autoscale = FALSE)
    pri_int = normal(prior_intercept_mean, prior_intercept_sd, autoscale = FALSE)
  }

  stan_glm(fml, data = d, family = binomial(link = "logit"),
           prior = pri, prior_intercept = pri_int,
           iter = iter, chains = chains, refresh = 0, seed = seed)
}

extract_post = function(fit) {
  draws = as.matrix(fit)
  cols  = c("(Intercept)", cab_cols)
  list(mean = colMeans(draws[, cols]),
       sd   = apply(draws[, cols], 2, sd))
}

predict_ovr = function(fits, newdata) {
  ps = sapply(fits, function(f) colMeans(posterior_epred(f, newdata = newdata)))
  ps = ps / rowSums(ps)
  colnames(ps) = names(fits)
  ps
}

evaluate_ovr = function(fits, newdata, label) {
  ps    = predict_ovr(fits, newdata)
  pred  = factor(colnames(ps)[max.col(ps)], levels = classes)
  truth = newdata$target_vector

  prob_df = as_tibble(ps); names(prob_df) = paste0(".pred_", names(prob_df))
  prob_cols = paste0(".pred_", classes)

  res = tibble(truth = truth, pred = pred,
               Conc. = newdata$Conc., batch = newdata$batch) |>
    bind_cols(prob_df)

  list(
    scenario = label,
    metrics  = tibble(
      scenario    = label,
      accuracy    = accuracy_vec(truth, pred),
      kap         = kap_vec(truth, pred),
      mn_log_loss = mn_log_loss_vec(truth, as.matrix(res[prob_cols])),
      roc_auc     = roc_auc_vec(truth, as.matrix(res[prob_cols]),
                                estimator = "hand_till")
    ),
    preds = mutate(res, scenario = label)
  )
}

# Build update_set helper
update_with_prior = function(update_data, prior_posts, label) {
  cat("\n[", label, "] fitting on n =", nrow(update_data), "...\n", sep = "")
  fits = setNames(lapply(classes, function(cls) {
    cat("   ", cls, "\n")
    p = prior_posts[[cls]]
    fit_binary(update_data, cls,
               prior_coefs          = p$mean[cab_cols],
               prior_sds            = pmax(p$sd[cab_cols], 1e-3),
               prior_intercept_mean = p$mean["(Intercept)"],
               prior_intercept_sd   = pmax(p$sd["(Intercept)"], 1e-3))
  }), classes)
  fits
}

# ============================================================================
# O: original-only (weak prior) — the deployable starting model
# ============================================================================
cat("\n[O] Fitting original-only (weak prior), 4 OVR models...\n")
fits_O = setNames(lapply(classes, function(cls) {
  cat("   ", cls, "\n")
  fit_binary(orig_z, cls)
}), classes)
posts_O = lapply(fits_O, extract_post)

# ============================================================================
# O+1, O+2, O+B: sequential updates from O posterior
# ============================================================================
fits_O1 = update_with_prior(b1_train,                          posts_O, "O+1")
fits_O2 = update_with_prior(b2_train,                          posts_O, "O+2")
fits_OB = update_with_prior(bind_rows(b1_train, b2_train),     posts_O, "O+B")

# ============================================================================
# W+B: control — weak prior, fit on both batches' train (no transfer from O)
# ============================================================================
cat("\n[W+B] Control — weak prior on batch1+batch2 train...\n")
fits_WB = setNames(lapply(classes, function(cls) {
  cat("   ", cls, "\n")
  fit_binary(bind_rows(b1_train, b2_train), cls)
}), classes)

# ----- Evaluate -----
cat("\nEvaluating on shared eval slice (n =", nrow(eval_set), ")...\n")

eval_O  = evaluate_ovr(fits_O,  eval_set, "O: original only")
eval_O1 = evaluate_ovr(fits_O1, eval_set, "O+1: original prior + batch1")
eval_O2 = evaluate_ovr(fits_O2, eval_set, "O+2: original prior + batch2")
eval_OB = evaluate_ovr(fits_OB, eval_set, "O+B: original prior + both batches")
eval_WB = evaluate_ovr(fits_WB, eval_set, "W+B: weak prior + both batches")

scenarios = list(eval_O, eval_O1, eval_O2, eval_OB, eval_WB)
all_metrics = bind_rows(lapply(scenarios, `[[`, "metrics"))
all_preds   = bind_rows(lapply(scenarios, `[[`, "preds"))

cat("\n===== Performance on combined eval slice =====\n")
print(all_metrics)

# Per-batch metric breakdown
per_batch_metrics = all_preds |>
  group_by(scenario, batch) |>
  summarise(
    n        = n(),
    accuracy = mean(pred == truth),
    .groups  = "drop"
  )
cat("\n===== Accuracy by batch =====\n")
print(per_batch_metrics)

# ----- Per-concentration accuracy by scenario -----
per_conc = all_preds |>
  mutate(Conc. = round(Conc., 4),
         correct = as.integer(pred == truth)) |>
  group_by(scenario, Conc.) |>
  summarise(num_correct     = sum(correct),
            total           = n(),
            percent_correct = round(num_correct / total * 100, 2),
            .groups         = "drop") |>
  arrange(scenario, Conc.)

cat("\n===== Per-concentration accuracy by scenario =====\n")
print(per_conc, n = Inf)

# Color palette
scen_colors = c(
  "O: original only"                   = "#7BA7AB",
  "O+1: original prior + batch1"       = "#266B6E",
  "O+2: original prior + batch2"       = "#1F4E50",
  "O+B: original prior + both batches" = "#0E2E2F",
  "W+B: weak prior + both batches"     = "#B5403F"
)
scen_order = names(scen_colors)

# ----- Metric chart -----
p_metrics = all_metrics |>
  mutate(scenario = factor(scenario, levels = scen_order)) |>
  pivot_longer(c(accuracy, kap, mn_log_loss, roc_auc),
               names_to = "metric", values_to = "value") |>
  ggplot(aes(scenario, value, fill = scenario)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_text(aes(label = round(value, 3)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = scen_colors) +
  labs(title = "Original screen as prior — sequential updates",
       subtitle = paste0("Combined held-out slice (b1+b2) | Conc. > ",
                         conc_threshold, " ng/ml"),
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x     = element_text(angle = 25, hjust = 1))
ggsave(file.path(out_dir, "metrics.png"), p_metrics,
       width = 11, height = 7, dpi = 150)

# ----- Per-conc bar chart -----
p_per_conc = per_conc |>
  mutate(scenario = factor(scenario, levels = scen_order),
         Conc.    = factor(Conc.)) |>
  ggplot(aes(Conc., percent_correct, fill = scenario)) +
  geom_col() +
  geom_text(aes(label = percent_correct), vjust = 1.5, color = "white", size = 2.7) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(values = scen_colors) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                     limits = c(0, 110)) +
  labs(title    = "Accuracy by concentration — original prior + updates",
       subtitle = paste0("Combined held-out slice (b1+b2) | Conc. > ",
                         conc_threshold, " ng/ml"),
       x        = "Concentration (ng/ml)",
       y        = "Percent Correct") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "percent_correct_by_conc.png"), p_per_conc,
       width = 9, height = 11, dpi = 150)

# ----- Coefficient shift: O → O+B per class -----
shift_df = bind_rows(lapply(classes, function(cls) {
  bind_rows(
    as_tibble(as.matrix(fits_O[[cls]])[, cab_cols]) |>
      pivot_longer(everything(), names_to = "cAb", values_to = "coef") |>
      mutate(scenario = "O: original only"),
    as_tibble(as.matrix(fits_OB[[cls]])[, cab_cols]) |>
      pivot_longer(everything(), names_to = "cAb", values_to = "coef") |>
      mutate(scenario = "O+B: original prior + both batches")
  ) |> mutate(class = cls)
}))

if (requireNamespace("ggridges", quietly = TRUE)) {
  p_shift = shift_df |>
    ggplot(aes(coef, cAb, fill = scenario)) +
    ggridges::geom_density_ridges(alpha = 0.6, scale = 0.9) +
    facet_wrap(~ class) +
    scale_fill_manual(values = c("O: original only" = "#7BA7AB",
                                 "O+B: original prior + both batches" = "#0E2E2F")) +
    labs(title = "Coefficient posteriors: original prior vs. fully-updated",
         x = "Coefficient (log-odds, standardized predictors)", y = NULL) +
    theme_minimal()
  ggsave(file.path(out_dir, "coef_shift.png"), p_shift,
         width = 11, height = 8, dpi = 150)
}

# ----- Save tables -----
write.csv(all_metrics,        file.path(out_dir, "metrics.csv"),       row.names = FALSE)
write.csv(per_conc,           file.path(out_dir, "per_conc.csv"),      row.names = FALSE)
write.csv(per_batch_metrics,  file.path(out_dir, "per_batch.csv"),     row.names = FALSE)

cat("\nWrote outputs to ", out_dir, "/\n", sep = "")
