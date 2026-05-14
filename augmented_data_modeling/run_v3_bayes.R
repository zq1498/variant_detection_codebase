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
dir.create("v3_outputs", showWarnings = FALSE)

# ----- Data -----
raw = read.csv("combined_augment_with_batch.csv") |>
  as_tibble() |>
  mutate(target_vector = factor(target_vector),
         batch         = factor(batch)) |>
  select(-cAb13)

conc_threshold = 0.9   # same threshold as the RF run

dat = raw |>
  filter(Conc. > conc_threshold) |>
  select(-Replicate)
# Conc. is kept on the frame so we can group test predictions by it later.
# It is NOT used as a predictor — the formula in fit_binary references cab_cols
# explicitly.

cab_cols = grep("^cAb", names(dat), value = TRUE)
classes  = levels(dat$target_vector)

cat("rows after threshold:", nrow(dat), "\n")
cat("by batch:\n"); print(table(dat$batch))
cat("by class:\n"); print(table(dat$target_vector))

# ----- Standardize predictors using batch-1 statistics -----
# This makes the Normal(0, 2.5) prior scale meaningful AND keeps the same scale
# across both batches so a posterior from batch 1 can be used as a prior for
# batch 2 without rescaling.
b1 = dat |> filter(batch == "1")
b2 = dat |> filter(batch == "2")

mu_b1 = colMeans(b1[cab_cols])
sd_b1 = apply(b1[cab_cols], 2, sd)
sd_b1[sd_b1 == 0] = 1

scale_x = function(df) {
  z = sweep(as.matrix(df[cab_cols]), 2, mu_b1, "-")
  z = sweep(z, 2, sd_b1, "/")
  as_tibble(z) |> mutate(target_vector = df$target_vector,
                         batch         = df$batch,
                         Conc.         = df$Conc.)
}

b1_z = scale_x(b1)
b2_z = scale_x(b2)

# ----- Split: hold out 25% of batch 2 as the common eval set -----
set.seed(2024)
b2_split = initial_split(b2_z, prop = 0.75, strata = target_vector)
b2_train = training(b2_split)
b2_eval  = testing(b2_split)

cat("\nbatch-2 train:", nrow(b2_train), " batch-2 eval:", nrow(b2_eval), "\n")

# ----- One-vs-rest helper -----
ovr_target = function(df, cls) {
  df |> mutate(y = as.integer(target_vector == cls))
}

# Fit a one-vs-rest stan_glm for one class.  prior_coefs/prior_sds, when
# supplied, override the default weakly-informative prior with the
# (potentially discounted) batch-1 posterior.
fit_binary = function(df, cls, prior_coefs = NULL, prior_sds = NULL,
                      prior_intercept_mean = NULL, prior_intercept_sd = NULL,
                      iter = 2000, chains = 2, seed = 2024) {
  d   = ovr_target(df, cls)
  fml = as.formula(paste("y ~", paste(cab_cols, collapse = " + ")))

  if (is.null(prior_coefs)) {
    pri = normal(0, 2.5, autoscale = FALSE)
    pri_int = normal(0, 5, autoscale = FALSE)
  } else {
    pri = normal(prior_coefs, prior_sds, autoscale = FALSE)
    pri_int = normal(prior_intercept_mean, prior_intercept_sd, autoscale = FALSE)
  }

  stan_glm(fml, data = d, family = binomial(link = "logit"),
           prior = pri, prior_intercept = pri_int,
           iter = iter, chains = chains, refresh = 0, seed = seed)
}

# Extract posterior mean and sd for the 11 cAb coefs + intercept
extract_post = function(fit) {
  draws = as.matrix(fit)
  cols  = c("(Intercept)", cab_cols)
  list(mean = colMeans(draws[, cols]),
       sd   = apply(draws[, cols], 2, sd))
}

# Predict P(class | x) by averaging posterior_epred over draws (one-vs-rest).
# Returns matrix [n_rows × n_classes].  Rows are renormalised so probs sum to 1
# (one-vs-rest fits don't share a softmax denominator).
predict_ovr = function(fits, newdata) {
  ps = sapply(fits, function(f) colMeans(posterior_epred(f, newdata = newdata)))
  ps = ps / rowSums(ps)
  colnames(ps) = names(fits)
  ps
}

evaluate_ovr = function(fits, newdata, label) {
  ps = predict_ovr(fits, newdata)
  pred = factor(colnames(ps)[max.col(ps)], levels = classes)
  truth = newdata$target_vector

  prob_df = as_tibble(ps)
  names(prob_df) = paste0(".pred_", names(prob_df))

  res = tibble(truth = truth, pred = pred, Conc. = newdata$Conc.) |>
    bind_cols(prob_df)

  prob_cols = paste0(".pred_", classes)

  acc = accuracy_vec(res$truth, res$pred)
  kap = kap_vec(res$truth, res$pred)
  ll  = mn_log_loss_vec(res$truth, as.matrix(res[prob_cols]))
  auc = roc_auc_vec(res$truth, as.matrix(res[prob_cols]),
                    estimator = "hand_till")

  list(scenario = label,
       metrics  = tibble(scenario = label,
                         accuracy = acc, kap = kap,
                         mn_log_loss = ll, roc_auc = auc),
       preds    = mutate(res, scenario = label))
}

# ============================================================================
# Scenario A: weak prior, fit on batch 1 only
# ============================================================================
cat("\n[A] Fitting batch-1-only baseline (weak priors), 4 OVR models...\n")
fits_A = setNames(lapply(classes, function(cls) {
  cat("   ", cls, "\n")
  fit_binary(b1_z, cls)
}), classes)

posts_A = lapply(fits_A, extract_post)

# ============================================================================
# Scenario B: weak prior, fit on batch-2 train only (no prior-knowledge transfer)
# ============================================================================
cat("\n[B] Fitting batch-2-only baseline (weak priors), 4 OVR models...\n")
fits_B = setNames(lapply(classes, function(cls) {
  cat("   ", cls, "\n")
  fit_binary(b2_train, cls)
}), classes)

# ============================================================================
# Scenario C: sequential update — batch-2 train with batch-1 posterior as prior
#   discount = 1.0 means trust batch-1 fully.  Lower values inflate prior sd
#   (less confidence in batch-1 transfer due to drift).
# ============================================================================
discount = 1.0
cat(sprintf("\n[C] Sequential update on batch-2 train, prior = batch-1 posterior (discount = %.2f)...\n", discount))
fits_C = setNames(lapply(classes, function(cls) {
  cat("   ", cls, "\n")
  p = posts_A[[cls]]
  fit_binary(b2_train, cls,
             prior_coefs          = p$mean[cab_cols],
             prior_sds            = pmax(p$sd[cab_cols] / discount, 1e-3),
             prior_intercept_mean = p$mean["(Intercept)"],
             prior_intercept_sd   = pmax(p$sd["(Intercept)"] / discount, 1e-3))
}), classes)

# ----- Evaluate all three on the common batch-2 eval slice -----
cat("\nEvaluating on shared batch-2 eval set (n =", nrow(b2_eval), ")...\n")
eval_A = evaluate_ovr(fits_A, b2_eval, "A: batch1 only (weak prior)")
eval_B = evaluate_ovr(fits_B, b2_eval, "B: batch2-train only (weak prior)")
eval_C = evaluate_ovr(fits_C, b2_eval, "C: batch2-train + batch1 prior")

all_metrics = bind_rows(eval_A$metrics, eval_B$metrics, eval_C$metrics)
cat("\n===== Bayesian one-vs-rest performance (held-out batch-2 eval) =====\n")
print(all_metrics)

# ----- Coefficient shift visualisation: A → C (selected class) -----
coef_shift = function(cls) {
  draws_A = as.matrix(fits_A[[cls]])[, cab_cols] |> as_tibble() |>
    pivot_longer(everything(), names_to = "cAb", values_to = "coef") |>
    mutate(scenario = "A: batch1 only")
  draws_C = as.matrix(fits_C[[cls]])[, cab_cols] |> as_tibble() |>
    pivot_longer(everything(), names_to = "cAb", values_to = "coef") |>
    mutate(scenario = "C: A-prior + batch2 update")
  bind_rows(draws_A, draws_C) |> mutate(class = cls)
}

shift_df = bind_rows(lapply(classes, coef_shift))

p_shift = shift_df |>
  ggplot(aes(coef, cAb, fill = scenario)) +
  ggridges::geom_density_ridges(alpha = 0.6, scale = 0.9) +
  facet_wrap(~ class) +
  scale_fill_manual(values = c("A: batch1 only" = "#7BA7AB",
                               "C: A-prior + batch2 update" = "#0e1372")) +
  labs(title = "Coefficient posteriors: prior (A) vs sequentially-updated (C)",
       x = "Coefficient (log-odds units, standardized predictors)", y = NULL) +
  theme_minimal()

# ggridges is optional — fall back if missing
if (!requireNamespace("ggridges", quietly = TRUE)) {
  p_shift = shift_df |>
    ggplot(aes(coef, fill = scenario)) +
    geom_density(alpha = 0.5) +
    facet_grid(cAb ~ class, scales = "free_y") +
    scale_fill_manual(values = c("A: batch1 only" = "#7BA7AB",
                                 "C: A-prior + batch2 update" = "#266B6E")) +
    labs(title = "Coefficient posteriors: prior (A) vs sequentially-updated (C)",
         x = "Coefficient (log-odds units)", y = NULL) +
    theme_minimal()
}

ggsave("v3_outputs/v3_bayes_coef_shift.png", p_shift, width = 11, height = 8, dpi = 150)

# ----- Side-by-side metric chart -----
p_metrics = all_metrics |>
  pivot_longer(c(accuracy, kap, mn_log_loss, roc_auc),
               names_to = "metric", values_to = "value") |>
  ggplot(aes(scenario, value, fill = scenario)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  geom_text(aes(label = round(value, 3)), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("A: batch1 only (weak prior)"      = "#7BA7AB",
                               "B: batch2-train only (weak prior)" = "#B5403F",
                               "C: batch2-train + batch1 prior"    = "#266B6E")) +
  labs(title = "Bayesian sequential update vs. baselines",
       subtitle = paste0("Held-out batch-2 eval | Conc. > ", conc_threshold, " ng/ml"),
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("v3_outputs/v3_bayes_metrics.png", p_metrics, width = 9, height = 6, dpi = 150)

# ----- Per-concentration accuracy by scenario -----
all_preds = bind_rows(eval_A$preds, eval_B$preds, eval_C$preds)

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

p_per_conc = per_conc |>
  mutate(Conc. = factor(Conc.)) |>
  ggplot(aes(Conc., percent_correct, fill = scenario)) +
  geom_col() +
  geom_text(aes(label = percent_correct), vjust = 1.5, color = "white", size = 3) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(values = c("A: batch1 only (weak prior)"      = "#7BA7AB",
                               "B: batch2-train only (weak prior)" = "#B5403F",
                               "C: batch2-train + batch1 prior"    = "#266B6E")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                     limits = c(0, 110)) +
  labs(title    = "Bayesian models — accuracy by concentration",
       subtitle = paste0("Held-out batch-2 eval | Conc. > ", conc_threshold, " ng/ml"),
       x        = "Concentration (ng/ml)",
       y        = "Percent Correct") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("v3_outputs/v3_bayes_percent_correct_by_conc.png", p_per_conc,
       width = 9, height = 8, dpi = 150)

# ----- Save metric tables -----
write.csv(all_metrics, "v3_outputs/v3_bayes_metrics.csv", row.names = FALSE)
write.csv(per_conc,    "v3_outputs/v3_bayes_per_conc.csv", row.names = FALSE)

cat("\nWrote v3_outputs/v3_bayes_metrics.png, v3_bayes_coef_shift.png,",
    "v3_bayes_percent_correct_by_conc.png, v3_bayes_metrics.csv, v3_bayes_per_conc.csv\n")
