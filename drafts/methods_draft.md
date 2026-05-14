# Methods (first draft)

> Items flagged with `[FILL]` need bench-side or citation specifics that aren't recoverable from the code.

## 2.1 Sandwich immunoassay and screening design

Recombinant SARS-CoV-2 spike protein from four variants — wild-type (Wuhan-Hu-1), Beta, Delta, and Omicron — was probed against an arrayed library of capture antibodies (cAbs) and detection antibodies (dAbs) in a fluorescence-based sandwich format. Each slide combinatorially paired a fixed dAb with a panel of cAbs, producing per-spot fluorescence intensities for every dAb–cAb pair. Variants were measured at eight spike-protein concentrations spanning four orders of magnitude (10, 1, 0.5, 0.25, 0.1, 0.05, 0.01, and 0 µg/mL), in three replicate spots per condition. The full screen comprises ~13 dAbs × ~29 cAbs ≈ 400 unique antibody pairs, evaluated under 8 (concentration) × 4 (variant) × 3 (replicate) = 96 conditions per dAb panel.

## 2.2 Image acquisition and per-spot quantification

Slides were scanned on a `[FILL: GenePix scanner model]` and analyzed with GenePix Pro to extract per-spot mean foreground (F635.Mean) and median local background (B635.Mean) at 635 nm. All downstream analysis was performed in R (v4.1.2) on the resulting GenePix Results (.gpr) tables.

For each spot we computed a background-subtracted intensity *I* = F635.Mean − B635.Mean, with negative or sub-unity values floored at 1 to permit subsequent log-transformation. To suppress hybridization artifacts and dust-attributable spikes within technical replicates, we applied a Tukey-style fence at three within-condition standard deviations: any replicate exceeding `mean ± 3·SD` of its sibling spots was replaced with the within-condition mean. The three QC-cleaned replicates were then averaged to a single intensity per dAb–cAb–concentration–variant condition, log10-transformed, and reshaped from long to wide format with one column per dAb–cAb pair (named `dabX_cabY`) and one row per (variant × concentration × replicate) observation. Variant identity was extracted from the slide label and recorded as the four-level outcome `target_vector`.

## 2.3 Pre-screen filtering of non-responsive pairs

Pairs with no dose response cannot contribute to variant discrimination and contribute only noise to downstream feature-selection algorithms. Each pair was therefore screened against a per-variant dose-response criterion: a pair was retained if, for *any* variant, its log10 intensity fell below 3 at the lowest concentrations *and* exceeded 4 at the highest concentrations. This corresponds to a base-noise ceiling of ~10³ relative fluorescence units and a dynamic-range floor of ~10⁴ — empirical thresholds for the assay platform. Pairs failing this check across all four variants were removed before any model fitting.

## 2.4 Per-dAb partitioning

Because the deployed device is constrained to a single dAb per chip for practical reasons (reagent cost, optical layout, washing protocol), all feature selection and modeling was performed within each dAb separately rather than over the joint pair space. The filtered table was therefore split by the `dab` prefix of each column name into 13 per-dAb subsets, each retaining only that dAb's cAb pairs together with the shared `Conc.`, `Replicates`, and `target_vector` columns. All subsequent analyses operate on each per-dAb subset independently and are aggregated only at reporting time.

## 2.5 Three-level nested holdout cross-validation scheme

A central methodological challenge in this dataset is that *n* per dAb is small (≤ 96 observations across 4 classes) while *p* is moderate (~29 features per dAb), and we must simultaneously estimate three quantities: a hyperparameter (e.g., LASSO penalty, RF `mtry`), the optimal number of antibody pairs to retain *N\**, and an unbiased estimate of generalization performance. Naïve combination of these on a single split risks selection bias of well-documented severity (Cawley & Talbot, 2010). We therefore adopted a three-level nested holdout scheme, applied identically across all six feature-selection methods (Section 2.6):

1. **Final test holdout (Level 1).** A 20% stratified holdout (`final_test`) was withheld from the full per-dAb table at the start of the run and not touched again until the final unbiased evaluation step.

2. **Feature-count validation holdout (Level 2).** Of the remaining 80% (`training_data`), a further stratified 20% was held out as a `feature_val_test` set used *only* to choose the optimal feature count *N\**. The remainder (`feature_val_train`, ≈ 64% of the original data) was used for hyperparameter tuning and feature ranking.

3. **Hyperparameter tuning (Level 3).** When the per-method default permitted (e.g., `cv.glmnet`, ranger OOB), tuning used the algorithm's internal cross-validation on `feature_val_train` directly. When external *k*-fold CV was used (k = 4–5), it was strata-balanced and confined to `feature_val_train`.

After the best hyperparameters were locked, features were ranked on the full `feature_val_train` at those settings; subsets of the top *N* features (with *N* drawn from {3, 4, 5, 6, 7, 8, 10} and capped at the per-method ceiling — number of nonzero LASSO coefficients, number of Boruta-confirmed pairs, etc.) were trained on `feature_val_train` and scored on `feature_val_test`. *N\** was selected as the smallest *N* maximizing accuracy on `feature_val_test`. The model was then re-trained on the full `training_data` (i.e. `feature_val_train` ∪ `feature_val_test`) at the locked hyperparameters and *N\** features — this re-merge is leakage-free because all selection decisions are now fixed — and finally evaluated *once* on the original `final_test`. To support engineering trade-off analysis, performance curves over *N* were additionally produced on the true `final_test` holdout by the same retrain-on-`training_data` procedure.

Stratification on `target_vector` was enforced at every split via `rsample::initial_split(strata = …)`.

## 2.6 Feature selection methods

Six feature-selection methods were embedded in the scheme above so that comparisons reflect intrinsic algorithmic differences rather than CV-protocol differences. All methods are implemented in `tidymodels` workflows. Outcomes are four-class (WT, Beta, Delta, Omicron).

**LASSO (multinomial logistic regression with L1 penalty).** Predictors were Box–Cox transformed and z-standardized within each recipe; the engine's internal standardization was disabled (`standardize = FALSE`) to prevent double-scaling. Multiclass support was via grouped multinomial regression (`glmnet::cv.glmnet(family = "multinomial", type.multinomial = "grouped")`) with 1-SE rule penalty selection. Pair importance was the mean absolute coefficient across the four variant classes. (`Embedded_FS/lasso_tidymodels.R`)

**Random forest (ranger backend).** Forests of 500 trees were grown without input transformation (RF is invariant to monotone transforms; only zero-variance predictors were removed). When external CV was disabled, `mtry` was selected on out-of-bag error over a grid of {2, ⌊√p⌋, ⌊p/3⌋, ⌊p/2⌋} with `min_n` fixed at 5; otherwise `mtry` × `min_n` were jointly tuned. Pair importance was the permutation importance averaged across trees. When evaluating subsets of *N* < `mtry`, `mtry` was capped at *N*. (`Embedded_FS/rf_tidymodels.R`)

**MARS (multivariate adaptive regression splines, multiclass via flexible discriminant analysis).** Implemented through `discrim::discrim_flexible()` with the `earth` engine — required because `parsnip::mars()` is binary-only, while `mda::fda` wraps `earth` in a multinomial flexible-discriminant framework. Hyperparameters `num_terms` ∈ {5, 10, 15, 20} and `prod_degree` ∈ {1, 2} were tuned by 5-fold CV; importance was earth's GCV-based variable importance, with features never selected receiving zero importance. (`Embedded_FS/mars_tidymodels.R`)

**PLS-DA (partial least squares discriminant analysis).** Multiclass PLS-DA via `plsmod::pls(engine = "mixOmics")`. Number of components was tuned in {2, 3, 5, 8, 10} (capped at min(*p*, *n*−1)) by 5-fold CV. Pair importance was the cumulative Variable Importance in Projection (VIP) score across all selected components. (`Embedded_FS/pls_tidymodels.R`)

**mRMR (minimum redundancy, maximum relevance).** A model-free filter method (`mRMRe::mRMR.classic`) that greedily selects features maximizing mutual information with the target while penalizing redundancy with already-selected features. Targets were integer-coded for the MI estimator; ranking is therefore independent of any classifier. Feature subsets were subsequently scored using a default RF (sqrt(*p*) `mtry`, 500 trees) for comparability with the other methods (Peng et al., 2005; De Jay et al., 2013). (`Embedded_FS/mrmr_tidymodels.R`)

**Boruta (all-relevant wrapper).** A shadow-feature significance test atop random forest: each feature is iteratively compared against the highest importance among permutation-shuffled "shadow" copies, accepting features whose importance significantly exceeds the shadow ceiling. We used `[FILL: maxRuns]` iterations at p < `[FILL: 0.01]`, resolved tentative features with `TentativeRoughFix`, and ranked confirmed features by mean shadow-relative importance. The downstream evaluation classifier was the same RF protocol as the standalone RF method. (`Embedded_FS/boruta_tidymodels.R`)

## 2.7 Multi-seed replication and uncertainty

To quantify run-to-run variability arising from the small data partitions, every (method × dAb) pipeline was replicated under three independent random seeds (1969, 1986, 2024), each producing fresh Level-1, Level-2, and tuning splits. Reported point estimates are the across-seed mean; uncertainty bands are ±1 standard error of the mean across the three seeds (or ±1 SE across the three-seeds × thirteen-dAbs pool when aggregating). Final-test confusion matrices were summed across seeds for per-class diagnostics.

## 2.8 Composite ranking across feature-selection methods

To identify antibody pairs that are robust to method choice, per-dAb feature rankings from each method were standardized to a common (dAb, feature, importance) schema and combined via Borda count aggregation (`Embedded_FS/composite_score.R`). Within each (dAb, method), features were rank-ordered from most to least important (rank 1 = most important; tied importances received the average of the tied ranks). Ranks were summed across methods to produce a per-(dAb, feature) composite score, lower scores indicating features consistently prioritized across methods. Method-specific normalizations: Boruta-rejected features were excluded prior to ranking; mRMR ranks (already ordinal) were negated to align with the higher-is-better convention; multi-seed rankings were averaged before Borda assignment.

## 2.9 Performance metrics

All models were assessed on three metrics: classification accuracy, Cohen's κ (which corrects accuracy for chance agreement under class imbalance), and macro-averaged F1 score (the mean of per-class F1, which weights minority and majority classes equally). The feature-count selection step used accuracy alone; reporting and figures use all three. `[FILL: any specific cite for kappa or macro-F1 if your committee will want it]`

## 2.10 Software and reproducibility

All analyses were performed in R 4.1.2 using the `tidymodels` framework (v `[FILL]`). Engine packages: `glmnet` (LASSO), `ranger` (random forest, Boruta evaluation), `earth` and `mda` via `discrim` (MARS/FDA), `mixOmics` via `plsmod` (PLS-DA), `mRMRe` (mRMR), and `Boruta`. Auxiliary packages: `tidyverse` for data wrangling, `psych` and `sda` for preliminary screen filtering, and `openxlsx` for results export. All seeds, split proportions, and hyperparameter grids are documented in the per-method scripts under `Embedded_FS/` in the project repository.

---

## Editing notes

- **§2.1** needs the bench-side specifics (instrument vendor, slide chemistry, blocking buffer, incubation times).
- **§2.2** sanity-check that the QC fence is `mean ± 3·SD` per (Slide × Block × Row) group at the bench-protocol level — that's what `tukey_fence = 3` does in code.
- **§2.6 (Boruta)** maxRuns and pValue are placeholders; fill in actual values used.
- **§2.8** drop if Borda composite isn't included in dissertation results.
- **Augmented data (Phase 2)** and the **stacking ensemble** experiments aren't included; fold in if needed.
