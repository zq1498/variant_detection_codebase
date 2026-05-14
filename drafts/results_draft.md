# Results (first draft)

> **Drafting notes (delete before final):**
> - Two throughlines to keep visible across sections:
>   1. **dAb18 / dAb19 narrative.** These are the candidates carried into the augmented-data experiments; introduce them by name as early as 2.2.3, keep them visible in the FS results (2.2.6–2.2.8), and close the loop in 2.2.9–2.2.12.
>   2. **Metrics matter for porting in silico → experimental performance.** Set up in 2.2.5 (why accuracy alone undersells the 4-class problem; why kappa and macro-F1 are co-reported), reinforce in 2.2.6 (how feature-count curves look different in accuracy vs kappa near the boundary), and close in 2.2.12 (SVM-RBF matches RF on accuracy but loses badly on log loss / ROC AUC — the metric you pick decides whether the model is "deployable").
> - "Show diverse results, not exhaustive" — pick representative dAbs for the FS plots (one top performer, one mid, one weak, plus dAb18/19). Full per-dAb tables can live in supplementary.

---

## 2.2.1 Problem statement and design constraints

- **Task.** Classify SARS-CoV-2 spike-protein variants (WT, Beta, Delta, Omicron) from a fluorescence-based sandwich immunoassay using the smallest antibody panel that still discriminates reliably.
- **Device-imposed constraints that shape every modeling choice that follows:**
  - One detection antibody (dAb) per device → models are built and evaluated **per-dAb**, not pooled across dAbs.
  - Multiplex capacity on the capture (cAb) side is small → the deployable panel is 3–10 dAb–cAb pairs.
  - Inference happens on a single sample at a time → no batch-level rescaling is available at deployment, which constrains the preprocessing pipeline (see 2.2.4) and the calibration discussion (2.2.12).
- **Operational definition of "success."** A pair subset is a candidate for downstream wet-lab validation only if (a) it classifies above chance with margin on held-out CV folds, (b) the gain over a smaller subset is meaningful in *kappa* / *macro-F1*, not just accuracy, and (c) the result is stable across reshuffles of the data (multi-seed).
- **Why this is non-trivial in silico, and why accuracy alone is the wrong primary readout.** 4 classes × 8 concentration levels × 3 replicates per dAb gives only ~96 observations per dAb (≈ 20 per held-out CV fold). Three structural problems with reporting accuracy as the headline metric in this regime, ordered by how much they bite on this specific problem:
  1. **Small test folds make accuracy noisy.** At ~20 rows per fold, one flipped prediction shifts accuracy by ~5pp; even at the larger augmented-data test set (n ≈ 57), it is ~1.75pp. A 2–3pp accuracy gap between two models can easily reflect which samples landed in which split rather than a real model difference. We mitigate this with multi-seed reporting (mean ± SE across 3–5 seeds), but the noise is still large enough that a single-seed result is not a reliable rank.
  2. **Chance baseline is 25%, not 50%.** A 4-class problem inflates accuracy values relative to binary classification; "80% accuracy" reads differently here than in most published benchmarks. Cohen's kappa is chance-corrected and so makes "is this model doing work above chance" structurally explicit — this matters even at infinite n.
  3. **Variants are balanced, but error structure is concentration-dependent.** A model can score respectably on accuracy by getting the high-concentration samples right and quietly failing on the low-concentration regime that the device will actually operate in. Macro-F1 (per-class average) and per-concentration accuracy panels surface this; raw accuracy hides it.
- **Throughline setup (metrics).** Together these motivate the metric choices in 2.2.5: we co-report accuracy, kappa, macro-F1, ROC AUC, and log loss because they disagree exactly at the boundary where deployment-readiness is decided, and that disagreement is more informative than any single number.

*Figures/tables needed:* schematic of the assay (existing); device-constraint diagram showing 1 dAb × N cAbs per panel.

---

## 2.2.2. Data set description

- **Panel layout.** ~13 detection antibodies × ~29 capture antibodies = ~400 unique dAb–cAb pair readouts; each dAb is analyzed as an independent dataset of cAb features.
- **Design.** 8 concentration levels × 4 variants × 3 replicates = 96 observations per dAb. Concentration sweep spans roughly [low_conc] – [high_conc] ng/mL (fill in exact range).
- **Sources of variation captured.** Variant (target signal), concentration (nuisance covariate that mostly modulates magnitude, not class identity, above the LoD), replicate (within-run technical variation), and — later, in the augmented data — batch / protocol drift (2.2.11).
- **Files.** Primary cleaned screen: `data/covid_screen_data_cleaned.csv`. Augmented batches (collected against dAb18 and dAb19 specifically) live in `data/augmented_data/`. *(verify path before finalizing)*

*Figures/tables needed:* table listing dAbs, n_cAbs each, concentration grid; sample-count breakdown by variant × concentration.

---

## 2.2.3. Screen results, EDA

- **Dynamic-range overview.** Across dAbs, fluorescence response increases monotonically with antigen concentration above the LoD but flattens or becomes noise-dominated below it — visible in concentration-response curves overlaid by variant.
- **Class separability is dAb-dependent, sharply.** A handful of dAbs separate at least one variant cleanly across most cAbs (dAb13, dAb17); others show essentially overlapping fluorescence distributions across variants at every concentration (dAb4, dAb11). This pre-figures the FS results in 2.2.6. WHERE WAS THIS CLASS SEPARABILITY SUCCESS DRAWN FROM? LIKE LITERALLY WHERE IS THE DATA
- **Where dAb18 and dAb19 sit in the raw EDA.** Both show clear variant-dependent response patterns across a non-trivial fraction of their cAb partners, but with cleaner separability at higher concentrations — establishing why the ≥ 0.9 ng/mL threshold becomes load-bearing in the augmented work (2.2.9).
- **PCA / unsupervised view.** First two PCs on the standardized per-dAb feature blocks show partial variant clustering for the strong dAbs and no obvious structure for the weak ones, consistent with the supervised results that follow.

*Figures/tables needed:*
- Concentration-response curves for one strong dAb (dAb13 or dAb17), one weak (dAb4 or dAb11), and **dAb18 + dAb19** (the latter two highlighted in the caption).
- PCA score plots for the same set.
- Boxplot or heatmap of variant-wise fluorescence dispersion across dAbs.

---

## 2.2.4. Data set preprocessing (brief, can be more detailed in methods)

- **Per-pair normalization.** Box-Cox + centering/scaling at the dAb level for methods that need it (LASSO, MARS, PLS); zero-variance filtering only for RF.
- **Crucially: no preprocessing leakage.** Each CV split divides the data into a *training portion* (used to fit) and a *held-out portion* (used to evaluate / score feature ranks). Preprocessing parameters — Box-Cox λ, feature means, feature standard deviations, kNN imputation neighbours — are estimated **only on the training portion of each split** and then applied to the held-out portion. They are *not* estimated once on the full dataset before splitting. This matters because the standard accidental shortcut (fitting Box-Cox/normalize globally, then splitting) lets the held-out rows influence their own normalization, which inflates evaluation metrics and biases the FS ranks toward features whose scale happened to align with the held-out set. Mechanically: tidymodels `recipes` are baked inside the resample machinery, so the leakage is excluded by construction rather than by discipline.
- **Pair naming fix.** A substring-match bug in the original preprocessor (e.g. dAb1 matching dAb13) is fixed via anchored matching so the pair list is correct.
- **Class label encoding.** Variant labels are 4-level factors; concentration is retained as a covariate where applicable but is *not* a predictor in the variant-classification models.
- Detailed protocol in Methods §X.

*Figures/tables needed:* none in the main text; a preprocessing diagram could go in methods.

---

## 2.2.5. Cross validation + safeguards on overfitting

- **Three-level CV.** Outer split for held-out evaluation; inner repeated k-fold for hyperparameter / feature-count tuning; an innermost layer for ranking-stability when applicable. This is the structural defense against the most common in-silico-vs-experimental gap: optimistic ranks from re-using the same data for FS and evaluation.
- **Multi-seed runs.** Three seeds (1969, 1986, 2024) for the screen-level FS work; five seeds (1969, 1986, 2024, 2025, 2026) for the augmented-data model comparison. We report mean ± SE across seeds, not the best seed. Per-seed variability is itself a result — large SE on a dAb is a soft signal that the dAb is borderline.
- **Metric set (this is the throughline lift).** For every CV fold we record:
  - **Accuracy** — intuitive but uninformative when one class dominates.
  - **Cohen's kappa** — chance-corrected; the metric that punishes "guess the majority class" models hardest.
  - **Macro-F1** — averages F1 across the four variants so a single hard-to-classify variant cannot hide behind the easy ones.
  - **ROC AUC (macro one-vs-rest)** — captures whether the model's probabilistic ranking is right even when the argmax is wrong. *Invariant to monotonic transforms of the predicted probabilities, so it is calibration-independent — this property is load-bearing in 2.2.12.*
  - **Log loss** — penalizes ranking *and* calibration; the metric that most directly reflects whether a calibrated confidence display will be honest at deployment.
- These metrics split into two groups by what they evaluate: accuracy / kappa / macro-F1 only use the predicted class (**argmax-based**); ROC AUC and log loss use the full predicted distribution (**probability-based**). They tend to agree on easy problems and disagree exactly at the boundary that matters for deployment. 2.2.12 makes this disagreement concrete on the four model families (SVM-RBF: competitive on every argmax metric, broken on both probability metrics).
- **Per-dAb framing.** Models and FS are always run per-dAb, never pooled. This is the right structural match for the device but it also means we cannot rescue a weak dAb by borrowing strength — every dAb stands or falls on its own ~96 observations.

*Figures/tables needed:* CV-design schematic (outer/inner/innermost); table listing the metric panel and what each one penalizes.

---

## 2.2.6. Feature selection results

> Largest section. Suggested sub-structure: (a) why we use a panel of FS methods, (b) per-method results with one or two representative dAb plots each, (c) cross-method comparison and the engineering sweet spot.

### 2.2.6.a Why a panel of FS methods
- No single FS method is consistently right at this n. LASSO favors sparsity but is hostile to correlated cAb features. RF impurity captures non-linear interactions but is biased toward high-cardinality features. MARS picks up monotone breakpoints. PLS-VIP rewards features aligned with the projection that maximizes class-explained variance. Boruta gates by shadow-feature significance. mRMR balances relevance against redundancy.
- Each is included for a specific reason; the composite score in 2.2.7 leverages disagreement between them.

### 2.2.6.b Per-method results (representative dAbs)

| Method     | What it ranks                                    | Strong dAbs (≥80% acc)           | Notable weak dAbs       |
|------------|--------------------------------------------------|----------------------------------|-------------------------|
| LASSO      | non-zero coefficients in 1-SE λ model            | dAb13, dAb17, **dAb18**, dAb16, **dAb19** | dAb4 (38%), dAb11 (53%) |
| RF (Gini)  | impurity-decrease importance                     | dAb5, dAb25, dAb13               | dAb4 (48%), dAb11 (52%) |
| MARS       | knot count + GCV-weighted importance             | *(fill in from `all_mars_runs_multiseed.rds`)* | *(fill in)*              |
| PLS-VIP    | projection-loading magnitude                     | *(fill in from `all_pls_runs_multiseed.rds`)*  | *(fill in)*              |
| Boruta     | confirmed > shadow                               | *(fill in from `all_boruta_runs_multiseed.rds`)* | *(fill in)*              |
| mRMR       | relevance − redundancy                           | *(fill in from `all_mrmr_runs_multiseed.rds`)*   | *(fill in)*              |

- **LASSO vs RF disagreement is informative, not noise.** RF wins dAb5 by ~12pp; LASSO wins dAb16 by ~20pp. This is the cleanest evidence that picking one FS method up front would have biased the final panel; the composite in 2.2.7 handles it explicitly.
- **dAb18 and dAb19 are in the top tier on both LASSO (82–83% acc) and the upper middle on RF (~78%) — and crucially with low across-seed SE.** Stable rank under reseeding is the property we want for a candidate going into expensive wet-lab augmentation.

### 2.2.6.c Performance vs. feature count — the engineering sweet spot
- Across LASSO and RF, mean accuracy + kappa + macro-F1 are essentially flat from 15 features down to ~6 per dAb, then degrade sharply between 5 and 3.
- **Kappa drops earlier and steeper than accuracy** — exactly the kind of metric disagreement flagged in 2.2.5. A panel chosen on the accuracy curve alone would still look "good" at 4 features; the kappa curve says it is not.
- Target panel size for deployment: **5–6 cAb pairs per dAb**.

*Figures/tables needed:*
- Aggregate `plot_perf_by_features` panel (accuracy / kappa / macro-F1 vs. feature count) for LASSO and RF.
- Per-dAb facet plot, with dAb18 and dAb19 panels emphasized.
- Stability heatmap: feature × seed selection frequency for dAb18 and dAb19.
- Cross-method top-feature overlap (Venn or upset plot) for one strong dAb and for dAb18, dAb19.

---

## 2.2.7. Composite score / Borda list / antibody subset results

- **Why a composite.** Six FS methods × 13 dAbs × 3 seeds gives a high-variance ranking surface. Borda-style aggregation across methods (and across seeds within a method) lets a feature that is "second-best everywhere" beat a feature that is "first in one method, absent in the others" — the failure mode that single-method FS most often produces at this n.
- **Construction.** For each dAb: (i) take the top-k features from each method × seed; (ii) score by inverse rank; (iii) sum across methods and seeds; (iv) report the dAb-level Borda list and the resulting candidate 5–6-feature subsets.
- **Antibody-level result.** A small number of dAbs survive consistently across the panel:
  - **Top stratum (high agreement, low SE, performance ≥ 80% acc with kappa keeping pace):** dAb13, dAb17, **dAb18**, **dAb19**.
  - **Method-dependent stratum:** dAb5 (RF only), dAb16 (LASSO only), dAb25 (RF only). Useful but more fragile.
  - **Weak stratum:** dAb4, dAb11 — drop.
- **Per-dAb subset result.** For dAb18 and dAb19 specifically, the top-5 Borda subset is stable across all three seeds, with the leading 3 pairs unchanged. Report the actual pair identities here once the table is regenerated.

*Figures/tables needed:*
- Borda-score heatmap (cAb features × dAbs).
- Bar chart of dAb-level composite scores, dAb18/19 highlighted.
- Pair-identity table for the recommended 5-feature subset per top-stratum dAb.

---

## 2.2.8. Final subset selection (short)

- We carry forward **dAb18 and dAb19** as the candidates for additional data collection.
- Rationale, tied to the criteria stated in 2.2.1:
  1. Stable top-tier performance across FS methods (LASSO top 5, RF upper-middle) and across seeds (small SE).
  2. Feature-count curve is flat from 15 → 6, supporting a 5-pair deployable panel.
  3. Practical wet-lab constraints (cost, availability, reagent compatibility) — fill in the project-specific reasons here.
- We acknowledge that dAb13 and dAb17 are also strong candidates on the screen; they were *not* down-selected for expensive augmentation for reasons covered in the discussion (resource allocation; the chosen pair gives broader epitope coverage; etc. — fill in).

*Figures/tables needed:* a small summary table — dAb · screen accuracy · kappa · panel size · selected for augmentation Y/N.

---

## 2.2.9. Data augmentation of final selected dAbs

- **Scope.** Augmentation runs were performed against dAb18 and dAb19 specifically. The augmented dataset extends the per-class sample count in the concentration regime where the deployed device will actually operate.
- **Concentration threshold.** Conc > 0.9 ng/mL applied **pre-split** so neither the CV folds nor the held-out test ever see sub-threshold rows. The threshold itself is a result: below 0.9 ng/mL, even the strongest dAb panel collapses to near-chance accuracy (acc 0.563 at threshold 0.0009 vs 0.825 at threshold 0.9 in single-seed RF baseline).
- **Panel used.** 10 cAb pairs (the dAb13 self-pair is dropped throughout, as in the screen panel).
- **Single-seed RF baseline on augmented data** (Phase 1): accuracy **0.825**, kappa **0.760**, ROC AUC **0.984** on the held-out test split. This is the in-silico "ceiling" we then ask the model-family comparison (2.2.12) to defend.

*Figures/tables needed:*
- Sample-count comparison: original screen vs. augmented per-class.
- RF baseline confusion matrix on the augmented test split.
- Accuracy / kappa vs. concentration threshold plot — supports the 0.9 ng/mL choice.

---

## 2.2.10. Alignment with original data

- **Question.** Can a model trained on the original screen alone be carried forward to the new (augmented) batches, or does it need updating?
- **Approach.** Bayesian sequential update with rstanarm — one-vs-rest, 4 binary fits per scenario. Original-screen posterior is used as a prior; the augmented batches (batch 1, batch 2, both) update it.
- **Result.** Sequential update consistently lifts performance over either single-batch fit alone. **Original-prior contribution decays as field data accumulates** — Original+Both ≈ Weak+Both ≈ 0.81 accuracy on the combined holdout, which is the correct Bayesian behavior and operationally means the prior matters most early in deployment.
- **Implication for the throughline.** Even when the in-silico model "transfers" technically, the prior carries diminishing weight; the *deployment workflow*, not the prior, is what makes the long-run model accurate.

*Figures/tables needed:*
- Accuracy / kappa table by scenario (orig only, batch-1 only, batch-2 only, sequential A, sequential B, sequential C, orig + batch-1, orig + batch-2, orig + both).
- Posterior-coefficient shrinkage plot for one representative one-vs-rest binary head.

---

## 2.2.11. Batch effects

- **Headline finding, worth stating bluntly.** A model trained on the original screen (cleaned CSV, dAb13 panel) and evaluated directly on the augmented assay performs **at chance** — accuracy 0.186, kappa −0.078, ROC AUC 0.494.
- **Why.** The original screen used different protocol / instrumentation than the final augmented assay format; this is expected divergence, not a bug. Confirmed by the protocol audit (see Methods §X).
- **Why this matters for the broader argument.** This is the single strongest piece of evidence in the dissertation for the metrics-vs-deployment throughline: a model can pass every in-silico check and still be *zero* useful on a new instrument generation. The defense against this is not a better algorithm — it is the in-field feedback loop (2.2.10) and OOD detection (2.2.13).
- **Predictor normalization.** Between-batch normalization (anchored on one source's mean/sd, applied to all) is required for the prior-as-prior interpretation in the sequential Bayesian work and is also what a real device would do — at inference you have one sample and cannot rescale by its own batch statistics.

*Figures/tables needed:*
- Confusion matrix: original-screen model → augmented test data.
- Per-feature mean shift (original vs. augmented), one panel per pair, ordered by shift magnitude.
- A "before/after sequential update" performance comparison bar chart.

---

## 2.2.12. Modeling and prediction results with new data

- **Setup.** Five seeds (1969, 1986, 2024, 2025, 2026); fresh `initial_split` + `vfold_cv` per seed; combined augmented data (both batches, threshold ≥ 0.9 ng/mL); identical 10-cAb panel and identical recipe across all four model families (no model is given a feature-count advantage).
- **Multi-seed results** (mean ± SE across 5 seeds, test set, n ≈ 57). Two metric groups, separated to make the throughline visible:

  |          | *Argmax-based (class prediction)* |             |              | *Probability-based (ranking + calibration)* |                     |
  |----------|-----------------------------------|-------------|--------------|---------------------------------------------|---------------------|
  | Model    | Accuracy        | Kappa         | Macro-F1     | ROC AUC                  | Log loss            |
  | **RF**   | 0.864 ± 0.011   | 0.812 ± 0.014 | 0.866 ± 0.017 | **0.975 ± 0.003**       | 0.471 ± 0.021       |
  | SVM-RBF  | 0.857 ± 0.015   | 0.802 ± 0.019 | 0.856 ± 0.017 | 0.910 ± 0.015           | 0.783 ± 0.055       |
  | ElasNet  | 0.853 ± 0.012   | 0.798 ± 0.016 | 0.848 ± 0.016 | 0.967 ± 0.008           | 0.473 ± 0.030       |
  | XGBoost  | 0.843 ± 0.017   | 0.784 ± 0.022 | 0.840 ± 0.019 | 0.969 ± 0.004           | **0.421 ± 0.028**   |

- **Single-seed result was a misleading three-way tie at ~0.86 acc.** Multi-seed splits the rank cleanly and consistently across seeds: RF > SVM-RBF > ElasNet > XGBoost on the three argmax-based metrics, with RF best on AUC and XGBoost best on log loss.

- **The metrics throughline payoff (the punch line of the chapter).** On all three **argmax-based** metrics — accuracy, kappa, macro-F1 — RF and SVM-RBF are within ~1pp of each other and within SE overlap. Read on those metrics alone, the two models are interchangeable. On the two **probability-based** metrics they separate decisively: ROC AUC by 6.5pp and log loss by 66%, both well outside SE.

  This is the concrete instance of the in-silico-vs-experimental-metric gap the chapter has been pointing at: *the metric you choose to evaluate on decides whether the model is "deployable."*

- **Why this isn't just a calibration artifact of SVM's Platt scaling** (anticipated reader objection). The two probability-based metrics differ in what they penalize:
  - **ROC AUC** is invariant to any monotonic transformation of the predicted probabilities. It depends only on the *ranking* of probabilities — not their magnitudes, not the decision-boundary width, not how Platt's sigmoid stretched them. If the SVM's only problem were that its probabilities were "compressed" or "wider," AUC would be unaffected.
  - **Log loss** penalizes ranking *and* magnitude calibration.
  - Both are worse for SVM-RBF. That means the underlying probabilistic *ranking* is genuinely less informative on this problem — not merely miscalibrated.
  - Mechanistically, SVM-RBF's class probabilities come from Platt scaling fitted on an inner CV split. At this training-set size (~225 rows) the sigmoid fit is unstable, especially when decision values cluster. RF's probabilities come from tree-vote proportions — a frequency estimate that is well-behaved at small n. This is a documented pattern, not specific to this assay.

- **Why log loss is the right deployment-time metric here, not just a training-time objective.** The device displays a variant label *and* a confidence number to the operator. The downstream action — accept, retest, request orthogonal confirmation — is conditional on that confidence. Log loss is precisely the metric that asks whether the displayed confidence is honest: a model with low log loss says "85%" when it is right ~85% of the time, while a model with high log loss says "85%" when it is right 60% of the time. This holds independently of what loss the model was trained on.

- **Best hyperparameters were stable across seeds** (RF: mtry=3, min_n ∈ {2,5}; ElasNet: penalty ≈ 5e-3, **mixture ≈ 0.25 — ridge-leaning, not pure LASSO**; SVM-RBF: cost=4, σ=0.05; XGBoost: depth=2, lr=0.05, mtry=5, min_n=5). Two side-observations:
  1. ElasNet *could* zero out features but chose not to: the winning mixture is ridge-leaning, meaning all 10 features are retained with shrinkage. This is concordant with the 2.2.6 conclusion — the FS pipeline has already trimmed the panel to the engineering sweet spot; there is no remaining redundancy for ElasNet to discover.
  2. XGBoost's log-loss win despite the lowest accuracy is the standard boosted-trees pattern; it is the model to reach for if probability calibration is the *primary* deciding criterion.

- **Why we recommend RF for deployment.** Best on argmax metrics, best AUC by a clear margin, second-best log loss (within SE of ElasNet's), and the smallest SE on AUC across seeds. RF additionally provides **tractable per-feature permutation/impurity importances**, which feed back into the panel-design loop (which pair, on which dAb, is doing the work?) in a way that SVM-RBF's RBF-kernelized decision function does not. For a device that may need ongoing panel revision as new variants emerge, that interpretability is itself a deployment advantage, not a nicety.

> **TODO — panel-size stress test.** The numbers above use the full 10-cAb panel. The device's design constraint targets a 5–6-pair panel (see 2.2.1, 2.2.6.c). Re-run the four-family multi-seed comparison while artificially constraining the antibody count — e.g. sweep n_pairs ∈ {3, 4, 5, 6, 8, 10}, selecting each subset from the composite-score rankings in 2.2.7 — to characterise *where* each model's performance degrades and *which metric degrades first*. The expected reportable result is a per-model degradation curve (one per metric) that tells us (i) whether RF holds its AUC/log-loss advantage at the design limit, (ii) whether one of the simpler models catches up at small panels because the kernel/tree complexity no longer pays off, and (iii) what the lower bound on n_pairs is before the probability-based metrics collapse — i.e. before the device can no longer honestly display confidence. This makes the in-silico → device handoff concrete in the section that immediately precedes the deployment discussion.

*Figures/tables needed:*
- The table above as a clean publication-ready table.
- Per-seed scatter (one dot per seed) for accuracy, kappa, AUC, log loss across the four models — visually demonstrates the SVM-RBF metric divergence.
- Reliability diagram / calibration curve for RF vs. SVM-RBF — the single most persuasive figure for the calibration argument.
- Per-class one-vs-rest ROC curves for the four families, dAb18 and dAb19 separately.

---

## 2.2.13. Deployment method and considerations

- **Concentration threshold as a deployment gate.** Reject any inference below the 0.9 ng/mL floor with an "out-of-range" rather than a confident class call; the 2.2.9 result is the operational justification.
- **Uncertainty-aware output.** The device should display calibrated class probabilities, not just the argmax. The choice of RF (or XGBoost, if probability calibration is paramount) over SVM-RBF is decided by 2.2.12.
- **OOD / batch-shift detection.** Motivated directly by 2.2.11: a new instrument generation can drive a screen-only model to chance. Two complementary safeguards:
  1. **In-field feedback loop** — periodic Bayesian update from collected field samples; the prior decays as field data accumulates (2.2.10), which is correct behavior, but it requires that field samples are actually collected and labeled.
  2. **Run-time OOD score** — a per-sample distance/likelihood check against the training distribution (e.g., Mahalanobis on the normalized 10-pair vector, or a one-class score) that flags samples falling outside the manifold the model was trained on.
- **Algorithm compatibility with new incoming data.** Between-batch normalization is anchored on training-time statistics; new samples are normalized to those statistics, not to their own batch. This matches real device operation (one sample at a time) and is consistent with the Bayesian-update setup.
- **What "deployment" looks like in practice (summary).**
  1. Run sample on device → 10-pair fluorescence vector (dAb18 or dAb19 panel).
  2. Normalize against training stats.
  3. Threshold check: ≥ 0.9 ng/mL inferred from total signal? If not → out-of-range.
  4. OOD score within tolerance? If not → flag for re-test / log for retraining.
  5. RF predicts class probabilities → display variant + calibrated confidence.
  6. Periodically: collected, confirmed field samples are appended to the next Bayesian update of the model.

*Figures/tables needed:*
- Deployment-flow diagram (the 6-step list above as a figure).
- Reliability diagram for the deployed RF (carries over from 2.2.12).
- An illustrative OOD-score histogram: in-distribution training samples vs. a held-out batch from 2.2.11.

> **TODO — leave-one-variant-out OOD demonstration.** We can give the OOD safeguard a concrete in-data demonstration without needing a new instrument-shift dataset, by treating one variant at a time as an unseen OOD class. Procedure: for each variant V in {WT, Beta, Delta, Omicron}, (i) remove all V rows from training; (ii) train a 3-class classifier on the remaining variants plus an OOD score / 4th "unknown" head (Mahalanobis distance, one-class SVM, or an evidential / Dirichlet-prior head — pick one and document); (iii) at test time, score the held-out three known variants *and* the held-out V samples; (iv) report (a) classification accuracy on the three knowns, (b) OOD-detection AUC on V-vs-known. Cycling V through all four variants shows the method works regardless of which variant is "novel" and produces a 4-panel figure that directly visualises the deployment-time scenario the chapter has been arguing for. This is the cleanest possible illustration of the OOD-safeguard claim in 2.2.13 and the most direct response to "what happens when a new variant shows up in the field."

---

### Cross-cutting figures/tables (suggest pulling into one place near the methods/results boundary)

- **Master FS-method table:** method · what it ranks · what metric it optimizes · why it's in the panel · which composite weight it gets. Resolves the open note at the bottom of the original outline.
- **Glossary of metrics:** accuracy, kappa, macro-F1, ROC AUC, log loss — one sentence each on what they reward and what they punish. Anchors the throughline.

### NOTE (carry-over from original outline): include the master FS-method table — handled above under "Cross-cutting figures/tables."
