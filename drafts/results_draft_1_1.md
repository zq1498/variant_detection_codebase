# 2.2. Results

> **Integration notes (delete before final).**
> - This file merges `results_draft.md` (skeleton + technical content, bullet voice) and `results_draft_1.md` (Zach's prose expansion). Prose voice from `results_draft_1.md` is the spine; technical content from the skeleton is folded in.
> - Substantive points to review:
>   - **2.2.1 reason-for-metrics paragraph rewritten** — original "large feature space + small n" rationale replaced with the three sharper reasons we agreed on (sampling variance from small folds; 4-class chance baseline; concentration-dependent error structure). Your prose voice preserved.
>   - **2.2.11 explanation kept as yours, not mine** — the smaller-array-→-less-cross-reactivity story is your domain claim; I added the specific chance-level numbers to attach evidence but did not override the mechanism.
>   - **2.2.5 "three-level CV" softened to two-level** — the code (LASSO/RF tidymodels) implements an outer split + inner k-fold; a true innermost ranking-stability layer is not implemented. Verify before finalising.
>   - **2.2.13 cut-off sentence** completed cautiously and tied to the leave-one-variant-out OOD TODO.
> - TODO call-outs (panel-size stress test, leave-one-variant-out OOD) preserved as blockquotes.
> - All `[VERIFY: ...]` markers are open questions where I am not confident; please check before submission.

---

## 2.2.1 Problem statement and design constraints

Our goal, ideally, is to identify a minimal antibody microarray panel that can reliably classify SARS-CoV-2 S-protein variants (WT, Beta, Delta, and Omicron). The design of the D4 assay itself imposes a number of practical design constraints (FIG X). First, only one dAb can be feasibly measured at once. As such, models must be evaluated per-dAb, and data cannot be pooled, nor antibody pairs drawn, across different dAbs. Second, while spatial multiplexing allows for straightforward quantification of many array spots, the physical space within the field-deployable cassette, the acquisition area of a POC-compatible fluorescent detector, and the need for replicate cAb spots decreases the number of unique cAbs that can be printed. As such, the optimal assay would be limited to between 3–10 unique cAbs. Moreover, since variant prediction must occur on a single sample at a time (i.e., the current patient sample), we cannot employ any form of batch-level rescaling at deployment. This limits our preprocessing pipeline as well as our ability to calibrate predictions (2.2.12).

Here, we consider a 'successful' antibody pair subset as one that (a) classifies above chance with margin on held-out CV folds, (b) provides a meaningful improvement in Cohen's κ or macro-F1, not just accuracy, over a smaller subset, and (c) yields stable results across reshuffles of the data (multi-seed).

While the D4 allows us to generate relatively large amounts of screening data across antibody pairs, the multi-class, multi-concentration nature of the problem creates a large theoretical sample space that we can only sparsely populate with training examples. Accounting for training samples that sit below the practical limit of quantitation, the overall n per dAb is comparatively small — roughly 96 observations, which translates to ~20 rows per held-out CV fold. In this regime, accuracy alone is a poor primary readout for three reasons that bite hardest *on this specific problem*:

1. **Small test folds amplify sampling variance into accuracy.** At ~20 rows per fold, a single flipped prediction is a ~5pp swing in accuracy; even on the larger augmented-data test split (n ≈ 57), it is ~1.75pp. A 2–3pp accuracy gap between two models or two seeds can easily reflect which samples happened to fall in which split rather than a true model difference. Multi-seed reporting (below) mitigates this, but it does not eliminate the fact that a single accuracy number is noisy.
2. **The chance baseline for a four-class problem is 25%, not 50%.** Accuracy values are accordingly inflated relative to binary-classification reporting; an 80% accuracy reads differently here than in most published benchmarks. Cohen's κ is chance-corrected and so makes the "is the model doing work above chance" question structurally explicit at any n.
3. **Variants are roughly balanced overall, but the error structure is concentration-dependent.** A model can score respectably on accuracy by getting the high-concentration samples right and quietly failing in the low-concentration regime where the deployed device will actually operate. Macro-F1 (per-class average) and per-concentration accuracy panels surface this; raw accuracy hides it.

To maximize the likelihood of successful experimental validation, we therefore comprehensively report accuracy, Cohen's κ, macro-F1, AUC-ROC, and log-loss, finding that the combination of these metrics yields more realistic evaluation of the models and can directly inform the engineering design process. The justification for *each* metric, and how they split into two interpretive groups, is laid out in 2.2.5.

*Figures/tables needed:* assay schematic (existing); device-constraint diagram showing 1 dAb × N cAbs per panel.

---

## 2.2.2. Data collection and description

In total, 13 detection antibodies were each screened against 29 unique capture antibodies for a total of over 400 unique dAb–cAb pair readouts. Each dAb is analyzed as an independent dataset of cAb features. For each of these, 8 concentrations of spike protein (ranging from 0 ng/mL to 500 ng/mL) were run in triplicate for each of the four variants, yielding a total of 96 observations per dAb. Evaluated together, numerous sources of variation must be accounted for, including the variant label (the target signal), the concentration (a covariate that mostly modulates signal magnitude rather than class identity above the LoD), replicate (within-run technical variation), and — as we will see later during data augmentation — batch drift (2.2.11).

The primary cleaned screen dataset is in `data/covid_screen_data_cleaned.csv`. The augmented batches (collected against dAb18 and dAb19 specifically) live alongside, in [VERIFY: confirm augmented-data path before finalizing].

*Figures/tables needed:* table listing dAbs with the number of cAbs each was screened against, the concentration grid, and a sample-count breakdown by variant × concentration.

---

## 2.2.3. Screen results and exploratory data analysis

First, we analyzed the dynamic range of each cAb:dAb pair, noting that fluorescence intensity increases monotonically with antigen concentration above the LoD but flattens or becomes dominated by noise below it. As a heuristic pre-filter, we discarded cAb:dAb pairs whose change in fluorescence intensity was less than [VERIFY: insert X] orders of magnitude across the entire concentration range. This removed pairs with poor dose-response behavior such as extremely high background or chronically low signal intensity.

Next, we investigated the relationship between class separability and choice of dAb. PCA on the standardized per-dAb feature blocks, together with unsupervised inspection of variant-stratified fluorescence distributions, indicates that some dAbs cleanly separate at least one variant across most of their cAb partners (dAb13, dAb17 are representative), while others show essentially overlapping fluorescence distributions across variants at every concentration (dAb4, dAb11 are representative). When decomposed to their principal components, the first two PCs of the standardized per-dAb feature blocks show partial variant clustering for 'strong' dAb candidates such as dAb13, whereas 'weak' dAbs like dAb4 show much less clear partitioning, generally with minimal obvious structure.

Taken together, this analysis indicates two things that motivate the supervised feature-selection work that follows: (1) certain dAbs are inherently more amenable to variant discrimination, despite each being interrogated by the same cAb array, and (2) promising dAb feature subsets contain significant redundant information among their cAbs, suggesting that smaller deployable panels are achievable without large losses in discriminative power. Encouraged, we shifted our focus to supervised learning methods to identify which exact pair subsets would be optimal for further characterization and validation.

*Figures/tables needed:*
- Concentration-response curves for one strong dAb (dAb13 or dAb17), one weak (dAb4 or dAb11), and **dAb18 + dAb19** (the latter two highlighted in the caption — these are the candidates we will carry into the augmented-data experiments).
- PCA score plots for the same set of dAbs.
- Boxplot or heatmap of variant-wise fluorescence dispersion across all dAbs.

---

## 2.2.4. Data set preprocessing

Briefly, per-dAb feature subsets were Box-Cox transformed to obtain approximately normal distributions, then centered and scaled for methods requiring normalization (LASSO, MARS, PLS). Zero-variance features were then removed. Concentration was retained as a label for performance visualization but was not included as a predictor, as this information would be unavailable during a clinical use case.

To prevent data leakage, all preprocessing parameters are estimated *inside* each cross-validation fold rather than on the full dataset before splitting. Concretely: each CV split divides the data into a *training portion* (used to fit) and a *held-out portion* (used to evaluate, and to score feature ranks). Box-Cox λ, feature means, feature standard deviations, and any kNN-imputation neighbours are estimated only on the training portion of each split and then applied to the held-out portion. They are *not* estimated once on the full dataset and reused. The standard accidental shortcut — fitting Box-Cox or normalize globally and then splitting — would let the held-out rows influence their own normalization, which inflates evaluation metrics and biases feature ranks toward features whose scale happened to align with the held-out set. In practice, this is enforced by `tidymodels` `recipes` being baked inside the resample machinery, so the leakage is excluded by construction rather than by discipline.

A small substring-matching bug in the original preprocessor (e.g. `dAb1` matching `dAb13`) was identified and fixed via anchored matching so that the pair list propagated to feature selection is correct.

Detailed protocol in Methods §X.

*Figures/tables needed:* none in the main text; a preprocessing diagram could go in methods.

---

## 2.2.5. Cross-validation method and safeguards against overfitting

In an effort to minimize the impact of overfitting our models to the small sample size during the extensive feature selection process, we instituted a nested cross-validation strategy. An outer split was held out for final evaluation of downselected feature sets, and an inner repeated k-fold cross-validation was employed for hyperparameter and feature-count tuning. [VERIFY: the original draft mentioned an innermost layer for ranking-stability; the implemented code uses two CV levels rather than three. Either describe the implemented two-level design, or, if a third level was intended and is to be added before submission, restore the three-level language here.] By partitioning the data this way, we can partly safeguard against over-optimistic feature-selection results stemming from reuse of the same data for both feature selection and evaluation.

Moreover, each feature selection and evaluation run is repeated three times (seeds 1969, 1986, 2024) for the screen-level feature selection and five times (additionally including 2025, 2026) for the augmented-data modeling in 2.2.12. We report the mean ± SE across these seeds, with large per-seed SEs for a given dAb serving as a soft indicator of potentially unreliable or stochastic feature selection.

For every CV fold we record the following metric panel:

***[insert metric set table here — one row per metric: name · what it rewards · what it penalizes · why it is in the panel.]***

- **Accuracy** — intuitive, but noisy at small fold sizes and uninformative about class-wise behavior or probabilistic quality.
- **Cohen's κ** — chance-corrected; the metric that punishes "guess the majority class" models structurally.
- **Macro-F1** — averages F1 across the four variants so a single hard-to-classify variant cannot hide behind the easy ones.
- **AUC-ROC (macro one-vs-rest)** — captures whether the model's probabilistic *ranking* is correct, even when the argmax is wrong. AUC-ROC is invariant to any monotonic transformation of the predicted probabilities, which makes it calibration-independent — a property that is load-bearing in 2.2.12, where it lets us distinguish "the model is miscalibrated" from "the model's probability ranking is genuinely worse."
- **Log loss** — penalizes ranking *and* magnitude calibration; the metric that most directly reflects whether a calibrated confidence display will be honest at deployment.

These metrics split into two groups based on what they are specifically evaluating: either the argmax of the variant probability vector, or the full predicted probability distribution from the model. Accuracy, κ, and macro-F1 can be considered **argmax-based** — they only consider the final predicted class during their calculation. By comparison, AUC-ROC and log loss are **probability-based** — they utilize the full predicted distribution. While often in agreement when evaluating 'easier' problems, such as binary classification or prediction with large quantities of training data, the two groups can notably disagree in equally many cases, especially in multi-class settings and at deployment. 2.2.12 demonstrates this disagreement concretely on the four model families considered: SVM-RBF is competitive on every argmax metric but breaks on both probability metrics.

Finally, all models and FS runs are conducted per-dAb, never pooled. This is the right structural match for the device — one dAb is interrogated at a time — but it also means we cannot rescue a weak dAb by borrowing strength across dAbs. Every dAb stands or falls on its own ~96 observations.

*Figures/tables needed:* CV-design schematic (outer / inner [/ innermost, if added]); the metric-panel table referenced above.

---

## 2.2.6. Feature selection results

### 2.2.6.a Why a panel of FS methods, and what each one tells us

No single feature-selection method is consistently right at this sample size, and each carries a different inductive bias that surfaces a different aspect of the cAb feature space. LASSO favours sparsity but is hostile to correlated cAb features, and so tends to pick one representative from a correlated cluster and drop the rest. Random Forest impurity / permutation importance captures non-linear interactions and combinations but is mildly biased toward high-cardinality features. MARS picks up monotone breakpoints and so highlights features with clear dose-response inflections. PLS-VIP rewards features aligned with the projection that maximises class-explained variance. Boruta gates by comparison against shadow-feature significance, so its output is more conservative and reflects "is this feature reliably better than a permuted version of itself." mRMR balances relevance against redundancy, so it tends to pick complementary sets rather than the *k* highest-scoring individuals.

Each is included because it surfaces a different kind of useful feature, and the composite score in 2.2.7 leverages their *disagreement* — features that survive across multiple methods are more likely to carry real signal than features that win on a single criterion.

### 2.2.6.b What the feature-selection results show

Across the panel of methods, three patterns emerge: dAbs that consistently rank near the top regardless of method (dAb13, dAb17, **dAb18**, **dAb19** in LASSO's top tier at ≥ 80% accuracy; dAb13 also strong in RF, with dAb18 and dAb19 in the upper-middle tier at ~78% RF accuracy); dAbs whose strength is method-dependent (dAb5 and dAb25 are strong under RF but not LASSO; dAb16 is strong under LASSO but not RF); and dAbs that are consistently weak across methods (dAb4 and dAb11, sitting at 38–53% accuracy).

| Method     | What it ranks                                    | Strong dAbs (≥ 80% acc)                          | Notable weak dAbs       |
|------------|--------------------------------------------------|--------------------------------------------------|-------------------------|
| LASSO      | non-zero coefficients in 1-SE λ model            | dAb13, dAb17, **dAb18**, dAb16, **dAb19**        | dAb4 (38%), dAb11 (53%) |
| RF (Gini / permutation) | impurity-decrease / permutation importance | dAb5, dAb25, dAb13                            | dAb4 (48%), dAb11 (52%) |
| MARS       | knot count + GCV-weighted importance             | [VERIFY: fill from `all_mars_runs_multiseed.rds`] | [VERIFY]                |
| PLS-VIP    | projection-loading magnitude                     | [VERIFY: fill from `all_pls_runs_multiseed.rds`]  | [VERIFY]                |
| Boruta     | confirmed > shadow                               | [VERIFY: fill from `all_boruta_runs_multiseed.rds`] | [VERIFY]              |
| mRMR       | relevance − redundancy                           | [VERIFY: fill from `all_mrmr_runs_multiseed.rds`]   | [VERIFY]              |

The disagreement between methods is itself informative rather than noise: RF wins dAb5 by roughly 12pp over LASSO, while LASSO wins dAb16 by roughly 20pp over RF. This is the clearest empirical justification for the composite score in 2.2.7 — choosing any single FS method up front would have biased the final panel toward the kinds of structures that method preferentially detects. Importantly for the narrative that follows, **dAb18 and dAb19 are in LASSO's top tier and in RF's upper-middle tier, with consistently low across-seed SE.** Stable rank under reseeding is exactly the property one wants for a candidate that is going into expensive wet-lab augmentation.

### 2.2.6.c Performance vs. feature count — balancing performance and engineering considerations

Across both LASSO and RF, mean accuracy, κ, and macro-F1 are essentially flat from 15 features down to about 6 per dAb, after which performance degrades sharply between 5 and 3 features. Notably, **κ drops earlier and more steeply than accuracy** — the same disagreement between argmax-based metrics that we previewed in 2.2.5. A deployable panel chosen on the accuracy curve alone would still look respectable at 4 features; the κ curve says that it is not. Taken together, the curves point to a target deployment panel size of **5–6 cAb pairs per dAb**, which is also comfortably within the 3–10-pair device constraint from 2.2.1.

*Figures/tables needed:*
- Aggregate `plot_perf_by_features` panel (accuracy / κ / macro-F1 vs. feature count) for LASSO and for RF.
- Per-dAb facet plot of the same, with the dAb18 and dAb19 panels emphasised.
- Stability heatmap: feature × seed selection frequency for dAb18 and dAb19.
- Cross-method top-feature overlap (Venn or upset plot) for one representative strong dAb and for dAb18 and dAb19.

---

## 2.2.7. Composite score / Borda list / antibody subset results

Across six distinct feature-selection methods and 13 dAbs, antibody-pair subset rankings can vary widely across algorithms, or even between seeds within a single method. By first normalizing, then aggregating the key rankings from each method, we can identify features that are broadly determined to be important. This helps to avoid cases wherein certain pair subsets are selected based on a single very strong result in one method despite poor results in other methods.

To construct our composite score for each dAb: (1) extract the top-k features from each method × seed; (2) score by inverse importance ranking; (3) sum across all methods and seeds; (4) report the final ranked list with candidate feature subsets.

Three strata of pairs emerge (Fig …). A **top stratum** contains pairs that perform well across many feature-selection methods (dAb13, dAb17, dAb18, dAb19); a **middle stratum** consists of pairs that are of high importance in a single, or very few, methods (dAb5 under RF; dAb16 under LASSO; dAb25 under RF) — these remain useful candidates but are more fragile in the sense of being method-dependent; and a **lower stratum** of pairs that are consistently selected against across all methods (dAb4, dAb11), which we drop from further consideration. For dAb18 and dAb19 specifically, strong importance ranking can be observed in all methods and across all seeds, and the top-5 Borda subset is stable across all three seeds with the leading three pairs unchanged. [VERIFY: regenerate the per-dAb Borda subset table and insert the actual pair identities for dAb18 and dAb19 here.]

*Figures/tables needed:*
- Borda-score heatmap (cAb features × dAbs).
- Bar chart of dAb-level composite scores, with dAb18 and dAb19 highlighted.
- Pair-identity table for the recommended 5-feature subset per top-stratum dAb.

---

## 2.2.8. Final subset selection for experimental validation

Having demonstrated strong performance across all feature-selection methods, as well as other favorable characteristics from a device fabrication perspective — high post-printing stability, consistent fluorescent labeling, low LoDs, and minimal cross-reactivity in complex media — dAb18 and dAb19 were selected for follow-up data augmentation. We acknowledge that dAb13 and dAb17 are also strong candidates on the basis of the screen alone; they were not down-selected for expensive augmentation for the engineering and reagent-related reasons noted above and discussed further in §[discussion].

To collect the augmented data, we again employed the D4 assay. Data was collected in two batches, utilizing different lots of [VERIFY: insert reagent — slide/polymer/cAb stock] during fabrication. Assays were performed after being desiccated for one month under vacuum at room temperature, and the batches themselves were produced approximately two months apart. This staggered fabrication and storage protocol was an intentional stress on the assay reproducibility, since any candidate for deployment must remain performant across reagent lots and shelf time.

*Figures/tables needed:* a small summary table — dAb · screen accuracy · κ · panel size · selected for augmentation Y/N · reason.

---

## 2.2.9. Data augmentation of final dAb subsets

Having identified and characterized the performance of promising cAb/dAb pair subsets, we next sought to obtain additional experimental data to better fill out our predefined feature space, focused specifically on the dAb18 and dAb19 panels and on the concentration regime in which the deployed device will actually operate.

A concentration threshold of > 0.9 ng/mL was applied **pre-split**, so that neither the CV folds nor the held-out test set ever see sub-threshold rows. The threshold itself is a result rather than an arbitrary cut: below 0.9 ng/mL, even the strongest dAb panel collapses to near-chance accuracy. A single-seed RF baseline on the augmented data at a threshold of 0.0009 ng/mL achieves only 0.563 accuracy; at the 0.9 ng/mL threshold the same pipeline achieves 0.825 accuracy, 0.760 κ, and 0.984 AUC-ROC on the held-out test split. The augmented panel itself consists of 10 cAb pairs (the dAb13 self-pair is dropped throughout, as in the screen).

These single-seed numbers establish the in-silico "ceiling" that the multi-seed, multi-model comparison in 2.2.12 then has to defend.

*Figures/tables needed:*
- Sample-count comparison: original screen vs. augmented, per class.
- RF baseline confusion matrix on the augmented test split.
- Accuracy / κ vs. concentration threshold — supports the 0.9 ng/mL choice.

---

## 2.2.10. Integration of new data with the original screen

In a production context, it is natural to question whether injecting additional data taken from entirely distinct production lots could adversely affect the robustness of the model. Two related questions follow: can a model trained on the original screen alone be carried forward to the new augmented batches and updated, or rather should it be replaced entirely by the more targeted data-augmentation batches?

To answer this, we performed a Bayesian sequential analysis using `rstanarm` — one-vs-rest, four binary fits per scenario — in which the original-screen posterior was used as a prior and the augmented batches (batch 1, batch 2, or both) provided the likelihood updates. The result is straightforward and concordant with the batch-effects analysis in 2.2.11 below: sequential updates that include augmented batches consistently lift performance over either single-batch fit alone, and the contribution of the original-screen prior decays as field data accumulates. Original + Both ≈ Weak prior + Both ≈ 0.81 accuracy on the combined holdout — which is correct Bayesian behaviour and operationally means the original-screen prior matters most early in deployment, before field data has had time to accumulate. The deployment workflow, not the prior, is what makes the long-run model accurate.

*Figures/tables needed:*
- Accuracy / κ table by scenario (orig only, batch-1 only, batch-2 only, sequential through each combination, orig + batch-1, orig + batch-2, orig + both).
- Posterior-coefficient shrinkage plot for one representative one-vs-rest binary head.

---

## 2.2.11. Batch effects

When comparing across the three sources of experimental data — original screen data, augmentation batch 1, and augmentation batch 2 — it is clear that incorporation of the original screen data, either as a Bayesian prior (2.2.10) or as an additional set of batch-adjusted samples to combine with the augmentation batches during training, does not yield improved predictive performance over what the augmented batches alone provide. A model trained on the original screen and evaluated directly on the augmented assay performs essentially at chance (accuracy 0.186, κ −0.078, AUC-ROC 0.494).

The most plausible reason for this is the difference in assay setup between the screen and the augmentation campaign — specifically, the drastically decreased size of the cAb array in the augmentation batches relative to the original screen. More highly multiplexed arrays naturally produce more sources of variation, principally in the form of cross-reactivity between cAbs. The original screen's apparent predictive power is therefore likely fed in part by this cross-reactivity signal, which is no longer available in the smaller augmented array. Despite producing seemingly strong predictive models for variant discrimination in isolation, the screen data cannot be meaningfully included as additional training material for our new, more targeted device — its predictive structure relies on features that the device cannot reproduce.

By contrast, the two augmentation batches *do* combine usefully despite being fabricated from different lots of polymerized slides, printed in different runs, and using different aliquots of the same cAbs and dAbs: including the second batch with the first produces significantly better models than either batch alone. This is the load-bearing positive result of this section — it tells us that the relevant batch effects across reagent lots and runs are within the regime that a single model can tolerate, while the screen-vs-augmentation gap is not.

*Figures/tables needed:*
- Confusion matrix: original-screen model → augmented test data (visualises the chance-level result).
- Per-feature mean shift (original vs. augmented), one panel per pair, ordered by shift magnitude — gives a mechanistic view consistent with the cross-reactivity hypothesis.
- Performance comparison bar chart: augmented batch 1 alone, batch 2 alone, both combined, and (for contrast) screen-as-prior + augmentation.

---

## 2.2.12. Modeling and prediction results with new data

Having established that the augmented batches form a self-consistent training set (2.2.11), we next compared four model families — Random Forest, XGBoost, Elastic Net, and SVM-RBF — on the combined augmented data. Each seed (1969, 1986, 2024, 2025, 2026) drew a fresh `initial_split` plus `vfold_cv`, with identical recipes and an identical 10-cAb panel across all four model families so that no model is given a feature-count advantage. Hyperparameters were tuned per seed by inner-fold cross-validation.

The multi-seed test-set results (mean ± SE across the five seeds, n ≈ 57 per seed) are organized below into the two metric groups introduced in 2.2.5:

|          | *Argmax-based (class prediction)* |             |              | *Probability-based (ranking + calibration)* |                     |
|----------|-----------------------------------|-------------|--------------|---------------------------------------------|---------------------|
| Model    | Accuracy        | Cohen's κ     | Macro-F1     | AUC-ROC                  | Log loss            |
| **RF**   | 0.864 ± 0.011   | 0.812 ± 0.014 | 0.866 ± 0.017 | **0.975 ± 0.003**       | 0.471 ± 0.021       |
| SVM-RBF  | 0.857 ± 0.015   | 0.802 ± 0.019 | 0.856 ± 0.017 | 0.910 ± 0.015           | 0.783 ± 0.055       |
| ElasNet  | 0.853 ± 0.012   | 0.798 ± 0.016 | 0.848 ± 0.016 | 0.967 ± 0.008           | 0.473 ± 0.030       |
| XGBoost  | 0.843 ± 0.017   | 0.784 ± 0.022 | 0.840 ± 0.019 | 0.969 ± 0.004           | **0.421 ± 0.028**   |

A single-seed run of this comparison had produced a misleading three-way tie at ~0.86 accuracy. The multi-seed view splits the rank cleanly and consistently: RF > SVM-RBF > ElasNet > XGBoost on all three argmax-based metrics, with RF best on AUC-ROC and XGBoost best on log loss.

The result that bears most directly on the throughline of this chapter is the contrast between the two metric groups. On all three argmax-based metrics — accuracy, κ, macro-F1 — RF and SVM-RBF sit within roughly 1pp of each other and well within SE overlap. Read on those metrics alone, the two models are interchangeable. On the two probability-based metrics, however, they separate decisively: AUC-ROC differs by 6.5pp and log loss by ~66%, both well outside SE. This is the concrete instance of the in-silico-vs-experimental-metric gap that the chapter has been pointing at: *the metric that one chooses to evaluate on directly decides whether a given model is "deployable" — even when two models look identical on the most commonly reported metric.*

A reasonable concern at this point is whether SVM-RBF's higher log loss is simply an artifact of Platt-scaled probability calibration — for instance, whether the SVM's decision boundary is "wider" or its probability outputs are "compressed," giving it a higher log loss without genuinely worse classification structure. The decisive test here is the AUC-ROC. AUC-ROC is invariant to any monotonic transformation of the predicted probabilities — it depends only on the *ranking* of probabilities, not their magnitudes, the decision-boundary width, or how Platt's sigmoid stretched them. If the SVM's only problem were that its probabilities were miscalibrated, the AUC-ROC would be unaffected. Both AUC-ROC and log loss are worse for SVM-RBF, which tells us the underlying probabilistic *ranking* is genuinely less informative on this problem, not merely miscalibrated. Mechanistically, SVM-RBF's class probabilities come from Platt scaling fitted on an inner CV split, and at our training-set size (~225 rows) the sigmoid fit is unstable, especially when decision values cluster. RF's probabilities come from tree-vote proportions — a frequency estimate that is well-behaved at small *n*. This is a documented pattern, not specific to this assay.

Log loss is also the right deployment-time metric here, not merely a training-time objective, because the device displays a variant label *and* a calibrated confidence number to the operator. The downstream action — accept the call, request a retest, or send for orthogonal confirmation — is conditional on that confidence. Log loss is precisely the metric that asks whether the displayed confidence is honest: a model with low log loss says "85%" when it is right roughly 85% of the time; a model with high log loss says "85%" when it is right roughly 60% of the time. This is independent of which loss the model was trained on.

The best hyperparameters were stable across seeds (RF: `mtry = 3`, `min_n ∈ {2, 5}`; ElasNet: `penalty ≈ 5e-3`, **`mixture ≈ 0.25` — ridge-leaning, not pure LASSO**; SVM-RBF: `cost = 4`, σ = 0.05; XGBoost: depth = 2, learning rate = 0.05, `mtry = 5`, `min_n = 5`). Two observations worth flagging: ElasNet could in principle have zeroed out features but chose not to — the winning mixture is ridge-leaning, meaning all 10 features were retained with shrinkage, which is concordant with the conclusion in 2.2.6 that the FS pipeline has already trimmed the panel to the engineering sweet spot. Second, XGBoost's log-loss win despite the lowest accuracy is the standard boosted-trees pattern; it is the model to reach for if probability calibration is the *primary* deciding criterion.

On balance, RF is the deployment recommendation. It is best on the argmax metrics, best on AUC-ROC by a clear margin, second-best on log loss (within SE of ElasNet's), and has the smallest SE on AUC-ROC across seeds. RF additionally provides **tractable per-feature permutation and impurity importances**, which feed back into the panel-design loop — which pair, on which dAb, is doing the work? — in a way that SVM-RBF's RBF-kernelised decision function cannot. For a device that may need ongoing panel revision as new variants emerge, that interpretability is itself a deployment advantage rather than a nicety.

> **TODO — panel-size stress test.** The numbers above use the full 10-cAb panel. The device's design constraint targets a 5–6-pair panel (see 2.2.1 and 2.2.6.c). Re-run the four-family multi-seed comparison while artificially constraining the antibody count — for example, sweeping `n_pairs ∈ {3, 4, 5, 6, 8, 10}` and selecting each subset from the composite-score rankings in 2.2.7 — to characterise *where* each model's performance degrades and *which metric degrades first*. The reportable result would be a per-model degradation curve per metric, which would tell us (i) whether RF retains its AUC-ROC and log-loss advantage at the device's design limit, (ii) whether one of the simpler models catches up at small panel sizes because the kernel/tree complexity no longer pays off, and (iii) what the lower bound on `n_pairs` is before the probability-based metrics collapse — i.e. before the device can no longer honestly display a confidence number. This makes the in-silico → device handoff concrete in the section immediately preceding the deployment discussion.

*Figures/tables needed:*
- The table above as a clean publication-ready table.
- Per-seed scatter (one dot per seed) for accuracy, κ, AUC-ROC, log loss across the four models — visually demonstrates the SVM-RBF metric divergence and the cross-seed rank stability.
- Reliability diagram / calibration curve for RF vs. SVM-RBF — the single most persuasive figure for the calibration argument.
- Per-class one-vs-rest ROC curves for the four families, with dAb18 and dAb19 panels reported separately.

---

## 2.2.13. Deployment method and considerations

A number of decisions must be made when deploying the final assay.

**Concentration thresholding as a deployment gate.** There is a clear degradation of predictive performance at S-protein concentrations below 1 ng/mL (2.2.9). In such cases, even if a positive signal is obtained, a variant determination should not be provided. This concentration threshold can be determined in one of two ways: using a known SARS-CoV-2 N-protein readout as a proxy for S-protein levels, or using the assay fluorescence intensity directly. The former provides the benefit of a lower LoD readout — a product of comparatively higher concentration of SARS-CoV-2 nucleocapsid protein in clinical samples — and a more robust capacity to diagnose SARS-CoV-2 itself, even in cases where a specific variant classification is not possible. The latter, by contrast, alleviates the need to include and optimize the assay for another antibody, at the cost of slightly diminished functionality at the low-concentration boundary.

**Uncertainty-aware output.** Additional confidence safeguards can be added by displaying the calibrated class probabilities rather than only the argmax. Predictions whose top-class probability falls below a predefined threshold can be flagged — either as assay failures, potentially emerging variants, or unknown anomalies. The choice of RF (or XGBoost, if probability calibration is paramount) over SVM-RBF for this purpose is determined directly by the calibration argument in 2.2.12.

**Algorithm compatibility with new incoming data.** Between-batch normalization is anchored on training-time statistics; new samples are normalized to those statistics, not to their own batch. This matches real device operation (one sample at a time, no batch context available) and is consistent with the Bayesian-update setup in 2.2.10.

**OOD and batch-shift detection.** Motivated directly by the batch-effects results in 2.2.11, the results of any deployed device should be collected and monitored for changes in the distribution of incoming data. If periodic Bayesian updates are then applied based on field-collected samples, the prior — i.e. the data collected during development and training — decays as field data accumulates. While this is the intended Bayesian behaviour, incoming updates must be inspected for potential assay artifacts that could bias the model toward noise or systemic error generated by in-field use rather than true signal. This risk can likely be mitigated by periodically repeating the in-house data-augmentation protocol against current reagent lots, so that the field-update stream is interleaved with controlled-condition data. In the case where an anomalous signal is the result of an emerging variant that is genuinely out of the training-set distribution, a per-sample OOD score (e.g. Mahalanobis distance on the normalized 10-pair vector, a one-class detector, or an evidential / Dirichlet-prior head) can flag the sample for review rather than producing a falsely confident variant call.

**What "deployment" looks like in practice (summary).**

1. Run sample on device → 10-pair fluorescence vector (dAb18 or dAb19 panel).
2. Normalize against training-time statistics.
3. Concentration-threshold check (N-protein proxy or fluorescence-intensity-based): if below 1 ng/mL, return out-of-range rather than a confident class.
4. OOD-score check: if outside training-distribution tolerance, flag for retest and log for the next retraining cycle.
5. RF predicts class probabilities → display variant label and calibrated confidence.
6. Periodically: confirmed field samples are appended to the next Bayesian update of the model, and in-house augmented batches are produced against current reagent lots to keep the field-update stream anchored.

*Figures/tables needed:*
- Deployment-flow diagram (the 6-step list above as a figure).
- Reliability diagram for the deployed RF (carries over from 2.2.12).
- An illustrative OOD-score histogram: in-distribution training samples vs. a held-out batch from 2.2.11.

> **TODO — leave-one-variant-out OOD demonstration.** The OOD safeguard above can be given a concrete in-data demonstration without needing a new instrument-shift dataset, by treating one variant at a time as an unseen OOD class. Procedure: for each variant V in {WT, Beta, Delta, Omicron}, (i) remove all V rows from training; (ii) train a 3-class classifier on the remaining variants together with an OOD-detection head (Mahalanobis distance, one-class SVM, or an evidential / Dirichlet-prior head — pick one and document); (iii) at test time, score the held-out three known variants *and* the held-out V samples; (iv) report (a) classification accuracy on the three known variants and (b) OOD-detection AUC for V-vs-known. Cycling V through all four variants demonstrates that the method works regardless of which variant is treated as novel, and produces a 4-panel figure that directly visualises the deployment-time scenario this chapter has been arguing for. This is the cleanest illustration of the OOD-safeguard claim and the most direct in-paper response to "what happens when a new variant shows up in the field."

---

### Cross-cutting figures/tables (suggest gathering near the methods/results boundary)

- **Master FS-method table:** method · what it ranks · what metric it optimises · why it is in the panel · which composite weight it receives. Resolves the open note that originally lived at the bottom of the outline.
- **Glossary of metrics:** accuracy, Cohen's κ, macro-F1, AUC-ROC, log loss — one sentence each on what they reward and what they penalize, with the argmax-vs-probability grouping made visible. Anchors the throughline for a reader who skims to the methods section first.
