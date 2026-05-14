# Thesis Figure Captions

Journal-style captions for every figure in `Thesis_figures/figs/`. Each entry
gives a bold title sentence, the methodological details a reader needs to
interpret the plot, and the takeaway. Drop directly into LaTeX or Word; adjust
figure numbering to match the chapter that hosts each plot.

---

### Figure 1 -- 01_cross_method_perf_grid.png

**Classifier performance as a function of selected feature count for six
feature-selection (FS) methods.** Mean test-set accuracy (light blue),
macro-averaged F1 score (medium blue), and Cohen's kappa (dark blue) are
plotted against the number of antibody pairs retained in the model, with
+/-1 standard error bars aggregated across three random seeds (1969, 1986,
2024) and all 13 detection antibodies. The x-axis is reversed so that fewer
pairs (the deployment-relevant regime) appear at right. All six FS methods
plateau between 6 and 10 pairs and degrade meaningfully below 5 pairs;
PLS-DA loses the most performance over this range, while Boruta and random
forest (RF) are most stable. The 5-6 pair regime defines the engineering
sweet spot for device deployment.

---

### Figure 2 -- 02_method_dab_acc_heatmap.png

**Final-test classification accuracy by feature-selection method and
detection antibody.** Each cell shows the mean test-set accuracy across three
random seeds at each method's auto-selected best number of features. Cell
shading saturates from white at chance (0.50) to dark blue at perfect
classification (1.00). Detection antibody dAb5 yields the highest accuracy
under every FS method (range 0.83-0.92), whereas dAb16 is the most
method-dependent antibody (LASSO 0.75 vs RF 0.62, a 13-point spread)
because RF systematically misclassifies Delta as Beta on dAb16's feature
subset (Fig. 4). The heatmap simultaneously identifies universally
informative dAbs and dAbs whose informativeness is contingent on the
selection algorithm.

---

### Figure 3 -- 03_feature_frequency_heatmap.png

**Selection frequency of capture antibodies (cAbs) as final-model features
across all FS runs.** Cells count, for each (detection antibody x capture
antibody) pair, the number of (FS method x seed) runs that retained that
pair in the final best-feature subset; the maximum possible count is 18
(6 methods x 3 seeds). White cells indicate pairs that were never selected
by any method or seed. cAb25 and cAb26 emerge as near-universal partners
across the eight top-performing detection antibodies (counts 13-18 in most
rows); cAb13, cAb14 and cAb15 form a secondary consensus cluster. dAb8 is
unusual in selecting cAb1 (17/18) as its top partner. Single-digit counts
identify pairs whose selection is method- or seed-specific and therefore
unreliable for the final assay design.

---

### Figure 4 -- 04_confusion_matrices_rf_top.png

**Random-forest test-set confusion matrices for the eight top detection
antibodies.** Each panel shows raw counts and row percentages of true
(rows) versus predicted (columns) variant calls, with the per-cell counts
summed over three random seeds at each detection antibody's auto-selected
best number of features. Cell shading denotes row-normalised proportion.
Detection antibodies dAb5, dAb6, dAb13 and dAb18 produce nearly diagonal
confusion matrices (>=80% per-class recall on all four variants), confirming
clean four-way separation. The off-diagonal mass on dAb16 -- 8 of 15 Delta
samples are called Beta -- explains its low overall accuracy in Fig. 2.
dAb19 achieves perfect Delta recall (15/15) at the cost of fragmented WT
performance (5/15 Beta).

---

### Figure 5 -- 05_per_class_metrics_rf.png

**Per-variant precision, recall and F1 score for random-forest classifiers
on each top detection antibody.** Bars show class-level metrics derived from
the confusion matrices in Fig. 4, with the value annotated above each bar.
The per-class view exposes failure modes that the macro-averaged accuracy
in Fig. 2 obscures. dAb16's Delta F1 is 0.55 (precision = 0.50, recall =
0.67) despite an overall accuracy of 0.62, identifying Delta detection as
the specific weakness of the dAb16 model. dAb13 achieves perfect Delta
recall (1.00) but only 0.79 Delta precision because Beta and Omicron leak
into the Delta column. WT and Omicron are the easiest classes overall
(median F1 ~0.83 across the eight panels).

---

### Figure 6 -- 06_dose_response_top_pairs.png

**Dose-response curves for the highest-LASSO-importance capture-antibody
partner of each top detection antibody.** Lines and ribbons show the mean
+/- one standard deviation of fluorescence across the three replicates at
each concentration; the x-axis is on a log10 scale. Variants are
distinguished by colour. Most top pairs only separate the four variants
above approximately 10 ng/mL, defining a practical lower bound on the
useful clinical concentration range. dAb13 x cAb26 is the only pair in
which Beta separates from the other three variants across the full
concentration range, accounting for its prominence in the LASSO-selected
features for Beta detection.

---

### Figure 7 -- 07_replicate_cv_ranking.png

**Replicate reproducibility ranking of antibody pairs.** Coefficient of
variation (CV = standard deviation / mean) was computed across the three
technical replicates at each (concentration x variant) cell, then averaged
across all cells per pair. Bars show the 25 most reproducible (left) and
25 least reproducible (right) pairs across all 13 detection antibodies and
their 29 capture partners; bar colour denotes the parent dAb. All pairs
fall below 6% CV, supporting the use of raw replicate-level fluorescence
readings as model inputs. The most reproducible pairs (e.g. dAb8 x cAb20,
dAb6 x cAb10) are not the same as the most discriminative pairs in
Figs. 3 and 12; reproducibility and class discriminative power are
distinct properties.

---

### Figure 8 -- 08_pca_per_dab.png

**Principal-component projections of the 29-column capture-antibody
fluorescence matrix for each top detection antibody, coloured by variant.**
For each detection antibody, principal-component analysis (PCA) was run on
the centred and scaled fluorescence matrix across all replicates and
concentrations; the first two principal components are plotted with 68%
confidence ellipses by variant. Variance explained by PC1 and PC2 is
annotated in each panel header. PC1 captures most of the variance
(typically >=80%) and aligns with concentration; variants separate primarily
along this axis with secondary structure on PC2. The figure demonstrates
that variant identity is recoverable from raw fluorescence readouts before
any supervised feature selection, motivating the per-dAb classifier
strategy.

---

### Figure 9 -- 09_augmented_vs_original.png

**Distribution overlays of original screen data versus the augmented
dataset for detection antibody dAb16.** Density curves of fluorescence
readings are shown for each (capture antibody x variant) combination, with
the original screen in blue and the augmented dataset in orange. The
augmented file maps onto dAb16 (11 overlapping capture antibodies). The
augmentation procedure preserves the modal fluorescence values for most
cells, but introduces a heavy low-fluorescence tail not present in the
original screen, particularly for the Beta and Delta variants. This
distributional drift should be characterised before relying on augmented
data for Phase 2 model validation.

---

### Figure 10 -- 10_roc_curves_rf_top.png

**One-versus-rest receiver operating characteristic (ROC) curves for
random-forest models on each top detection antibody.** Each panel shows
sensitivity (true positive rate) plotted against 1 - specificity (false
positive rate) for each variant treated as the positive class. Curves are
constructed by stacking probability predictions across three seeds on each
model's held-out final-test split. The dashed grey diagonal denotes
chance performance. Curves cluster in the upper-left corner for dAb5,
dAb13 and dAb18, indicating near-perfect class separation at deployable
feature counts. dAb16's curves are visibly looser, mirroring its
classification accuracy and Delta confusion in Figs. 2 and 4. Omicron is
the easiest variant to discriminate across all detection antibodies.

---

### Figure 11 -- 11_calibration_rf.png

**Probability calibration of the pooled random-forest classifier across
top detection antibodies.** Predicted probabilities were binned in 0.1
increments (bin sample size encoded by point area, minimum 5
observations per bin) and the observed positive rate within each bin is
plotted against the mean predicted probability for each variant. The
dashed grey diagonal denotes perfect calibration. WT predictions sit
slightly above the diagonal at moderate probabilities and Beta sits
slightly below in the 0.2-0.4 range, but all four variants track close
to perfect calibration. The pattern is consistent with the global
temperature-scaling result in Fig. 11b.

---

### Figure 11b -- 11b_temp_scale_T_k.png

**Per-class temperature parameters fitted independently for each (seed x
detection antibody) random-forest model.** Box plots show the
distribution of the four learned scalar temperatures T_k (one per
variant), with each point representing one of 40 fits (5 seeds x 8 top
detection antibodies). The y-axis is on a log scale; the dashed line at
T = 1 marks no scaling. Temperatures were fitted by minimising negative
log-likelihood on five-fold cross-validated out-of-fold predictions
within each fit's training partition. All four classes cluster tightly
around T = 1, indicating that per-instance calibration adjustments are
small. WT and Beta lean marginally above 1.0 (mild flattening), Delta
and Omicron marginally below 1.0 (mild sharpening), consistent with the
qualitative pattern in Fig. 11.

---

### Figure 11c -- 11c_temp_scale_ece_before_after.png

**Expected Calibration Error (ECE) before and after per-instance
temperature scaling, for each of the 40 (seed x detection antibody)
fits.** Each line connects the uncalibrated and calibrated ECE for one
fit; line colour denotes the detection antibody. Per-instance
calibration on this sample size yields no systematic improvement: the
mean ECE moves from 0.184 +/- 0.060 (uncalibrated) to 0.196 +/- 0.048
(calibrated), with calibration helping roughly half of the fits and
hurting the other half. Pooling out-of-fold predictions across all 40
fits and refitting one global per-class T_k (Section X.X) reduces ECE
from 0.078 to 0.077, confirming that the random-forest classifier is
already near-calibrated and that further calibration on this dataset
provides negligible benefit.

---

### Figure 12 -- 12_lasso_pair_importance.png

**Top eight capture-antibody partners per top detection antibody ranked by
mean LASSO importance across three seeds.** Bars show the mean absolute
LASSO coefficient (averaged across the four variant classes) for each
capture-antibody partner; the numeric value is annotated to the right of
each bar. Bar fill encodes how many of three seeds retained the pair in
the final best-feature subset (light = 0/3, dark = 3/3). Error bars show
+/-1 standard deviation across seeds. dAb16 x cAb14 is the single most
important pair (mean |coefficient| = 2.43); dAb18 x cAb25, dAb19 x cAb26
and dAb13 x cAb26 are the next strongest individual contributors. The
fill encoding identifies the stable feature core of each model; for
example dAb18's top three partners (cAb25, cAb6, cAb13) all appear in
3/3 seeds, while dAb19's top features are more variable.

---

### Figure 13 -- 13_beeswarm_top_dabs.png

**Per-variant fluorescence distributions for the eight top detection
antibodies, with all capture-antibody partners pooled.** Each point
represents one fluorescence reading from a single (capture-antibody
partner x replicate x concentration) measurement for the indicated
detection antibody. Points are positioned along a quasi-random horizontal
jitter and coloured by variant; black crossbars mark per-variant medians.
The spread of medians across variants foreshadows classification
difficulty: detection antibodies whose median fluorescence is well
separated by variant (dAb5, dAb13) yield high per-dAb accuracy in Fig. 2,
while dAb16's compressed median spread is consistent with its weaker
downstream classification.

---

### Figure 14 -- 14_performance_summary.png

**Summary of feature-selection results across all six FS methods. (A)**
Mean classifier performance (accuracy, macro-F1, kappa) versus the number
of selected antibody pairs, by FS method, with +/-1 standard error bars
across three seeds and all 13 detection antibodies. **(B)** Mean test
accuracy by FS method and detection antibody at each method's
auto-selected best number of features. Together, panels A and B identify
both the engineering sweet spot for the deployable assay (5-6 antibody
pairs) and the detection antibodies that drive classification quality
(dAb5, dAb13, dAb18, dAb19). Reproduced from Figs. 1 and 2.
