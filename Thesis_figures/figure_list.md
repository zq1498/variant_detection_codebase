# List of Figures

For thesis front matter. Title sentences are the bolded leads from
`figure_captions.md`, trimmed to one line. Renumber to match final
chapter ordering.

| # | File | Title |
|---|------|-------|
| 1 | `figs/01_cross_method_perf_grid.png` | Classifier performance as a function of selected feature count for six feature-selection methods. |
| 2 | `figs/02_method_dab_acc_heatmap.png` | Final-test classification accuracy by feature-selection method and detection antibody. |
| 3 | `figs/03_feature_frequency_heatmap.png` | Selection frequency of capture antibodies as final-model features across all FS runs. |
| 4 | `figs/04_confusion_matrices_rf_top.png` | Random-forest test-set confusion matrices for the eight top detection antibodies. |
| 5 | `figs/05_per_class_metrics_rf.png` | Per-variant precision, recall and F1 score for random-forest classifiers on each top detection antibody. |
| 6 | `figs/06_dose_response_top_pairs.png` | Dose-response curves for the highest-LASSO-importance capture-antibody partner of each top detection antibody. |
| 7 | `figs/07_replicate_cv_ranking.png` | Replicate reproducibility ranking of antibody pairs. |
| 8 | `figs/08_pca_per_dab.png` | Principal-component projections of the 29-column capture-antibody fluorescence matrix for each top detection antibody. |
| 9 | `figs/09_augmented_vs_original.png` | Distribution overlays of original screen data versus the augmented dataset for detection antibody dAb16. |
| 10 | `figs/10_roc_curves_rf_top.png` | One-versus-rest receiver operating characteristic curves for random-forest models on each top detection antibody. |
| 11 | `figs/11_calibration_rf.png` | Probability calibration of the pooled random-forest classifier across top detection antibodies. |
| 11b | `figs/11b_temp_scale_T_k.png` | Per-class temperature parameters fitted independently for each (seed x detection antibody) random-forest model. |
| 11c | `figs/11c_temp_scale_ece_before_after.png` | Expected calibration error before and after per-instance temperature scaling, for each of the 40 (seed x detection antibody) fits. |
| 12 | `figs/12_lasso_pair_importance.png` | Top eight capture-antibody partners per top detection antibody ranked by mean LASSO importance across three seeds. |
| 13 | `figs/13_beeswarm_top_dabs.png` | Per-variant fluorescence distributions for the eight top detection antibodies, with all capture-antibody partners pooled. |
| 14 | `figs/14_performance_summary.png` | Summary of feature-selection results across all six FS methods. |

---

## Suggested chapter placement

Map figures onto the standard thesis structure to check coverage:

- **Methods / Materials** -- Fig. 7 (assay reproducibility justifies raw-replicate inputs)
- **Results -- raw screen** -- Figs. 13, 6, 8 (what the data looks like; pre-classification structure)
- **Results -- model performance** -- Figs. 14 (= 1 + 2), 3, 12 (chapter-opener summary; selected features)
- **Results -- error analysis** -- Figs. 4, 5, 10, 11, 11b, 11c (where models fail; discrimination quality; calibration story)
- **Discussion -- caveats and future work** -- Fig. 9 (augmentation drift, motivates Phase 2 validation)
