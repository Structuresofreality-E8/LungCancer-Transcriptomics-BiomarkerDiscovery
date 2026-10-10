# Model card

Task: classify a resected lung specimen as LUAD tumor or adjacent nonmalignant tissue.
Training resource: open TCGA-LUAD STAR counts with FPKM and annotations.
External resource: GSE40419 RPKM, 87 tumors / 77 adjacent normals.

Development protocol: 5 outer patient-grouped folds, 3 inner patient-grouped folds;
variance filtering, ANOVA selection (10/25/50), standardization and class-weighted
elastic-net logistic regression. C = 0.01/0.1/1; l1_ratio = 0.5/0.9. Inner ROC-AUC
selects a candidate. Final fitting uses all TCGA samples. Threshold: 0.5.

Every prediction requires expression for the saved fixed gene universe. Rank each
sample over that universe with average ties, then apply the fitted pipeline.
Do not rank only the selected genes, independently scale the GEO cohort, zero-fill
missing selected genes, or substitute raw counts for length-normalized expression.
Future gene panels require a separately developed normalization scheme and validation.

Association analysis: integer unstranded counts; complete TCGA pairs; patient +
condition design. GEO tests: paired Wilcoxon of log2(RPKM+1); BH over the selected
signature. These quantities are on different scales. Logistic coefficients are not
DESeq2 log fold changes, and they are not causal effects.

Artifact `locked_model.joblib` was produced by Python serialization. Load only a
trusted copy, with the matching source module and environment. The exporter never
loads this artifact and never calls model.fit.

Limitations: spectrum and selection bias; surgical specimens; geography and protocol
shift; sparse demographic controls; potential cell-composition effects; incomplete
confounder measurements; conditional bootstrap uncertainty; uncertain panel stability.
No claim of clinical deployment, diagnostic approval, early blood detection, prognosis,
causality or treatment benefit. See actual run's LIMITATIONS.md for further details.
