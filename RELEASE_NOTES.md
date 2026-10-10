## v1.0.0 — LUAD transcriptomic reproducibility study

First portfolio release of a retrospective TCGA-to-GSE40419
lung adenocarcinoma tissue-classification analysis.

### Results
- TCGA: 576 samples from 518 patients.
- External GSE40419: 164 samples from 87 patients.
- External ROC AUC: 0.9927.
- Conditional patient-bootstrap AUC 95% interval: 0.9799–0.9997.
- External confusion matrix: TN=77, FP=0, FN=7, TP=80.
- TCGA QC-only baseline ROC AUC: 0.6685.

### Included
Source, Colab notebook, environment records, input manifests, saved model,
result tables, figures, uncertainty estimates, software checks and
a Streamlit results explorer.

### Interpretation
All 25 selected genes meet the implemented paired association-replication
criteria, but gene selection is not uniformly stable. Model scores correlate
with composition proxies. The rank transformation requires a fixed universe
of 16,937 genes; this is not a standalone 25-gene assay.

This is a research tissue classifier, not a validated screening test.
Saved artifacts cannot certify the complete history of external-data
inspection. See RUN_HISTORY.md and reports/SCIENTIFIC_AUDIT.md.
