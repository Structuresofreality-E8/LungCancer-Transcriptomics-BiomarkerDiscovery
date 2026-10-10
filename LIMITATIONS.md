# Interpretation and next evidence required

This experiment classifies resected LUAD tissue versus adjacent nonmalignant
lung tissue. Adjacent normal is not a healthy screening population. Tissue
discrimination does not establish blood detectability or preclinical detection.

Mitigated: patient leakage (all folds grouped), global supervised feature
selection (fold-local), cross-study units (per-sample gene ranks), duplicate
aliquots (fixed lexical choice), DE pairing (complete-pair design), multiplicity
(BH), and external test tuning (locked model). Ranks do not remove all platform,
annotation, gene-length, cellular-composition, or population differences.

Patient bootstrap intervals describe test-patient sampling uncertainty at a
fixed trained model; nested OOF intervals do not capture full training-process
uncertainty. Average precision depends on prevalence. Class weighting means
scores and Brier values are not calibrated clinical probabilities. Threshold
0.5 was fixed before external testing. Bootstrap all model training stages or
repeat nested CV for a stronger, more expensive uncertainty analysis.

Signature stability is a descriptive frequency across outer training folds;
it is not a probability that a gene is causal. Replication p-values are BH
adjusted over the locked selected signature only. External log2(RPKM+1)
paired differences are not directly comparable to DESeq2 log2 fold changes.
GSEA uses gene-set permutations; gene correlation can affect significance.

Bulk immune/stromal/epithelial marker averages are confounding diagnostics,
not estimated cell fractions or purity. Confirm using pathology-derived purity,
validated deconvolution, and independent single-cell or spatial data. Do not
interpret signature genes as tumor-cell-intrinsic or causal without that work.

Age, smoking, library preparation, and medication may be incompletely captured.
GDC tissue source site is a center proxy, not a sequencing batch. Sparse
subgroup estimates are omitted, not proof of fairness. Basic STAR summary
metrics are not FASTQ QC; raw-read processing is required for read-level checks.

Next study: lock the signature, assay and threshold; test blinded in additional
independent tissue cohorts including benign nodules, inflammatory disease,
other cancers and metastatic lesions. For a blood-screening claim, assay the
markers in blood, assess detection limits and tissue-to-blood transport, then
validate in a prospective intended-use screening cohort with healthy controls,
benign disease, appropriate prevalence and stage distribution. Measure
sensitivity at prespecified high specificity, predictive values, calibration,
and clinical benefit. Neither TCGA nor GSE40419 can supply that evidence.

Stage-I tissue sensitivity is descriptive and does not demonstrate early
detection. Pathway/disease annotation is contextual evidence, not clinical
validation or novelty. No causal, diagnostic, prognostic, or treatment claim.
