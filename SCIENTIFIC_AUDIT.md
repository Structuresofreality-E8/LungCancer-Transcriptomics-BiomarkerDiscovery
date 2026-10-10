# LUAD saved-run audit — 10 October 2026

## Assessment

Suitable for a transparent GitHub portfolio release after adding real authorship and a chosen code license. The saved outputs support strong cross-cohort tissue discrimination. They do not establish a clinically validated assay or a novel mechanism. A research manuscript needs a sharper robustness question and additional controls.

## Verified findings

| Evaluation | Patients | Samples | ROC AUC | Sensitivity | Specificity |
| --- | --- | --- | --- | --- | --- |
| TCGA nested out-of-fold | 518 | 576 | 0.99849 | 512/517 | 58/59 |
| GSE40419 external | 87 | 164 | 0.99269 | 80/87 | 77/77 |

TCGA contains 58 complete tumor/normal pairs; GEO contains 77. No repeated patient-condition records were found in prediction tables. All observations belonging to each TCGA patient share an outer fold. Metrics and confusion matrices were independently recalculated from saved predictions. All 78 original release hashes matched. Four supplied software invariant tests passed. These checks do not constitute a rerun of downloaded expression data or verification of the original human analysis history.

External AUC patient-bootstrap 95% interval: 0.97988–0.99971, conditional on the fixed fitted model. External sensitivity exact 95% interval: 84.12–96.70%; specificity: 95.32–100%. The zero-width specificity bootstrap interval is an artifact of observing no false positives; use the exact interval alongside it.

## What changes the interpretation

**Technical summaries contain signal.** The TCGA QC-only baseline AUC is 0.66849 (patient-bootstrap 95% interval 0.60246–0.73946). Expression exceeds this baseline by 0.33000 AUC (paired interval 0.25940–0.39422). Thus simple library summaries alone do not reproduce expression performance, but technical confounding has not been ruled out. QC average precision 0.943 should be interpreted against the high TCGA tumor prevalence of 0.898.

**Cellular composition remains a serious alternative explanation.** Within GEO tumors, score correlates negatively with the immune proxy (Spearman rho −0.592; BH q=4.75e−9) and positively with the epithelial proxy (rho 0.550; q=6.96e−8). Within GEO normals, score correlates positively with the stromal proxy (rho 0.664; q=2.79e−10). TCGA tumor associations are weaker (absolute rho about 0.11–0.13). These are marker-rank correlations, not measured cell fractions. Shared expression input and biological programs can also induce correlation. The findings justify purity/reference-based follow-up; they do not demonstrate that composition fully explains the classifier.

**Association replication is stronger than selection stability.** All 25 selected genes meet the implemented paired direction/FDR replication criteria. ANGPTL7 and ADAMTS8 appear in the final model but in zero of five outer-fold models. Their replicated differential expression does not make them indispensable predictors. Five-fold frequency is descriptive and coarse, not a confidence probability.

**A sparse estimator is not a 25-gene assay.** The 25 coefficients operate on percentile ranks requiring the fixed 16,937-gene measurement universe. Measuring only those 25 genes would change the inputs.

**Pathway evidence is primarily proliferation-related.** There are 462 Reactome terms at reported GSEA q<0.05. Representative positive terms include Cell Cycle Checkpoints (NES 3.315), G2/M Checkpoints (3.225) and DNA Replication (3.060). Leading genes include CDC25C, CCNB1, CCNB2, CHEK1, CDK1 and replication machinery. Many terms share genes and hierarchy; this is not 462 independent discoveries. Reported q=0 values are finite-permutation software estimates, not zero error probabilities. Pathways were tested in TCGA only and have not been independently replicated as pathways in GEO.

## GitHub presentation

1. Extract the ZIP. Create a repository named LungCancer-Transcriptomics-BiomarkerDiscovery. Upload the contents of the extracted folder, so README.md is at the root; do not upload only the ZIP.
2. Enter your real details in CITATION.cff.example and save the completed file as CITATION.cff. Select a code license and add LICENSE. Complete RUN_HISTORY.md honestly, including any prior GEO-informed changes.
3. Feature the external confusion matrix, validation curves and pathway plot. Retain uncertainty and the composition caveat near the headline performance.
4. Use docs/demo_script.md for a 90-second walkthrough; run `streamlit run app.py` locally after installing dependencies. A recording has not been created by this audit.
5. Tag the reviewed snapshot v1.0.0 and pin the repository to your profile. Link the repository and a short demonstration from your CV or portfolio. A versioned archive/DOI is useful for citation but is not peer review.

Suggested portfolio wording: “Built and evaluated a reproducible LUAD tissue-expression classifier using patient-grouped nested validation and an independent RNA-seq cohort. External ROC AUC 0.993; quantified uncertainty, technical baseline performance, selection instability and composition-related limitations.” Use this only if you can explain and defend the implementation and your contribution.

## Publication decision

Public-data secondary analysis can support a research paper; the absence of new wet-lab data is not an automatic disqualification. This result alone is not enough to claim a new clinically useful biomarker panel. The strongest paper framing is: **Cross-cohort transfer, selection stability and composition sensitivity of LUAD transcriptomic classifiers**.

Before submission: compare appropriate simple expression baselines under identical grouped validation; assess repeated-resampling stability; investigate composition with suitable independent purity/reference evidence; document novelty against related work; specify a new untouched cohort before evaluating any revised model. A third cohort is a recommended strengthening step, not a universal journal requirement. Mechanistic claims need appropriate experimental evidence; the present narrow association/classification claim does not automatically require wet-lab experiments.

Conditional journal options:
- PLOS ONE: plausible for a rigorously justified reproducibility/secondary-analysis study; its criteria explicitly discuss replication studies. https://journals.plos.org/plosone/s/criteria-for-publication
- PeerJ: relevant biological research scope; fit depends on the final contribution and rigor. https://peerj.com/about/aims-and-scope/
- Scientific Reports: broad original-research scope; consider once robustness and contribution are convincingly developed. https://www.nature.com/srep/about/aims

These are scope-based options, not acceptance predictions. Review current author requirements and charges before choosing. No submission or external publication was made.

## Provenance and limits of this audit

Pipeline, notebook, fitted model, original prediction tables, scientific results and figures were preserved byte-for-byte. Documentation was updated only. The original manifest is retained as reports/ORIGINAL_RELEASE_MANIFEST.json; the new root manifest covers this reviewed package. No retraining, threshold adjustment or GEO optimization occurred in this audit. Saved files cannot establish that GEO was untouched in every earlier human decision. Raw input checksum records are present, but bulk raw downloads were not supplied and were not independently rehashed. A bioinformatician should review the full study before manuscript submission.
