# Biological interpretation: what is supported now

The saved signature associates PYCR1 and GOLM1 positively with tumor classification,
and GPM6A negatively. The saved paired tests also show matching paired expression directions
across cohorts. These observations are associations; bulk RNA cannot identify the cell
of origin or the cause of the change.

## PYCR1: existing biological context

Gao et al. (2020) reported that PYCR1 knockdown reduced proliferation, migration and
invasion in LUAD cell lines and investigated JAK/STAT signaling (PMID 32133692;
https://doi.org/10.1002/mc.23174). The present positive association is consistent with
prior cancer biology. It does not reproduce that mechanism or establish a new target.

## GOLM1: prior lung adenocarcinoma association

Liu, Chen and Zhang (2018) analyzed public cohorts and reported elevated GOLM1 and
prognostic associations in LUAD (https://doi.org/10.1177/1073274818778001).
The current classification result is compatible with existing expression evidence.
No survival analysis was performed here. The prior study may reuse TCGA/GEO resources,
so it is literature context, not an additional independent validation cohort.

## GPM6A: a previously studied candidate

A 2025 multi-omics LUAD study included GPM6A among genes downregulated in tumor tissue
(https://doi.org/10.1186/s12885-025-14943-x). Its appearance here is not evidence that
this project discovered the gene's association. Tissue composition, stromal variation
and differentiation are hypotheses to investigate, not demonstrated mechanisms.

## Pathways and composition from the completed run

TCGA enrichment is dominated by positive cell-cycle and DNA-replication signals: Cell Cycle Checkpoints NES 3.315, G2/M Checkpoints 3.225 and DNA Replication 3.060. This is consistent with proliferation-associated transcription, not proof of a novel mechanism. Related terms share many leading-edge genes; 462 significant terms do not represent 462 independent discoveries. Reported zero permutation q-values must not be interpreted as certainty. These pathways have not been validated as pathways in GEO.

GEO tumor scores correlate negatively with immune marker ranks (rho −0.592) and positively with epithelial ranks (rho 0.550); GEO normal scores correlate positively with stromal ranks (rho 0.664). All three have BH q<1e−7 in the saved diagnostics. TCGA tumor correlations are weaker. Composition and shared biological programs therefore remain plausible contributors. Marker averages are not cell fractions; bulk RNA alone cannot assign these signals to malignant cells. See ../reports/SCIENTIFIC_AUDIT.md for numerical context and follow-up priorities.

## Candidate selection

ANGPTL7 and ADAMTS8 were selected by the final model but in none of the five outer-fold
models. Do not delete them now to improve the reported external result. Document this
instability. Any stability-filtered panel is a new model version.

This is a targeted contextual literature check, not a systematic novelty review.
Before submission, examine original papers, publication updates and overlap in datasets.
