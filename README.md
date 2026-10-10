# LUAD transcriptomic transfer and reproducibility

A retrospective study of whether a sparse gene-expression classifier transfers
from TCGA lung adenocarcinoma tissue to the independent GSE40419 cohort.

Status: saved-run audit completed; portfolio release candidate. Independent scientific review is recommended before journal submission.

**Research scope:** resected tumor versus adjacent nonmalignant lung tissue.
The outcome is tissue classification, not screening, prognosis, treatment response,
or blood-based early detection.

<!-- RESULTS_START -->
## Results verified against saved predictions

| Evaluation | n | roc_auc | average_precision | balanced_accuracy | sensitivity | specificity |
| --- | --- | --- | --- | --- | --- | --- |
| TCGA nested OOF | 576 | 0.998492 | 0.999835 | 0.98669 | 0.990329 | 0.983051 |
| GEO locked external | 164 | 0.992685 | 0.994681 | 0.95977 | 0.91954 | 1 |
| TCGA QC-only baseline | 576 | 0.668492 | 0.942995 | 0.642658 | 0.60735 | 0.677966 |

External confusion matrix: TN=77, FP=0, FN=7, TP=80.
TCGA includes 518 patients; GEO includes 87 patients.
25/25 selected genes satisfy the implemented paired association-replication criteria.

See [uncertainty, baseline, pathway and composition audit](reports/evidence_review.md).
The classifier uses 25 genes after percentile normalization across 16937 genes;
it is not a stand-alone 25-gene assay.

<!-- RESULTS_END -->

## Audit findings

External AUC 95% patient-bootstrap interval: **0.9799–0.9997**. Exact external sensitivity interval: **84.12–96.70%**; specificity: **95.32–100%**. QC-only AUC is **0.6685**. GEO composition proxies correlate with scores, and two final genes were absent from all five outer-fold signatures. These limits are part of the result.

Read the [completed scientific audit and publication assessment](reports/SCIENTIFIC_AUDIT.md).

![External confusion matrix](figures/external_confusion_matrix.png)

## Why this project is useful

- Patient-grouped 5-fold outer / 3-fold inner cross-validation.
- Fold-local variance filtering, ANOVA feature ranking, scaling and elastic-net fitting.
- Primary differential expression on complete TCGA pairs with patient + condition design.
- Fixed gene availability across cohorts, length-normalized expression and per-sample ranks.
- A model and threshold locked before external scoring.
- Gene association replication, selection stability, uncertainty and confounding diagnostics.

This is an application of established methods. No new learning algorithm or previously
unknown causal biomarker is claimed. The estimator is sparse, but rank normalization
requires the entire fixed shared gene universe, not only the 25 coefficient genes.

## Reproduce

1. Upload `LUAD_Biomarker_Discovery_Colab.ipynb` to Colab and select CPU.
2. Install its pinned dependencies. Restart the Python session after installation
   if scientific libraries were already imported; then resume without reinstalling.
3. Run the source-writing cell, optional Drive cell, and pipeline cell sequentially.
4. Wait for completion before inspecting or downloading results.

Local execution (Python 3.12 tested for the software smoke test):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export LUAD_ROOT="$PWD/run_outputs"
python -c "import luad_pipeline as pipeline; pipeline.main()"
python -m unittest discover -s tests -v
```

No GPU or GDC access token is used. The full eligible public cohort is downloaded.
The GDC manifest pins file identities/checksums; the actual run's package snapshot
is in `requirements-lock.txt` after export. Colab availability and download speeds vary.

## Inspect the completed release

- `reports/evidence_review.md`: saved-file checks, uncertainty, baseline, pathways and composition.
- `results/locked_signature.csv`: standardized coefficients and descriptive selection frequency.
- `results/external_paired_gene_replication.csv`: external paired gene tests.
- `results/metrics.json`: original run metrics and saved patient-bootstrap intervals.
- `figures/`: original analysis figures plus an external confusion matrix.
- `docs/biology_interpretation.md`: literature context, with limits on interpretation.
- `docs/model_card.md`: preprocessing, intended use, and failure modes.
- `docs/showcase_and_release.md`: GitHub and portfolio presentation.

After completing the export, run `streamlit run app.py` from the repository root.
The Streamlit app examines saved results; it does not rerun model fitting.
See `docs/demo_script.md` for a short recording outline.

## Important limits

The two cohorts differ in protocols, population and expression processing. Ranks do
not remove every study effect. Adjacent normal tissue is not a healthy screening
population. High discrimination may reflect tissue/cell composition. Marker proxies
are not deconvolution or measured purity. Scores are not calibrated clinical risk.
There are only 77 external normal samples; zero observed errors is not perfect population
specificity. Bootstrap uncertainty is conditional on the fitted models.

Keep this evaluation frozen. Record any previous GEO-informed iteration honestly.
Changes motivated by these results define a new model version requiring another
untouched external cohort. Saved artifacts cannot prove the absence of human tuning.

## Data, attribution and authorship

[TCGA/GDC](https://portal.gdc.cancer.gov/),
[GSE40419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40419),
[PyDESeq2](https://pydeseq2.readthedocs.io/en/v0.5.2/),
[Reactome](https://reactome.org/), and [Open Targets](https://platform.opentargets.org/).
Full references and reuse notes: `docs/references.md` and `docs/data_and_authorship.md`.

Initial pipeline code and repository packaging were drafted with ChatGPT/Codex assistance.
The maintainer is responsible for execution, verification, interpretation and the final
scientific claims. No hospital collaboration or newly collected clinical data is implied.
Before public release, enter the maintainer's real authorship details and choose a code
license. Do not use an institution as an affiliation without a valid relationship.


## Maintainer and citation

Maintained by [Tarun Malik](https://github.com/Structuresofreality-E8).
Citation metadata: [CITATION.cff](CITATION.cff).
Code license: [MIT](LICENSE). Third-party datasets, annotations and
dependencies retain their own terms; this license does not relicense them.

## Featured validation figure

![Internal and external validation curves](figures/validation_curves.png)

External ROC AUC: **0.9927**; 95% conditional patient-bootstrap interval:
**0.9799–0.9997**. Read the uncertainty and composition limitations above.

## Demonstration

Run `streamlit run app.py` after installing the project dependencies.
A recording outline is provided in [docs/demo_script.md](docs/demo_script.md).
