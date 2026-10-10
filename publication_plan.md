# Publication and showcase assessment

Updated using the actual completed run on 10 October 2026. Read [the completed audit](../reports/SCIENTIFIC_AUDIT.md) for verified metrics, baseline and composition interpretation. GitHub portfolio release is appropriate after authorship/license completion; manuscript submission requires further scientific development.

## Scientific contribution worth developing

Working title: **Cross-cohort transfer and selection stability of transcriptomic
classifiers for lung adenocarcinoma**.

Question: does excellent tissue discrimination coexist with unstable gene selection,
and how much does preprocessing or cellular composition affect transfer? This is a
testable reproducibility/robustness question. It is more defensible than claiming that
another combination of known lung-cancer-associated genes is a new diagnostic tool.

The available literature already connects PYCR1 to LUAD functional experiments, GOLM1
to LUAD expression/prognosis, and GPM6A to existing LUAD models. Reusing TCGA or GEO in a
prior paper does not constitute independent corroboration of the present result.
See the biological note and references inside the kit. This was a targeted search,
not a systematic novelty review; the exact proposed contribution must still be compared
against closely related studies.

## Work needed before submission

1. **Extend the saved-run audit.** Verify cohort accounting, complete patient pairs, repeated
   samples, original data checksums, metadata missingness, module status and source/model
   provenance. Confirm that gene selection, scaling and hyperparameter search stayed
   inside the inner folds. Document all prior inspection of GEO results.
2. **Strengthen controls.** Report the already saved QC-only baseline. In a new development
   analysis, compare sensible simple expression baselines and normalization choices on
   TCGA under the same patient-grouped protocol. Repeat nested CV to study selection
   stability; treat this as exploratory development after the original evaluation.
   If using permutation tests, preserve the paired/block design rather than shuffling
   specimen labels indiscriminately across patients.
3. **Fix a new evaluation plan before touching another test cohort.** Choose the biological
   scope, protocol, missing-gene policy and endpoint first. Then evaluate a truly untouched
   suitable cohort. A third cohort is a recommended strengthening step, not a universal
   journal rule. A model revised using GSE40419 cannot use it as an untouched test again.
4. **Match controls to the claim.** For LUAD specificity, include other lung cancer
   histologies and appropriate benign/inflammatory lung samples. For a narrower
   tumor-versus-adjacent-normal claim, state that limitation directly. Reprocessing raw
   reads consistently or testing a suitable orthogonal platform can probe protocol effects.
5. **Resolve the biological interpretation.** Read actual leading-edge pathways; examine
   composition/purity; distinguish known associations from unexpected evidence. Functional
   experiments are appropriate for mechanistic claims but not a universal requirement
   for every computational paper. Avoid adding a large unrelated analysis merely to make
   the paper longer.
6. **Write and archive.** Include sample-flow accounting, exact methods, locked evaluation,
   confusion matrices, uncertainty, failures and limitations. Share source, dependencies,
   necessary derived data and a versioned artifact. Cite original data generators.

Do not change the reported classifier now to remove ANGPTL7 or ADAMTS8. Their instability
is a finding to report. A revised panel needs a new version and evaluation.

## Journal shortlist: possible fit after strengthening

These are scope-based options, not acceptance predictions or endorsements of the current
manuscript. Fees, waivers and specific policies should be checked before submission.

| Journal | Potential fit | Why the present version may be rejected |
|---|---|---|
| PLOS ONE | A technically rigorous secondary analysis or explicitly justified reproducibility study | A derivative signature without a clear scientific rationale, incomplete controls or overstated conclusions |
| PeerJ (Life & Environment; bioinformatics/genomics area) | A well-supported computational biology research article | Thin biological insight, unresolved validation concerns or insufficiently developed research question |
| Scientific Reports | A broader, fully developed original study of cross-cohort robustness | A routine two-cohort classifier with little additional evidence or interpretation |

PLOS ONE explicitly expects justification for studies that replicate or closely resemble
previous work. Its emphasis on scientific validity does not make every executable
analysis publishable. PeerJ includes bioinformatics in scope; Scientific Reports covers
original natural-science research. These scope statements do not guarantee fit for a
specific submission or waive methodological requirements.

BMC Bioinformatics is not my first choice for the unchanged pipeline: its scope emphasizes
computational algorithms, software, models and tools. A genuinely reusable and evaluated
methodological contribution could change that assessment. PLOS Computational Biology
expects novelty and biological significance; the present evidence does not support
targeting it. This is an assessment of fit, not a prediction that editors will reject.

## GitHub versus a research paper

A repository demonstrates software, reproducibility and scientific judgment. A journal
article must communicate a defensible research contribution and survive peer review.
An archived release with a DOI is citable but not peer-reviewed. A preprint is likewise
not an accepted paper. PLOS explicitly permits preprints; check other journals' current
policies for the exact planned disclosure.

For hiring, prioritize a clean repository, a concise result/limitations case study,
actual tests and a short explanation in your own words. Present the external sample
counts and false negatives. Explain AI assistance and what you personally verified.
Do not borrow hospital or university affiliation for work conducted independently.

## Primary sources

- PLOS ONE criteria: https://journals.plos.org/plosone/s/criteria-for-publication
- PLOS ONE data policy: https://journals.plos.org/plosone/s/data-availability
- PeerJ scope: https://peerj.com/about/aims-and-scope/
- Scientific Reports scope: https://www.nature.com/srep/about/aims
- BMC Bioinformatics scope: https://bmcbioinformatics.biomedcentral.com/
- PLOS Computational Biology criteria: https://journals.plos.org/ploscompbiol/s/submission-guidelines
- Public-data research guidance: https://doi.org/10.1371/journal.pcbi.1010749
- PYCR1 original study: https://pubmed.ncbi.nlm.nih.gov/32133692/
- GOLM1 original study: https://doi.org/10.1177/1073274818778001
- LUAD multi-omics study including GPM6A: https://doi.org/10.1186/s12885-025-14943-x
- GitHub import instructions: https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github
- Zenodo: https://help.zenodo.org/docs/get-started/quickstart/
