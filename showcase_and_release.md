# GitHub and portfolio release

1. In the completed Colab session, run the supplied exporter. It copies real outputs,
   checks metrics against saved predictions, adds uncertainty and baseline reports,
   and creates a lightweight repository ZIP. It performs no retraining.
2. Download and extract the resulting LUAD_GitHub_Ready ZIP. Review README.md and
   reports/evidence_review.md. Resolve failed modules; do not describe them as complete.
3. Create a GitHub repository named `luad-transcriptomic-transfer-validation`.
   GitHub Desktop is the simplest upload route: add the extracted folder as a local
   repository, review the staged files, commit, and publish. Upload source, docs,
   results, figures and manifests; exclude bulk expression files and secrets.
4. Add a concise About description: "Patient-grouped transcriptomic modeling with
   independent LUAD tissue-cohort evaluation and reproducibility diagnostics."
   Suggested topics: bioinformatics, rna-seq, lung-adenocarcinoma, reproducible-research,
   scikit-learn, pydeseq2. Pin the repository on your profile.
5. Run the tests. Verify reproduction in a clean environment. Include a screenshot
   of the results explorer and a 60–90 second walkthrough using docs/demo_script.md.
6. Fill authorship/citation details, select a license for your code, and preserve
   upstream data/software licenses. Tag a release such as v0.1.0 after the audit.
7. Optionally archive that version in Zenodo and put its assigned DOI in the README
   and citation metadata. A DOI is an archival identifier, not peer review.

Suggested README image order: external confusion matrix, validation curves, then
signature/replication figure. Show four decimal places in the summary table so 0.9985
is not presented as perfect performance. Always show the external sample counts.

Portfolio presentation: link the repo, a short case study, and the demo. Explain the
biological question, the hardest methodological decision, the external result, and
what remains unresolved. Disclose AI assistance and describe only work you can defend.
Do not claim a new algorithm, hospital patient analysis or clinical diagnostic tool.

GitHub instructions: https://docs.github.com/en/migrations/importing-source-code/using-the-command-line-to-import-source-code/adding-locally-hosted-code-to-github
Zenodo instructions: https://help.zenodo.org/docs/get-started/quickstart/

Public code/results sharing is distinct from journal publication. PLOS explicitly
allows preprints, but check the exact target journal's current prior-publication rules
before posting a manuscript. Never submit the same manuscript simultaneously to journals.
