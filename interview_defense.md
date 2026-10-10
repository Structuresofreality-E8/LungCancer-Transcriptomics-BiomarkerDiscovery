# Questions to be able to answer

1. Why not feed RPKM to PyDESeq2? It models count sampling with sample size factors;
   length/library-normalized continuous expression is not the original count process.
2. Why group by patient? Related tumor/normal specimens share patient-specific biology;
   putting them on opposite sides of a split can produce optimistic evaluation.
3. Why nested CV? Inner folds choose features/hyperparameters; outer patients estimate
   the performance of that entire procedure. Global supervised selection leaks labels.
4. Why ranks? They reduce per-sample scale differences between length-normalized inputs;
   they neither remove every study bias nor support inference from only 25 assayed genes.
5. Why can a replicated gene be unstable? Correlated genes can substitute for one another
   under elastic-net selection. Paired association tests ask a different question from
   whether that exact predictor is repeatedly selected in a multivariable model.
6. Why isn't 0.9927 AUC 99.27% diagnostic accuracy? AUC measures ranking across thresholds;
   performance at threshold 0.5 is summarized by the separate confusion matrix.
7. Why isn't 100% specificity certain? Only 77 external negatives were observed; report an
   interval that accounts for finite sample size as well as the conditional bootstrap.
8. What is the QC-only baseline for? To reveal whether library summary information predicts
   labels. It is a warning diagnostic, not proof that the expression model is confounded
   or proof of its absence when baseline discrimination is low.
9. What was AI-assisted? Initial implementation and packaging. Describe your actual
   verification, biological interpretation, debugging and project decisions honestly.
10. What would falsify your claim? Failure to transfer to an untouched suitable cohort,
    a dominant technical explanation, or lack of specificity to the claimed tissue task.
