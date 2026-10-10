# Preserve the external evaluation

v0.1.0 represents the already completed model and its first reported GEO evaluation.
Preserve model bytes, signature, fixed gene universe, code, configuration, predictions,
source identifiers and checksums. The exporter creates SHA-256 records after collecting
these artifacts; these are post-run provenance, not a prospective registration.

Do not optimize against the seven missed GEO tumors, drop low-stability genes to improve
GEO scores, adjust the threshold or change normalization and keep the same validation
claim. Any such change creates v0.2.0, with GEO as development evidence and a newly
untouched cohort required for final evaluation.

Descriptive audit of the current seven false negatives is allowed if clearly labeled
post hoc. It cannot retrospectively become prespecified validation. Add a RUN_HISTORY.md
with actual reruns, fixes and which outcomes were inspected; do not invent a clean history.
