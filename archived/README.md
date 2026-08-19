# Archived workflow components

This directory preserves the former master orchestration and annotation stack
for historical reference:

- `prokaryont.sh` and `run_all.sh`
- `06_annotate_assess.sh`
- `07a_predict_genes.sh` through `07e_reconcile_merge.sh`
- their annotation helper programs under `utils/`

These files are no longer part of the active ProkaryONT workflow and their
relative paths have not been adapted for execution from this directory. The
active pipeline is intentionally step-by-step: invoke the numbered scripts in
the repository root separately and review each stage before starting the next.
