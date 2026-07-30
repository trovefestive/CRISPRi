# Dolcetto library representation QC

Run the analysis from the project root on the external SSD:

```bash
Rscript analysis/run_dolcetto_qc.R --root /Volumes/Subash/CRISPRi --out-dir analysis/results
```

The script reads only R1 for guide counting, streams the compressed FASTQs in chunks, preserves exact matches as the authoritative counts, and reports one-substitution matches separately as diagnostics.

For a limited validation run:

```bash
Rscript analysis/run_dolcetto_qc.R --root /Volumes/Subash/CRISPRi --out-dir analysis/test_results --max-reads 100000 --skip-md5
```
