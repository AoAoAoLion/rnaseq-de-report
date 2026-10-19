# Bulk RNA-seq differential expression, end to end

A small, reproducible pipeline that turns a gene-by-sample count matrix into a delivered
analysis: QC, PCA, DESeq2 differential expression, publication-style figures and a
self-contained HTML report with a methods paragraph.

Built for standard contract work: one command in, one report out, with the scripts included
so the analysis can be re-run by whoever receives it.

```bash
make setup      # create the virtualenv and install dependencies
make example    # generate the synthetic example dataset
make run        # run the pipeline, writes results/
open results/report.html
```

## What it produces

| Output | Content |
|---|---|
| `results/report.html` | self-contained report: summary table, QC, PCA, volcano, MA, heatmap, methods text |
| `results/figures/` | every figure as a 150 dpi PNG |
| `results/tables/deseq2_results.csv` | full DESeq2 result table (log2FC, SE, stat, p, padj, baseMean) |
| `results/tables/qc_per_sample.csv` | library size and detected genes per sample |
| `results/summary.json` | machine-readable summary of the run |

## Input

Two CSVs.

`counts.csv` — genes in rows, samples in columns, integer counts:

```
gene,ctrl_1,ctrl_2,trt_1,trt_2
GENE00001,1523,1610,44,51
```

`metadata.csv` — one row per sample, indexed by sample name, with a condition column:

```
sample,condition
ctrl_1,control
trt_1,treated
```

If your data starts as FASTQ, quantify first (featureCounts, HTSeq or salmon/tximport) and
feed the resulting matrix here. That step is deliberately outside this repo: reference genome
and annotation version are study-specific decisions that should be documented per project.

## Pipeline steps

1. **Load and validate** — sample names in the matrix and sample sheet must match; the
   pipeline refuses to run on a silent mismatch.
2. **QC** — library size and detected genes per sample, plotted and tabulated. Low-count
   genes are filtered (total count < 10). Samples are never dropped silently.
3. **PCA** — on log2 CPM, centred per gene. This is exploratory only; it is not used for
   the differential expression model.
4. **Differential expression** — DESeq2 (PyDESeq2), median-of-ratios normalisation,
   Benjamini-Hochberg adjustment. The contrast is passed explicitly (`test` vs `ref`),
   so the sign of log2FC is always interpretable.
5. **Figures** — volcano, MA, and a heatmap of the top 40 genes by adjusted p-value among
   those passing padj < 0.05 and |log2FC| >= 1.
6. **Enrichment** — Enrichr (KEGG 2021 Human, GO Biological Process 2023) on the significant
   gene list. Requires network access and real gene symbols; skipped with a clear message
   otherwise.
7. **Report** — HTML rendered from the results, including a methods paragraph that can be
   adapted for a manuscript.

Thresholds (0.05, |log2FC| >= 1) are defaults, visible in the code, not hidden in a config
file three directories away.

## Example data are synthetic, and say so

`data/make_example_data.py` draws counts from a negative-binomial model with a fixed seed
and 800 spiked-in DE genes. The reason is reproducibility: the demo runs offline, in under a
minute, on any machine, and it comes with a truth set, so the pipeline's own behaviour can be
checked: `--truth data/example/truth.csv` reports how many spiked genes were recovered.

Current run on the example data: 800 true DE genes, 585 called at padj < 0.05 (540 of them
true, 45 additional calls), and 331 genes passing the additional fold-change filter. That is
the expected shape for 6 vs 6 samples at ~24M reads: good sensitivity, honest false positives.

Real datasets behave in their own ways (batch effects, outliers, zero inflation). Point
`--counts` and `--metadata` at your own files to run the identical code on real data.

## What this pipeline deliberately does not do

- No alignment or quantification. Those are study-specific and belong upstream.
- No automatic outlier removal. Removing samples changes the answer, so it is a decision
  documented in the report, not a silent preprocessing step.
- No batch correction by default. If a batch covariate exists, add it to the design
  (`--condition-col` is one factor; multi-factor designs are a small change to the model
  formula, not a hidden default).
- No enrichment when the input uses placeholder gene IDs. Fabricating a gene-set result is
  worse than skipping it.

## Reproducibility

- Fixed random seed in the data generator; counts are byte-identical across runs.
- Versions pinned in `requirements.txt`; the report records the PyDESeq2 version used.
- Every figure is regenerated from the saved tables — nothing is hand-edited.

## Tests

```bash
make test
```

The test suite checks the data generator's determinism, the result table's schema, that the
pipeline actually detects the spiked-in signal, and that the report is written with figures.

## Layout

```
data/make_example_data.py   synthetic dataset generator (documented, seeded)
src/pipeline.py             the pipeline and the report template
templates/                  (report template lives in src/pipeline.py to keep it single-file)
results/                    generated output
tests/test_pipeline.py      pytest suite
```

## License

MIT.
