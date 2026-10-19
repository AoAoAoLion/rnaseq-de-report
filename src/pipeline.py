#!/usr/bin/env python3
"""Bulk RNA-seq differential expression pipeline (standard delivery).

Input:  a gene x sample count matrix + a sample sheet with the condition column.
Output: QC tables, PCA / volcano / MA / heatmap figures, DESeq2 results and a
        self-contained HTML report.

Usage:
    python src/pipeline.py --counts data/example/counts.csv \
        --metadata data/example/metadata.csv --outdir results
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from jinja2 import Template
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sns.set_theme(style="whitegrid", context="notebook")
PALETTE = {"control": "#4C72B0", "treated": "#C44E52"}


def log(msg: str) -> None:
    print(f"[pipeline] {msg}", flush=True)


def load(counts_path: Path, meta_path: Path, condition_col: str):
    counts = pd.read_csv(counts_path, index_col=0)
    meta = pd.read_csv(meta_path, index_col=0)
    missing = set(counts.columns) - set(meta.index)
    if missing:
        raise SystemExit(f"samples in counts but not in metadata: {sorted(missing)}")
    meta = meta.loc[counts.columns]
    if condition_col not in meta.columns:
        raise SystemExit(f"column '{condition_col}' not in metadata: {list(meta.columns)}")
    return counts, meta


def qc(counts: pd.DataFrame, meta: pd.DataFrame, condition_col: str, outdir: Path) -> dict:
    lib = counts.sum()
    detected = (counts > 0).sum()
    keep = counts.sum(axis=1) >= 10
    filtered = counts.loc[keep]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    order = lib.sort_values().index
    colors = [PALETTE.get(c, "#777777") for c in meta.loc[order, condition_col]]
    axes[0].bar(range(len(order)), lib[order] / 1e6, color=colors)
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels(order, rotation=45, ha="right", fontsize=9)
    axes[0].set_ylabel("library size (millions of reads)")
    axes[0].set_title("Library size per sample")

    axes[1].bar(range(len(order)), detected[order], color=colors)
    axes[1].set_xticks(range(len(order)))
    axes[1].set_xticklabels(order, rotation=45, ha="right", fontsize=9)
    axes[1].set_ylabel("genes detected (count > 0)")
    axes[1].set_title("Detected genes per sample")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "01_qc_library.png", dpi=150)
    plt.close(fig)

    pd.DataFrame({"library_size": lib, "genes_detected": detected}).to_csv(
        outdir / "tables" / "qc_per_sample.csv"
    )
    log(f"QC: {counts.shape[0]} genes, {counts.shape[1]} samples; "
        f"{counts.shape[0] - filtered.shape[0]} genes filtered (total count < 10)")
    return {"genes_input": int(counts.shape[0]), "genes_after_filter": int(filtered.shape[0])}


def pca_figure(counts: pd.DataFrame, meta: pd.DataFrame, condition_col: str, outdir: Path) -> None:
    logcpm = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1)
    x = logcpm.sub(logcpm.mean(axis=1), axis=0)
    u, s, vt = np.linalg.svd(x.values, full_matrices=False)
    var = (s ** 2) / np.sum(s ** 2) * 100
    pcs = vt[:2].T

    fig, ax = plt.subplots(figsize=(7, 5.5))
    groups = meta[condition_col].astype(str)
    for grp in groups.unique():
        m = (groups == grp).values
        ax.scatter(pcs[m, 0], pcs[m, 1], s=70, label=grp, color=PALETTE.get(grp, None), edgecolor="k")
    # label a few representative samples per group instead of all of them
    dist = np.hypot(pcs[:, 0], pcs[:, 1])
    for grp in groups.unique():
        idx = np.where((groups == grp).values)[0]
        for i in idx[np.argsort(dist[idx])[-2:]]:
            ax.annotate(counts.columns[i], (pcs[i, 0], pcs[i, 1]), fontsize=8,
                        xytext=(5, 5), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% variance)")
    ax.set_title("PCA on log2 CPM")
    ax.legend(frameon=False, fontsize=9, loc="lower right")
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "02_pca.png", dpi=150)
    plt.close(fig)


def differential_expression(counts, meta, condition_col, ref_level, test_level, outdir: Path, cpus: int):
    dds = DeseqDataSet(
        counts=counts.T,          # pydeseq2 expects samples x genes
        metadata=meta,
        design=f"~{condition_col}",
        refit_cooks=True,
        n_cpus=cpus,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=[condition_col, test_level, ref_level],
        n_cpus=cpus,
        quiet=True,
    )
    stats.summary()
    res = stats.results_df.copy()
    res.index.name = "gene"
    res = res.sort_values("padj", na_position="last")
    res.to_csv(outdir / "tables" / "deseq2_results.csv")
    return dds, res


def plots(res: pd.DataFrame, dds, counts, meta, condition_col, ref_level, test_level, outdir: Path):
    sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() >= 1)]

    fig, ax = plt.subplots(figsize=(9, 7))
    ns = res[(res["padj"] >= 0.05) | (res["log2FoldChange"].abs() < 1)]
    ax.scatter(ns["log2FoldChange"], -np.log10(ns["padj"]), s=6, c="#B0B7C3", label="not significant")
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]
    ax.scatter(up["log2FoldChange"], -np.log10(up["padj"]), s=8, c="#C44E52", label=f"up ({len(up)})")
    ax.scatter(down["log2FoldChange"], -np.log10(down["padj"]), s=8, c="#4C72B0", label=f"down ({len(down)})")
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=1)
    ax.axvline(1, ls="--", c="k", lw=1)
    ax.axvline(-1, ls="--", c="k", lw=1)
    ax.set_xlabel(f"log2 fold change ({test_level} / {ref_level})")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Volcano plot")
    ax.legend(frameon=False, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "03_volcano.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(res["baseMean"], res["log2FoldChange"], s=6, c="#B0B7C3")
    ax.scatter(sig["baseMean"], sig["log2FoldChange"], s=8, c="#C44E52")
    ax.set_xscale("log")
    ax.set_xlabel("mean normalised count")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot")
    ax.axhline(0, ls="--", c="k", lw=1)
    fig.tight_layout()
    fig.savefig(outdir / "figures" / "04_ma.png", dpi=150)
    plt.close(fig)

    top = sig.reindex(sig["padj"].sort_values().index).head(40).index
    if len(top):
        norm = dds.layers["normed_counts"].T if "normed_counts" in dds.layers else counts
        mat = pd.DataFrame(norm, index=counts.index, columns=counts.columns).loc[top]
        mat = np.log2(mat.div(mat.sum(axis=0), axis=1) * 1e6 + 1)
        mat = mat.sub(mat.mean(axis=1), axis=0)
        cg = sns.clustermap(mat, cmap="RdBu_r", center=0, figsize=(9, 11),
                            col_colors=meta[condition_col].map(PALETTE).values,
                            yticklabels=True, xticklabels=True, cbar_pos=(0.02, 0.83, 0.03, 0.12))
        cg.fig.suptitle(f"Top 40 DE genes (padj < 0.05, |log2FC| >= 1)", y=1.02)
        cg.fig.savefig(outdir / "figures" / "05_heatmap_top_de.png", dpi=150, bbox_inches="tight")
        plt.close(cg.fig)
    return sig


def enrichment(sig: pd.DataFrame, outdir: Path) -> pd.DataFrame | None:
    """Gene-set enrichment through Enrichr. Silently skipped without internet."""
    try:
        import gseapy as gp

        genes = sig.index.tolist()[:2000]
        if len(genes) < 10:
            return None
        enr = gp.enrichr(gene_list=genes, gene_sets=["KEGG_2021_Human", "GO_Biological_Process_2023"],
                         outdir=None, cutoff=1.0)
        df = enr.results.sort_values("P-value").head(30)
        if df.empty:
            log("enrichment skipped: no gene-set overlap (input gene IDs are not standard symbols)")
            return None
        df.to_csv(outdir / "tables" / "enrichment.csv", index=False)
        fig, ax = plt.subplots(figsize=(10, 7))
        top = df.head(15).iloc[::-1]
        ax.barh(top["Term"].str.slice(0, 60), -np.log10(top["P-value"]), color="#55A868")
        ax.set_xlabel("-log10 p-value")
        ax.set_title("Top enriched terms")
        fig.tight_layout()
        fig.savefig(outdir / "figures" / "06_enrichment.png", dpi=150)
        plt.close(fig)
        log(f"enrichment: {len(df)} terms written")
        return df
    except Exception as exc:  # noqa: BLE001
        log(f"enrichment skipped ({type(exc).__name__}: {str(exc)[:80]})")
        return None


REPORT = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Bulk RNA-seq differential expression report</title>
<style>
 body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;max-width:980px;margin:40px auto;padding:0 20px;color:#1b2430;line-height:1.6}
 h1{border-bottom:3px solid #2ec4b6;padding-bottom:10px}
 h2{margin-top:38px;border-bottom:1px solid #e3e8ef;padding-bottom:6px}
 table{border-collapse:collapse;width:100%;font-size:14px}
 th,td{border:1px solid #e3e8ef;padding:8px 10px;text-align:left}
 th{background:#f6f8fb}
 img{max-width:100%;border:1px solid #e3e8ef;border-radius:6px;margin:10px 0}
 code{background:#f2f4f7;padding:2px 5px;border-radius:3px;font-size:13px}
 .note{background:#fff8e6;border-left:4px solid #f0b429;padding:12px 16px;margin:16px 0}
 .meta{color:#5b6b7c;font-size:14px}
</style></head><body>
<h1>Bulk RNA-seq differential expression report</h1>
<p class="meta">Contrast: <code>{{ test }}</code> vs <code>{{ ref }}</code> &middot; {{ n_samples }} samples &middot;
{{ genes_after_filter }} genes after filtering &middot; generated {{ date }}</p>

<h2>1. Summary</h2>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Samples analysed</td><td>{{ n_samples }}</td></tr>
<tr><td>Genes in input matrix</td><td>{{ genes_input }}</td></tr>
<tr><td>Genes after low-count filter (total &lt; 10)</td><td>{{ genes_after_filter }}</td></tr>
<tr><td>Genes with padj &lt; 0.05</td><td>{{ n_padj }}</td></tr>
<tr><td>Genes with padj &lt; 0.05 and |log2FC| &ge; 1</td><td>{{ n_sig }}</td></tr>
<tr><td>Up / down</td><td>{{ n_up }} / {{ n_down }}</td></tr>
{% if recall %}<tr><td>Recall against truth set (simulated data)</td><td>{{ recall }}</td></tr>{% endif %}
</table>

{% if recall_note %}<div class="note">{{ recall_note }}</div>{% endif %}

<h2>2. Sample quality</h2>
<img src="figures/01_qc_library.png" alt="library size and detected genes">
<img src="figures/02_pca.png" alt="PCA">
<p>Library sizes and detected-gene counts are reported per sample. Samples were not excluded;
any sample that deviates strongly is flagged here rather than silently dropped.</p>

<h2>3. Differential expression</h2>
<p>Model: <code>~{{ condition_col }}</code>, DESeq2 with median-of-ratios normalisation
(implementation: PyDESeq2 {{ pydeseq2_version }}). Genes are called at padj &lt; 0.05
(Benjamini-Hochberg) and |log2 fold change| &ge; 1 for the figures.</p>
<img src="figures/03_volcano.png" alt="volcano">
<img src="figures/04_ma.png" alt="MA plot">
{% if heatmap %}<img src="figures/05_heatmap_top_de.png" alt="heatmap">{% endif %}

<h2>4. Enrichment</h2>
{% if enrichment %}
<img src="figures/06_enrichment.png" alt="enrichment">
<table><tr><th>Term</th><th>Overlap</th><th>P-value</th><th>Adjusted P</th></tr>
{% for r in enrichment %}<tr><td>{{ r.Term }}</td><td>{{ r.Overlap }}</td><td>{{ '%.2e'|format(r['P-value']) }}</td><td>{{ '%.2e'|format(r['Adjusted P-value']) }}</td></tr>{% endfor %}
</table>
{% else %}<p>Enrichment was not run for this report (no network access or too few genes).</p>{% endif %}

<h2>5. Files</h2>
<ul>
<li><code>tables/qc_per_sample.csv</code> &mdash; library size and detected genes per sample</li>
<li><code>tables/deseq2_results.csv</code> &mdash; full DESeq2 result table</li>
<li><code>tables/enrichment.csv</code> &mdash; enrichment terms (if run)</li>
<li><code>src/pipeline.py</code> &mdash; the script that produced this report</li>
</ul>

<h2>6. Methods text (adapt for your manuscript)</h2>
<p>Read counts were filtered to genes with a total count of at least 10 across samples.
Differential expression was tested with DESeq2 using a median-of-ratios normalisation and
the Wald test for the contrast {{ test }} versus {{ ref }}. P-values were adjusted for
multiple testing with the Benjamini-Hochberg procedure; genes with adjusted p &lt; 0.05 were
considered differentially expressed. Figures show the full result set; the heatmap displays
the top 40 genes ranked by adjusted p-value among those passing both the adjusted p-value and
the |log2 fold change| &ge; 1 threshold.</p>
</body></html>
"""


def write_report(outdir: Path, ctx: dict, enrichment_df) -> None:
    enrichment_rows = enrichment_df.to_dict("records") if enrichment_df is not None else None
    html = Template(REPORT).render(**ctx, enrichment=enrichment_rows)
    (outdir / "report.html").write_text(html)
    log(f"report written: {outdir/'report.html'}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", default="data/example/counts.csv")
    ap.add_argument("--metadata", default="data/example/metadata.csv")
    ap.add_argument("--condition-col", default="condition")
    ap.add_argument("--ref", default="control")
    ap.add_argument("--test", default="treated")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--cpus", type=int, default=2)
    ap.add_argument("--truth", default=None, help="optional truth.csv for a recall check")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    for stale in (outdir / "figures").glob("*.png"):
        stale.unlink()
    (outdir / "tables").mkdir(parents=True, exist_ok=True)

    counts, meta = load(Path(args.counts), Path(args.metadata), args.condition_col)
    info = qc(counts, meta, args.condition_col, outdir)

    filtered = counts.loc[counts.sum(axis=1) >= 10]
    pca_figure(filtered, meta, args.condition_col, outdir)
    dds, res = differential_expression(filtered, meta, args.condition_col, args.ref, args.test,
                                       outdir, args.cpus)
    sig = plots(res, dds, filtered, meta, args.condition_col, args.ref, args.test, outdir)
    enr = enrichment(sig, outdir)

    recall = None
    recall_note = None
    if args.truth and Path(args.truth).exists():
        truth = pd.read_csv(args.truth).set_index("gene")
        called = set(res.index[res["padj"] < 0.05])
        true_de = set(truth.index[truth["is_de"]])
        hit = len(called & true_de)
        recall = f"{hit}/{len(true_de)} ({100*hit/len(true_de):.1f}%)"
        recall_note = (f"This is simulated data with a known truth set: {len(true_de)} genes were "
                       f"spiked in. {recall} were recovered at padj&lt;0.05, with "
                       f"{len(called - true_de)} additional calls outside the truth set. "
                       f"Applying the additional |log2FC|&ge;1 filter narrows the list to {len(sig)} genes.")

    import pydeseq2

    ctx = dict(
        test=args.test, ref=args.ref, n_samples=int(counts.shape[1]),
        genes_input=info["genes_input"], genes_after_filter=info["genes_after_filter"],
        n_padj=int((res["padj"] < 0.05).sum()), n_sig=int(len(sig)),
        n_up=int((sig["log2FoldChange"] > 0).sum()), n_down=int((sig["log2FoldChange"] < 0).sum()),
        condition_col=args.condition_col, pydeseq2_version=pydeseq2.__version__,
        date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"), recall=recall,
        recall_note=recall_note, heatmap=len(sig) > 0,
    )
    write_report(outdir, ctx, enr)

    summary = {k: ctx[k] for k in
               ["test", "ref", "n_samples", "genes_input", "genes_after_filter",
                "n_padj", "n_sig", "n_up", "n_down", "recall"]}
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"done: {ctx['n_padj']} genes padj<0.05, {ctx['n_sig']} passing fold-change filter")


if __name__ == "__main__":
    main()
