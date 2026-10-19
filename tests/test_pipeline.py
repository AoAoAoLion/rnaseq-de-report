"""Tests for the RNA-seq pipeline.

The suite is deliberately end-to-end on a small synthetic dataset: it checks that the
pipeline runs, that the result table has the expected schema, that the spiked-in signal is
actually recovered, and that the report is written with its figures.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
PY = sys.executable


@pytest.fixture(scope="module")
def small_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("run")
    subprocess.run(
        [PY, str(REPO / "data" / "make_example_data.py"), "--outdir", str(tmp / "data"),
         "--genes", "3000", "--per-group", "4", "--n-de", "300"],
        check=True, capture_output=True,
    )
    out = tmp / "results"
    proc = subprocess.run(
        [PY, str(REPO / "src" / "pipeline.py"),
         "--counts", str(tmp / "data" / "counts.csv"),
         "--metadata", str(tmp / "data" / "metadata.csv"),
         "--outdir", str(out),
         "--truth", str(tmp / "data" / "truth.csv"),
         "--cpus", "2"],
        check=True, capture_output=True, text=True,
    )
    return tmp, out, proc.stdout


def test_generator_is_deterministic(tmp_path):
    for d in ("a", "b"):
        subprocess.run([PY, str(REPO / "data" / "make_example_data.py"),
                        "--outdir", str(tmp_path / d), "--genes", "500", "--per-group", "3"],
                       check=True, capture_output=True)
    a = (tmp_path / "a" / "counts.csv").read_text()
    b = (tmp_path / "b" / "counts.csv").read_text()
    assert a == b, "simulated counts must be reproducible across runs"


def test_results_schema(small_run):
    _, out, _ = small_run
    res = pd.read_csv(out / "tables" / "deseq2_results.csv", index_col=0)
    for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
        assert col in res.columns, f"missing column {col}"
    assert len(res) > 1000
    assert res["padj"].notna().sum() > 0


def test_detects_spiked_signal(small_run):
    tmp, out, _ = small_run
    res = pd.read_csv(out / "tables" / "deseq2_results.csv", index_col=0)
    truth = pd.read_csv(tmp / "data" / "truth.csv").set_index("gene")
    called = set(res.index[res["padj"] < 0.05])
    true_de = set(truth.index[truth["is_de"]])
    recall = len(called & true_de) / len(true_de)
    assert recall > 0.4, f"recall {recall:.2f} too low for 300 spiked genes"
    precision = len(called & true_de) / max(len(called), 1)
    assert precision > 0.5, f"precision {precision:.2f} suggests inflated false positives"


def test_report_and_figures(small_run):
    _, out, _ = small_run
    html = (out / "report.html").read_text()
    assert "Bulk RNA-seq differential expression report" in html
    assert "Methods text" in html
    figs = sorted(p.name for p in (out / "figures").glob("*.png"))
    assert {"01_qc_library.png", "02_pca.png", "03_volcano.png", "04_ma.png"} <= set(figs)


def test_summary_written(small_run):
    _, out, stdout = small_run
    import json
    summary = json.loads((out / "summary.json").read_text())
    assert summary["n_samples"] == 8
    assert summary["n_padj"] > 0
    assert "recall" in summary and summary["recall"]
    assert "done" in stdout
