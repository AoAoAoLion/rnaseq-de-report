#!/usr/bin/env python3
"""Generate a realistic synthetic bulk RNA-seq count matrix for the demo.

Why synthetic: the demo must run end to end on any machine with no external
downloads. Counts follow a negative-binomial / multinomial model with
library-size variation, gene-specific dispersion and a known set of
differentially expressed genes, so the pipeline can be validated against a
truth set.

Replace data/example/counts.csv and metadata.csv with your own files to run
the same pipeline on real data.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

GENES = 20000
PER_GROUP = 6
LIB_SIZE_MEAN = 24_000_000
LIB_SIZE_SD = 2_000_000
N_DE = 800
SEED = 20261018


def simulate(outdir: Path, genes: int = GENES, per_group: int = PER_GROUP,
             n_de: int = N_DE) -> None:
    rng = np.random.default_rng(SEED)
    outdir.mkdir(parents=True, exist_ok=True)

    samples = [f"ctrl_{i}" for i in range(1, per_group + 1)] + [
        f"trt_{i}" for i in range(1, per_group + 1)
    ]
    condition = ["control"] * per_group + ["treated"] * per_group

    # baseline gene means: log-normal, a handful of very high expressers
    mu = rng.lognormal(mean=5.2, sigma=1.7, size=genes)

    # gene-wise dispersion (NB): small for high-expression genes
    dispersion = 0.02 + 0.3 / np.sqrt(mu + 1)

    # spike in DE genes with a multiplicative effect
    de_idx = rng.choice(genes, size=min(n_de, genes // 4), replace=False)
    log2fc = np.zeros(genes)
    log2fc[de_idx] = np.clip(rng.normal(0, 1.1, size=len(de_idx)), -3.5, 3.5)
    weak = (np.abs(log2fc[de_idx]) < 0.6)
    log2fc[de_idx[weak]] *= 1.8

    lib = rng.normal(LIB_SIZE_MEAN, LIB_SIZE_SD, size=len(samples)).clip(16e6, 34e6)

    counts = np.zeros((genes, len(samples)), dtype=np.int64)
    for j, _ in enumerate(samples):
        base = mu * (2.0 ** log2fc) if condition[j] == "treated" else mu
        # gamma-distributed gene rates -> negative-binomial counts after multinomial sampling,
        # which keeps library size exactly as drawn and avoids composition artefacts
        rate = rng.gamma(shape=1.0 / dispersion, scale=base * dispersion)
        p = rate / rate.sum()
        counts[:, j] = rng.multinomial(int(lib[j]), p)

    gene_ids = [f"GENE{i:05d}" for i in range(1, genes + 1)]
    df_counts = pd.DataFrame(counts, index=gene_ids, columns=samples)
    df_meta = pd.DataFrame({"sample": samples, "condition": condition}).set_index("sample")

    df_counts.to_csv(outdir / "counts.csv")
    df_meta.to_csv(outdir / "metadata.csv")
    pd.DataFrame({"gene": gene_ids, "truth_log2fc": log2fc, "is_de": log2fc != 0}).to_csv(
        outdir / "truth.csv", index=False
    )
    print(f"wrote {outdir/'counts.csv'} {df_counts.shape}, "
          f"{outdir/'metadata.csv'}, {outdir/'truth.csv'} ({len(de_idx)} DE genes)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/example")
    ap.add_argument("--genes", type=int, default=GENES)
    ap.add_argument("--per-group", type=int, default=PER_GROUP)
    ap.add_argument("--n-de", type=int, default=N_DE)
    a = ap.parse_args()
    simulate(Path(a.outdir), genes=a.genes, per_group=a.per_group, n_de=a.n_de)
