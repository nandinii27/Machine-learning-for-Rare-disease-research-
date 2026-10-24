"""
Estimate mosaic aneuploidy burden per chromosome from a single-sample VCF.

Two orthogonal signals are computed from the same call set:

  BAF deviation
      At a heterozygous SNP in a purely diploid sample, the alternate allele
      fraction is 0.5. If a fraction f of cells carry a trisomy for that
      chromosome, het sites split into branches at (1+f)/(2+f) and 1/(2+f);
      for monosomy, (1-f)/(2-f) and 1/(2-f). In both cases the absolute
      deviation from 0.5 is f / (2 * (2 +/- f)), which is approximately f/4
      for small f. The deviation is therefore an estimator of mosaic fraction
      but is symmetric: it cannot distinguish gain from loss.

  Relative read depth
      A chromosome at average copy number 2+f carries relative coverage
      1 + f/2; at 2-f, coverage 1 - f/2. Noisier than BAF because of GC and
      mappability bias, but signed, so it resolves the direction that BAF
      leaves ambiguous.

Interpretation requires both. A chromosome with elevated BAF deviation and
elevated depth indicates mosaic gain; elevated deviation with depressed depth
indicates mosaic loss; elevated deviation with flat depth is more consistent
with copy-neutral loss of heterozygosity or with sample contamination than
with aneuploidy.

This is a screen, not a diagnostic. Established tools for mosaic chromosomal
alterations (MoChA and comparable callers) model phase, segment within
chromosomes, and fit explicit mixture likelihoods. This script computes
chromosome-level summary statistics only, and will miss segmental and
arm-level events entirely.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from statistics import median

import pandas as pd


AUTOSOMES = [str(i) for i in range(1, 23)]

MIN_DP = 15
MAX_DP = 100
MIN_GQ = 30
VAF_WINDOW = (0.15, 0.85)


@dataclass
class ChromAccumulator:
    vafs: list[float] = field(default_factory=list)
    depths: list[int] = field(default_factory=list)


def _format_index(fmt: str, key: str) -> int | None:
    parts = fmt.split(":")
    return parts.index(key) if key in parts else None


def collect_het_sites(vcf_path: str, max_records: int | None = None) -> dict[str, ChromAccumulator]:
    """Stream the VCF and accumulate per-chromosome het-SNP VAFs and depths.

    Filters applied, each with a reason:
      PASS only          - non-PASS calls have unreliable allele depths
      biallelic SNV      - indels and multiallelics distort AD
      GT heterozygous    - BAF is only defined at het sites
      DP and GQ bounds   - low depth inflates VAF variance; very high depth
                           marks collapsed repeats where AD is meaningless
      VAF window         - excludes likely homozygous miscalls and somatic
                           or artefactual low-fraction calls
    """
    acc: dict[str, ChromAccumulator] = {c: ChromAccumulator() for c in AUTOSOMES}
    seen = 0

    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue

            chrom = cols[0].removeprefix("chr")
            if chrom not in acc:
                continue

            if cols[6] != "PASS":
                continue

            ref, alt = cols[3], cols[4]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue

            fmt, sample = cols[8], cols[9]
            gt_i = _format_index(fmt, "GT")
            ad_i = _format_index(fmt, "AD")
            dp_i = _format_index(fmt, "DP")
            gq_i = _format_index(fmt, "GQ")
            if None in (gt_i, ad_i, dp_i, gq_i):
                continue

            vals = sample.split(":")
            if len(vals) <= max(gt_i, ad_i, dp_i, gq_i):
                continue

            gt = vals[gt_i].replace("|", "/")
            if gt not in ("0/1", "1/0"):
                continue

            try:
                dp = int(vals[dp_i])
                gq = int(vals[gq_i])
                ad_ref, ad_alt = (int(x) for x in vals[ad_i].split(",")[:2])
            except ValueError:
                continue

            if not (MIN_DP <= dp <= MAX_DP) or gq < MIN_GQ:
                continue

            total = ad_ref + ad_alt
            if total == 0:
                continue

            vaf = ad_alt / total
            if not (VAF_WINDOW[0] <= vaf <= VAF_WINDOW[1]):
                continue

            acc[chrom].vafs.append(vaf)
            acc[chrom].depths.append(dp)

            seen += 1
            if max_records is not None and seen >= max_records:
                break

    return acc


def summarise(acc: dict[str, ChromAccumulator], min_sites: int = 2000) -> pd.DataFrame:
    """Per-chromosome BAF deviation and relative depth, with mosaic estimates."""
    rows = []
    for chrom in AUTOSOMES:
        a = acc[chrom]
        if len(a.vafs) < min_sites:
            continue
        devs = [abs(v - 0.5) for v in a.vafs]
        rows.append({
            "chrom": chrom,
            "n_het_sites": len(a.vafs),
            "mean_abs_baf_dev": sum(devs) / len(devs),
            "median_depth": median(a.depths),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    # Genome-wide baselines. The median across chromosomes is used rather than
    # the mean so that a genuinely affected chromosome does not inflate its own
    # reference point.
    baseline_dev = df["mean_abs_baf_dev"].median()
    baseline_depth = df["median_depth"].median()

    df["excess_baf_dev"] = df["mean_abs_baf_dev"] - baseline_dev
    df["relative_depth"] = df["median_depth"] / baseline_depth

    # d ~= f/4 for small f, so f ~= 4d. Valid only where the excess is real;
    # negative excess is reported as zero rather than a negative fraction.
    df["est_mosaic_fraction_baf"] = (4 * df["excess_baf_dev"]).clip(lower=0.0)

    # Depth-based estimate, signed: relative coverage 1 + f/2 for gain.
    df["est_mosaic_fraction_depth"] = 2 * (df["relative_depth"] - 1.0)

    def call(row):
        if row["est_mosaic_fraction_baf"] < 0.02:
            return "no evidence"
        if row["est_mosaic_fraction_depth"] > 0.02:
            return "possible mosaic gain"
        if row["est_mosaic_fraction_depth"] < -0.02:
            return "possible mosaic loss"
        return "BAF shift, depth flat - check LOH or contamination"

    df["interpretation"] = df.apply(call, axis=1)
    return df.sort_values("est_mosaic_fraction_baf", ascending=False)


def plot_baf_distributions(acc: dict[str, ChromAccumulator], chroms: list[str]):
    """Overlaid VAF histograms. A mosaic chromosome shows a widened or
    bimodal distribution relative to the diploid baseline.

    Returns the figure; does not write to disk.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 5))
    for c in chroms:
        vafs = acc[c].vafs
        if not vafs:
            continue
        ax.hist(vafs, bins=120, histtype="step", density=True, label=f"chr{c}")
    ax.axvline(0.5, linestyle="--", linewidth=0.8)
    ax.set_xlabel("alternate allele fraction at heterozygous SNPs")
    ax.set_ylabel("density")
    ax.legend()
    fig.tight_layout()
    plt.show()
    return fig


if __name__ == "__main__":
    VCF = "./data/WGS_EX2312012_HGWCNDSX7.vcf.gz"

    acc = collect_het_sites(VCF)
    summary = summarise(acc)

    pd.set_option("display.width", 160)
    pd.set_option("display.max_columns", 20)
    print(summary.to_string(index=False))

    top = summary.head(3)["chrom"].tolist()
    baseline = summary.sort_values("est_mosaic_fraction_baf").head(1)["chrom"].tolist()
    plot_baf_distributions(acc, top + baseline)
