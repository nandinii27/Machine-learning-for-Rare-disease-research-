"""
Value-of-information analysis for N-of-1 variant interpretation.

Two decision problems are modelled:

  1. DIAGNOSTIC  - which causal configuration explains the proband, given the
     candidate set recovered from the provided call set. Computed from the
     pre-submission evidence state only; leaderboard feedback is deliberately
     not conditioned on, because the purpose is to show what the pipeline
     would recommend in an undiagnosed case where no oracle exists.

  2. MECHANISTIC - how much residual BubR1 function the missense allele
     retains. This remains open after the diagnostic question closes, and it
     is the assumption any therapeutic hypothesis rests on.

For each candidate experiment e with discrete outcomes r, the expected
reduction in Shannon entropy over the hypothesis set C is

    VoI(e) = H(C) - sum_r P(r) * H(C | e = r)

with P(r) = sum_h P(h) P(r | h) and P(h | e = r) by Bayes. VoI is reported in
bits, and also normalised by cost, since the ranking that matters clinically
is information per unit of cost and delay, not information alone.

All priors and likelihoods below are explicit elicited subjective estimates,
stated so they can be challenged and revised. They are not fitted to data and
are not claimed to be calibrated. The rationale for each is given inline.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log2


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Experiment:
    """A candidate experiment.

    likelihood[hypothesis][outcome] = P(outcome | hypothesis).
    Each hypothesis row must sum to 1.
    """
    name: str
    description: str
    outcomes: tuple[str, ...]
    likelihood: dict[str, tuple[float, ...]]
    turnaround_days: float
    relative_cost: float

    def validate(self, hypotheses: dict[str, float]) -> None:
        for h in hypotheses:
            if h not in self.likelihood:
                raise ValueError(f"{self.name}: no likelihood row for {h}")
            row = self.likelihood[h]
            if len(row) != len(self.outcomes):
                raise ValueError(f"{self.name}/{h}: row length != outcome count")
            if abs(sum(row) - 1.0) > 1e-9:
                raise ValueError(f"{self.name}/{h}: row sums to {sum(row)}, not 1")


def entropy(dist: dict[str, float]) -> float:
    return -sum(p * log2(p) for p in dist.values() if p > 0)


def outcome_marginal(prior: dict[str, float], exp: Experiment) -> tuple[float, ...]:
    return tuple(
        sum(prior[h] * exp.likelihood[h][i] for h in prior)
        for i in range(len(exp.outcomes))
    )


def posterior(prior: dict[str, float], exp: Experiment, outcome_index: int) -> dict[str, float]:
    joint = {h: prior[h] * exp.likelihood[h][outcome_index] for h in prior}
    z = sum(joint.values())
    if z == 0:
        raise ValueError(f"{exp.name}: outcome {exp.outcomes[outcome_index]} has zero marginal")
    return {h: v / z for h, v in joint.items()}


def expected_posterior_entropy(prior: dict[str, float], exp: Experiment) -> float:
    marg = outcome_marginal(prior, exp)
    return sum(
        p_r * entropy(posterior(prior, exp, i))
        for i, p_r in enumerate(marg)
        if p_r > 0
    )


def value_of_information(prior: dict[str, float], exp: Experiment) -> float:
    return entropy(prior) - expected_posterior_entropy(prior, exp)


def rank_experiments(
    prior: dict[str, float],
    experiments: list[Experiment],
) -> list[dict[str, float | str]]:
    h0 = entropy(prior)
    rows = []
    for e in experiments:
        e.validate(prior)
        voi = value_of_information(prior, e)
        rows.append({
            "experiment": e.name,
            "voi_bits": voi,
            "frac_resolved": voi / h0 if h0 > 0 else 0.0,
            "days": e.turnaround_days,
            "cost": e.relative_cost,
            "bits_per_cost": voi / e.relative_cost,
            "bits_per_day": voi / e.turnaround_days,
        })
    return sorted(rows, key=lambda r: -r["voi_bits"])


def decision_boundary(
    prior: dict[str, float],
    exp: Experiment,
    hypothesis: str,
) -> list[tuple[str, float, float]]:
    """Posterior on one hypothesis under each outcome.

    This is the 'what would change my mind' output: it names the specific
    result that would move a hypothesis from favoured to disfavoured.
    """
    marg = outcome_marginal(prior, exp)
    return [
        (exp.outcomes[i], marg[i], posterior(prior, exp, i)[hypothesis])
        for i in range(len(exp.outcomes))
    ]


# ---------------------------------------------------------------------------
# Problem 1: diagnostic
# ---------------------------------------------------------------------------
#
# Priors reflect the pre-submission evidence state: a PASS, well-covered
# nonsense allele at gnomAD AF 1e-4 in exon 17/23, plus a gnomAD-absent
# missense in exon 23/23 with concordant damaging predictions, against a
# three-gene prior in which CEP57 carried no rare variation and TRIP13's two
# gene-body candidates were shown to be in cis by read-backed phasing.
#
# H3 is given substantial weight because second alleles in BUB1B compound
# heterozygotes are frequently non-coding or structural, and the analysis used
# the provided call set without realignment.

DIAGNOSTIC_PRIOR: dict[str, float] = {
    "H1_nonsense_plus_missense": 0.55,
    "H2_nonsense_plus_deep_intronic": 0.15,
    "H3_nonsense_plus_uncalled": 0.20,
    "H4_cis_or_non_BUB1B": 0.10,
}

PARENTAL_SEQUENCING = Experiment(
    name="parental_sequencing",
    description="Targeted sequencing of both parents at the two candidate loci.",
    outcomes=("trans_confirmed", "cis_detected", "uninformative"),
    # Under H1 the two alleles are genuinely in trans. Under H2 and H3 the
    # missense is incidental, so its phase relative to the nonsense allele is
    # close to a coin flip - which is why this experiment separates H1 from H4
    # sharply but from H2/H3 only weakly.
    likelihood={
        "H1_nonsense_plus_missense": (0.93, 0.02, 0.05),
        "H2_nonsense_plus_deep_intronic": (0.48, 0.47, 0.05),
        "H3_nonsense_plus_uncalled": (0.48, 0.47, 0.05),
        "H4_cis_or_non_BUB1B": (0.10, 0.85, 0.05),
    },
    turnaround_days=21,
    relative_cost=1.0,
)

RNA_SEQ = Experiment(
    name="rna_seq_fibroblast",
    description="RNA-seq of patient fibroblasts: pseudoexon activation, NMD, allelic expression.",
    outcomes=("aberrant_splicing", "normal_splicing", "inconclusive"),
    likelihood={
        "H1_nonsense_plus_missense": (0.08, 0.82, 0.10),
        "H2_nonsense_plus_deep_intronic": (0.85, 0.05, 0.10),
        "H3_nonsense_plus_uncalled": (0.35, 0.55, 0.10),
        "H4_cis_or_non_BUB1B": (0.15, 0.75, 0.10),
    },
    turnaround_days=28,
    relative_cost=2.0,
)

LONG_READ_WGS = Experiment(
    name="long_read_wgs",
    description="Long-read WGS: structural and mobile-element variants, plus phase across 10.9 kb.",
    outcomes=("sv_or_mei_found", "no_sv_phase_trans", "no_sv_phase_cis"),
    likelihood={
        "H1_nonsense_plus_missense": (0.05, 0.90, 0.05),
        "H2_nonsense_plus_deep_intronic": (0.07, 0.46, 0.47),
        "H3_nonsense_plus_uncalled": (0.80, 0.12, 0.08),
        "H4_cis_or_non_BUB1B": (0.10, 0.12, 0.78),
    },
    turnaround_days=35,
    relative_cost=4.0,
)

CHECKPOINT_ASSAY = Experiment(
    name="checkpoint_functional_assay",
    description="Mitotic checkpoint competence in patient cells after microtubule perturbation.",
    # Confirms that the checkpoint is defective, but does not localise which
    # allele is responsible - the near-identical rows for H1 to H3 are the
    # model correctly reporting that this experiment answers a different
    # question from the one currently in doubt.
    outcomes=("severe_defect", "moderate_defect", "near_normal"),
    likelihood={
        "H1_nonsense_plus_missense": (0.55, 0.40, 0.05),
        "H2_nonsense_plus_deep_intronic": (0.55, 0.40, 0.05),
        "H3_nonsense_plus_uncalled": (0.55, 0.40, 0.05),
        "H4_cis_or_non_BUB1B": (0.20, 0.35, 0.45),
    },
    turnaround_days=60,
    relative_cost=6.0,
)

DIAGNOSTIC_EXPERIMENTS = [
    PARENTAL_SEQUENCING,
    RNA_SEQ,
    LONG_READ_WGS,
    CHECKPOINT_ASSAY,
]


# ---------------------------------------------------------------------------
# Problem 2: mechanistic
# ---------------------------------------------------------------------------
#
# BubR1 is a pseudokinase: the C-terminal domain has no catalytic activity and
# functions as a scaffold, recruiting PP2A-B56 to the kinetochore for
# checkpoint silencing and chromosome alignment. p.Asn1002Lys therefore
# degrades an interaction surface rather than an active site, and the defect
# is expected to be graded rather than binary.
#
# With one null allele, the severity of the phenotype is set by how much
# function the missense allele retains. That quantity is unmeasured, and every
# therapeutic hypothesis is conditional on it.

MECHANISTIC_PRIOR: dict[str, float] = {
    "M1_severe_hypomorph": 0.40,   # under ~20 percent residual scaffold function
    "M2_moderate_hypomorph": 0.40,  # roughly 20 to 60 percent
    "M3_mild": 0.20,                # over 60 percent; phenotype driven mainly by the null allele
}

IN_SILICO_DDG = Experiment(
    name="in_silico_stability",
    description="Structure-based stability and interface-disruption prediction for p.Asn1002Lys.",
    outcomes=("strongly_destabilising", "mildly_destabilising", "neutral"),
    likelihood={
        "M1_severe_hypomorph": (0.70, 0.25, 0.05),
        "M2_moderate_hypomorph": (0.30, 0.55, 0.15),
        "M3_mild": (0.10, 0.35, 0.55),
    },
    turnaround_days=0.1,
    relative_cost=0.01,
)

PP2A_COIP = Experiment(
    name="pp2a_b56_binding",
    description="Co-immunoprecipitation quantifying PP2A-B56 recruitment by the variant protein.",
    outcomes=("low_binding", "partial_binding", "normal_binding"),
    likelihood={
        "M1_severe_hypomorph": (0.85, 0.12, 0.03),
        "M2_moderate_hypomorph": (0.20, 0.70, 0.10),
        "M3_mild": (0.05, 0.25, 0.70),
    },
    turnaround_days=14,
    relative_cost=3.0,
)

MITOTIC_TIMING = Experiment(
    name="mitotic_timing",
    description="Live-cell imaging of nuclear envelope breakdown to anaphase onset in patient cells.",
    outcomes=("shortened", "intermediate", "normal"),
    likelihood={
        "M1_severe_hypomorph": (0.80, 0.17, 0.03),
        "M2_moderate_hypomorph": (0.30, 0.60, 0.10),
        "M3_mild": (0.08, 0.32, 0.60),
    },
    turnaround_days=21,
    relative_cost=4.0,
)

ANEUPLOIDY_RATE = Experiment(
    name="aneuploidy_rate",
    description="Karyotyping or single-cell sequencing to quantify the fraction of aneuploid cells.",
    outcomes=("high", "moderate", "low"),
    likelihood={
        "M1_severe_hypomorph": (0.70, 0.25, 0.05),
        "M2_moderate_hypomorph": (0.40, 0.50, 0.10),
        "M3_mild": (0.15, 0.45, 0.40),
    },
    turnaround_days=10,
    relative_cost=2.0,
)

MECHANISTIC_EXPERIMENTS = [
    IN_SILICO_DDG,
    PP2A_COIP,
    MITOTIC_TIMING,
    ANEUPLOIDY_RATE,
]


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_ranking(title: str, prior: dict[str, float], experiments: list[Experiment]) -> None:
    print(f"\n{title}")
    print("=" * len(title))
    print(f"Prior entropy: {entropy(prior):.3f} bits over {len(prior)} hypotheses\n")
    for h, p in sorted(prior.items(), key=lambda kv: -kv[1]):
        print(f"  {p:5.2f}  {h}")
    print()
    header = f"{'experiment':<30}{'bits':>8}{'frac':>8}{'days':>8}{'cost':>8}{'bits/cost':>12}"
    print(header)
    print("-" * len(header))
    for r in rank_experiments(prior, experiments):
        print(
            f"{r['experiment']:<30}{r['voi_bits']:>8.3f}{r['frac_resolved']:>8.2f}"
            f"{r['days']:>8.1f}{r['cost']:>8.2f}{r['bits_per_cost']:>12.2f}"
        )


def print_decision_boundary(
    prior: dict[str, float],
    exp: Experiment,
    hypothesis: str,
) -> None:
    print(f"\nWhat would change the verdict on {hypothesis}")
    print(f"via {exp.name} (prior {prior[hypothesis]:.2f})")
    print("-" * 60)
    for outcome, p_outcome, post in decision_boundary(prior, exp, hypothesis):
        direction = "up" if post > prior[hypothesis] else "down"
        print(f"  {outcome:<24} p={p_outcome:.2f}  ->  posterior {post:.2f}  ({direction})")


if __name__ == "__main__":
    print_ranking(
        "Diagnostic decision problem (pre-submission evidence state)",
        DIAGNOSTIC_PRIOR,
        DIAGNOSTIC_EXPERIMENTS,
    )
    print_decision_boundary(DIAGNOSTIC_PRIOR, RNA_SEQ, "H1_nonsense_plus_missense")
    print_decision_boundary(DIAGNOSTIC_PRIOR, LONG_READ_WGS, "H3_nonsense_plus_uncalled")

    print_ranking(
        "Mechanistic decision problem (open)",
        MECHANISTIC_PRIOR,
        MECHANISTIC_EXPERIMENTS,
    )
    print_decision_boundary(MECHANISTIC_PRIOR, PP2A_COIP, "M1_severe_hypomorph")
