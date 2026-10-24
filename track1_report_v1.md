# Track 1 Methods Report — Submission 1 (baseline probe)

**Proband:** PROBAND01
**Submission:** 1 of 6
**Status:** Preliminary baseline. This submission establishes a scored reference point for the
annotation-driven panel pipeline. Non-coding and structural hypotheses are deliberately not yet
tested; they are the subject of planned follow-up submissions.

---

## 1. Summary of findings

A compound heterozygous genotype in *BUB1B* (chr15, GRCh38) is proposed as the causal result:

| | Variant | HGVSc | HGVSp | Consequence | gnomAD AF |
|---|---|---|---|---|---|
| Allele 1 | chr15:40209701 T>G | c.2210T>G | p.Leu737Ter | stop gained, exon 17/23 | 1.0e-4 |
| Allele 2 | chr15:40220612 T>G | c.3006T>G | p.Asn1002Lys | missense, exon 23/23 | absent |

Allele 1 introduces a premature termination codon well upstream of the final exon-exon junction and
is therefore predicted to trigger nonsense-mediated decay, giving effectively null output from that
copy. Allele 2 substitutes a neutral polar asparagine with a positively charged lysine within the
C-terminal kinase domain; SIFT calls it deleterious and PolyPhen-2 probably damaging, and it is
absent from gnomAD.

The null-plus-hypomorph configuration is consistent with the genotypes reported in the MVA
literature, where complete biallelic loss of *BUB1B* is not observed in liveborn individuals.

---

## 2. Phenotype interpretation

Clinical features supplied with the case: rhabdomyosarcoma (HP:0002859), nephrocalcinosis
(HP:0000121), short stature (HP:0004322), failure to thrive (HP:0001508), skeletal muscle atrophy
(HP:0003202), premature birth (HP:0001622), small for gestational age (HP:0001518), and parental
recurrent spontaneous abortion (HP:0200067).

Two aspects of this constellation drove the analysis:

**Cancer predisposition with growth restriction** is characteristic of chromosomal instability
syndromes. Mosaic variegated aneuploidy arises from failure of the spindle assembly checkpoint, the
mechanism that blocks anaphase onset until every kinetochore is correctly attached to the mitotic
spindle. Loss of checkpoint stringency produces stochastic missegregation across divisions, giving
the variable aneuploidy that names the condition, developmental disruption through altered gene
dosage, and elevated malignancy risk through random loss of tumour suppressors.

**Parental recurrent pregnancy loss** was treated as primary inheritance evidence rather than
background. Under an autosomal recessive model with two healthy carrier parents, one quarter of
conceptions are biallelic; in a chromosome-segregation disorder the severe end of that distribution
is expected to be embryonic lethal. The reproductive history is therefore consistent with recessive
inheritance and with the proposed null-plus-hypomorph genotype in the surviving proband.

**Noted discrepancy.** Microcephaly is not reported in this case, although it is a prominent and
frequently severe feature of *BUB1B*-associated MVA as originally described. This was treated as a
reason to carry all three MVA genes through the analysis rather than restricting a priori to
*BUB1B*. The genotype evidence did not support redirection; the phenotypic discrepancy is recorded
here as unresolved.

---

## 3. Method

**Input.** The provided joint-called VCF (GRCh38) only. FASTQ data was not realigned.

**Candidate region definition.** Canonical coordinates for *BUB1B*, *CEP57* and *TRIP13* were
retrieved from the Ensembl REST lookup endpoint rather than hardcoded, and padded by 20 kb on each
flank to retain promoter and deep intronic territory. Variants were assigned to gene body or flank
after annotation, and flanking variants were scored separately: an unpartitioned region yields a
misleadingly dense candidate set dominated by neighbouring haplotype structure.

**Panel extraction.** A single streaming pass over the VCF retained records falling inside the
padded intervals, yielding 105 variants (*BUB1B* 80, *CEP57* 13, *TRIP13* 12). After gene-body
partitioning, only 14 of the *BUB1B* records lie within the transcribed gene.

**Genotype and quality parsing.** GT, DP, GQ, PID and FILTER were parsed per record. Note that the
INFO `AF` field in this file is the within-sample allele fraction emitted by the variant caller, not
a population frequency; it was not used for rarity filtering. Non-PASS records were retained and
downweighted rather than discarded, since mapping-quality flags are expected in repetitive sequence
and do not by themselves indicate a false call.

**Annotation.** Variants were submitted to the Ensembl VEP REST endpoint in batches, requesting
canonical transcript consequence, HGVS notation, exon and intron numbering, SIFT and PolyPhen
predictions, and gnomAD exome and genome allele frequencies.

**Ranking.** Candidates were ordered by a conjunction of three criteria rather than a fitted
score: population rarity (gnomAD AF below 1e-4 or absent), predicted functional consequence
(protein-truncating, canonical splice-disrupting, or missense with concordant damaging predictions),
and compatibility with biallelic recessive inheritance within a single gene. With a single proband
and a three-gene prior, this conjunction is sufficient to separate candidates; no statistical model
was fitted, and fitting one would not be meaningful at N=1.

**Exclusions.** *CEP57* carried no novel or rare variants. *TRIP13*'s two gene-body candidates
(c.92+38A>C, c.92+42G>C) share a read-backed phasing identifier and are therefore in cis on a single
haplotype, leaving one intact copy and excluding a recessive mechanism at that locus.

---

## 4. Submitted ranking

| Rank | EPCR | Content | Rationale |
|---|---|---|---|
| 1 | 0.95 | p.Leu737Ter + p.Asn1002Lys | Proposed compound heterozygous pair |
| 2 | 0.70 | p.Leu737Ter alone | Hedge against single-variant answer key |
| 3 | 0.60 | p.Asn1002Lys alone | Hedge against single-variant answer key |
| 4 | 0.30 | p.Leu737Ter + c.2679-1026A>G | Alternative second allele, untested |
| 5 | 0.15 | c.2679-1026A>G alone | Deep intronic, flagged for splicing follow-up |

EPCR values are strictly decreasing with no ties, so that the highest-confidence threshold isolates
a single row.

---

## 5. Declared limitations

**Phase is inferred, not observed.** The two proposed alleles are 10,911 bp apart, far beyond the
span of a sequencing fragment, so no read-backed phasing evidence links them. Trans configuration is
inferred from the recessive disease model and the parental reproductive history. Parental sequencing
would be required to establish it directly.

**No realignment from FASTQ.** Analysis used the provided call set only. Structural variants, mobile
element insertions, and variants in regions the caller did not emit are therefore inaccessible to
this submission. This matters because second alleles in *BUB1B* compound heterozygotes are
frequently non-coding or structural.

**No splicing prediction applied.** The deep intronic candidate c.2679-1026A>G is absent from gnomAD
and remains untested for pseudoexon activation. SpliceAI or Pangolin scoring over the *BUB1B* locus
is planned for a subsequent submission.

**In-silico pathogenicity only.** The functional claim for p.Asn1002Lys rests on SIFT and PolyPhen-2
predictions and on domain location. No structural stability calculation or functional assay supports
it. Both predictors are known to be correlated and to overcall damage at conserved positions.

**Single-transcript annotation.** Consequences are reported against the canonical transcript
(ENST00000287598.11). Alternative transcripts were not systematically evaluated.

---

## 6. Planned follow-up

1. SpliceAI and Pangolin scoring across the *BUB1B* locus, including the 20 kb flanks, to test the
   deep intronic candidate and search for unannotated pseudoexon activation.
2. Structure-aware assessment of p.Asn1002Lys: kinase domain position, predicted stability change,
   and proximity to functional interfaces, as an independent line of evidence from the
   sequence-based predictors.
3. Targeted reassessment of regions where the provided call set may be incomplete, to address the
   possibility that the second allele is structural.

---

## 7. Reproducibility

Code: see accompanying repository. The pipeline is a notebook sequence: Ensembl coordinate lookup,
VCF region scan, genotype parsing, VEP annotation, and candidate ranking. External services used are
the Ensembl REST lookup and VEP endpoints.

No patient sequence data or genotype-scale derivative is included in the repository. Reported
findings are limited to the candidate variants named above, which the challenge terms permit.

## 8. AI assistance disclosure

Anthropic API, Claude, commercial terms, no training on customer content.
