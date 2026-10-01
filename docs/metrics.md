# bbglab/deepCSA: Computed metrics and how to interpret them

This document explains **what deepCSA computes** and **how each output should be interpreted**. It is organised by analysis layer, following the order in which the pipeline processes your data. For the file layout of every output directory see [Output](output.md), and for the mathematical details of individual methods see [Tools](tools.md).

> **Reading order suggestion:** if you are new to the pipeline, read the [Interpretation guide](#interpretation-guide-where-to-start) at the bottom of this document first — it maps common biological questions to the outputs you should look at.

## Table of contents

- [1. Depth and panel definition](#1-depth-and-panel-definition)
- [2. Mutation filtering and somatic calling](#2-mutation-filtering-and-somatic-calling)
- [3. Mutation burden: mutation density](#3-mutation-burden-mutation-density)
- [4. Mutational processes: profiles and signatures](#4-mutational-processes-profiles-and-signatures)
- [5. Positive selection](#5-positive-selection)
- [6. Clonal structure: mutated genomes and cells](#6-clonal-structure-mutated-genomes-and-cells)
- [7. Interindividual variability: regressions](#7-interindividual-variability-regressions)
- [8. Quality control](#8-quality-control)
- [Interpretation guide: where to start](#interpretation-guide-where-to-start)

---

## 1. Depth and panel definition

**What is computed.** Per-position sequencing depth for every sample, and the set of genomic positions ("panels") that are well covered enough to be analysed.

**Where to find it.**

| Output | Content |
|---|---|
| `depths/individual/` | Per-position depth table, one column per sample. |
| `depths/summary/` | Average depth per sample/gene for three region definitions: `exons`, `exons_cons` (exons within the consensus panel), `all_cons` (all well-covered regions). |
| `depths/plots_per_group/` | Depth distribution plots per sample and per group. |
| `regions/consensuspanels/` | Cohort-level consensus panel (positions covered ≥ `consensus_panel_min_depth` in ≥ 80% of samples). |
| `regions/samplepanels/` | Per-sample panels for each region type (exons, introns, protein-affecting, non-protein-affecting, synonymous). |
| `regions/expandedregions/` | Subgenic / domain / exon expansions used by omega and dNdScv. |

**How to interpret it.**

- **Depth is the foundation of every downstream metric.** Mutation density, mutability, and selection estimates are all normalised by depth, so a sample with poor or uneven coverage will produce unreliable values. Always check `depths/summary/` and the depth plots before interpreting anything else.
- The **consensus panel** defines the genomic space in which cohort-level comparisons are made. If a gene is absent from the consensus panel it will not be included in cohort-level selection metrics.
- Note that the depth reported per mutation (used for VAF) is **N-discounted** (N bases are not counted), while the depth used for density normalisation is not. This is expected and documented in [Output — Depth analysis](output.md#depth-analysis).

## 2. Mutation filtering and somatic calling

**What is computed.** Each input mutation is annotated (VEP), then classified and filtered at two levels:

- **Sample-level filters** — e.g. VAF distortion (`VAF_AM / VAF > vaf_distortion_threshold`), low depth, N-rich context, lack of pileup support.
- **Cohort-level filters** — e.g. `other_sample_SNP` (the "somatic" mutation is a common SNP in another sample), `repetitive_variant` (seen in too many samples), `not_covered`, `not_in_exons`.

Mutations are labelled **germline** (present in multiple samples, likely inherited) or **somatic** (private to a sample).

**Where to find it.**

| Output | Content |
|---|---|
| `mutations/germline_somatic/` | All calls with their germline/somatic label. |
| `mutations/clean_somatic/` | Somatic calls after all filtering — **this is the main mutation table used downstream**. |
| `mutations/clean_germline_somatic/` | Cleaned germline + somatic calls. |
| `processing_files/flagged_positions/` | Positions flagged by the cohort-level filters. |
| `plots/mutations_summary/` | PDFs with per-sample/per-gene mutation counts, filter statistics, and the distribution of mutations per number of mutated reads (overall and per mutation type). |

**How to interpret it.**

- The **number of mutations removed by each filter** is shown in the `mutations_summary` PDFs. A large fraction removed by `other_sample_SNP` usually indicates shared germline variants that were not filtered upstream; a large fraction removed by depth/VAF filters may indicate low-quality calling.
- Downstream analyses (density, profiles, selection) use the **clean somatic** set. If your biological question concerns germline variation, use `clean_germline_somatic` instead.
- The **needle plots** (`plots/needle_plots/`) show, per sample and per gene, the VAF of each mutation against the gene's depth — useful to visually inspect clonal substructure (a mutation present in a subset of cells shows a lower VAF than the dominant clone).

## 3. Mutation burden: mutation density

**What is computed.** The number of mutations per megabase of sequenced DNA, per sample, per gene, and per consequence-type group (all types, SNVs, indels, protein-affecting, non-protein-affecting, synonymous). Two flavours are produced:

- **Flat density** (`mutdensity/`): `N_MUTS / DEPTH × 10⁶`, i.e. mutations per Mb of sequenced DNA. Also reported per mutated read (`MUTREADSDENSITY_MB`).
- **Adjusted density** (`mutdensity_adjusted/`): corrects the flat density for the **trinucleotide composition** of the analysed sites, so that two regions with the same mutational process but different triplet content become comparable. See [Tools — Adjusted mutation density](tools.md#adjusted-mutation-density) for the full derivation.

**Where to find it.**

| Output | Content |
|---|---|
| `mutdensity/individual_vals/` | Flat density per sample/gene/region. |
| `mutdensity_adjusted/individual_vals/` | Trinucleotide-adjusted density per sample/gene/region. |
| `qc/mutdensityqc/` | QC plots of density per sample and per gene. |

**How to interpret it.**

- **Flat density** answers "how many mutations per Mb did this sample accumulate in this region?". It is the right metric to compare **overall mutational burden** between samples (e.g. smokers vs non-smokers).
- **Adjusted density** answers "how many mutations per Mb, given the mutational process acting on these sites?". Use it when comparing **different genomic regions or consequence classes** (e.g. missense vs synonymous sites), because it removes the confounding effect of trinucleotide context.
- A sample whose density collapses at low depth is a red flag — check `qc/metrics_vs_depth/` (see [Quality control](#8-quality-control)).

## 4. Mutational processes: profiles and signatures

**What is computed.**

- **Mutational profiles** (`mutational_profile/`): the probability of each of the 96 SBS trinucleotide contexts (e.g. `C>T` in a `TCA` context), computed in up to four normalisation conditions: all regions, exons only, non-protein-affecting regions, and introns/intergenic regions. Profiles are computed per sample and for the cohort, with optional Bayesian shrinkage for low-burden samples.
- **Mutational signatures**:
  - `signatures/sigprofilerassignment/` — assignment of known COSMIC signatures to the cohort (and per sample), with activity tables and plots.
  - `signatures/signatures_hdp/` + `signatures/hdp_decomposition_spa/` — de novo signature extraction with a Hierarchical Dirichlet Process, followed by reassignment of the extracted signatures to samples.
  - `signatures/sigprofilerassignment_indels/` — the same assignment for indels.
  - `processing_files/mutations_matrix/` — per-sample SBS count matrices, ready for external tools such as MSA.

**How to interpret it.**

- The **profile** is a fingerprint of the mutational processes acting in your samples. Compare profiles across groups (e.g. by smoking status) to identify which processes differ.
- **Signature assignment** decomposes each sample's profile into known processes (e.g. SBS4/SBS5 for smoking, SBS1 for age-related deamination). The activity table tells you the relative contribution of each signature.
- **De novo extraction (HDP)** is useful when your samples are expected to harbour processes not in COSMIC (e.g. in non-cancer tissues or unusual exposures).
- **Caveat:** signature estimates are unstable in low-burden samples. Check the profile stability files (`*.profile_stability.tsv`) and the `qc/mutational_profiles_comparison/` clustermap before over-interpreting a single sample.

## 5. Positive selection

This is the core of deepCSA: a battery of complementary metrics that test whether mutations in a gene (or sub-region) are **enriched relative to the neutral expectation**. No single metric is sufficient on its own — they make different assumptions, and concordance between them strengthens a selection call.

### 5.1 Omega (dN/dS)

**What is computed.** A dN/dS-style ratio per gene (and per subgenic region) comparing the observed number of non-synonymous mutations to the expected number given the mutational profile and the number of synonymous sites. Two variants:

- `selection/omega/` — uses **per-sample** mutational profiles and synonymous rates.
- `selection/omegagloballoc/` — uses a **global cohort** profile and synonymous rates, which stabilises estimates in low-burden samples.

An optional covariate model (`selection/omega_covariates/`) regresses omega against sample-level covariates.

**Key columns** (in `estimator/output_mle.<sample>.tsv` and the cohort tables):

| Column | Meaning |
|---|---|
| `dnds` | Point estimate of the selection ratio. |
| `lower` / `upper` | Confidence interval. |
| `pvalue` | P-value for `dnds > 1`. |
| `pvalue_adj` | Benjamini-Hochberg adjusted p-value (corrected separately for cohort, groups, and per-sample sets, and for gene-level vs subgenic regions). |

**How to interpret it.**

- `dnds ≈ 1`: neutral evolution (mutations accumulate as expected).
- `dnds > 1` with a significant `pvalue_adj`: **positive selection** — non-synonymous mutations are enriched.
- `dnds < 1`: purifying selection or insufficient power.
- Use `omega` for sample-specific signals and `omegagloballoc` for conservative cohort-level estimates (see [Tools — Interpreting outputs](tools.md#interpreting-outputs-sanity-checks-and-key-metrics)).
- Genes with very few mutations produce wide confidence intervals; treat non-significant estimates in low-burden samples as "no evidence" rather than "no selection".

### 5.2 Site comparison

**What is computed.** For each site (or amino-acid residue / residue change), the observed number of mutations is compared to the expected number from the mutability model.

**Where to find it.** `selection/sitecomparison/` — eight combinations of background model (`single`, `multi`, `glocsingle`, `glocmulti`) and counting scheme (`count_single`, `count_multi`).

**Key columns:** `OBSERVED_MUTS`, `EXPECTED_MUTS`, `OBS/EXP` (enrichment ratio), `p_value` (Poisson).

**How to interpret it.**

- This is the **finest-resolution** selection metric: it tells you *which specific sites* are under selection, not just which genes.
- Recommended starting points: `bckg_single_count_single` for cohort-level reporting; `bckg_single_count_multi` / `bckg_multi_count_multi` when you want to account for multiple occurrences of the same mutation (e.g. recurrent hotspot mutations).
- A site with `OBS/EXP >> 1` and a small `p_value` is a candidate **driver mutation / hotspot**.

### 5.3 dNdScv

**What is computed.** The [dNdScv](https://github.com/im3sanger/dndscv) R package, run with a per-run reference CDS built dynamically from your panel. It models mutation rates with covariates (epigenomic features, trinucleotide context, ...) and estimates dN/dS per gene.

**Where to find it.**

| Output | Content |
|---|---|
| `selection/dndscv/cv/` | `*.cv.tsv` — per-gene dN/dS with covariate-adjusted estimates. |
| `selection/dndscv/persample/` | `*.globaldnds.tsv` — per-sample global dN/dS. |
| `selection/dndscv/local/` | `*.loc.tsv` — local (per-gene) dN/dS. |

**How to interpret it.** dNdScv is the most covariate-aware dN/dS estimator in the pipeline; it is the best choice when you suspect that mutation rates vary systematically across the genome (e.g. by replication timing or chromatin state). Compare its per-gene calls with omega: concordant genes are the most robust selection candidates.

### 5.4 dN/dS proxy

**What is computed.** A fast per-gene ratio of non-synonymous vs synonymous **adjusted mutation densities** (`selection/dndsproxy/*.gene_mutdensities_n_dnds.tsv`).

**How to interpret it.** A quick sanity check only — it has no significance testing. Use it to get a first impression of which genes look selected before running the full omega/dNdScv analyses.

### 5.5 OncodriveFML

**What is computed.** A functional-impact-based selection score per gene that combines the excess of protein-affecting mutations with CADD scores of the affected residues.

**Where to find it.** `selection/oncodrivefml/`.

**How to interpret it.** OncodriveFML is sensitive to *which residues* are hit, not just how many mutations occur. A gene with a modest omega but a high OncodriveFML score is likely accumulating a few highly deleterious mutations (e.g. in a critical domain).

### 5.6 Oncodrive3D

**What is computed.** Selection analysis in 3D protein space: mutations are mapped onto protein structures and clustered, testing whether selected mutations cluster in specific structural regions.

**Where to find it.** `selection/oncodrive3d/run/` (per-sample results) and `plots/selection/oncodrive3d/chimerax/` (3D visualisations).

**How to interpret it.** Complements the sequence-based metrics: it can detect selection acting on structurally clustered residues that are far apart in the sequence. The ChimeraX plots let you visually inspect the clustering.

### 5.7 Indel selection

**What is computed.** Selection analysis restricted to indels (`indels = true`).

**How to interpret it.** Indels are rarer and harder to call than SNVs, so indel selection calls should be interpreted with more caution than SNV-based ones.

## 6. Clonal structure: mutated genomes and cells

**What is computed.** From the VAF distribution of each mutation, deepCSA estimates how many **mutated genomes** (and, with `mutated_cells_vaf`, how many **mutated cells**) carry each mutation, and summarises the clonal architecture per sample.

**Where to find it.** Subdirectories under `selection/` and `mutations/` (see [Output — Additional clonal structure metrics](output.md#additional-clonal-structure-metrics)).

**How to interpret it.**

- A mutation with VAF ≈ 50% in a diploid sample is likely present in **all** cells (clonal); a mutation with VAF ≈ 25% is present in about **half** the cells (subclonal).
- The distribution of mutated-genome counts per sample describes its **clonal architecture**: a sample dominated by one high-VAF mutation has a simple architecture, while a sample with many intermediate-VAF mutations has a complex, multi-clone architecture.
- These estimates assume a known ploidy and no copy-number distortion; interpret them as approximations.

## 7. Interindividual variability: regressions

**What is computed.** Univariate and multivariate linear regressions between clonal-structure/selection metrics (mutation density, omega, ...) and sample-level covariates from your feature groups table (e.g. age, sex, smoking status).

**Where to find it.** `regressions/` (model tables and plots).

**How to interpret it.**

- Each regression row gives a coefficient, standard error, and p-value for one covariate's effect on one metric.
- A significant positive coefficient for "smoking" on mutation density means smokers accumulate more mutations per Mb.
- Use the **group-level** regressions (from `features_groups_list`) to test whether a covariate's effect differs between subgroups.
- Regressions are only as good as the covariates you provide — see [File formatting — Feature groups](file_formatting.md#feature-groups).

## 8. Quality control

The `qc/` directory collects all quality-control views. **Run through this list before trusting any biological conclusion.**

| QC output | What it checks | What to look for |
|---|---|---|
| `qc/metrics_vs_depth/` | Whether mutation density and omega estimates are confounded by sequencing depth. | Samples/genes whose metric values drift with depth should be treated with caution (or excluded). |
| `qc/mutdensityqc/` | Density distributions per sample and per gene. | Outlier samples with extreme burden. |
| `qc/trinucleotide_proportions/` | Trinucleotide composition of the panel vs the genome. | Large deviations indicate a biased capture. |
| `qc/mutational_profiles_comparison/` | Similarity of mutational profiles across samples (clustermap). | Samples that cluster apart may have different mutational processes or quality issues. |
| `qc/mutationspecific/` | Per-mutation QC (VAF distributions, read support). | Mutations with unusual VAF or low read support. |
| `qc/omega_flagged/` | Genes/samples with unstable omega estimates. | Flagged entries should not be reported as selection calls. |
| `qc/evaluate_omega_globalloc/` | Agreement between per-sample and global-loc omega. | Large disagreement indicates low-burden or profile-unstable samples. |
| `qc/contamination/` | Cross-sample contamination estimates. | High contamination invalidates the somatic/germline classification for the affected samples. |

---

## Interpretation guide: where to start

A quick map from biological question to the outputs you should look at:

| Question | Start here | Then check |
|---|---|---|
| Is my data of good quality? | `depths/summary/`, `qc/` (all) | `plots/mutations_summary/` (filter stats) |
| How many mutations does each sample carry? | `mutdensity/individual_vals/` | `plots/mutations_summary/` (per-sample counts) |
| Which mutational processes are active? | `mutational_profile/`, `signatures/sigprofilerassignment/` | `qc/mutational_profiles_comparison/` |
| Which genes show positive selection? | `selection/omega/` + `selection/omegagloballoc/` | `selection/dndscv/cv/`, `selection/oncodrivefml/` |
| Which specific sites/residues are selected? | `selection/sitecomparison/` | `plots/selection/oncodrive3d/chimerax/` |
| What is the clonal architecture of each sample? | Mutated-genomes outputs under `selection/` | `plots/needle_plots/` |
| Does a covariate (age, smoking, ...) explain variation? | `regressions/` | `plots/interindividual_variability/` |
| Are two groups of samples different? | Group-level density/omega tables (from `features_groups_list`) | `plots/selection_summary/` |

### A note on statistical significance

- Selection metrics report **Benjamini-Hochberg adjusted** p-values (`pvalue_adj`). Use `pvalue_adj < 0.05` (or your chosen FDR threshold) as the primary significance criterion, not the raw `pvalue`.
- Multiple selection metrics are computed on the same data; **concordance** between omega, dNdScv, OncodriveFML, and site comparison is the strongest evidence for a true selection signal.
- Low-burden samples (few mutations) have low power: non-significant results there mean "not enough evidence", not "no selection".
