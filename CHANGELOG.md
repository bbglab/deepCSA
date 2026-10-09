# bbglab/deepCSA: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of bbglab/deepCSA, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- New output `panel_exons_protein_intervals.tsv` in the `DNA2PROTEINMAPPING` step (published to `regions/annotations/`), reporting for each exon of the panel transcripts the interval of protein coordinates it covers, using the same structure as the domains info file (`Ens_Transcr_ID`, `Begin`, `End`, `NAME`, `GENE`, `DOMAIN_ID`).
- Saturation kinetics curves (`COMPUTE_SATURATION_KINETICS` step, published to `plots/saturation_kinetics/`). For each group and each (resolution, impact) combination — genomic/residue × protein_affecting/nonsense/truncating/missense/synonymous — the step produces:
  - `{group}.curves/{sites}_{impact}_empirical.pdf` — empirical discovery index curves (proportion of mutated sites vs sequencing depth) for all genes, obtained by downsampling the observed mutations with Bernoulli replicates.
  - `{group}.curves/{sites}_{impact}_theoretical_empirical.pdf` — the same empirical curves overlaid with the theoretical neutral saturation curve derived from the per-site relative mutability and the synonymous mutation rate.
  - `{group}.curves/{sites}_{impact}_slopes.pdf` — per-gene comparison of the rate of change (Δ proportion / Δ log10 depth) of the empirical curve against the theoretical neutral curve over identical depth intervals.
  - `{group}_mutations_{sites}_rates.{impact}.tsv` — per-gene/per-site unique-mutation probabilities at each subsampling depth.
  - `{group}_slopes_{sites}.{impact}.tsv` — per-gene interval slopes (empirical, theoretical and their ratio) for cross-run comparison.
  - Requires both `--omega` and one of the mutability-driven analyses (`--oncodrivefml`, `--oncodriveclustl` or `--oncodrive3d`) to be enabled, since it consumes the omega preprocessing mutability table and the relative mutability per site.

### `Fixed`

### `Dependencies`

### `Deprecated`
