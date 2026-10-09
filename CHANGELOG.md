# bbglab/deepCSA: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of bbglab/deepCSA, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- New output `panel_exons_protein_intervals.tsv` in the `DNA2PROTEINMAPPING` step (published to `regions/annotations/`), reporting for each exon of the panel transcripts the interval of protein coordinates it covers, using the same structure as the domains info file (`Ens_Transcr_ID`, `Begin`, `End`, `NAME`, `GENE`, `DOMAIN_ID`).

### `Fixed`

### `Dependencies`

### `Deprecated`
