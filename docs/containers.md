# Containers used in deepCSA

This document lists all the container images used across the deepCSA pipeline,
where they are declared, and where their build recipes can be found when
available.

Containers are declared in two ways:

1. **Per process** with the `container` directive inside each module's
   `main.nf` file.
2. **Per label** in `conf/modules.config`, using `withLabel` scopes. Processes
   carrying the label inherit the container, which avoids repeating the same
   image in dozens of modules.

> **Note on registries**: the default registry for Docker and Apptainer is set
> to `quay.io` in `nextflow.config` (`docker.registry` /
> `apptainer.registry`). Images declared with an explicit `docker.io/` prefix
> are pulled from Docker Hub regardless of this default; images without a
> registry prefix (e.g. `biocontainers/...`) resolve to the configured default
> registry.

---

## bbglab containers

Images maintained by the BBGLab. The Dockerfile recipes for the `bbglab/*`
images are tracked in the lab repository
[bbglab/containers-recipes](https://github.com/bbglab/containers-recipes).

| Container | URL | Used by | Build recipe |
|-----------|-----|---------|--------------|
| `bbglab/deepcsa-core:0.2.0` | https://hub.docker.com/r/bbglab/deepcsa-core | All processes labelled `deepcsa_core` (~63 local Python processes: mutation processing, QC plots, omega post-processing, saturation, etc.) | [bbglab/containers-recipes](https://github.com/bbglab/containers-recipes) |
| `bbglab/deepcsa_bed:latest` | https://hub.docker.com/r/bbglab/deepcsa-bed | Panel creation: `CREATEPANELS:CONSENSUS`, `CREATEPANELS:CAPTURED` (`modules/local/createpanels/`) | [bbglab/containers-recipes](https://github.com/bbglab/containers-recipes) |
| `bbglab/omega:0.2.1` | https://hub.docker.com/r/bbglab/omega | Omega estimation: `OMEGA:PREPROCESS`, `OMEGA:MUTABILITIES`, `OMEGA:ESTIMATOR` (`modules/local/bbgtools/omega/`) | [bbglab/containers-recipes](https://github.com/bbglab/containers-recipes) |

Declared in `conf/modules.config` (`deepcsa_core` label) and in the
corresponding module files.

---

## Collaborator / lab member containers

Images built and maintained by individual lab members on their Docker Hub
accounts. Build recipes are generally kept in personal repositories; where a
public recipe is known it is linked.

<!-- TODO start transferring all these images to the bbglab dockerhub and also make sure that we pair them with the right receipes from the bbglab/container-receipes repo -->

| Container | URL | Used by | Build recipe |
|-----------|-----|---------|--------------|
| `ferriolcalvet/oncodriveclustl:latest` | https://hub.docker.com/r/ferriolcalvet/oncodriveclustl | `ONCODRIVECLUSTL` (`modules/local/bbgtools/oncodriveclustl/main.nf`) | [bbglab/oncodriveclustl](https://github.com/bbglab/oncodriveclustl) |
| `ferriolcalvet/oncodrivefml:latest` | https://hub.docker.com/r/ferriolcalvet/oncodrivefml | `ONCODRIVEFML` (`modules/local/bbgtools/oncodrivefml/main.nf`) | [bbglab/oncodrivefml](https://github.com/bbglab/oncodrivefml) |
| `ferriolcalvet/dnds:latest` | https://hub.docker.com/r/ferriolcalvet/dnds | `DNDS:BUILDREF`, `DNDS:RUN` (`modules/local/dnds/`) | — |
| `ferriolcalvet/omegacovariates:v0.1.0` | https://hub.docker.com/r/ferriolcalvet/omegacovariates | `OMEGACOVARIATES:RUN` (`modules/local/omega_covariates/run/main.nf`) | — |
| `ferriolcalvet/runningr:v1` | https://hub.docker.com/r/ferriolcalvet/runningr | `WGSCALEDMUTDENSITY` (`modules/local/mut_density/wgscaled/main.nf`) | — |
| `ferriolcalvet/saturation:v0.1.0` | https://hub.docker.com/r/ferriolcalvet/saturation | `SATURATION_KINETICS:COMPUTE` (`modules/local/saturation_kinetics/compute/main.nf`) | — |
| `ferriolcalvet/hdp_wrapper` | https://hub.docker.com/r/ferriolcalvet/hdp_wrapper | HDP signature pipeline: `CHAINSAMPLING`, `COMPARESIGS`, `NORMALIZESIGS`, `PREPAREINPUT`, `PROCESSRESULTS` (`modules/local/signatures/hdp/`) | — |
| `ferriolcalvet/msighdp:latest` | https://hub.docker.com/r/ferriolcalvet/msighdp | `MSIGHDP` (`modules/local/signatures/msighdp/main.nf`) | — |
| `ferriolcalvet/musical:latest` | https://hub.docker.com/r/ferriolcalvet/musical | `MUSICAL` (`modules/local/signatures/musical/main.nf`) | — |
| `ferriolcalvet/sigprofiler_assignment:1.1.3` | https://hub.docker.com/r/ferriolcalvet/sigprofiler_assignment | `SIGPROFILER:ASSIGNMENT:COSMICFIT`, `SIGPROFILER:ASSIGNMENT:DECOMPOSEFIT` | — |
| `ferriolcalvet/sigprofilermatrixgenerator:1.3.5` | https://hub.docker.com/r/ferriolcalvet/sigprofilermatrixgenerator | `SIGPROFILER:MATRIXGENERATOR` | — |
| `ferriolcalvet/sigprofilerassignment` | https://hub.docker.com/r/ferriolcalvet/sigprofilerassignment | `SIGPROFILER:EXTRACTOR` | — |
| `spellegrini87/oncodrive3d:1.0.9-light` | https://hub.docker.com/r/spellegrini87/oncodrive3d | `ONCODRIVE3D:RUN`, `ONCODRIVE3D:PLOT` (`modules/local/bbgtools/oncodrive3d/`) | [bbglab/oncodrive3d](https://github.com/bbglab/oncodrive3d) |
| `spellegrini87/oncodrive3d:1.0.9-chimerax` | https://hub.docker.com/r/spellegrini87/oncodrive3d | `ONCODRIVE3D:PLOT_CHIMERAX` | [bbglab/oncodrive3d](https://github.com/bbglab/oncodrive3d) |
| `rblancomi/bbgregressions:dev` | https://hub.docker.com/r/rblancomi/bbgregressions | All processes labelled `bbgregressions` (regression analyses) | — |
| `ferranmuinos/test_mutated_genomes` | https://hub.docker.com/r/ferranmuinos/test_mutated_genomes | `MUTATED_CELLS_FROM_VAF`, `MUTATED_GENOMES_FROM_VAF` (`modules/local/mutated_*_from_vaf/`) | — |
| `axelrosendahlhuber/expected_mutrate:latest` | https://hub.docker.com/r/axelrosendahlhuber/expected_mutrate | `MUTATED_CELLS_EXPECTED` (`modules/local/mutated_cells_expected/main.nf`) | — |

> **Warning**: several of these images use the `latest` tag or have no tag at
> all, which makes runs non-reproducible. See the note on pinning at the end
> of this document.

---

## nf-core / biocontainers containers

Standard community images used by the nf-core-style modules. These are
well-established images published by the [Biocontainers](https://biocontainers.pro)
project; their recipes are generated from the corresponding Bioconda packages
in [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes).

| Container | URL | Used by | Build recipe |
|-----------|-----|---------|--------------|
| `biocontainers/multiqc:1.20--pyhdfd78af_0` (docker) / `biocontainers/multiqc:1.18--pyhdfd78af_0` (singularity) | https://quay.io/repository/biocontainers/multiqc | `MULTIQC` (`modules/nf-core/multiqc/main.nf`); also `DUMPSOFTWAREVERSIONS` (requires pyyaml) | [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/multiqc) |
| `biocontainers/tabix:1.11--hdfd78af_0` | https://quay.io/repository/biocontainers/tabix | `TABIX:BGZIPTABIX`, `TABIX:BGZIPTABIXQUERY` (`modules/nf-core/tabix/`) | [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/tabix) |
| `biocontainers/samtools:1.18--h50ea8bc_1` | https://quay.io/repository/biocontainers/samtools | `COMPUTEDEPTHS` (`modules/local/computedepths/main.nf`) | [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/samtools) |
| `biocontainers/pybedtools:0.9.1--py38he0f268d_0` | https://quay.io/repository/biocontainers/pybedtools | `CREATEPANELS:CUSTOMBEDFILE` (`modules/local/createpanels/custombedfile/main.nf`) | [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/pybedtools) |
| `biocontainers/python:3.8.3` | https://quay.io/repository/biocontainers/python | `CHECK_SAMPLESHEET` (`modules/local/samplesheet_check.nf`) | [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/python) |

Several of these modules use the conditional pattern:

```groovy
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'biocontainers/tool:version' :
    'docker://biocontainers/tool:version' }"
```

which selects a Docker Hub-prefixed reference when running under
Singularity/Apptainer and the plain reference otherwise.

---

## Reproducibility notes

- **Unpinned images**: `ferriolcalvet/oncodriveclustl`, `oncodrivefml`,
  `dnds`, `msighdp`, `musical`, `hdp_wrapper`, `sigprofilerassignment`,
  `ferranmuinos/test_mutated_genomes`, `axelrosendahlhuber/expected_mutrate`
  and both `bbglab/deepcsa_bed` usages resolve to `latest` or no tag. Pinning
  them to immutable versions (as already done for `bbglab/omega:0.2.1`,
  `sigprofilermatrixgenerator:1.3.5` or `oncodrive3d:1.0.9-light`) is
  recommended for reproducibility.
- **`rblancomi/bbgregressions:dev`** is a development tag; the config itself
  carries a FIXME asking for an update to a bbglab registry image.
- **Registry mirror**: if a strict reproducibility setup is required, all
  external images should be mirrored to an internal registry and
  `docker.registry` / `apptainer.registry` in `nextflow.config` pointed at it.
