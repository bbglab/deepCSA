# VAF vs Depth Plots

## Overview

The VAF vs Depth plotting module generates scatter plots showing the relationship between sequencing depth and various quantitative metrics. These plots help identify patterns in mutation data and assess whether mutations follow expected depth-dependent patterns.

## Features

The module creates plots for:

1. **VAF (Variant Allele Frequency) vs Depth per Site**
   - Scatter plot of individual mutations
   - Overlay of hyperbolic curves (N/depth for N=1,2,3...)
   - Counts of mutations following each curve

2. **Mutation Density vs Average Depth per Gene**
   - Separate plots for different mutation types:
     - All mutations
     - SNV only
     - Insertions/Deletions
     - Protein-affecting vs non-protein-affecting
   - Reference curves showing expected densities

3. **Omega (dN/dS) vs Average Depth per Gene**
   - Color-coded by statistical significance (p < 0.05)
   - Separate plots for missense and truncating mutations
   - Reference line at neutral selection (ω=1)

4. **OncodriveFML Score vs Average Depth per Gene**
   - Color-coded by q-value significance
   - Identifies genes with functional bias

## Hyperbolic Curves

The hyperbolic curves represent the expected VAF for N mutant reads at various depths:

```
VAF = N / depth
```

Where:
- N = number of mutant reads (1, 2, 3, ...)
- depth = total sequencing depth at that position

These curves help identify:
- Low-frequency mutations (N=1,2,3)
- Technical artifacts vs true biological signals
- Sequencing depth effects on variant detection

## Usage

### Command Line

```bash
plot_vaf_depth.py \
  --sample_name SAMPLE_ID \
  --maf_file mutations.maf \
  --mutdensity_file mutdensity.tsv \
  --depth_file depth_per_gene.tsv \
  --omega_file omega_results.tsv \
  --oncodrivefml_file oncodrivefml_results.tsv \
  --output_prefix output_name \
  --max_n 10
```

### In Nextflow Pipeline

```groovy
include { PLOT_VAF_DEPTH } from './modules/local/plot/vaf_depth/main'

workflow {
    // Prepare input channels
    maf_ch = Channel.fromPath(params.maf)
        .map { file -> [[id: 'sample1'], file, null, null, null] }
    depth_ch = Channel.fromPath(params.depth)
        .map { file -> [[id: 'sample1'], file] }
    
    // Run plotting
    PLOT_VAF_DEPTH(maf_ch, depth_ch)
}
```

## Input Files

### MAF File (Required for VAF plots)
Tab-separated file with columns:
- `CHROM`: Chromosome
- `POS`: Position
- `DEPTH`: Total sequencing depth
- `VAF`: Variant allele frequency (0-1)
- `ALT_DEPTH`: Number of reads supporting alternate allele
- `GENE`: Gene name

### Mutation Density File (Required for density plots)
Tab-separated file with columns:
- `SAMPLE_ID`: Sample identifier
- `GENE`: Gene name
- `MUTTYPES`: Type of mutations (e.g., 'all_types', 'SNV', 'DELETION')
- `MUTDENSITY_MB`: Mutation density per megabase
- `N_MUTS`: Number of mutations

### Depth Per Gene File (Required for gene-level plots)
Tab-separated file with columns:
- `GENE`: Gene name
- `SAMPLE_ID`: Sample identifier
- `MEAN_GENE_DEPTH`: Average sequencing depth across gene

### Omega File (Optional)
Tab-separated file with columns:
- `gene`: Gene name
- `impact`: Impact type ('missense' or 'truncating')
- `dnds`: dN/dS ratio (omega)
- `pvalue`: Statistical significance
- `mutations`: Number of mutations

### OncodriveFML File (Optional)
Tab-separated file with columns:
- `SYMBOL`: Gene symbol
- `SCORE`: OncodriveFML score
- `PVALUE`: P-value
- `QVALUE`: FDR-corrected q-value

## Output

A single PDF file containing all generated plots:
- `{output_prefix}.vaf_depth_plots.pdf`

Each plot includes:
- Scatter points for data
- Hyperbolic curves (for VAF plots)
- Reference lines (for selection plots)
- Annotation boxes with counts and statistics
- Proper axis labels and titles

## Interpretation

### VAF vs Depth Plots

**Mutations on low N curves (N=1,2,3)**
- May represent rare variants or technical artifacts
- Higher sequencing depth improves detection

**Mutations deviating from curves**
- Could indicate:
  - Copy number alterations
  - Subclonal populations
  - Technical issues with depth calculation

### Mutation Density Plots

**High density at low depth**
- May indicate:
  - Poor sequencing quality in specific genes
  - Regions prone to sequencing errors

**Consistent density across depths**
- Suggests uniform mutation processes

### Selection Plots (Omega)

**ω > 1** (red points if p < 0.05)
- Positive selection
- More non-synonymous than expected

**ω < 1**
- Purifying selection
- Fewer non-synonymous than expected

**ω ≈ 1**
- Neutral evolution

## Parameters

- `--sample_name`: Sample identifier (required)
- `--maf_file`: Path to MAF file (optional)
- `--mutdensity_file`: Path to mutation density file (optional)
- `--depth_file`: Path to depth per gene file (optional)
- `--omega_file`: Path to omega results (optional)
- `--oncodrivefml_file`: Path to OncodriveFML results (optional)
- `--output_prefix`: Output file prefix (required)
- `--max_n`: Maximum N for hyperbolic curves (default: 10)

## Example

```bash
# Generate all available plots
plot_vaf_depth.py \
  --sample_name bladder_sample_01 \
  --maf_file results/bladder_sample_01.maf \
  --mutdensity_file results/bladder_sample_01.mutdensities.tsv \
  --depth_file results/bladder_sample_01.depth_per_gene.tsv \
  --omega_file results/bladder_sample_01.omega.tsv \
  --oncodrivefml_file results/bladder_sample_01.oncodrivefml.tsv \
  --output_prefix bladder_sample_01 \
  --max_n 15
```

This will create `bladder_sample_01.vaf_depth_plots.pdf` containing all plots.

## Notes

- At least one input file (MAF, mutation density, omega, or OncodriveFML) must be provided
- The script will generate only the plots for which input files are provided
- Depth per gene file is required for mutation density, omega, and OncodriveFML plots
- All plots use consistent styling from `utils_plot.py`
- Hyperbolic curves are semi-transparent to avoid obscuring data points
