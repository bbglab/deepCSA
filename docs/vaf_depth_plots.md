# VAF and Selection Metrics vs Depth Plots

## Overview

The depth relationship plotting modules generate scatter plots showing the relationship between sequencing depth and various quantitative metrics. These plots help identify patterns in mutation data and assess whether mutations follow expected depth-dependent patterns.

## Features

The implementation consists of two independent modules:

### Module 1: VAF and Mutation Density (with hyperbolic curves)

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

### Module 2: Selection Metrics (without hyperbolic curves)

3. **Omega (dN/dS) vs Average Depth per Gene**
   - Color-coded by statistical significance (p < 0.05)
   - Separate plots for missense and truncating mutations
   - Reference line at neutral selection (ω=1)
   - **No hyperbolic curves**

4. **OncodriveFML Score vs Average Depth per Gene**
   - Color-coded by q-value significance
   - Identifies genes with functional bias
   - **No hyperbolic curves**

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

**For VAF and Mutation Density (with hyperbolic curves):**
```bash
plot_vaf_mutdensity_depth.py \
  --sample_name SAMPLE_ID \
  --maf_file mutations.maf \
  --mutdensity_file mutdensity.tsv \
  --depth_file depth_per_gene.tsv \
  --output_prefix output_name \
  --max_n 10
```

**For Selection Metrics (without hyperbolic curves):**
```bash
plot_selection_depth.py \
  --sample_name SAMPLE_ID \
  --omega_file omega_results.tsv \
  --oncodrivefml_file oncodrivefml_results.tsv \
  --depth_file depth_per_gene.tsv \
  --output_prefix output_name
```

### In Nextflow Pipeline

The modules are integrated via the `PLOT_DEPTH_RELATIONSHIPS` subworkflow. Enable with:

```yaml
# In nextflow.config or params file
plot_depth_relationships = true
```

The workflow automatically runs the appropriate modules based on available data:
- Runs VAF/mutation density module if MAF or mutation density files exist
- Runs selection metrics module if omega or OncodriveFML results exist
- Each module operates independently

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

**For VAF and mutation density analysis:**
```bash
plot_vaf_mutdensity_depth.py \
  --sample_name bladder_sample_01 \
  --maf_file results/bladder_sample_01.maf \
  --mutdensity_file results/bladder_sample_01.mutdensities.tsv \
  --depth_file results/bladder_sample_01.depth_per_gene.tsv \
  --output_prefix bladder_sample_01 \
  --max_n 15
```
Creates: `bladder_sample_01.vaf_mutdensity_depth.pdf`

**For selection metrics analysis:**
```bash
plot_selection_depth.py \
  --sample_name bladder_sample_01 \
  --omega_file results/bladder_sample_01.omega.tsv \
  --oncodrivefml_file results/bladder_sample_01.oncodrivefml.tsv \
  --depth_file results/bladder_sample_01.depth_per_gene.tsv \
  --output_prefix bladder_sample_01
```
Creates: `bladder_sample_01.selection_depth.pdf`

## Notes

- The modules are independent and can run separately
- **VAF/mutation density module**: At least MAF or mutation density file required; generates plots with hyperbolic curves
- **Selection metrics module**: At least omega or OncodriveFML file required; generates plots WITHOUT hyperbolic curves
- Depth per gene file is required for mutation density, omega, and OncodriveFML plots
- All plots use consistent styling from `utils_plot.py`
- Hyperbolic curves appear only on VAF and mutation density plots, not on selection metrics plots
