#!/bin/bash
# Example script for running VAF vs depth plots
#
# This script demonstrates how to use the plot_vaf_depth.py tool
# with different combinations of input files.

set -e

SAMPLE_NAME="example_sample"
OUTPUT_DIR="./results"
mkdir -p "$OUTPUT_DIR"

# Example 1: Plot only VAF vs depth (requires MAF file only)
echo "Example 1: VAF vs depth per site"
plot_vaf_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --maf_file data/sample.maf \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_vaf_only" \
  --max_n 10

# Example 2: Plot mutation density vs depth (requires mutation density + depth files)
echo "Example 2: Mutation density vs depth"
plot_vaf_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --mutdensity_file data/sample.mutdensities.tsv \
  --depth_file data/sample.depth_per_gene.tsv \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_mutdensity" \
  --max_n 10

# Example 3: Plot all metrics together
echo "Example 3: All metrics combined"
plot_vaf_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --maf_file data/sample.maf \
  --mutdensity_file data/sample.mutdensities.tsv \
  --depth_file data/sample.depth_per_gene.tsv \
  --omega_file data/sample.omega.tsv \
  --oncodrivefml_file data/sample.oncodrivefml.tsv \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_all" \
  --max_n 15

echo "All plots generated successfully!"
echo "Output files:"
ls -lh "$OUTPUT_DIR"/*.vaf_depth_plots.pdf
