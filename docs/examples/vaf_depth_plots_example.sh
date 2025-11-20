#!/bin/bash
# Example script for running VAF and selection metric vs depth plots
#
# This script demonstrates how to use the independent plotting modules
# with different combinations of input files.

set -e

SAMPLE_NAME="example_sample"
OUTPUT_DIR="./results"
mkdir -p "$OUTPUT_DIR"

# Example 1: Plot only VAF vs depth (requires MAF file only)
echo "Example 1: VAF vs depth per site (with hyperbolic curves)"
plot_vaf_mutdensity_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --maf_file data/sample.maf \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_vaf_only" \
  --max_n 10

# Example 2: Plot mutation density vs depth (requires mutation density + depth files)
echo "Example 2: Mutation density vs depth (with hyperbolic curves)"
plot_vaf_mutdensity_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --mutdensity_file data/sample.mutdensities.tsv \
  --depth_file data/sample.depth_per_gene.tsv \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_mutdensity" \
  --max_n 10

# Example 3: Plot VAF and mutation density together
echo "Example 3: VAF and mutation density combined (with hyperbolic curves)"
plot_vaf_mutdensity_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --maf_file data/sample.maf \
  --mutdensity_file data/sample.mutdensities.tsv \
  --depth_file data/sample.depth_per_gene.tsv \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_vaf_mutdens" \
  --max_n 15

# Example 4: Plot selection metrics (omega and OncodriveFML)
echo "Example 4: Selection metrics vs depth (WITHOUT hyperbolic curves)"
plot_selection_depth.py \
  --sample_name "$SAMPLE_NAME" \
  --omega_file data/sample.omega.tsv \
  --oncodrivefml_file data/sample.oncodrivefml.tsv \
  --depth_file data/sample.depth_per_gene.tsv \
  --output_prefix "$OUTPUT_DIR/${SAMPLE_NAME}_selection"

echo "All plots generated successfully!"
echo "Output files:"
ls -lh "$OUTPUT_DIR"/*.vaf_mutdensity_depth.pdf "$OUTPUT_DIR"/*.selection_depth.pdf
