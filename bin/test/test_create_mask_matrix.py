#!/usr/bin/env python3
"""
Unit tests for create_mask_matrix.py

Tests the functionality of creating position × sample mask matrices from BED files.
"""

import os
import sys
import tempfile
import shutil
import unittest
from pathlib import Path
import pandas as pd
import gzip

# Add the bin directory to the path to import the module
sys.path.insert(0, str(Path(__file__).parent.parent))
from create_mask_matrix import create_mask_matrix


class TestCreateMaskMatrix(unittest.TestCase):
    """Test suite for create_mask_matrix.py"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        os.chdir(self.test_dir)
        
    def tearDown(self):
        """Clean up test fixtures"""
        os.chdir(self.original_dir)
        shutil.rmtree(self.test_dir)

    def create_mock_bed_file(self, sample_name, positions):
        """
        Create a mock BED file with specified positions.
        
        Parameters
        ----------
        sample_name : str
            Name of the sample
        positions : list of tuples
            List of (chrom, start, end, filter) tuples
        """
        bed_file = f"{sample_name}.flagged-pos.bed"
        with open(bed_file, 'w') as f:
            for chrom, start, end, filter_val in positions:
                f.write(f"{chrom}\t{start}\t{end}\t{filter_val}\n")
        return bed_file

    def read_mask_matrix(self, filename="flagged_positions.mask.tsv.gz"):
        """Read and return the mask matrix"""
        return pd.read_csv(filename, sep="\t", compression='gzip')

    def test_single_sample_single_position(self):
        """Test creating mask matrix with one sample and one position"""
        bed_files = [
            self.create_mock_bed_file("sample1", [("chr1", 100, 100, "n_rich")])
        ]
        
        create_mask_matrix(bed_files)
        
        # Check output exists
        self.assertTrue(Path("flagged_positions.mask.tsv.gz").exists())
        
        # Check content
        mask_df = self.read_mask_matrix()
        self.assertEqual(len(mask_df), 1)
        self.assertEqual(mask_df.iloc[0]['CHROM'], 'chr1')
        self.assertEqual(mask_df.iloc[0]['POS'], 100)
        self.assertEqual(mask_df.iloc[0]['sample1'], 0)  # Should be masked

    def test_multiple_samples(self):
        """Test creating mask matrix with multiple samples"""
        bed_files = [
            self.create_mock_bed_file("sample1", [("chr1", 100, 100, "n_rich")]),
            self.create_mock_bed_file("sample2", [("chr1", 200, 200, "nanoseq_snp")]),
            self.create_mock_bed_file("sample3", [("chr2", 300, 300, "nanoseq_noise")])
        ]
        
        create_mask_matrix(bed_files)
        
        mask_df = self.read_mask_matrix()
        
        # Should have 3 positions total
        self.assertEqual(len(mask_df), 3)
        
        # Check that each sample has correct columns
        self.assertIn('sample1', mask_df.columns)
        self.assertIn('sample2', mask_df.columns)
        self.assertIn('sample3', mask_df.columns)
        
        # Check masking pattern: position should be 0 for its sample, 1 for others
        row1 = mask_df[(mask_df['CHROM'] == 'chr1') & (mask_df['POS'] == 100)].iloc[0]
        self.assertEqual(row1['sample1'], 0)  # Masked
        self.assertEqual(row1['sample2'], 1)  # Not masked
        self.assertEqual(row1['sample3'], 1)  # Not masked

    def test_overlapping_positions(self):
        """Test when multiple samples have the same position flagged"""
        bed_files = [
            self.create_mock_bed_file("sample1", [("chr1", 100, 100, "n_rich")]),
            self.create_mock_bed_file("sample2", [("chr1", 100, 100, "nanoseq_snp")]),
        ]
        
        create_mask_matrix(bed_files)
        
        mask_df = self.read_mask_matrix()
        
        # Should have only 1 position
        self.assertEqual(len(mask_df), 1)
        
        # Both samples should have 0 for this position
        row = mask_df.iloc[0]
        self.assertEqual(row['CHROM'], 'chr1')
        self.assertEqual(row['POS'], 100)
        self.assertEqual(row['sample1'], 0)
        self.assertEqual(row['sample2'], 0)

    def test_range_expansion(self):
        """Test that ranges are expanded to individual positions (1-based inclusive)"""
        bed_files = [
            self.create_mock_bed_file("sample1", [("chr1", 100, 102, "nanoseq_noise")])
        ]
        
        create_mask_matrix(bed_files)
        
        mask_df = self.read_mask_matrix()
        
        # Should have 3 positions: 100, 101, 102 (1-based inclusive)
        self.assertEqual(len(mask_df), 3)
        positions = sorted(mask_df['POS'].tolist())
        self.assertEqual(positions, [100, 101, 102])
        
        # All should be masked for sample1
        for pos in [100, 101, 102]:
            row = mask_df[mask_df['POS'] == pos].iloc[0]
            self.assertEqual(row['sample1'], 0)

    def test_empty_bed_file(self):
        """Test handling of empty BED file"""
        # Create an empty BED file
        bed_file = "sample1.flagged-pos.bed"
        with open(bed_file, 'w') as f:
            pass  # Empty file
        
        bed_files = [bed_file]
        create_mask_matrix(bed_files)
        
        # Should create empty matrix with just CHROM and POS columns
        mask_df = self.read_mask_matrix()
        self.assertEqual(len(mask_df), 0)
        self.assertIn('CHROM', mask_df.columns)
        self.assertIn('POS', mask_df.columns)

    def test_multiple_chromosomes(self):
        """Test handling of multiple chromosomes"""
        bed_files = [
            self.create_mock_bed_file("sample1", [
                ("chr1", 100, 100, "n_rich"),
                ("chr2", 200, 200, "nanoseq_snp"),
                ("chrX", 300, 300, "nanoseq_noise")
            ])
        ]
        
        create_mask_matrix(bed_files)
        
        mask_df = self.read_mask_matrix()
        
        # Should have 3 positions on different chromosomes
        self.assertEqual(len(mask_df), 3)
        chroms = sorted(mask_df['CHROM'].unique())
        self.assertIn('chr1', chroms)
        self.assertIn('chr2', chroms)
        self.assertIn('chrX', chroms)

    def test_sorting(self):
        """Test that output is sorted by CHROM and POS"""
        bed_files = [
            self.create_mock_bed_file("sample1", [
                ("chr2", 200, 200, "nanoseq_noise"),
                ("chr1", 150, 150, "nanoseq_snp"),
                ("chr1", 100, 100, "n_rich"),
            ])
        ]
        
        create_mask_matrix(bed_files)
        
        mask_df = self.read_mask_matrix()
        
        # Check sorting
        chroms = mask_df['CHROM'].tolist()
        positions = mask_df['POS'].tolist()
        
        # chr1 should come before chr2
        chr1_indices = [i for i, c in enumerate(chroms) if c == 'chr1']
        chr2_indices = [i for i, c in enumerate(chroms) if c == 'chr2']
        
        if chr1_indices and chr2_indices:
            self.assertTrue(max(chr1_indices) < min(chr2_indices))
        
        # Within chr1, positions should be sorted
        chr1_positions = [positions[i] for i in chr1_indices]
        self.assertEqual(chr1_positions, sorted(chr1_positions))

    def test_complex_scenario(self):
        """Test a complex scenario with multiple samples, overlaps, and ranges"""
        bed_files = [
            self.create_mock_bed_file("sample1", [
                ("chr1", 100, 102, "n_rich"),  # 3 positions
                ("chr2", 200, 200, "nanoseq_snp"),  # 1 position
            ]),
            self.create_mock_bed_file("sample2", [
                ("chr1", 101, 101, "nanoseq_noise"),  # Overlaps with sample1
                ("chr3", 300, 301, "n_rich"),  # 2 positions
            ]),
            self.create_mock_bed_file("sample3", [
                ("chr2", 200, 200, "nanoseq_snp"),  # Overlaps with sample1
            ])
        ]
        
        create_mask_matrix(bed_files)
        
        mask_df = self.read_mask_matrix()
        
        # Total unique positions: chr1:100,101,102 + chr2:200 + chr3:300,301 = 6
        self.assertEqual(len(mask_df), 6)
        
        # Check specific position: chr1:101 (masked in sample1 and sample2)
        row_101 = mask_df[(mask_df['CHROM'] == 'chr1') & (mask_df['POS'] == 101)].iloc[0]
        self.assertEqual(row_101['sample1'], 0)
        self.assertEqual(row_101['sample2'], 0)
        self.assertEqual(row_101['sample3'], 1)
        
        # Check specific position: chr2:200 (masked in sample1 and sample3)
        row_200 = mask_df[(mask_df['CHROM'] == 'chr2') & (mask_df['POS'] == 200)].iloc[0]
        self.assertEqual(row_200['sample1'], 0)
        self.assertEqual(row_200['sample2'], 1)
        self.assertEqual(row_200['sample3'], 0)

    def test_sample_name_extraction(self):
        """Test that sample names are correctly extracted from filenames"""
        # Create BED file with complex name
        bed_file = "complex_sample_name-123.flagged-pos.bed"
        with open(bed_file, 'w') as f:
            f.write("chr1\t100\t100\tn_rich\n")
        
        create_mask_matrix([bed_file])
        
        mask_df = self.read_mask_matrix()
        
        # Sample name should be extracted correctly
        self.assertIn('complex_sample_name-123', mask_df.columns)

    def test_all_samples_empty(self):
        """Test when all samples have empty BED files"""
        bed_files = [
            self.create_mock_bed_file("sample1", []),
            self.create_mock_bed_file("sample2", []),
        ]
        
        create_mask_matrix(bed_files)
        
        # Should create empty matrix
        mask_df = self.read_mask_matrix()
        self.assertEqual(len(mask_df), 0)


if __name__ == '__main__':
    unittest.main()
