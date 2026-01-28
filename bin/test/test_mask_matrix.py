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
from merge_annotation_depths import apply_mask_matrix


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


class TestApplyMaskMatrix(unittest.TestCase):
    """Test suite for apply_mask_matrix function"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        os.chdir(self.test_dir)
        
    def tearDown(self):
        """Clean up test fixtures"""
        os.chdir(self.original_dir)
        shutil.rmtree(self.test_dir)

    def create_mock_depths(self, data):
        """
        Create a mock annotated depths dataframe.
        
        Parameters
        ----------
        data : dict
            Dictionary with columns as keys and lists as values
        
        Returns
        -------
        pd.DataFrame
            Mock depths dataframe
        """
        return pd.DataFrame(data)

    def create_mock_mask_matrix(self, data, filename="mask_matrix.tsv.gz"):
        """
        Create a mock mask matrix file.
        
        Parameters
        ----------
        data : dict
            Dictionary with columns as keys and lists as values
        filename : str
            Output filename
        """
        mask_df = pd.DataFrame(data)
        mask_df.to_csv(filename, sep="\t", index=False, compression='gzip')
        return filename

    def test_basic_masking(self):
        """Test basic masking: positions with 0 in mask should have depth=0"""
        # Create mock depths
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1', 'chr1'],
            'POS': [100, 200, 300],
            'CONTEXT': ['ACA', 'TCG', 'GGC'],
            'sample1': [50, 60, 70],
            'sample2': [40, 50, 60]
        })
        
        # Create mask matrix: sample1 position 200 should be masked
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1', 'chr1'],
            'POS': [100, 200, 300],
            'sample1': [1, 0, 1],  # Position 200 masked
            'sample2': [1, 1, 1]   # All kept
        })
        
        # Apply masking
        result = apply_mask_matrix(depths, mask_file)
        
        # Verify results
        self.assertEqual(result.loc[result['POS'] == 100, 'sample1'].values[0], 50)  # Kept
        self.assertEqual(result.loc[result['POS'] == 200, 'sample1'].values[0], 0)   # Masked
        self.assertEqual(result.loc[result['POS'] == 300, 'sample1'].values[0], 70)  # Kept
        
        # sample2 should remain unchanged
        self.assertEqual(result.loc[result['POS'] == 100, 'sample2'].values[0], 40)
        self.assertEqual(result.loc[result['POS'] == 200, 'sample2'].values[0], 50)
        self.assertEqual(result.loc[result['POS'] == 300, 'sample2'].values[0], 60)

    def test_multiple_samples_masked(self):
        """Test masking across multiple samples at the same position"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'CONTEXT': ['ACA', 'TCG'],
            'sample1': [50, 60],
            'sample2': [40, 50],
            'sample3': [30, 40]
        })
        
        # Mask chr1:100 for sample1 and sample2
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'sample1': [0, 1],  # 100 masked
            'sample2': [0, 1],  # 100 masked
            'sample3': [1, 1]   # All kept
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # Position 100 should be 0 for sample1 and sample2, kept for sample3
        self.assertEqual(result.loc[result['POS'] == 100, 'sample1'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 100, 'sample2'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 100, 'sample3'].values[0], 30)
        
        # Position 200 should be kept for all
        self.assertEqual(result.loc[result['POS'] == 200, 'sample1'].values[0], 60)
        self.assertEqual(result.loc[result['POS'] == 200, 'sample2'].values[0], 50)
        self.assertEqual(result.loc[result['POS'] == 200, 'sample3'].values[0], 40)

    def test_realistic_example(self):
        """Test with realistic example from user's data"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
            'POS': [1071610, 11059842, 19851697, 23495229, 23495334, 26729792],
            'CONTEXT': ['ACA', 'TCG', 'GGC', 'TAT', 'CGC', 'AAA'],
            'P19_0001_BDO_01': [45, 50, 55, 60, 65, 70],
            'P19_0001_BTR_01': [35, 40, 45, 50, 55, 60],
            'P19_0004_BTR_01': [25, 30, 35, 40, 45, 50],
            'P19_0005_BTR_01': [15, 20, 25, 30, 35, 40]
        })
        
        # Create mask based on user's example
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
            'POS': [1071610, 11059842, 19851697, 23495229, 23495334, 26729792],
            'P19_0001_BDO_01': [1, 0, 1, 1, 1, 0],
            'P19_0001_BTR_01': [1, 0, 1, 1, 1, 1],
            'P19_0004_BTR_01': [0, 0, 1, 0, 1, 1],
            'P19_0005_BTR_01': [1, 0, 1, 1, 1, 1]
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # Check specific positions
        # chr1:1071610 - P19_0004_BTR_01 should be masked (0)
        self.assertEqual(result.loc[result['POS'] == 1071610, 'P19_0004_BTR_01'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 1071610, 'P19_0001_BDO_01'].values[0], 45)
        
        # chr1:11059842 - All samples should be masked except (P19_0001_BTR_01 and others have 0)
        self.assertEqual(result.loc[result['POS'] == 11059842, 'P19_0001_BDO_01'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 11059842, 'P19_0001_BTR_01'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 11059842, 'P19_0004_BTR_01'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 11059842, 'P19_0005_BTR_01'].values[0], 0)
        
        # chr1:23495229 - P19_0004_BTR_01 should be masked
        self.assertEqual(result.loc[result['POS'] == 23495229, 'P19_0004_BTR_01'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 23495229, 'P19_0001_BDO_01'].values[0], 60)
        
        # chr1:26729792 - P19_0001_BDO_01 should be masked
        self.assertEqual(result.loc[result['POS'] == 26729792, 'P19_0001_BDO_01'].values[0], 0)
        self.assertEqual(result.loc[result['POS'] == 26729792, 'P19_0001_BTR_01'].values[0], 60)

    def test_all_positions_kept(self):
        """Test when all positions have mask=1 (nothing masked)"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'CONTEXT': ['ACA', 'TCG'],
            'sample1': [50, 60],
            'sample2': [40, 50]
        })
        
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'sample1': [1, 1],
            'sample2': [1, 1]
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # All depths should remain unchanged
        pd.testing.assert_frame_equal(result[['CHROM', 'POS', 'sample1', 'sample2']], 
                                      depths[['CHROM', 'POS', 'sample1', 'sample2']])

    def test_all_positions_masked(self):
        """Test when all positions have mask=0 (everything masked)"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'CONTEXT': ['ACA', 'TCG'],
            'sample1': [50, 60],
            'sample2': [40, 50]
        })
        
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'sample1': [0, 0],
            'sample2': [0, 0]
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # All depths should be 0
        self.assertTrue((result['sample1'] == 0).all())
        self.assertTrue((result['sample2'] == 0).all())

    def test_empty_mask_matrix(self):
        """Test handling of empty mask matrix"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1'],
            'POS': [100],
            'CONTEXT': ['ACA'],
            'sample1': [50]
        })
        
        # Create empty mask matrix
        mask_file = self.create_mock_mask_matrix({
            'CHROM': [],
            'POS': []
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # Should return original depths unchanged
        pd.testing.assert_frame_equal(result, depths)

    def test_partial_overlap(self):
        """Test when mask matrix has only some positions from depths"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1', 'chr1'],
            'POS': [100, 200, 300],
            'CONTEXT': ['ACA', 'TCG', 'GGC'],
            'sample1': [50, 60, 70]
        })
        
        # Mask only has positions 100 and 200
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'sample1': [0, 1]
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # Position 100 should be masked
        self.assertEqual(result.loc[result['POS'] == 100, 'sample1'].values[0], 0)
        # Position 200 should be kept
        self.assertEqual(result.loc[result['POS'] == 200, 'sample1'].values[0], 60)
        # Position 300 not in mask, should remain unchanged
        self.assertEqual(result.loc[result['POS'] == 300, 'sample1'].values[0], 70)

    def test_structure_preservation(self):
        """Test that CHROM, POS, and CONTEXT columns are preserved"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr2'],
            'POS': [100, 200],
            'CONTEXT': ['ACA', 'TCG'],
            'sample1': [50, 60]
        })
        
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr2'],
            'POS': [100, 200],
            'sample1': [0, 1]
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # Check structure is preserved
        self.assertListEqual(list(result.columns), ['CHROM', 'POS', 'CONTEXT', 'sample1'])
        self.assertListEqual(result['CHROM'].tolist(), ['chr1', 'chr2'])
        self.assertListEqual(result['POS'].tolist(), [100, 200])
        self.assertListEqual(result['CONTEXT'].tolist(), ['ACA', 'TCG'])

    def test_zero_depth_remains_zero(self):
        """Test that positions already at depth=0 remain at 0 regardless of mask"""
        depths = self.create_mock_depths({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'CONTEXT': ['ACA', 'TCG'],
            'sample1': [0, 60]  # Already 0
        })
        
        mask_file = self.create_mock_mask_matrix({
            'CHROM': ['chr1', 'chr1'],
            'POS': [100, 200],
            'sample1': [1, 0]  # Keep 100, mask 200
        })
        
        result = apply_mask_matrix(depths, mask_file)
        
        # Position 100 was already 0, should remain 0 even with mask=1
        self.assertEqual(result.loc[result['POS'] == 100, 'sample1'].values[0], 0)
        # Position 200 should be masked to 0
        self.assertEqual(result.loc[result['POS'] == 200, 'sample1'].values[0], 0)

if __name__ == '__main__':
    unittest.main()
