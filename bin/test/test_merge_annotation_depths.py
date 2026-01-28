#!/usr/bin/env python3
"""
Unit tests for merge_annotation_depths.py

Tests the functionality of applying mask matrices to depth data.
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
from merge_annotation_depths import apply_mask_matrix


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
