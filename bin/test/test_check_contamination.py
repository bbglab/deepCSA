#!/usr/bin/env python3
"""
Unit tests for check_contamination.py

Tests the functionality of checking contamination in genomic data.
"""

import sys
import unittest
from pathlib import Path
import pandas as pd
import numpy as np

# Add the bin directory to the path to import the module
sys.path.insert(0, str(Path(__file__).parent.parent))
from check_contamination import compute_shared_variants

class TestComputeSharedVariants(unittest.TestCase):
    """Test suite for compute_shared_variants function."""

    def _somatic_df(self):
        """Return a minimal somatic variants DataFrame."""
        return pd.DataFrame({
            'SAMPLE_ID': ['s1', 's1', 's2'],
            'MUT_ID':   ['mut1', 'mut2', 'mut3'],
        })

    def _germline_df(self):
        """Return a minimal germline variants DataFrame."""
        return pd.DataFrame({
            'SAMPLE_ID': ['s1', 's1', 's2'],
            'MUT_ID':   ['mut1', 'mut4', 'mut5'],
        })

    # ---- basic behaviour ------------------------------------------------

    def test_returns_dataframe(self):
        result = compute_shared_variants(self._somatic_df(), self._germline_df())
        self.assertIsInstance(result, pd.DataFrame)

    def test_dtype_is_int(self):
        result = compute_shared_variants(self._somatic_df(), self._germline_df())
        self.assertEqual(result.dtypes['s1'], np.int64)

    def test_shared_mutations_count(self):
        """s1 shares mut1 with s1 → cell [s1, s1] == 1."""
        result = compute_shared_variants(self._somatic_df(), self._germline_df())
        self.assertEqual(result.loc['s1', 's1'], 1)

    def test_no_shared_mutations(self):
        """s2 has mut3 which s1/s2 don't have → cell [s2, s1] == 0."""
        result = compute_shared_variants(self._somatic_df(), self._germline_df())
        self.assertEqual(result.loc['s2', 's1'], 0)

    def test_multiple_shared_mutations(self):
        som = pd.DataFrame({'SAMPLE_ID': ['s1'], 'MUT_ID': ['mut1']})
        ger = pd.DataFrame({
            'SAMPLE_ID': ['s1', 's1', 's1'],
            'MUT_ID':    ['mut1', 'mut1', 'mut2'],   # mut1 appears twice
        })
        result = compute_shared_variants(som, ger)
        # mut1 is in both sets → intersection size is 1 (set semantics)
        self.assertEqual(result.loc['s1', 's1'], 1)

    def test_empty_somatic(self):
        result = compute_shared_variants(
            pd.DataFrame({'SAMPLE_ID': [], 'MUT_ID': []}),
            self._germline_df(),
        )
        self.assertEqual(len(result), 0)

    def test_empty_germline(self):
        result = compute_shared_variants(
            self._somatic_df(),
            pd.DataFrame({'SAMPLE_ID': [], 'MUT_ID': []}),
        )
        self.assertEqual(len(result.columns), 0)

    def test_both_empty(self):
        result = compute_shared_variants(
            pd.DataFrame({'SAMPLE_ID': [], 'MUT_ID': []}),
            pd.DataFrame({'SAMPLE_ID': [], 'MUT_ID': []}),
        )
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (0, 0))

    def test_samples_sorted(self):
        """Sample IDs should appear in sorted order in index/columns."""
        som = pd.DataFrame({'SAMPLE_ID': ['z', 'a'], 'MUT_ID': ['m1', 'm2']})
        ger = pd.DataFrame({'SAMPLE_ID': ['y', 'b'], 'MUT_ID': ['m3', 'm4']})
        result = compute_shared_variants(som, ger)
        self.assertEqual(list(result.index), ['a', 'z'])
        self.assertEqual(list(result.columns), ['b', 'y'])
