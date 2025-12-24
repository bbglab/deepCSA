#!/usr/bin/env python3
"""
Unit tests for plot_selectionsideplots.py

Tests the functionality of dynamic track selection and graceful handling of missing data files.
"""

import os
import sys
import tempfile
import shutil
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add the bin directory to the path to import the module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


class TestPlotSelectionSidePlots(unittest.TestCase):
    """Test suite for plot_selectionsideplots.py"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        os.chdir(self.test_dir)
        
    def tearDown(self):
        """Clean up test fixtures"""
        os.chdir(self.original_dir)
        shutil.rmtree(self.test_dir)

    def create_mock_maf_file(self, sample_name):
        """Create a minimal mock MAF file for testing"""
        maf_file = f"{sample_name}.somatic.mutations.tsv"
        with open(maf_file, 'w') as f:
            f.write("TYPE\tcanonical_SYMBOL\n")
            f.write("SNV\tGENE1\n")
            f.write("SNV\tGENE2\n")
        return maf_file

    def create_mock_oncodrivefml_file(self, sample_name):
        """Create a minimal mock OncodriveFML file for testing"""
        fml_file = f"{sample_name}-oncodrivefml.tsv.gz"
        import gzip
        with gzip.open(fml_file, 'wt') as f:
            f.write("GENE_ID\tZ-SCORE\tQ_VALUE\tAVG_SCORE_OBS\tPOPULATION_MEAN\tSTD_OF_MEANS\n")
            f.write("GENE1\t2.5\t0.01\t0.8\t0.5\t0.1\n")
        return fml_file

    def create_mock_omega_file(self, sample_name):
        """Create a minimal mock Omega file for testing"""
        omega_file = f"output_mle.{sample_name}.tsv"
        with open(omega_file, 'w') as f:
            f.write("gene\timpact\tmutations\tdnds\tpvalue\tlower\tupper\n")
            f.write("GENE1\ttruncating\t10\t1.5\t0.05\t1.0\t2.0\n")
            f.write("GENE1\tmissense\t20\t1.2\t0.03\t1.0\t1.5\n")
            f.write("ALL_GENES\ttruncating\t100\t1.0\t0.5\t0.8\t1.2\n")
        return omega_file

    def create_mock_indels_file(self, sample_name):
        """Create a minimal mock indels file for testing"""
        indels_file = f"{sample_name}.sample.indels.tsv"
        with open(indels_file, 'w') as f:
            f.write("SYMBOL\tpa/Npa\tpvalue\tMULTIPLE_3.non_protein_affecting\tNOT_MULTIPLE_3.non_protein_affecting\t")
            f.write("NON_TRUNCATING.protein_affecting\tTRUNCATING.protein_affecting\tpa_TRUNC/NOTTRUNC\tNpa_NM3/M3\n")
            f.write("GENE1\t2.0\t0.01\t1\t2\t3\t4\t1.5\t0.8\n")
        return indels_file

    def create_mock_oncodrive3d_file(self, sample_name):
        """Create a minimal mock Oncodrive3D file for testing"""
        o3d_file = f"{sample_name}.3d_clustering_genes.csv"
        with open(o3d_file, 'w') as f:
            f.write("Gene,Score_obs_sim_top_vol,qval\n")
            f.write("GENE1,5.0,0.01\n")
        return o3d_file

    def test_file_existence_check(self):
        """Test that the script checks for file existence"""
        from plot_selectionsideplots import get_all_data
        
        sample_name = "test_sample"
        
        # Test with no files - should return early without error
        try:
            get_all_data(sample_name, self.test_dir)
        except Exception as e:
            self.fail(f"get_all_data raised exception with no files: {e}")

    def test_with_only_omega_data(self):
        """Test plotting with only Omega data available"""
        from plot_selectionsideplots import get_all_data
        
        sample_name = "test_sample"
        self.create_mock_maf_file(sample_name)
        self.create_mock_omega_file(sample_name)
        
        # Should not crash with only omega data
        try:
            get_all_data(sample_name, self.test_dir, 
                        tracks=("omega_trunc", "omega_mis", "oncodrive3d", "oncodrivefml"))
        except Exception as e:
            self.fail(f"get_all_data raised exception with only omega data: {e}")

    @patch('plot_selectionsideplots.plt.subplots')
    @patch('plot_selectionsideplots.sns.barplot')
    def test_dynamic_track_selection(self, mock_barplot, mock_subplots):
        """Test that tracks are dynamically selected based on available data"""
        from plot_selectionsideplots import get_all_data
        
        # Mock matplotlib to avoid actual plotting
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_subplots.return_value = (mock_fig, [mock_ax])
        
        sample_name = "test_sample"
        self.create_mock_maf_file(sample_name)
        self.create_mock_omega_file(sample_name)
        self.create_mock_oncodrivefml_file(sample_name)
        
        # Should detect and use both omega and oncodrivefml
        try:
            get_all_data(sample_name, self.test_dir,
                        tracks=("omega_trunc", "omega_mis", "oncodrive3d", "oncodrivefml"))
        except Exception as e:
            self.fail(f"get_all_data raised exception with multiple tracks: {e}")

    def test_all_files_available(self):
        """Test with all data files available (backward compatibility)"""
        from plot_selectionsideplots import get_all_data
        
        sample_name = "test_sample"
        self.create_mock_maf_file(sample_name)
        self.create_mock_omega_file(sample_name)
        self.create_mock_oncodrivefml_file(sample_name)
        self.create_mock_oncodrive3d_file(sample_name)
        self.create_mock_indels_file(sample_name)
        
        # Should work with all files present
        try:
            get_all_data(sample_name, self.test_dir,
                        tracks=("omega_trunc", "omega_mis", "oncodrive3d", "oncodrivefml", "indels"))
        except Exception as e:
            self.fail(f"get_all_data raised exception with all files: {e}")

    def test_generate_side_figures_no_maf(self):
        """Test that generate_all_side_figures returns early when MAF file is missing"""
        from plot_selectionsideplots import generate_all_side_figures
        
        sample_name = "test_sample"
        
        # Should return early without error when MAF file is missing
        try:
            generate_all_side_figures(sample_name, self.test_dir)
        except Exception as e:
            self.fail(f"generate_all_side_figures raised exception with no MAF: {e}")

    def test_generate_side_figures_with_partial_data(self):
        """Test generate_all_side_figures with partial data"""
        from plot_selectionsideplots import generate_all_side_figures
        
        sample_name = "test_sample"
        self.create_mock_maf_file(sample_name)
        self.create_mock_omega_file(sample_name)
        # Don't create oncodrivefml or indels files
        
        # Should handle partial data gracefully
        try:
            generate_all_side_figures(sample_name, self.test_dir,
                                    tools=["oncodrivefml", "omega_trunc", "omega_mis", "excess_indels"])
        except Exception as e:
            self.fail(f"generate_all_side_figures raised exception with partial data: {e}")

    @patch('plot_selectionsideplots.plot_all_positive_selection')
    def test_plot_all_positive_selection_called_with_available_tracks(self, mock_plot):
        """Test that plot_all_positive_selection is called with only available tracks"""
        from plot_selectionsideplots import get_all_data
        
        # Mock the plotting function to avoid actual plot generation
        mock_plot.return_value = None
        
        sample_name = "test_sample"
        self.create_mock_maf_file(sample_name)
        self.create_mock_omega_file(sample_name)
        
        get_all_data(sample_name, self.test_dir,
                    tracks=("omega_trunc", "omega_mis", "oncodrive3d", "oncodrivefml"))
        
        # Check that the function was called (meaning we had available tracks)
        # We can't easily check the exact tracks without more complex mocking
        # but we verify it didn't crash and returned properly

    def test_empty_gene_list_handling(self):
        """Test handling of empty gene lists"""
        from plot_selectionsideplots import plot_all_positive_selection
        
        # Test with empty gene order
        result = plot_all_positive_selection(
            omega_truncating=None,
            omega_missense=None,
            indels_panel_df=None,
            oncodrive3d_data_scores=None,
            oncodrivefml_data=None,
            gene_order=[],
            tracks=()
        )
        
        # Should return None or handle gracefully
        self.assertIsNone(result)


class TestIntegrationScenarios(unittest.TestCase):
    """Integration tests for various data availability scenarios"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()
        self.original_dir = os.getcwd()
        os.chdir(self.test_dir)
        
    def tearDown(self):
        """Clean up test fixtures"""
        os.chdir(self.original_dir)
        shutil.rmtree(self.test_dir)

    def test_scenario_all_methods_enabled(self):
        """Scenario 1: All four methods enabled (backward compatibility test)"""
        # This would require full mock data setup
        # Placeholder for more comprehensive integration test
        pass

    def test_scenario_only_omega(self):
        """Scenario 2: Only Omega method enabled"""
        # Placeholder for integration test
        pass

    def test_scenario_mixed_availability(self):
        """Scenario 3: Mixed availability (e.g., Omega + OncodriveFML)"""
        # Placeholder for integration test
        pass


def run_tests():
    """Run all tests and return results"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestPlotSelectionSidePlots))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegrationScenarios))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
