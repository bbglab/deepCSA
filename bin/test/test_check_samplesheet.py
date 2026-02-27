#!/usr/bin/env python3
"""
Unit tests for check_samplesheet.py sample name validation

Tests the security enhancements to prevent shell injection attacks.
"""

import sys
import tempfile
import unittest
from pathlib import Path

# Add bin to path
sys.path.insert(0, str(Path(__file__).parent))
from check_samplesheet import check_samplesheet


class TestSampleNameValidation(unittest.TestCase):
    """Test suite for sample name validation and sanitization"""

    def setUp(self):
        """Set up test fixtures"""
        self.test_dir = tempfile.mkdtemp()

    def test_valid_alphanumeric(self):
        """Test that valid alphanumeric sample names pass"""
        input_file = Path(self.test_dir) / "valid1.csv"
        output_file = Path(self.test_dir) / "valid1_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample123,test.vcf.gz,test.bam\n")
        
        try:
            check_samplesheet(input_file, output_file)
            self.assertTrue(output_file.exists())
        except SystemExit:
            self.fail("Valid alphanumeric sample name was rejected")

    def test_valid_underscore(self):
        """Test that sample names with underscores pass"""
        input_file = Path(self.test_dir) / "valid2.csv"
        output_file = Path(self.test_dir) / "valid2_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample_test_123,test.vcf.gz,test.bam\n")
        
        try:
            check_samplesheet(input_file, output_file)
            self.assertTrue(output_file.exists())
        except SystemExit:
            self.fail("Sample name with underscores was rejected")

    def test_valid_hyphen(self):
        """Test that sample names with hyphens pass"""
        input_file = Path(self.test_dir) / "valid3.csv"
        output_file = Path(self.test_dir) / "valid3_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample-test-123,test.vcf.gz,test.bam\n")
        
        try:
            check_samplesheet(input_file, output_file)
            self.assertTrue(output_file.exists())
        except SystemExit:
            self.fail("Sample name with hyphens was rejected")

    def test_valid_dots(self):
        """Test that sample names with dots pass"""
        input_file = Path(self.test_dir) / "valid4.csv"
        output_file = Path(self.test_dir) / "valid4_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample.test.123,test.vcf.gz,test.bam\n")
        
        try:
            check_samplesheet(input_file, output_file)
            self.assertTrue(output_file.exists())
        except SystemExit:
            self.fail("Sample name with dots was rejected")

    def test_space_replacement(self):
        """Test that spaces are replaced with underscores"""
        input_file = Path(self.test_dir) / "space.csv"
        output_file = Path(self.test_dir) / "space_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample test,test.vcf.gz,test.bam\n")
        
        try:
            check_samplesheet(input_file, output_file)
            self.assertTrue(output_file.exists())
            # Verify space was replaced
            with open(output_file, 'r') as f:
                content = f.read()
                self.assertIn("sample_test", content)
                self.assertNotIn("sample test", content)
        except SystemExit:
            self.fail("Sample name with space was not properly sanitized")

    def test_reject_semicolon(self):
        """Test that semicolons are rejected (shell command separator)"""
        input_file = Path(self.test_dir) / "invalid1.csv"
        output_file = Path(self.test_dir) / "invalid1_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample;rm -rf /,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_dollar(self):
        """Test that dollar signs are rejected (command substitution)"""
        input_file = Path(self.test_dir) / "invalid2.csv"
        output_file = Path(self.test_dir) / "invalid2_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample$(whoami),test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_backtick(self):
        """Test that backticks are rejected (command substitution)"""
        input_file = Path(self.test_dir) / "invalid3.csv"
        output_file = Path(self.test_dir) / "invalid3_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample`whoami`,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_pipe(self):
        """Test that pipes are rejected (command chaining)"""
        input_file = Path(self.test_dir) / "invalid4.csv"
        output_file = Path(self.test_dir) / "invalid4_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample|cat /etc/passwd,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_ampersand(self):
        """Test that ampersands are rejected (background execution)"""
        input_file = Path(self.test_dir) / "invalid5.csv"
        output_file = Path(self.test_dir) / "invalid5_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample&whoami,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_redirect(self):
        """Test that redirect operators are rejected"""
        input_file = Path(self.test_dir) / "invalid6.csv"
        output_file = Path(self.test_dir) / "invalid6_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample>output.txt,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_parentheses(self):
        """Test that parentheses are rejected (subshell)"""
        input_file = Path(self.test_dir) / "invalid7.csv"
        output_file = Path(self.test_dir) / "invalid7_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample(test),test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_hyphen_start(self):
        """Test that sample names starting with hyphen are rejected"""
        input_file = Path(self.test_dir) / "invalid8.csv"
        output_file = Path(self.test_dir) / "invalid8_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("-sample1,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_slash(self):
        """Test that slashes are rejected (path traversal)"""
        input_file = Path(self.test_dir) / "invalid9.csv"
        output_file = Path(self.test_dir) / "invalid9_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write("sample/../../etc/passwd,test.vcf.gz,test.bam\n")
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)

    def test_reject_quotes(self):
        """Test that quotes are rejected"""
        input_file = Path(self.test_dir) / "invalid10.csv"
        output_file = Path(self.test_dir) / "invalid10_out.csv"
        
        with open(input_file, 'w') as f:
            f.write("sample,vcf,bam\n")
            f.write('sample"test,test.vcf.gz,test.bam\n')
        
        with self.assertRaises(SystemExit):
            check_samplesheet(input_file, output_file)


def run_tests():
    """Run all tests and return results"""
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    suite.addTests(loader.loadTestsFromTestCase(TestSampleNameValidation))
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
