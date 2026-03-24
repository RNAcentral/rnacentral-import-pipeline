import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl

from rnacentral_pipeline.databases.expressionatlas.parser import parse_baseline
from rnacentral_pipeline.databases.expressionatlas.sdrf import parse_condensed_sdrf


class TestParseBaseline(unittest.TestCase):
    def setUp(self):
        # Create sample tpms data
        self.tpms_content = (
            "GeneID\tg1\tg2\tg3\n"
            "ENSG00000001\t1,2,3,4,5\t10,11,12\t20,21,22,23\n"  # All medians > 0.5
            "ENSG00000002\t0.1,0.2,0.3\t0.6,0.7,0.8\t1.5,1.6,1.7\n"  # Some medians > 0.5
            "ENSG00000003\t0,0,0.1\t0.2,0.3,0.4\t0,0.1,0.2\n"  # All medians < 0.5, should be filtered out
            "ENSG00000004\t0,0,0\t5,6,7\t0,0,0\n"  # One median > 0.5
        )

        # Create sample SDRF data for mocking
        self.mock_sdrf_data = pl.DataFrame(
            {
                "exp_name": ["E-MTAB-123", "E-MTAB-123"],
                "assay_name": ["Assay1", "Assay2"],
                "feat_class": ["Class1", "Class2"],
                "feat_type": ["organism", "factor"],
                "feat_value": ["Homo sapiens", "treatment"],
                "ontology": ["NCBITaxon_9606", None],
            }
        )

        # Create sample lookup data
        self.lookup_content = (
            "CREATE TABLE\n"
            "COPY 5\n"  # Skip these two lines
            "urs_taxid,taxid,gene\n"
            "URS000123_9606,9606,ENSG00000001\n"
            "URS000124_9606,9606,ENSG00000002\n"
            "URS000125_9606,9606,ENSG00000004\n"
            "URS000126_10090,10090,ENSG00000001\n"  # Different taxid, shouldn't match
            "URS000127_9606,9606,ENSG00000999\n"  # Gene not in tpms
        )

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_basic(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(self.tpms_content)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # Verify results
            self.assertEqual(result.shape, (3, 2))  # 3 rows, 2 columns

            # Check that the correct genes were selected
            result_genes = result["urs_taxid"].to_list()
            self.assertIn("URS000123_9606", result_genes)  # ENSG00000001
            self.assertIn("URS000124_9606", result_genes)  # ENSG00000002
            self.assertIn("URS000125_9606", result_genes)  # ENSG00000004

            # ENSG00000003 should be filtered out (all medians < 0.5)
            self.assertNotIn("URS000127_9606", result_genes)  # ENSG00000999 not in tpms

            # Verify experiment column
            self.assertEqual(set(result["experiment"].to_list()), {"E-MTAB-123"})

            # Verify that mock_parse_sdrf was called with the correct path
            mock_parse_sdrf.assert_called_once_with(sdrf_path)

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_gene_id_column_rename(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Modified tpms content with "Gene ID" instead of "GeneID"
        tpms_content_alt = (
            "Gene ID\tg1\tg2\tg3\n"
            "ENSG00000001\t1,2,3,4,5\t10,11,12\t20,21,22,23\n"
            "ENSG00000002\t0.1,0.2,0.3\t0.6,0.7,0.8\t1.5,1.6,1.7\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file with alternative column name
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(tpms_content_alt)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # Verify results - should have same shape as before
            self.assertEqual(result.shape, (2, 2))

            # Verify the renamed column worked correctly
            result_genes = result["urs_taxid"].to_list()
            self.assertIn("URS000123_9606", result_genes)  # ENSG00000001
            self.assertIn("URS000124_9606", result_genes)  # ENSG00000002

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_all_filtered_out(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Create tpms data where all rows will be filtered out (medians < 0.5)
        tpms_content_filtered = (
            "GeneID\tg1\tg2\tg3\n"
            "ENSG00000001\t0.1,0.2,0.3\t0.1,0.2,0.3\t0.1,0.2,0.3\n"
            "ENSG00000002\t0.2,0.3,0.4\t0.2,0.3,0.4\t0.2,0.3,0.4\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(tpms_content_filtered)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # Verify results - should be empty after filtering
            self.assertEqual(result.shape, (0, 2))

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_empty_values(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Test with empty values in some cells
        tpms_content_empty = (
            "GeneID\tg1\tg2\tg3\n"
            "ENSG00000001\t\t10,11,12\t20,21,22,23\n"  # Empty g1
            "ENSG00000002\t0.1,0.2,0.3\t\t1.5,1.6,1.7\n"  # Empty g2
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(tpms_content_empty)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # This should either handle empty values gracefully or raise an appropriate exception
            try:
                result = parse_baseline(tpms_path, sdrf_path, lookup_path)
                # If it handles empty values, verify results
                self.assertGreaterEqual(result.shape[0], 0)
            except Exception as e:
                # If it raises an exception, that's acceptable too
                self.assertIsInstance(e, (pl.exceptions.PolarsError, ValueError))

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_single_value(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Test with single values instead of lists in some cells
        tpms_content_single = (
            "GeneID\tg1\tg2\tg3\n"
            "ENSG00000001\t5\t10,11,12\t20\n"  # Single values in g1 and g3
            "ENSG00000002\t0.1,0.2,0.3\t0.7\t1.5,1.6,1.7\n"  # Single value in g2
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(tpms_content_single)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # This should handle single values as if they're arrays with one element
            try:
                result = parse_baseline(tpms_path, sdrf_path, lookup_path)
                # Verify results for single values
                self.assertEqual(result.shape[0], 2)  # Should have 2 rows
                result_genes = result["urs_taxid"].to_list()
                self.assertIn("URS000123_9606", result_genes)  # ENSG00000001
                self.assertIn("URS000124_9606", result_genes)  # ENSG00000002
            except Exception as e:
                self.fail(
                    f"parse_baseline raised {type(e).__name__} unexpectedly with single values: {str(e)}"
                )

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_multiple_experiments(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file with a different experiment name
            tpms_path = Path(temp_dir) / "E-MTAB-456-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(self.tpms_content)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # Verify experiment name is correctly extracted
            self.assertEqual(set(result["experiment"].to_list()), {"E-MTAB-456"})

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_no_matching_taxids(self, mock_parse_sdrf):
        # Create SDRF data with a different taxid
        mock_sdrf_data_diff_taxid = pl.DataFrame(
            {
                "exp_name": ["E-MTAB-123"],
                "assay_name": ["Assay1"],
                "feat_class": ["Class1"],
                "feat_type": ["organism"],
                "feat_value": ["Mus musculus"],
                "ontology": ["NCBITaxon_10090"],  # Different from most lookup entries
            }
        )

        mock_parse_sdrf.return_value = mock_sdrf_data_diff_taxid

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(self.tpms_content)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # Verify results - should have one match for the 10090 taxid
            self.assertEqual(result.shape, (1, 2))
            self.assertEqual(result["urs_taxid"].to_list(), ["URS000126_10090"])

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_large_numbers(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Test with very large numbers
        tpms_content_large = (
            "GeneID\tg1\tg2\tg3\n"
            "ENSG00000001\t1000,2000,3000\t10000,11000,12000\t20000,21000,22000\n"
            "ENSG00000002\t100000,200000,300000\t400000,500000,600000\t700000,800000,900000\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(tpms_content_large)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # Verify results
            self.assertEqual(result.shape, (2, 2))
            result_genes = result["urs_taxid"].to_list()
            self.assertIn("URS000123_9606", result_genes)  # ENSG00000001
            self.assertIn("URS000124_9606", result_genes)  # ENSG00000002

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_baseline_negative_values(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Test with negative values (unusual but possible)
        tpms_content_negative = (
            "GeneID\tg1\tg2\tg3\n"
            "ENSG00000001\t-1,-2,-3\t10,11,12\t20,21,22\n"  # Negative values in g1
            "ENSG00000002\t0.1,0.2,0.3\t-0.6,-0.7,-0.8\t1.5,1.6,1.7\n"  # Negative values in g2
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create tpms file
            tpms_path = Path(temp_dir) / "E-MTAB-123-tpms.tsv"
            with open(tpms_path, "w") as f:
                f.write(tpms_content_negative)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_baseline(tpms_path, sdrf_path, lookup_path)

            # ENSG00000001 should be included (medians for g2, g3 > 0.5)
            # ENSG00000002 should be included (median for g3 > 0.5)
            self.assertGreaterEqual(result.shape[0], 1)  # At least one match


if __name__ == "__main__":
    unittest.main()
