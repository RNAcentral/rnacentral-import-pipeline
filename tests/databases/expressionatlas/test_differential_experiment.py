import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import polars as pl
import pytest

from rnacentral_pipeline.databases.expressionatlas.parser import parse_differential
from rnacentral_pipeline.databases.expressionatlas.sdrf import parse_condensed_sdrf


class TestParseDifferential(unittest.TestCase):
    def setUp(self):
        # Create sample analytics data
        self.analytics_content = (
            "GeneID\tlog2foldchange\tp-value\n"
            "ENSG00000001\t2.5\t0.01\n"
            "ENSG00000002\t1.2\t0.04\n"
            "ENSG00000003\t0.5\t0.02\n"  # Will be filtered out (log2fc < 1)
            "ENSG00000004\t-1.5\t0.03\n"
            "ENSG00000005\t3.0\tNA\n"  # Will be filtered out (null p-value)
            "ENSG00000006\t1.8\t0.06\n"  # Will be filtered out (p-value > 0.05)
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
            "URS000127_9606,9606,ENSG00000999\n"  # Gene not in analytics
        )

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_basic(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create analytics file
            analytics_path = Path(temp_dir) / "E-MTAB-123-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(self.analytics_content)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_differential(analytics_path, sdrf_path, lookup_path)
            # Verify results
            self.assertEqual(result.shape, (3, 2))  # 2 rows, 2 columns

            # Check that the correct genes were selected
            result_genes = result["urs_taxid"].to_list()
            self.assertIn("URS000123_9606", result_genes)  # ENSG00000001
            self.assertIn("URS000125_9606", result_genes)  # ENSG00000004

            # ENSG00000002 should be included (it passes filters)
            self.assertIn("URS000124_9606", result_genes)

            # Verify experiment column
            self.assertEqual(set(result["experiment"].to_list()), {"E-MTAB-123"})

            # Verify that mock_parse_sdrf was called with the correct path
            mock_parse_sdrf.assert_called_once_with(sdrf_path)

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_gene_id_column_rename(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Modified analytics content with "Gene ID" instead of "GeneID"
        analytics_content_alt = (
            "Gene ID\tlog2foldchange\tp-value\n"
            "ENSG00000001\t2.5\t0.01\n"
            "ENSG00000002\t1.2\t0.04\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create analytics file with alternative column name
            analytics_path = Path(temp_dir) / "E-MTAB-123-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(analytics_content_alt)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_differential(analytics_path, sdrf_path, lookup_path)

            # Verify results - should have same shape as before
            self.assertEqual(result.shape, (2, 2))

            # Verify the renamed column worked correctly
            result_genes = result["urs_taxid"].to_list()
            self.assertIn("URS000123_9606", result_genes)  # ENSG00000001
            self.assertIn("URS000124_9606", result_genes)  # ENSG00000002

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_no_matching_genes(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Create analytics with genes that don't match the lookup
        analytics_content_no_match = (
            "GeneID\tlog2foldchange\tp-value\n"
            "ENSG99999001\t2.5\t0.01\n"
            "ENSG99999002\t1.2\t0.04\n"
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create analytics file
            analytics_path = Path(temp_dir) / "E-MTAB-123-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(analytics_content_no_match)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_differential(analytics_path, sdrf_path, lookup_path)

            # Verify results - should be empty
            self.assertEqual(result.shape, (0, 2))

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_no_matching_taxids(self, mock_parse_sdrf):
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
            # Create analytics file
            analytics_path = Path(temp_dir) / "E-MTAB-123-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(self.analytics_content)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_differential(analytics_path, sdrf_path, lookup_path)

            # Verify results - should have one match for the 10090 taxid
            self.assertEqual(result.shape, (1, 2))
            self.assertEqual(result["urs_taxid"].to_list(), ["URS000126_10090"])

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_multiple_experiments(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create analytics file with a different experiment name
            analytics_path = Path(temp_dir) / "E-MTAB-456-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(self.analytics_content)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_differential(analytics_path, sdrf_path, lookup_path)

            # Verify experiment name is correctly extracted
            self.assertEqual(set(result["experiment"].to_list()), {"E-MTAB-456"})

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_empty_analytics(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Create empty analytics data (just header)
        analytics_content_empty = "GeneID\tlog2foldchange\tp-value\n"

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create analytics file
            analytics_path = Path(temp_dir) / "E-MTAB-123-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(analytics_content_empty)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            with pytest.raises(ValueError):
                result = parse_differential(analytics_path, sdrf_path, lookup_path)

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_parse_differential_all_filtered_out(self, mock_parse_sdrf):
        mock_parse_sdrf.return_value = self.mock_sdrf_data

        # Create analytics data where all rows will be filtered out
        analytics_content_filtered = (
            "GeneID\tlog2foldchange\tp-value\n"
            "ENSG00000001\t0.5\t0.01\n"  # log2fc < 1
            "ENSG00000002\t1.2\t0.06\n"  # p-value > 0.05
            "ENSG00000003\t2.0\tNA\n"  # null p-value
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create analytics file
            analytics_path = Path(temp_dir) / "E-MTAB-123-analytics.tsv"
            with open(analytics_path, "w") as f:
                f.write(analytics_content_filtered)

            # Create lookup file
            lookup_path = Path(temp_dir) / "lookup.csv"
            with open(lookup_path, "w") as f:
                f.write(self.lookup_content)

            # Create mock SDRF path
            sdrf_path = Path(temp_dir) / "condensed-sdrf.tsv"

            # Run the function
            result = parse_differential(analytics_path, sdrf_path, lookup_path)

            # Verify results - should be empty after filtering
            self.assertEqual(result.shape, (0, 2))


if __name__ == "__main__":
    unittest.main()
