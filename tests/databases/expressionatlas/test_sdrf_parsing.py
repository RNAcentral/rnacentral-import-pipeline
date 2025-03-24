import io
import logging
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import mock_open, patch

import polars as pl
import pytest

# Import the module containing the parse_condensed_sdrf function
from rnacentral_pipeline.databases.expressionatlas.sdrf import parse_condensed_sdrf


class TestParseCondensedSdrf(unittest.TestCase):
    def setUp(self):
        # Set up logging to capture log messages
        self.log_capture = io.StringIO()
        self.log_handler = logging.StreamHandler(self.log_capture)
        logging.getLogger().addHandler(self.log_handler)
        logging.getLogger().setLevel(logging.INFO)

    def tearDown(self):
        # Remove the log handler after each test
        logging.getLogger().removeHandler(self.log_handler)

    def test_standard_seven_column_format(self):
        # Test with standard 7-column format
        content = (
            "E-MTAB-123\t\tAssay1\tClass1\torganism\tHomo sapiens\tNCBITaxon_9606\n"
            "E-MTAB-123\t\tAssay2\tClass2\torganism\tMus musculus\tNCBITaxon_10090\n"
            "E-MTAB-123\t\tAssay3\tClass3\tfactor\ttreatment\t"
        )

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_file.write(content)
            temp_path = temp_file.name

        try:
            result = parse_condensed_sdrf(temp_path)

            # Check if the DataFrame has the expected structure
            self.assertEqual(result.shape, (3, 6))
            self.assertEqual(
                list(result.columns),
                [
                    "exp_name",
                    "assay_name",
                    "feat_class",
                    "feat_type",
                    "feat_value",
                    "ontology",
                ],
            )

            # Check specific values
            self.assertEqual(
                result["exp_name"].to_list(), ["E-MTAB-123", "E-MTAB-123", "E-MTAB-123"]
            )
            self.assertEqual(
                result["assay_name"].to_list(), ["Assay1", "Assay2", "Assay3"]
            )
            self.assertEqual(
                result["feat_class"].to_list(), ["Class1", "Class2", "Class3"]
            )
            self.assertEqual(
                result["feat_type"].to_list(), ["organism", "organism", "factor"]
            )
            self.assertEqual(
                result["feat_value"].to_list(),
                ["Homo sapiens", "Mus musculus", "treatment"],
            )
            self.assertEqual(
                result["ontology"].to_list(),
                ["NCBITaxon_9606", "NCBITaxon_10090", None],
            )

            # Check log message
            self.assertIn(
                f"Loading sdrf data from {temp_path}", self.log_capture.getvalue()
            )
        finally:
            os.unlink(temp_path)

    def test_unusual_six_column_format(self):
        # Test with unusual 6-column format
        content = (
            "E-MTAB-456\tAssay1\tClass1\torganism\tHomo sapiens\tNCBITaxon_9606\n"
            "E-MTAB-456\tAssay2\tClass2\torganism\tMus musculus\tNCBITaxon_10090\n"
            "E-MTAB-456\tAssay3\tClass3\tfactor\ttreatment\t"
        )

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_file.write(content)
            temp_path = temp_file.name

        try:
            result = parse_condensed_sdrf(temp_path)

            # Check if the DataFrame has the expected structure
            self.assertEqual(result.shape, (3, 6))

            # Check specific values for the unusual format
            self.assertEqual(
                result["exp_name"].to_list(), ["E-MTAB-456", "E-MTAB-456", "E-MTAB-456"]
            )
            self.assertEqual(
                result["assay_name"].to_list(), ["Assay1", "Assay2", "Assay3"]
            )
            self.assertEqual(
                result["feat_class"].to_list(), ["Class1", "Class2", "Class3"]
            )
            self.assertEqual(
                result["feat_type"].to_list(), ["organism", "organism", "factor"]
            )
            self.assertEqual(
                result["feat_value"].to_list(),
                ["Homo sapiens", "Mus musculus", "treatment"],
            )
            self.assertEqual(
                result["ontology"].to_list(),
                ["NCBITaxon_9606", "NCBITaxon_10090", None],
            )

            # Check warning log message for unusual format
            self.assertIn(
                "Unusual sdrf parsing with 6 columns, not 7 for experiment E-MTAB-456",
                self.log_capture.getvalue(),
            )
        finally:
            os.unlink(temp_path)

    def test_mixed_row_length(self):
        # Test with mixed row lengths (some with ontology, some without)
        content = (
            "E-MTAB-789\t\tAssay1\tClass1\torganism\tHomo sapiens\tNCBITaxon_9606\n"
            "E-MTAB-789\t\tAssay2\tClass2\torganism\tMus musculus\n"  # Missing ontology
            "E-MTAB-789\t\tAssay3\tClass3\tfactor\ttreatment\tONTO_12345"
        )

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_file.write(content)
            temp_path = temp_file.name

        try:
            result = parse_condensed_sdrf(temp_path)

            # Check the ontology column which should handle missing values
            self.assertEqual(
                result["ontology"].to_list(), ["NCBITaxon_9606", None, "ONTO_12345"]
            )
        finally:
            os.unlink(temp_path)

    def test_empty_file(self):
        # Test with an empty file
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_path = temp_file.name

        try:
            with pytest.raises(ValueError):
                _ = parse_condensed_sdrf(temp_path)

            # # Should return an empty DataFrame with the correct columns
            # self.assertEqual(result.shape, (0, 6))
            # self.assertEqual(list(result.columns), ["exp_name", "assay_name", "feat_class", "feat_type", "feat_value", "ontology"])
        finally:
            os.unlink(temp_path)

    def test_file_with_headers(self):
        # Test with a file that includes headers
        content = (
            "Experiment\tEmpty\tAssay\tClass\tType\tValue\tOntology\n"
            "E-MTAB-123\t\tAssay1\tClass1\torganism\tHomo sapiens\tNCBITaxon_9606\n"
            "E-MTAB-123\t\tAssay2\tClass2\torganism\tMus musculus\tNCBITaxon_10090"
        )

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_file.write(content)
            temp_path = temp_file.name

        try:
            result = parse_condensed_sdrf(temp_path)

            # Check that headers are treated as data
            self.assertEqual(result.shape, (3, 6))
            self.assertEqual(result["exp_name"].to_list()[0], "Experiment")
            self.assertEqual(result["assay_name"].to_list()[0], "Assay")
        finally:
            os.unlink(temp_path)

    ## I don't think this can happen, but I may need to test against it if we get weird stuff going on
    # def test_file_with_tabs_in_values(self):
    #     # Test handling of extra tabs within field values (simulating malformed data)
    #     content = (
    #         "E-MTAB-123\t\tAssay1\tClass1\torganism\tHomo\tsapiens\tNCBITaxon_9606\n"  # Extra tab in value
    #         "E-MTAB-123\t\tAssay2\tClass2\torganism\tMus musculus\tNCBITaxon_10090"
    #     )

    #     with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_file:
    #         temp_file.write(content)
    #         temp_path = temp_file.name

    #     try:
    #         result = parse_condensed_sdrf(temp_path)
    #         print(result)
    #         # The first row will be parsed differently due to the extra tab
    #         self.assertEqual(result["feat_value"].to_list()[0], "Homo")
    #         self.assertEqual(result["ontology"].to_list()[0], "sapiens")
    #     finally:
    #         os.unlink(temp_path)

    def test_file_not_found(self):
        # Test with a non-existent file
        non_existent_path = "/path/to/nonexistent/file.tsv"

        with self.assertRaises(FileNotFoundError):
            parse_condensed_sdrf(non_existent_path)

    def test_path_object_vs_string(self):
        # Test that the function works with both string and Path objects
        content = "E-MTAB-123\t\tAssay1\tClass1\torganism\tHomo sapiens\tNCBITaxon_9606"

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_file.write(content)
            temp_path = temp_file.name

        try:
            # Test with string path
            result_str = parse_condensed_sdrf(temp_path)

            # Test with Path object
            result_path = parse_condensed_sdrf(Path(temp_path))

            # Both should return equivalent DataFrames
            self.assertTrue(result_str.equals(result_path))
        finally:
            os.unlink(temp_path)


if __name__ == "__main__":
    unittest.main()
