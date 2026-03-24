import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import polars as pl
import pytest

from rnacentral_pipeline.databases.expressionatlas import sdrf

# Import the module where find_all_taxids is located
# from your_module import find_all_taxids, sdrf
from rnacentral_pipeline.databases.expressionatlas.helpers import find_all_taxids


class TestFindAllTaxids(unittest.TestCase):
    def setUp(self):
        # Create sample data that mimics the structure returned by sdrf.parse_condensed_sdrf
        self.sample_data1 = pl.DataFrame(
            {
                "feat_type": ["organism", "other", "organism"],
                "ontology": ["NCBITaxon_9606", "something_else", "NCBITaxon_10090"],
            }
        )

        self.sample_data2 = pl.DataFrame(
            {
                "feat_type": ["organism", "other", "organism"],
                "ontology": ["NCBITaxon_10090", "another_thing", "NCBITaxon_7227"],
            }
        )

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_empty_directory(self, mock_parse):
        # Test with a directory that has no matching files
        with tempfile.TemporaryDirectory() as temp_dir:
            with pytest.raises(FileNotFoundError):
                _ = find_all_taxids(temp_dir)
            mock_parse.assert_not_called()

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_single_file(self, mock_parse):
        # Test with a directory containing a single sdrf file
        mock_parse.return_value = self.sample_data1

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create a mock sdrf file
            file_path = os.path.join(temp_dir, "test-condensed-sdrf.tsv")
            with open(file_path, "w") as f:
                f.write("dummy content")

            # Patch Path.rglob to return our mock file
            with patch("pathlib.Path.rglob", return_value=[Path(file_path)]):
                result = find_all_taxids(temp_dir)
                self.assertEqual(result, [9606, 10090])
                mock_parse.assert_called_once()

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_multiple_files(self, mock_parse):
        # Test with a directory containing multiple sdrf files
        mock_parse.side_effect = [self.sample_data1, self.sample_data2]

        with tempfile.TemporaryDirectory() as temp_dir:
            # Create mock sdrf files
            file_path1 = os.path.join(temp_dir, "test1-condensed-sdrf.tsv")
            file_path2 = os.path.join(temp_dir, "test2-condensed-sdrf.tsv")

            for path in [file_path1, file_path2]:
                with open(path, "w") as f:
                    f.write("dummy content")

            # Patch Path.rglob to return our mock files
            with patch(
                "pathlib.Path.rglob", return_value=[Path(file_path1), Path(file_path2)]
            ):
                result = find_all_taxids(temp_dir)
                self.assertEqual(sorted(result), [7227, 9606, 10090])
                self.assertEqual(mock_parse.call_count, 2)

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_no_organism_entries(self, mock_parse):
        # Test when there are no organism entries in the data
        mock_parse.return_value = pl.DataFrame(
            {
                "feat_type": ["gene", "other", "protein"],
                "ontology": ["value1", "value2", "value3"],
            }
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, "test-condensed-sdrf.tsv")
            with open(file_path, "w") as f:
                f.write("dummy content")

            with patch("pathlib.Path.rglob", return_value=[Path(file_path)]):
                result = find_all_taxids(temp_dir)
                self.assertEqual(result, [])

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_invalid_taxid_format(self, mock_parse):
        # Test handling of invalid taxid formats
        mock_parse.return_value = pl.DataFrame(
            {
                "feat_type": ["organism", "organism"],
                "ontology": ["NCBITaxon_invalid", "NCBITaxon_123abc"],
            }
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, "test-condensed-sdrf.tsv")
            with open(file_path, "w") as f:
                f.write("dummy content")

            with patch("pathlib.Path.rglob", return_value=[Path(file_path)]):
                # This should raise an exception because casting to Int64 will fail
                with self.assertRaises(Exception):
                    find_all_taxids(temp_dir)

    @patch("rnacentral_pipeline.databases.expressionatlas.sdrf.parse_condensed_sdrf")
    def test_different_ontology_format(self, mock_parse):
        # Test with different ontology format
        mock_parse.return_value = pl.DataFrame(
            {
                "feat_type": ["organism", "organism"],
                "ontology": ["OtherPrefix_9606", "NotNCBI_10090"],
            }
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            file_path = os.path.join(temp_dir, "test-condensed-sdrf.tsv")
            with open(file_path, "w") as f:
                f.write("dummy content")

            with patch("pathlib.Path.rglob", return_value=[Path(file_path)]):
                # The function should throw an exception on unexpected taxonomies
                with pytest.raises(ValueError):
                    _ = find_all_taxids(temp_dir)


if __name__ == "__main__":
    unittest.main()
