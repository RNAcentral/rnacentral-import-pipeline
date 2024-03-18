# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import typing as ty
import xml.etree.ElementTree as ET
from pathlib import Path


def parse(file_path):
    """
    Parse a config file into a python dictionary
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    config_dict = {
        "experimentType": root.attrib["experimentType"].strip(),
        "r_data": root.attrib["r_data"],
    }

    analytics_dict = {}
    analytics = root.find("analytics")
    if analytics.find("array_design") is not None:
        analytics_dict["array_design"] = analytics.find("array_design").text.strip()
    else:
        analytics_dict["array_design"] = None

    assay_groups_dict = {}
    assay_groups = analytics.find("assay_groups")
    for assay_group in assay_groups.findall("assay_group"):
        group_id = assay_group.attrib["id"]
        assay_list = [assay.text for assay in assay_group.findall("assay")]
        assay_groups_dict[group_id] = assay_list

    analytics_dict["assay_groups"] = assay_groups_dict

    contrasts_dict = {}
    contrasts = analytics.find("contrasts")
    if contrasts is not None:
        for contrast in contrasts.findall("contrast"):
            contrast_id = contrast.attrib["id"]
            contrast_dict = {
                "name": contrast.find("name").text,
                "reference_assay_group": contrast.find("reference_assay_group").text,
                "test_assay_group": contrast.find("test_assay_group").text,
            }
            contrasts_dict[contrast_id] = contrast_dict

        analytics_dict["contrasts"] = contrasts_dict

    config_dict["analytics"] = analytics_dict

    return config_dict


def extract_experiment_name(file_path):
    """
    Experiment config filenames look like this: E-ATMX-6-configuration.xml
    This should give the experiment name E-ATMX-6
    """
    basename = Path(file_path).name
    parts = basename.split("-")
    experiment_name = "-".join(parts[:3])
    return experiment_name


def build_exp_config_lookup(file_list: ty.List[Path]):
    """
    file_list should be a list of Path objects
    """
    lookup_dict = {}
    for file in file_list:
        exp_name = extract_experiment_name(file)
        config_dict = parse(file)
        lookup_dict[exp_name] = config_dict
    return lookup_dict
