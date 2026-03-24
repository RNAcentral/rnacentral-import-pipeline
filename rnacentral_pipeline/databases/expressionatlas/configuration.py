# This module handles the parsing of the configuration file
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from lxml import etree


@dataclass
class Contrast:
    id: str
    name: str
    ref_group: str  # Mapped from reference_assay_group
    test_group: str  # Mapped from test_assay_group


@dataclass
class Contrasts:
    contrast: List[Contrast] = field(default_factory=list)


@dataclass
class Assay:
    id: str
    technical_replicate_id: Optional[str] = None


@dataclass
class AssayGroup:
    id: str
    label: Optional[str] = None
    assays: List[str] = field(default_factory=list)


@dataclass
class AssayGroups:
    assay_group: List[AssayGroup] = field(default_factory=list)


@dataclass
class Analytics:
    assay_groups: AssayGroups
    array_design: Optional[str] = None
    contrasts: Optional[Contrasts] = None


@dataclass
class Config:
    exp_type: str  # Mapped from experimentType
    r_data: Optional[str] = None  # Added to match the example XML
    analytics: List[Analytics] = field(default_factory=list)


def parse_config(file_path):
    """Parse the XML configuration file into a Config object."""
    try:
        with open(file_path, "r") as file:
            xml_content = file.read()

        root = etree.fromstring(xml_content.encode("utf-8"))

        # Create Config object - the root is 'configuration' with experimentType attribute
        config = Config(
            exp_type=root.get("experimentType"), r_data=root.get("r_data"), analytics=[]
        )

        # Parse analytics
        for analytics_elem in root.findall("analytics"):
            assay_groups = AssayGroups(assay_group=[])

            # Parse assay groups - in the example, they're directly under 'assay_groups'
            for assay_group_elem in analytics_elem.findall(
                "./assay_groups/assay_group"
            ):
                assays = []

                # Process each assay element and extract technical_replicate_id if present
                for assay_elem in assay_group_elem.findall("assay"):
                    assays.append(assay_elem.text)

                assay_group = AssayGroup(
                    id=assay_group_elem.get("id"),
                    label=assay_group_elem.get("label"),
                    assays=assays,
                )
                assay_groups.assay_group.append(assay_group)

            # The example doesn't have contrasts, but keeping the code for completeness
            contrasts_elem = analytics_elem.find("contrasts")
            contrasts = None
            if contrasts_elem is not None:
                contrasts = Contrasts(contrast=[])
                for contrast_elem in contrasts_elem.findall("contrast"):
                    contrast = Contrast(
                        id=contrast_elem.get("id"),
                        name=contrast_elem.find("name").text,
                        ref_group=contrast_elem.find("reference_assay_group").text
                        if contrast_elem.find("reference_assay_group") is not None
                        else contrast_elem.find("ref_group").text,
                        test_group=contrast_elem.find("test_assay_group").text
                        if contrast_elem.find("test_assay_group") is not None
                        else contrast_elem.find("test_group").text,
                    )
                    contrasts.contrast.append(contrast)

            # The example doesn't have array_design, but keeping the code for completeness
            array_design_elem = analytics_elem.find("arrayDesign")
            array_design = (
                array_design_elem.text if array_design_elem is not None else None
            )

            analytics = Analytics(
                assay_groups=assay_groups,
                array_design=array_design,
                contrasts=contrasts,
            )
            config.analytics.append(analytics)

        return config

    except Exception as e:
        raise Exception(f"Error parsing configuration file: {e}")
