#!/usr/bin/env python3

"""
This is the main entry point to various parts of the RNAcentral pipeline. This
is mean to provide a single point for the processing of data from external
databases. In order for loading to be efficient is is handled by pgloader
externally.
"""

import logging

from rnacentral_pipeline.cli import cli

# Stays on stderr: several commands write their data to stdout (`output`
# defaults to `-`), so log lines there would corrupt the payload.
logging.basicConfig(
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
cli()
