#!/usr/bin/env python

"""
This is the main entry point to various parts of the RNAcentral pipeline. This
is mean to provide a single point for the processing of data from external
databases. In order for loading to be efficient is is handled by pgloader
externally.
"""

import logging

from rnacentral_pipeline.cli import cli

logging.basicConfig()
cli()
