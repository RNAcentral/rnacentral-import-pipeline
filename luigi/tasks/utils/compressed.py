# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

import os
import gzip
import tarfile

from luigi.local_target import atomic_file


def expand_tar(filename, directory):
    """
    This will expand a tarball (compressed or otherwise) in to the given
    directory. It will provide a list of the paths that are were produced.
    """

    with tarfile.open(filename, 'r:*') as tar:
        tar.extractall(directory)
        return [os.path.join(directory, f) for f in tar.getnames()]


def expand_gzip(filename, directory):
    """
    This will expand a gzipped file into the given directory. It will return a
    list with one element, the path to the file that was expanded.
    """

    basename = os.path.basename(filename)
    output_filename = os.path.join(directory, basename.replace('.gz', ''))
    with atomic_file(output_filename) as out:
        with gzip.open(filename, 'r') as raw:
            out.write(raw.read())
        return [output_filename]


def expand(filename, directory):
    """
    Expand the given file into the given directory. This can expand tarball's
    (compressed or not) and gzip files. Gzip files must be named ending in
    '.gz'. If the given file is neither of those then a ValueError will be
    thrown.
    """

    if tarfile.is_tarfile(filename):
        return expand_tar(filename, directory)

    if filename.endswith('.gz'):
        return expand_gzip(filename, directory)

    raise ValueError("Unknown compressed file type: " + filename)
