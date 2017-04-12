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

from urlparse import urlparse

from luigi.parameter import Parameter
from luigi.format import Gzip as GzipFormat
from luigi.contrib.ftp import RemoteTarget as FtpTarget
from luigi.local_target import LocalTarget


class GenericFileParameter(Parameter):
    def serialize(self, value):
        return value

    def parse(self, value):
        try:
            self.as_target(value)
            return value
        except:
            raise ValueError("Unknown type of file %s" % value)

    def as_target(self, path, **kwargs):
        parsed = urlparse(path)
        fmt = None
        if path.endswith('.gz'):
            fmt = GzipFormat
        kwargs['format'] = fmt

        if not parsed.scheme or parsed.scheme == 'file':
            return LocalTarget(path=os.path.normalize(path), **kwargs)
        if parsed.scheme == 'ftp':
            return FtpTarget(parsed.path, parsed.netloc, **kwargs)

        raise ValueError("Cannot build target for %s" % path)


class CommaGenericFileParameter(Parameter):
    def seralize(self, value):
        return ','.join(sorted(value))

    def parse(self, value):
        if not value:
            raise ValueError("Must provide files to deduplicate")

        try:
            targets = self.as_targets(value)
            return [path.strip() for path in value.split(',')]
        except:
            raise ValueError("Unhandable type of files %s" % value)

    def as_targets(self, value, **kwargs):
        parameter = GenericFileParameter()
        paths = [p.strip() for p in value.split(',')]
        return [parameter.as_target(path) for path in paths]


class PathParameter(Parameter):
    def seralize(self, value):
        return value

    def parse(self, value):
        return os.path.normalie(value)
