"""
Copyright [2009-2014] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Usage:

python path/to/this/file.py <Greengenes/Noncode>
    --local-scheduler
    --destination /path/to/output/files
    --input-folder /path/to/input/folder/
"""

import luigi

from glob import glob
from json_parser_greengenes import JsonParserGreengenes
from json_parser_noncode import JsonParserNoncode


class Noncode(luigi.Task): # pylint: disable=W0232
    """
    Luigi task for processing all Json files found in an input folder.
    """
    input_folder = luigi.Parameter()
    destination = luigi.Parameter(default='/tmp')
    test = luigi.BoolParameter(default=False, significant=False)

    def requires(self):
        return [JsonParserNoncode(input_file=f, destination=self.destination, test=self.test) for f in glob(self.input_folder + '*.json')]


class Greengenes(luigi.Task): # pylint: disable=W0232
    """
    Luigi task for processing all Json files found in an input folder.
    """
    input_folder = luigi.Parameter()
    destination = luigi.Parameter(default='/tmp')
    test = luigi.BoolParameter(default=False, significant=False)

    def requires(self):
        return [JsonParserGreengenes(input_file=f, destination=self.destination, test=self.test) for f in glob(self.input_folder + '*.json')]


# main entry point
if __name__ == '__main__':
    luigi.run()
