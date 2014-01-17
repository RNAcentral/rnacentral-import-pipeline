# Copyright [2009-2014] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
	Create sql scripts for updating the rnc_accessions table
	for immediate data visualisation on the site.
	Makes it easy to import tmRNA Website accessions at any time,
	even if the data gets overwritten.

    external id is the tmRNA external id
    optional id is the 2-piece tmRNA mate id
"""


def analyze_file(filename, filehandle):
	f = open(filename)
	for line in f:
		data = line.split('\t')
		try:
			query = "UPDATE rnc_accessions SET external_id = '%s' WHERE parent_ac='%s' and database = 'TMRNA_WEB' and is_composite = 'Y';\n" % (data[1].strip(), data[0])
			filehandle.write(query)
		except Exception, e:
			print 'Problem with line %s' % line
	f.close()

def analyze_two_piece_tmrna_file(filename, filehandle):
	f = open(filename)

	for line in f:
		data = line.split('\t')
		try:
			query = "UPDATE rnc_accessions SET external_id = '%s', optional_id = '%s' WHERE parent_ac='%s' and database = 'TMRNA_WEB' and is_composite = 'Y';\n" % (data[1], data[2].strip(), data[0])
			filehandle.write(query)
		except Exception, e:
			print 'Problem with line %s' % line
	f.close()


outputfile = open('update_tmrna_accessions.sql', 'w')

# import tmRNA_Web gene file
analyze_file('tmRNAdb_20131122_DNA.tsv', outputfile)

# import tmRNA_Web 1-piece tmRNA file
analyze_file('tmRNAdb_20131122_RNA1P.tsv', outputfile)

# import tmRNA_Web 2-piece file
analyze_two_piece_tmrna_file('tmRNAdb_20131122_RNA2P.tsv', outputfile)

outputfile.write('COMMIT;')

outputfile.close()
