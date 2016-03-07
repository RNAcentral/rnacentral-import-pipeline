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
"""

import hashlib


class RNAcentralEntry(object):
    """
    A class that represents an RNAcentral entry and can prepare the data for
    writing out to csv files.
    """
    def __init__(self, **kwargs):
        self.accession = ''
        self.allele = ''
        self.anticodon = ''
        self.assembly_info = []
        self.chromosome = ''
        self.common_name = ''
        self.database = ''
        self.db_xrefs = ''
        self.description = ''
        self.division = ''
        self.experiment = ''
        self.feature_location_end = ''
        self.feature_location_start = ''
        self.feature_type = ''
        self.function = ''
        self.gene = ''
        self.gene_synonym = ''
        self.inference = ''
        self.is_composite = ''
        self.keywords = ''
        self.lineage = ''
        self.locus_tag = ''
        self.map = ''
        self.mol_type = ''
        self.ncbi_tax_id = ''
        self.ncrna_class = ''
        self.non_coding_id = ''
        self.note = ''
        self.old_locus_tag = ''
        self.operon = ''
        self.optional_id = ''
        self.ordinal = ''
        self.organelle = ''
        self.parent_accession = ''
        self.primary_id = ''
        self.product = ''
        self.project = ''
        self.pseudogene = ''
        self.references = []
        self.sequence = '',
        self.seq_version = ''
        self.species = ''
        self.standard_name = ''
        self.__dict__.update(kwargs)

    def CRC64digest(self, input_string):
        """
        Python re-implementation of SWISS::CRC64
        Adapted from:
        http://code.activestate.com/recipes/259177-crc64-calculate-the-cyclic-redundancy-check/
        """
        POLY64REVh = 0xd8000000L
        CRCTableh = [0] * 256
        CRCTablel = [0] * 256
        isInitialized = False
        crcl = 0
        crch = 0
        if isInitialized is not True:
            isInitialized = True
            for i in xrange(256):
                partl = i
                parth = 0L
                for _ in xrange(8):
                    rflag = partl & 1L
                    partl >>= 1L
                    if parth & 1:
                        partl |= (1L << 31L)
                    parth >>= 1L
                    if rflag:
                        parth ^= POLY64REVh
                CRCTableh[i] = parth
                CRCTablel[i] = partl

        for item in input_string:
            shr = 0L
            shr = (crch & 0xFF) << 24
            temp1h = crch >> 8L
            temp1l = (crcl >> 8L) | shr
            tableindex = (crcl ^ ord(item)) & 0xFF

            crch = temp1h ^ CRCTableh[tableindex]
            crcl = temp1l ^ CRCTablel[tableindex]
        return "%08X%08X" % (crch, crcl)

    def get_md5(self, input_string):
        """
        Get md5 hash based on a string.
        """
        md5 = hashlib.md5()
        md5.update(input_string)
        result = md5.hexdigest()
        if len(result) != 32:
            class MD5Exception(Exception):
                pass
            raise MD5Exception('MD5 value is too short! %s' % result)
        return result

    def md5(self):
        """
        Convenience method for getting md5 of the RNA sequence.
        """
        return self.get_md5(self.sequence)

    def crc64(self):
        """
        Convenience method for getting crc64 of the RNA sequence.
        """
        return self.CRC64digest(self.sequence)

    def is_valid(self):
        """
        Check that the entry meets minimum quality requirements.
        """
        length = len(self.sequence)
        if length < 10:
            if verbose:
                print 'WARN: {accession} is too short (length = {length})'.format(accession=self.accession, length=length)
            return False

        if length > 1000000:
            if verbose:
                print 'WARN: {accession} is too long (length = {length})'.format(accession=self.accession, length=length)
            return False

        mandatory_fields = ['accession', 'sequence', 'primary_id', 'ncbi_tax_id']
        for field in mandatory_fields:
            if not getattr(self, field):
                if verbose:
                    print 'WARN: Attribute {field} is missing in {accession}'.format(field=field, accession=self.accession)
                return False
        return True

    def format_references(self):
        """
        Format the entry for writing out into csv files.
        """
        lines = []
        for reference in self.references:
            lines.append([
                self.get_md5(''.join([reference['authors'], reference['location'], reference['title']])),
                self.accession,
                reference['authors'],
                reference['location'],
                reference['title'],
                reference['pmid'],
                reference['doi'],
            ])
        return lines

    def format_sequence_line(self):
        """
        Format the entry for writing out into csv files.
        """
        return [
            self.crc64(),
            len(self.sequence),
            self.sequence,
            self.database,
            self.accession,
            self.optional_id,
            self.seq_version,
            self.ncbi_tax_id,
            self.md5(),
        ]

    def format_genomic_locations(self):
        """
        Format the entry for writing out into csv files.
        """
        lines = []
        for exon in self.assembly_info:
            if 'complement' in exon:
                complement = -1
            else:
                complement = 1
            lines.append([
                self.accession,
                '.'.join([self.parent_accession, self.seq_version]),
                exon['primary_start'],
                exon['primary_end'],
                complement,
            ])
        return lines

    def format_ac_line(self):
        """
        Format the entry for writing out into csv files.
        """
        return  [
            self.accession,
            self.parent_accession,
            self.seq_version,
            self.feature_location_start,
            self.feature_location_end,
            self.feature_type,
            self.ordinal,
            self.is_composite,
            self.non_coding_id,
            self.database,
            self.primary_id,
            self.optional_id,
            self.project,
            self.division,
            self.keywords,
            self.description,
            self.species,
            self.common_name,
            self.organelle,
            self.lineage,
            self.allele,
            self.anticodon,
            self.chromosome,
            self.experiment,
            self.function,
            self.gene,
            self.gene_synonym,
            self.inference,
            self.locus_tag,
            self.map,
            self.mol_type,
            self.ncrna_class,
            self.note,
            self.old_locus_tag,
            self.operon,
            self.product,
            self.pseudogene,
            self.standard_name,
            self.db_xrefs,
        ]
