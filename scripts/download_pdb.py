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

Download RNA sequences from PDB
and prepare the data for importing into RNAcentral.

Algorithm outline:

* Get RNA chains from all RNA-containing PDB files
* Filter out mRNAs
* Format the data as in the rest of the data import pipeline.
"""

import csv
import hashlib
import json
import re
import requests # pip install requests
import xml.etree.ElementTree as ET


def get_rna_containing_pdb_ids():
    """
    Get PDB ids of all RNA-containing 3D structures
    using the RCSB PDB REST API.
    """
    query = """
    <orgPdbQuery>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <containsProtein>I</containsProtein>
    <containsDna>I</containsDna>
    <containsRna>Y</containsRna>
    <containsHybrid>I</containsHybrid>
    </orgPdbQuery>
    """
    url = 'http://www.rcsb.org/pdb/rest/search'

    request = requests.post(url, data=query,
        # without this header the request is redirected incorrectly
        headers={'content-type': 'application/x-www-form-urlencoded'})

    if request.status_code == 200:
        pdb_ids = request.text.rstrip().split('\n')
    else:
        pdb_ids = None

    return pdb_ids


def get_custom_report(pdb_ids, fields):
    """
    Get custom report about PDB files in a tabular format.
    """
    url = 'http://www.rcsb.org/pdb/rest/customReport.csv'
    data = {
        'pdbids': ','.join(pdb_ids),
        'customReportColumns': ','.join(fields),
        'format': 'csv',
    }
    request = requests.post(url, data=data)

    if request.status_code == 200:
        return request.text.split('<br />')
    else:
        return None


def get_chain_descriptions(pdb_ids):
    """
    Get per-chain information about each RNA sequence.
    Return a dictionary that looks like this:
    {
        '1S72_A': {'structureId': '1S72', 'chainId': '0' etc}
    }
    """
    fields = [
        'structureId',
        'chainId',
        'structureTitle',
        'experimentalTechnique',
        'releaseDate',
        'ndbId',
        'emdbId',
        'classification',
        'entityId',
        'sequence',
        'chainLength',
        'db_id',
        'db_name',
        'entityMacromoleculeType',
        'source',
        'taxonomyId',
        'compound',
        'resolution',
    ]

    def filter_out_chains(report):
        """
        * read custom report
        * skip non-rna chains
        * skip chains shorter than 10 nucleotides
        * skip mRNAs
        * return data for RNA-containing chains
        """
        reader = csv.DictReader(report, delimiter=',', quotechar='"')
        data = dict()
        min_seq_length = 10
        for row in reader:
            # skip proteins
            if 'RNA' not in row['entityMacromoleculeType']:
                continue
            # skip short sequences
            if int(row['chainLength']) < min_seq_length:
                continue
            # skip mRNAs
            if re.search('mRNA', row['compound'], re.IGNORECASE) and \
               not re.search('tmRNA', row['compound'], re.IGNORECASE):
                continue

            # use entityId to ensure that the id is unique when chainIds
            # are only different in case ('A' and 'a')
            accession = '{structureId}_{chainId}_{entityId}'.format(**row)
            # if no taxonomy id, use that of the synthetic construct
            if row['taxonomyId'] == '':
                row['taxonomyId'] = '32630' # synthetic construct
            data[accession] = row
        return data

    report = get_custom_report(pdb_ids, fields)
    data = filter_out_chains(report)
    return data


def get_literature_references(pdb_ids):
    """
    Get literature citations for each PDB file.
    Return a dictionary that looks like this:
    {
        '1S72': {'structureId': '1S72', etc}
    }
    """
    fields = [
        'structureId',
        'citationAuthor',
        'firstPage',
        'lastPage',
        'journalName',
        'title',
        'volumeId',
        'publicationYear',
        'pubmedId',
        'pmc',
        'doi',
    ]
    report = get_custom_report(pdb_ids, fields)

    reader = csv.DictReader(report, delimiter=',', quotechar='"')
    data = dict()
    for row in reader:
        data[row['structureId']] = row
    return data


def get_full_lineages(chains):
    """
    Based on taxid retrieve full taxonomic lineage
    and common name using the ENA API.
    """
    taxids = dict()
    url = 'http://www.ebi.ac.uk/ena/data/view/Taxon:{taxid}&display=xml'

    for accession, data in chains.iteritems():
        taxid = data['taxonomyId']
        # skip if this taxid was processed
        if taxid in data:
            continue
        taxids[taxid] = dict()

        request = requests.get(url.format(taxid=taxid))

        root = ET.fromstring(request.text)

        organism = root.findall('taxon')[0]
        taxids[taxid]['species'] = organism.get('scientificName')
        taxids[taxid]['common_name'] = organism.get('commonName') or ''

        lineage = []
        for taxon in root.findall('./taxon/lineage/taxon'):
            if 'hidden' not in taxon.attrib:
                lineage.append(taxon.get('scientificName'))
        lineage.reverse()
        lineage.append(taxids[taxid]['species'])
        taxids[taxid]['lineage'] = '; '.join(lineage)

    return taxids


def get_md5(input_string):
    """
    Get MD5 hash of a string.
    """
    md5 = hashlib.md5()
    md5.update(input_string)
    return md5.hexdigest()


def write_sequence_files(chains):
    """
    Write sequence files to be imported into the RNAcentral database.
    """
    filename_short = 'pdb_short.csv'
    filename_long = 'pdb_long.csv'
    f_short = open(filename_short, 'w')
    f_long = open(filename_long, 'w')
    max_seq_length = 4000

    for accession, data in chains.iteritems():
        sequence = data['sequence'].replace('U', 'T').upper()
        fields = {
            'crc64': CRC64digest(sequence),
            'length': data['chainLength'],
            'seq': sequence,
            'dbid': 'PDB',
            'accession': accession,
            'optional_id': '',
            'version': 1,
            'taxid': data['taxonomyId'],
            'md5': get_md5(sequence),
        }
        line = '{crc64},{length},{seq},{dbid},{accession},' + \
               '{optional_id},{version},{taxid},{md5}'
        line = line.format(**fields)
        if int(data['chainLength']) < max_seq_length:
            f_short.write(line + "\n")
        else:
            f_long.write(line + "\n")

    f_short.close()
    f_long.close()


def sanitize(input_string):
    """
    Remove leading and trailing quotes
    and scape the remaining double quotes.
    """
    input_string = re.sub('^"+', '', input_string) # leading quotes
    input_string = re.sub('"+$', '', input_string) # trailing quotes
    input_string = re.sub('"', '""', input_string) # escape double quotes
    return input_string


def write_literature_references(chains, refs):
    """
    For each RNA chain write out literature citations
    to be imported into the RNAcentral database.
    """
    def get_journal_name():
        """
        Format the journal name as a single string.
        Example:
        Genome Res. 22(9):1775-1789(2012).
        """
        return '{journalName} {volumeId}:{firstPage}-{lastPage}({publicationYear}).'.format(**ref)

    def get_unique_id():
        """
        authors + location + title
        This unique string is used during the database update
        to identify distinct references.
        """
        return '{citationAuthor}{citationAuthor}{title}'.format(**ref)

    filename = 'pdb_refs.csv'
    f = open(filename, 'w')

    for accession, data in chains.iteritems():
        ref = refs[data['structureId']] # primary reference of this PDB id
        ref['journal'] = get_journal_name()
        fields = {
            'md5': get_md5(get_unique_id()),
            'accession': accession,
            'authors': sanitize(ref['citationAuthor']),
            'location': sanitize(ref['journal']),
            'title': sanitize(ref['title']),
            'pmid': sanitize(ref['pubmedId']),
            'doi': sanitize(ref['doi']),
        }
        line = '"{md5}","{accession}","{authors}","{location}",' + \
               '"{title}","{pmid}","{doi}"'
        line = line.format(**fields)
        f.write(line + "\n")

    f.close()


def write_accession_info(chains, species):
    """
    For each RNA chain write out the annotation information.
    """
    def get_description():
        """
        Get description of the RNA chain.
        Truncate compound if it is too long.
        """
        max_length = 80
        compound =  data['compound'][:max_length] + \
                    (data['compound'][max_length:] and '...')
        line = compound + ' from {source} (PDB {structureId}, chain {chainId})'.format(**data)
        return line

    def get_ncrna_type():
        """
        INSDC RNA types include:
        * tRNA
        * tmRNA
        * rRNA
        * precursor_RNA
        * misc_RNA
        * ncRNA (subtypes http://www.insdc.org/rna_vocab.html)
        """
        data['feature_name'] = 'misc_RNA'
        data['ncrna_class'] = ''
        compound = data['compound'].upper()
        ribosomes = ['5S', '5.8S', '16S', '18S', '23S', '28S', '30S', '40S', '60S', '80S']

        # tRNA
        if 'TRNA' in compound:
            data['feature_name'] = 'tRNA'
        # tmRNA
        if 'TMRNA' in compound:
            data['feature_name'] = 'tmRNA'
        # rRNA
        elif [x for x in ribosomes if x in compound]:
            data['feature_name'] = 'rRNA'
        # SRP
        elif 'SRP' in compound:
            data['feature_name'] = 'ncRNA'
            data['ncrna_class'] = 'SRP_RNA'
        # Ribozyme
        elif 'RIBOZYME' in compound and 'HAMMERHEAD' not in compound:
            data['feature_name'] = 'ncRNA'
            data['ncrna_class'] = 'ribozyme'
        # Hammerhead ribozyme
        elif 'RIBOZYME' in compound and 'HAMMERHEAD' in compound:
            data['feature_name'] = 'ncRNA'
            data['ncrna_class'] = 'hammerhead_ribozyme'
        # snRNA
        elif 'SNRNA' in compound:
            data['feature_name'] = 'ncRNA'
            data['ncrna_class'] = 'snRNA'
        # snoRNA
        elif 'SNORNA' in compound:
            data['ncrna_class'] = 'snoRNA'
            data['feature_name'] = 'ncRNA'

    def get_db_xref():
        """
        Put NDB and EMDB xrefs in the db_xref field.
        """
        line = ''
        if data['ndbId']:
            line += 'NDB:%s' % data['ndbId']
        if data['db_name']:
            line += ' %s:%s' % (data['db_name'], data['db_id'])
        return sanitize(line.strip())

    def get_structured_note():
        """
        Store additional information in the structured note.
        """
        notes = {}
        fields = [
            'structureTitle',
            'experimentalTechnique',
            'resolution',
            'releaseDate',
        ]
        for field in fields:
            if field in data and data[field]:
                notes[field] = data[field]
        line = json.dumps(notes)
        return sanitize(line)

    filename = 'pdb_ac_info.csv'
    f = open(filename, 'w')

    for accession, data in chains.iteritems():
        get_ncrna_type()
        taxon = species[data['taxonomyId']]
        if data['source'] == '':
            data['source'] = taxon['species']
        fields = {
            'accession': accession,
            'parent_ac': data['structureId'],
            'seq_version': data['entityId'],
            'feature_start': 1,
            'feature_end': data['chainLength'],
            'feature_name': data['feature_name'],
            'ordinal': 1,
            'is_composite': 'N',
            'non_coding_id': '',
            'database': 'PDBE',
            'external_id': data['structureId'],
            'optional_id': data['chainId'],
            'project': '',
            'division': '',
            'keywords': '',
            'description': get_description(),
            'species': taxon['species'],
            'common_name': taxon['common_name'],
            'organelle': '',
            'classification': taxon['lineage'],
            'allele': '',
            'anticodon': '',
            'chromosome': '',
            'experiment': '',
            'function': '',
            'gene': '',
            'gene_synonym': '',
            'inference': '',
            'locus_tag': '',
            'map': '',
            'mol_type': '',
            'ncrna_class': data['ncrna_class'],
            'note': get_structured_note(),
            'old_locus_tag': '',
            'operon': '',
            'product': data['compound'],
            'pseudogene': '',
            'standard_name': '',
            'db_xref': get_db_xref(),
        }
        line = '"{accession}","{parent_ac}","{seq_version}","{feature_start}",' + \
               '"{feature_end}","{feature_name}","{ordinal}","{is_composite}",' + \
               '"{non_coding_id}","{database}","{external_id}","{optional_id}",' + \
               '"{project}","{division}","{keywords}","{description}","{species}",' + \
               '"{common_name}","{organelle}","{classification}","{allele}",' + \
               '"{anticodon}","{chromosome}","{experiment}","{function}","{gene}",' + \
               '"{gene_synonym}","{inference}","{locus_tag}","{map}","{mol_type}",' + \
               '"{ncrna_class}","{note}","{old_locus_tag}","{operon}","{product}",' + \
               '"{pseudogene}","{standard_name}","{db_xref}"'
        line = line.format(**fields)
        f.write(line + "\n")

    f.close()


def CRC64(aString):
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
            for j in xrange(8):
                rflag = partl & 1L
                partl >>= 1L
                if parth & 1:
                    partl |= (1L << 31L)
                parth >>= 1L
                if rflag:
                    parth ^= POLY64REVh
            CRCTableh[i] = parth
            CRCTablel[i] = partl

    for item in aString:
        shr = 0L
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8L
        temp1l = (crcl >> 8L) | shr
        tableindex = (crcl ^ ord(item)) & 0xFF

        crch = temp1h ^ CRCTableh[tableindex]
        crcl = temp1l ^ CRCTablel[tableindex]
    return (crch, crcl)


def CRC64digest(aString):
    """
    Get string-formatted CRC64.
    """
    return "%08X%08X" % (CRC64(aString))


def main():
    """
    Main program.
    """

    pdb_ids = get_rna_containing_pdb_ids()
    print len(pdb_ids)

    # pdb_ids = ['1S72', '1J5E']
    # pdb_ids = pdb_ids[:10]

    chains = get_chain_descriptions(pdb_ids)
    refs = get_literature_references(pdb_ids)
    species = get_full_lineages(chains)

    # import sys
    # sys.exit()

    write_sequence_files(chains)
    write_literature_references(chains, refs)
    write_accession_info(chains, species)


# main entry point
if __name__ == "__main__":
    main()
