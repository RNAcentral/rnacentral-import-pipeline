# CIRCpedia V3 Import Module

This module implements the import pipeline for CIRCpedia V3 circular RNA data into RNAcentral.

## Overview

CIRCpedia V3 is a comprehensive circular RNA database containing:
- 2.6 million circular RNAs
- Data from 20 species
- Expression profiles from 2350 NGS and 63 TGS datasets
- Community-recommended nomenclature
- Genomic coordinates and annotations

**Reference**: Zhai SN, Zhang YY, Chen MH, Fu ZC, Chen LL, Ma XK, Yang L. "CIRCpedia v3: an interactive database for circular RNA characterization and functional exploration." Nucleic Acids Research, 2025. DOI: 10.1093/nar/gkaf1039

## Data Format

The module expects two files:
1. **Annotation file**: Tab-delimited text (.txt) with circRNA metadata
2. **Sequence file**: FASTA format (.fa) with circRNA sequences

### Annotation File (TSV)

#### Required Columns

- **circID**: CIRCpedia circular RNA identifier (e.g., "hsa_circ_0001")
- **species**: Species name in binomial nomenclature (e.g., "Homo sapiens")
- **Location**: Genomic location with strand in format "chr:start-end(strand)" (e.g., "chr1:12345-67890(+)")

#### Optional Columns

- **gene_Refseq**: RefSeq gene identifier
- **gene_Ensembl**: Ensembl gene identifier
- **circname**: Circular RNA name (e.g., "circ-TP53(1-5)")
- **length**: Length of circular RNA
- **subcell_location**: Subcellular localization (e.g., "cytoplasm", "nucleus")
- **editing_site**: RNA editing sites
- **DIS3_signal**: DIS3 degradation signals
- **Orthology**: Orthology information across species
- **TGS**: Third-generation sequencing support
- **transcript_Ensembl**: Ensembl transcript ID
- **transcript_Refseq**: RefSeq transcript ID

### Sequence File (FASTA)

FASTA file with sequences keyed by circID:

```fasta
>hsa_circ_0001
ATCGATCGATCGATCGATCGATCG...
>hsa_circ_0002
GCTAGCTAGCTAGCTAGCTAGCTA...
```

### Example Data

See `example_annotation.txt` and `example_sequences.fa` for sample data files.

## Architecture

The implementation follows the RNAcentral import pipeline architecture:

### Components

1. **parser.py**: Main parsing logic
   - Loads sequences from FASTA file
   - Reads TSV annotation files using Polars for efficient processing
   - Converts circRNA data to RNAcentral Entry format
   - Handles taxonomy lookups
   - Processes genomic coordinates with combined location field

2. **helpers.py**: Helper functions
   - Combined location field parsing (chr:start-end(strand))
   - Exon position parsing
   - Taxonomy ID lookup with fallback mapping
   - URL and accession generation with CIRCPEDIA: prefix
   - Note data construction from TSV fields

3. **CLI (cli/circpedia.py)**: Command-line interface
   - Integrates with the `rnac` command system
   - Command: `rnac circpedia parse <taxonomy> <annotation_file> <fasta_file> <output>`

4. **Nextflow (workflows/databases/circpedia.nf)**: Pipeline integration
   - Fetches both annotation and FASTA files from configured sources
   - Runs parser with taxonomy context
   - Outputs standardized CSV files for database loading

## Usage

### Command Line

```bash
# Parse CIRCpedia TSV annotation and FASTA sequence files
rnac circpedia parse context.db circpedia_annotation.txt circpedia_sequences.fa output_dir/ --assembly GRCh38
```

### Nextflow

```groovy
// In parse-databases.nf
include { circpedia } from './databases/circpedia'

workflow {
  build_context | set { context }
  circpedia(context) | set { data }
}
```

### Configuration

In `config/databases.config`:

```groovy
circpedia {
  run = false                    // Set to true to enable
  needs_taxonomy = true          // Requires taxonomy database
  process.directives.memory = 8.GB
  remote {
    annotation = '/path/to/circpedia_v3_annotation.txt'
    fasta = '/path/to/circpedia_v3_sequences.fa'
  }
  assembly = 'GRCh38'           // Genome assembly ID
}
```

## Output

The parser generates standard RNAcentral CSV files:

- **accessions.csv**: Sequence accessions and metadata
- **seq_short.csv** / **seq_long.csv**: RNA sequences
- **references.csv**: Publication references
- **regions.csv**: Genomic regions and coordinates
- Additional files as defined by EntryWriter

## Testing

Unit tests are located in `tests/databases/circpedia/`:

- **helpers_test.py**: Tests for helper functions
- **parser_test.py**: Tests for parser logic

Tests use minimal mock data and avoid dependencies on:
- Network resources
- External databases
- Large test data files

Run tests with:
```bash
pytest tests/databases/circpedia/
```

## Implementation Notes

### Technology Choices

- **Polars**: Used for efficient TSV processing (preferred over pandas per requirements)
- **psycopg**: Would be used for database queries if needed
- **SqliteDict**: Used for taxonomy lookups

### Circular RNA Representation

- RNA type: SO:0000593 (circular RNA)
- Database: CIRCPEDIA
- Accessions: Use CIRCpedia's own IDs with CIRCPEDIA: prefix (e.g., "CIRCPEDIA:hsa_circ_0001_1:100-200")
- URLs: Direct links to circRNA detail pages (e.g., https://bits.fudan.edu.cn/circpediav3/circrna/hsa_circ_0001)
- Coordinate system: 1-based, fully-closed (same as GFF/GTF)

### Sequence Handling

- Sequences are loaded from separate FASTA file
- Sequences must be DNA (T not U) - automatically converted if needed
- If a circID is missing from FASTA, the entry is skipped with a warning
- Back-splicing junction sequences may need special handling in production

### Species Support

The module includes fallback taxonomy mapping for common species:
- Human (9606)
- Mouse (10090)
- Rat (10116)
- Zebrafish (7955)
- Fruit fly (7227)
- C. elegans (6239)
- And 14 more common model organisms

## Future Enhancements

Potential improvements for production use:

1. **Sequence Extraction**: Integrate genome sequence extraction for circRNAs
2. **Back-splice Junction**: Handle circular junction sequences
3. **Expression Data**: Expand expression profile support
4. **Alternative Splicing**: Track alternative back-splicing events
5. **Conservation**: Import conservation scores
6. **Validation**: Add data quality checks
7. **Batch Processing**: Handle large datasets in chunks

## Data Source

CIRCpedia V3 is available at: https://bits.fudan.edu.cn/circpediav3

For data download options and format details, consult the CIRCpedia V3 documentation.

## License

Copyright [2009-2025] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0
