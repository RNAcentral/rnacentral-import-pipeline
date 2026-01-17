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

The module expects CSV files with the following columns:

### Required Columns

- **circid**: CIRCpedia circular RNA identifier (e.g., "hsa_circ_0001")
- **species**: Species name in binomial nomenclature (e.g., "Homo sapiens")
- **location**: Genomic location in format "chr:start-end" (e.g., "chr1:12345-67890")

### Optional Columns

- **strand**: DNA strand ('+' or '-')
- **gene**: Host gene symbol
- **sequence**: RNA sequence (if not provided, placeholder will be used)
- **exon_positions** or **exonstart_exonend**: Exon coordinates in format "start1-end1,start2-end2,..."
- **fpm**: Expression level in Fragments Per Million
- **cell_line** or **tissue**: Cell line or tissue type
- **sequencing_type**: Type of sequencing used
- **conservation**: Conservation information
- **rnase_r_enrichment**: RNase R enrichment fold change

### Example CSV

```csv
circid,species,location,strand,gene,sequence,exon_positions,fpm,cell_line
hsa_circ_0001,Homo sapiens,chr1:1000-2000,+,TP53,ATCGATCG...,1000-1200;1400-2000,10.5,HeLa
hsa_circ_0002,Homo sapiens,chr2:3000-4000,-,MYC,GCTAGCTA...,3000-3500;3600-4000,5.2,K562
mmu_circ_0001,Mus musculus,chr3:5000-6000,+,Trp53,TTTTAAAA...,5000-5500;5700-6000,8.7,MEF
```

## Architecture

The implementation follows the RNAcentral import pipeline architecture:

### Components

1. **parser.py**: Main parsing logic
   - Reads CSV files using Polars for efficient processing
   - Converts circRNA data to RNAcentral Entry format
   - Handles taxonomy lookups
   - Processes genomic coordinates and exon positions

2. **helpers.py**: Helper functions
   - Genomic location parsing
   - Exon position parsing
   - Taxonomy ID lookup with fallback mapping
   - URL and accession generation
   - Note data construction

3. **CLI (cli/circpedia.py)**: Command-line interface
   - Integrates with the `rnac` command system
   - Command: `rnac circpedia parse <taxonomy> <csv_file> <output>`

4. **Nextflow (workflows/databases/circpedia.nf)**: Pipeline integration
   - Fetches data from configured source
   - Runs parser with taxonomy context
   - Outputs standardized CSV files for database loading

## Usage

### Command Line

```bash
# Parse CIRCpedia CSV data
rnac circpedia parse context.db circpedia_v3.csv output_dir/ --assembly GRCh38
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
  remote = '/path/to/circpedia_v3.csv'
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

- **Polars**: Used for efficient CSV processing (preferred over pandas per requirements)
- **psycopg**: Would be used for database queries if needed
- **SqliteDict**: Used for taxonomy lookups

### Circular RNA Representation

- RNA type: SO:0000593 (circular RNA)
- Database: CIRCPEDIA
- Accessions: Generated using MD5 hash of ID + location for uniqueness
- URLs: Link to CIRCpedia V3 search interface

### Sequence Handling

If sequences are not provided in the CSV:
- Placeholder sequences are used (marked with warning)
- In production, sequences should be extracted from genome assemblies
- Back-splicing junction sequences may need special handling

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
