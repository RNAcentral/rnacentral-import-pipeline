# CIRCpedia V3 Import Implementation

## Summary

This document summarizes the implementation of CIRCpedia V3 circular RNA data import for the RNAcentral pipeline.

## What is CIRCpedia V3?

CIRCpedia V3 is a comprehensive circular RNA database published in Nucleic Acids Research (2025):
- **2.6 million circular RNAs** from 20 species
- Expression data from **2350 NGS and 63 TGS datasets**
- Community-recommended nomenclature
- Comprehensive genomic annotations
- Available at: https://bits.fudan.edu.cn/circpediav3

**Reference**: Zhai SN, Zhang YY, Chen MH, Fu ZC, Chen LL, Ma XK, Yang L. DOI: 10.1093/nar/gkaf1039

## Implementation Overview

A complete import pipeline has been implemented following RNAcentral's architecture patterns, including:

### 1. Python Parser Module
**Location**: `rnacentral_pipeline/databases/circpedia/`

Files created:
- `__init__.py` - Module initialization
- `parser.py` - Main TSV/FASTA parsing logic using Polars
- `helpers.py` - Helper functions for data processing
- `README.md` - Comprehensive documentation
- `example_annotation.txt` - Sample TSV annotation file
- `example_sequences.fa` - Sample FASTA sequence file

**Key Features**:
- Parses TSV annotation files with circular RNA metadata
- Loads sequences from FASTA files
- Converts to RNAcentral Entry format
- Handles taxonomy lookups with fallback mappings
- Processes combined genomic location fields (chr:start-end(strand))
- Uses **Polars** for efficient dataframe operations (as requested)
- Generates proper RNA type (SO:0000593 for circular RNA)
- Uses CIRCpedia's own IDs with CIRCPEDIA: prefix for accessions
- Creates direct links to circRNA detail pages

### 2. CLI Integration
**Location**: `rnacentral_pipeline/cli/circpedia.py`

Command:
```bash
rnac circpedia parse <taxonomy> <annotation_file> <fasta_file> <output> [--assembly ASSEMBLY]
```

Registered in: `rnacentral_pipeline/cli/__init__.py`

### 3. Nextflow Workflow
**Location**: `workflows/databases/circpedia.nf`

Workflow processes:
- `fetch_data`: Downloads CIRCpedia annotation and FASTA files from configured sources
- `parse_data`: Runs parser with taxonomy context and both input files
- Outputs standard CSV files for RNAcentral loading

Integrated into: `workflows/parse-databases.nf`

### 4. Configuration
**Location**: `config/databases.config`

Configuration added:
```groovy
circpedia {
  run = false                    // Disabled by default
  needs_taxonomy = true          // Requires taxonomy database
  process.directives.memory = 8.GB
  remote {
    annotation = '/nfs/production/agb/rnacentral/provided-data/circpedia/circpedia_v3_annotation.txt'
    fasta = '/nfs/production/agb/rnacentral/provided-data/circpedia/circpedia_v3_sequences.fa'
  }
  assembly = 'GRCh38'
}
```

### 5. Unit Tests
**Location**: `tests/databases/circpedia/`

Files created:
- `__init__.py`
- `helpers_test.py` - Tests for helper functions (35 test cases)
- `parser_test.py` - Tests for parser logic (13 test cases)

**Test Design**:
- Uses minimal mock data
- No network dependencies
- No RNAcentral database dependencies
- No large test data files
- All tests can run in isolation

## Data Flow

```
CIRCpedia Annotation (TSV) + Sequences (FASTA)
    ↓
Nextflow: fetch_data
    ↓
Nextflow: parse_data
    ↓
Python: parser.parse()
    ↓
1. Load sequences from FASTA into dictionary
2. Read TSV with Polars
3. For each row:
   - Parse combined location field (chr:start-end(strand))
   - Lookup taxonomy ID
   - Get sequence from FASTA dictionary
   - Create Entry object with CIRCPEDIA: prefix
    ↓
EntryWriter
    ↓
Standard RNAcentral CSV files:
  - accessions.csv
  - seq_short.csv / seq_long.csv
  - regions.csv
  - references.csv
  - etc.
```

## Expected Data Format

### Annotation File (TSV)

#### Required Columns
- `circID` - CIRCpedia circular RNA ID
- `species` - Species name (e.g., "Homo sapiens")
- `Location` - Combined genomic location with strand (e.g., "V:15874634-15876408(-)")

#### Optional Columns
- `gene_Refseq` - RefSeq gene identifier
- `gene_Ensembl` - Ensembl gene identifier
- `circname` - Circular RNA name
- `length` - Length of circular RNA
- `subcell_location` - Subcellular localization
- `editing_site` - RNA editing sites
- `DIS3_signal` - DIS3 degradation signals
- `Orthology` - Orthology information
- `TGS` - Third-generation sequencing support
- `transcript_Ensembl` - Ensembl transcript ID
- `transcript_Refseq` - RefSeq transcript ID

### Sequence File (FASTA)

Standard FASTA format with headers matching circIDs:
```
>circID_1
ATCGATCG...
>circID_2
GCTAGCTA...
```

## Technology Stack

Following requirements:
- ✅ **Polars** for dataframe operations (instead of pandas)
- ✅ **psycopg** ready for database queries (if needed)
- ✅ Standard RNAcentral patterns (Entry, EntryWriter, etc.)
- ✅ Nextflow for pipeline orchestration
- ✅ Click for CLI interface

## Testing

Run tests:
```bash
# All circpedia tests
pytest tests/databases/circpedia/

# Specific test files
pytest tests/databases/circpedia/helpers_test.py
pytest tests/databases/circpedia/parser_test.py

# With coverage
pytest --cov=rnacentral_pipeline.databases.circpedia tests/databases/circpedia/
```

## Files Created/Modified

### New Files (12)
```
rnacentral_pipeline/databases/circpedia/
  ├── __init__.py
  ├── parser.py
  ├── helpers.py
  ├── README.md
  ├── example_annotation.txt
  └── example_sequences.fa

rnacentral_pipeline/cli/
  └── circpedia.py

workflows/databases/
  └── circpedia.nf

tests/databases/circpedia/
  ├── __init__.py
  ├── helpers_test.py
  └── parser_test.py

CIRCPEDIA_IMPLEMENTATION.md
```

### Modified Files (3)
```
rnacentral_pipeline/cli/__init__.py
  - Added circpedia import
  - Registered circpedia CLI command

workflows/parse-databases.nf
  - Added circpedia workflow include
  - Added circpedia to workflow mix

config/databases.config
  - Added circpedia configuration
```

## Usage Instructions

### For Development/Testing

1. **Prepare test data**:
   ```bash
   # Use example data
   cp rnacentral_pipeline/databases/circpedia/example_annotation.txt test_annotation.txt
   cp rnacentral_pipeline/databases/circpedia/example_sequences.fa test_sequences.fa
   ```

2. **Build taxonomy context** (if not already available):
   ```bash
   # This is normally done by the pipeline
   wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
   tar xzf new_taxdump.tar.gz
   mkdir taxdump
   mv *.dmp taxdump
   rnac context build taxdump context.db
   ```

3. **Parse data**:
   ```bash
   rnac circpedia parse context.db test_annotation.txt test_sequences.fa output/ --assembly GRCh38
   ```

### For Production

1. **Enable in config**:
   Edit `config/databases.config`:
   ```groovy
   circpedia {
     run = true  // Enable
     remote {
       annotation = '/path/to/circpedia_v3_annotation.txt'
       fasta = '/path/to/circpedia_v3_sequences.fa'
     }
   }
   ```

2. **Run pipeline**:
   ```bash
   nextflow run main.nf -profile standard
   ```

## Implementation Notes

### Circular RNA Specifics
- **RNA Type**: SO:0000593 (circular RNA from Sequence Ontology)
- **Database**: CIRCPEDIA (uppercase as per RNAcentral convention)
- **Accessions**: Use CIRCpedia's own IDs with CIRCPEDIA: prefix (e.g., "CIRCPEDIA:hsa_circ_0001_1:100-200")
- **URLs**: Direct links to circRNA detail pages (https://bits.fudan.edu.cn/circpediav3/circrna/{circ_id})
- **Coordinate System**: 1-based, fully-closed (same as GFF/GTF format)

### Sequence Handling
- Sequences loaded from separate FASTA file
- Sequences keyed by circID
- DNA sequences (T not U) - automatically converted if needed
- If circID missing from FASTA: entry skipped with warning
- Production may need special handling for back-splice junctions

### Species Support
Includes fallback taxonomy mapping for 20 common species:
- Human, Mouse, Rat, Zebrafish, Fruit fly, C. elegans
- Yeast, Arabidopsis, Chicken, Chimp, Macaque, Dog
- Cow, Pig, Xenopus, Medaka, Fugu, Sea urchin
- Sea squirt, Mosquito

### Performance Considerations
- Uses Polars for efficient CSV processing
- Memory: 8GB configured (adjustable based on data size)
- Batch processing with progress logging every 10,000 rows
- Taxonomy database: SqliteDict for fast lookups

## Future Enhancements

Potential production improvements:

1. **Sequence Extraction**
   - Integrate with genome assemblies
   - Extract sequences from coordinates
   - Handle back-splice junctions

2. **Expression Data**
   - Expand expression profile support
   - Store TPM/FPKM values
   - Link to expression atlas

3. **Alternative Splicing**
   - Track alternative back-splicing events
   - Multiple isoforms per circRNA

4. **Conservation**
   - Import conservation scores
   - Cross-species mappings

5. **Data Quality**
   - Validation checks
   - Duplicate detection
   - Coordinate verification

6. **Performance**
   - Chunked processing for large files
   - Parallel parsing
   - Optimized memory usage

## Sources

Research sources used:

- [CIRCpedia v3 Paper - NAR 2025](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaf1039/8296757)
- [CIRCpedia v2 Documentation](https://pmc.ncbi.nlm.nih.gov/articles/PMC6203687/)
- [circAtlas Database Information](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02018-y)

## License

Copyright [2009-2025] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0

---

**Implementation Status**: ✅ Complete (Local only - not pushed to remote as requested)

**Branch**: `claude/add-circpedia-v3-m3I0I`

**Ready for**: Code review, testing with real data, integration into pipeline
