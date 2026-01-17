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
- `parser.py` - Main CSV parsing logic using Polars
- `helpers.py` - Helper functions for data processing
- `README.md` - Comprehensive documentation
- `example_data.csv` - Sample data format

**Key Features**:
- Parses CSV files with circular RNA data
- Converts to RNAcentral Entry format
- Handles taxonomy lookups with fallback mappings
- Processes genomic coordinates and exon positions
- Uses **Polars** for efficient dataframe operations (as requested)
- Generates proper RNA type (SO:0000593 for circular RNA)

### 2. CLI Integration
**Location**: `rnacentral_pipeline/cli/circpedia.py`

Command:
```bash
rnac circpedia parse <taxonomy> <csv_file> <output> [--assembly ASSEMBLY]
```

Registered in: `rnacentral_pipeline/cli/__init__.py`

### 3. Nextflow Workflow
**Location**: `workflows/databases/circpedia.nf`

Workflow processes:
- `fetch_data`: Downloads CIRCpedia data from configured source
- `parse_data`: Runs parser with taxonomy context
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
  remote = '/nfs/production/agb/rnacentral/provided-data/circpedia/circpedia_v3.csv'
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
CIRCpedia CSV
    ↓
Nextflow: fetch_data
    ↓
Nextflow: parse_data
    ↓
Python: parser.parse()
    ↓
For each row:
  - Parse genomic location
  - Lookup taxonomy ID
  - Parse exon positions
  - Create Entry object
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

## Expected CSV Format

### Required Columns
- `circid` - CIRCpedia circular RNA ID
- `species` - Species name (e.g., "Homo sapiens")
- `location` - Genomic location (e.g., "chr1:12345-67890")

### Optional Columns
- `strand` - DNA strand ('+' or '-')
- `gene` - Host gene symbol
- `sequence` - RNA sequence
- `exon_positions` - Exon coordinates
- `fpm` - Expression level
- `cell_line` - Cell line/tissue
- `sequencing_type` - Sequencing type
- `conservation` - Conservation info
- `rnase_r_enrichment` - RNase R enrichment

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

### New Files (11)
```
rnacentral_pipeline/databases/circpedia/
  ├── __init__.py
  ├── parser.py
  ├── helpers.py
  ├── README.md
  └── example_data.csv

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
   cp rnacentral_pipeline/databases/circpedia/example_data.csv test_data.csv
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
   rnac circpedia parse context.db test_data.csv output/ --assembly GRCh38
   ```

### For Production

1. **Enable in config**:
   Edit `config/databases.config`:
   ```groovy
   circpedia {
     run = true  // Enable
     remote = '/path/to/circpedia_v3.csv'  // Set data path
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
- **Accessions**: Generated using MD5 hash of circID + location for uniqueness
- **URLs**: Point to CIRCpedia V3 search interface

### Sequence Handling
- If sequences provided in CSV: used directly
- If sequences missing: placeholder used with warning
  - Production should extract from genome assemblies
  - Special handling may be needed for back-splice junctions

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
