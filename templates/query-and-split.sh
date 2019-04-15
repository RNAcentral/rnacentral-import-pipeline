#!/usr/bin/env bash

set -euo pipefail

psql -v ON_ERROR_STOP=1 ${variables} -f "$query" "$PGDATABASE" > raw.json
json2fasta.py raw.json rnacentral.fasta
seqkit shuffle --two-pass rnacentral.fasta > shuffled.fasta
seqkit split --two-pass --by-size ${chunk_size} --out-dir 'parts/' shuffled.fasta