process active {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
  path(query)
  path('template.txt')

  output:
  path('rnacentral_active.fasta.gz'), emit: 'active'
  path('example.txt'), emit: 'example'
  path('readme.txt'), emit: 'readme'

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" | json2fasta.py - rnacentral_active.fasta
  head rnacentral_active.fasta > example.txt
  gzip rnacentral_active.fasta
  cp template.txt readme.txt
  """
}

process inactive {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
  path(query)

  output:
  path('rnacentral_inactive.fasta.gz')

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" | json2fasta.py - - | gzip > rnacentral_inactive.fasta.gz
  """
}

process species_specific {
  input:
  path(query)

  output:
  path('rnacentral_species_specific_ids.fasta')

  script:
  """
  set -euo pipefail

  export PYTHONIOENCODING=utf8
  psql -q -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > raw.json
  json2fasta.py raw.json rnacentral_species_specific_ids.fasta
  """
}

process create_ssi {
  publishDir "${params.export.ftp.publish}/sequences/.internal/", mode: 'copy'
  memory '8GB'

  input:
  path('rnacentral_species_specific_ids.fasta')

  output:
  path('rnacentral_species_specific_ids.fasta.ssi')

  script:
  """
  esl-sfetch --index rnacentral_species_specific_ids.fasta
  """
}

process compress_species_fasta {
  publishDir "${params.export.ftp.publish}/sequences/", mode: 'copy'

  input:
  path('rnacentral_species_specific_ids.fasta')

  output:
  path('rnacentral_species_specific_ids.fasta.gz')

  script:
  """
  gzip < rnacentral_species_specific_ids.fasta > rnacentral_species_specific_ids.fasta.gz
  """
}


process find_dbs {
  input:
  path(query)

  output:
  path('dbs.txt')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > dbs.txt
  """
}

process database_specific {
  tag { db }
  maxForks params.export.ftp.sequences.by_database.max_forks
  publishDir "${params.export.ftp.publish}/sequences/by-database", mode: 'copy'

  input:
  tuple val(db), val(full_descr), path(query)

  output:
  file('*.fasta.gz')

  script:
  """
  export PYTHONIOENCODING=utf8
  psql -v ON_ERROR_STOP=1 -f "$query" -v db='%${full_descr}%' "\$PGDATABASE" > raw.json
  json2fasta.py raw.json ${db.toLowerCase().replaceAll(' ', '_').replace('/', '_')}.fasta
  gzip ${db.toLowerCase().replaceAll(' ', '_').replace('/', '_')}.fasta
  """
}

process by_species_readme {
  publishDir "${params.export.ftp.publish}/sequences/by-species/", mode: 'copy'

  input:
  path('template.txt')

  output:
  path('readme.txt')

  script:
  """
  cp template.txt readme.txt
  """
}

process find_species {
  input:
  path(query)

  output:
  path('species.csv')

  script:
  """
  psql -v ON_ERROR_STOP=1 -f "$query" "\$PGDATABASE" > species.csv
  """
}

process split_by_species {
  cpus 4
  memory '16 GB'
  time '12h'
  publishDir "${params.export.ftp.publish}/sequences/by-species", mode: 'copy'

  input:
  tuple path(species_csv), path(fasta)

  output:
  path('**/*.fasta.gz')

  script:
  """
  set -euo pipefail

  # Split the species-specific FASTA into one (uncompressed) file per species,
  # organised into kingdom subdirectories. Compression is done afterwards in
  # parallel, so this stays single-threaded and cheap.
  python3 - << 'PYEOF'
import csv, os, re, subprocess
from collections import Counter

KINGDOMS = [
    ('bacteria',      'Bacteria'),
    ('archaea',       'Archaea'),
    ('viruses',       'Viruses'),
    ('plants',        'Viridiplantae'),
    ('fungi',         'Fungi'),
    ('vertebrates',   'Vertebrata'),
    ('metazoa',       'Metazoa'),
    ('environmental', 'unclassified entries'),
]

def get_kingdom(lineage):
    for name, marker in KINGDOMS:
        if marker in lineage:
            return name
    return 'protists'

def safe_name(name):
    return re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')

# Load species from the single taxonomy query
rows = []
with open('${species_csv}') as f:
    for taxid, name, lineage in csv.reader(f):
        rows.append((taxid, get_kingdom(lineage), safe_name(name)))

# Detect (kingdom, name) collisions so we can disambiguate just those with
# the taxid, while leaving the common case as a clean species name.
collisions = {key for key, n in Counter((k, n) for _, k, n in rows).items() if n > 1}

species = {}
for taxid, k, name in rows:
    os.makedirs(k, exist_ok=True)
    fname = f'{name}_{taxid}' if (k, name) in collisions else name
    species[taxid] = (k, fname)

# Stream FASTA into sort, annotating each line with its record's taxid.
# Sequences are ordered by URS (not taxid), so we sort to group by taxid,
# allowing a single sequential pass with O(1) open file handles. Keep sort's
# scratch + memory in the work dir so it doesn't fill a small /tmp.
sort_proc = subprocess.Popen(
    ['sort', '-T', '.', '-S', '4G', '--stable', '-k1,1', '-t', '\t', '-o', 'sorted.tsv'],
    stdin=subprocess.PIPE, text=True
)
with open('${fasta}') as fin:
    taxid = ''
    for line in fin:
        if line.startswith('>'):
            taxid = line.split()[0].rsplit('_', 1)[-1]
        sort_proc.stdin.write(taxid + '\t' + line)
sort_proc.stdin.close()
if sort_proc.wait() != 0:
    raise SystemExit('sort failed')

# Single sequential pass: open one file at a time, switching on taxid change.
current_taxid = None
current_handle = None

with open('sorted.tsv') as fin:
    for raw in fin:
        taxid, line = raw.split('\t', 1)
        if taxid != current_taxid:
            if current_handle:
                current_handle.close()
                current_handle = None
            current_taxid = taxid
            if taxid in species:
                k, name = species[taxid]
                current_handle = open(f'{k}/{name}.fasta', 'w')
        if current_handle:
            current_handle.write(line)

if current_handle:
    current_handle.close()

os.unlink('sorted.tsv')
PYEOF

  # Compress every per-species file in parallel with stock gzip.
  find . -mindepth 2 -name '*.fasta' -print0 | xargs -0 -r -P ${task.cpus} gzip
  """
}

workflow fasta_export {
  channel.fromPath('files/ftp-export/sequences/active.sql') | set { active_sql }
  channel.fromPath('files/ftp-export/sequences/readme.txt') | set { readme }
  active(active_sql, readme)

  channel.fromPath('files/ftp-export/sequences/inactive.sql') | inactive
  channel.fromPath('files/ftp-export/sequences/species-specific.sql') \
  | species_specific \
  | (compress_species_fasta & create_ssi)

  channel.fromPath('files/ftp-export/sequences/databases.sql') \
  | find_dbs \
  | splitCsv \
  | combine(channel.fromPath('files/ftp-export/sequences/database-specific.sql')) \
  | database_specific

  channel.fromPath('files/ftp-export/sequences/species.sql') \
  | find_species \
  | combine(species_specific.out) \
  | split_by_species

  channel.fromPath('files/ftp-export/sequences/by-species-readme.txt') | by_species_readme
}
