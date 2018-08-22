===================================================================
RNAcentral Genomic Coordinates Data
===================================================================

This directory contains genomic coordinates for a subset of RNAcentral ids
where such mapping is available.

* Bed
Format description:
http://www.ensembl.org/info/website/upload/bed.html
http://genome.ucsc.edu/FAQ/FAQformat.html

* Gff2
Format description:
http://www.sanger.ac.uk/resources/software/gff/spec.html

* Gff3
Format description:
http://www.sequenceontology.org/gff3.shtml

* track_hub/
UCSC-style track hub description:
https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html

Track hub folder structure:
- genomes.txt [list of annotated genomes]
- hub.txt [track hub description]
- hg38 [human GRCh38 assembly]
-- rnacentral.BigBed [bigBed binary data file]
-- rnacentral.html []
-- trackDb.txt [track description]
