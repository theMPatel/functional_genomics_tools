# Functional genomics tools

Functional genomics tools actively used at the Centers for Disease Control and Prevention for foodborne surveillance tasks.

## Description
This suite of tools is used for detecting biological traits of interest in bacteria from whole genome sequence data. These tools are specifically built to operate as a job in a High Performance Compute environment.


## Dependencies

  1. Python 2.7
  2. [NCBI Blast 2.7.1](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  3. [MUMmer 3.x](https://github.com/garviz/MUMmer)
  4. [Samtools 1.8](https://github.com/samtools/samtools)
  5. [BWA 0.7.17](https://github.com/lh3/bwa)
  6. [BowTie 2.3.4.1](https://github.com/BenLangmead/bowtie2)
  7. [SeqSero 1.0](https://github.com/denglab/SeqSero)

#### Citations
Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]

Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

Zhang S, Yin Y, Jones MB, Zhang Z, Deatherage Kaiser BL, Dinsmore BA, Fitzgerald C, Fields PI, Deng X.
Salmonella serotype determination utilizing high-throughput genome sequencing data.
J Clin Microbiol. 2015 May;53(5):1685-92.PMID:25762776
