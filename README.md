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

## Usage

### HPC Implementation

Before delving into the usage of this module, you need to first understand the directory heirarchy:

```
          ce
          ├── bin
          │   └── opensrc_algos
          │       ├── config.json
          │       ├── cewrapper
          │       │   ├── __init__.py
          |       |   |
          |       |   ....
          │       ├── genotyping
          │       │   ├── init.py
          |       |   |
          |       |   ....
          │       └── tools
          │           ├── init.py
          |           |
          |           ....
          ...
```

The opensrc_algos contains all of the **_Python_** modules that can be run as a job. Execution flow:
  1. User requests certain job which is handled by _cewrapper_;
  2. _cewrapper_ translates the job name into specific modules within the opensrc_algos package and executes;
  3. All job specific execution is handled by the modules in charge, which return with/without error to the _cewrapper_;
  4. _cewrapper_ notifies user of job success or failure;
  
### Commandline Usage

These tools were not developed with the user in mind, meaning it is not necessarily straight forward to run this as a user sitting at a terminal. That said, with software **_ANYTHING IS POSSIBLE!_**

#### Basic argument requirements:

##### cewrapper.py

The cewrapper needs a handful of basic arguments passed to the commandline (or as a file) in order to properly set up logging and the job environment. An important thing to note is that extra arguments that are not defined below can and _should_ be passed along to the cewrapper. These extra arguments are passed along to the job being called. **These are the bare minimum required arguments.**

```
        [--nThreads]        The number of threads to use for a requested job
        [--localdir]        A local working directory
        [--tempdir]         A temporary shared directory
        [--resultsdir]      A place for a job to put its results
        [--shareddir]       Datafiles shared by multiple jobs
        [--toolsdirs]       A great place for external dependecies
        [--algorithm]       The job to load and execute
```

#### genotyping.py

This is the handler script for the functional genomics pipeline. It will load specific modules and execute different parts of the functional genotyping algorithm depending on the genus of the organism being analyzed.

```
        [--organism]        Genus of the species needed to be analyzed
        [--query]           A path to an assembled genome to be analyzed
        [--query-reads]     Un-assembled read files to be analyzed
```

**NOTE**: _You do not call this module yourself, you let the cewrapper set up the environment and do the module calling_

### Citations

I would like to thank the below scientists, programmers, colleagues in arms for building incredible tools that save the world on a daily basis. Incredible work should never go unrecognized!

Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]

Li H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943]

Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. Nature Methods. 2012, 9:357-359.

Zhang S, Yin Y, Jones MB, Zhang Z, Deatherage Kaiser BL, Dinsmore BA, Fitzgerald C, Fields PI, Deng X.
Salmonella serotype determination utilizing high-throughput genome sequencing data.
J Clin Microbiol. 2015 May;53(5):1685-92.PMID:25762776

Identification of acquired antimicrobial resistance genes.
Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV.
J Antimicrob Chemother. 2012 Jul 10.
PMID: 22782487         doi: 10.1093/jac/dks261

Real-time whole-genome sequencing for routine typing, surveillance, and outbreak detection of verotoxigenic Escherichia coli.
Joensen KG, Scheutz F, Lund O, Hasman H, Kaas RS, Nielsen EM, Aarestrup FM.
J. Clin. Micobiol. 2014. 52(5): 1501-1510.

PlasmidFinder and pMLST: in silico detection and typing of plasmids.
Carattoli A, Zankari E, Garcia-Fernandez A, Voldby Larsen M, Lund O, Villa L, Aarestrup FM, Hasman H.
Antimicrob. Agents Chemother. 2014. April 28th.

Joensen, K. G., A. M. Tetzschner, A. Iguchi, F. M. Aarestrup, and F. Scheutz. 2015. Rapid and easy in silico serotyping of Escherichia coli using whole genome sequencing (WGS) data. J.Clin.Microbiol. 53(8):2410-2426. doi:JCM.00008-15 [pii];10.1128/JCM.00008-15 [doi]
