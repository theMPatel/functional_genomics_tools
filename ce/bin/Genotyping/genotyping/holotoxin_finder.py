###################################################################
#
# Assembly based holotoxin finder
#
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

from tools.environment import (
    log_message,
    log_error,
    log_progress,
    log_ephemeral,
    log_algo_version,
    write_results
)

from tools.dbinfo import (
    DbInfo,
    SequenceInfo
)

from tools.tools import (
    is_fasta,
    parse_fasta
)

from .ab_detection import (
    mutation_detector
)

import os
import json
from functools import partial
from collections import namedtuple, defaultdict

"""
 1769  for f in ./*; do mkdir "${f%.*}" && python ~/tools/all_executables/bowtie2-build $f "${f%.*}/index"; done
 1770  ll
 1771  for f in ./*; do   if [ -f "$f" ]; then    bowtie2 -x "${f%.*}/index" -p 4 --reorder --local --sensitive-local --no-unal --all -1 ../2010C-4541-M947-14-038-Jun11_S5_L001_R1_001.fastq -2 ../2010C-4541-M947-14-038-Jun11_S5_L001_R2_001.fastq > "${f%.*}.sam";   fi; done
 1772  history
 1773  ll
 1774  mv stx1_consensus.bam stx1_consensus.sam
 1775  ll
 1776  samtools view -b -S -o ./stx1_consensus.bam ./stx1_consensus.sam
 1777  samtools sort stx1_consensus.bam stx1_consensus.sorted
 1778  samtools sort -T stx1inter -o stx1_consensus.sorted stx1_consensus.bam
 1779  samtools index stx1_consensus.sorted
 1780  samtools mpileup -aa -f holotoxins.fasta --ff 1540 -s h.sorted > h.pup
 """