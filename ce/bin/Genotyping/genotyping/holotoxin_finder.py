###################################################################
#
# Assembly based holotoxin finder
#
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import json
from functools import partial
from collections import namedtuple, defaultdict

from tools.environment import (
    log_message,
    log_error,
    log_progress,
    log_ephemeral,
    log_algo_version,
    write_results,
    full_path,
    valid_dir
)

from tools.dbinfo import (
    DbInfo,
    SequenceInfo
)

from tools.tools import (
    is_fasta,
    parse_fasta,
    popen
)

deletions = re.compile(r'-[0-9]+[ACGTNacgtn]+')
insertions = re.compile(r'\+[0-9]+[ACGTNacgtn]+')

def bowtie_index(reference, env):

    if not isinstance(reference, basestring) or not \
        os.path.exists(reference):

        raise RuntimeError('Invalid reference file provided: {}'.foramt(
            str(reference)))

    reference_name = os.path.basename(reference).split('.')[0]
    index_dir = os.path.join(env.localdir, reference_name, 'index')

    valid_dir(index_dir)
    
    bowtie_index_path = os.path.join(env.toolsdir, 'bowtie2', 'bowtie2-build')

    if not os.path.exists(bowtie_index_path):
        raise RuntimeError('Missing Bowtie2 indexer')

    cmd_args = [
        os.path.abspath(os.path.realpath(sys.executable)),
        bowtie_index_path,
        reference,
        index_dir
    ]

    return_code, out, err = popen(cmd_args)

    if return_code:
        log_error(err.strip())
        raise RuntimeError('Error building bowtie2 index')

    return index_dir

def paired_bowtie2(read_files, env, index_path='', reference= ''):

    if not isinstance(read_files, list) or not \
        all(os.path.exists(x) for x in read_files):

        raise RuntimeError('Invalid read files provided')

    bowtie2_path = os.path.join(env.toolsdir, 'bowtie2', 'bowtie2')

    if not os.path.exists(bowtie2_path):
        raise RuntimeError('Missing bowtie2 executable')

    cmd_args = [
        os.path.realpath(bow)
    ]


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