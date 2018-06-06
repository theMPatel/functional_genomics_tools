###################################################################
#
# Bowtie runtime interface
#
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import subprocess as sp
from .tools import (
    popen,
    parse_paired_files
)
from .environment import (
    log_message,
    valid_dir,
    full_path
)

def bowtie_index(reference, env, name=''):

    if not isinstance(reference, basestring) or not \
        os.path.exists(reference):

        raise RuntimeError('Invalid reference file provided: {}'.foramt(
            str(reference)))

    if name:
        reference_name = name
    else:
        reference_name = os.path.basename(reference).split('.')[0]
    
    index_dir = full_path(
        os.path.join(env.localdir,'bowtie', reference_name, 'index')
    )

    valid_dir(index_dir)

    bowtie_index_path = full_path(
        os.path.join(env.toolsdir, 'bowtie2', 'bowtie2-build')
    )

    if not os.path.exists(bowtie_index_path):
        raise RuntimeError('Missing Bowtie2 indexer')

    cmd_args = [
        full_path(sys.executable),
        bowtie_index_path,
        reference,
        index_dir
    ]

    log_message('Running bowtie2 index args: {}'.format(
        ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args)
    
    if return_code:
        log_error(err.strip())
        raise RuntimeError('Error building bowtie2 index')

    return index_dir

def paired_bowtie2(read_files, env, index_path='', reference= ''):

    if not isinstance(read_files, list) or not \
        all(os.path.exists(x) for x in read_files):

        raise RuntimeError('Invalid read files provided')

    bowtie2_path = full_path(
        os.path.join(env.toolsdir, 'bowtie2', 'bowtie2')
    )

    if not os.path.exists(bowtie2_path):
        raise RuntimeError('Missing bowtie2 executable')

    if not index_path:
        if not reference:
            raise RuntimeError('Bowtie index or reference provided!')

        index_path = bowtie_index(reference, env)

    output_dir = full_path(
        os.path.join(env.localdir, 'bowtie')
    )

    output_filename = os.path.join(output_dir, 'output.sam')

    valid_dir(output_dir)

    # This function is expecting that the read files are in the following order:
    # read_files = [R1, R2]
    # Make sure that this is the case
    cmd_args = [
        bowtie2_path,
        '-x',
        index_path,
        '-p', str(env.threads-1),
        '--reorder',
        '--local',
        '--sensitive-local',
        '--no-unal',
        '--all',
        '-1', read_files[0],
        '-2', read_files[1],
        '-S', output_filename
    ]

    log_message('Running bowtie2 args: {}'.format(
        ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args, cwd=output_dir)

    if return_code:
        log_error(err.strip())
        raise RuntimeError('Error running bowtie2 alignment')

    if not os.path.exists(output_filename):
        raise RuntimeError('Missing output file from bowtie2')

    return output_filename

def sam_view(sam_file, env):

    if not os.path.exists(sam_file) or not \
        sam_file.endswith('.sam'):

        raise RuntimeError('Invalid file provided for bam conversion')

    # Run samtools view to get the bam file
    bam_dir, sam_name = os.path.split(sam_file)
    bam_name = sam_name.split('.')[0] + '.bam'

    bam_path = full_path(os.path.join(bam_dir, bam_name))

    samtools_path = full_path(
        os.path.join(
            env.toolsdir,
            'samtools'
        )
    )

    if not os.path.exists(samtools_path):
        raise RuntimeError('Missing samtools executable')

    cmd_args = [
        samtools_path,
        'view',
        '-b', '-S',
        '-o', bam_path,
        sam_file
    ]

    log_message('Running samtools view args: {}'.format(
        ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args)

    if return_code:
        log_error(err)
        raise RuntimeError('Error running samtools view command')

    if not os.path.exists(bam_path):
        raise RuntimeError('Missing bam results file')

    return bam_path

def bam_sort(bam_file, env):

    if not isinstance(bam_file, basestring) or not \
        os.path.exists(bam_file):
        raise RuntimeError('Invalid bam file provided')

    bam_dir, bam_name = os.path.split(bam_file)
    sorted_name = bam_name.split('.')[0] + '.sorted'

    bam_sorted_path = full_path(
        os.path.join(bam_dir, sorted_name))

    samtools_path = full_path(
        os.path.join(
            env.toolsdir,
            'samtools'
        )
    )

    if not os.path.exists(samtools_path):
        raise RuntimeError('Missing samtools executable')

    itermediate_name = bam_name.split('.')[0] + '_intermediate'

    cmd_args = [
        samtools_path,
        'sort',
        '-T', itermediate_name,
        '-o', bam_sorted_path,
        bam_file
    ]

    log_message('Running samtools sort args: {}'.format(
        ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args)

    if return_code:
        log_error(err.strip())
        raise RuntimeError('Error running samtools sort')

    if not os.path.exists(bam_sorted_path):
        raise RuntimeError('Missing sorted bam file')

    return bam_sorted_path

def pile_up_sam(bam_sorted, reference, env):

    if not isinstance(bam_sorted, basestring) or not \
        os.path.exists(bam_sorted):

        raise RuntimeError('Invalid sorted bam file')

    if not isinstance(reference, basestring) or not \
        os.path.exists(reference):

        raise RuntimeError('Invalid reference file')

    # We want to make sure that the sorted bam file has been indexed
    sorted_index_name = os.path.basename(bam_sorted).split('.')[0]

    index_path = bowtie_index(bam_sorted, env, name=sorted_index_name)
    parent_dir = index_path

    # This will be the bowtie2 dir that is created
    # in the local dir
    for _ in range(2):
        parent_dir = os.path.dirname(parent_dir)

    pileup_name = os.path.split(reference)[1].split('.')[0] + '.pup'

    samtools_path = full_path(
        os.path.join(
            env.toolsdir,
            'samtools'
        )
    )

    if not os.path.exists(samtools_path):
        raise RuntimeError('Missing samtools executable')

    cmd_args = [
        samtools_path,
        'mpileup',
        '-aa',
        '-f', reference,
        # This will push through secondary alignments
        '--ff', '1540',
        '-s', bam_sorted,
        '-o', pileup_name

    ]

    log_message('Running samtools mpileup args: {}'.format(
        ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args, cwd=parent_dir)

    if return_code:
        log_error(err.strip())
        raise RuntimeError('Error running mpileup')

    pileup_path = full_path(os.path.join(parent_dir, pileup_name))

    if not os.path.exists(pileup_path):
        raise RuntimeError('Missing pileup_path')

    return pileup_path
