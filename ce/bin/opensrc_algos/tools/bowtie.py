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
        os.path.join(env.toolsdir, 'bowtie2_dir', 'bowtie2-build')
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
        os.path.join(env.toolsdir, 'bowtie2_dir', 'bowtie2')
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
