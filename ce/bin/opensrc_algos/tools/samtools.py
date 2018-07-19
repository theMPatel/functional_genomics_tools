###################################################################
#
# Samtools runtime interface
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
    log_error,
    log_warning,
    valid_dir,
    full_path
)

class SamtoolsFilter(object):
    """
    Use this class to create filter objects for sam files to pass to
    sam view.

    These are the default alignment flags:

    0x1     PAIRED  paired-end (or multiple-segment) sequencing technology
    0x2     PROPER_PAIR each segment properly aligned according to the aligner
    0x4     UNMAP   segment unmapped
    0x8     MUNMAP  next segment in the template unmapped
    0x10    REVERSE SEQ is reverse complemented
    0x20    MREVERSE    SEQ of the next segment in the template is reverse complemented
    0x40    READ1   the first segment in the template
    0x80    READ2   the last segment in the template
    0x100   SECONDARY   secondary alignment
    0x200   QCFAIL  not passing quality controls
    0x400   DUP PCR or optical duplicate
    0x800   SUPPLEMENTARY   supplementary alignment
    """

    flags = {
        'paired'        : 0x1,
        'proper_pair'   : 0x2,
        'unmapped'      : 0x4,
        'mate_unmapped' : 0x8,
        'reversed'      : 0x10,
        'mate_reversed' : 0x20,
        'read_1'        : 0x40,
        'read_2'        : 0x80,
        'secondary'     : 0x100,
        'qc_fail'       : 0x200,
        'duplicate'     : 0x400,
        'supplementary' : 0x800
    }

    def __init__(self, filter_flags=None, require_flags=None):
        self._require_flag = 0x0
        self._filter_flag = 0x0

        if filter_flags is not None:
            if isinstance(filter_flags, (set, frozenset)):
                filter_flags = list(filter_flags)

            elif not isinstance(filter_flags, list):
                raise TypeError('Filter flags must be list type!')

            self.filter(*filter_flags)

        if require_flags is not None:
            if isinstance(require_flags, (set, frozenset)):
                require_flags = list(require_flags)

            elif not isinstance(require_flags, list):
                raise TypeError('Require flags must be list type!')

            self.require(*require_flags)

    def require(self, *args):
        for flag in set(args):

            if flag not in SamtoolsFilter.flags:
                raise RuntimeError('Inappropriate flag '
                    'provided: {}'.format(flag))

            self._require_flag |= SamtoolsFilter.flags[flag]

    def filter(self, *args):
        for flag in set(args):

            if flag not in SamtoolsFilter.flags:
                raise RuntimeError('Inappropriate flag '
                    'provided: {}'.format(flag))

            self._filter_flag |= SamtoolsFilter.flags[flag]

    def reset(self, require=False, flter=False):
        if require:
            self._require_flag = 0x0

        if flter:
            self._filter_flag = 0x0

    @property
    def require_flag(self):
        return str(self._require_flag)
    
    @property
    def filter_flag(self):
        return str(self._filter_flag)
    
def sam_view(sam_file, env, kwargs=None, *args):
    """
    This function takes care of the input/output for SAM file
    just provide the args and keyword args you want. Returns
    the filepath of the resulting file
    """

    if not os.path.exists(sam_file):
        raise RuntimeError('Invalid file provided for SAMTools view')

    # Run samtools view to get the bam file
    out_dir, sam_name = os.path.split(sam_file)
    new_name = sam_name.split('.')[0] + '.out'

    out_path = full_path(os.path.join(out_dir, new_name))

    samtools_path = full_path(
        os.path.join(
            env.toolsdir,
            'samtools'
        )
    )

    if not os.path.exists(samtools_path):
        raise RuntimeError('Missing samtools executable')

    cmd_args = [samtools_path, 'view']

    # TODO: Replace with new functions that validate args/kwargs and add
    #       to cmd_args in-place

    for arg in args:
        if not isinstance(arg, basestring) or not \
            arg.startswith('-'):

            log_warning('Invalid argument provided for sam view: '
                '{}, skipping..'.format(arg))

            continue

        cmd_args.append(arg)

    for arg, value in kwargs.iteritems():

        if not isinstance(arg, basestring) or not \
            arg.startswith('-'):

            log_warning('Invalid argument provided for sam view:'
                '{} -> {}, skipping...'.format(arg, value))

            continue

        cmd_args.extend([arg, value])
    
    cmd_args += ['-o', out_path, sam_file]

    log_message('Running samtools view args: {}'.format(
        ' '.join(cmd_args)))

    return_code, out, err = popen(cmd_args)

    if return_code:
        log_error(err)
        raise RuntimeError('Error running samtools view command')

    if not os.path.exists(out_path):
        raise RuntimeError('Missing bam results file')

    return out_path

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
            'all_tools',
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
            'all_tools',
            'samtools'
        )
    )

    if not os.path.exists(samtools_path):
        raise RuntimeError('Missing samtools executable')

    default_filter = SamtoolsFilter()
    
    default_filter.filter(
        'unmapped',
        'qc_fail',
        'duplicate'
    )

    cmd_args = [
        samtools_path,
        'mpileup',
        '-aa',
        '-f', reference,
        # This will push through secondary alignments
        '--ff', default_filter.filter_flag,
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