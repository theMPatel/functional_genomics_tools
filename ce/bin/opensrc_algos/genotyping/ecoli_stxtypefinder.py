###################################################################
#
# Detects the various stx gene variants for Escherichia
#
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import json
from collections import namedtuple

from tools.environment import (
    log_algo_version,
    log_message,
    write_results,
    valid_dir,
)

from tools.tools import (
    parse_fasta,
    is_fasta
)

from tools.dbinfo import (
    DbInfo
)

from .rb_detection import (
    build_sequences
)

from tools.bowtie import (
    bowtie_index,
    paired_bowtie2,
)

from tools.samtools import (
    sam_view,
    bam_sort,
    pile_up_sam
)

STXTarget = namedtuple('STXTarget', [
    'locus',
    'allele',
    'accession',
    'sequence'
])

def sequence_parser(header, sequence, sep='_'):
    info = header.split(sep)

    while len(info) < 4:
        info.append('')

    return STXTarget(
        locus=info[0]+info[1],
        allele=info[2],
        accession=info[3],
        sequence=sequence
    )

def prepare_pileup(settings, env, dbinfo):

    log_message('Dumping reference sequences...')

    out_file = os.path.join(env.localdir, 'bowtie', 'stxrefs.fasta')
    
    dbinfo.export_sequences(out_file)

    log_message('Success!')

    log_message('Indexing reference files...')

    index_file = bowtie_index(out_file, env)

    log_message('Success!')

    log_message('Mapping reads against references...')

    sam_file = paired_bowtie2(settings.query_reads, env, index_path=index_file)

    log_message('Success!')

    log_message('Sorting sam file and creating pileup...')

    bam_path = sam_view(
        sam_file,
        env,
        # Output to bam
        '-b',

        # Input is SAM format
        '-S'
    )

    bam_sorted = bam_sort(bam_path, env)

    pileup_path = pile_up_sam(bam_sorted, out_file, env)

    log_message('Success!')

    return pileup_path

def main(settings, env):
    
    log_message('Beginning reads based STX detection algorithm')

    # Write the version number of the database and algorithm
    version_path = settings['version']

    # Set the initial version information
    log_algo_version(
        algo_version = None,
        settings = settings,
        env = env
    )

    # Check just to make sure we have reads
    if len(settings.query_reads) < 2:
        raise RuntimeError('No read files provided or improperly paired reads')

    # Get the database path
    database_path = env.get_sharedpath(settings.database)

    log_message('Database path found at: {}'.format(
        database_path))

    log_message('Loading stx gene sequences and associated'
        ' information')

    sequence_database = DbInfo(
        database_path, seq_parser = sequence_parser)

    log_message('Succesfully loaded sequences and metadata')
    log_message('Indexing references and executing mapping')

    pileup_path = prepare_pileup(settings, env, sequence_database)

    log_message('Searching alignments for stx subtypes')

    found_sequences = build_sequences(
        pileup_path, settings.min_ambiguity, settings.min_coverage, settings.max_complexity)

    results_out = sequence_database.results_parser(found_sequences, f=results_parser)

    log_message('Found {} unique subtypes!'.format(
        sum(results_out['results'].values())
        )
    )
    
    log_message('Outputting results..')

    write_results('ecoli.stxfinder_genotypes.json', json.dumps(results_out))

    log_message('Successfully ran STX detection algorithm!')

def results_parser(dbinfo, results):

    results_out = {
        'results': {},
        'extra': []
    }

    for info in dbinfo.sequences.itervalues():
        gene_name = info.locus
        results_out['results'][gene_name] = False

    for seq in results:
        reference_info = dbinfo.sequences.get(seq.ref.split('|')[0], None)
        
        if reference_info is None:
            raise RuntimeError('Missing reference')

        results_out['results'][reference_info.locus] = True

        results_out['extra'].append({
            'locus': reference_info.locus,
            'allele': reference_info.allele,
            'accession': reference_info.accession,
            'sequence': ''.join(seq.get_fragment()),
            'coverage': str(seq.coverage),
            'complexity': str(seq.complexity)
        })

    return results_out