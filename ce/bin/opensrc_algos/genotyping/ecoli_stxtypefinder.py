###################################################################
#
# Detects the various stx gene variants for Escherichia
#
# Author: Milan Patel, with some direction from key mentors ;)
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
    time_now,
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
    build_consensus
)

from tools.bowtie import (
    bowtie_index,
    paired_bowtie2,
    sam_view,
    bam_sort,
    pile_up_sam
)

STXTarget = namedtuple('STXTarget', [
    'locus'
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

    log_message('Dumping reference sequences...', 3)

    out_file = os.path.join(env.localdir, 'bowtie', 'stxrefs.fasta')
    
    dbinfo.export_sequences(out_file)

    log_message('Success!', 4)

    log_message('Indexing reference files...', 3)

    index_file = bowtie_index(out_file, env)

    log_message('Success!', 4)

    log_message('Mapping reads against references...', 3)

    sam_file = paired_bowtie2(settings.query_reads, env, index_path=index_file)

    log_message('Success!', 4)

    log_message('Sorting sam file and creating pileup...', 3)

    bam_path = sam_view(sam_file, env)

    bam_sorted = bam_sort(bam_path, env)

    pileup_path = pile_up_sam(bam_sorted, out_file, env)

    log_message('Success!', 4)

    return pile_up_path

def main(settings, env):
    
    log_message('Beginning reads based STX detection algorithm', 1)

    # Write the version number of the database and algorithm
    version_path = settings['version']

    # Set the initial version information
    log_algo_version(
        algo_version = None,
        settings = settings,
        env = env
    )

    # Check just to make sure we have reads

    if not len(settings.query_reads) != 2:
        raise RuntimeError('No read files provided or improperly paired reads')

    # Get the database path
    database_path = env.get_sharedpath(settings.database)

    log_message('Database path found at: {}'.format(
        database_path), 1)

    log_message('Loading stx gene sequences and associated'
        ' information', 2)

    sequence_database = DbInfo(
        database_path, seq_parser = sequence_parser)

    log_message('Succesfully loaded sequences and metadata', 3)
    log_message('Indexing references and executing mapping'2)

    pileup_path = prepare_pileup(settings, env, dbinfo)

    log_message('Searching alignments for stx subtypes', 2)

    found_sequences = build_consensus(pileup_path, settings.min_ambiguity)
    results = [f for f in found_sequences if not f.complexity]

    results_out = sequence_database.results_parser(found_sequences, f=results_parser)

    log_message('Found {} unique subtypes!'.format(sum(results_out['results'].values()), 3))
    log_message('Outputting results..', 2)

    write_results(ecoli.stxfinder_genotypes.json, json.dumps(results_out))

    log_message('Successfully ran STX detection algorithm!', 1)

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

        results_out['results'][gene_name] = True

        extra.append([
            'locus': gene_name,
            'allele': allele,
            'accession': reference_info.accession,
            'sequence': ''.join(seq.get_fragment())
        ])

    return results_out