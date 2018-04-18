###################################################################
#
# Caller script for the presence/absence ecoli serotype workflow
# 
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

# Set the version here. It makes more sense to set it as a variable
# here because if you edit the file, you can edit the version at the
# the same time

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

from .ab_detection import (
    presence_detector
)

import json
from functools import partial
from collections import namedtuple, defaultdict

def sequence_parser(header, sequence, sep = '_'):
    parts = header.split(sep)

    while len(parts) < 4:
        parts.append('')

    return SequenceInfo(
        locus = parts[0],
        allele = parts[1],
        accession = parts[2],
        sequence = sequence,
        # The O/H types information
        other = parts[3].split('/')
    )


def main(settings, env):

    log_message('Starting running presence/absence ecoli serotype algorithm', 1)

    # Write the version number of the database and algorithm
    log_algo_version(
        algo_version = None,
        settings = settings,
        env = env
    )

    # Get the database path
    database_path = env.get_sharedpath(settings.database)

    # Log it
    log_message('Database path found at: {}'.format(
        database_path), 2)

    # Load the resistance sequences:
    log_message('Loading serotype sequences and associated'
        ' information...', 3)

    # Load it
    sequence_database = DbInfo(
        database_path, seq_parser = sequence_parser)

    # Loading the sequences
    log_message('Successfully loaded sequences', 3)

    # We were successful in running the algorithm
    log_message('Running serotype algorithm...', 2)

    # Run the presence detector
    results = presence_detector(
        sequence_database,
        settings.query,
        settings.cached_query,
        settings.percent_identity,
        settings.min_relative_coverage,
        settings.min_merge_overlap,
        settings.search_fragments,
        env
    )

    # Get the results we want to write
    results_out = sequence_database.results_parser(results, f=results_parser)

    # Write the results out
    log_message('Writing results out...', 2)
    write_results('ecoli.serotype.json', json.dumps(results_out))

    # Success!
    log_message('Successfully ran serotype algorithm!', 3)


def results_parser(dbinfo, results):

    sequences = dbinfo.sequences
        
    results_out = {
        'results': {
            'O': [],
            'H': []
        },
        'extra': []
    }
    
    for result, geno_regions in results.iteritems():

        # The sequence information for this result
        sequence_info = sequences[result]

        # The locus of this sequence
        locus = sequence_info.locus

        # The allele of the found genotype, if there is one
        allele = sequence_info.allele

        # the O/Htypes:
        otypes = sequence_info.other

        for otype in otypes:

            if otype.lower().startswith('o'):
                results_out['results']['O'].append(otype)

            elif otype.lower().startswith('h'):
                results_out['results']['H'].append(otype)

        # This will be the extra information that
        # will get dumped into the entry logs
        # no need to make parsing the results
        # more difficult than it has to be
        hit_information = [
            {
                    'locus': locus,
                    'identity': geno_region.identity,
                    'coverage': geno_region.coverage,
                    'allele' : allele,
                    'otypes' : otypes,
                    'hits': [
                                {
                                'contig_id': hit.query_id,
                                'query_start': hit.query_start,
                                'query_stop': hit.query_stop,
                                'reference_start': hit.reference_start,
                                'reference_stop': hit.reference_stop,
                                'full_match' : hit.full_match,
                                } for hit in geno_region.locations
                            ]
            } for geno_region in geno_regions
        ]

        # Add it to the results out
        results_out['extra'].extend(hit_information)

    return results_out