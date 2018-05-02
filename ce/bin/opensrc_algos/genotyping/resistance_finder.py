###################################################################
#
# Caller script for the presence/absence resistance workflow
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
    sequence_parser,
    notes_parser
)

from .ab_detection import (
    presence_detector,
)

import os
import json
import importlib
from functools import partial
from collections import namedtuple, defaultdict

def main(settings, env):

    # Log the inital message
    log_message('Starting running presence/absence resistance algorithm', 1)

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
        database_path))

    # Make the partial function that will parse the resistance information
    # *** IMPORTANT ***
    resistance_seq_parser = partial(sequence_parser, sep='_')
    
    # Load the resistance sequences:
    log_message('Loading resistance sequences and associated'
        ' information...', 1)
    
    # Load it
    sequence_database = DbInfo(
        database_path, seq_parser = resistance_seq_parser)

    # Loading the sequences
    log_message('Successfully loaded sequences', 2)

    # We were successful in running the algorithm
    log_message('Running resistance algorithm...', 1)

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
    log_message('Writing results out...', 1)
    write_results('resistance.json', json.dumps(results_out))

    # Success!
    log_message('Successfully ran resistance algorithm!', 2)

    # If point finder is in the available modules for this organism
    # run it

    if hasattr(settings, 'mutation_finder'):
        mutation_finder_env = env.copy()
        env.localdir = os.path.join(
            os.path.dirname(env.localdir),
            'mutation_finder'
        )

        mutation_finder = importlib.import_module('.mutation_finder', __package__)

        settings.mutation_finder.query = settings.query
        mutation_finder.main(settings.mutation_finder, mutation_finder_env)

def results_parser(dbinfo, results):

    sequences = dbinfo.sequences
    notes = dbinfo.notes

    results_out = {
        'results': {},
        'extra': []
    }

    present = set()

    # For each of the result, get the information we need
    for result, geno_regions in results.iteritems():

        # The sequence information for this result
        sequence_info = sequences[result]
        
        # The locus of this sequence
        locus = sequence_info.locus

        # The allele of the found genotype, if there is one
        allele = sequence_info.allele
        
        # The specifc locus information
        if notes:
            notes_info = notes[locus]


        # Updating for better output for surveillance
        # gene_name = '_'.join([locus, allele])
        gene_name = allele

        # We have both the locus and the resistance conferred
        # results_out['results'][gene_name] = True
        present.add(gene_name)

        # This will be the extra information that
        # will get dumped into the entry logs
        # no need to make parsing the results
        # more difficult than it has to be
        hit_information = [
            {
                    'locus': locus,
                    'identity': geno_region.identity,
                    'coverage': geno_region.coverage,
                    'allele': allele,
                    'antibiotic' : notes_info.antibiotic,
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

    results_out['results'].update((gene_name, True) for gene_name in present)

    for sequence_info in dbinfo.sequences.itervalues():

        locus = sequence_info.locus
        allele = sequence_info.allele
        # gene_name = '_'.join([locus, allele])
        gene_name = locus
        
        if gene_name in present:
            continue

        results_out['results'][gene_name] = False

    return results_out
