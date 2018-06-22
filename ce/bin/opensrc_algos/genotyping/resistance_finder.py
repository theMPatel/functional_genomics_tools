###################################################################
#
# Caller script for the presence/absence resistance workflow
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

from tools.environment import (
    log_message,
    log_error,
    log_progress,
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
import re
import json
import importlib
import mutation_finder
from functools import partial
from collections import namedtuple, defaultdict


def main(settings, env):

    # Log the inital message
    log_message('Starting running presence/absence resistance algorithm')

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
        ' information...')
    
    # Load it
    sequence_database = DbInfo(
        database_path, seq_parser = resistance_seq_parser)

    # Loading the sequences
    log_message('Successfully loaded sequences')

    # We were successful in running the algorithm
    log_message('Running resistance algorithm...')

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
    results_out, notes_out = sequence_database.results_parser(results, f=results_parser)

    # Write the results out
    log_message('Writing results out...')

    write_results('resistance.json', json.dumps(results_out))
    
    # If point finder is in the available modules for this organism
    # run it
    mutation_finder_antibios = {
        'results' : {}, 
        'extra' : []
    }

    if hasattr(settings, 'mutation_finder'):
        mutation_finder_env = env.copy()
        env.localdir = os.path.join(
            os.path.dirname(env.localdir),
            'mutation_finder'
        )

        settings.mutation_finder.query = settings.query
        mutation_finder_antibios = mutation_finder.main(settings.mutation_finder, mutation_finder_env)

    # Update the antiobiotic information from pointfinder if there are results
    for key, value in mutation_finder_antibios['results'].iteritems():

        if key in notes_out['results']:
            notes_out['results'][key] = notes_out['results'][key] or value

        else:
            notes_out['results'][key] = value

    # Update the hit information, this will now include information
    # from both resfinder and pointfinder now.
    notes_out['extra'].extend(mutation_finder_antibios['extra'])
    
    write_results('resistance.antibios.json', json.dumps(notes_out))
    
    # Success!
    log_message('Successfully ran resistance algorithm!')

gene_parser = re.compile(r'(^[a-zA-Z]*)(?:-\w*)?(?:-[\d]+$)')
re_comp_type = type(gene_parser)
def parse_res_genes(string, parser=gene_parser):
    # The purpose of this function is to return a condensed
    # name for resistance genes in the format
    # xxxx-123123 where the numbers are used as a counter
    # in the resistance finder database

    if isinstance(parser, re_comp_type):
        match = parser.search(string)

        if match:
            return match.group(1)
        
        return string

    else:
        return string

def results_parser(dbinfo, results):

    sequences = dbinfo.sequences
    notes = dbinfo.notes

    results_out = {
        'results': {},
        'extra': []
    }

    notes_out = {
        'results' : {},
        'extra' : []
    }

    present_genes = set()
    present_antibios = set()

    # For each of the result, get the information we need
    for result, geno_regions in results.iteritems():

        # The sequence information for this result
        sequence_info = sequences[result]
        
        # The locus of this sequence
        locus = sequence_info.locus

        # The allele of the found genotype, if there is one
        allele = sequence_info.allele
        
        # The specific locus information
        antibiotic = ''
        if notes:
            notes_info = notes.get(locus, None)

            if notes_info is not None:
                antibiotic = notes_info.antibiotic

        if antibiotic:
            if isinstance(antibiotic, list):
                present_antibios.update(antibiotic)
            else:
                present_antibios.add(antibiotic)

        # Updating for better output for surveillance
        # gene_name = '_'.join([locus, allele])
        gene_name = parse_res_genes(locus)

        # We have both the locus and the resistance conferred
        # results_out['results'][gene_name] = True
        present_genes.add(gene_name)

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
                    'antibiotic' : antibiotic,
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
        notes_out['extra'].extend(hit_information)

    results_out['results'].update((gene_name, True) for gene_name in present_genes)
    notes_out['results'].update((antibio, True) for antibio in present_antibios)

    for sequence_info in dbinfo.sequences.itervalues():

        locus = sequence_info.locus
        allele = sequence_info.allele

        gene_name = parse_res_genes(locus)

        if gene_name in present_genes:
            continue

        results_out['results'][gene_name] = False

    for note in notes.itervalues():

        antibio = note.antibiotic

        if isinstance(antibio, list):
            for a in antibio:
                if a in present_antibios:
                    continue

                notes_out['results'][a] = False

        elif antibio in present_antibios:
            continue

        else:
            notes_out['results'][antibio] = False

    return results_out, notes_out
