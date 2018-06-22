###################################################################
#
# Caller script for the presence/absence ecoli pathotype workflow
# 
# Author: Milan Patel
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

_PATHOTYPES = {
    'STEC': "Shiga toxin-producing Escherichia coli",
    'ETEC': "Enterotoxigenic Escherichia coli",
    'EPEC': "Enteropathogenic Escherichia coli",
    'EPEC (typical)': "Enteropathogenic Escherichia coli (typical)",
    'EAEC': "Enteroaggregative Escherichia coli",
    'EIEC/Shigella': "Enteroinvasive Escherichia coli",
    'DAEC': "Diffusely adherent Escherichia coli",
    'STEC/EAEC': "Hybrid shiga toxin-producing / enteroaggregative Escherichia coli",
    }

def sequence_parser(header, sequence, sep = ':'):
    # etec.fasta has has a pipe in the references file
    # that's why we have to do this
    parts = header.split('|')[0].split(sep)

    while len(parts) < 4:
        parts.append('')

    return SequenceInfo(
        locus = parts[0],
        allele = parts[1],
        accession = parts[2],
        sequence = sequence,
        other = ''
    )

def main(settings, env):

    log_message('Starting running presence/absence ecoli pathotype algorithm')

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

    # Load the resistance sequences:
    log_message('Loading pathotype sequences and associated'
        ' information...')

    # Load it
    sequence_database = DbInfo(
        database_path, seq_parser = sequence_parser)

    # Loading the sequences
    log_message('Successfully loaded sequences')

    # We were successful in running the algorithm
    log_message('Running pathotype algorithm...')

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
    results_out = sequence_database.results_parser(results)

    final_pathotypes = {
        'results' : set(),
        'extras' : []
    }

    present_genes = [gene for gene, boolean in \
        results_out['results'].iteritems() if boolean]

    final_pathotypes['results'].update(get_pathotypes(present_genes))

    for typ in final_pathotypes['results']:
        final_pathotypes['extras'].extend(_PATHOTYPES[typ] for typ in final_pathotypes['results'])

    final_pathotypes['results'] = list(final_pathotypes['results'])

    # Write the results out
    log_message('Writing results out...')

    write_results('ecoli.pathotype_genotypes.json', json.dumps(results_out))

    write_results('ecoli.pathotype_pathotypes.json', json.dumps(final_pathotypes))

    # Success!
    log_message('Successfully ran ecoli pathotype algorithm!')

def get_pathotypes(loci):

    if any(locus.startswith('ipaH') for locus in loci) or 'ipaD' in loci:
        return ['EIEC/Shigella']

    if any(locus.startswith('stx') for locus in loci) and \
        ('aaiC' in loci or 'aggR' in loci or 'aatA' in loci or 'aap' in loci):
        return ['STEC/EAEC']

    if any(locus.startswith('stx') for locus in loci):
        return ['STEC']

    if 'ltcA' in loci or 'sta1' in loci or 'stb' in loci:
        return ['ETEC']

    if 'eae' in loci:
        if ('bfpA' in loci or 'eaf' in loci):
            return ['EPEC (typical)']
        else:
            return ['EPEC']

    if 'aaiC' in loci or 'aggR' in loci or 'aatA' in loci or 'aap' in loci:
        return ['EAEC']

    if 'daaC' in loci:
        return ['DAEC']

    return []

