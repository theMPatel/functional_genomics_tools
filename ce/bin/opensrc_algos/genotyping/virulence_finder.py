###################################################################
#
# Caller script for the presence/absence ecoli virulence finder 
# workflow
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
    SequenceInfo,
    sequence_parser
)

from .ab_detection import (
    presence_detector
)

import json
from functools import partial
from collections import namedtuple, defaultdict

def main(settings, env):

    log_message('Starting running presence/absence virulence finder algorithm', 1)

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
    log_message('Loading virulence finder sequences and associated'
        ' information...', 3)


    virulence_seq_parser = partial(sequence_parser, sep=':')

    # Load it
    sequence_database = DbInfo(
        database_path, seq_parser = virulence_seq_parser)

    # Loading the sequences
    log_message('Successfully loaded sequences', 3)

    # We were successful in running the algorithm
    log_message('Running virulence finder algorithm...', 2)

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

    # Write the results out
    log_message('Writing results out...', 2)

    write_results('virulence.json', json.dumps(results_out))

    # Success!
    log_message('Successfully ran virulence finder algorithm!', 2)
