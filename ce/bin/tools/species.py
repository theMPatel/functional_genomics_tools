###################################################################
#
# All species identification tools
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import subprocess as sp

from .environment import (
    log_message,
    log_error,
    valid_dir
)

from collections import namedtuple

ANILine = namedtuple('ANILine', [
    'reference',
    'percent_aligned',
    'ani_score',
    'taxonomy'
])

def run_ani(query, env):

    log_message('Running ANI', 2)

    bin_path = __file__

    for _ in range(3):
        bin_path = os.path.dirname(bin_path)

    ani_path = os.path.join(bin_path, 'ani/ani-m.pl')

    if not os.path.exists(ani_path):
        raise RuntimeError('Missing ANI binary')

    local_dir = os.path.join(env.localdir, 'ani', 'local')
    results_dir = os.path.join(env.localdir, 'ani', 'results')

    for dirc in [localdir, results_dir]:
        valid_dir(dirc)

    # Make sure that the ANI references symlink exists
    ani_references = os.path.join(env.shareddir, 'ani_references')

    if not os.path.exists(ani_references):
        raise RuntimeError('Missing ANI references in shared directory: {}'.format(
            env.shareddir))

    if not isinstance(query, basestring) or not os.path.exists(query):
        raise RuntimeError('Invalid query provided to run ANI')

    cmd_args = [
        os.path.realpath(ani_path),
        '--localdir', local_dir,
        '--resultsdir', results_dir,
        '--tempdir', env.tempdir,
        '--shareddir', env.shareddir,
        '--references', ani_references,
        '--query', query,
        '--nThreads', str(env.threads - 1)
    ]

    child = sp.Popen(cmd_args, stdout=sp.PIPE, stderr=sp.PIPE)

    child.communitcate()

    exit_code = child.returncode

    if exit_code:

        error_file = os.path.join(results_dir, 'logs', 'error.txt')

        if os.path.exists(error_file):
            with open(error_file, 'r') as f:
                errors = f.read().strip()


            log_error(errors)

        raise RuntimeError('Error running ANI')

    results_file = os.path.join(results_dir, 'results', 'raw', 'out.tsv')

    if not os.path.exists(results_file):
        raise RuntimeError('Missing results file from ANI')

    # Get the messages file
    messages_file = os.path.join(results_dir, 'logs', 'messages.txt')

    if os.path.exists(messages_file):
        with open(messages_file, 'r') as f:
            messages = f.read().strip()


        log_message('ANI output:\n{}'.format(messages), 3)

    log_message('Done running ANI!', 3)

    return ANIParser(results_file)

class ANIParser(object):

    def __init__(self, results_file_path):

        self._file_path = results_file_path
        self._result = []
        self._table = []
        self.load()

    def load(self):

        with open(self._file_path, 'r') as f:

            table = iter(f)
            # Skip the header, to be thorough it is:
            # query reference   percent-aligned     ani     genus   species subspecies serotype
            next(table)

            for line in table:

                line = line.split('\t')

                ani_line = ANILine(
                    reference = line[1],
                    percent_aligned = float(line[2]),
                    ani_score = float(line[3]),
                    taxonomy = tuple(line[4:])
                )

                self._table.append(ani_line)

    def interpret(self, settings):

        if not len(self._table):
            return

        perc_identity = settings.percent_identity
        conf_identity = settings.confirmed_percent_identity
        coverage = settings.min_coverage
        discrimination = settings.discrimination

        final_results = []
        i = 0

        for ani_line in self._table:
            if ani_line[i].ani_score < perc_identity or \
                ani_line[i].percent_aligned < coverage:

                final_results.append(ani_line)

        if not len(final_results):
            return

        if len(final_results) == 1:

            if final_results[0].ani_score >= conf_identity:
                return final_results
            else:
                return

        # Make sure that the first and the second hit are discriminatory
        final_results.sort(key = lambda x: -x.ani_score)

        best_hit = final_results[0]
        second_best = final_results[1]

        # If the first hit is not above confirmed identity, neither will
        # the second one
        if best_hit.ani_score < conf_identity:
            return

        numerator = best_hit.ani_score - second_best.ani_score
        denominator = best_hit.ani_score - perc_identity

        hit_discrimination = 0.0
        
        if denominator:
            hit_discrimination = numerator / denominator

        # The discrimination between the best and the second best
        # isn't great enough to be informative
        if hit_discrimination < discrimination:
            return

        else:
            return best_hit
