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

    for dirc in [local_dir, results_dir]:
        valid_dir(dirc)

    # Make sure that the ANI references symlink exists
    ani_references = os.path.join(env.shareddir, 'ani_references', 'species.tsv')

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

    child = sp.Popen(
        cmd_args,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
        cwd=os.path.join(env.localdir, 'ani')
    )

    child.communicate()

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

    results = ANIParser()

    results.load(path=results_file)

    return results

class ANIParser(object):

    def __init__(self):
        self._result = []
        self._table = []

    def load(self, path=None, flobj=None):

        if flobj is not None:
            self._load(flobj)

        elif path is not None:
            with open(path, 'r') as flobj:
                self._load(flobj)

        else:
            raise RuntimeError('Nothing to load!')

    def _load(self, flobj):

        table = iter(flobj)
        # Skip the header, to be thorough it is:
        # query reference   percent-aligned     ani     genus   species subspecies serotype
        next(table)

        for line in table:

            line = line.strip().split('\t')

            ani_line = ANILine(
                reference = line[1],
                percent_aligned = float(line[2]),
                ani_score = float(line[3]),
                taxonomy = tuple(line[4:])
            )

            self._table.append(ani_line)

    def build_output_table(self, results):

        final_table = ['\t'.join(['Taxonomy', '% Aligned', 'ANI Score'])]

        for line in results:
            final_table.append('\t'.join(map(str, [
                line.taxonomy[2],
                line.percent_aligned,
                line.ani_score
            ])))

        return final_table

    def log_table(self, table, depth=2):

        for line in table:
            log_message(line, depth)

    def interpret(self, settings):

        if not len(self._table):
            return

        perc_identity = settings.percent_identity * 100
        coverage = settings.min_coverage * 100
        discrimination = settings.discrimination * 100

        final_results = []

        for ani_line in self._table:
            if ani_line.ani_score > perc_identity and \
                ani_line.percent_aligned > coverage:

                final_results.append(ani_line)

        if not len(final_results):
            return

        if len(final_results) == 1:
            return final_results[0]

        # Make sure that the first and the second hit are discriminatory
        final_results.sort(key = lambda x: -x.ani_score)

        log_message('ANI Results:', 2)
        pretty_table = self.build_output_table(final_results)
        self.log_table(pretty_table, depth=3)

        best_hit = final_results[0]
        second_best = final_results[1]

        if best_hit.ani_score >= second_best.ani_score + discrimination:
            return best_hit

        return

if __name__ == '__main__':

    trial = [
        ANILine(
            reference='IV_SalTDH2012K0845.fasta',
            percent_aligned=74.9,
            ani_score=95.08,
            taxonomy=('Salmonella', 'enterica', 'IV', '\n')
        ),

        ANILine(
            reference='II_62_3163_GLE7372A19.031473.circlator.quiver.031513.fasta',
            percent_aligned=73.6,
            ani_score=95.56,
            taxonomy=('Salmonella',
            'enterica',
            'II',
            '\n')),

        ANILine(
            reference='V_04_0440_GLE7360A3.031221.circlator.quiver.031490.fasta',
            percent_aligned=73.0,
            ani_score=90.3,
            taxonomy=('Salmonella',
            'bongori',
            '',
            '\n')),

        ANILine(
            reference='IIIa_2014K_1020_GLE7360A13.031225.circlator.quiver.031499.fasta',
            percent_aligned=73.25,
            ani_score=93.66,
            taxonomy=('Salmonella',
            'enterica',
            'IIIa',
            '\n')),

        ANILine(
            reference='I_SalJFX2010K2370c1Infantis.fa',
            percent_aligned=87.19,
            ani_score=98.8,
            taxonomy=('Salmonella',
            'enterica',
            'I',
            '\n')),

        ANILine(
            reference='VI_68_4603_GLE7360A6.031228.circlator.quiver.031493.fasta',
            percent_aligned=78.54,
            ani_score=95.84,
            taxonomy=('Salmonella',
            'enterica',
            'VI',
            '\n')),

        ANILine(
            reference='II_2011K_1440_GLE7372A2.031537.circlator.quiver.031618.fasta',
            percent_aligned=80.29,
            ani_score=96.25,
            taxonomy=('Salmonella',
            'enterica',
            'II',
            '\n')),

        ANILine(
            reference='II_08_0466_GLE7360A7.031230.canu.quiver.031612.fasta',
            percent_aligned=78.8,
            ani_score=94.85,
            taxonomy=('Salmonella',
            'enterica',
            'II',
            '\n')),

        ANILine(
            reference='IIIb_2015K_1072_GLE7367A18.031335.circlator.quiver.031597.fasta',
            percent_aligned=75.38,
            ani_score=95.5,
            taxonomy=('Salmonella',
            'enterica',
            'IIIb',
            '\n')),

        ANILine(
            reference='IV_64_2439_GLE7360A5.031226.circlator.quiver.031492.fasta',
            percent_aligned=72.97,
            ani_score=93.91,
            taxonomy=('Salmonella',
            'enterica',
            'IV',
            '\n'))
    ]

    import StringIO

    flobj = StringIO.StringIO(
        'query\treference\tpercent-aligned\tani\tgenus\tspecies\tsubspecies\tserotype\n')

    test_ani = ANIParser()
    test_ani.load(flobj=flobj)
    test_ani._table = trial

    settings = namedtuple('settings', [
        'percent_identity',
        'min_coverage',
        'discrimination'
    ])

    s = settings(
        percent_identity=0.8,
        min_coverage=0.7,
        discrimination=0.02
    )

    a = test_ani.interpret(s)
