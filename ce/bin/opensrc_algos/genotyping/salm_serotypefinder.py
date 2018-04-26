###################################################################
#
# Caller script for SeqSero
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

# Set the version here. This should be the seq_sero version number!

__version__ = '1.0'

from tools.environment import (
    log_algo_version,
    log_message,
    log_error,
    log_progress,
    log_ephemeral,
    write_results,
    valid_dir
)

from tools.species import (
    run_ani
)

from tools.fancy_tools import (
    DecisionTree
)

import os
import sys
import csv
import json
import subprocess as sp
from functools import partial
from itertools import combinations
from collections import namedtuple, defaultdict

import pdb

def validate(**kwargs):

    if 'query_reads' in kwargs:

        if len(kwargs['query_reads']) != 2:
            raise ValueError('No read files provided')

        for file in kwargs['query_reads']:
            if not os.path.exists(file):
                raise IOError(
                    'File not found: {}'.format(file))

def set_environ_vars(env):

    BWA_PATH = os.path.join(env.toolsdir, 'bwa')

    if not os.path.exists(BWA_PATH):
        raise RuntimeError('Missing BWA, a SeqSero dependency')

    SAMTOOLS_PATH = os.path.join(env.toolsdir, 'samtools')

    if not os.path.exists(SAMTOOLS_PATH):
        raise RuntimeError('Missing Samtools, a SeqSero dependency')


    os.environ['bwa'] = BWA_PATH
    os.environ['samtools'] = SAMTOOLS_PATH

    hts_lib_path = os.path.join(
        os.path.dirname(
            os.path.dirname(
                os.readlink(os.environ['samtools'])
            ),
        ),
        'lib'
    )

    os.environ['LD_LIBRARY_PATH'] = '{}:{}'.format(
        os.environ['LD_LIBRARY_PATH'],
        hts_lib_path
    )

def parse_results(local_path):
    log_message('Retrieving SeqSero results...', 2)

    results = {
        'formula' : '',
        'serotype': ''
    }
    results_path = [local_path]

    for folder in os.listdir(local_path):

        if folder.startswith('SeqSero_result_'):
            results_path.extend([folder, 'Seqsero_result.txt'])
            break

    if len(results_path) != 3 or not \
        os.path.exists(os.path.join(*results_path)):

        raise RuntimeError('Could not find results file: {}'.format(str(results_path)))

    try:
        with open(os.path.join(*results_path), 'r') as f:
            content = f.read()

        for line in content.split('\n'):

            if line.startswith('Predicted antigenic profile:'):
                results['formula'] = line.split('\t')[1]

    except:
        raise RuntimeError('Could not read SeqSero results file...')

    return results

def assembly_run_seqsero(settings, env):
    # 
    # Seq_sero options:
    #   -m <int> (input data type, 
    #       '1' for interleaved paired-end reads ,
    #       '2' for separated paired-end reads, 
    #       '3' for single reads,
    #       '4' for genome assembly) 
    # 
    #   -i <file> (/path/to/input/file) 
    # 
    #   -b <string> (algorithms for bwa mapping; 
    #                       'mem' for mem, 
    #                       'sam' for samse/sampe;
    #                       default=sam; optional)
    # 

    raise NotImplementedError('Assembly based SeqSero not ready yet')

    # log_message('Running SeqSero {}'.format(__version__), 2)

    # local_dir = os.path.join(env.localdir, 'seq_sero')

    # valid_dir(local_dir)

    # cmd_args = [
    #     'python',
    #     os.path.join(env.toolsdir, 'SeqSero', 'SeqSero.py'),
    #     '-m', str(4),
    #     '-i', settings.query
    # ]

    # child = sp.Popen(cmd_args, cwd=local_dir, stdout=sp.PIPE, stderr=sp.PIPE)

    # stdout, stderr = child.communicate()

    # exit_code = child.returncode

    # log_ephemeral(stdout.stip())

    # if exit_code:
    #     log_error(stderr.strip())
    #     raise RuntimeError('Error running SeqSero')

    # log_message('Done running SeqSero', 3)

    # return local_dir

def reads_run_seqsero(settings, env):
    # 
    # Seq_sero options:
    #   -m <int> (input data type, 
    #       '1' for interleaved paired-end reads ,
    #       '2' for separated paired-end reads, 
    #       '3' for single reads,
    #       '4' for genome assembly) 
    # 
    #   -i <file> (/path/to/input/file) 
    # 
    #   -b <string> (algorithms for bwa mapping; 
    #                       'mem' for mem, 
    #                       'sam' for samse/sampe;
    #                       default=sam; optional)
    # 

    log_message('Running SeqSero {}'.format(__version__), 2)

    local_dir = os.path.join(env.localdir, 'seq_sero')

    valid_dir(local_dir)

    # The -m option below is for paired-end reads
    cmd_args = [
        sys.executable,
        os.path.join(env.toolsdir, 'SeqSero', 'SeqSero.py'),
        '-m', str(2),
        '-i', settings.query_reads[0],
        settings.query_reads[1]
    ]

    child = sp.Popen(cmd_args, cwd=local_dir, stdout=sp.PIPE, stderr=sp.PIPE)

    stdout, stderr = child.communicate()

    exit_code = child.returncode

    log_message(stdout.strip())
    log_message(stderr.strip())

    if exit_code:
        log_error(stderr.strip())
        raise RuntimeError('Error running SeqSero')

    log_message('Done running SeqSero', 3)

    return local_dir

def interpret_insilicopcr(env, ssp_antigenic, sslookup):
    
    insilicopcr_path = os.path.join(env.resultsdir, 'insilicopcr.json')

    if not os.path.exists(insilicopcr_path):
        return None

    else:
        with open(insilicopcr_path, 'r') as f:
            pcr_results = json.load(f)

    branch = {}
    # There are 6 typi so you need to split on the - to flatten
    for key, value in pcr_results['results'].iteritems():
        new_key = key.split('-')[0]

        if new_key in branch:
            branch[new_key] = branch[new_key] or value

        else:
            branch[new_key] = value

    branch = [[key, value] for key, value in branch.iteritems()]
    branch.sort(key=lambda leg: leg[0])

    serotype = sslookup.contingency_table[ssp_antigenic].recurse(branch)

    return serotype

def interpret_results(results, sslookup, settings, env):

    formula = results['formula']
    serotype = ''

    if __debug__:
        pdb.set_trace()

    if settings.ani_value is None:
        
        ani_results = run_ani(settings.query, env)
        best = ani_results.interpret(settings)
        
        if best is not None:
            ani_value = best.taxonomy[2].strip()

        else:
            return results

    results['formula'] = '{spp} {formula}'.format(
        spp=ani_value,
        formula=formula
    )

    if sslookup.lookup_table.get(ani_value, False) and \
        sslookup.lookup_table[ani_value].get(formula, False):

        serotype = sslookup.lookup_table[ani_value][formula]

    else:
        serotype = 'Needs further review'

    if serotype != 'to genotyper':
        results['serotype'] = serotype

    else:
        results['serotype'] = serotype

    return results

    # For now let's just get basic serotpying working before
    # we try to do insilicopcr
    ssp_antigenic = ' '.join([ani_value, antigenic_f])
    serotype = interpret_insilicopcr(env, ssp_antigenic, sslookup)

    if serotype is None:
        return results

    else:
        results['serotype'] = serotype

    return results

def main(settings, env):

    log_message('Starting running CDC-SeqSero serotype prediction algorithm', 1)

    log_algo_version(
        algo_version = __version__,
        settings = settings,
        env = env
    )

    # Double check to make sure that we have read files
    validate(
        query_reads=settings.query_reads
    )

    # Get the path of the lookup table
    database_path = env.get_sharedpath(settings.database)

    # Create the lookup table object
    sslookup = SeqSeroLookup(database_path)

    # Set the environment variables for SeqSero:
    # set_environ_vars(env)

    # Run Seq Sero
    local_dir = reads_run_seqsero(settings, env)

    # Parse the results
    results = parse_results(local_dir)

    # Interpret the results and update them if need be
    results_out = interpret_results(results, sslookup, settings, env)

    results_out = {
        'results': results,
        'extra': []
    }

    write_results('salmonella.serotype.json', json.dumps(results_out))

class SeqSeroLookup(object):
    # This is going to be a fancy wrapper for a couple 
    # dictionaries

    def __init__(self, db_path):

        if not os.path.exists(db_path) or not os.path.isdir(db_path):
            raise RuntimeError('Not a valid path for Salmonella serotyping')

        self.lookup_table = defaultdict(dict)
        self.contingency_table = {}
        self.unresolved_formulas = set()
        self._db_path = db_path

        self.load()

    def load(self):

        lookup_table_files = [
            'lookup_table.csv',
            'contingency_table.csv'
        ]

        for file in lookup_table_files:
            
            file_path = os.path.join(self._db_path, file)

            if not os.path.exists(file_path):
                raise RuntimeError('Missing: {}'.format(file))

        lookup_path = os.path.join(self._db_path, 'lookup_table.csv')

        with open(lookup_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)

            for row in reader:

                subspecies = row[0]
                seq_sero_form = row[1]
                final_serotype = row[3]

                self.lookup_table[subspecies][seq_sero_form] = final_serotype
                self.unresolved_formulas.add(seq_sero_form)

        contingency_path = os.path.join(self._db_path, 'contingency_table.csv')

        #self.build_contingency(contingency_path)

    def build_contingency(self, contingency_path):

        temp_table = defaultdict(list)

        with open(contingency_path, 'r') as f:
            reader = csv.DictReader(f)

            for row in reader:

                new_dict = { key:value for key, value in row.iteritems() \
                    if value in ('+', '-') or key == 'BN Serotype' }

                temp_table[row['BN Formula']].append(new_dict)

        for key, value in temp_table.iteritems():

            headers = set()
            headers.update(x for dct in value for x in dct.keys())
            headers = list(headers - set(['BN Report']))

            for i in range(len(headers)):

                headers[i].extend(branch[headers[i][0]] == '+' \
                    for branch in value)

            headers.sort(key=lambda header: header[0])

            headers.append(['BN Report'])
            headers[-1].extend(branch[headers[-1][0]] for branch in value)

            self.contingency_table[key] = DecisionTree()

            for i in range(1, len(headers[0])):
                new_branch = [[sub_branch[0], sub_branch[i]] for sub_branch in headers]
                self.contingency_table[key].update(new_branch)

