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
    write_results,
    valid_dir
)

from tools.species import (
    run_ani
)

from tools.tools import (
    check_b64encoded
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

def parse_results(local_path, results):
    log_message('Retrieving SeqSero results...')

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

            log_message(line)

            if line.startswith('Predicted antigenic profile:'):
                results['formula'] = line.split('\t')[1]

    except:
        raise RuntimeError('Could not read SeqSero results file...')

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

    # log_message('Running SeqSero {}'.format(__version__))

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

    # log_message(stdout.stip())

    # if exit_code:
    #     log_error(stderr.strip())
    #     raise RuntimeError('Error running SeqSero')

    # log_message('Done running SeqSero')

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

    log_message('Running SeqSero {}'.format(__version__))

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

    if exit_code:
        log_error(stderr.strip())
        raise RuntimeError('Error running SeqSero')

    log_message('Done running SeqSero')

    return local_dir

def interpret_insilicopcr(env, ssp_antigenic, sslookup):
    
    # First make sure this is something that we have information for
    if not sslookup.c_table.get(ssp_antigenic, False):
        return None

    insilicopcr_path = os.path.join(env.resultsdir, 'insilicopcr.json')

    if not os.path.exists(insilicopcr_path):
        return None

    else:
        with open(insilicopcr_path, 'r') as f:
            data = f.read()

        if check_b64encoded(data):
            pcr_results = json.loads(base64.b64decode(data))
        else:
            pcr_results = json.loads(data)

    lines = sslookup.c_table[ssp_antigenic]
    condensed_view = pcr_results.get('results', {})
    
    if not condensed_view:
        return None
    
    serotype = ''

    for line in lines:
        if serotype:
            break

        sslkp_arr = [k for k in line.iterkeys() if k != 'BN Serotype']
        sslkp_arr.sort()

        sslkp_vals = [line[k] for k in sslkp_arr]
        pcr_vals = [condensed_view[k] for k in sslkp_arr]

        found = True

        for lkp_val, pcr_val in zip(sslkp_vals, pcr_vals):

            if lkp_val == pcr_val:
                continue
            else:
                found = False
                break

        if found:
            serotype = line.get('BN Serotype')

    return serotype

def interpret_results(results, sslookup, settings, env):

    formula = results['formula']
    ani_value = settings.ani_value
    serotype = ''

    if ani_value == 'S. bongori':
        results['serotype'] = 'Needs further review'
        return results
    
    if __debug__:
        pdb.set_trace()

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

        # For now let's just get basic serotpying working before
        # we try to do insilicopcr
        serotype = interpret_insilicopcr(env, results['formula'], sslookup)

        if not serotype:
            results['serotype'] = 'Needs further review'
            return results

        else:
            results['serotype'] = serotype

    return results

def main(settings, env):

    log_message('Starting running CDC-SeqSero serotype prediction algorithm')

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

    results = {
        'formula' : '',
        'serotype' : ''
    }

    if not settings.ani_value:
        ani_results = run_ani(settings.query, env)
        best = ani_results.interpret(settings)

        if best is not None:
            settings.ani_value = best.taxonomy[2].strip()

        else:
            results['serotype'] = 'Needs further review'

    if settings.ani_value is not None:

        # Run Seq Sero
        local_dir = reads_run_seqsero(settings, env)

        # Parse the results
        parse_results(local_dir, results)

        # Interpret the results and update them if need be
        interpret_results(results, sslookup, settings, env)

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

        lookup_table_files = {
            'lookup' : 'lookup_table.csv',
            'contingency' : 'contingency_table.tsv'
        }

        for file in lookup_table_files.itervalues():
            
            file_path = os.path.join(self._db_path, file)

            if not os.path.exists(file_path):
                raise RuntimeError('Missing: {}'.format(file))

        lookup_path = os.path.join(self._db_path, lookup_table_files['lookup'])

        with open(lookup_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)

            for row in reader:

                subspecies = row[0]
                seq_sero_form = row[1]
                final_serotype = row[3]

                self.lookup_table[subspecies][seq_sero_form] = final_serotype
                self.unresolved_formulas.add(seq_sero_form)

        contingency_path = os.path.join(self._db_path, lookup_table_files['contingency'])

        self.build_contingency(contingency_path)

    def build_contingency(self, contingency_path):

        self._table = defaultdict(list)

        with open(contingency_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')

            for row in reader:

                new_dict = { key:value for key, value in row.iteritems() \
                    if value in ('+', '-') or key == 'BN Serotype' }

                for key, value in new_dict.iteritems():
                    if key == 'BN Serotype':
                        continue

                    new_dict[key] = value == '+'

                self._table[row['BN Formula']].append(new_dict)

    @property
    def c_table(self):
        return self._table
    
