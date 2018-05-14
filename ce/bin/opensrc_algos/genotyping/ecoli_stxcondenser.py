###################################################################
#
# Caller script for the stx prediction algorithm
# 
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import re
import json

from tools.environment import (
    log_message,
    log_error,
    log_progress,
    log_ephemeral,
    log_algo_version,
    write_results
)

flnames = {

    'stxfinder' : 'ecoli.stxfinder_genotypes.json',
    'virulence' : 'virulence.json',
    'insilicopcr' : 'insilicopcr.json',
    'pathotyper' : 'ecoli.pathotype_genotypes.json'
}

stx_matcher = re.compile(r'^stx(?P<xtype>\d)(?:[AB])?(?P<subtype>[a-f])?$')

def main(settings, env):

    # We have all the results calculated since this is designed to
    # be the last genotyper run in the set
    results_dir = env.resultsdir

    # Lets parse the data for the experiment type first
    expr_out = {
        'results': {},
        'extra': []
    }

    fld_out = {
        'results' : [],
        'extra' : []
    }

    # stx types
    needed_chars = {
        'eaeA' : False,
        'ehxA' : False,
        'stx1' : False,
        'stx2' : False
    }

    # stx subtypes
    subtypes = {
        'stx1': range(ord('a'), ord('f')+1),
        'stx2': range(ord('a'), ord('f')+1)
    }

    stx_subtypes = {}

    # Set up the needed characters for the subtypes
    # easier than typing it out ;)
    for gene, r in subtypes.iteritems():
        for c in r:
            stx_subtypes[gene+chr(c)] = False

    # Coalesce each of the file names if they are present
    for file in flnames.itervalues():

        flpath = os.path.join(results_dir, file)

        if not os.path.exists(flpath):
            continue

        with open(flpath, 'r') as f:
            mod_results = json.load(f)

        mod_brief_results = mod_results.get('results', {})

        if not mod_brief_results:
            continue

        for char, present in mod_brief_results.iteritems():

            if not present:
                continue

            match_obj = stx_matcher.search(char)

            if not match_obj:
                
                if char in needed_chars:
                    needed_chars[char] = present

                continue

            xtype = match_obj.group('xtype')
            subtype = match_obj.group('subtype')
            final_type = 'stx'

            if xtype and final_type + xtype in needed_chars:
                final_type += xtype
                needed_chars[final_type] = present

            if subtype and final_type + subtype in stx_subtypes:
                final_type += subtype

                stx_subtypes[final_type] = present

    expr_out['results'].update(needed_chars.iteritems())
    
    for key, value in stx_subtypes.iteritems():

        if value:

            fld_out['results'].append(key)

    write_results('stx_condenser_expr.json', json.dumps(expr_out))
    write_results('stx_condenser_flds.json', json.dumps(fld_out))