###################################################################
#
# Caller script for insilico pcr
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
    write_results,
    valid_dir,
    log_algo_version
)

from tools.dbinfo import (
    DbInfo,
    SequenceInfo
)

from tools.tools import (
    is_fasta,
    parse_fasta,
    check_mismatches,
    reverse_complement
)

from tools.align import (
    BLASTSettings,
    align_blast,
    align_blast_nodb,
    GenotypeHit,
    create_blastdb
)

import os
import json
from functools import partial
from itertools import combinations
from collections import namedtuple, defaultdict

PcrTarget = namedtuple('PcrTarget', [
    'primer_id',
    'locus',
    'forward_primer_id',
    'forward_sequence',
    'reverse_primer_id',
    'reverse_sequence',
    'length'
])

def sequence_parser(line):
    parts = line.split('\t')

    while len(parts) < 7:
        parts.append('')

    return PcrTarget(
        primer_id = parts[0],
        locus = parts[1],
        forward_primer_id = parts[2],
        forward_sequence = parts[3],
        reverse_primer_id = parts[4],
        reverse_sequence = parts[5],
        length = parts[6]
    )

def main(settings, env):

    # Log the initial message
    log_message('Starting running insilico PCR algorithm', 1)
    
    # Set the initial version information
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

    # Let's load it all
    log_message('Loading insilico PCR targets and associated'
        ' information', 3)

    sequence_database = DbInfo(
        database_path, seq_parser = sequence_parser)

    log_message('Successfully loaded PCR targets and associated'
        ' information', 3)

    log_message('Running insilico PCR...', 2)

    log_message('Exporting PCR target database...', 3)

    # Make the path
    reference_dir = os.path.join(env.localdir, 'blastdb')
    
    # Check to make sure that its a real dir
    valid_dir(reference_dir)

    reference_path = os.path.join(reference_dir, 'references.fasta')

    # Export the reference sequences for blast database creation
    sequence_database.export_sequences(reference_path)

    log_message('Successfully exported reference database...', 4)
    
    # Create the path to the blast database
    blast_db_path = os.path.join(env.localdir, 'blastdb', 'db.fasta')

    log_message('Creating blast database...', 3)

    # Create the blast database
    create_blastdb(reference_path, blast_db_path, env)
    
    # Log that we were successful
    log_message('Successfully created blast database!', 4)   

    blast_settings = BLASTSettings(
        task = 'blastn-short',
        identity = settings.percent_identity,
        relative_minlen = 0,
        absolute_minlen = 0,
        include_sequences = True
    )

    log_message('BLASTing query genome against PCR targets', 3)

    # Run the alignment
    results = align_blast(
        settings.query,
        blast_db_path,
        blast_settings,
        env
    )

    log_message('Successfully BLASTed query genome against PCR targets', 4)

    log_message('Searching for ideal primer pairs...', 3)

    interpretations = find_targets(
        sequence_database,
        settings.cached_query,
        results,
        settings.max_mismatches,
        settings.max_nonIUPAC,
        settings.max_length_deviation,
        settings.percent_identity,
        env
    )

    if not len(interpretations):
        return

    results = sequence_database.results_parser(interpretations)

    log_message('Writing results...', 3)

    write_results('insilicopcr.json', json.dumps(results))

    log_message('Success!', 4)

def find_targets(sequence_database, cached_query, blast_results, max_mismatches,
    max_nonIUPAC, max_length_deviation, percent_identity, env):
    
    regions = defaultdict(list)

    for hit in blast_results.hits:

        # Filter the hits that have more than zero insertions/deletions
        if hit.num_gap_opens > 0:
            continue

        contig = cached_query.get(hit.query_id, None)
        if contig is None:
            raise RuntimeError('Could not find contig information')

        real_reference_seq = sequence_database.sequences[hit.reference_id]

        if hit.forward:

            hit.query_start = hit.query_start - hit.reference_start
            hit.query_stop = hit.query_stop + hit.reference_len - hit.reference_stop - 1

            pre_n = 'N' * (-1 * min(0, hit.query_start))
            post_n = 'N' * (-1 * min(0, len(contig)-1 - hit.query_stop))
            
            # It might be the case that the primer is hanging off the beginning
            # of a contig XOR the end.
            hit.query_seq = pre_n + \
                contig[max(0, hit.query_start):min(hit.query_stop, len(contig)-1)+1] + \
                post_n

            hit.reference_start = 0
            hit.reference_stop = hit.reference_len - 1

            hit.reference_seq = real_reference_seq

            hit.num_mismatches = check_mismatches(
                hit.query_seq, hit.reference_seq)

            hit.identity = 1.0 - \
                (float(hit.num_mismatches) / float(hit.reference_len))

        else:

            hit.query_start = hit.query_start - \
                hit.reference_len + hit.reference_stop + 1

            hit.query_stop = hit.query_stop + hit.reference_start

            pre_n = 'N' * (-1 * min(0, len(contig)-1 - hit.query_stop))
            post_n = 'N' * (-1 * min(0, hit.query_start))

            new_q_seq = post_n + \
                contig[max(0, hit.query_start):min(hit.query_stop, len(contig)-1)+1] + \
                pre_n

            hit.query_seq = reverse_complement(new_q_seq)

            hit.reference_start = 0
            hit.reference_stop = hit.reference_len - 1
            hit.reference_seq = real_reference_seq

            hit.num_mismatches = check_mismatches(
                hit.query_seq, hit.reference_seq)

            hit.identity = 1.0 - \
                (float(hit.num_mismatches) / float(hit.reference_len))

        # We exported the reference database as a string of
        # primerid&fwdid/revid
        primer_id = hit.reference_id.split('&')[0]

        # If the primer pairs exist, they should be part of the same
        # region
        regions[primer_id].append(hit)

    log_message('Found {} potential primer pair regions'.format(
        len(regions)), 4)

    if not len(regions):
        return None

    log_message('Filtering hits...', 4)

    # Keep only the best hits
    best = {}

    for region, hits in regions.iteritems():

        # Bare minimum you need two hits because they are a pair
        if len(hits) < 2:
            continue

        to_remove = []

        for i, hit in enumerate(hits):

            # Trust me, this is the fastest way
            # to do this (in python)
            counts = {'A':0, 'T':0, 'C':0, 'G':0}
            for nuc in hit.query_seq:
                if nuc in counts:
                    counts[nuc] += 1
                else:
                    counts[nuc] = 1

            # Check to make sure that the amount of
            # non ACTG is less than our threshold
            nonIUPAC = len(hit.query_seq) - (counts['A'] + counts['C'] + \
                counts['T'] + counts['G'])

            if nonIUPAC > max_nonIUPAC or hit.num_mismatches > max_mismatches:
                to_remove.append(i)
                continue

        # remove the poor quality hits
        for i in reversed(to_remove):
            del hits[i]

        # if we have a pair then keep those
        if len(hits) >= 2:
            best[region] = hits

    log_message('After filtering, retained {} PCR target regions'.format(
        len(best)), 4)

    if not len(best):
        return None

    final_results = defaultdict(list)

    for target, best_hits in best.iteritems():

        target_info = sequence_database.targets[target]

        for hit1, hit2 in combinations(best_hits, 2):

            # The two hits should be on the same contig
            if hit1.query_id != hit2.query_id:
                continue

            # These should be different from one another
            if hit1.forward == hit2.forward:
                continue

            # Get the length of the PCR product
            start = min(hit1.query_start, hit2.query_start)
            stop = max(hit1.query_stop, hit2.query_stop)

            length =  float(stop - start + 1)
            target_length = float(target_info.length)

            # If the length of the amplicon is not within (+/-)
            # of the deviation, continue
            if not (1.0 - max_length_deviation) * target_length <= length <= \
                (1.0 + max_length_deviation) * target_length:

                continue

            contig = cached_query.get(hit1.query_id, '')

            if not contig:
                raise RuntimeError('This should not have happened'
                    ' this late in the game')

            # Start is the start of the earliest primer
            # while stop is the end of the other
            sequence = contig[start:stop+1]

            forward = None
            reverse = None

            if hit1.reference_id.split('&')[1] == target_info.forward_primer_id:
                forward = hit1
                reverse = hit2
            else:
                forward = hit2
                reverse = hit1

            result = {
                'primer_id' : target_info.primer_id,
                'locus': target_info.locus,
                'forward_primer_id': forward.reference_id.split('&')[1],
                'forward_sequence': forward.query_seq,
                'reverse_primer_id': reverse.reference_id.split('&')[1],
                'reverse_sequence': reverse.query_seq,
                'forward_mismatch': forward.num_mismatches,
                'reverse_mismatch': reverse.num_mismatches,
                'expected_len' : target_info.length,
                'actual_len': length,
                'contig' : forward.query_id,
                'amplicon' : sequence
            }

            final_results[target].append(result)

    log_message('Retained {} PCR targets after analysis'.format(
        len(final_results)), 4)

    return final_results

class DbInfo(DbInfo):

    @property
    def targets(self):
        return self._targets

    @property
    def sequences(self):
        return self._sequences

    def load_database(self, dirpath, seq_parser, note_parser):

        self._targets = {}
        
        sequence_counts = defaultdict(dict)
        allele_id_template = '{}&{}'

        primer_path = os.path.join(dirpath, 'primers.txt')

        if not os.path.exists(primer_path):
            raise RuntimeError('Missing primers to run insilico PCR!')

        with open(primer_path, 'r') as f:

            for line in f:

                line = line.strip()

                if not line or line[0] == '#':
                    continue

                primer_info = seq_parser(line)

                forward_id = allele_id_template.format(
                    primer_info.primer_id, primer_info.forward_primer_id)

                assert forward_id not in self._sequences

                self._sequences[forward_id] = primer_info.forward_sequence

                reverse_id = allele_id_template.format(
                    primer_info.primer_id, primer_info.reverse_primer_id)

                assert reverse_id not in self._sequences

                self._sequences[reverse_id] = primer_info.reverse_sequence

                self._targets[primer_info.primer_id] = primer_info

    def export_sequences(self, filepath):

        valid_dir(os.path.dirname(filepath))

        ostr = '>{}|{}\n{}\n'

        with open(filepath, 'w') as f:

            for primer_id, sequence in self._sequences.iteritems():

                f.write(
                
                    ostr.format(
                        primer_id,
                        str(len(sequence)),
                        sequence
                
                    )
                )

    def results_parser(self, results):

        results_out = {
            'results' : {},
            'extra': [result for targets in \
                results.itervalues() for result in targets]
        }

        for target in self.targets:

            # This is an update since we don't need to know whether it was
            # stx2d-3 or stx2d-1
            pretty_name = target.rsplit(
                '-',
                1
            )[0]

            results_out['results'][pretty_name] = target in results

        return results_out