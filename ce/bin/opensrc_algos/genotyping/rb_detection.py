###################################################################
#
# Reads based detection tools
#
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import re
import os
import sys
from itertools import izip, izip_longest, product
from collections import defaultdict, namedtuple

from tools.tools import (
    counter
)

_deletions = re.compile(r'-([0-9]+)([ACGTNacgtn]+)')
_insertions = re.compile(r'\+([0-9]+)([ACGTNacgtn]+)')
_substitutions = re.compile(r'[ACGTNacgtn]')
_remove = re.compile(r'[$<>]')
_start_read = re.compile(r'[\^]')
_reference = re.compile(r'[.,]')
_asterisk = re.compile(r'[\*]')

class ConsensusSequence(object):

    def __init__(self):
        self._start = -1
        self._stop = -1
        self._seq = []
        self._count = -1
        self._ambiguous = {}
        self._ref = ''

    def initialize(self, ref, start, stop, nuc, count, ref_call):

        self._start = start
        self._stop = stop

        if not isinstance(nuc, ConsensusPosition):
            cp = ConsensusPosition(stop, nuc, count, ref_call)

        else:
            cp = nuc

        if cp.ambiguous:
            self._ambiguous[stop] = cp

        self._seq.append(cp)
        self._count = count
        self._ref = ref

    def add_nuc(self, ref, pos, nuc, count, ref_call):

        if self._start == -1:
            return self.initialize(ref, pos, pos, nuc, count, ref_call)

        assert self.ref == ref

        if not isinstance(nuc, ConsensusPosition):
            cp = ConsensusPosition(pos, nuc, count, ref_call)

        else:
            cp = nuc

        if not pos == self._stop+1:
            raise RuntimeError('Adding to consensus sequence needs to be contiguous')

        if cp.ambiguous:
            self._ambiguous[pos-1] = cp

        self._stop += 1
        self._seq.append(cp)
        self._count += count

    @staticmethod
    def merge(first, second):

        if type(first) != type(second):
            raise RuntimeError('Type mismatch: {} vs. {}'.format(
                type(first), type(second)))

        if first.start > second.start:
            return ConsensusSequence.merge(second, first)


        if first.start <= second.end <= first.end:
            raise RuntimeError('Cannot merge overlapping sequences')

        # Get the gap between the two sequences in the case
        # that we have only parts of the gene covered
        gap_count = max(0, second.start - first.end - 1)

        while gap_count:
            first.add_nuc(first.end+1, '-', 0)
            gap_count -= 1

        # Append the sequence to the other
        for i, nuc in enumerate(second.seq, second.start):
            first.add_nuc(i, nuc, 0)

        # Add the coverage value so we can
        # get the coverage later
        first.count += second.count

        return first

    def flatten(self):

        rm = []

        for k, v in self.ambiguous.iteritems():
            if not v.flatten():
                rm.append(k)

        for k in rm:
            del self.ambiguous[k]

        self.count = sum(p.count for p in self)

    def get_fragment(self, start=None, stop=None):

        if start is None:
            start = self.start

        if stop is None:
            stop = self.stop

        for i in range(start-1, stop):
            if i+1 in self.ambiguous:
                raise StopIteration

            if self.seq[i].nuc == '-':
                continue
            else:
                yield self.seq[i]

    @property
    def complexity(self):
        return len(self.ambiguous)

    @property
    def ambiguous(self):
        return self._ambiguous

    @property
    def ref(self):
        return self._ref

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def seq(self):
        return self._seq

    @property
    def count(self):
        return self._count

    @count.setter
    def count(self, val):
        if not isinstance(val, (int, float)):
            raise RuntimeError("Count must be int or float! Got"
                " {} instead".format(val))

    @property
    def coverage(self):
        return int(float(self.count) / float(self.stop - self.start + 1))

    def __iter__(self):
        return iter(self._seq)

    def __len__(self):
        return len(self._seq)

    def __getitem__(self, index):
        return self._seq[index]

class ConsensusPosition(object):

    ambiguity_thresh = 1.0 / 3.0

    def __init__(self, pos, nuc, count, ref_call):

        self._pos = -1
        self._nuc = None
        self._count = -1
        self._ref_call = ''
        self._ambiguous = False
        self._analyze = True

        self.initialize(pos, nuc, count, ref_call)

    def initialize(self, pos, nuc, count, ref_call):

        self.pos = pos
        self.count = count
        self.ref_call = ref_call

        if isinstance(nuc, list):
            self.ambiguous = True
            self.nuc = nuc

        elif isinstance(nuc, basestring):
            self.nuc = nuc.upper()

        else:
            raise ValueError('Nucleotide position must be of type string')

    def flatten(self):

        if not self.ambiguous or not self.analyze:
            return self.ambiguous

        total = float(len(self.nuc))
        to_keep = set()
        counts = counter(self.nuc)

        for k, c in counts.iteritems():

            if float(c) / total >= ConsensusPosition.ambiguity_thresh:
                to_keep.add(k)

        if len(to_keep) == 1:
            self.nuc = to_keep.pop()
            self.ambiguous = not self._nuc == self._ref_call
            self.count = counts[self.nuc]
            

        elif self.ref_call in to_keep:
            self.nuc = self.ref_call
            self.ambiguous = False
            self.count = counts[self.nuc]
        
        else:
            i = 0
            while i < len(self.nuc):
                if self.nuc[i] in to_keep:
                    i += 1
                else:
                    # remove the sequence
                    removed = self.nuc.pop(i)

                    # remove the count of nucleotides
                    self.count -= len(removed)

        return self.ambiguous

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, val):
        if isinstance(val, int):
            self._pos = val

        else:
            raise ValueError('Nucleotide position should be of type int;'
                ' got {} instead'.format(type(val)))

    @property
    def nuc(self):
        return self._nuc

    @nuc.setter
    def nuc(self, val):
        if isinstance(val, basestring):
            self._nuc = val

        elif isinstance(val, list):
            if not self.ambiguous:
                raise ValueError('You must specify that this position'
                    ' is an ambiguous position before setting ambiguous values')

            self._nuc = val

        else:
            raise ValueError('Nucleotide position must be string/list type;'
                ' got {} instead'.format(type(val)))

    @property
    def count(self):
        return self._count

    @count.setter
    def count(self, val):
        if isinstance(val, int):
            self._count = val

        else:
            raise ValueError('Nucleotide count must be of type int;'
                ' got {} instead'.format(type(val)))

    @property
    def ambiguous(self):
        return self._ambiguous

    @ambiguous.setter
    def ambiguous(self, val):

        if isinstance(val, bool):
            self._ambiguous = val

        else:
            raise ValueError('Ambiguity must be of type bool;'
                ' got {} instead'.format(type(val)))

    @property
    def analyze(self):
        return self._analyze

    @analyze.setter
    def analyze(self, val):

        if isinstance(val, bool):
            self._analyze = val
        else:
            raise ValueError('Analyze attribute must be of type bool;'
                ' got {} instead'.format(type(val)))

    @property
    def ref_call(self):
        return self._ref_call

    @ref_call.setter
    def ref_call(self, val):
        if isinstance(val, basestring):
            self._ref_call = val.upper()
        else:
            raise RuntimeError('Reference nucleotide must be string type'
                ' got {} instead'.format(type(val)))
    
def pileup_iterator(flname=None, flobj=None):

    if flobj:
        for line in _pileup_iterator(flobj):
            yield line

    elif flname:

        if not isinstance(flname, basestring) or not \
            os.path.exists(flname):

            raise RuntimeError('Invalid pileup file')

        with open(flname, 'r') as f:

            for line in _pileup_iterator(f):
                yield line

    else:
        raise RuntimeError('No fileobj or filename provided for pileup scanning')

def _pileup_iterator(flobj):

    for line in flobj:
        split_line = line.strip().split('\t')

        while len(line) < 7:
            split_line.append('')

        yield split_line

def base_action(read_calls, index, final, match):
    return match.end()

def delete_action(read_calls, index, final, match):

    # Just to cover our bases to make sure
    # that the logic called the proper
    # function
    assert read_calls[index] == '-'

    # Add the number of deletions to the last sequence
    # We will need it to know how much to expand the
    # window by
    final[-1] = (final[-1], int(match.group(1)))

    return match.end()

def insert_action(read_calls, index, final, match):

    # Just to cover our bases to make sure
    # that the logic called the proper
    # function
    assert read_calls[index] == '+'

    final[-1] += match.group(2).upper()

    return match.end()

def _process_iter_re(read_calls, match_obj, index=0, final=None):

    if final is None:
        final = []

    upto = match_obj.start()

    while index < upto:
        final.append(read_calls[index])
        index += 1

    return index, final

def process_iter_re(read_calls, matcher=None, match_objs=None, action=None):

    # the purpose of this function is to parse out the pile_up line
    # into discrete pieces:
    # E.g.
    #
    # ,  +1t  ,+1t  .+1T  .  .+1T  .+1T  .+1T  .  .+1T  .+1T
    # |    |             |  |
    # ^^^^^^             ^^^^
    #    |                  There are no insertions here
    # These two are associated
    #

    if action is None or not callable(action):
        action = base_action

    index = 0
    final = []

    if match_objs:

        for match in match_objs:
            index, final = _process_iter_re(read_calls, match, index=index, final=final)
            index = action(read_calls, index, final, match)

        while index < len(read_calls):
            final.append(read_calls[index])
            index += 1

    elif matcher:

        for match in matcher.finditer(read_calls):

            index, final = _process_iter_re(read_calls, match, index=index, final=final)
            index = action(read_calls, index, final, match)

        while index < len(read_calls):
            final.append(read_calls[index])
            index +=1

    return final

def process_line(line):

    # All the information associated to a line
    ref = line[0]
    position = int(line[1])
    reference_call = line[2]
    read_count = int(line[3])
    read_calls = line[4]
    quality_scores = map(lambda x: ord(x)-33, line[5])
    mapping_quality = map(lambda x: ord(x)-33, line[6])

    ins = list(_insertions.finditer(read_calls))
    dels = list(_deletions.finditer(read_calls))

    if ins or dels:
        # Insertions are very straight forward,
        # Let's assume that the nucleotide position is a kmer of k > x
        # where x is the length of the insertion. The read count
        # is going to be the number of positions where there is an
        # insertion

        # Haven't handled the case where there are insertions and deletions
        assert bool(ins) ^ bool(dels)

        if ins:
            extracted = process_iter_re(read_calls, match_objs=ins, action=insert_action)

            extracted_with_ref_calls = map(
                lambda x: _reference.sub(reference_call, x), extracted)

        else:
            extracted = process_iter_re(read_calls, match_objs=dels, action=delete_action)
            extracted_with_ref_calls = []
            
            for e in extracted:
                if isinstance(e, tuple):
                    extracted_with_ref_calls.append(
                        (e[0].replace('.', reference_call).replace(',', reference_call), e[1])
                    )

                else:
                    extracted_with_ref_calls.append(
                        (e.replace(',', reference_call).replace('.', reference_call), 0)
                    )

        return bool(dels), ref, position, read_count, extracted_with_ref_calls, reference_call

    read_calls_list = []

    # Remove the dollar signs and other things
    # that don't tell us things we actually want to know
    read_calls = _remove.sub(r'', read_calls)

    i = 0
    for match in _start_read.finditer(read_calls):

        if not match:
            continue

        group = match.group()

        if not group:
            continue

        if isinstance(group, tuple):
            raise RuntimeError("We were only matching one character!")

        upto = match.start()

        while i < upto:
            read_calls_list.append(read_calls[i])
            i+=1

        i += 2

    while i < len(read_calls):
        read_calls_list.append(read_calls[i])
        i += 1

    # Make sure this assertion passes, because otherwise
    # we did something wrong
    #assert len(read_calls_list) == read_count

    read_calls = ''.join(read_calls_list)
    read_calls = _asterisk.sub(r'-', read_calls)

    if _substitutions.search(read_calls) or '-' in read_calls:
        # Replace all of the calls that are forward or reverse
        # matches with the reference call
        read_calls = _reference.sub(reference_call, read_calls).upper()
        return False, ref, position, read_count, list(read_calls), reference_call

    else:
        # If there is no other call than the reference
        # we can just return the call
        return False, ref, position, read_count, reference_call, reference_call

def detect_deletion_window(pile, ref, position, read_count, read_calls,
    current_consensus, ref_call,):

    max_jump = max(read_calls, key=lambda x: x[1])[1]

    stack = []

    # In the weird off chance that the deletion is
    # at the end of the sequence but samtools  gives us a deletion
    # larger than the remainder of the data stream

    try:
        for _ in range(max_jump):
            stack.append(next(pile))
    except StopIteration:
        max_jump = len(stack)

    processed_stack = [process_line(line) for line in stack]
    all_lines = [read_calls]
    all_lines.extend(line[4] for line in processed_stack)

    #
    # all_lines represents a window: [
    #   POS 108 -> [nuc, nuc, ... , nuc, nuc, nuc]
    #   POS 109 -> [nuc, nuc, ... , nuc, nuc]
    #   POS 110 -> [nuc, nuc, ... , nuc]
    #   POS 111 -> [nuc, nuc, ... , nuc, nuc]
    # ]
    #
    #
    # The first list in all_lines is the read_calls this function was called with
    # so the first position prior to the deletion window
    # and the last list being the last position in the deletion window
    #

    sub_seq = []

    for window in izip_longest(*all_lines, fillvalue='-'):
        # If bowtie says that this is a deletion window
        # Then it better be consistent
        # Namely, this happens:
        #   POS 108 -> .-3TTA
        #   POS 109 -> A
        #   POS 110 -> *
        #   POS 111 -> T
        #
        # Clearly this is stupid, mpileup says there's a deletion here
        # but that is not reflected in subsequent reads
        #
        # The opposite will also happen where:
        #
        #   POS 108 ->  .
        #   POS 109 ->  *
        #   POS 110 ->  *
        #   POS 111 ->  *
        #
        # Clearly a deletion event of 3 was not flagged
        # and yet all reads afterwards are deleted

        if all(x=='-' for x in window[1:]):
            sub_seq.append(tuple([window[0][0]] + list(window[1:])))

        elif all(x!='-' for x in window[1:]):
            sub_seq.append(tuple([window[0][0]] + list(window[1:])))

    # We want to get the counts of the deletion window
    # meaning how many reads support a certain substring
    counts = counter(sub_seq)
    to_remove = []
    for seq, count in counts.iteritems():
        if float(count) / len(sub_seq) < ConsensusPosition.ambiguity_thresh:
            to_remove.append(seq)

    for rm in to_remove:
        del counts[rm]

    # If we reduce the complexity down
    # to just one substring, hooray!
    if len(counts) == 1:
        final_seqs = counts.keys()[0]
        count = counts.values()[0]
        current_consensus.add_nuc(ref, position, final_seqs[0], count, ref_call)


        for i, next_nuc in enumerate(final_seqs[1:]):
            current_consensus.add_nuc(
                ref,
                processed_stack[i][2],
                next_nuc,
                0 if next_nuc=='-' else count,
                processed_stack[i][5]
        )

        return 0

    else:

        i = len(sub_seq)-1
        while i >= 0:
            if sub_seq[i] not in counts:
                sub_seq.pop(i)
            i-=1

        final_seqs = map(list, zip(*sub_seq))

        first_cp = ConsensusPosition(
            position,
            final_seqs[0],
            sum(x!='-' for x in final_seqs[0]),
            processed_stack[i][5]
        )

        # We've already done the frequency analysis
        # there's no longer a need to do it
        #first_cp.analyze = False

        current_consensus.add_nuc(
            ref,
            first_cp.pos,
            first_cp,
            first_cp.count,
            ref_call
        )

        for i, next_nuc in enumerate(final_seqs[1:]):

            cp = ConsensusPosition(
                processed_stack[i][2],
                next_nuc,
                sum(x!='-' for x in next_nuc),
                processed_stack[i][5]
            )

            #cp.analyze = False

            current_consensus.add_nuc(
                processed_stack[i][1],
                cp.pos,
                cp,
                cp.count,
                processed_stack[i][5]
            )

        return 0

def build_consensus(pileup_file, ambiguity_threshold):

    reference = None
    current_consensus = ConsensusSequence()
    position_offset = 0
    pile = pileup_iterator(flname=pileup_file)
    ConsensusPosition.ambiguity_thresh = ambiguity_threshold
    
    while pile:

        next_line = next(pile)
        dels, ref, position, read_count, read_calls, reference_call = process_line(next_line)

        if position < current_consensus.stop:
            # We've started on the next reference
            # in the pileup file

            yield reference, current_consensus

            current_consensus = ConsensusSequence()
            position_offset = 0

        reference = ref

        if not read_count:
            current_consensus.add_nuc(reference, position+position_offset, '-', read_count, reference_call)
            continue

        if dels:
            position_offset += detect_deletion_window(
                pile,
                ref, position,
                read_count,
                read_calls,
                current_consensus,
                reference_call
            )
            continue

        current_consensus.add_nuc(reference, position+position_offset, read_calls, read_count, reference_call)
