###################################################################
#
# Read based holotoxin finder
#
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import re
import os
import sys
from copy import deepcopy
from collections import defaultdict, namedtuple, Counter

_deletions = re.compile(r'-([0-9]+)([ACGTNacgtn]+)')
_insertions = re.compile(r'\+([0-9]+)([ACGTNacgtn]+)')
_substitutions = re.compile(r'[ACGTNacgtn]')
_remove = re.compile(r'[$<>\*]')
_start_read = re.compile(r'[\^]')
_reference = re.compile(r'[.,]')

class ConsensusSequence(object):

    def __init__(self):
        self._start = -1
        self._stop = -1
        self._seq = []
        self._count = -1
        self._ambiguous = {}
        self._ref = ''

    def initialize(self, ref, start, stop, nuc, count):

        self._start = start
        self._stop = stop

        if not isinstance(nuc, ConsensusPosition):
            cp = ConsensusPosition(stop, nuc, count)

        else:
            cp = nuc

        if cp.ambiguous:
            self._ambiguous[stop] = cp

        self._seq.append(cp)
        self._count = count
        self._ref = ref

    def add_nuc(self, ref, pos, nuc, count):

        if self._start == -1:
            return self.initialize(ref, pos, pos, nuc, count)

        assert self.ref == ref

        if not isinstance(nuc, ConsensusPosition):
            cp = ConsensusPosition(pos, nuc, count)

        else:
            cp = nuc

        if not pos == self._stop+1:
            raise RuntimeError('Adding to consensus sequence needs to be contiguous')

        if cp.ambiguous:
            self._ambiguous[pos] = cp

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

    @property
    def coverage(self):
        return float(self.count) / float(self.stop - self.start + 1)

class ConsensusPosition(object):

    def __init__(self, pos, nuc, count):

        self._pos = -1
        self._nuc = None
        self._count = -1
        self._ambiguous = False

        self.initialize(pos, nuc, count)

    def initialize(self, pos, nuc, count):

        self.pos = pos
        self.count = count

        if isinstance(nuc, list):
            self.ambiguous = True
            self.nuc = nuc

        elif isinstance(nuc, basestring):
            self.nuc = nuc.upper()

        else:
            raise ValueError('Nucleotide position must be of type string')

    @property
    def pos(self):
        return self._pos

    @pos.setter
    def pos(self, val):
        if isinstance(val, int):
            self._pos = val

        else:
            raise ValueError('Nucleotide position should be of type int')

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
            raise ValueError('Nucleotide position must be string/list type')

    @property
    def count(self):
        return self._count
    
    @count.setter
    def count(self, val):
        if isinstance(val, int):
            self._count = val

        else:
            raise ValueError('Nucleotide count must be of type int')

    @property
    def ambiguous(self):
        return self._ambiguous
    
    @ambiguous.setter
    def ambiguous(self, val):

        if isinstance(val, bool):
            self._ambiguous = val

        else:
            raise ValueError('Ambiguity must be of type bool')    
    
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

    # Here we are safe because the logic will replace
    # asterisks with gaps
    return match.end()

def insert_action(read_calls, index, final, match):
    
    # Just to cover our bases to make sure
    # that the logic called the proper
    # function
    assert read_calls[index] == '+'

    final[-1] += match.group(2).upper()

    return match.end()

def _process_iter_re(read_calls, match_obj, index=0, final=None):

    if not final:
        final = []

    upto = match_obj.start()

    while index < upto:
        final.append(read_calls[index])
        i+=1

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
            index, final = _process_iter_re(read_calls, match, index, final)
            index = action(read_calls, index, final, match)

        while index < len(read_calls):
            final.append(read_calls[index])
            index += 1

    elif matcher:
        
        for match in matcher.finditer(read_calls):

            index, final = _process_iter_re(read_calls, match, index, final)
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

        if ins:
            extracted = process_iter_re(read_calls, match_objs=ins, action=insert_action)

        else:
            extracted = process_iter_re(read_calls, match_objs=dels, action=delete_action)

        return True, ref, position, read_count, read_calls

    read_calls_list = []

    # Remove the dollar signs and other things
    # that don't tell us things we actually want to know
    read_calls = _remove.sub(r'', read_calls)
    read_calls = _reference.sub(reference_call, read_calls).upper()
    
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
    assert len(read_calls_list) == read_count

    read_calls = ''.join(read_calls_list)

    if _substitutions.match(read_calls):
        # Replace all of the calls that are forward or reverse
        # matches with the reference call

        return False, ref, position, read_count, read_calls


    else:
        # If there is no other call than the reference
        # we can just return the call
        return False, ref, position, read_count, reference_call

def build_consensus(pileup_file):

    reference = None
    consensus = []
    current_consensus = ConsensusSequence()
    
    pile = pileup_iterator(flname=pileup_file)

    while pile:

        next_line = next(pile)
        indel, ref, position, read_count, read_calls = process_line(next_line)

        if position < current_consensus.stop:
            # We've started on the next reference
            # in the pileup file
            
            consensus.append(current_consensus)
            
            yield reference, consensus

            consensus = []
            current_consensus = ConsensusSequence()

        reference = ref

        if not read_count:
            current_consensus.add_nuc(reference, position, '-', read_count)
            continue

        if indel:
            # This is because of insertions
            # or deletions
            raise NotImplementedError('Indel analysis not yet implemented')

        current_consensus.add_nuc(reference, position, read_calls, read_count)


if __name__ == '__main__':

    seqs = []
    pileup_path = os.path.join(os.getcwd(), 'h.pup')

    for ref, seqs in build_consensus(pileup_path):

        print(ref)
        print(seqs)