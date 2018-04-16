###################################################################
#
# Detects the various stx gene variants for Escherichia
# You can use this for basically anything you want to
# implement into the Calculation Engine
#
# Author: Hannes Pouseele
# Contact: hannes_pouseele@applied-math.com
# Version 1.0
#
###################################################################

import os
import json

from tools.environment import (
    log_algo_version,
    log_message,
    write_results,
    time_now,
    valid_dir,
)

from tools.tools import (
    parse_fasta,
    is_fasta
)

from tools.dbinfo import (
    DbInfo,
    LocusInfo,
    notes_parser
)

from .rb_detection import Main2 as ReadFindPresenceAbsence
from collections import namedtuple

def ValidateInput(settings):
    pass

SeqInfo = namedtuple('SeqInfo', ['locus', 'allele', 'accession', 'sequence', 'other'])
StxSeqInfo = namedtuple('StxSeqInfo', SeqInfo._fields + ('stx', 'variant'))

def ParserForSeq(header, seq, separator=':'):
    parts = header.split('|')[0].split(separator)
    while len(parts) < 2:
        parts.append('')

    stx = parts[0][:4]
    variant = parts[0][-1]
    return StxSeqInfo(parts[0], '', parts[1], seq, '{}:{}'.format(stx, variant), stx, variant)

class StxDbInfo(object):
    def __init__(self, folder=None, seqParser=ParserForSeq, notesParser=notes_parser, acceptFunc=None):
        self._references = {}
        self._notes = {}
        self._sequences = {}

        if folder is not None:
            self.Load(folder, seqParser, notesParser
                      , acceptFunc)

    def Load(self, folder, seqParser, notesParser, acceptFunc):

        def MyAcceptFunc(seqInfo):
            return True

        flName = os.path.join(folder, "stx1-all.fasta")
        self.LoadFile(flName, seqParser, notesParser, MyAcceptFunc)
        flName = os.path.join(folder, "stx2-no-f-all.fasta")
        self.LoadFile(flName, seqParser, notesParser, MyAcceptFunc)
        flName = os.path.join(folder, "stx2-f-all.fasta")
        self.LoadFile(flName, seqParser, notesParser, MyAcceptFunc)

        seqs = parse_fasta(os.path.join(folder, 'stx_consensus.fasta'))
        self._references = {
            k: SeqInfo(k, '', '', v, '') for k, v in seqs.iteritems()
        }

    def LoadFile(self, flName, seqParser, notesParser, acceptFunc):
        seqs = parse_fasta(flName)
        for header, seq in seqs.iteritems():
            info = seqParser(header, seq)
            if acceptFunc and not acceptFunc(info):
                continue

            alleleId = '{}_{}'.format(info.locus, info.allele) if info.allele else info.locus
            idx = 2
            while alleleId in self._sequences:
                alleleId = '{}_{}-{}'.format(info.locus, info.allele, idx)
                idx += 1
            if alleleId not in self._sequences:
                self._sequences[alleleId] = info

    def GetRefSeqId(self, seqId):
        if 'stx2' in seqId:
            if 'f' not in seqId:
                return 'consensus-stx2-no-f'
            else:
                return 'consensus-stx2-f'
        else:
            return 'consensus-stx1'

    def GetSeq(self, seqId):
        return self._sequences[seqId].sequence.replace('-', '')

    def GetAlignedSeq(self, seqId):
        return self._sequences[seqId].sequence

    def GetAlignedRefSeq(self, seqId):
        return self._sequences[seqId].sequence

    def GetRefSeq(self, seqId):
        return self._sequences[seqId].sequence

    def Allele(self, seqId):
        return self._sequences[seqId]

    def Locus(self, locus):
        return self._notes.get(locus, None)

    def ExportToFasta(self, flName):
        valid_dir(os.path.split(flName)[0])
        with open(flName, 'w') as oFile:
            for seqId, seqInfo in self._sequences.iteritems():
                oFile.write('>{}|{}\n{}\n'.format(seqId, len(seqInfo.sequence), seqInfo.sequence))

    def ExportRefsToFasta(self, flName):
        valid_dir(os.path.split(flName)[0])
        with open(flName, 'w') as oFile:
            for seqId, seqInfo in self._sequences.iteritems():
                oFile.write('>{}\n{}\n'.format(seqId, seqInfo.sequence))

    @property
    def RefSequences(self):
        return self._references

    @property
    def Sequences(self):
        return self._sequences

    def ExportRefsToFasta(self, flName):
        valid_dir(os.path.split(flName)[0])
        with open(flName, 'w') as oFile:
            for seqId, seqInfo in self._references.iteritems():
                oFile.write('>{}\n{}\n'.format(seqId, seqInfo.sequence))

def main(settings, env):
    # write basic logging info

    log_message('Starting running reads based stx subtyping algorithm', 1)

    log_algo_version(
        algo_version=None,
        settings=settings,
        env=env
    )

    # perform input data validation
    ValidateInput(settings)

    # load database info
    database_path = env.get_sharedpath(settings.database)

    log_message('Database path found at: {}'.format(
        database_path), 2
    )

    dbInfo = StxDbInfo(database_path, ParserForSeq)

    # Just in case that the reads do not get pushed from the client
    results = {}

    # run the finder
    log_message("Determining fixed-position base calls...")
    if 'query_reads' in settings and settings.query_reads:

        mappingSettings = {
            'Bowtie2': {
                '--reorder': None,  # to avoid bowtie2 from crashing when wrongly doing multithreaded output
                '--local': None,  # perform local mapping
                '--sensitive-local': None,  # default setting: sensitive
                # '--rdg': '3,1', #make it gappy-er
                # '--score-min': 'L,90,0', #adjust the minimum scoring function to f(readlen) = A + B*readlen
                '--no-unal': None,  # avoid writing unaligned reads in the sam file
                '--all': None,  # report all alignments
            },
            'BAMCoveragePrinter': {
                'output_raw_coverages': None,
                'count_multiple_alignments': None,
            },
            # 'samtools' : {}
            'LoadCoverage': True,
        }

        results = ReadFindPresenceAbsence(
            dbInfo,
            settings.query_reads,
            settings.reads.percent_identity,
            settings.reads.min_coverage,
            settings.reads.min_relative_coverage,
            mappingSettings,
            env.threads,
            env
        )


    fixed_results = {
        'results': {},
        'extra' : []
    }

    found_genotypes = set(g['seqInfo'][0] for g in results)
    all_genotypes = {g: False for g in dbInfo.Sequences.keys()}

    for g in found_genotypes:
        all_genotypes[g] = True

    fixed_results['results'] = all_genotypes
    fixed_results['extra'] = results

    write_results('ecoli.stxfinder_genotypes.json', json.dumps(fixed_results))

    # # parse position information into stx types
    # log_message("Performing stx calling on fixed-position base calls...")
    # stxTypeResults = {}
    # if 'stx-aligned_1' in results:
    #     stxTypeResults = {
    #         tpe: sum(any(c in results['stx-aligned_1'][p] for c in call[p]) for p in call)
    #         for tpe, call in dbInfo.Calls.iteritems()
    #     }

    # write_results('StxTypeFinder_stxtypes.json', json.dumps(stxTypeResults))


if __name__ == '__main__':
    from collections import namedtuple
    Settings = namedtuple('Settings', [
        'localDir',
        'resultsDir',
        'sharedDir',
        'toolsDir',
        'query',
        'percentIdentity',
        'minRelCoverage',
        'minOverlapForMerge',
        'searchForFragments'
    ])

    settings = Settings(
        r'c:\users\hannes\temp\test_genotyping\local',
        r'c:\users\hannes\temp\test_genotyping',
        r'C:\tfs\Internal\Test\Hannes\scripts\genotyping\am_databases',
        r'C:\Program Files\Applied Maths\BioNumerics 7.6 bleeding edge\SequenceData\Blast',
        # r'c:\users\hannes\temp\test_genotyping\Escherichia_coli_O157H7_Sakai.fasta.gz',
        # r'c:\users\hannes\temp\test_genotyping\Shigella_boydii_965_58.fasta.gz',
        r'c:\users\hannes\temp\test_genotyping\query.fasta',
        0.9,
        0.6,
        0.9,
        True
    )

    Main(settings)





