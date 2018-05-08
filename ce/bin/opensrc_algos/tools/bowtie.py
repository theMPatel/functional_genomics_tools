# -*- coding: utf-8 -*-

import os
import subprocess as sp
from .tools import (
    parse_paired_files
)
from .environment import (
    log_message,
    valid_dir
)

def DoBowtie2Index(fastaFlName, indexFolder, env):
    indexFlName = os.path.join(indexFolder, 'index')

    cmdLine = [
        'python',
        os.path.join(env.toolsdir, 'bowtie2', 'bowtie2-build'),
        fastaFlName,
        indexFlName
    ]

    log_message("Bowtie index cmdline: " + ' '.join(cmdLine), 1)
    sp.call(cmdLine, shell=True, cwd=indexFolder)

    log_message("Bowtie index is done!", 1)
    return indexFlName


def DoBowtie2(indexFl, env, nThreads, singleReadFl=None, pairedReadFl1=None, pairedReadFl2=None, settings=None):
    outputFlName = os.path.join(env.localdir, "output.sam")

    cmdLine = [
        'perl',
        os.path.join(env.toolsdir, 'bowtie2', 'bowtie2'),
        '-x', os.path.normpath(indexFl),
        '-p', str(nThreads),
        '--reorder'
    ]
    if 'Bowtie2' in settings:
        for k, v in settings['Bowtie2'].iteritems():
            cmdLine.extend([k, str(v)] if v else [k])

    if singleReadFl:
        cmdLine += ['-U', singleReadFl.replace('"', '')]
    if pairedReadFl1 and pairedReadFl2:
        cmdLine += ['-1', pairedReadFl1.replace('"', ''), '-2', pairedReadFl2.replace('"', '')]
    cmdLine += ['-S', outputFlName]

    log_message("Bowtie cmdline: " + ' '.join(cmdLine), 1)
    #print(os.getenv('PATH'))
    
    child = sp.Popen(cmdLine, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, env={'PATH' : os.getenv('PATH')})

    out, err = child.communicate()

    if err:
        log_message(err.strip())

    #sp.call(cmdLine, shell=True)

    # the following set of commands "just" use the BAMCoveragePrinter to create the data needed to construct the coverage vector
    if 'BAMCoveragePrinter' in settings:
        cmdLine = [
            os.path.join(env.toolsdir, 'BAMCoveragePrinter', 'BAMCoveragePrinter'),
            os.path.normpath(outputFlName),
        ]
        for k, v in settings['BAMCoveragePrinter'].iteritems():
            cmdLine.extend([k, str(v)] if v else [k])
        log_message("BAMCoveragePrinter: " + ' '.join(cmdLine), 1)
        #sp.call(cmdLine, cwd=env.localdir, shell=True)

        outputFlName_coverages = os.path.join(env.localdir, 'raw_coverages.txt')
        return outputFlName, outputFlName_coverages

    # the following set of commands depend on samtools to create the data needed to construct the coverage vector
    if 'samtools' in settings:
        outputFlName_bam = outputFlName.replace('.sam', '.bam')
        cmdLine = [
            os.path.join(env.toolsdir, 'samtools', 'samtools'),
            'view',
            '-bS', outputFlName,
            '>', outputFlName_bam
        ]
        log_message("samtools cmdline: " + ' '.join(cmdLine), 1)
        sp.call(cmdLine, shell=True)

        outputFlName_bam_sorted = outputFlName.replace('.sam', '_sorted')
        cmdLine = [
            os.path.join(env.toolsdir, 'samtools', 'samtools'),
            'sort',
            outputFlName_bam,
            outputFlName_bam_sorted
        ]
        log_message("samtools cmdline: " + ' '.join(cmdLine), 1)
        sp.call(cmdLine, shell=True)
        outputFlName_bam_sorted = outputFlName_bam_sorted + '.bam'

        outputFlName_mpileup = outputFlName.replace('.sam', '.mpileup')
        cmdLine = [
            os.path.join(env.toolsdir, 'samtools', 'samtools'),
            'mpileup',
            outputFlName_bam_sorted,
            '>', outputFlName_mpileup
        ]
        log_message("samtools cmdline: " + ' '.join(cmdLine), 1)
        sp.call(cmdLine, shell=True)

        return outputFlName, outputFlName_mpileup


class PositionCoverage(object):
    def __init__(self):
        for k in ['nAFwd', 'nCFwd', 'nGFwd', 'nTFwd', 'nGapFwd',
                  'nARev', 'nCRev', 'nGRev', 'nTRev', 'nGapRev']:
            self.__setattr__(k, 0)

def CallsToCov(calls):
    cov = PositionCoverage()
    insertions = []
    deletions = []
    skipCount = 0
    for i, b in enumerate(calls):
        if skipCount:
            skipCount -= 1
            continue
        if b == '^':
            skipCount = 1
            continue
        if b == 'A':
            cov.nAFwd += 1
            continue
        if b == 'C':
            cov.nCFwd += 1
            continue
        if b == 'G':
            cov.nGFwd += 1
            continue
        if b == 'T':
            cov.nTFwd += 1
            continue
        if b == 'a':
            cov.nARev += 1
            continue
        if b == 'c':
            cov.nCRev += 1
            continue
        if b == 'g':
            cov.nGRev += 1
            continue
        if b == 't':
            cov.nTRev += 1
            continue
        if b == '+':  # an insertions
            skipCountStr = ''.join(ele for ele in calls[i + 1:] if ele.isdigit())
            insertSize = int(skipCountStr)
            skipCount = len(skipCountStr) + insertSize
            insertions.append(calls[i + 1:i + 1 + insertSize])
        if b == '-':  # a deletion
            skipCountStr = ''.join(ele for ele in calls[i + 1:] if ele.isdigit())
            insertSize = int(skipCountStr)
            skipCount = len(skipCountStr) + insertSize
            deletions.append(calls[i + 1:i + 1 + insertSize])
            if calls[i + 1].isupper():
                cov.nGapFwd += 1
            else:
                cov.nGapRev += 1

    return cov, insertions, deletions


class CoverageVector(object):
    def __init__(self, l):
        self._fwd_A = [0]*l
        self._fwd_C = [0]*l
        self._fwd_G = [0]*l
        self._fwd_T = [0]*l
        self._fwd_Tot = [0]*l
        self._rev_A = [0]*l
        self._rev_C = [0]*l
        self._rev_G = [0]*l
        self._rev_T = [0]*l
        self._rev_Tot = [0]*l

    def Add(self, pos, cov):
        self._fwd_A[pos] += cov.nAFwd
        self._fwd_C[pos] += cov.nCFwd
        self._fwd_G[pos] += cov.nGFwd
        self._fwd_T[pos] += cov.nTFwd
        self._fwd_Tot[pos] += cov.nAFwd + cov.nCFwd + cov.nGFwd + cov.nTFwd + cov.nGapFwd
        self._rev_A[pos] += cov.nARev
        self._rev_C[pos] += cov.nCRev
        self._rev_G[pos] += cov.nGRev
        self._rev_T[pos] += cov.nTRev
        self._rev_Tot[pos] += cov.nARev + cov.nCRev + cov.nGRev + cov.nTRev + cov.nGapRev

    def GetCovA(self, pos):
        return self._fwd_A[pos] + self._rev_A[pos]

    def GetCovC(self, pos):
        return self._fwd_C[pos] + self._rev_C[pos]

    def GetCovG(self, pos):
        return self._fwd_G[pos] + self._rev_G[pos]

    def GetCovT(self, pos):
        return self._fwd_T[pos] + self._rev_T[pos]

    def GetCovGap(self, pos):
        return self.GetCovTot(pos) - (self.GetCovA(pos) + self.GetCovC(pos) + self.GetCovG(pos) + self.GetCovT(pos))

    def GetCovTot(self, pos):
        return self._fwd_Tot[pos] + self._rev_Tot[pos]

    def GetCovAFwd(self, pos):
        return self._fwd_A[pos]

    def GetCovCFwd(self, pos):
        return self._fwd_C[pos]

    def GetCovGFwd(self, pos):
        return self._fwd_G[pos]

    def GetCovTFwd(self, pos):
        return self._fwd_T[pos]

    def GetCovGapFwd(self, pos):
        return self.GetCovTotFwd(pos) - (
                    self.GetCovAFwd(pos) + self.GetCovCFwd(pos) + self.GetCovGFwd(pos) + self.GetCovTFwd(pos))

    def GetCovTotFwd(self, pos):
        return self._fwd_Tot[pos]

    def GetCovARev(self, pos):
        return self._rev_A[pos]

    def GetCovCRev(self, pos):
        return self._rev_C[pos]

    def GetCovGRev(self, pos):
        return self._rev_G[pos]

    def GetCovTRev(self, pos):
        return self._rev_T[pos]

    def GetCovGapRev(self, pos):
        return self.GetCovTotRev(pos) - (
                    self.GetCovARev(pos) + self.GetCovCRev(pos) + self.GetCovGRev(pos) + self.GetCovTRev(pos))

    def GetCovTotRev(self, pos):
        return self._rev_Tot[pos]

    def GetCov(self, pos, base):
        if base == 'A':
            return self.GetCovA(pos)
        if base == 'C':
            return self.GetCovC(pos)
        if base == 'G':
            return self.GetCovG(pos)
        if base == 'T':
            return self.GetCovT(pos)
        if base == '-':
            return self.GetCovGap(pos)
        return 0

    def GetCovFwd(self, pos, base):
        if base == 'A':
            return self.GetCovAFwd(pos)
        if base == 'C':
            return self.GetCovCFwd(pos)
        if base == 'G':
            return self.GetCovGFwd(pos)
        if base == 'T':
            return self.GetCovTFwd(pos)
        if base == '-':
            return self.GetCovGapFwd(pos)
        return 0

    def GetCovRev(self, pos, base):
        if base == 'A':
            return self.GetCovARev(pos)
        if base == 'C':
            return self.GetCovCRev(pos)
        if base == 'G':
            return self.GetCovGRev(pos)
        if base == 'T':
            return self.GetCovTRev(pos)
        if base == '-':
            return self.GetCovGapRev(pos)
        return 0

    def GetRelCov(self, pos, base):
        if self.GetCovTot(pos) == 0:
            return 0.0
        return float(self.GetCov(pos, base)) / float(self.GetCovTot(pos))

    def IsAmbiguous(self, pos, minCov):
        return sum(self.GetCov(pos, b) > minCov for b in 'ACGT') > 1


class Coverage(object):

    def __init__(self, references):
        self._coverages = {
            k: CoverageVector(len(v.sequence)) for k, v in references.iteritems()
        }

    def LoadRawCoverages(self, flName):
        for line in open(flName):
            line = line.strip()
            if not line:  # an empty line
                continue
            if line[0] == '#':  # a comment line
                continue

            if line[0:3] == '>>>':  # a contig id and a length
                contigId, _ = line[3:].split('\t')
                continue

            parts = map(int, line.split('\t'))
            position = parts[0]
            posCov = PositionCoverage()
            posCov.nAFwd = parts[1]
            posCov.nCFwd = parts[2]
            posCov.nGFwd = parts[3]
            posCov.nTFwd = parts[4]
            posCov.nGapFwd = parts[9] - sum(parts[1:5])
            posCov.nARev = parts[5]
            posCov.nCRev = parts[6]
            posCov.nGRev = parts[7]
            posCov.nTRev = parts[8]
            posCov.nGapRev = parts[10] - sum(parts[5:9])

            self._coverages[contigId].Add(position, posCov)

    def LoadMPilupData(self, mPileUpFlName):

        for line in open(mPileUpFlName):
            line = line.strip()
            if not line:  # an empty line
                continue
            if line[0] == '@':  # a comment line
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            contigId = parts[0]
            position = int(parts[1]) - 1
            refCall = parts[2]
            calls = parts[4].replace('.', refCall.upper()).replace(',', refCall.lower())
            posCovs, insertions, deletions = CallsToCov(calls)
            self._coverages[contigId].Add(position, posCovs)

    def GetCoverage(self, contigId):
        # print contigId
        return self._coverages.get(contigId, CoverageVector(0))

def RunBowtie2(dbInfo, readFiles, mappingSettings, nThreads, env):
    # create the database
    log_message("Exporting search data")
    refsFlName = os.path.join(env.localdir, 'refs.fasta')
    dbInfo.ExportRefsToFasta(refsFlName)

    # create the index
    indexDir = os.path.join(env.localdir, 'index')
    valid_dir(indexDir)
    indexFl = DoBowtie2Index(refsFlName, indexDir, env)

    # establish the pairedness of the files
    singleFiles, pairedFiles = readFiles
    # report the files we are going to use
    log_message("Found the following single-end files:")
    for fl in singleFiles:
        log_message(fl, 1)
    log_message("Found the following paired-end files:")
    for fl in pairedFiles:
        log_message("{}, {}".format(fl[0], fl[1]), 1)

    # perform the mapping
    log_message("Performing mapping...")

    # using the following index file
    # run the contamination detection
    bowtie2Results = Coverage(dbInfo.RefSequences)
    log_message("Going through single-end files...")
    for fl in singleFiles:
        log_message("Handling file {} ...".format(fl), 1)
        samFlName, mPilupFlName = DoBowtie2(indexFl, env=env, nThreads=nThreads, singleReadFl=fl,
                                            settings=mappingSettings)
        log_message("Performing read assignment ...", 1)
        if mappingSettings['LoadCoverage']:
            if 'BAMCoveragePrinter' in mappingSettings and mPilupFlName:
                bowtie2Results.LoadRawCoverages(mPilupFlName)
            elif 'samtools' in mappingSettings:
                bowtie2Results.LoadMPilupData(mPilupFlName)
    log_message("Going through paired-end files...")
    for fl in pairedFiles:
        log_message("Handling files {} and {}...".format(fl[0], fl[1]), 1)
        samFlName, mPilupFlName = DoBowtie2(indexFl, env=env, nThreads=nThreads, pairedReadFl1=fl[0],
                                            pairedReadFl2=fl[1], settings=mappingSettings)
        log_message("Performing read assignment ...", 1)
        if mappingSettings['LoadCoverage']:
            if 'BAMCoveragePrinter' in mappingSettings and mPilupFlName:
                bowtie2Results.LoadRawCoverages(mPilupFlName)
            elif 'samtools' in mappingSettings:
                bowtie2Results.LoadMPilupData(mPilupFlName)
    log_message("Mapping done!")

    return bowtie2Results

