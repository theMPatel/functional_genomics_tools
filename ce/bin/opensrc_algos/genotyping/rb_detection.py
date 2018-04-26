###################################################################
#
# Environment specifics for the genotyping algorithm
# You can use this for basically anything you want to
# implement into the Calculation Engine
#
# Author: Milan Patel, with some direction from key mentors ;)
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
from tools.bowtie import RunBowtie2

def CompareCov(cov, seq, start, stop, reflen, minCov):
    return (1.0 * sum(cov.GetCovTot(i) >= minCov for i in xrange(start, stop) if seq[i] in 'ACGT')) / (
                1.0 * (stop - start))

def CompareSeq(cov, seq, start, stop, reflen, minCov):
    return (1.0 * sum(cov.GetCov(i, seq[i]) >= minCov for i in xrange(start, stop) if seq[i] in 'ACGT')) / (
                1.0 * (stop - start))

def Main(dbInfo, readFiles, percentIdentity, minCovDepth, minCovRef, mappingSettings, nThreads, env):
    bowtie2Results = RunBowtie2(dbInfo, readFiles, mappingSettings, nThreads, env)

    # determine coverages
    results = []
    for seqId, seqInfo in dbInfo.Sequences.iteritems():

        refSeqId = dbInfo.GetRefSeqId(seqId)
        alignedRefSeq = dbInfo.GetAlignedRefSeq(refSeqId)
        refSeq = dbInfo.GetRefSeq(refSeqId)

        bowtie2Result = bowtie2Results.GetCoverage(refSeqId)

        assert bowtie2Result is not None, "Hmm are we missing a reference?"

        start = alignedRefSeq.find('-')
        if start == -1: start = 0
        stop = alignedRefSeq.rfind('-')
        if stop == -1: stop = len(alignedRefSeq)
        refLen = len(alignedRefSeq) - alignedRefSeq.count('-')
        covScore = CompareCov(bowtie2Result, alignedRefSeq, start, stop, refLen, minCovDepth)
        seqIdScore = CompareSeq(bowtie2Result, alignedRefSeq, start, stop, refLen, minCovDepth)
        isPresent = covScore >= minCovRef
        isCloseEnough = seqIdScore >= percentIdentity
        if isPresent and isCloseEnough:
            results.append(
                {'seqId': seqIdScore, 'cov': covScore, 'seqInfo': seqInfo}
            )

            # report this out
    return results

# this will find the desired sequences in an exhaustive way
def Main2(dbInfo, readFiles, percentIdentity, minCovDepth, minCovRef, mappingSettings, nThreads, env):
    bowtie2Results = RunBowtie2(dbInfo, readFiles, mappingSettings, nThreads, env)

    # determine coverages
    results = []

    seqDone = {}
    for cntRefs, (seqId, seqInfo) in enumerate(dbInfo.Sequences.iteritems()):

        if cntRefs % 10 == 0:
            print cntRefs

        alignedSeq = dbInfo.GetAlignedSeq(seqId)
        ungappedSequence = dbInfo.GetSeq(seqId)

        # mapping is done against the reference sequences
        # so we need to translate the seqId into the right refSeqId
        # to be able to pick up the coverage vector
        refSeqId = dbInfo.GetRefSeqId(seqId)
        bowtie2Result = bowtie2Results.GetCoverage(refSeqId)

        assert bowtie2Result is not None, "Hmm are we missing a reference?"

        sampleLen = len(alignedSeq)
        covInfo = [
            {b: bowtie2Result.GetCov(i, b) for b in 'ACGT-' if bowtie2Result.GetCov(i, b) >= minCovDepth}
            for i in xrange(sampleLen)
        ]

        # try finding the sequence in the coverage
        # dynamic programming kind of walk-through, not very efficient though
        toCheck = [(0, 0, 0, [])]
        found = []
        minPenalty = 10  # max number of mismatches allowed
        iterCount = 0
        maxDiag = 5
        posDone = set()
        while toCheck:
            # dummy check to avoid infinite loops
            if iterCount > 1.1 * (2 * maxDiag + 1) * sampleLen:
                raise RuntimeError("Problem with species identification - infinite loop")

            # handle a stepstone iteration in the path
            sampleIdx, refIdx, penalty, covlist = toCheck.pop()
            posDescr = '{}|{}'.format(sampleIdx, refIdx)
            if posDescr in posDone: continue

            posDone.add(posDescr)

            iterCount += 1
            if refIdx + 1 == len(ungappedSequence):  # we are at the end of the sequence
                minPenalty = min(minPenalty,
                                 penalty)  # adjust min penalty to force an early stop for non-promising paths
                found.append([penalty, covlist])
            else:
                # check if we have a base
                # if seqInfo._IsReal(ungappedSequence[refIdx], sampleIdx): #slower
                if ungappedSequence[refIdx] in covInfo[sampleIdx]:
                    # if so: consume it
                    if sampleIdx + 1 < sampleLen:
                        toCheck.append((sampleIdx + 1, refIdx + 1, penalty,
                                        covlist + [covInfo[sampleIdx][ungappedSequence[refIdx]]]))
                else:
                    # if not: add a mismatch
                    if sampleIdx + 1 < sampleLen and (penalty + 1 <= minPenalty):
                        toCheck.append((sampleIdx + 1, refIdx + 1, penalty + 1, covlist + [-1]))

                # add a gap if that is possible
                if '-' in covInfo[sampleIdx]:
                    if sampleIdx + 1 < sampleLen:  # and abs(refIdxMap[sampleIdx]-refIdx)<maxDiag:
                        toCheck.append((sampleIdx + 1, refIdx, penalty, covlist + [covInfo[sampleIdx]['-']]))
        if found:
            # determine best-scoring path
            minIdx = min(enumerate(found), key=lambda x: x[1])[0]
            # grab coverage
            sampleCov = found[minIdx][1]
            # calculate abundance
            diffPosns = [i for i in xrange(len(sampleCov)) if bowtie2Result.IsAmbiguous(i, minCovDepth)]
            abundance = sum(sampleCov[i] for i in diffPosns) / len(diffPosns) if len(diffPosns) > 0 else -1
            # add prediction
            results.append(
                {'penalty': found[minIdx][0], 'abundance': abundance, 'seqInfo': seqInfo}
            )
            # keep track of sequence already done
            seqDone[ungappedSequence] = abundance

            # report this out
    return results

