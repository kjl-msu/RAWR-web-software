#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2020 12 - Wei Wang <weiwang.msu@gmail.com>
# This file is used for MSA and tree support estimation.

import os
import warnings
import numpy as np
import seqs
import itertools
import argparse
import MSA_support_csv_2_jalview_sequence_annotation

def calculateMSASupport(alndata, sampleDir, samplenum):
    """
    calculate MSA support.
    """
    seqnum, seqlen = alndata.shape
    pairInOneCol = int(seqnum * (seqnum - 1) / 2)
    totalSamplePairs = np.array([0 for x in range(seqlen * pairInOneCol)])
    totalPositivePairs = np.array([0 for x in range(seqlen * pairInOneCol)])

    for i in range(1, samplenum + 1):
        # print(i)
        sampleAlnFile = sampleDir + "/" + str(i) + ".aln.fasta"
        sampleIdxFile = sampleDir + "/" + str(i) + ".index"
        if os.path.isfile(sampleAlnFile) and os.path.isfile(sampleIdxFile):
            sampleAlnData = seqs.getAlnData(sampleAlnFile)
            sampleIndex = np.loadtxt(sampleIdxFile)
            sampleSeqData = alndata.iloc[:, sampleIndex]
            # print(sampleSeqData)
            sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
            sampleSeqIndexData = seqs.sampleSeqDataToSampleSeqIndexData(
                sampleIndex, sampleSeqData)
            samplePairs = seqs.countPairs(sampleSeqIndexData, alndata.shape[0],
                                          alndata.shape[1])
            sampleAlnIndexData = seqs.sampleAlnDataToSampleAlnIndexData(
                sampleIndex, sampleSeqData, sampleAlnData)
            positivePairs = seqs.countPairs(sampleAlnIndexData,
                                            alndata.shape[0], alndata.shape[1])
            totalSamplePairs = totalSamplePairs + samplePairs
            totalPositivePairs = totalPositivePairs + positivePairs
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in divide")
        support = totalPositivePairs / totalSamplePairs             # numpy handles this innately, suppress the warning
    support = np.where(np.isnan(support), 0, support)
    return support


def getValidPair(alndata):
    seqnum, seqlen = alndata.shape
    pairInOneCol = int(seqnum * (seqnum - 1) / 2)
    validPair = [0 for x in range(seqlen * pairInOneCol)]
    for cidx in range(seqlen):
        col = np.array(alndata.iloc[:, cidx].tolist())
        notGap = np.where(col != "-")[0]
        for i in range(len(notGap)):
            for j in range(i + 1, len(notGap)):
                ridx1 = notGap[i]
                ridx2 = notGap[j]
                pairIdx = int(cidx * pairInOneCol + (2 * seqnum -
                                                     ridx1 - 1) * ridx1 / 2 + (ridx2 - ridx1 - 1))
                validPair[pairIdx] = 1
    return validPair


def writeSupport(alndata, support, supportfile):
    """
    write MSA support value to csv file.
    """
    seqnum, seqlen = alndata.shape
    pairInOneCol = int(seqnum * (seqnum - 1) / 2)
    pairlist = list(itertools.combinations([i for i in range(seqnum)], 2))
    validPair = getValidPair(alndata)
    with open(supportfile, "w") as outf:
        outf.write("columnIndex,rowIndex1,rowIndex2,supportValue\n")
        for pairIdx, s in enumerate(support):
            if validPair[pairIdx] == 1:
                cidx, ridx = divmod(pairIdx, pairInOneCol)
                outf.write(",".join(
                    [str(x) for x in [cidx, pairlist[ridx][0], pairlist[ridx][1], s]]) + "\n")


def main():

    parser = argparse.ArgumentParser(
        description='Parameters for calculate support for MSA.')
    parser.add_argument(
        '--alnfile',
        '-a',
        help='fasta file of estimated alignment.',
        required=True)
    parser.add_argument(
        '--sampleDir',
        '-d',
        help='directory for sampled alignments.',
        required=True)
    parser.add_argument(
        '--sampleNum',
        '-n',
        help='the total number of samples.',
        required=False)
    parser.add_argument(
        '--supportFile',
        '-s',
        help='the support output file.',
        required=False)
    args = parser.parse_args()

    alnfile = args.alnfile
    sampleDir = args.sampleDir
    samplenum = int(args.sampleNum)
    supportfile = args.supportFile

    alndata = seqs.getAlnData(alnfile)
    support = calculateMSASupport(alndata, sampleDir, samplenum)
    writeSupport(alndata, support, supportfile)
    MSA_support_csv_2_jalview_sequence_annotation.support_csv_2_jalview(supportfile,"pink")

    return


if __name__ == '__main__':
    main()
