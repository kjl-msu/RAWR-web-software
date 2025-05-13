#!/usr/bin/env python
# -*- coding: utf-8 -*-

import seqs
import sampleSeq
import numpy as np

def calPairSupport(alndata, sampleIndex, sampleSeqData, sampleAlnData):
    sampleSeqIndexData = seqs.sampleSeqDataToSampleSeqIndexData(
        sampleIndex, sampleSeqData)
    samplePairs = seqs.countPairs(sampleSeqIndexData, alndata.shape[0],
                                  alndata.shape[1])
    sampleAlnIndexData = seqs.sampleAlnDataToSampleAlnIndexData(
        sampleIndex, sampleSeqData, sampleAlnData)
    positivePairs = seqs.countPairs(sampleAlnIndexData, alndata.shape[0],
                                    alndata.shape[1])
    pairSupport = np.array(positivePairs) / np.array(samplePairs)
    pairSupport = np.where(np.isnan(pairSupport), 0, pairSupport)
    return pairSupport, samplePairs


def calColumnSupport(alndata, sampleIndex, sampleSeqData, sampleAlnData):
    sampleSeqIndexData = seqs.sampleSeqDataToSampleSeqIndexData(
        sampleIndex, sampleSeqData)
    sampleColumnPairs = seqs.countColumnPairs(sampleSeqIndexData,
                                              alndata.shape[1])
    sampleAlnIndexData = seqs.sampleAlnDataToSampleAlnIndexData(
        sampleIndex, sampleSeqData, sampleAlnData)
    positiveColumnPairs = seqs.countColumnPairs(sampleAlnIndexData,
                                                alndata.shape[1])
    columnSupport = np.array(positiveColumnPairs) / np.array(sampleColumnPairs)
    columnSupport = np.where(np.isnan(columnSupport), 0, columnSupport)
    return columnSupport, sampleColumnPairs


def pairTruePositiveRate(trueAlnData, alndata):
    trueIndexData, estiIndexData = seqs.trueAlnDataToTrueAlnIndexData(
        trueAlnData, alndata)
    truePairs = seqs.countPairs(trueIndexData, alndata.shape[0],
                                alndata.shape[1])
    estiPairs = seqs.countPairs(estiIndexData, alndata.shape[0],
                                alndata.shape[1])
    return truePairs, estiPairs


def colTruePositiveRate(trueAlnData, alndata):
    trueIndexData, estiIndexData = seqs.trueAlnDataToTrueAlnIndexData(
        trueAlnData, alndata)
    trueColumnPairs = seqs.countColumnPairs(trueIndexData, alndata.shape[1])
    estiColumnPairs = seqs.countColumnPairs(estiIndexData, alndata.shape[1])
    truePositiveForColumn = np.array(trueColumnPairs) / np.array(
        estiColumnPairs)
    truePositiveForColumn = np.where(
        np.isnan(truePositiveForColumn), 0, truePositiveForColumn)
    return truePositiveForColumn


def sampleAndGetSampleAlnData(alndata, samplenum, sampleMethod, mu):
    """
    Use Bootstrap with normal distribution to sample sequences from the original alignment
    Align sequences and calculate the support value for sampled sites.
    :param alndata: pandas dataframe, input alignment
            row_index: seq id
            column_index: 0 to length of input alignment
            content: a t c g and - of input alignment at each position
    :param mu: int, mean of the normal distribution.
    :param samplenum: the index of sampled file.
    :return: defaultdict(int), support value for column index correspond to input alignment.
    """
    import subprocess
    # sample seq data
    if sampleMethod == 'bootstrapNormal':
        sampleIndex, sampleSeqData = sampleSeq.bootstrapNormal(
            alndata, mu, 100)
    elif sampleMethod == 'bootstrap':
        sampleIndex, sampleSeqData = sampleSeq.bootstrap(alndata)
    print("sampleIndex:", sampleIndex)
    # write seq data to outfile
    sampleIdxFile = 'sample' + str(samplenum) + '.index'
    sampleSeqFile = 'sample' + str(samplenum) + '.seq.fasta'
    sampleAlnFile = 'sample' + str(samplenum) + '.aln.fasta'
    with open(sampleSeqFile, 'w') as outf:
        for i in sampleSeqData.index:
            outf.write('>' + str(sampleSeqData.loc[i].name) + '\n')
            outf.write(''.join(sampleSeqData.loc[i]).replace('-', '') + '\n')
    np.save(sampleIdxFile, np.array(sampleIndex))
    # generate alignment
    cmd = 'mafft ' + sampleSeqFile + ' > ' + sampleAlnFile
    with open("NUL", "w") as fh:
        p = subprocess.Popen(cmd, shell=True, stdout=fh, stderr=fh)
        output = p.communicate()[0]
    sampleAlnData = seqs.getAlnData(sampleAlnFile)
    return sampleAlnData


# def sampleAndGetColSupport(alndata, samplenum, sampleMethod, mu):
#     """
#     Use Bootstrap with normal distribution to sample sequences from the original alignment
#     Align sequences and calculate the support value for sampled sites.
#     :param alndata: pandas dataframe, input alignment
#             row_index: seq id
#             column_index: 0 to length of input alignment
#             content: a t c g and - of input alignment at each position
#     :param mu: int, mean of the normal distribution.
#     :param samplenum: the index of sampled file.
#     :return: defaultdict(int), support value for column index correspond to input alignment.
#     """
#     import subprocess
#     # sample seq data
#     if sampleMethod == 'bootstrapNormal':
#         sampleIndex, sampleSeqData = sampleSeq.bootstrapNormal(alndata, mu, 100)
#     elif sampleMethod == 'bootstrap':
#         sampleIndex, sampleSeqData = sampleSeq.bootstrap(alndata)
#     print("sampleIndex:", sampleIndex)
#     # write seq data to outfile
#     sampleIdxFile = 'sample' + str(samplenum) + '.index'
#     sampleSeqFile = 'sample' + str(samplenum) + '.seq.fasta'
#     sampleAlnFile = 'sample' + str(samplenum) + '.aln.fasta'
#     with open(sampleSeqFile, 'w') as outf:
#         for i in sampleSeqData.index:
#             outf.write('>' + sampleSeqData.loc[i].name + '\n')
#             outf.write(''.join(sampleSeqData.loc[i]).replace('-', '') + '\n')
#     with open(sampleIdxFile, 'w') as outf:
#         outstr = ','.join([str(x) for x in sampleIndex])
#         outf.write(outstr + '\n')
#     # generate alignment
#     cmd = 'mafft ' + sampleSeqFile + ' > ' + sampleAlnFile
#     with open("NUL", "w") as fh:
#         p = subprocess.Popen(cmd, shell=True, stdout=fh, stderr=fh)
#         output = p.communicate()[0]
#     sampleAlnData = seqs.getAlnData(sampleAlnFile)
#     columnSupport, sampleColumnPairs=calColumnSupport(alndata, sampleIndex, sampleSeqData, sampleAlnData)
#     return columnSupport, sampleColumnPairs


def main():
    import argparse
    import seqs
    import numpy as np
    import os
    from os import listdir
    from os.path import isfile, join

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
    args = parser.parse_args()

    alnfile = args.alnfile
    sampleDir = args.sampleDir
    if not args.sampleNum:
        sampleNum = max([
            int(f.strip().split(".")[0][6:]) for f in listdir("samples")
            if f[:6] == "sample"
        ])
    else:
        sampleNum = int(args.sampleNum)
    alndata = seqs.getAlnData(alnfile)
    seqnum, seqlen = alndata.shape
    pairInOneCol = int(seqnum * (seqnum - 1) / 2)
    totalSamplePairs = np.array([0 for x in range(seqlen * pairInOneCol)])
    totalPositivePairs = np.array([0 for x in range(seqlen * pairInOneCol)])

    for i in range(1, sampleNum + 1):
        print(i)
        sampleAlnFile = sampleDir + "/sample" + str(i) + ".aln.fasta"
        sampleIdxFile = sampleDir + "/sample" + str(i) + ".index.npy"
        if os.path.isfile(sampleAlnFile) and os.path.isfile(sampleIdxFile):
            sampleAlnData = seqs.getAlnData(sampleAlnFile)
            sampleIndex = seqs.getIndexData(sampleIdxFile)
            sampleSeqData = alndata.iloc[:, sampleIndex]
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
        warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
        support = totalPositivePairs / totalSamplePairs             # numpy handles this innately, suppress the warning
    support = np.where(np.isnan(support), 0, support)

if __name__ == "__main__":
    main()
