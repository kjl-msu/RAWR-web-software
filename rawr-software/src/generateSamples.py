#!/usr/bin/env python
import numpy as np
import os
import time
from src import sampleSeq, seqs
import subprocess


def writeSampleSeqAndIndex(sampleSeqData, sampleIndex, sampleDir, n):
    """
    Write sampled sequences and index to output file.
    :param sampleSeqData: pandas.DataFrame, sequence data.
    :param sampleIndex: numpy.Array, index data.
    :param outputPath: string, output prefix.
    :param n: int, sample number.
    :return:
    """
    with open(sampleDir + str(n) + ".seq.fasta", 'w') as outf:
        for i in sampleSeqData.index:
            outf.write('>' + sampleSeqData.loc[i].name + '\n')
            outf.write(''.join(sampleSeqData.loc[i]).replace('-', '') +
                       '\n')
    with open(sampleDir + str(n) + ".index", "w") as outf:
        for idx in sampleIndex:
            outf.write(str(idx) + "\n")


def generateSampleSeq(alnData, parameters, trigger=None, currValue=0):
    """
    Generate SERES or RAWR resampled sequences.
    :return: None.
    """
    print("Generate sampled sequences.")
    if parameters["algorithm"] == "RAWR":
        reverseRate = parameters["reverseRate"]
        samplenum = parameters["sampleNum"]
        sampleDir = parameters["outputDir"] + "/samples/"
        if not os.path.exists(sampleDir):
            os.mkdir(sampleDir)

        for n in range(1, samplenum + 1):
            sampleIndex, sampleSeqData = sampleSeq.rawrSample(
                alnData, float(reverseRate))
            writeSampleSeqAndIndex(sampleSeqData, sampleIndex, sampleDir, n)
            if trigger and n % 10 == 1:
                currValue += 1
                trigger.emit(currValue)

    elif parameters["algorithm"] == "SERES":
        anchorLen = parameters["anchorLen"]
        anchorNum = parameters["anchorNum"]
        reverseRate = parameters["reverseRate"]
        samplenum = parameters["sampleNum"]
        sampleDir = parameters["outputDir"] + "/samples/"
        if not os.path.exists(sampleDir):
            os.mkdir(sampleDir)
        barrier = sampleSeq.getAnchor(
            alnData, int(anchorLen), int(anchorNum))
        for n in range(1, samplenum + 1):
            sampleIndex, sampleSeqData = sampleSeq.seresSample(
                alnData, int(anchorLen), int(anchorNum), float(reverseRate), barrier)
            writeSampleSeqAndIndex(sampleSeqData, sampleIndex, sampleDir, n)
            if trigger and n % 10 == 1:
                currValue += 1
                trigger.emit(currValue)

    return currValue


