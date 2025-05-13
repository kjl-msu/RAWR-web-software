#!/usr/bin/env python
import numpy as np
import os
import time
from src import sampleSeq, seqs
import subprocess
import multiprocessing
from joblib import Parallel, delayed

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


def generateSampleSeq(alnData, parameters, num_cores,trigger=None, currValue=0):
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

        if trigger:
            trigger_flag = True
        else:
            trigger_flag = False
        mytuple = []
        for i in range(1, samplenum + 1):
            mytuple.append((alnData, float(reverseRate), sampleDir, i, currValue, trigger_flag))
        results = Parallel(n_jobs=num_cores)(delayed(generateSampleSeq_rawr)(i) for i in mytuple)
        currValue = sum(results)
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

        if trigger:
            trigger_flag = True
        else:
            trigger_flag = False
        mytuple = []
        for i in range(1, samplenum + 1):
            mytuple.append((alnData, int(anchorLen), int(anchorNum), float(reverseRate), barrier, sampleDir, i, currValue, trigger_flag))
        results = Parallel(n_jobs=num_cores)(delayed(generateSampleSeq_seres)(i) for i in mytuple)
        currValue = sum(results)
        trigger.emit(currValue)

    return currValue

def generateSampleSeq_rawr(mytuple):
    alnData, reverseRate, sampleDir, i, currValue, trigger = mytuple
    sampleIndex, sampleSeqData = sampleSeq.rawrSample(alnData, reverseRate)
    writeSampleSeqAndIndex(sampleSeqData, sampleIndex, sampleDir, i)
    if trigger and i % 10 == 1:
        currValue += 1
    return currValue

def generateSampleSeq_seres(mytuple):
    alnData, anchorLen, anchorNum, reverseRate, barrier, sampleDir, i, currValue, trigger = mytuple
    sampleIndex, sampleSeqData = sampleSeq.seresSample(alnData, anchorLen, anchorNum, reverseRate, barrier)
    writeSampleSeqAndIndex(sampleSeqData, sampleIndex, sampleDir, i)
    if trigger and i % 10 == 1:
        currValue += 1
    return currValue
