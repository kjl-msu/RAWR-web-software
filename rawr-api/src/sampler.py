#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2020 12 - Wei Wang <weiwang.msu@gmail.com>
# This file is used for sequence resampling.

from sklearn.preprocessing import normalize
from collections import Counter
import itertools
from . import seqs
import os
import pandas as pd
import numpy as np
from pathlib import Path
import sys

class sampler(object):
    def __init__(self, alnFile, outputDir, reverseRate=0.1, samplenum=10):
        self.alnFile = os.path.normpath(alnFile)
        self.outputDir = os.path.normpath(outputDir)
        self.sampleDir = os.path.join(self.outputDir,"samples")
        if not os.path.exists(self.sampleDir):
            Path(self.sampleDir).mkdir(parents=True, exist_ok=True)
        self.reverseRate = reverseRate
        self.samplenum = samplenum
        self.alndata = seqs.getAlnData(self.alnFile)

    def validSample(self, sampleSeqData):
        """
        Exam the validation of sampled sequences.
        Valid sampled sequences means there is no empty sequence.

        :param sampleSeqData: pandas dataframe, sampled input alignment
                        row_index: seq id;
                        column_index: 0 to length of sampled seq.
        :return:
            True: if resampled sequences contains no empty sequences.
            False: if there is empty sequences in sampled sequences.
        """
        for i in sampleSeqData.index:
            if len(''.join(sampleSeqData.loc[i]).replace('-', '')) == 0:
                print("All empty seq:", i)
                return False
        return True

    def writeSampleSeqAndIndex(self, sampleSeqData, sampleIndex, n):
        """
        Write sampled sequences and index to output file.
        :param sampleSeqData: pandas.DataFrame, sequence data.
        :param sampleIndex: numpy.Array, index data.
        :param outputPath: string, output prefix.
        :param n: int, sample number.
        :return:
        """
        with open(os.path.join(self.sampleDir,str(n) + ".seq.fasta"), 'w') as outf:
            for i in sampleSeqData.index:
                outf.write('>' + sampleSeqData.loc[i].name + '\n')
                outf.write(''.join(sampleSeqData.loc[i]).replace('-', '') +'\n')
        with open(os.path.join(self.sampleDir,str(n) + ".index"), "w") as outf:
            for idx in sampleIndex:
                outf.write(str(idx) + "\n")

    def sample(self):
        sampleIndex = []
        sampleSeqData = pd.DataFrame()
        return sampleIndex, sampleSeqData

    def sampleSeqs(self):
        """
        Generate SERES or RAWR resampled sequences.
        :return: None.
        """
        print("Start generating", self.samplenum, "sampled sequences.")
        for n in range(1, self.samplenum + 1):
            sampleIndex, sampleSeqData = self.sample()
            self.writeSampleSeqAndIndex(sampleSeqData, sampleIndex, n)
        return


class rawrSampler(sampler):
    def __init__(self, alnFile, outputDir, reverseRate=0.1, samplenum=10):
        sampler.__init__(self, alnFile, outputDir, reverseRate, samplenum)

    def sample(self, start=-1):
        isvalid = False
        startIndex, endIndex = self.alndata.columns[0], self.alndata.columns[-1]
        # print(startIndex, endIndex)
        seqlen = self.alndata.shape[1]
        directionName = {1: "r", -1: "l"}
        while isvalid == False:
            if start == -1:
                currIndex = np.random.choice(seqlen, 1)[0] + startIndex
            else:
                currIndex = start
            currStart = currIndex
            direction = np.random.choice([-1, 1])  # -1: left, 1: right
            reverseList = np.random.choice(
                [-1, 1], seqlen + 1,
                p=[self.reverseRate,
                   1 - self.reverseRate])  # -1 turn over, 1 stay the same direction
            sampleIndex = []
            while len(sampleIndex) < seqlen:
                if reverseList[len(sampleIndex) - 1] == -1:
                    # print(currStart, currIndex, directionName[direction])
                    currStart = currIndex
                sampleIndex.append(currIndex)
                currIndex += direction
                if currIndex < startIndex:
                    currIndex, direction = startIndex, 1
                elif currIndex > endIndex:
                    currIndex, direction = endIndex, -1
                else:
                    direction *= reverseList[len(sampleIndex) - 1]
            sampleSeqData = self.alndata.loc[:, sampleIndex]
            isvalid = self.validSample(sampleSeqData)

        sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
        return sampleIndex, sampleSeqData


class seresSampler(sampler):
    def __init__(self, alnFile, outputDir, reverseRate=0.1,
                 samplenum=10, anchorLen=-1, anchorNum=-1):
        sampler.__init__(self, alnFile, outputDir, reverseRate, samplenum)
        self.anchorLen = 5 if anchorLen == -1 else anchorLen
        self.anchorNum = self.alndata.shape[1] // 20 if anchorNum == - \
            1 else anchorNum
        self.barrier = self.getAnchor()

    def distance(self, s1, s2):
        if len(s1) != len(s2):
            raise ValueError("Undefined for sequences of unequal length")
        l1 = np.array(list(s1))
        l2 = np.array(list(s2))

        f1 = np.invert(np.isin(l1, '-'))
        l1 = l1[f1]
        l2 = l2[f1]
        f2 = np.invert(np.isin(l2, '-'))
        l1 = l1[f2]
        l2 = l2[f2]
        return (l1 != l2).mean()

    def getColSimilarity(self, col):
        totalpair = list(itertools.combinations(col, 2))
        count = Counter(totalpair)
        return sum(
            [count[x] for x in count if (x[0] == x[1] and x[0] != '-')])

    def similarity(self):
        sim = []
        for c in self.alndata.columns:
            sim.append(self.getColSimilarity(self.alndata[c]))
        sim = np.array(sim)
        normsim = normalize(sim[:, None], norm='max', axis=0)
        return normsim.flatten()

    def getAnchor(self):
        normsim = self.similarity()
        seqnum, seqlen = self.alndata.shape
        totalScore = {
            i: sum(normsim[(i - self.anchorLen + 1):(i + 1)])
            for i in range(self.anchorLen - 1, seqlen)
        }
        anchorPool = [
            x[0] for x in sorted(
                totalScore.items(), key=lambda d: d[1], reverse=True)
        ]

        minDis = int(max(seqlen / (2 * (self.anchorNum + 1)), self.anchorLen))
        barrier = [0]
        minDis = int(max(seqlen / (2 * (self.anchorNum + 1)), self.anchorLen))
        for i in range(0, 0 + minDis + 1):
            if i in anchorPool:
                anchorPool.remove(i)
        for i in range(seqlen - minDis - 1, seqlen):
            if i in anchorPool:
                anchorPool.remove(i)
        for i in range(self.anchorNum):
            bestAnchorPos = anchorPool[0]
            barrier.append(bestAnchorPos - self.anchorLen + 1)
            barrier.append(bestAnchorPos)
            # delete sites around selected anchor from anchorpool
            for j in range(bestAnchorPos - minDis, bestAnchorPos + minDis + 1):
                if j in anchorPool:
                    anchorPool.remove(j)
        barrier.append(seqlen - 1)
        barrier = sorted(barrier)
        return barrier

    def sample(self):
        seqnum, seqlen = self.alndata.shape
        reverse = [
            np.random.choice(
                np.arange(2), p=[
                    (1 - self.reverseRate), self.reverseRate])
            for x in range(seqlen)
        ]
        # sampleValid = False

        # while sampleValid == False:
        # bi = 0  # 20190521 find bug here, if bi=0, then the sampleValid is
        # always false and it stuck here.
        bi = np.random.choice(np.arange(len(self.barrier)))
        # 20200103 find bug here, if not __init__ial sampleIndex=[] here, the
        # sampleIndex will increase inf__init__ely.
        sampleIndex = []
        prevDirection = np.random.choice(np.arange(2))
        for turnover in reverse:
            if bi == 0:
                direction = 1
            elif bi == (len(self.barrier) - 1):
                direction = 0
            else:
                direction = abs(turnover - prevDirection)
            if bi == (len(self.barrier) - 1) or direction == 0:  # go backward
                temIndex = [x for x in range(
                    self.barrier[bi], self.barrier[bi - 1], -1)]
                # print(
                # bi, direction, barrier[bi], barrier[bi - 1], 'l', sep=', ')
                bi -= 1
            elif bi == 0 or direction == 1:
                temIndex = [x for x in range(
                    self.barrier[bi], self.barrier[bi + 1])]
                # print(
                # bi, direction, barrier[bi], barrier[bi + 1], 'r', sep=', ')
                bi += 1
            sampleIndex.extend(temIndex)
            prevDirection = direction
            if len(sampleIndex) >= seqlen:
                break
        sampleSeqData = self.alndata.iloc[:, sampleIndex]
        # sampleValid = validSample(sampleSeqData)

        sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
        return sampleIndex, sampleSeqData
