#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import itertools
from collections import Counter
from sklearn.preprocessing import normalize


def validSample(sampleSeqData):
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


def rawrSample(alndata, reverseRate, start=-1):
    isvalid = False
    startIndex, endIndex = alndata.columns[0], alndata.columns[-1]
    # print(startIndex, endIndex)
    seqlen = alndata.shape[1]
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
            p=[reverseRate,
                1 - reverseRate])  # -1 turn over, 1 stay the same direction
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
        sampleSeqData = alndata.loc[:, sampleIndex]
        isvalid = validSample(sampleSeqData)

    sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
    return sampleIndex, sampleSeqData


def distance(s1, s2):
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


def getColSimilarity(col):
    totalpair = list(itertools.combinations(col, 2))
    count = Counter(totalpair)
    return sum(
        [count[x] for x in count if (x[0] == x[1] and x[0] != '-')])


def similarity(alndata):
    sim = []
    for c in alndata.columns:
        sim.append(getColSimilarity(alndata[c]))
    sim = np.array(sim)
    normsim = normalize(sim[:, None], norm='max', axis=0)
    return normsim.flatten()


def getAnchor(alndata, anchorlen = "default", anchornum = "default"):
    normsim = similarity(alndata)
    seqnum, seqlen = alndata.shape
    if anchorlen == "default":
        anchorlen = 5
    if anchornum == "default":
        anchornum = seqlen // 20
    totalScore = {
        i: sum(normsim[(i - anchorlen + 1):(i + 1)])
        for i in range(anchorlen - 1, seqlen)
    }
    anchorPool = [
        x[0] for x in sorted(
            totalScore.items(), key=lambda d: d[1], reverse=True)
    ]

    minDis = int(max(seqlen / (2 * (anchornum + 1)), anchorlen))
    barrier = [0]
    minDis = int(max(seqlen / (2 * (anchornum + 1)), anchorlen))
    for i in range(0, 0 + minDis + 1):
        if i in anchorPool:
            anchorPool.remove(i)
    for i in range(seqlen - minDis - 1, seqlen):
        if i in anchorPool:
            anchorPool.remove(i)
    for i in range(anchornum):
        bestAnchorPos = anchorPool[0]
        barrier.append(bestAnchorPos - anchorlen + 1)
        barrier.append(bestAnchorPos)
        # delete sites around selected anchor from anchorpool
        for j in range(bestAnchorPos - minDis, bestAnchorPos + minDis + 1):
            if j in anchorPool:
                anchorPool.remove(j)
    barrier.append(seqlen - 1)
    barrier = sorted(barrier)
    return barrier


def seresSample(alndata, anchorlen, anchornum, reverserate, barrier):
    seqnum, seqlen = alndata.shape
    if anchornum == -1:
        print('No specific anchor number.')
        anchornum = int(seqlen / 20) - 1
        print('Anchor number: ' + str(anchornum + 2))
    else:
        print('User set anchor number: ' + str(anchornum))
        anchornum = int(anchornum) - 2

    # barrier = getAnchor()
    reverse = [
        np.random.choice(np.arange(2), p=[(1 - reverserate), reverserate])
        for x in range(seqlen)
    ]
    # sampleValid = False

    # while sampleValid == False:
    # bi = 0  # 20190521 find bug here, if bi=0, then the sampleValid is always false and it stuck here.
    bi = np.random.choice(np.arange(len(barrier)))
    # 20200103 find bug here, if not initial sampleIndex=[] here, the sampleIndex will increase infinitely.
    sampleIndex = []
    prevDirection = np.random.choice(np.arange(2))
    for turnover in reverse:
        if bi == 0:
            direction = 1
        elif bi == (len(barrier) - 1):
            direction = 0
        else:
            direction = abs(turnover - prevDirection)
        if bi == (len(barrier) - 1) or direction == 0:  # go backward
            temIndex = [x for x in range(barrier[bi], barrier[bi - 1], -1)]
            # print(
            #     bi, direction, barrier[bi], barrier[bi - 1], 'l', sep=', ')
            bi -= 1
        elif bi == 0 or direction == 1:
            temIndex = [x for x in range(barrier[bi], barrier[bi + 1])]
            # print(
            #     bi, direction, barrier[bi], barrier[bi + 1], 'r', sep=', ')
            bi += 1
        sampleIndex.extend(temIndex)
        prevDirection = direction
        if len(sampleIndex) >= seqlen:
            break
    sampleSeqData = alndata.iloc[:, sampleIndex]
    # sampleValid = validSample(sampleSeqData)

    sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
    return sampleIndex, sampleSeqData
