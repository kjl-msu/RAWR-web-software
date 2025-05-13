#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is used for MSA and tree support estimation.

import os
import sys
import numpy as np
from src import seqs
from Bio import AlignIO
import subprocess
import itertools
from sys import platform
import warnings
from time import sleep
import ete3

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(
        sys, '_MEIPASS', os.path.dirname(
            os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

def runShellCmd(cmd):
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # we're not using wait() or communicate() so that stdout of slow processes can be shown to users to indicate what the software is processing. wait() comes with the danger of blocking due to memory restrictions and communicate() has the problem where if the command fails to run properly, you wouldn't know and the software continues to the next section without throwing any messages.
    sleep(0.1)
    rc = p.poll()
    while rc is None or ( rc == 0 and not p.stderr):      # only allow programs to exit with valid return code. MAFFT for some reason will return rc = 0 when not finished, so we force it to wait 
        while True:
            line = p.stdout.readline()
            if not line:   # terminates with EOF
                break
            print(line.decode())
            sys.stdout.flush()   # may not be necessary as we are reading from stdout which prevents it from filling up to buffer size and blocking 
        rc = p.poll()   
        sleep(0.1)
    if hasattr(p.stderr, 'read'):
        #print("p.stderr is a BufferedReader")
        print()
    elif rc != 0 and p.stderr:
        print("Program did not exit normally. ERROR MESSAGE: ", p.stderr)
    p.kill()

def estimateSampleAln(parameters, trigger=None, currValue=0):
    print("Estimate alignments. ")
    if platform == "linux" or platform == "linux2":
        # linux
        pathOfBinary = resource_path(os.path.join("mafft-linux","mafft.bat"))
    elif platform == "darwin":
        # OS X
        pathOfBinary = resource_path(os.path.join("mafft-mac","mafft.bat"))
    elif platform == "win32":
        # Windows
        pathOfBinary = resource_path(os.path.join("mafft-win","mafft.bat"))
    else:
        print("I don't believe your system is supported by this software.")
    print(pathOfBinary)
    samplenum = parameters["sampleNum"]
    sampleDir = os.path.join(parameters["outputDir"],"samples")
    for n in range(1, samplenum + 1):
        sampleSeqFile = os.path.join(sampleDir,str(n) + ".seq.fasta")
        sampleAlnFile = os.path.join(sampleDir,str(n) + ".aln.fasta")
        cmd = pathOfBinary + " " + sampleSeqFile + " > " + sampleAlnFile
        runShellCmd(cmd)
        if trigger and n % 10 == 1:
            currValue += 2
            trigger.emit(currValue)
    return currValue


def estimateSampleTree(parameters, trigger=None, currValue=0):
    print("Estimate trees.")
    if platform == "linux" or platform == "linux2":
        # linux
        pathOfBinary = resource_path(os.path.join("raxmlHPC","raxmlHPC_debian"))
        mv = "mv "
        rm = "rm "
        cat = "cat "
        cp = "cp "
    elif platform == "darwin":
        # OS X
        pathOfBinary = resource_path(os.path.join("raxmlHPC","raxmlHPC_darwin"))
        mv = "mv "
        rm = "rm "
        cat = "cat "
        cp = "cp "
    elif platform == "win32":
        # Windows
        pathOfBinary = resource_path(os.path.join("raxmlHPC","raxmlHPC_windows.exe"))
        mv = "move "
        rm = "del "
        cat = "type "
        cp = "copy "
    else:
        print("I don't believe your system is supported by this software.")
    print(pathOfBinary)
    samplenum = parameters["sampleNum"]
    sampleDir = os.path.join(parameters["outputDir"],"samples")
    rng = np.random.uniform(1, 1000000, 2)

    # remove residual RAxML files in this folder
    cmd = rm + os.path.join(sampleDir,"RAxML_info.*")
    runShellCmd(cmd)    

    for n in range(1, samplenum + 1):
        sampleAlnFileFasta = os.path.join(sampleDir,str(n) + ".aln.fasta")
        cmd = pathOfBinary + " -f a -s " + sampleAlnFileFasta + " -n " + str(n) + " -m GTRCAT -p "+str(rng[0])+" -x "+str(rng[1])+ " -# 10 -w " + sampleDir

        runShellCmd(cmd)
        cmd = mv + os.path.join(sampleDir,"RAxML_bestTree." + str(n)) + " " + os.path.join(sampleDir,str(n) + ".tree")
        runShellCmd(cmd)
        if trigger and n % 5 == 1:
            currValue += 3
            trigger.emit(currValue)
    return currValue


def calTreeSupport(inputTree, parameters):
    print("Calculate tree support.")
    outputDir = parameters["outputDir"]
    if platform == "linux" or platform == "linux2":
        # linux
        pathOfBinary = resource_path(os.path.join("raxmlHPC","raxmlHPC_debian"))
        mv = "mv "
        rm = "rm "
        cat = "cat "
        cp = "cp "
    elif platform == "darwin":
        # OS X
        pathOfBinary = resource_path(os.path.join("raxmlHPC","raxmlHPC_darwin"))
        mv = "mv "
        rm = "rm "
        cat = "cat "
        cp = "cp "
    elif platform == "win32":
        # Windows
        pathOfBinary = resource_path(os.path.join("raxmlHPC","raxmlHPC_windows.exe"))
        mv = "move "
        rm = "del "
        cat = "type "
        cp = "copy "
    else:
        print("I don't believe your system is supported by this software.")
    print(pathOfBinary)
    # remove residual RAxML files in this folder
    cmd = rm + os.path.join(outputDir,"RAxML_info.*")
    runShellCmd(cmd)

    cmd = cat + os.path.join(outputDir,"samples","*.tree") + " > " + os.path.join(outputDir,"sample.trees")
    runShellCmd(cmd)
    with open(os.path.join(outputDir,"input.tree"), "w") as outf:
        outf.write(inputTree.write() + "\n")

    cmd = pathOfBinary + " -f b -m GTRGAMMA -t " + os.path.join(outputDir,"input.tree") + " -z " + os.path.join(outputDir,"sample.trees") + " -n support -w " + outputDir
    runShellCmd(cmd)

    cmd = cp + os.path.join(outputDir,"RAxML_bipartitions.support") + " " + os.path.join(outputDir, "tree.support.txt")
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE) # we're not using wait() or communicate() so that stdout of slow processes >
    sleep(0.05)

def calculateMSASupport(alndata, parameters, trigger=None, currValue=0):
    """
    calculate MSA support.
    """
    print("Calculate MSA support.")
    samplenum = parameters["sampleNum"]
    sampleDir = os.path.join(parameters["outputDir"],"samples")
    seqnum, seqlen = alndata.shape
    pairInOneCol = int(seqnum * (seqnum - 1) / 2)
    totalSamplePairs = np.array([0 for x in range(seqlen * pairInOneCol)])
    totalPositivePairs = np.array([0 for x in range(seqlen * pairInOneCol)])

    for i in range(1, samplenum + 1):
        print(i)
        sampleAlnFile = os.path.join(sampleDir,str(i) + ".aln.fasta")
        sampleIdxFile = os.path.join(sampleDir,str(i) + ".index")
        if os.path.isfile(sampleAlnFile) and os.path.isfile(sampleIdxFile):
            sampleAlnData = seqs.getAlnData(sampleAlnFile)
            sampleIndex = np.loadtxt(sampleIdxFile)
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
        if trigger and i % 5 == 1:
            currValue += 3
            trigger.emit(currValue)
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
