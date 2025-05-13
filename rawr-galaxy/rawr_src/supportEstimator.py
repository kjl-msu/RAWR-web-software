#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2020 12 - Wei Wang <weiwang.msu@gmail.com>
# This file is used for MSA and tree support estimation.

import os
import sys
import numpy as np
from Bio import AlignIO
import subprocess
import itertools
from time import sleep
import warnings
from . import seqs
from . import MSA_support_csv_2_jalview_sequence_annotation

class supportEstimator(object):

    def __init__(self, alnFile, outputDir, samplenum=10):
        self.outputDir = os.path.normpath(outputDir)
        self.sampleDir = os.path.join(self.outputDir, "samples")
        self.alnFile = os.path.normpath(alnFile)
        self.alndata = seqs.getAlnData(alnFile)
        self.samplenum = int(samplenum)
        self.basePath = getattr(
        sys, '_MEIPASS', os.path.dirname(
            os.path.abspath(__file__)))
        if sys.platform == "linux" or sys.platform == "linux2":
            # linux
            self.raxmlPath = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))), "raxmlHPC","raxmlHPC_debian")
            self.mafftPath = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))), "mafft-linux","mafft.bat")
            self.mv = "mv "
            self.rm = "rm "
            self.cat = "cat "
        elif sys.platform == "darwin":
            # OS X
            self.raxmlPath = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))), "raxmlHPC","raxmlHPC_darwin")
            self.mafftPath = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))), "mafft-mac","mafft.bat")
            self.mv = "mv "
            self.rm = "rm "
            self.cat = "cat "
        elif sys.platform == "win32":
            # Windows
            self.raxmlPath = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))), "raxmlHPC","raxmlHPC_windows.exe")
            self.mafftPath = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))), "mafft-win","mafft.bat")
            self.mv = "move "
            self.rm = "del "
            self.cat = "type "

    def runShellCmd(self, cmd):
        
        ### Alternative subprocess running in commented code below. This is useful for developers because you will see longer subprocesses print to terminal, and it will alert you if the subprocess ends abnormally for any reason. ###
        print(cmd)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # we're not using wait() or communicate() so that stdout of slow processes can be shown to users to indicate the software is still working. wait() comes with the danger of blocking due to memory restrictions and communicate() has the problem where if the command fails to run properly, you wouldn't know and the software continues to the next section without throwing any messages.
        sleep(1)
        rc = p.poll()
        while rc != 0:      # only allow programs to exit with valid return code or EOF 
            while True:
                line = p.stdout.readline()
                if not line:   # terminates with EOF
                    break
                print(line.decode())
                sys.stdout.flush()   # may not be necessary as we are reading from stdout which prevents it from filling up to buffer size and blocking 
            rc = p.poll()   
        assert rc == 0
        """
        ### This subprocess command is nice and silent, but if anything goes awry the program will attempt to continue instead of throwing a subprocess error message ###
        with open("NUL", "w") as fh:
            p = subprocess.Popen(cmd, shell=True, stdout=fh, stderr=fh)
            output = p.communicate()[0]
        return output
        """

    def estimateSampleAln(self):
        print("Estimate alignments. ")
        # pathOfBinary = os.path.join(self.basePath, "mafft-mac/mafft.bat")
        for n in range(1, self.samplenum + 1):
            sampleSeqFile = os.path.join(self.sampleDir,str(n) + ".seq.fasta")
            sampleAlnFile = os.path.join(self.sampleDir,str(n) + ".aln.fasta")
            cmd = self.mafftPath + " " + sampleSeqFile + " > " + sampleAlnFile
            #print(cmd)
            self.runShellCmd(cmd)

class msaSupportEstimator(supportEstimator):
    def __init__(self, alnFile, outputDir, samplenum=10, your_mafftPath=None):
        supportEstimator.__init__(self, alnFile, outputDir, samplenum)
        if your_mafftPath is not None:
            self.mafftPath = your_mafftPath

    def getValidPair(self):
        seqnum, seqlen = self.alndata.shape
        pairInOneCol = int(seqnum * (seqnum - 1) / 2)
        validPair = [0 for x in range(seqlen * pairInOneCol)]
        for cidx in range(seqlen):
            col = np.array(self.alndata.iloc[:, cidx].tolist())
            notGap = np.where(col != "-")[0]
            for i in range(len(notGap)):
                for j in range(i + 1, len(notGap)):
                    ridx1 = notGap[i]
                    ridx2 = notGap[j]
                    pairIdx = int(cidx * pairInOneCol + (2 * seqnum -
                                                         ridx1 - 1) * ridx1 / 2 + (ridx2 - ridx1 - 1))
                    validPair[pairIdx] = 1
        return validPair


    def calculateSupport(self):
        """
        calculate MSA support.
        """
        print("Calculate MSA support.")
        self.estimateSampleAln()
        seqnum, seqlen = self.alndata.shape
        pairInOneCol = int(seqnum * (seqnum - 1) / 2)
        totalSamplePairs = np.array([0 for x in range(seqlen * pairInOneCol)])
        totalPositivePairs = np.array([0 for x in range(seqlen * pairInOneCol)])

        for i in range(1, self.samplenum + 1):
            print(i)
            sampleAlnFile = os.path.join(self.sampleDir, str(i) + ".aln.fasta")
            sampleIdxFile = os.path.join(self.sampleDir, str(i) + ".index")
            if os.path.isfile(sampleAlnFile) and os.path.isfile(sampleIdxFile):
                sampleAlnData = seqs.getAlnData(sampleAlnFile)
                sampleIndex = np.loadtxt(sampleIdxFile)
                sampleSeqData = self.alndata.iloc[:, sampleIndex]
                sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
                sampleSeqIndexData = seqs.sampleSeqDataToSampleSeqIndexData(
                    sampleIndex, sampleSeqData)
                samplePairs = seqs.countPairs(sampleSeqIndexData, self.alndata.shape[0],
                                              self.alndata.shape[1])
                sampleAlnIndexData = seqs.sampleAlnDataToSampleAlnIndexData(
                    sampleIndex, sampleSeqData, sampleAlnData)
                positivePairs = seqs.countPairs(sampleAlnIndexData,
                                                self.alndata.shape[0], self.alndata.shape[1])
                totalSamplePairs = totalSamplePairs + samplePairs
                totalPositivePairs = totalPositivePairs + positivePairs
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
            support = totalPositivePairs / totalSamplePairs             # numpy handles this innately, suppress the warning
        support = np.where(np.isnan(support), 0, support)
        self.writeMSASupport(support)
        support_outfile = os.path.join(self.outputDir, "MSA.support.csv")
        print("MSA support is written to", support_outfile)

        #JalView
        MSA_support_csv_2_jalview_sequence_annotation.support_csv_2_jalview(support_outfile,"pink")
        print("JalView format file for MSA support is written to ", os.path.join(self.outputDir,"MSA.support.csv_jalview_annotation.txt"))
        return support



    def writeMSASupport(self, support):
        """
        write MSA support value to csv file.
        """
        supportfile = os.path.join(self.outputDir,"MSA.support.csv")
        seqnum, seqlen = self.alndata.shape
        pairInOneCol = int(seqnum * (seqnum - 1) / 2)
        pairlist = list(itertools.combinations([i for i in range(seqnum)], 2))
        validPair = self.getValidPair()
        with open(supportfile, "w") as outf:
            outf.write("columnIndex,rowIndex1,rowIndex2,supportValue\n")
            for pairIdx, s in enumerate(support):
                if validPair[pairIdx] == 1:
                    cidx, ridx = divmod(pairIdx, pairInOneCol)
                    outf.write(",".join(
                        [str(x) for x in [cidx, pairlist[ridx][0], pairlist[ridx][1], s]]) + "\n")



class treeSupportEstimator(supportEstimator):
    def __init__(self, alnFile, outputDir, inputTreeFile, samplenum=10, your_mafftPath=None, your_raxmlPath=None):
        supportEstimator.__init__(self, alnFile, outputDir, samplenum)
        self.inputTreeFile = os.path.normpath(inputTreeFile)
        if your_mafftPath is not None:
            self.mafftPath = your_mafftPath
        if your_raxmlPath is not None:
            self.raxmlPath = your_raxmlPath

    def estimateSampleTree(self):
        print("Estimate trees.")
        sampleDir = os.path.join(self.outputDir,"samples")
        rng = np.random.uniform(1, 1000000, 2)
        for n in range(1, self.samplenum + 1):
            sampleAlnFile = os.path.join(sampleDir,str(n) + ".aln.fasta")
            cmd = self.raxmlPath + " -f a -s " + sampleAlnFile + " -n " + \
                str(n) + " -m GTRCAT -p "+str(rng[0])+" -x "+str(rng[1])+ " -# 10 -w " + sampleDir

            self.runShellCmd(cmd)
            cmd = self.mv + os.path.join(sampleDir,"RAxML_bestTree."+ str(n)) + " " + os.path.join(sampleDir,str(n) + ".tree")
            self.runShellCmd(cmd)
        cmd = self.rm + os.path.join(sampleDir,"RAxML_*")
        self.runShellCmd(cmd)


    def calculateSupport(self):
        print("Calculate tree support.")
        self.estimateSampleAln()
        self.estimateSampleTree()
        cmd = self.cat + os.path.join(self.outputDir,"samples","*.tree") + " > " + os.path.join(self.outputDir,"sample.trees")
        self.runShellCmd(cmd)
        cmd = self.raxmlPath + " -f b -m GTRGAMMA -t " + self.inputTreeFile + " -z " + os.path.join(self.outputDir,"sample.trees") + " -n support -w " + self.outputDir
        self.runShellCmd(cmd)
        cmd = self.mv + os.path.join(self.outputDir,"RAxML_bipartitions.support") +" " + os.path.join(self.outputDir,"tree.support.txt")
        self.runShellCmd(cmd)
        cmd = self.rm + os.path.join(self.outputDir,"RAxML_*")
        self.runShellCmd(cmd)
        print("Tree support is written to", os.path.join(self.outputDir,"tree.support.txt"))


