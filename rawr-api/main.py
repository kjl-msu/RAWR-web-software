# coding:utf-8
# Available under GPLv3 license
import os
from src import sampler, supportEstimator
from src.MSA_support_csv_2_jalview_sequence_annotation import support_csv_2_jalview

if __name__ == '__main__':
        """
        This file explores RAWR v.1 functions and parameters.
        """
        workDir = os.path.abspath(os.path.dirname(__file__))

        print("### example of using RAWR for MSA, simple default options ###")
        ### example of using RAWR for MSA, simple default options ###
        # 1) set up input and output directories
        alnFile = os.path.join(workDir, "example","dataset-10taxa","alignment.fasta")
        outputDir1 = os.path.join(workDir, "outputDir_rawr_msa")
        # 2) RAWR resample the alignment file into /samples/ directory
        rawr1 = sampler.rawrSampler(alnFile, outputDir1) 
        rawr1.sampleSeqs()
        # 3) use the resampled files to do multiple sequence alignment support estimation
        rawr1_msa = supportEstimator.msaSupportEstimator(alnFile, outputDir1) 
        rawr1_msa.calculateSupport()

        print()
        print("### example of converting MSA output to JalView annotation ###")
        #### example of converting RAWR output to JalView annotation ####
        rawrMSAoutput=os.path.join(outputDir1,"MSA.support.csv")
        support_csv_2_jalview(rawrMSAoutput,"pink")
        print("Jalview annotation file at:" + str(os.path.join(outputDir1,"MSA.support.csv_jalview_annotation.txt")))

        print()
        print("### example of using RAWR for phylogenetic tree support with more options ###")
        ### example of using RAWR for phylogenetic tree support with more options ###
        """
        rawrSampler(alnFile, outputDir, reverseRate=0.1, samplenum=10)
        can input reverseRate (0,1) not including 0 or 1.
        can input any samplenum >= 2
        """
        # 1) set up input and output directories
        alnFile = os.path.join(workDir, "example","dataset-10taxa","alignment.fasta")
        treeFile = os.path.join(workDir, "example","dataset-10taxa","infer.tree")
        outputDir2 = os.path.join(workDir, "outputDir_rawr_tree")
        reverseRate = 0.1
        samplenum = 3
        # 2) RAWR resample the alignment file into /samples/ directory
        rawr2 = sampler.rawrSampler(alnFile, outputDir2,reverseRate,samplenum) 
        rawr2.sampleSeqs()
        # 3) use the resampled files to do phylogenetic tree support estimation
        rawr2_tree = supportEstimator.treeSupportEstimator(alnFile, outputDir2, treeFile,samplenum) 
        rawr2_tree.calculateSupport()

        print()
        print("### example of using SERES for MSA with more options ###")
        ### example of using SERES for MSA with more options ###
        """
        seresSampler(alnFile, outputDir, reverseRate=0.1, samplenum=10, anchorLen=5, anchorNum=20)
        can input reverseRate (0,1) not including 0 or 1.
        can input any samplenum >= 2
        can input any anchorLen >= 1
        can input any anchorNum >= 1
        """
        # 1) set up input and output directories
        alnFile = os.path.join(workDir,"example","dataset-10taxa","alignment.fasta")
        outputDir3 = os.path.join(workDir, "outputDir_seres_msa")
        reverseRate = 0.1
        samplenum = 3
        anchorLen = 5
        anchorNum = 20
        # 2) RAWR resample the alignment file into /samples/ directory
        seres1 = sampler.seresSampler(alnFile, outputDir3,reverseRate,samplenum,anchorLen,anchorNum)  # note: if you want to use SERES algorithm, simply change from rawrSampler to seresSampler here.
        seres1.sampleSeqs()
        # 3) use the resampled files to do multiple sequence alignment support estimation
        seres1_msa = supportEstimator.msaSupportEstimator(alnFile, outputDir3,samplenum) 
        seres1_msa.calculateSupport()

        print()
        print("### example of converting MSA output to JalView annotation ###")
        #### example of converting RAWR output to JalView annotation ####
        seresMSAoutput=os.path.join(outputDir3,"MSA.support.csv")
        support_csv_2_jalview(seresMSAoutput,"pink")
        print("Jalview annotation file at:" + str(os.path.join(outputDir3,"MSA.support.csv_jalview_annotation.txt")))

        print()
        print("### example of using SERES for phylogenetic tree support, simple default options ###")
        ### example of using SERES for phylogenetic tree support, simple default options ###
        # 1) set up input and output directories
        alnFile = os.path.join(workDir, "example","dataset-10taxa","alignment.fasta")
        treeFile = os.path.join(workDir, "example","dataset-10taxa","infer.tree")
        outputDir4 = os.path.join(workDir, "outputDir_seres_tree")
        # 2) RAWR resample the alignment file into /samples/ directory
        seres2 = sampler.seresSampler(alnFile, outputDir4)
        seres2.sampleSeqs()
        # 3) use the resampled files to do phylogenetic tree support estimation
        seres2_tree = supportEstimator.treeSupportEstimator(alnFile, outputDir4, treeFile)
        seres2_tree.calculateSupport()

        print()
        print("### If you want to use your own RAxML and MAFFT instead of the standalone MAFFT version 7.487 and RAxML version 8.2.12 we packaged into /src/, add the path when you call supportEstimator ###")
        ### If you want to use your own RAxML and MAFFT instead of the standalone MAFFT version 7.487 and RAxML version 8.2.12 we packaged into /src/, add the path when you call supportEstimator ###
        """
        # 3) use the resampled files to do multiple sequence alignment support estimation
        rawr1_msa = supportEstimator.msaSupportEstimator(alnFile, outputDir1,samplenum, "/path/to/mafft") 
        supportEstimator.msaSupportEstimator.calculateSupport(rawr1_msa) # alternative way of calling support functions

        # 3) use the resampled files to do phylogenetic tree support estimation
        rawr2_tree = supportEstimator.treeSupportEstimator(alnFile, outputDir2,treeFile,samplenum,"/path/to/mafft","/path/to/raxml")
        supportEstimator.treeSupportEstimator.calculateSupport(rawr2_tree) # alternative way of calling support functions
        """
