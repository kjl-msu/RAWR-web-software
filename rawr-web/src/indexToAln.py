#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import seqs


def indexToSeq(alndata, sampleIndex, indexfile):
    sampleSeqData = alndata.iloc[:, sampleIndex]
    sampleSeqData.columns = [x for x in range(sampleSeqData.shape[1])]
    sampleSeqFile = indexfile.strip().split("index")[:-1][0] + "seq.fasta"
    with open(sampleSeqFile, 'w') as outf:
        for i in sampleSeqData.index:
            outf.write('>' + sampleSeqData.loc[i].name + '\n')
            outf.write(''.join(sampleSeqData.loc[i]).replace('-', '') + '\n')
    return sampleSeqFile


def seqToAln(sampleSeqFile):
    import subprocess
    sampleAlnFile = sampleSeqFile.strip().split("seq")[:-1][0] + "aln.fasta"
    cmd = 'mafft ' + sampleSeqFile + ' > ' + sampleAlnFile
    with open("NUL", "w") as fh:
        p = subprocess.Popen(cmd, shell=True, stdout=fh, stderr=fh)
        output = p.communicate()[0]
    return sampleAlnFile


def main():
    parser = argparse.ArgumentParser(
        description='Parameters for calculate the running AUC.')
    parser.add_argument(
        '--indexfile',
        '-i',
        help='index file of sampled sequences.',
        required=True)
    parser.add_argument(
        '--alnfile',
        '-a',
        help='fasta file of estimated alignment.',
        required=True)
    args = parser.parse_args()
    indexfile = args.indexfile
    alnfile = args.alnfile
    alndata = seqs.getAlnData(alnfile)
    sampleIndex = seqs.getIndexData(indexfile)
    sampleSeqFile = indexToSeq(alndata, sampleIndex, indexfile)
    sampleAlnFile = seqToAln(sampleSeqFile)
    print(sampleSeqFile, sampleAlnFile)


if __name__ == "__main__":
    main()
