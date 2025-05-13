#!/usr/bin/env bash
RANDOM=$$
workDir=${1%/}
sampleNum=$2
scriptDir=${3%/}
# Detect the platform 
OS="`uname`"
case $OS in
  'Linux')
	echo "Linux"
    raxml=${scriptDir}/raxmlHPC/raxmlHPC_debian
	mafft=${scriptDir}/mafft-linux/mafft.bat
    ;;
  'Darwin') 
	echo "Mac"
    raxml=${scriptDir}/raxmlHPC/raxmlHPC_darwin
	mafft=${scriptDir}/mafft-mac/mafft.bat
    ;;

  *) 
    echo "Your operating system is not compatible with this program."
	;;
esac

cd $workDir/samples
for ((i=1;i<=$sampleNum;i++))
do
    seqFile=$i.seq.fasta
    alnFile=$i.aln.fasta
    $mafft $seqFile > $alnFile
done

echo "python $scriptDir/calMSASupport.py --alnfile $workDir/alignment.fasta --sampleDir $workDir/samples/ --sampleNum $sampleNum --supportFile $workDir/MSA_Support.csv"
python $scriptDir/calMSASupport.py --alnfile $workDir/alignment.fasta --sampleDir $workDir/samples/ --sampleNum $sampleNum --supportFile $workDir/MSA_Support.csv
