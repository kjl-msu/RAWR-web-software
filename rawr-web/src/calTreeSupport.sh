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
    treeFile=RAxML_bestTree.$i
	$mafft $seqFile > $alnFile
    $raxml -f a -s $alnFile -n $i -m GTRCAT -p $RANDOM -x $RANDOM  -# 10

    mv RAxML_bestTree.$i $i.tree
done
rm RAxML_*

cd $workDir

cat samples/*.tree > sample.trees
$raxml -f b -m GTRGAMMA -t input.tree -z sample.trees -n support
mv RAxML_bipartitions.support tree.support.txt
rm RAxML_*

echo "python $scriptDir/plotTree.py --treeFile tree.support.txt --figFile tree.support.png"
python $scriptDir/plotTree.py --treeFile tree.support.txt --figFile tree.support.png

