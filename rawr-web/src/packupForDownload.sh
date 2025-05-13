#!/usr/bin/env bash
RANDOM=$$
rootDir=${1%/}
unixTime=$2
supportFile=$3

uploadDir=$rootDir/upload
resultDir=$rootDir/result
staticDir=$rootDir/static

cd $uploadDir/$unixTime
tar czvf samples.tar.gz samples/
rm -r samples/
if [ -f "$uploadDir/$unixTime/tree.support.png" ]
then
    cp $uploadDir/$unixTime/tree.support.png $staticDir/$unixTime.tree.support.png
fi
cp $uploadDir/$unixTime/$supportFile $resultDir/$unixTime.$supportFile
cp $uploadDir/$unixTime/${supportFile}_jalview_annotation.txt $resultDir/$unixTime.${supportFile}_jalview_annotation.txt
cp $uploadDir/$unixTime/samples.tar.gz $resultDir/$unixTime.samples.tar.gz
cp $uploadDir/$unixTime/tree.support.png $resultDir/$unixTime.tree.support.png


