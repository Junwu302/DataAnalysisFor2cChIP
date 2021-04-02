#!/bin/bash
for bam in `ls *.sorted.bam`
do
  echo $bam
  name=${bam%.sorted.bam}
  java -Xmx40g -jar /media/DISK3TB/Softwares/Picard/picard.jar MarkDuplicates INPUT=$bam OUTPUT=${name}.rmdup.bam METRICS_FILE=test.metrics REMOVE_DUPLICATES=TRUE
done
