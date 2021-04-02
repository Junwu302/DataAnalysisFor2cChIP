#!/bin/bash
for var in `ls *_R1.fastq.gz`
do
  name=${var%_R1.fastq.gz}
  echo $name
  java -jar /media/DISK3TB/Softwares/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 6 \
  ${name}_R1.fastq.gz ${name}_R2.fastq.gz \
  ${name}_1.qc.fq ${name}_1.se.fq.gz \
  ${name}_2.qc.fq ${name}_2.se.fq.gz \
  ILLUMINACLIP:/media/DISK3TB/Softwares/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 2>&1 | tee ${name}.log
  rm ${name}_1.se.fq.gz ${name}_2.se.fq.gz
  bowtie2 -p 16 --no-unal --no-mixed --no-discordant -x /media/DISK3TB/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -1 ${name}_1.qc.fq -2 ${name}_2.qc.fq -S $name.sam 2>&1 | tee -a ${name}.log
  samtools view -bS -o ${name}.bam $name.sam
  samtools sort -o ${name}.sorted.bam ${name}.bam
  samtools index ${name}.sorted.bam
  rm ${name}.sam ${name}.bam
  java -Xmx40g -jar /media/DISK3TB/Softwares/Picard/picard.jar MarkDuplicates INPUT=${name}.sorted.bam OUTPUT=${name}.rmdup.bam METRICS_FILE=${name}.metrics REMOVE_DUPLICATES=TRUE
  samtools index ${name}.rmdup.bam
  N=10000000
  n="$(samtools view -c  ${name}.rmdup.bam  2>&1 )"
  sf=$(printf "%.5f" `echo "scale=5;$N/$n"|bc`)
  bamCoverage --scaleFactor $sf -b ${name}.rmdup.bam -e 300 -bs 100 --smoothLength 1000 -o ${name}.10M.bw
done