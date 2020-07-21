#!/bin/bash

while read line
do
R1=`echo $line | awk '{print $2}'`;
R2=`echo $line | awk '{print $3}'`;
ID=`echo $line | awk '{print $1}'`;
REF="/projects/tewhey-lab/projects/COVID/reference_files/NC_045512.fa"
IDX="/projects/tewhey-lab/projects/COVID/reference_files/hg_38_gencode34_idx"

cp slurm/header.txt slurm/slurm.${ID}.runCMD.sh
echo "java -jar /projects/tewher/bin/Trimmomatic-0.38/trimmomatic-0.38.jar PE $R1 $R2 ../reads/${ID}.R1.trimmed.fastq ../reads/${ID}.R1.unmated.fastq ../reads/${ID}.R2.trimmed.fastq ../reads/${ID}.R2.unmated.fastq ILLUMINACLIP:/projects/tewher/bin/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE MINLEN:25 &> ${ID}.trim.log" >> slurm/slurm.${ID}.runCMD.sh

echo "bwa aln -t 4 $IDX ../reads/${ID}.R1.trimmed.fastq > ../mapping/${ID}.R1.sai" >> slurm/slurm.${ID}.runCMD.sh
echo "bwa aln -t 4 $IDX ../reads/${ID}.R2.trimmed.fastq > ../mapping/${ID}.R2.sai" >> slurm/slurm.${ID}.runCMD.sh
echo "bwa sampe $IDX ../mapping/${ID}.R1.sai ../mapping/${ID}.R2.sai $R1 $R2 > ../mapping/${ID}_human.sam" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools fastq -1 ../mapping/${ID}_filtered.R1.fastq -2 ../mapping/${ID}_filtered.R2.fastq -0 /dev/null -n -f 0x4 ../mapping/${ID}_human.sam" >> slurm/slurm.${ID}.runCMD.sh

echo "bwa mem -k 12 -B 1 -t 12 $REF ../mapping/${ID}_filtered.R1.fastq ../mapping/${ID}_filtered.R2.fastq | samtools view -u - | samtools sort -O BAM -o ../mapping/${ID}.bam -" >> slurm/slurm.${ID}.runCMD.sh
echo "/projects/tewhey-lab/rtewhey/COVID/bin/samtools/samtools ampliconclip --both-ends --strand  --filter-len 20  --no-excluded -b /projects/tewhey-lab/projects/COVID/reference_files/artic_primers_v3.bed ../mapping/${ID}.bam |samtools view -u - | samtools sort -O BAM -o ../mapping/${ID}.clipped.bam" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools index ../mapping/${ID}.clipped.bam" >> slurm/slurm.${ID}.runCMD.sh

echo "samtools mpileup -A -d 0 -Q 0 -B ../mapping/${ID}.clipped.bam | ivar consensus -p ../mapping/${ID}.consensus" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools mpileup -A -d 0 -Q 0 --reference $REF ../mapping/${ID}.clipped.bam | ivar variants -g ${REF%%.fa}.gff -r $REF -p mapping/${ID}.consensus -t 0.05" >> slurm/slurm.${ID}.runCMD.sh
done < samples.list

for i in `ls slurm/*runCMD.sh`; do echo $i; sbatch $i; done