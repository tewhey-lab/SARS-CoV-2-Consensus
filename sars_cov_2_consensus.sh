#!/bin/bash

while read line
do
R1=`echo $line | awk '{print $2}'`;
R2=`echo $line | awk '{print $3}'`;
ID=`echo $line | awk '{print $1}'`;
REF="/projects/tewhey-lab/projects/COVID/reference_files/NC_045512.fa"
IDX="/projects/tewhey-lab/projects/COVID/reference_files/hg38_gencode34_cov2_spike96_idx2"

echo $ID

cp slurm/header.txt slurm/slurm.${ID}.runCMD.sh
echo "java -jar /projects/tewher/bin/Trimmomatic-0.38/trimmomatic-0.38.jar PE $R1 $R2 reads/${ID}.R1.trimmed.fastq reads/${ID}.R1.unmated.fastq reads/${ID}.R2.trimmed.fastq reads/${ID}.R2.unmated.fastq ILLUMINACLIP:/projects/tewher/bin/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE MINLEN:25 &> reads/${ID}.trim.log" >> slurm/slurm.${ID}.runCMD.sh

echo "bwa mem -t 4 $IDX reads/${ID}.R1.trimmed.fastq reads/${ID}.R2.trimmed.fastq > mapping/${ID}_human.sam" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools view -F 256 mapping/${ID}_human.sam | awk '{ if (\$3==\"*\" || \$3==\"NC_045512.2\") {print \$1} }' > mapping/${ID}_to_keep.txt" >> slurm/slurm.${ID}.runCMD.sh
echo "picard FilterSamReads -Xmx3g I=mapping/${ID}_human.sam O=mapping/${ID}_human_filtered.sam READ_LIST_FILE=mapping/${ID}_to_keep.txt FILTER=includeReadList" >>slurm/slurm.${ID}.runCMD.sh
echo "samtools fastq -1 mapping/${ID}_filtered.R1.fq -2 mapping/${ID}_filtered.R2.fq -0 /dev/null -s /dev/null -n mapping/${ID}_human_filtered.sam" >> slurm/slurm.${ID}.runCMD.sh

# echo "bwa mem -k 12 -B 1 -t 12 $REF mapping/${ID}_filtered.R1.fq mapping/${ID}_filtered.R2.fq | samtools view -u - | samtools sort -O BAM -o mapping/${ID}.bam -" >> slurm/slurm.${ID}.runCMD.sh
echo "minimap2 -K 20M -x sr -a $REF mapping/${ID}_filtered.R1.fq mapping/${ID}_filtered.R2.fq | samtools view -u -h -F 4 - | samtools sort -O BAM -o mapping/${ID}.bam -" >> slurm/slurm.${ID}.runCMD.sh

##Remove soft-clipped reads
echo "samtools view mapping/${ID}.bam | awk '(\$6 ~ /H|S/){print \$1}' | sort -k1,1 | uniq > mapping/${ID}_names.txt" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools view mapping/${ID}.bam | sort -k1,1 > mapping/${ID}_names_tmp.sam" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools view -H mapping/${ID}.bam > mapping/${ID}_ClipRemoved.sam" >> slurm/slurm.${ID}.runCMD.sh
echo "join -t $'\t' -v 1 -1 1 -2 1 mapping/${ID}_names_tmp.sam mapping/${ID}_names.txt >> mapping/${ID}_ClipRemoved.sam" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools view -bS mapping/${ID}_ClipRemoved.sam > mapping/${ID}_ClipRemoved.bam" >> slurm/slurm.${ID}.runCMD.sh
echo "rm mapping/${ID}_names_tmp.sam mapping/${ID}_ClipRemoved.sam" >> slurm/slurm.${ID}.runCMD.sh

echo "/projects/tewhey-lab/rtewhey/COVID/bin/samtools/samtools ampliconclip --both-ends --tolerance 1 --filter-len 20  --no-excluded -b /projects/tewhey-lab/projects/COVID/reference_files/artic_primers_v3.bed mapping/${ID}_ClipRemoved.bam |samtools view -u - | samtools sort -O BAM -o mapping/${ID}.clipped.bam" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools index mapping/${ID}.clipped.bam" >> slurm/slurm.${ID}.runCMD.sh

echo "samtools view -b mapping/${ID}.clipped.bam | genomeCoverageBed -d -ibam stdin > QC/${ID}.hist" >> slurm/slurm.${ID}.runCMD.sh
echo "Rscript /projects/tewhey-lab/projects/COVID/scripts/coverage_hist.R QC/${ID} QC/${ID}.hist" >> slurm/slurm.${ID}.runCMD.sh

echo "samtools mpileup -A -d 0 -Q 0 -B mapping/${ID}.clipped.bam | ivar consensus -t 0 -p working_consensus/${ID}.consensus" >> slurm/slurm.${ID}.runCMD.sh
echo "samtools mpileup -A -d 0 -Q 0 --reference $REF mapping/${ID}.clipped.bam | ivar variants -g ${REF%%.fa}.gff -r $REF -p working_consensus/${ID}.consensus -t 0.05" >> slurm/slurm.${ID}.runCMD.sh
echo "python /projects/tewhey-lab/projects/COVID/scripts/trim_ends_fasta.py working_consensus/${ID}.consensus.fa 80 working_consensus/${ID}.consensus.trimmed" >> slurm/slurm.${ID}.runCMD.sh
done < 20210303_samples_list.txt

for i in `ls slurm/*runCMD.sh`; do echo $i; sbatch $i; done
