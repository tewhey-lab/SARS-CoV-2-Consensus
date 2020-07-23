# SARS-CoV-2 Consensus
 Tewhey Lab pipeline for consensus calling of SARS-CoV-2 samples

*sars_cov_2_consensus.sh*

Creates individual bash script for calling consensus for each sample provided. The format for providing samples should be document with three columns, with the first column being listed as the sample name, the second column should be the full path to the read 1 fastq file, and the third column should be the full path to the read 2 fastq file.

The user should have the trimmomatic jar available to them and change that directory to where it lives.
The user should also have modules for `java`, `bwa`, `samtools` and `ivar` available. It is also important to note that at this point `ampliconclip` is not in the version 1.9 conda release of samtools. It is available through the git repository should be cloned and the path modified.

The IDX file is a bwa index file generated from a single file comprised of the human genome and transcriptome references.
The REF file is the SARS-CoV2 reference fasta.
