# SARS-CoV-2 Consensus
 Tewhey Lab pipeline for consensus calling of SARS-CoV-2 samples

*sars_cov_2_consensus.sh*

Creates individual bash script for calling consensus for each sample provided.

The user should have the trimmomatic jar available to them and change that directory to where it lives.
The user should also have modules for java, bwa, samtools and ivar available.

The IDX file is a bwa index file generated from a single file comprised of the human genome and transcriptome references.
The REF file is the SARS-CoV2 reference fasta.
