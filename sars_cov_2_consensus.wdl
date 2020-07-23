# WDL for mapping and getting consensus for covid samples

workflow CovidMap {
    Array[File] read_1_raw #Array of full paths to read 1 fastq files
    Array[File] read_2_raw #Array of full paths to read 2 fastq files
    Array[String] id_out #Array of sample IDs, should be in the same order as fastq files
    Array[Pair[File,File]] to_trim = zip(read_1_raw,read_2_raw)
    Array[Pair[String,Pair[File,File]]] all_info = zip(id_out,to_trim)
    File primer_bed
    File reference_fasta
    File reference_gff
    String IDX #Path, including prefix, of index files to be used for bwa
    String results_directory #Path to where you want overall results
    String script_directory #Path to coverage_hist.R
    String trimmomatic_jar #Full path to trimmomatic jar location
    String samtools_git #Full path to cloned samtools git repository e.x. /projects/tewhey-lab/rtewhey/COVID/bin/samtools/samtools

    scatter (sample in all_info) {
      Pair[File,File] reads = sample.right

      call Trim { input:
                       read_1_raw = reads.left,
                       read_2_raw = reads.right,
                       trimmomatic_jar = trimmomatic_jar,
                       id_out = sample.left
                     }
      call MapHuman { input:
                            trimmed_1 = Trim.outR1,
                            trimmed_2 = Trim.outR2,
                            IDX = IDX,
                            id_out = sample.left
                          }
      call MapBam { input:
                          filtered_1 = MapHuman.R1_filtered,
                          filtered_2 = MapHuman.R2_filtered,
                          reference_fasta = reference_fasta,
                          id_out = sample.left
                      }
      call ClipPrimer { input:
                              primer_bed = primer_bed,
                              full_bam = MapBam.out,
                              id_out = sample.left,
                              samtools_git = samtools_git
                            }
      call Index { input:
                        clipped_bam = ClipPrimer.out,
                        id_out = sample.left
                      }
      call Consensus { input:
                            clipped_bam = ClipPrimer.out,
                            id_out = sample.left
                          }
      call ConsensusGFF { input:
                                clipped_bam = ClipPrimer.out,
                                reference_fasta = reference_fasta,
                                reference_gff = reference_gff,
                                id_out = sample.left
                              }
      call get_depth { input:
                          bam=ClipPrimer.out,
                          id_out=sample.left
                      }
      call make_hist { input:
                          depth_hist=get_depth.out,
                          id_out=sample.left,
                          working_directory=script_directory
                      }
      call ConveneFiles { input:
                                full_bam = MapBam.out,
                                clipped_bam = ClipPrimer.out,
                                index_bam = Index.out,
                                consensus_fasta = Consensus.out1,
                                consensus_quality = Consensus.out2,
                                consensus_tsv = ConsensusGFF.out,
                                qc_plot = make_hist.out,
                                results_directory = results_directory
                              }
    }
  }

task Trim {
    File read_1_raw
    File read_2_raw
    String trimmomatic_jar
    String id_out
    command {
        java -jar ${trimmomatic_jar} PE ${read_1_raw} ${read_2_raw} ${id_out}_R1_001.trimmed.fastq ${id_out}_R1_001.unmated.fastq ${id_out}_R2_001.trimmed.fastq ${id_out}_R2_001.unmated.fastq ILLUMINACLIP:/projects/tewher/bin/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10:2:TRUE MINLEN:25
      }
    runtime { runtime_minutes: 15
              requested_memory_mb_per_core: 900
              cpus: 1}
    output {
        File outR1 = "${id_out}_R1_001.trimmed.fastq"
        File outR2 = "${id_out}_R2_001.trimmed.fastq"
      }
  }

task MapHuman {
  File trimmed_1
  File trimmed_2
  String IDX
  String id_out
  command {
    bwa aln -t 4 ${IDX} ${trimmed_1} > ${id_out}.R1.sai
    bwa aln -t 4 ${IDX} ${trimmed_2} > ${id_out}.R2.sai
    bwa sampe ${IDX} ${id_out}.R1.sai ${id_out}.R2.sai ${trimmed_1} ${trimmed_2} > ${id_out}_human.sam
    samtools fastq -1 ${id_out}_filtered.R1.fastq -2 ${id_out}_filtered.R2.fastq -0 /dev/null -s /dev/null -n -f 0x4 ${id_out}_human.sam
   }
  runtime { runtime_minutes: 15
             requested_memory_mb_per_core: 900
             cpus: 1}
  output {
    File R1_filtered = "${id_out}_filtered.R1.fastq"
    File R2_filtered = "${id_out}_filtered.R2.fastq"
   }
 }

task MapBam {
  File filtered_1
  File filtered_2
  File reference_fasta
  String id_out
  command {
      bwa mem -k 12 -B 1 -t 12 ${reference_fasta} ${filtered_1} ${filtered_2} | samtools view -u - | samtools sort -O BAM -o ${id_out}.bam -
    }
  runtime { runtime_minutes: 15
            requested_memory_mb_per_core: 900
            cpus: 1}
  output {
      File out = "${id_out}.bam"
    }
  }

task ClipPrimer {
  File primer_bed
  File full_bam
  String id_out
  String samtools_git
  command {
      ${samtools_git} ampliconclip --both-ends --strand --filter-len 20 --no-excluded -b ${primer_bed} ${full_bam} | samtools view -u - | samtools sort -O BAM -o ${id_out}.clipped.bam
    }
  runtime { runtime_minutes: 15
            requested_memory_mb_per_core: 900
            cpus: 1}
  output {
      File out = "${id_out}.clipped.bam"
    }
  }

task Index {
    File clipped_bam
    String id_out
    command {
        samtools index ${clipped_bam}
      }
    runtime { runtime_minutes: 15
              requested_memory_mb_per_core: 900
              cpus: 1}
    output {
        File out = "${id_out}.clipped.bam.bai"
      }
  }

task Consensus {
    File clipped_bam
    String id_out
    command {
        samtools mpileup -A -d 0 -Q 0 -B ${clipped_bam} | ivar consensus -p ${id_out}.consensus
      }
    runtime { runtime_minutes: 15
              requested_memory_mb_per_core: 900
              cpus: 1}
    output {
        File out1 = "${id_out}.consensus.fa"
        File out2 = "${id_out}.consensus.qual.txt"
      }
  }

task ConsensusGFF {
    File clipped_bam
    File reference_fasta
    File reference_gff
    String id_out
    command {
        samtools mpileup -A -d 0 -Q 0 --reference ${reference_fasta} ${clipped_bam} | ivar variants -g ${reference_gff} -r ${reference_fasta} -p ${id_out}.consensus -t 0.05
      }
    runtime { runtime_minutes: 15
              requested_memory_mb_per_core: 900
              cpus: 1}
    output {
        File out = "${id_out}.consensus.tsv"
      }
  }

task ConveneFiles {
    File full_bam
    File clipped_bam
    File index_bam
    File consensus_fasta
    File consensus_quality
    File consensus_tsv
    File qc_plot
    String results_directory
    command {
        rsync -v ${full_bam} ${clipped_bam} ${index_bam} ${consensus_fasta} ${consensus_quality} ${consensus_tsv} ${qc_plot} ${results_directory}/
      }
  }

task get_depth {
    File bam
    String id_out
    command {
        samtools view -b ${bam} | genomeCoverageBed -d -ibam stdin > ${id_out}.hist
      }
    runtime { runtime_minutes: 15
              requested_memory_mb_per_core: 900
              cpus: 1}
    output {
        File out="${id_out}.hist"
      }
  }
task make_hist {
    File depth_hist
    String id_out
    String working_directory
    command {
        module load R/3.6.0
        Rscript ${working_directory}/coverage_hist.R ${id_out} ${depth_hist}
      }
    runtime { runtime_minutes: 15
              requested_memory_mb_per_core: 900
              cpus: 1}
    output {
        File out="${id_out}_log_depth.png"
      }
  }
