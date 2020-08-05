import sys
argv = sys.argv

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

fastafile = argv[1]
to_trim = int(argv[2])
trimmed_name = argv[3]

with open("%s.fa" % trimmed_name, "w") as trimmed_fasta:
    with open("%s_length_check.txt" % trimmed_name, "w") as trim_length:
        with open(fastafile, "rt") as handle:
          for record in SeqIO.parse(handle, "fasta"):
              seq_only = record.seq
              seq_only = str(seq_only)

              seq_len = len(seq_only)
              end_trim = seq_len-to_trim
              new_seq = seq_only[to_trim:end_trim]
              record.seq = Seq(new_seq)
              new_len = len(new_seq)

              SeqIO.write(record, trimmed_fasta, "fasta")
              trim_length.write("%s\t%d\t%d\t%d\n" % (record.name, seq_len, seq_len-new_len, new_len))
