from Bio import SeqIO

# seq = ">seq1\n"+"CCTACGACTGAACAGACTCTCCTACGACTG\n"+">seq2\n"+"CCTACGACTTAACAGACTCTCCTACGACTT"

# for record in SeqIO.parse(open('fasta.txt'), "fasta"):
#     print(record)
#         # Add this record to our list
# dna = "CCTACGACTG"     
# comp = dna.translate(str.maketrans("AGCT", "TCGA"))
# lcomp = list(comp)
# # lcomp.reverse()
# print("".join(lcomp))

import Orange as o
learner = o.classification.NaiveBayesLearner()