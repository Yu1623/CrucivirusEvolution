import re
import sys
from readFasta import readFasta
def writeFasta(genes, title):
    file = open("%s.txt" % (title), 'w')
    for gene in genes:
        file.write(">" + gene + "\n" + genes[gene] + "\n")
    file.close()

def combineFasta(fileName1, fileName2, title):
    file = open("%s.txt" % (title), 'w')
    file1 = readFasta(fileName1)
    file2 = readFasta(fileName2)
    for gene in file1:
        file.write(">" + gene + "\n" + file1[gene] + "\n")
    for gene in file2:
        file.write(">" + gene + "\n" + file2[gene] + "\n")

combineFasta("305CRESS_Reps.fasta", "DNARepRev.txt", "DNARepCombine")
