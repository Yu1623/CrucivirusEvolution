import re
import sys
def writeFasta(genes, title):
    file = open("%s.txt" % (title), 'w')
    for gene in genes:
        file.write(">" + gene + "\n" + genes[gene] + "\n")
    file.close()