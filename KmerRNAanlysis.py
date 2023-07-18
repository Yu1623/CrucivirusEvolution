from readFasta import readFasta
from writeFasta import writeFasta
import matplotlib.pyplot as plt
import numpy as np

'''
Make RNA virus genome file (removing CP and reverse complements)
'''

genes = readFasta('1626RNAviruses.fasta')
cpGenes = readFasta('1526RNA_CPs.fasta')
newGenes = {}
for gene in genes:
    geneSeq = genes[gene]
    removeSeq = []
    for cpGene in cpGenes:
        cpSeq = cpGenes[cpGene]
        if (cpSeq in geneSeq):
            removeSeq += [cpSeq]
    for seq in removeSeq:
        geneSeq = geneSeq.replace(seq, '')
    if (gene in newGenes.keys()):
        newGenes[gene] += geneSeq
    else:
        newGenes[gene] = geneSeq
writeFasta(newGenes, "newRNA")
    