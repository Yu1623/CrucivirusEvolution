from readFasta import *
from kmerAnalysis import countGenome
import matplotlib.pyplot as plt
import numpy as np

def revComplement(fileName, geneName):
    file = readFasta(fileName)
    gene = file[geneName]
    newGene = gene.replace("A", "N").replace("T", "A").replace("N", "T").replace("G", "L").replace("C", "G").replace("L", "C")
    return newGene

def revGenome(fileName):
    file = readFasta(fileName)
    newGenome = {}
    for gene in file:
        newGene = revComplement(fileName, gene)
        if (gene in newGenome.keys()):
            newGenome[gene] += newGene
        else:
            newGenome[gene] = newGene
    return newGenome

def countGenomeVisualrevCompare(fileName1, fileName2, step):
    countsArray = countGenome(fileName1, step)
    revcountsArray = countGenome(fileName2, step)
    counts = []
    frequencies = []
    for count in countsArray:
        counts += [count]
        frequency = countsArray[count]
        frequencies += [frequency]
    for count in revcountsArray:
        counts += [count]
        frequency = revcountsArray[count]
        frequencies += [frequency]
    plt.scatter(frequencies, counts)
    plt.xlabel("%s-mer frequency" % (step))
    plt.ylabel("Count")
    plt.show()

countGenomeVisualrevCompare('855Reps.fasta', 'rev855Reps.txt', 8)