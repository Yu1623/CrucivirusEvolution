from readFasta import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

'''
Input: file name and gene name
Output: codons in the gene
'''
def codonSequence(fileName, geneName):
    file = readFasta(fileName)
    gene = file[geneName]
    numCodons = int(len(gene) - 2)
    codons = []
    for codon in range(numCodons):
        codonSeq = gene[codon:codon + 3]
        codons += [codonSeq]
    return codons

'''
Input: file name and gene name
Output: the percentage of nna, nnt, nng, nnc
'''
def codonEndingFrequency(fileName, geneName):
    codons = codonSequence(fileName, geneName)
    endingA = 0
    endingT = 0
    endingG = 0
    endingC = 0

    for codon in codons:
        codonEnding = codon[-1]
        if (codonEnding == 'A'):
            endingA += 1
        elif (codonEnding == 'T'):
            endingT += 1
        elif (codonEnding == 'G'):
            endingG += 1
        elif (codonEnding == 'C'):
            endingC += 1
    percentA = endingA / len(codons)
    percentT = endingT / len(codons)
    percentG = endingG / len(codons)
    percentC = endingC / len(codons)
    return percentA, percentT, percentG, percentC

'''
Input: file name
Output: the arrays of percentages of A, T, G, and C
'''
def codonPreference(fileName):
    file = readFasta(fileName)
    percentagesA = []
    percentagesT = []
    percentagesG = []
    percentagesC = []
    for gene in file:
        percentA, percentT, percentG, percentC = codonEndingFrequency(fileName, gene)
        percentagesA += [percentA]
        percentagesT += [percentT]
        percentagesG += [percentG]
        percentagesC += [percentC]
    return percentagesA, percentagesT, percentagesG, percentagesC
    
'''
Box plot of percentages nna, nnt, nng, and nnc
'''
def codonVisual(fileName):
    percentA, percentT, percentG, percentC = codonPreference(fileName)
    percentNucleotide = [percentA, percentT, percentG, percentC]

    figure = plt.figure(figsize = (10, 7))
    axis = figure.add_subplot(111)
    codonPlot = axis.boxplot(percentNucleotide, patch_artist = True, vert = 0)
    axis.set_yticklabels(['Adenine', 'Thymine', 'Guamine', 'Cytosine'])
    axis.set_xlabel("Percentage of Codon ending Usage")
    plt.title("%s Percentages of nucleotides as endings of codons" % (fileName))
    axis.get_xaxis().tick_bottom()
    axis.get_yaxis().tick_left()
    plt.show()

#codonVisual('22Reps.fasta')

def codonCompareVisual(fileName1, fileName2, label1, label2):
    percentA1, percentT1, percentG1, percentC1 = codonPreference(fileName1)
    percentA2, percentT2, percentG2, percentC2 = codonPreference(fileName2)
    datasetA = [percentA1, percentA2]
    datasetT = [percentT1, percentT2]
    datasetG = [percentG1, percentG2]
    datasetC = [percentC1, percentC2]
    labels = [label1, label2]

    figure, (axis1, axis2, axis3, axis4) = plt.subplots(nrows = 1, ncols = 4, figsize = (9, 4), sharey = True)
    boxplot1 = axis1.boxplot(datasetA, vert = True, patch_artist = True, labels = labels)
    axis1.set_title('nna')
    boxplot2 = axis2.boxplot(datasetT, vert = True, patch_artist = True, labels = labels)
    axis2.set_title('nnt')
    boxplot3 = axis3.boxplot(datasetG, vert = True, patch_artist = True, labels = labels)
    axis3.set_title('nng')
    boxplot4 = axis4.boxplot(datasetC, vert = True, patch_artist = True, labels = labels)
    axis4.set_title('nnc')

    plt.tight_layout()

    colors = ['lightblue', 'gray']
    for boxplot in (boxplot1, boxplot2, boxplot3, boxplot4):
        for patch, color in zip(boxplot['boxes'], colors):
            patch.set_facecolor(color)

    plt.show()