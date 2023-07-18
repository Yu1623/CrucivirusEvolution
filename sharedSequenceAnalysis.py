from readFasta import readFasta
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from kmerAnalysis import kmerGenomeCount

def sharedKmersCombinations(fileName1, fileName2, step):
    kmerSequence1 = kmerGenomeCount(fileName1, step)
    kmerSequence2 = kmerGenomeCount(fileName2, step)
    counts = 0
    kmerGenes = 0
    matchingRatio = 0
    sharedKmerSeq = []
    for kmer1 in kmerSequence1:
        for kmer2 in kmerSequence2:
            kmerGenes += 1
            if (kmer1 == kmer2):
                counts += 1
                sharedKmerSeq += [kmer1]
    matchingRatio = counts / kmerGenes
    
    return counts, matchingRatio, sharedKmerSeq

#counts, matchingRatio, sharedKmerSeq = sharedKmersCombinations('newRNA.txt', '879CPs.fasta', 7)
#print(counts, ", ", matchingRatio)

def sharedKmers(fileName1, fileName2, step):
    kmerSequence1 = kmerGenomeCount(fileName1, step)
    kmerSequence2 = kmerGenomeCount(fileName2, step)
    counts = 0
    kmerGenes = 0
    matchingRatio = 0
    for kmer1 in kmerSequence1:
        for kmer2 in kmerSequence2:
            kmerGenes += kmerSequence1[kmer1] + kmerSequence2[kmer2]
            if (kmer1 == kmer2):
                counts += kmerSequence1[kmer1] + kmerSequence2[kmer2]
    matchingRatio = counts / kmerGenes
    return counts, matchingRatio
#counts, matchingRatio = sharedKmers('newRNA.txt', '316CRESS.fasta', 7)
#print(counts, ", ", matchingRatio)

def compareSharedKmer(RNAfile, DNAfile, CPfile, Repfile, step):
    counts1, matchingRatio1 = sharedKmers(RNAfile, DNAfile, step)
    counts2, matchingRatio2 = sharedKmers(CPfile, Repfile, step)
    counts3, matchingRatio3 = sharedKmers(RNAfile, CPfile, step)

    figure, axis = plt.subplots(nrows = 1, ncols = 1)
    gene = ['RNA and DNA', 'CP and Rep', 'RNA and CP']
    ratio = [matchingRatio1, matchingRatio2, matchingRatio3]

    def addValues(x, y):
        for i in range(len(x)):
            plt.text(i, y[i], y[i])
        
    axis.bar(gene, ratio, color = 'lightblue', width = 0.20)
    addValues(gene, ratio)

    plt.show()

#compareSharedKmer('newRNA.txt', '316CRESS.fasta', '855Reps.fasta', '879CPs.fasta', 7)

def interactCompareKmer(RNAfile, DNAfile, Repfile, CPfile, step):
    counts1, ratio1 = sharedKmers(RNAfile, DNAfile, step)
    counts2, ratio2 = sharedKmers(CPfile, Repfile, step)
    counts3, ratio3 = sharedKmers(RNAfile, CPfile, step)

    gene = ['RNA and DNA', 'CP and Rep', 'RNA and CP']
    ratio = [ratio1, ratio2, ratio3]

    figure, axis = plt.subplots(nrows = 1, ncols = 1)

    axisSlider = plt.axes([0.2, 0.01, 0.65, 0.02])
    slider = Slider(axisSlider, 'sval', 2, 8, valinit = step)

    def update(val):
        sVal = slider.val
        newCounts1, newRatio1 = sharedKmers(RNAfile, DNAfile, step)
        newCounts2, newRatio2 = sharedKmers(CPfile, Repfile, step)
        newCounts3, newRatio3 = sharedKmers(RNAfile, CPfile, step)
        ratio = [newRatio1, newRatio2, newRatio3]
        axis.clear()
        axis.bar(gene, ratio, color = 'lightgreen', width = 0.20)
        figure.canvas.draw()
    
    slider.on_changed(update)
    axis.bar(gene, ratio, color = 'lightgreen', width = 0.20)

    plt.show()    

#interactCompareKmer('newRNA.txt', '316CRESS.fasta', '855Reps.fasta', '879CPs.fasta', 7)

def compareDifference(Repfile, CPfile):
    Repfile = readFasta(Repfile)
    CPfile = readFasta(CPfile)
    CPGenes = []
    RepGenes = []
    totalGenes = 0
    percentageOfDifference  = 0
    for gene in CPfile:
        CPGenes += [CPfile[gene]]
    for gene in Repfile:
        RepGenes += [Repfile[gene]]
    differences = 0
    minGenes = []
    if (len(CPGenes) > len(RepGenes)):
        minGenes = RepGenes
    else:
        minGenes = CPGenes
    for i in range(len(minGenes)):
        CPGene = CPGenes[i]
        RepGene = RepGenes[i]
        minGene = ""
        if (len(CPGene) > len(RepGene)):
            minGene = RepGene
        else:
            minGene = CPGene
        for i in range(len(minGene)):
            if (CPGene[i] != RepGene[i]):
                differences += 1
            totalGenes += 1
    percentageOfDifference = differences / totalGenes
    return percentageOfDifference
print("Cruci: ", compareDifference('855Reps.fasta', '879CPs.fasta'), "  RNA: ", compareDifference('1514RdRPs.fasta', '1526RNA_CPs.fasta'), " DNA: ", compareDifference('305CRESS_Reps.fasta', '270CRESS_CPs.fasta'))
                
        
            