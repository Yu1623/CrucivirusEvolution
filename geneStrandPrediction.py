from readFasta import readFasta
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display
from tabulate import tabulate

def geneStrandPrediction(fileName, geneName):
    genes = readFasta(fileName)
    gene = genes[geneName]
    numCodons = int(len(gene)/3)
    codonSequence = []
    for i in range(numCodons):
        codon = gene[i * 3 : i * 3 + 3]
        codonSequence += [codon]
    nna = 0
    nnt = 0
    nng = 0
    nnc = 0
    for codon in codonSequence:
        codonEnding = codon[-1]
        if (codonEnding == "A"):
            nna += 1
        elif (codonEnding == "T"):
            nnt += 1
        elif (codonEnding == "C"):
            nnc += 1
        elif (codonEnding == "G"):
            nng += 1
    percenta = nna / numCodons
    percentt = nnt / numCodons
    percentg = nng / numCodons
    percentc = nnc / numCodons
    bias = ""
    if (percenta > percentt):
        if (percentg > percentc):
            if (percenta > percentg):
                bias = "A"
            else:
                bias = "G"
        else:
            if (percenta > percentc):
                bias = "A"
            else:
                bias = "C"
    else:
        if (percentg > percentc):
            if (percentt > percentg):
                bias = "T"
            else:
                bias = "G"
        else:
            if (percentt > percentc):
                bias = "T"
            else:
                bias = "C"
    return bias

def geneStrandsPrediction(fileName):
    genes = readFasta(fileName)
    biases = {}
    for geneName in genes:
        bias = geneStrandPrediction(fileName, geneName)
        biases[geneName] = bias
    return biases

def genomeSensePrediction(CPfileName, RepfileName):
    CPbiases = geneStrandsPrediction(CPfileName)
    Repbiases = geneStrandsPrediction(RepfileName)  
    CPgenes = []
    Repgenes = []
    genomeNames = [] 
    GenomeSense = []
    ambisense = 0
    for CPgeneName in CPbiases:
        for RepgeneName in Repbiases:
            if (CPgeneName.split("_-_")[0] == RepgeneName.split("_-_")[0]):
                genomeSense = ""
                CPbias = CPbiases[CPgeneName]
                Repbias = Repbiases[RepgeneName]
                if (CPbias == Repbias):
                    genomeSense = "unisense"
                else:
                    genomeSense = "ambisense"
                    ambisense += 1
                CPgenes += [CPbiases[CPgeneName]]
                Repgenes += [Repbiases[RepgeneName]]
                GenomeSense += [genomeSense]
                genomeNames += [CPgeneName.split("_-_")[0]]
                percentAmbisense = ambisense / len(CPbiases)
    return percentAmbisense, CPgenes, Repgenes, GenomeSense, genomeNames
#print(geneStrandsPrediction("22CPs.fasta"))
#print(geneStrandsPrediction("22Reps.fasta"))
#print(genomeSensePrediction("22CPs.fasta", "22Reps.fasta"))

percentAmbisense, CPgenes, Repgenes, GenomeSense, genomeNames = genomeSensePrediction("879CPs.fasta", "855Reps.fasta")

rows = []
for genome in genomeNames:
    rows += [genome]
data = np.array([[column for column in range(4)] for row in range(len(rows))], dtype = object)
for i in range(len(data[:, 0])):
    data[i, 0] = rows[i]
for i in range(len(data[:, 1])):
    data[i, 1] = CPgenes[i]
for i in range(len(data[:, 2])):
    data[i, 2] = Repgenes[i]
for i in range(len(data[:, 3])):
    data[i, 3] = GenomeSense[i]
columns = ['Genome name', 'CP', 'Rep', 'Genome']
dict = {'Genome Name: ': genomeNames, 'CP': CPgenes, 'Rep': Repgenes, 'Genome': GenomeSense}

df = pd.DataFrame(dict)
pd.DataFrame(dict).to_csv('GenomeSensePrediction.csv')
print(df.to_markdown())
