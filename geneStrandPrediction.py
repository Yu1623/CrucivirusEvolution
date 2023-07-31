'''
This program predicts the direction, represented in codon-ending bias, of CP and Rep genes. 
The direction of CP and Rep genes can help predict the sense of the genome. The name of the genome, the direction of CP, the direction of Rep, and the sense of genome can be downloaded as a csv file and seen in the terminal.
'''

from readFasta import readFasta
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display
from tabulate import tabulate

'''
Input: file name and the gene name that makes the path to the gene that the program will predict the orientation of.
Output: the codon ending bias, which indicates the orientation, and the percent of codons ending in a and the percent ending in t
'''

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
    return bias, percenta, percentt

'''
Input: file name of the genes that will be predicted
Output: an associative array giving each gene and their corresponding codon-ending preference.
'''

def geneStrandsPrediction(fileName):
    genes = readFasta(fileName)
    biases = {}
    for geneName in genes:
        bias, percenta, percentt = geneStrandPrediction(fileName, geneName)
        biases[geneName] = bias
    return biases

'''
Input: CP file name and Rep file name
Output: an associative array that stores the genome names and the sense of the genomes.
'''

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

def genomeSenseTable(CPfile, Repfile):
    CPgenes, Repgenes, GenomeSense, genomeNames, CPending, Repending = genomeSensePrediction(CPfile, Repfile)
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
    dict = {'Genome Name': genomeNames, 'CP': CPgenes, 'Rep': Repgenes, 'Genome': GenomeSense}
    df = pd.DataFrame(dict)
    df.to_markdown()
    return dict

'''
Input: file name of the csv file that has the actual sense of the genomes.
Output: a matrix with the information from the csv file.
'''

def readCSV(fileName):
    f = open(fileName, 'r')
    lines = f.readlines()
    names = []
    sequenceNames = []
    lengths = []
    directions = []
    documentNames = []
    types = []
    sources = []
    for i in range(len(lines)):
        line = lines[i]
        if (i == 0):
            continue
        if ("Rep" in line.split(",")[0] or "CP" in line.split(",")[0] or "Capsid" in line.split(",")[0]):
            name = line.split(",")[0]
            sequenceName = line.split(",")[1]
            length = line.split(",")[2]
            direction = line.split(",")[3]
            documentName = line.split(",")[4]
            type = line.split(",")[5]
            source = line.split(",")[6]
            names += [name]
            sequenceNames += [sequenceName]
            lengths += [length]
            directions += [direction]
            documentNames += [documentName]
            types += [type]
            sources += [source]
    dict = {"Name": names, "Sequence Name": sequenceNames, "Length": lengths, "Direction": directions, "Document Name": documentNames, "Type": types, "Source": sources}
    df = pd.DataFrame(dict)
    return dict

'''
Checks for the accuracy of the predictions of the genome sense and outputs the accuracy and an array of genomes that are have wrong predictions.
'''

def checkPrediction(csvfileName, CPfileName, RepfileName):
    dict1 = readCSV(csvfileName)
    dict2 = genomeSenseTable(CPfileName, RepfileName)

    commonSense = 0
    senseTotal = 0
    incorrectPrediction = []
    for i in range(len(dict1["Name"])):
        if ("CP" in dict1["Name"][i] or "Capsid" in dict1["Name"][i]):
            for l in range(len(dict1["Name"])):
                if ("Rep" in dict1["Name"][l]):
                    if (dict1["Sequence Name"][i] == dict1["Sequence Name"][l]):
                        directionCP = dict1["Direction"][i]
                        directionRep = dict1["Direction"][l]
                        if (directionCP == directionRep):
                            genomeSense = "unisense"
                        else:
                            genomeSense = "ambisense"
                        for m in range(len(dict2["Genome Name"])):
                            if (dict2["Genome Name"][m] == dict1["Sequence Name"][l]):
                                genomeSensePrediction = dict2["Genome"][m]
                                senseTotal += 1
                                if (genomeSense == genomeSensePrediction):
                                    commonSense += 1
                                else:
                                    incorrectPrediction += [dict1["Sequence Name"][l]]
    accuracy = commonSense/senseTotal
    return accuracy, incorrectPrediction

