'''
This program analyzes the kmer rank between CP and Rep gene of the same genome. 
Please refer to the KmerRankCPRepVisual file to see the latest version of kmer rank analysis.
'''
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
import numpy as np

def readRanking(fileName):
    f = open(fileName, 'r')
    kmerRank = {}
    lines = f.readlines()
    for i in range(len(lines)):
        line = lines[i]
        geneName = line.split(": ")[0]
        gene = line.split(": ")[1]
        kmerRank[geneName] = gene
    return kmerRank

kmerRanks = ["kmerRankCruci.txt", "kmerRankDNA.txt", "kmerRankRNA.txt"]
percentsCommon = {}
geneNames = {}
for kmerRank in kmerRanks:
    kmerrank = readRanking(kmerRank)
    percentCommon = []
    geneName = []
    for kmer in kmerrank:
        percentCommon += [kmerrank[kmer]]
        geneName += [kmer]
    df = pd.DataFrame(percentCommon)
    fig, axis = plt.subplots(nrows=1, ncols=1, figsize=(50, 42))
    temp = hierarchy.linkage(df, method = 'ward')
    axis.set_title(kmerRank)
    axis.set_xlabel('Genomes')
    dn = hierarchy.dendrogram(temp, ax=axis, above_threshold_color='green', color_threshold=0.7, labels=geneName, leaf_font_size=2)
    plt.show()
    percentsCommon[kmerRank] = df
    geneNames[kmerRank] = geneName

fig, (axis1, axis2, axis3) = plt.subplots(nrows = 1, ncols = 3, figsize = (25, 20))
temp1 = hierarchy.linkage(percentsCommon['kmerRankCruci.txt'], method = 'ward')
axis1.set_title('Cruci')
axis1.set_xlabel('Genomes')
dn1 = hierarchy.dendrogram(temp1, ax=axis1, above_threshold_color='green', color_threshold=0.7)
axis1.set_ylim(0, 0.00078)

temp2 = hierarchy.linkage(percentsCommon["kmerRankRNA.txt"], method = 'ward')
axis2.set_title('RNA')
axis2.set_xlabel('Genomes')
dn2 = hierarchy.dendrogram(temp2, ax=axis2, above_threshold_color='blue', color_threshold=0.7)
axis2.set_ylim(0, 0.00078)

temp3 = hierarchy.linkage(percentsCommon["kmerRankDNA.txt"], method = 'ward')
axis3.set_title('DNA')
axis3.set_xlabel('Genomes')
dn3 = hierarchy.dendrogram(temp3, ax=axis3, above_threshold_color='gray', color_threshold=0.7)
axis3.set_ylim(0, 0.00078)

plt.show()


