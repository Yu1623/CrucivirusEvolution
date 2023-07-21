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
kmerRank = readRanking("kmerRankCruci.txt")
percentCommon = []
geneName = []
for kmer in kmerRank:
    percentCommon += [kmerRank[kmer]]
    geneName += [kmer]
df = pd.DataFrame(percentCommon)
print(df)

temp = hierarchy.linkage(df, method = 'ward')
plt.figure()
dn = hierarchy.dendrogram(temp, above_threshold_color='green', color_threshold=0.7, labels = geneName)
plt.show()


