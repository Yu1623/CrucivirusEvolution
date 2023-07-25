'''
This file analyze kmer rank for the similarity of CP and Rep genes of genomes.
The program makes a line graph that visualizes the frequency of the similarity percent between CP and Rep genes.
'''

import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import pandas as pd

def readKmerRanking(fileName):
    f = open(fileName, 'r')
    kmerRank = {}
    lines = f.readlines()
    for line in lines:
        genome = line.split(": ")[0]
        percentSimilar = line.split(": ")[1]
        kmerRank[genome] = percentSimilar
    return kmerRank

def kmerRankShaded(fileName):
    kmerRank = readKmerRanking(fileName)
    percentsCounts = {}
    for genome in kmerRank:
        percentSimilar = round(float(kmerRank[genome]), 7)
        if (percentSimilar in percentsCounts.keys()):
            percentsCounts[percentSimilar] += 1
        else:
            percentsCounts[percentSimilar] = 1
    newpercentsCounts = sorted(list(percentsCounts))
    updatedpercentsCounts = {}
    for percentSimilar1 in newpercentsCounts:
        for percentSimilar2 in percentsCounts:
            if (percentSimilar1 == percentSimilar2):
                updatedpercentsCounts[percentSimilar1] = percentsCounts[percentSimilar2]
    percentsSimilar = []
    countsSimilar = []
    for percentSimilar in updatedpercentsCounts:
        percentsSimilar += [percentSimilar]
        countsSimilar += [percentsCounts[percentSimilar]]
    plt.plot(percentsSimilar, countsSimilar)
    plt.grid()
    #plt.fill_between(x=percentsSimilar, y1=countsSimilar, where=(-1<percentsSimilar)&(percentsSimilar>1), color='b', alpha=0.2)
    plt.xlabel("Percent Similar between CP and Rep")
    plt.ylabel("Count")
    plt.title("The Count of each percent similar between CP and Rep of Cruciviruses")
    plt.show()
    return percentsSimilar, countsSimilar
percentsSimilar, countsSimilar = kmerRankShaded("kmerRankCruci.txt")
print(percentsSimilar)
print(countsSimilar)
