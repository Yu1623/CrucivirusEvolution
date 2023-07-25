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

'''
x = [6.34e-05, 6.15e-05, 6.23e-05, 6.89e-05, 6.23e-05, 6.47e-05, 6.34e-05, 6.82e-05, 6.79e-05, 6.21e-05, 6.38e-05, 6.17e-05, 6.54e-05, 6.59e-05]
y = [1, 2, 4, 1, 2, 1, 1, 1, 1, 2, 4, 2, 1, 1]
plt.plot(x, y)
plt.show()
'''    