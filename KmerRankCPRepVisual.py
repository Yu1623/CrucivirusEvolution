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


def kmerRankShadedCompare(Crucifile, RNAfile, DNAfile):
    kmerRankCruci = readKmerRanking(Crucifile)
    kmerRankRNA = readKmerRanking(RNAfile)
    kmerRankDNA = readKmerRanking(DNAfile)
    percentsCountsCruci = {}
    percentsCountsRNA = {}
    percentsCountsDNA = {}
    for kmer in kmerRankCruci:
        percentSimilar = round(float(kmerRankCruci[kmer]), 7)
        if (percentSimilar in percentsCountsCruci.keys()):
            percentsCountsCruci[percentSimilar] += 1
        else:
            percentsCountsCruci[percentSimilar] = 1
    newpercentsCountsCruci = sorted(list(percentsCountsCruci))
    updatedpercentsCountsCruci = {}
    for percentSimilar1 in newpercentsCountsCruci:
        for percentSimilar2 in percentsCountsCruci:
            if (percentSimilar1 == percentSimilar2):
                updatedpercentsCountsCruci[percentSimilar1] = percentsCountsCruci[percentSimilar2]
    percentsSimilarCruci = []
    countsSimilarCruci = []
    for percentSimilar in updatedpercentsCountsCruci:
        percentsSimilarCruci += [percentSimilar]
        countsSimilarCruci += [updatedpercentsCountsCruci[percentSimilar]]

    for kmer in kmerRankRNA:
        percentSimilar = round(float(kmerRankRNA[kmer]), 7)
        if (percentSimilar in kmerRankRNA.keys()):
            percentsCountsRNA[percentSimilar] += 1
        else:
            percentsCountsRNA[percentSimilar] = 1
    newpercentsCountsRNA = sorted(list(percentsCountsRNA))
    updatedpercentsCountsRNA = {}
    for percentSimilar1 in newpercentsCountsRNA:
        for percentSimilar2 in percentsCountsRNA:
            if (percentSimilar1 == percentSimilar2):
                updatedpercentsCountsRNA[percentSimilar1] = percentsCountsRNA[percentSimilar2]
    percentsSimilarRNA = []
    countsSimilarRNA = []
    for percentSimilar in updatedpercentsCountsRNA:
        percentsSimilarRNA += [percentSimilar]
        countsSimilarRNA += [updatedpercentsCountsRNA[percentSimilar]]

    for kmer in kmerRankDNA:
        percentSimilar = round(float(kmerRankDNA[kmer]), 7)
        if (percentSimilar in percentsCountsDNA.keys()):
            percentsCountsDNA[percentSimilar] += 1
        else:
            percentsCountsDNA[percentSimilar] = 1
    newpercentsCountsDNA = sorted(list(percentsCountsDNA))
    updatedpercentsCountsDNA = {}
    for percentSimilar1 in newpercentsCountsDNA:
        for percentSimilar2 in percentsCountsDNA:
            if (percentSimilar1 == percentSimilar2):
                updatedpercentsCountsDNA[percentSimilar1] = percentsCountsDNA[percentSimilar2]
    percentsSimilarDNA = []
    countsSimilarDNA = []
    for percentSimilar in updatedpercentsCountsDNA:
        percentsSimilarDNA += [percentSimilar]
        countsSimilarDNA += [updatedpercentsCountsDNA[percentSimilar]]
        
    plt.plot(percentsSimilarCruci, countsSimilarCruci)
    plt.plot(percentsSimilarRNA, countsSimilarRNA)
    plt.plot(percentsSimilarDNA, countsSimilarDNA)
    plt.xlabel("Percent Similar between CP and Rep")
    plt.ylabel("Count")
    plt.title("The Count of each percent similar between CP and Rep")
    plt.legend(["Cruci", "RNA", "DNA"])
    plt.show()
kmerRankShadedCompare("kmerRankCruci.txt", "kmerRankRNA.txt", "kmerRankDNA.txt")
