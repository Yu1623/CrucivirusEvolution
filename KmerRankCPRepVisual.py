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
    return percentsSimilar, countsSimilar

def kmerCPRepVisual(fileName):
    percentsSimilar, countsSimilar = kmerRankShaded(fileName)
    plt.plot(percentsSimilar, countsSimilar)
    plt.xlabel("Percent Similar between CP and Rep")
    plt.ylabel("Count")
    plt.title("The Count of percents similar between CP and Rep")
    plt.show()

def kmerCPRepCompareVisual(Crucifile, RNAfile, DNAfile):
    percentsSimilarCruci, countsSimilarCruci = kmerRankShaded(Crucifile)
    percentsSimilarRNA, countsSimilarRNA = kmerRankShaded(RNAfile)
    percentsSimilarDNA, countsSimilarDNA = kmerRankShaded(DNAfile)
    plt.plot(percentsSimilarCruci, countsSimilarCruci, alpha = 0.7)
    plt.plot(percentsSimilarRNA, countsSimilarRNA, alpha = 0.7)
    plt.plot(percentsSimilarDNA, countsSimilarDNA, alpha = 0.7)
    plt.xlabel("Percent Similar between CP and Rep")
    plt.ylabel("Count")
    plt.title("Comparison of Count of Similarity Percentage")
    plt.legend(["Cruci", "RNA", "DNA"])
    plt.show()
