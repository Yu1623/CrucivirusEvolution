'''
Purpose: Prompts for a file that has the rank of percents of similar kmers between CP and Rep and makes graphs to visualize the rank.
The file can be rechieved from CPRepKmerRank.py.
'''

import numpy as np
import matplotlib.pyplot as plt
import re
import sys
import pandas as pd

'''
This program converts the file of kmer rank into an associative array that is used for the programs below.
'''
def readKmerRanking(fileName):
    f = open(fileName, 'r')
    kmerRank = {}
    lines = f.readlines()
    for line in lines:
        genome = line.split(": ")[0]
        percentSimilar = line.split(": ")[1]
        kmerRank[genome] = percentSimilar
    return kmerRank

'''
This program creates an array of the percents of similar kmers between CP and Rep and an array of the count of each percent.
'''

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

'''
This program creates a line graph with percents of similar kmers between CP and Rep against the counts of the percents.
'''

def kmerCPRepVisual(fileName):
    percentsSimilar, countsSimilar = kmerRankShaded(fileName)
    plt.plot(percentsSimilar, countsSimilar)
    plt.xlabel("Percent Similar between CP and Rep")
    plt.ylabel("Count")
    plt.title("The Count of percents similar between CP and Rep")
    plt.show()

'''
This program combines the three line graphs of crucivirus, RNA virus, and DNA virus into one to compare their percents of similarity between CP and Rep.
'''

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
