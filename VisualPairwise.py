'''
Purpose: Prompts for CP and Rep files and the kmer length, or step, and make heatmaps and visuals of the pairwise comparison percentages among the CP genes and Rep genes separately.
'''

from readFasta import readFasta
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy

'''
This program gets the percents of similar kmers with length 'step' between every possible pair of genes in the file, the labels, or the names of the genes, and the gene sequences.
'''

def pairwise(fileName, step):
    file = readFasta(fileName)
    genes = []
    for gene in file:
        genes += [file[gene]]
    percents = [[column for column in range(len(genes))] for row in range(len(genes))]
    labels = [[column for column in range(len(genes))] for row in range(len(genes))]
    for i in range(len(percents)):
        for l in percents[i]:
            percents[i][l] = 0
    for i in range(len(labels)):
        for l in labels[i]:
            labels[i][l] = ""
    for i in range(len(genes)):
        for l in range(len(genes)):
            if (l <= i):
                print(i, l)
                gene1 = genes[i]
                gene2 = genes[l]
                geneName1 = ""
                geneName2 = ""
                for geneName in file:
                    if (file[geneName] == gene1):
                        geneName1 = geneName
                    if (file[geneName] == gene2):
                        geneName2 = geneName

                print(labels[i][l])
                labels[i][l] = "%s, %s" % (geneName1, geneName2)
                kmers1 = []
                kmers2 = []
                similar = 0
                totalkmer = 0
                numkmers1 = int(len(gene1) - (step - 1))
                for kmer in range(numkmers1):
                    kmer1 = gene1[kmer : kmer + step]
                    kmers1 += [kmer1]
                numkmers2 = int(len(gene2) - (step - 1))
                for kmer in range(numkmers2):
                    kmer2 = gene2[kmer : kmer + step]
                    kmers2 += [kmer2]
                for kmer1 in kmers1:
                    for kmer2 in kmers2:
                        totalkmer += 1
                        if (kmer1 == kmer2):
                            similar += 1
                percent = similar / totalkmer
                percents[i][l] = percent
    return percents,labels, genes

'''
This program creates heatmaps and dendrograms that shows the pairwise comparison among the CP genes and among the Rep genes. 
Please note that the program is very slow in running through the genes for large datasets, such as datasets with over a hundred genes.
'''

def pairwiseVisual(CPfileName, RepfileName, step):
    percents, labels, genes = pairwise("22CPs.fasta", 7)
    df1 = pd.DataFrame(percents)
    labelsy = []
    labelsx = []
    for label in labels[:][-1]:
        print(label)
        labelsy += [label.split(", ")[1].split("_-_")[0]]
    print("\n")
    for label in labels[-1][:]:
        print(label)
        labelsx += [label.split(", ")[1].split("_-_")[0]]
    sns.heatmap(df1, vmin = 6.22e-05, vmax = 0.0008637, cmap = "terrain", xticklabels=labelsx, yticklabels=labelsy)
    plt.title("Heat map of the Pairwise Comparison of CP Genes")
    plt.show()
    
    percents, labels, genes = pairwise("22Reps.fasta", 7)
    df2 = pd.DataFrame(percents)
    labelsy = []
    labelsx = []
    for label in labels[:][-1]:
        labelsy += [label.split(", ")[1].split("_-_")[0]]
    for label in labels[-1][:]:
        labelsx += [label.split(", ")[1].split("_-_")[0]]
    sns.heatmap(df2, vmin=6e-05, vmax = 0.000864, cmap = 'terrain', xticklabels=labelsx, yticklabels=labelsy)
    plt.title("Heat map of the Pairwise Comparison of Rep Genes")
    plt.show()

    temp = hierarchy.linkage(df1, method = 'ward')
    dn1 = hierarchy.dendrogram(temp, above_threshold_color='green', color_threshold=0.7, leaf_rotation = 45, labels=labelsx)
    plt.title("Capsid-Protein Genes Pairwise Comparison")
    plt.show()
    
    temp = hierarchy.linkage(df2, method = 'ward')
    hierarchy.set_link_color_palette(['#174b20', ' #82a19d '])
    dn2 = hierarchy.dendrogram(temp, above_threshold_color='blue', color_threshold=0.7, leaf_rotation = 45, labels = labelsx)
    plt.title("Replication-Protein Genes Pairwise Comparison")
    plt.show()
