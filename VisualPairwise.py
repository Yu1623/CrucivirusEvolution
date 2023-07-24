from readFasta import readFasta
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster import hierarchy

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
file = readFasta("22CPs.fasta")
step = 7
kmers1 = []
kmers2 = []
totalkmer = 0
difference = 0
gene1 = file["CruV-319_-_Putative_CruV_Capsid_Protein"]
gene2 = file["CruV-319_-_Putative_CruV_Capsid_Protein"]
if (gene1 == gene2):
    print(True)
if (len(gene1) < len(gene2)):
    mingene = gene1
else:
    mingene = gene2
numkmers = int(len(mingene) - (step - 1))
for kmer in range(numkmers):
    kmer1 = gene1[kmer : kmer + step]
    kmer2 = gene2[kmer : kmer + step]
    print(kmer1, kmer2)
    kmers1 += [kmer1]
    kmers2 += [kmer2]
    if (kmer1 == kmer2):
        print(True)
for kmer1 in kmers1:
    for kmer2 in kmers2:
        totalkmer += 1
        if (kmer1 != kmer2):
            difference += 1
percent = difference / totalkmer

print(percent)

'''


percents, labels, genes = pairwise("879CPs.fasta", 7)
print(percents)

df1 = pd.DataFrame(percents)

'''
for y in range(df.shape[0]):
    for x in range(df.shape[1]):
        plt.text(x + 0.5, y + 0.5, labels[y][x])
        print(labels[y][x])
'''

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
plt.show()

percents, labels, genes = pairwise("855Reps.fasta", 7)
df2 = pd.DataFrame(percents)
labelsy = []
labelsx = []
for label in labels[:][-1]:
    labelsy += [label.split(", ")[1].split("_-_")[0]]
for label in labels[-1][:]:
    labelsx += [label.split(", ")[1].split("_-_")[0]]
sns.heatmap(df2, vmin=6e-05, vmax = 0.000864, cmap = 'terrain', xticklabels=labelsx, yticklabels=labelsy)
plt.show()

temp = hierarchy.linkage(df1, method = 'ward')
dn1 = hierarchy.dendrogram(temp, above_threshold_color='lightgreen', color_threshold=0.2, labels=labelsx)
temp = hierarchy.linkage(df2, method = 'ward')
hierarchy.set_link_color_palette(['#174b20', ' #82a19d '])
dn2 = hierarchy.dendrogram(temp, above_threshold_color='lightblue', color_threshold=0.2)

plt.show()










'''
def pairwise(fileName, step):
    file = readFasta(fileName)
    genes = []
    labels = []
    for gene in file:
        genes += [file[gene]]
        labels += [gene.split("_-_")[0]]
    percents = [[column for column in range(len(genes))] for row in range(len(genes))]
    for i in range(len(percents)):
        for l in percents[i]:
            percents[i][l] = 0
    for i in range(len(genes)):
        for l in range(len(genes)):
            if (l <= i):
                print(i, l)
                gene1 = genes[i]
                gene2 = genes[l]
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
'''
file = readFasta("22CPs.fasta")
step = 7
kmers1 = []
kmers2 = []
totalkmer = 0
difference = 0
gene1 = file["CruV-319_-_Putative_CruV_Capsid_Protein"]
gene2 = file["CruV-319_-_Putative_CruV_Capsid_Protein"]
if (gene1 == gene2):
    print(True)
if (len(gene1) < len(gene2)):
    mingene = gene1
else:
    mingene = gene2
numkmers = int(len(mingene) - (step - 1))
for kmer in range(numkmers):
    kmer1 = gene1[kmer : kmer + step]
    kmer2 = gene2[kmer : kmer + step]
    print(kmer1, kmer2)
    kmers1 += [kmer1]
    kmers2 += [kmer2]
    if (kmer1 == kmer2):
        print(True)
for kmer1 in kmers1:
    for kmer2 in kmers2:
        totalkmer += 1
        if (kmer1 != kmer2):
            difference += 1
percent = difference / totalkmer

print(percent)
'''

'''


percents, labels, genes = pairwise("22CPs.fasta", 7)
print(min(percents), max(percents))
print(labels)
df = pd.DataFrame(percents, columns=labels, rows = labels)
sns.heatmap(df, vmin = 6.22e-05, vmax = 0.0008637, cmap = "terrain", xticklabels=labels, yticklabels=labels)


temp = hierarchy.linkage(df, method = 'ward')
plt.figure()
dn = hierarchy.dendrogram(temp, above_threshold_color='lightgreen', color_threshold=0.2, labels = labels)
plt.show()
'''

