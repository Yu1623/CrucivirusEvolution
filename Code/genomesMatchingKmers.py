'''
Purpose: create a line plot of the shared kmer sequences and their count between two genomes
'''

from readFasta import readFasta
import matplotlib.pyplot as plt

def main():
    genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_88", "AlphaS_KT948075", 7)
    genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_88", "AlphaS_KF471057", 7)
    genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruCGE_296", "Bacil_MH617605", 7)
    genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_87", "Bacil_AB193315", 7)


'''
Input: names of two files and names of two genes, length of kmer
Output: arrays of kmer sequences
'''

def kmerSequence(fileName1, fileName2, genomeName1, genomeName2, step):
    file1 = readFasta(fileName1)
    file2 = readFasta(fileName2)
    genome1 = file1[genomeName1]
    genome2 = file2[genomeName2]
    numKmers1 = int(len(genome1) - (step - 1))
    kmers1 = []
    for kmer in range(numKmers1):
        kmers1 += [genome1[kmer : kmer + step]]
    numKmers2 = int(len(genome2) - (step - 1))
    kmers2 = []
    for kmer in range(numKmers2):
        kmers2 += [genome2[kmer : kmer + step]]
    return kmers1, kmers2

'''
Input: names of the two files, names of the two genomes, length of kmer
Output: associative array with the shared kmer sequences and the count
'''

def genomeskmerMatch(fileName1, fileName2, genomeName1, genomeName2, step):
    kmers1, kmers2 = kmerSequence(fileName1, fileName2, genomeName1, genomeName2, step)
    kmerCounts = {}
    for kmer1 in kmers1:
        if (kmer1 not in kmerCounts.keys()):
            for kmer12 in kmers1:
                if (kmer1 == kmer12):
                    if (kmer1 in kmerCounts.keys()):
                        kmerCounts[kmer1] += 1
                    else:
                        kmerCounts[kmer1] = 1
            for kmer2 in kmers2:
                if (kmer1 == kmer2):
                    if (kmer1 in kmerCounts.keys()):
                        kmerCounts[kmer1] += 1
                    else:
                        kmerCounts[kmer1] = 1
    return kmerCounts
#print(genomeskmerMatch("885crucis.fasta", "Cruci_CruV_244", "Cruci_CruV_336", 7)) 

'''
Input: names of the two files, names of the two genomes, length of kmer
Output: line graph with the shared kmer sequences and the count
'''

def genomeskmerMatchVisual(fileName1, fileName2, genomeName1, genomeName2, step):
    kmerCounts = genomeskmerMatch(fileName1, fileName2, genomeName1, genomeName2, step)
    sequences = []
    counts = []
    for kmer in kmerCounts:
        sequences += [kmer]
        counts += [kmerCounts[kmer]]
    print(sequences, counts)
    plt.bar(sequences, counts, width=1.0)
    totalCount = 0
    for count in counts:
        totalCount += count
    averageCount = totalCount / len(sequences)
    plt.axhline(y=averageCount, color = 'r', linestyle = '-')
    plt.xlabel("Kmer Sequence")
    plt.ylabel("Counts")
    plt.title("Comparison between %s and %s" % (genomeName1, genomeName2))
    plt.xticks(rotation=90, fontsize=0.02)
    #plt.ylim(0, 8)
    plt.show()
    print(kmerCounts)


main()
