from readFasta import readFasta
import matplotlib.pyplot as plt

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
genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_88", "AlphaS_KT948075", 7)
    


