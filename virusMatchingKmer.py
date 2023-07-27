from readFasta import readFasta
import seaborn as sns
import matplotlib.pyplot as plt

def virusKmerMatch(DNAfileName, CrucifileName, step):
    DNAfile = readFasta(DNAfileName)
    Crucifile = readFasta(CrucifileName)
    DNAs = []
    DNAnames = []
    newDNAs = []
    newDNAnames = []
    Crucis = []
    Crucinames = []
    newCrucis = []
    newCrucinames = []
    for genome in DNAfile:
        DNAs += [DNAfile[genome]]
        DNAnames += [genome]
    for genome in Crucifile:
        Crucis += [Crucifile[genome]]
        Crucinames += [genome]
    
    for i in range(len(DNAnames)):
        if (i <= 21):
            newDNAs += [DNAs[i]]
            newDNAnames += [DNAnames[i]]
    for i in range(len(Crucinames)):
        if (i <= 21):
            newCrucis += [Crucis[i]]
            newCrucinames += [Crucinames[i]]
    DNAs = newDNAs
    DNAnames = newDNAnames
    Crucis = newCrucis
    Crucinames = newCrucinames
        
    percentsMatch = [[column for column in range(22)] for row in range(22)]
    for i in range(len(DNAs)):
        for l in range(len(Crucis)):
            if ((i <= 21) and (l <= 21)):
                kmersDNA = []
                kmersCruci = []
                matchingSequence = 0
                totalSequence = 0
                numkmersDNA = int(len(DNAs[i]) - (step - 1))
                for kmer in range(numkmersDNA):
                    kmersDNA += [DNAs[i][kmer : kmer + step]]
                numkmersCruci = int(len(Crucis[l]) - (step - 1))
                for kmer in range(numkmersCruci):
                    kmersCruci += [Crucis[l][kmer : kmer + step]]

                for kmerDNA in kmersDNA:
                    for kmerCruci in kmersCruci:
                        if (kmerDNA == kmerCruci):
                            matchingSequence += 1
                        totalSequence += 1
                print(matchingSequence, totalSequence)
                percentMatch = matchingSequence / totalSequence
                percentsMatch[i][l] = percentMatch
                print(i, l, percentMatch)
    return percentsMatch, DNAnames, Crucinames

def virusKmerVisual(DNAfileName, CrucifileName, step):
    percentsMatched, DNAnames, Crucinames = virusKmerMatch(DNAfileName, CrucifileName, step)
    sns.heatmap(data = percentsMatched, xticklabels = Crucinames, yticklabels = DNAnames)
    plt.xlabel("Cruci")
    plt.ylabel("DNA")
    plt.title("Comparison of kmer sequences between DNA and Cruci genomes")
    plt.show()

virusKmerVisual("316CRESS.fasta", "885crucis.fasta", 7)

    


