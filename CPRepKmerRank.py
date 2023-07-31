from readFasta import readFasta

def main():
    writeRank("879CPs.fasta", "855Reps.fasta", 7, "Cruci")
    writeRank("1526RNA_CPs.fasta", "1514RdRPs.fasta", 7, "RNA")
    writeRank("270CRESS_CPs.fasta", "305CRESS_Reps.fasta", 7, "DNA")

def CPRepKmerComparisonRank(CPfileName, RepfileName, step):
    CPfile = readFasta(CPfileName)
    Repfile = readFasta(RepfileName)
    kmerRank = {}
    for i in range(len(CPfile)):
        CPgeneName = list(CPfile)[i]
        CPgenename = CPgeneName.split("_-_")[0]
        print(CPgenename)
        for l in range(len(Repfile)):
            RepgeneName = list(Repfile)[l]
            Repgenename = RepgeneName.split("_-_")[0]
            print(Repgenename)
            if (CPgenename == Repgenename):
                geneName = CPgeneName
                print(geneName)
                CPgene = CPfile[CPgeneName]
                Repgene = Repfile[RepgeneName]
                totalkmers = 0
                similarity = 0
                numKmer1 = int(len(CPgene) - (step - 1))
                numKmer2 = int(len(Repgene) - (step - 1))
                RepgeneKmer = []
                CPgeneKmer = []
                for kmer in range(numKmer1):
                    CPgeneKmer += [CPgene[kmer: kmer + step]]
                for kmer in range(numKmer2):
                    RepgeneKmer += [Repgene[kmer: kmer + step]]
                for kmer1 in RepgeneKmer:
                    for kmer2 in CPgeneKmer:
                        totalkmers += 1
                        if (kmer1 == kmer2):
                            similarity += 1
                percentCommon = similarity / totalkmers
                kmerRank[geneName] = percentCommon
    updatedKmerRank = sorted(kmerRank, key = kmerRank.get, reverse = True)
    newKmerRank = {}
    for gene in updatedKmerRank:
        newKmerRank[gene] = kmerRank[gene]
    return newKmerRank


def writeRank(CPfileName, RepfileName, step, name):
    kmerRank = CPRepKmerComparisonRank(CPfileName, RepfileName, step)
    file = open("kmerRank%s.txt" % (name), 'w')
    for geneName in kmerRank:
        file.write("%s: %s\n" % (geneName.split('_-_')[0], kmerRank[geneName]))
    file.close()

main()
