'''
Makes a txt.file containing the rank of genomes in a file from the genome with the highest percent of similar CP and Rep kmer sequences to the genome with the lowest percent of similarity.
'''

from readFasta import readFasta

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

