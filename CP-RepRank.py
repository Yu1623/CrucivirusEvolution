from readFasta import readFasta

def CPRepRank(CPfile, Repfile):
    CPfile = readFasta(CPfile)
    Repfile = readFasta(Repfile)
    CPgenes = []
    Repgenes = []
    for geneID in CPfile:
        gene = CPfile[geneID]
        CPgenes += [gene]
    for geneID in Repfile:
        gene = Repfile[geneID]
        Repgenes += [gene]
    geneRank = {}
    minGenes = {}
    if (len(CPgenes) > len(Repgenes)):
        minGenes = Repgenes
    else:
        minGenes = CPgenes
    for i in range(len(minGenes)):
        CPgene = CPgenes[i]
        Repgene = Repgenes[i]
        mingene = ""
        similarSequence = 0
        if (len(CPgene) > len(Repgene)):
            mingene = Repgene
        else:
            mingene = CPgene
        for i in range(len(mingene)):
            if(CPgene[i] == Repgene[i]):
                similarSequence += 1
        for geneID in CPfile:
            if (CPfile[geneID] == CPgene):
                if (geneID in geneRank.keys()):
                    geneRank[geneID] += similarSequence
                else:
                    geneRank[geneID] = similarSequence
    updatedgeneRank = sorted(geneRank, key = geneRank.get, reverse = True)
    return updatedgeneRank

rankingGene = CPRepRank('879CPs.fasta', '855Reps.fasta')
print(rankingGene)
