'''
Purpose: find the genomes that may show horizontal genetic transfer through unexpected genome sense
'''

from readFasta import readFasta
import geneStrandPrediction

'''
Input: file name and genome name of the gene that the program finds
Output: the gene that is in the genome
'''
def findGene(fileName, genomeName):
    file = readFasta(fileName)
    gene = ""
    for geneName in file:
        if (genomeName in geneName):
            gene = geneName
    return gene

'''
Identify the genomes that have unexpected genome sense
'''

def genomeSenseSurprising(csvfileName, CPfileName, RepfileName):
    accuracy, incorrectPrediction = geneStrandPrediction.checkPrediction(csvfileName, CPfileName, RepfileName)
    percentsSimilar = 0
    genomeSimilar = []
    for genome in incorrectPrediction:
        geneNameCP = findGene("879CPs.fasta", genome)
        geneNameRep = findGene("855Reps.fasta", genome)
        CPbias, CPpercenta, CPpercentt = geneStrandPrediction.geneStrandPrediction("879CPs.fasta", geneNameCP)
        Repbias, Reppercenta, Reppercentt = geneStrandPrediction.geneStrandPrediction("855Reps.fasta", geneNameRep)
        if ((abs(CPpercenta - CPpercentt) > 0.05) and (abs(Reppercenta - Reppercentt) > 0.05)):
            #print("CP %s percenta: %s, percentt: %s \nRep %s percenta: %s, percentt: %s" % (geneNameCP, CPpercenta, CPpercentt, geneNameRep, Reppercenta, Reppercentt))
            percentsSimilar += 1
            genomeSimilar += [genome]
    return percentsSimilar, genomeSimilar
