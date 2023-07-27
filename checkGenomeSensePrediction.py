from readFasta import readFasta
import geneStrandPrediction

def findGene(fileName, genomeName):
    file = readFasta(fileName)
    gene = ""
    for geneName in file:
        if (genomeName in geneName):
            gene = geneName
    return gene

accuracy, incorrectPrediction = geneStrandPrediction.checkPrediction()
percentsSimilar = 0
genomeSimilar = []
for genome in incorrectPrediction:
    geneNameCP = findGene("879CPs.fasta", genome)
    geneNameRep = findGene("855Reps.fasta", genome)
    CPbias, CPpercenta, CPpercentt = geneStrandPrediction.geneStrandPrediction("879CPs.fasta", geneNameCP)
    Repbias, Reppercenta, Reppercentt = geneStrandPrediction.geneStrandPrediction("855Reps.fasta", geneNameRep)
    if ((abs(CPpercenta - CPpercentt) > 0.05) or (abs(Reppercenta - Reppercentt) > 0.05)):
        #print("CP %s percenta: %s, percentt: %s \nRep %s percenta: %s, percentt: %s" % (geneNameCP, CPpercenta, CPpercentt, geneNameRep, Reppercenta, Reppercentt))
        percentsSimilar += 1
        genomeSimilar += [genome]
print(percentsSimilar, genomeSimilar)