import pandas as pd
from readFasta import readFasta
from geneStrandPrediction import genomeSenseTable

def readCSV(fileName):
    f = open(fileName, 'r')
    lines = f.readlines()
    names = []
    sequenceNames = []
    lengths = []
    directions = []
    documentNames = []
    types = []
    sources = []
    for i in range(len(lines)):
        line = lines[i]
        if (i == 0):
            continue
        if ("Rep" in line.split(",")[0] or "CP" in line.split(",")[0] or "Capsid" in line.split(",")[0]):
            name = line.split(",")[0]
            sequenceName = line.split(",")[1]
            length = line.split(",")[2]
            direction = line.split(",")[3]
            documentName = line.split(",")[4]
            type = line.split(",")[5]
            source = line.split(",")[6]
            names += [name]
            sequenceNames += [sequenceName]
            lengths += [length]
            directions += [direction]
            documentNames += [documentName]
            types += [type]
            sources += [source]
    dict = {"Name": names, "Sequence Name": sequenceNames, "Length": lengths, "Direction": directions, "Document Name": documentNames, "Type": types, "Source": sources}
    df = pd.DataFrame(dict)
    #df.to_csv("crucis875.csv")
    print(df)
    return dict

dict1 = readCSV("875crucisAnnotations.csv")
dict2 = genomeSenseTable("879CPs.fasta", "855Reps.fasta")

commonSense = 0
senseTotal = 0
for i in range(len(dict1["Name"])):
    if ("CP" in dict1["Name"][i] or "Capsid" in dict1["Name"][i]):
        for l in range(len(dict1["Name"])):
            if ("Rep" in dict1["Name"][l]):
                if (dict1["Sequence Name"][i] == dict1["Sequence Name"][l]):
                    print(i, l)
                    print(dict1["Sequence Name"][i], dict1["Sequence Name"][l])
                    directionCP = dict1["Direction"][i]
                    directionRep = dict1["Direction"][l]
                    if (directionCP == directionRep):
                        genomeSense = "unisense"
                    else:
                        genomeSense = "ambisense"
                    print(dict1["Sequence Name"][l], genomeSense)
                    for m in range(len(dict2["Genome Name"])):
                        if (dict2["Genome Name"][m] == dict1["Sequence Name"][l]):
                            genomeSensePrediction = dict2["Genome"][m]
                            senseTotal += 1
                            print(genomeSense, genomeSensePrediction)
                            if (genomeSense == genomeSensePrediction):
                                commonSense += 1
accuracy = commonSense/senseTotal
print(accuracy)