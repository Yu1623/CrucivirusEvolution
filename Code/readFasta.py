import re
import sys

def readFasta(fileName):
    f = open(fileName, 'r')
    lines = f.readlines()
    hre = re.compile('>(\S+)')
    lre = re.compile('^(\S+)$')
    genes = {}

    for line in lines:
        outh = hre.search(line)
        if outh:
            id = outh.group(1)
        else:
            outl = lre.search(line)
            if (id in genes.keys()):
                genes[id] += outl.group(1)
            else:
                genes[id] = outl.group(1)
    return genes
