'''
Graphs of kmer sequences and the number of times they appear in a genome. 
Compare kmer sequences of genomes
'''

from readFasta import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

def kmerSequence(fileName, geneName, step):
    file = readFasta(fileName)
    gene = file[geneName]
    numKmer = int(len(gene) / step)
    kmers = []
    for kmer in range(numKmer):
        kmerSeq = gene[kmer: kmer + step]
        kmers += [kmerSeq]
    return kmers

def kmerCount(fileName, geneName, step):
    kmers = kmerSequence(fileName, geneName, step)
    kmerSequenceCount = {}
    for kmer in kmers:
        if (kmer in kmerSequenceCount.keys()):
            kmerSequenceCount[kmer] += 1
        else:
            kmerSequenceCount[kmer] = 1
    return kmerSequenceCount        

def countArray(fileName, geneName, step):
    kmerSequenceCount = kmerCount(fileName, geneName, step)
    countsArray = {}
    for kmer in kmerSequenceCount:
        count = kmerSequenceCount[kmer]
        if (count in countsArray.keys()):
            countsArray[count] += 1
        else:
            countsArray[count] = 1
    return countsArray

def countVisual(fileName, geneName, step):
    countsArray = countArray(fileName, geneName, step)
    counts = []
    frequencies = []
    for count in countsArray:
        counts += [count]
        frequency = countsArray[count]
        frequencies += [frequency]
    plt.scatter(frequencies, counts)
    plt.xlabel("%s-mer frequency" % (step))
    plt.ylabel("Count")
    plt.show()

def kmerGenomeCount(fileName, step):
    file = readFasta(fileName)
    kmerSequenceArray = {}
    for gene in file:
        kmers = kmerSequence(fileName, gene, step)
        for kmer in kmers:
            if (kmer in kmerSequenceArray.keys()):
                kmerSequenceArray[kmer] += 1
            else:
                kmerSequenceArray[kmer] = 1
    return kmerSequenceArray


def countGenome(fileName, step):
    kmerSequenceArray = kmerGenomeCount(fileName, step)
    countsArray = {}
    for kmer in kmerSequenceArray:
        count = kmerSequenceArray[kmer]
        if (count in countsArray.keys()):
            countsArray[count] += 1
        else:
            countsArray[count] = 1
    return countsArray

def findKmerWithCount(fileName, step, count):
    kmerArray = kmerGenomeCount(fileName, step)
    kmerSequence = []
    for kmer in kmerArray:
        if (kmerArray[kmer] == count):
            kmerSequence += [kmer]
    return kmerSequence

def countGenomeVisual(fileName, step):
    countsArray = countGenome(fileName, step)
    counts = []
    frequencies = []
    for count in countsArray:
        counts += [count]
        frequency = countsArray[count]
        frequencies += [frequency]
    color = np.random.randint(1, 5)
    colors = len(frequencies) * [color]
    norm = plt.Normalize(1, 4)
    cmap = plt.cm.PiYG
    figure, axis = plt.subplots()
    scatter = axis.scatter(frequencies, counts)
    annotation = axis.annotate(text = '', xy = (0, 0), xytext = (15, 15), textcoords = 'offset points', bbox = {'boxstyle': 'round', 'fc': 'w'}, arrowprops = {'arrowstyle': '->'})
    annotation.set_visible(False)

    def motion_hover(event):
        annotation_visibility = annotation.get_visible()
        if event.inaxes == axis:
            is_contained, annotation_index = scatter.contains(event)
            if is_contained:
                data_point_location = scatter.get_offsets()[annotation_index['ind'][0]]
                annotation.xy = data_point_location
                kmerSequence = findKmerWithCount(fileName, step, data_point_location[1])
                kmerSeqString = ''
                for kmer in kmerSequence:
                    kmerSeqString += kmer + ' '
                coordinate_label = '({0:.2f}, {1:.2f})'.format(data_point_location[0], data_point_location[1])
                sequence_label = kmerSeqString
                annotation.set_text(coordinate_label + '\n' + sequence_label)
                annotation.get_bbox_patch().set_facecolor(cmap(norm(colors[annotation_index['ind'][0]])))
                annotation.set_alpha(0.4)

                annotation.set_visible(True)
                figure.canvas.draw_idle()
        else:
            if annotation_visibility:
                annotation.set_visible(False)
                figure.canvas.draw_idle()
    figure.canvas.mpl_connect('motion_notify_event', motion_hover)
    plt.xlabel("%s-mer frequency" % (step))
    plt.ylabel("Count")
    plt.show()

def countCompareGenomeVisual(fileName1, fileName2, step):
    countsArray1 = countGenome(fileName1, step)
    countsArray2 = countGenome(fileName2, step)
    counts1 = []
    frequencies1 = []
    counts2 = []
    frequencies2 = []
    for count in countsArray1:
        frequency = countsArray1[count]
        counts1 += [count]
        frequencies1 += [frequency]
    for count in countsArray2:
        frequency = countsArray2[count]
        counts2 += [count]
        frequencies2 += [frequency]
    figure, (axis1, axis2) = plt.subplots(nrows = 1, ncols = 2, sharey = True)
    axis1.scatter(frequencies1, counts1)
    axis2.scatter(frequencies2, counts2)
    axis1.set_title(fileName1)
    axis2.set_title(fileName2)
    plt.suptitle('%s-mer frequency' % (step))
    plt.ylabel("Count")
    plt.show()

def compareGenomeCountSameChart(fileName1, fileName2, step):
    countsArray1 = countGenome(fileName1, step)
    countsArray2 = countGenome(fileName2, step)
    counts1 = []
    frequencies1 = []
    counts2 = []
    frequencies2 = []
    for count in countsArray1:
        counts1 += [count]
        frequencies1 += [countsArray1[count]]
    for count in countsArray2:
        counts2 += [count]
        frequencies2 += [countsArray2[count]]
    plt.scatter(frequencies1, counts1)
    plt.scatter(frequencies2, counts2)
    plt.legend(['CP', 'Rep'])
    plt.xlabel("Frequency")
    plt.ylabel("Count")
    plt.show()

def genomeCountDotCompare(fileName1, fileName2, step):
    kmerCount1 = kmerGenomeCount(fileName1, step)
    kmerCount2 = kmerGenomeCount(fileName2, step)
    counts1 = []
    kmers1 = []
    counts2 = []
    kmers2 = []
    for kmer in kmerCount1:
        counts1 += [kmerCount1[kmer]]
        kmers1 += [kmer]
    for kmer in kmerCount2:
        counts2 += [kmerCount2[kmer]]
        kmers2 += [kmer]
    
    figure = go.Figure()
    figure.add_trace(go.Scatter(x = counts1, y = kmers1, name = 'Counts of Kmer Sequences CP', marker = dict(color = 'rgba(156, 165, 196, 0.94)', line_color = 'rgba(156, 165, 196, 1.2)')))
    figure.update_traces(mode = 'markers', marker = dict(line_width = 1, symbol = 'circle', size = 16))
    figure.add_trace(go.Scatter(x = counts2, y = kmers2, name = 'Counts of Kmer Sequence Rep', marker = dict(color = 'rgb(88, 123, 204)', line_color = 'rgba(217, 217, 217, 1.2)')))
    figure.update_traces(mode = 'markers', marker = dict(line_width = 1, symbol = 'circle', size = 10))
    figure.update_layout(title = "The number of times each kmer sequence appear in the gene", xaxis=dict(showgrid = False, showline = True, linecolor = 'rgb(102, 102, 101)', tickfont_color = 'rgb(102, 102, 101)', showticklabels = True, dtick = 'outside', tickcolor = 'rgb(102, 102, 101)'), margin = dict(l = 140, r = 40, b = 50, t = 80), legend = dict(font_size = 10, yanchor = 'middle', xanchor = 'right'), width = 800, height = 600, paper_bgcolor = 'white', plot_bgcolor = "white", hovermode = 'closest')
    figure.show()

def countMaxMin(fileName, step):
    kmerArray = kmerGenomeCount(fileName, step)
    minCount = 10000
    maxCount = -1
    for kmer in kmerArray:
        count = kmerArray[kmer]
        if (count <= minCount):
            minCount = count
        if (count >= maxCount):
            maxCount = count
    minKmerFrequency = []
    maxKmerFrequency = []
    for kmer in kmerArray:
        if (kmerArray[kmer] == minCount):
            minKmerFrequency += [kmer]
        if (kmerArray[kmer] == maxCount):
            maxKmerFrequency += [kmer]
    return maxKmerFrequency, maxCount
