import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from pylab import savefig
from kmerAnalysis import countGenome
import pickle
fileName = '22CPs.fasta' #data set that is to be processed
step = 4
countsArray = countGenome(fileName, step)
frequencies = []
counts = []
for count in countsArray:
    counts += [count]
    frequency = countsArray[count]
    frequencies += [frequency]
figure, axis = plt.subplots()
plt.xlabel("%s-mer frequency" % (step))
plt.ylabel("Count")
scatter = axis.scatter(frequencies, counts)

axisSlider = plt.axes([0.2, 0.01, 0.65, 0.02])
slider = Slider(axisSlider, 'sVal', 2, 8, valinit = step, valstep = 1)

def update(val):
    sVal = slider.val
    newCountsArray = countGenome(fileName, sVal)
    frequencies = []
    counts = []
    for count in newCountsArray:
        counts += [count]
        frequency = newCountsArray[count]
        frequencies += [frequency]
    axis.clear()
    axis.set_xlabel("%s-mer frequency" % (sVal))
    axis.scatter(frequencies, counts)
    figure.canvas.draw()

slider.on_changed(update)
#axis.figure.savefig('/home/yuxuan/Summer Internship/ChartShowingCoordinates5-merfrequency.png')
#pickle.dump(scatter, open(r"/home/yuxuan/Summer Internship/file.pickle", 'wb'))
#plt.show()

#ax = pickle.load(open(r'/home/yuxuan/Summer Internship/file.pickle', 'rb'))
#plt.show()