import matplotlib.pyplot as plt
import numpy as np

x = [1, 2, 3, 4, 5, 6]
y = [2, 4, 8, 16, 32, 64]

scatter = plt.scatter(x, y)

plt.show()

for i in range(len(x)):
    print(i)