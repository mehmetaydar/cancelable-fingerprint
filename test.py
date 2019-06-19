import numpy as np
#from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt1
from matplotlib import pyplot as plt2

data = np.array([
    [1, 2],
    [2, 3],
    [3, 6],
])
data2 = np.array([
    [1, 2],
    [3, 4],
    [6, 9],
])
x, y = data.T
x2, y2 = data2.T

plt1.scatter(x,y,s=1, c='black')
plt1.axis('off')
plt1.savefig('test1.png')
plt1.show()
plt1.close()

plt2.scatter(x2,y2)
plt2.savefig('test2.png')

#plt.show()