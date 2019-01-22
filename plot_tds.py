import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.loadtxt('tdsbr.dat'))
    fig.savefig('plot.pdf')
