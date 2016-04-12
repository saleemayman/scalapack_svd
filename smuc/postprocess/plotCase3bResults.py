import glob, os, sys
import numpy as np
from random import*
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors


def main(argv):
    logFileName = sys.argv[1]
    plotTitle = sys.argv[2]
    print logFileName

    timeSVD = np.genfromtxt(logFileName, delimiter = ',')

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.plot(timeSVD[:, 1], timeSVD[:, 3], color='r', marker='o')
    ax1.set_xlim([0, max(timeSVD[:, 1])])
    ax1.set_ylim([min(timeSVD[:, 3]) - 20, max(timeSVD[:, 3]) + 10])
    ax1.set_xlabel('Cores per Executor')
    ax1.set_ylabel('Time to compute SVD [sec]')
    ax1.set_title('Execution times')

    ax1.annotate('Workers=1\nExecutors=1\nPartitions = same as number \nof cores per executor',
                 xy=(-210, -100), xycoords='axes pixels',
                 bbox=dict(boxstyle='square', fc='yellow', alpha=0.3))

    # speed-up graph
    ax2 = fig.add_subplot(122)
    serialTime = timeSVD[0, 3]
    parallelTimes = timeSVD[:, 3]
    speedUp = serialTime/parallelTimes
    ax2.plot(timeSVD[:, 1], speedUp, color='g', marker='o')
    ax2.set_xlim([0, max(timeSVD[:, 1])])
    ax2.set_ylim([0, max(speedUp) + 1])
    ax2.set_xlabel('Cores per Executor')
    ax2.set_ylabel('Speed-up')
    ax2.set_title(plotTitle)
    ax2.annotate('Workers=1\nExecutors=1\nPartitions = same as number \nof cores per executor',
                 xy=(-210, 10), xycoords='axes pixels',
                 bbox=dict(boxstyle='square', fc='yellow', alpha=0.3))

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1])