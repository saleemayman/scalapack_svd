import glob, os, sys
import numpy as np
from random import*
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colors

def get_rand_color(val):
    h,s,v = random()*6, 0.5, 243.2
    colors = []
    for i in range(val):
        h += 3.75#3.708
        tmp = ((v, v-v*s*abs(1-h%2), v-v*s)*3)[5**int(h)/3%3::int(h)%2+1][:3]
        colors.append('#' + '%02x' *3%tmp)
        if i%5/4:
            s += 0.1
            v -= 51.2
    return colors

def main(argv):
    logFileName = sys.argv[1]
    plotTitle = sys.argv[2]
    print logFileName
    
    # get all csv file from current directory
    resultsCSV = sorted(glob.glob(logFileName + '*_results.csv'))

    # read the individual csv files for logFileName
    nFiles = len(resultsCSV)
    colors = ['r', 'g', 'b', 'k', 'c', 'y']

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # plot the time-line data
    for i in range(0, nFiles):
        timeSVD = np.genfromtxt(resultsCSV[i], delimiter = ', ')
        serialTime = timeSVD[0, 3]
        parallelTimes = timeSVD[:, 3]
        speedUp = serialTime/parallelTimes
        ax1.plot(timeSVD[:, 1], timeSVD[:, 3], color=colors[i], marker='o', label=str(int(timeSVD[0, 2])))
        ax2.plot(timeSVD[:, 1], speedUp, color=colors[i], marker='o', label=str(int(timeSVD[0, 2])))

    ax1.set_xlim([0, 65])
    ax1.set_ylim([0, 210])
    ax1.xaxis.set_ticks(range(0, 65, 8))
    ax1.set_xlabel('MPI Processes (cores)')
    ax1.set_ylabel('Time to compute SVD [sec]')
    # re-arrange the legend (sort both labels and handles by labels)
    handles, labels = ax1.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: int(t[0])))
    ax1.legend(handles, labels, loc='best', title='Blocking factor')
    # ax1.set_title(plotTitle)

    # the speed-up plot
    ax2.set_xlim([0, 65])
    ax2.xaxis.set_ticks(range(0, 65, 8))
    ax2.set_xlabel('MPI Processes (cores)')
    ax2.set_ylabel('Speed-up')
    handles, labels = ax2.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: int(t[0])))
    # ax2.legend(handles, labels, loc='best', title='Blocking factor')
    # ax2.set_title('Speed-up')

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1])
