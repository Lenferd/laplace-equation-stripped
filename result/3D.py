from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pylab
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import os
import sys


def build_3D_function(ax, title, X, Y, Z, stride):
    ax.plot_surface(X, Y, Z, alpha=0.6, linewidth=0.9, edgecolors='g',
                    rstride=stride, cstride=stride, antialiased=True, shade=False, cmap=cm.coolwarm)

    plt.title(title)


def main():
    if len(sys.argv) < 3:
        print("Not enough arguments")
        exit(-1)
    else:
        file1 = sys.argv[1]
        file2 = sys.argv[2]

    path_setting = os.path.join(os.pardir, "initial", "settings.ini")
    with open(path_setting) as file:
        setting = {line.split('=')[0]: float(line.split('=')[1]) for line in file}

    N = int(setting['DIM'])

    XSTART = setting['XSTART']
    XEND = setting['XEND']
    YSTART = setting['YSTART']
    YEND = setting['YEND']

    X = np.linspace(XSTART, XEND, N, dtype=float)
    Y = np.linspace(YSTART, YEND, N, dtype=float)
    X, Y = np.meshgrid(X, Y)

    fig = pylab.figure('PLOTS')

    expected_Z = np.sin(math.pi*X)*np.exp(-math.pi*Y)

    ZS = np.loadtxt(file1, unpack=True)
    ZP = np.loadtxt(file2, unpack=True)

    stride = N // 20
    print("размер шага равен {:d}".format(stride))

    ax1 = fig.add_subplot(121, projection='3d')
    build_3D_function(ax1, 'INIT', X, Y, ZS, stride)

    ax2 = fig.add_subplot(122, projection='3d')
    build_3D_function(ax2, 'RESULT', X, Y, ZP, stride)

    plt.show()


if __name__ == '__main__':
    main()
