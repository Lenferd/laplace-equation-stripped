import os
import re
import sys
import numpy as np

from matplotlib.ticker import FuncFormatter

pathSerial = "SerialResult.txt"
pathParallel = "ParallelResult.txt"

serial = np.loadtxt(pathSerial, dtype=float)
parallel = np.loadtxt(pathParallel, dtype=float)

serial = np.reshape(serial, len(serial)*len(serial[0]))
parallel = np.reshape(parallel, len(parallel)*len(parallel[0]))

absoluteFault = max([abs(xi-xj) for xi, xj in zip(serial, parallel)])

print("Абсолютная ошибка:\t%.15f" % absoluteFault)
