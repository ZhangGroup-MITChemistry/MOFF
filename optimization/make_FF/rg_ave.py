import sys
import os
import math
import numpy


filename=sys.argv[1]

rg=numpy.loadtxt(filename,skiprows=501)

print(numpy.mean(rg[:,1]))
print(numpy.std(rg[:,1]))
