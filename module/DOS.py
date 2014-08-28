#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from math import *
#import numpy
#from numpy import *
#from numpy import dot
#import csv
import re
#from matplotlib import *
#from pylab import plot, show, ylim, xlim, yticks, xlabel, ylabel
#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
#import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties
#import matplotlib
import os


def readAtomNum(line):
    meta = u"\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*"
    atomNumReader = re.compile(meta)
    m = atomNumReader.match(line)
    atomNum = int(m.group(1))
    return atomNum


def readParameter(line):
    meta = u"\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*"
    paraReader = re.compile(meta)
    m = paraReader.match(line)
    energyNum = int(m.group(3))
    fermiEnergy = float(m.group(4))
    return (energyNum, fermiEnergy)

def readDos(line, labelorbit, fermienergy, dosData):
    dosReader = re.compile("\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)\s*")
    m = dosReader.match(line)
    dosData["Energy"].append(float(m.group(1))-fermiEnergy)
    for i in range(2, 20):
        dosData[labelOrbit[i-2]].append(float(m.group(i)))
    return dosData

def getMaxDos(dosdata, outlabels):
    maximum=[]
    for i in range(0,len(outlabels)):
        maximum.append(max(dosdata[outlabels[i]]))
    return max(maximum)

def writeData(dosdata, outlabels, fileout):
    out = open(fileout, 'w')
    out.write("# ")
    for list in outlabels:
        out.write("%s\t"%(list))
    out.write("\n")
    for i in range(0,energyNum):
        for list in outlabels:
            out.write("%f\t"%(dosdata[list][i]))
        out.write("\n")


labelOrbit = ["s","s*","py","py*","pz","pz*","px","px*","dxy","dxy*","dyz","dyz*","dz2","dz2*","dxz","dxz*","dx2","dx2*"]

f = open("data/DOSCAR_polarized", "r")
line = f.readline()
f.close

atomNum = readAtomNum(line)
for i in range(0,5):
    line = f.readline()
energyNum, fermiEnergy = readParameter(line)

#prepare dosData box
dosData = []
for i in range(0,atomNum):
    dosData.append({})
    dosData[i].update({"Energy":[]})
    for j in range(0, 18):
        dosData[i].update({labelOrbit[j]:[]})

#skip lines
for i in range(0,energyNum):
    line = f.readline()
line = f.readline()
#skip lines

for i in range(0,atomNum):
    for j in range(0,energyNum):
        line = f.readline()
        readDos(line, labelOrbit, fermiEnergy,dosData[i])
    line = f.readline()

sumDos = {}
sumDos.update({"Energy":[]})
for i in range(0, 18):
    sumDos.update({labelOrbit[i-2]:[0]*energyNum})
sumDos["Energy"]=dosData[0]["Energy"]

for i in range(0, atomNum):
    for j in range(0,18):
        for k in range(0,energyNum):
            sumDos[labelOrbit[j]][k] += dosData[i][labelOrbit[j]][k]

def makeDOSdata(plotlistOrbital):
    outLabelsDown = ["Energy"]
    for list in plotlistOrbital:
        outLabelsDown.append("%s*"%(list))
    writeData(sumDos, outLabelsDown, "dos_down.dat")

    outLabelsUp = ["Energy"]
    for list in plotlistOrbital:
        outLabelsUp.append("%s"%(list))
    writeData(sumDos, outLabelsUp, "dos_up.dat")

maximum=getMaxDos(sumDos,labelOrbit)

print "Ef =%f"%(fermiEnergy)


