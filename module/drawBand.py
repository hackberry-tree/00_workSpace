 #!/usr/bin/python
 # -*- coding: utf-8 -*-
#import sys
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
import DOS

#plotlistOrbital = ["s","py","pz","px","dxy","dyz","dz2","dxz","dx2"]
orbitalLists = {"s":[1,"",""],"py":[2,"",""],"pz":[3,"",""],"px":[4,"",""],"dxy":[5,"green","dxy"],"dyz":[6,"salmon","dyz"],"dz2":[7,"purple","d3z^{/=7 2}-r^{/=7 2}"],"dxz":[8,"red","dxz"],"dx2":[9,"blue","dx^{/=7 2}-y^{/=7 2}"]}

def getMeta(line):
    metaMatcher = re.compile("#[^:]+:([^#]+)#[^:]+:([^#]+)#[^:]+:(.+)$")
    m = metaMatcher.match(line)
    kpointsCount = int(m.group(1))
    bandsCount = int(m.group(2))
    ionsCount = int(m.group(3))
    return (kpointsCount, bandsCount, ionsCount)

def readKPoint(stream, bandsCount, ionsCount):
    kpointHeaderMatcher = re.compile("\s*k-point[^:]+:\s*([^\s]+)\s*([^\s]+)\s*([^\s]+)")
    line = stream.readline()
    while (not kpointHeaderMatcher.match(line)):
        line = stream.readline()
    kpointHeaderLine = line
    m = kpointHeaderMatcher.match(kpointHeaderLine)
    kpoint_x = float(m.group(1))
    kpoint_y = float(m.group(2))
    kpoint_z = float(m.group(3))

    bandHeaderMatcher = re.compile(".+energy\s+([^#]+)# occ\.\s(.+)$")
    bands = []
    for i in range(0,bandsCount):
        line = stream.readline()
        while (not bandHeaderMatcher.match(line)):
            line = stream.readline()
        bandHeaderLine = line
        energy = float(bandHeaderMatcher.match(line).group(1))

        orbitalHeaderMatcher = re.compile("ion\s+s.+tot$")
        line = stream.readline()
        while (not orbitalHeaderMatcher.match(line)):
            line = stream.readline()
        orbitalHeaderLine = line
        orbitals = []
#各原子の合計のみ習得
        for j in range(0,ionsCount):
            line = stream.readline()
        line = stream.readline()
        orbitalAssignmentText = line.split()
        orbitals.append(orbitalAssignmentText)
#原子それぞれに分離
#    for j in range(0,ionsCount):
#      line = stream.readline()
#      orbitalAssignmentText = line.split()
#      orbitals.append(orbitalAssignmentText)
        bands.append({ "energy": energy,"orbitals": orbitals })
    return { "kpoint_x": kpoint_x, "kpoint_y": kpoint_y, "kpoint_z": kpoint_z, "bands": bands }

def readFile(procar):
    f = open(procar, "r")
    f.readline() # throwaway
    metaLine = f.readline() # metadata
    kpointsCount, bandsCount, ionsCount = getMeta(metaLine)
    kpoints = []
    for i in range(0,kpointsCount*2):
        kpoint = readKPoint(f, bandsCount, ionsCount)
        kpoints.append(kpoint)
    return (kpoints, bandsCount, ionsCount)


def writef(fname,lout):
    out = open(fname,'w')
    for line in lout:
        out.write(line)
    out.close

def changeGPloader(maxdos, outputfile, gploader):
    global orbitalLists
    global plotlistOrbital
    f = open(gploader, 'r')
    lines = f.readlines()
    f.close
    lines[0] = "maxdos = %f\n"%(maxdos)
    lines[1] = "color1 = \'%s\'\n"%(orbitalLists[plotlistOrbital[0]][1])
    lines[2] = "color2 = \'%s\'\n"%(orbitalLists[plotlistOrbital[1]][1])
    lines[3] = "label1 = \'{/Times-Italic %s}\'\n"%(orbitalLists[plotlistOrbital[0]][2])
    lines[4] = "label2 = \'{/Times-Italic %s}\'\n"%(orbitalLists[plotlistOrbital[1]][2])
    lines[5] = "set output \'%s\'\n"%(outputfile)
    out = open(gploader, 'w')
    for line in lines:
        out.write(line)

kpoints, bandsCount, ionsCount = readFile("data/PROCAR_band")
kpointsCount = len(kpoints)

#print kpoints[685]
#print kpoints[686]


orbitalNames = ["s","py","pz","px","dxy","dyz","dz2","dxz","dx2"]

#print "band,ion,s,py,pz,px,dxy,dyz,dz2,dxz,dx2,tot,energy,kpoint_x,kpoint_y,kpoint_z"


kpointNums = []
Ef = DOS.fermiEnergy


orbEnergies=[]
orbNum=0
for i in range(0, len(orbitalNames)+1):
    orbEnergies.append([])
    for j in range(0, bandsCount):
        orbEnergies[orbNum].append([])
    orbNum += 1

kpointNum = 1
for kpoint in kpoints:
    bandNum = 0
    for band in kpoint["bands"]:
        for ion in band["orbitals"]:
            character = {}
            for i in range(1,10):
                character.update({orbitalNames[i-1]:ion[i]})
        orbNum=0
        orbEnergies[orbNum][bandNum].append(band["energy"]-Ef)
        orbNum += 1
        for orb in orbitalNames:
            orbEnergies[orbNum][bandNum].append(float(character[orb]))
            #else :
            #  orbEnergies[orbNum][bandNum].append([band["energy"]-Ef,0])
            orbNum += 1
        bandNum += 1
        kpointNums.append(kpointNum)
    kpointNum += 1

def deGenerate(orbenergies):
    degenXY = []
    for i in range(0,bandsCount):
        degenXY.append([])
        for j in range(0,kpointsCount):
            degenXY[i].append((orbenergies[orbitalLists["dyz"][0]][i][j]+orbEnergies[orbitalLists["dxz"][0]][i][j])/1.414)
    orbenergies[orbitalLists["dyz"][0]] = degenXY
    orbenergies[orbitalLists["dxz"][0]] = degenXY

deGenerate(orbEnergies)
#print orbEnergies[orbitalLists["dyz"][0]][5][5]

def strageData(plotlistOrbital):
    strage_data_up = []
    strage_data_up.append("#\tenergy\t")
    for list in plotlistOrbital:
        strage_data_up[0] += "%s\t"%(list)
    strage_data_up[0] += "\n"
    for i in range(0,bandsCount):
        for j in range(0,kpointsCount/2):
            Energies = "%f\t"%(orbEnergies[0][i][j])
            for list in plotlistOrbital:
                Energies += "%f\t"%(orbEnergies[orbitalLists[list][0]][i][j])
            strage_data_up.append('%d\t%s\n'%(j,Energies))
        dummy = "-100.00000\t"
        for k in range(0, len(plotlistOrbital)):
            dummy += '%f\t'%(0)
        strage_data_up.append('%d\t%s\n'%(len(orbEnergies[0][i]),dummy))
    writef('band_up.dat',strage_data_up)

    strage_data_down=[]
    strage_data_down.append("#\tenergy\t")
    for list in plotlistOrbital:
        strage_data_down[0] += "%s\t"%(list)
    strage_data_down[0] += "\n"
    for i in range(0,bandsCount):
        for j in range(kpointsCount/2,kpointsCount):
            Energies = "%f\t"%(orbEnergies[0][i][j])
            for list in plotlistOrbital:
                Energies += "%f\t"%(orbEnergies[orbitalLists[list][0]][i][j])
            strage_data_down.append('%d\t%s\n'%(j-kpointsCount/2,Energies))
        dummy = "-100.00000  "
        for k in range(0, len(plotlistOrbital)):
            dummy += '%f\t'%(0)
        strage_data_down.append('%d\t%s\n'%(len(orbEnergies[0][i]),dummy))
    writef('band_down.dat',strage_data_down)

def epsTopdf(plotname):
    os.system("gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s.pdf %s.eps"%(plotname,plotname))
    os.system("rm %s.eps"%(plotname))

def drawBand(plotlistOrbital):
    DOS.makeDOSdata(plotlistOrbital)
    strageData(plotlistOrbital)
    os.system("mkdir %s_%s/"%(plotlistOrbital[0],plotlistOrbital[1]))
    plotNames = {"dual_plot":"00_dual","up_plot":"01_up","down_plot":"02_down"}
    for key in plotNames:
        changeGPloader(DOS.maximum, "%s_%s/%s.eps"%(plotlistOrbital[0],plotlistOrbital[1],plotNames[key]), key)
        os.system("gnuplot %s"%(key))
    for key in plotNames:
        epsTopdf("%s_%s/%s"%(plotlistOrbital[0],plotlistOrbital[1],plotNames[key]))

plotlistOrbital = ["dxy","dx2"]
print plotlistOrbital
drawBand(plotlistOrbital)

plotlistOrbital = ["dxy","dyz"]
print plotlistOrbital
drawBand(plotlistOrbital)

plotlistOrbital = ["dxz","dyz"]
print plotlistOrbital
drawBand(plotlistOrbital)

plotlistOrbital = ["dz2","dxz"]
print plotlistOrbital
drawBand(plotlistOrbital)

plotlistOrbital = ["dx2","dxz"]
print plotlistOrbital
drawBand(plotlistOrbital)


#os.system("open -a Preview *.eps")




#rcParams['lines.linewidth'] = 4.0
#rcParams['font.weight'] = "heavy"
#rcParams['axes.linewidth'] = 4.0
#rcParams['xtick.major.size'] = 0.0
#rcParams['xtick.minor.size'] = 0.0
#rcParams['ytick.major.size'] = 0.0
#rcParams['ytick.minor.size'] = 0.0

#for i in range(0,bandsCount):
#  plot(range(0,len(energies[i])),energies[i], 'bo')
#xlim(0,60)
#ylim(-10,20)
#plt.title("Band Structure")
#xlabel("K-Points")
#ylabel("Energy(eV)")
#plt.show()
