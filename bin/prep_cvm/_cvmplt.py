#!/usr/bin/python
# -*- coding: utf-8 -*-
import math, os, re, copy, commands, sys

def plt():
    CompList = ["A","B","C","D"]
    convlist = readConvergence()[1:]
    for list in convlist:
        os.chdir(list[0])
        makePlotdata(CompList)
        os.chdir('..')
    makeViewplt(convlist,2,3)
    os.system("gnuplot view.plt")


def makeViewplt(convlist,x,y):
    lines='''#!/usr/bin/gnuplot
        set terminal x11
        set xlabel ""
        set ylabel ""
        set title ""
        set pointsize 2
        set xzeroaxis
        set xtics autofreq
        set key outside
        plot 0 notitle'''
    for list in convlist:
        lines += ', '
        lines += '\'%s/plot.dat\' u %s:%s title \'%s:%s\''%(list[0],x,y,list[0],list[1])
    lines += "\n\n pause -1"
    writeFile("view.plt",lines)


def makePlotdata(complist):
    temp = valGet("f",2,"cv1.txt")
    comp = compGet("c","cv1.txt")
    ene = valGet("f",5,"cv1.txt")
    if len(temp) != len(comp) or len(temp) != len(ene):
        print "data list is not same !"
    else:
        #print "data list is OK !"
        pass
    compnum = len(comp[0])
    outLines = "#temp\t"
    for i in range(0,compnum):
        outLines += "comp_%s\t"%complist[i]
    outLines += "energy\n"
    for i in range(0,len(temp)):
        outLines += "%s\t"%temp[i]
        for val in comp[i]:
            outLines += "%s\t"%val
        energy = float(ene[i])*1.380662e-23*6.02e23/1000
        outLines += "%s\n"%energy
    writeFile("plot.dat",outLines)

def readConvergence():
    lines = readFile("convergence.txt")
    list = [["#RMAX","error"]]
    for i in range(0,len(lines)):
        if lines[i].split()[1] == list[-1][1]:
            pass
        else:
            list.append(lines[i].split())
    return list

def compGet(key,file):
    list = []
    lines = grepFile(key,file)
    for line in lines:
        list.append(line.split()[1:])
    return list

def valGet(key,pos,file):
    list = []
    lines = grepFile(key,file)
    for line in lines:
        if len(line.split()) < pos:
            print "There"
            os.system("pwd")
        else:
            list.append(line.split()[pos])
    return list

def grepFile(key,file):
    c = commands
    c.g = c.getoutput
    all_line = c.g('grep -n \"%s\" %s'%(key,file))
    lines = all_line.split('\n')
    return lines

def readFile(filename):
	f = open(filename,'r')
	lines = f.readlines()
	return lines

def writeFile(fname,lout):
	out = open(fname,'w')
	for line in lout:
		out.write(line)
	out.close

if __name__ == "__main__":
    plt()
