#!/usr/bin/python
# -*- coding: utf-8 -*-
import os, re, copy, commands, sys
OrigDir = commands.getoutput('which sequence_cvm.py | sed "s/sequence_cvm.py//g"')
sys.path.append(OrigDir)
import _cvmplt

#######################
Rmin = 5.0
Rmax = 10.0
Step = 0.2
DefFile = "def.txt"
DefLine = " " + commands.getoutput("""grep TO def.txt| grep RMAX | sed "s/.*RMAX/RMAX/g" | sed "s/\s.*//g" """)
DefCng = " RMAX=%s"
#######################

def main():
    global Rmin, Rmax, Step
    if len(sys.argv) == 1:
        rmin = raw_input("default R minimum 5.0\ninput : ")
        if rmin != '':
            Rmin = float(rmin)
        rmax = raw_input("default R maximum 10.0\ninput : ")
        if rmax != '':
            Rmax = float(rmax)
        step = raw_input("default step 0.2\ninput : ")
        if step != '':
            Step = float(step)
        print "Rmin:%s, Rmax:%s, step:%s"%(Rmin, Rmax, Step)
        _ = raw_input()
        run(Rmin, Step, Rmax)

    else:
        if sys.argv[1] == '-plot' or sys.argv[1] == '--plot':
            _cvmplt.plt()
        elif sys.argv[1] == '-def' or sys.argv[1] == '--def':
            run(Rmin, Step, Rmax)
        else:
            print "sequence_cvm.py has no " + sys.argv[1]

def run(rmin, step, rmax):
    dir = initialize()
    os.chdir(dir)
    cycleCVM(rmin, step, rmax)


def initialize():
    i = commands.getoutput("""ls -F | grep / | grep seq | wc -l""")
    i = int(i) + 1
    dir = 'seq%s'%i
    os.system('mkdir %s'%dir)
    os.system('mkdir -p %s/originals'%dir)
    os.system('cp def.txt IN.CVM energies.txt *.str %s/originals/'%dir)
    return dir

def cycleCVM(ini,stp,fin):
    c = commands
    c.g = c.getoutput
    pwd = c.g('pwd')
    dir = ini
    judge = 'OK'
    while judge == 'OK':
        os.system('mkdir -p %s'%dir)
        os.system('cp originals/* %s'%dir)
        os.system('cat originals/%s | sed -e "s/%s/%s/g" > %s/%s'%(DefFile,DefLine,DefCng%dir,dir,DefFile))
        os.chdir('%s'%dir)
        os.system('sh autorun_cvm.sh')
        err = c.g('grep "RELATIVE PREDICTIVE ERROR:" log.txt | sed -e "s/RELATIVE PREDICTIVE ERROR: //g"')
        os.chdir(pwd)
        if err == '':
            judge = 'NG'
        else:
            os.system('echo %s %s >> convergence.txt'%(dir,err))
        dir += stp
        if dir > fin:
            judge = 'NG'

def readFile(fname):
    f = open(fname,'r')
    lines = f.readlines()
    f.close
    return lines

def writeFile(fname,lout):
    out = open(fname,'w')
    for line in lout:
        out.write(line)
    out.close

if __name__ == "__main__":
    main()
