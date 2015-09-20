#!/home/enoki/.conda/envs/py3k/bin/python
# -*- coding: utf-8 -*-
from montecarlo import *
import sys


__date__ = "Nov 17 2014"

def main():
    montecarlo = {'--Cr': montecarloCr, '--Ti': montecarloTi,
                  '--AlCu': montecarloAlCu}[sys.argv[1]]
    montecarlo(sys.argv[2], sys.argv[3], sys.argv[4])
    #montecarlo(100, 30)
    #montecarlo(200, 30)
    #montecarlo(300, 30)

def montecarloAlCu(T, size, conc):
    conc = float(conc)
    fcc = FCCXtal.from_pickle_ecis("cluster.pickle",
                                   arrange='random', conc=conc, size=int(size))
    mc = MonteCarlo(fcc, T=int(T))

    fname = "POSCAR_T" + str(T) + "_init"
    fcc.make_poscar(fname)

    for i in range(10000):
        # mc = MonteCarlo(fcc, T=300)
        mc.loop_fcc_micro_single(100)
        # ene = mc.loop_fcc_micro(10)
        fname = "POSCAR_T" + str(T) + "_" + str(i)
        fcc.make_poscar(fname)
        ene = fcc.get_energy()
        lines = "{0}\n".format(ene)
        fname = "T" + str(T) + ".txt"
        with open(fname, 'a') as wfile:
            wfile.write(lines)
        fname = "T" + str(T) + ".pickle"
        fcc.save_cell(fname)

def montecarloCr(T, size, conc):
    conc = float(conc)
    nacl = NaClXtal.from_pickle_ecis(
        "cluster_int.pickle", conc_int=conc, conc_sub=1-conc, size=int(size))
    mc = MonteCarlo(nacl, T=int(T))
    # int -569.921196 sub -240.030646
    #nacl.sub_ecis[0] = -600.
    #nacl.int_ecis[0] = -500

    fname = "POSCAR_T" + str(T) + "_init"
    nacl.make_poscar(fname)
    for i in range(300):
        ene = mc.loop_pairflip_nacl(10)
        fname = "POSCAR_T" + str(T) + "_" + str(i)
        nacl.make_poscar(fname)

        lines = ""
        for l in ene:
            lines += "{0[0]}\t{0[1]}\t{0[2]}\n".format([str(x) for x in l])
        fname = "T" + str(T) + ".txt"
        with open(fname, 'a') as wfile:
            wfile.write(lines)

    fname = "T" + str(T) + ".pickle"
    nacl.save_cell(fname)

def montecarloTi(T, size, conc):
    conc = float(conc)
    nacl = NaClXtal.from_pickle_ecis(
        "cluster_int_Ti.pickle", conc_int=conc, conc_sub=1-conc, size=int(size))
    mc = MonteCarlo(nacl, T=int(T))
    # int -569.921196 sub -240.030646
    #nacl.sub_ecis[0] = -600.
    #nacl.int_ecis[0] = -500

    fname = "POSCAR_T" + str(T) + "_init"
    nacl.make_poscar(fname)
    for i in range(300):
        ene = mc.loop_pairflip_nacl(10)
        fname = "POSCAR_T" + str(T) + "_" + str(i)
        nacl.make_poscar(fname)
        lines = ""
        for l in ene:
            lines += "{0[0]}\t{0[1]}\t{0[2]}\n".format([str(x) for x in l])
        fname = "T" + str(T) + ".txt"
        with open(fname, 'a') as wfile:
            wfile.write(lines)

    fname = "T" + str(T) + ".pickle"
    nacl.save_cell(fname)


if __name__ == '__main__':
    main()
