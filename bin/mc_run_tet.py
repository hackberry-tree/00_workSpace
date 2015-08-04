#!/usr/bin/env python
# -*- coding: utf-8 -*-
from mc_fct import *

__date__ = "Nov 17 2014"

def main():
    montecarlo(30)

def montecarlo(T):
    cell = FaceCenterTetragonal(30, arrange='random', conc=0.99)
    ecis = np.array(
        [-11.370633, -611.080423, 688.864922, 271.55793, -10.2773419,
         -59.7960943, 0, 195.368747, 209.502264, 108.438473, 182.923418,
         -663.075602, -786.354929, -164.142113, 68.6509158, 648.47156, 0,
         -44.4990322]) * 1/1000.
    ecis[1] = -1300 * 1/1000.
    moncal = MonteCarlo(ecis, cell, T)

    fname = "T" + str(T) + "_init"
    cell.make_poscar(fname)
    for i in range(30):
        conc_ene = moncal._iterationTO_reserved_atoms(10)
        fname = "POSCAR_T" + str(T) + "_" + str(i)
        cell.make_poscar(fname)

    lines = ""
    for l in conc_ene:
        lines += "{0[0]}\t{0[1]}\n".format([str(x) for x in l])
    fname = "T" + str(T) + ".txt"
    with open(fname, 'w') as wfile:
         wfile.write(lines)
    fname = "T" + str(T) + ".pickle"
    cell.save_cell(fname)    

if __name__ == '__main__':
    main()
