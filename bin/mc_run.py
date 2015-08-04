#!/usr/bin/env python
# -*- coding: utf-8 -*-
from mc import *

__date__ = "Nov 17 2014"

def main():
    montecarlo(300)

def montecarlo(T):
    cell = FaceCenterCubic(30, arrange='random', conc=0.98)
    ecis = np.array([-112.688724565625, 92.0961460062501, 119.291397787500,
                -0.909217246875000, -72.9952017250000, -24.1141423375000,
                0, -10.0314454625000 ,1.01281285312500, 2.02562570625000,
                0.337604284375000]) * 1/1000.
    ecis[1] = -62.566 * 1/1000
    moncal = MonteCarlo(ecis, cell, T)
    #de = moncal.delta_e()
    #print(de)
    #print(moncal.transition_probability(de))
    fname = "T" + str(T) + "_init"
    cell.make_poscar(fname)
    for i in range(30):
        conc_ene = moncal._iterationTO_reserved_atoms(10)
        fname = "T" + str(T) + "_" + str(i)
        cell.make_poscar(fname)

    lines = ""
    for l in conc_ene:
        lines += "{0[0]}\t{0[1]}\n".format([str(x) for x in l])
    fname = "T" + str(T) + ".txt"
    with open(fname, 'w') as wfile:
         wfile.write(lines)

    print(cell.get_orderparam(1))
    print(cell.get_orderparam(2))
    print(cell.get_orderparam(3))
    print(cell.get_orderparam(4))
    print(cell.get_orderparam(6))
    print(cell.get_orderparam(11))
    fname = "T" + str(T)
    moncal.cell.save_cell(fname)
    # moncal.iteration(10)
    # print(cell.get_orderparam(1))
    # print(cell.get_orderparam(2))
    # print(cell.get_orderparam(3))
    # print(cell.get_orderparam(4))
    # print(cell.get_orderparam(6))
    # print(cell.get_orderparam(12))
    # moncal.iteration(10)
    # print(cell.get_orderparam(1))
    # print(cell.get_orderparam(2))
    # print(cell.get_orderparam(3))
    # print(cell.get_orderparam(4))
    # print(cell.get_orderparam(6))
    # print(cell.get_orderparam(12))

    #de = moncal.delta_e()
    #print(moncal.transition_probability(de))

if __name__ == '__main__':
    main()
