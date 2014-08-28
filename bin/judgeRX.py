#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Total Energyに対して10回緩和した場合のエネルギー差が1e-6以下の場合に緩和終了
それ以上は緩和を繰り返す
nCoresをexportしておく必要がある
"""
import math
import shutil
import vaspy
import os
from commopy import Bash

VASP_CMD = 'mpirun -n $nCores /opt/vasp5/vasp.5.2/vasp'
NUM = 10
PATH = '.'


def main():
    """main"""
    relax(PATH, NUM, VASP_CMD)


def judge_oszicar(path):
    """
    エネルギー差を確認して以下の場合 Trueをreturn
    1. NUMより計算回数が小さい
    2. 最初と最後の状態のエネルギー差が小さい
    """
    fname = os.path.join(path, 'OSZICAR')
    oszicar = vaspy.Oszicar(fname)
    ene_ini = oszicar.results[0]['energy']
    ene_fin = oszicar.results[-1]['energy']
    nsw_fin = oszicar.results[-1]['nsw_num']
    incar_path = os.path.join(path, 'INCAR')
    nsw = vaspy.IncarReadWriteMixin.read_incar(incar_path)['nsw']
    if nsw_fin == nsw:
        return False

    delta = 2 * (ene_ini - ene_fin) / (ene_ini + ene_fin)
    delta = math.fabs(delta)
    if delta < 1e-6:
        print("Difference of energy between initial and end is {0}"
              .format(delta))
        print("Enough small. Lattice relaxation is finished.")
        return True
    else:
        print("Difference of energy between initial and end is {0}"
              .format(delta))
        print("Not enough small.")

        return False


def relax(path, num, vaspcmd):
    """
    judgeの結果がTrueの場合終了
    それ以外はOSZICAR, OUTCARを保管してVASPを再実行
    """
    for i in range(0, num):
        judge = judge_oszicar(path)
        if judge:
            break
        else:
            shutil.copyfile('OUTCAR', "out.relax.{0}".format(i))
            shutil.copyfile('OSZICAR', "osz.relax.{0}".format(i))
            shutil.copyfile('CONTCAR', 'POSCAR')
            print("Start re-calculation.")
            Bash.execute(vaspcmd)

if __name__ == '__main__':
    main()
