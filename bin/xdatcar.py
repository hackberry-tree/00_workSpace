#!/home/enoki/.conda/envs/py3k/bin/python
# -*- coding: utf-8 -*-
import os
from vdos import Atoms

def main():
    make_poscars()


def make_poscars(dst="."):
    path = os.path.join(dst, "XDATCAR")
    atoms = Atoms.from_xdatcar(path, start=-100, end=-1)
    atoms.make_poscars(100, dst=dst)

if __name__ == '__main__':
    main()
