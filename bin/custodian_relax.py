#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
custodianを使用して構造緩和を行う
2回の構造緩和計算を実行
"""
from custodian.custodian import Custodian
from custodian.vasp.handlers import UnconvergedErrorHandler
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob
#from mycustodian import myVaspErrorHandler

print("start custodian_relax.py")
print("running...")

handlers = [VaspErrorHandler(), UnconvergedErrorHandler()]
vasp_cmd = ['mpirun', '-n', '$nCores', '/opt/vasp5n/vasp.5.2/vasp']
jobs = [VaspJob(vasp_cmd, final=False, suffix=".relax1", auto_npar=False),
        VaspJob(vasp_cmd, final=True, backup=False, suffix="", gzipped=False,
                auto_npar=False,
                settings_override=[{"file": "CONTCAR",
                                    "action": {"_file_copy":
                                               {"dest": "POSCAR"}}}])]
c = Custodian(handlers, jobs, max_errors=15)
c.run()

print("finish custodian_relax.py")
