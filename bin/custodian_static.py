#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
staticのDFT計算を行う
"""
from custodian.custodian import Custodian
from custodian.vasp.handlers import UnconvergedErrorHandler
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob
#from mycustodian import myVaspErrorHandler

print("run custodian_static.py")
print("running...")


handlers = [VaspErrorHandler()]
vasp_cmd = ['mpirun', '-n', '$nCores', '/opt/vasp5n/vasp.5.2/vasp']
jobs = [VaspJob(vasp_cmd, final=True, backup=False, suffix="",
                auto_npar=False, gzipped=False)]
c = Custodian(handlers, jobs, max_errors=3)
c.run()

print("finish custodian_static.py")

