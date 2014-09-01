#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
custodianを使用して構造緩和を行う
2回の構造緩和計算を実行
"""
from custodian.custodian import Custodian
from custodian.vasp.handlers import UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob
from mycustodian import myVaspErrorHandler

handlers = [myVaspErrorHandler(), UnconvergedErrorHandler()]
vasp_cmd = ['mpirun', '-n', '$nCores', '/opt/vasp5n/vasp.5.2/vasp']
jobs = [VaspJob(vasp_cmd, final=True, backup=False, suffix="",
       auto_npar=False, gzipped=False)]
c = Custodian(handlers, jobs, max_errors=15)
c.run()


