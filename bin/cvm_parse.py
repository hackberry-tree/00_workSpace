#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
cvm 関連の perser
"""
from parse_cvm import CVMLogEnth

log = CVMLogEnth.from_file('log.txt')
log.print_all_phase_energy()
