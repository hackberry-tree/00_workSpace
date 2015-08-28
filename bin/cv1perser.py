#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cv1parser
"""
from parse_cvm import CVMCv1Parser
import pylab

cv1s = CVMCv1Parser.from_file('cv1.txt')
cv1s.data[0]['data'].set_free_energy_00_from_logtxt('log.txt')
pylab.plot(cv1s.data[0]['data']['comp1'], cv1s.data[0]['data']['g00'])
pylab.show()

pylab.close()
