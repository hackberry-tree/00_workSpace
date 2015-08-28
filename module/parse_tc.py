#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""
thermocalc 関連の perther
"""
import re
import numpy as np
import pylab
from fractions import Fraction
from commopy import DataBox
from fitting_analysis import FitData


class ExpFile(DataBox):
    """
    .exp file を取り扱う
    """
    def __init__(self, data):
        DataBox.__init__(self, data)
        self.output_keys = []

    @classmethod
    def from_file(cls, fname):
        with open(fname, 'r') as rfile:
            lines = rfile.readlines()
        meta_head = re.compile('^BLOCK .*')
        meta_tail = re.compile('^BLOCKEND*')
        data = []
        for i in range(len(lines)):
            if meta_head.match(lines[i]):
                i += 2
                while not meta_tail.match(lines[i]):
                    data.append({'x': float(lines[i].split()[0]),
                                 'y': float(lines[i].split()[1])})
                    i += 1
        data.sort(key=lambda x: x['x'])
        return cls(data)

    def get_y_interporated(self, xval):
        if xval < self['x'].min() or xval > self['x'].max():
            print('x range is out')
            exit()
        i = np.where(self['x'] == xval)
        if i[0].any():
            return self['y'][i[0][0]]
        else:
            x = self['x'] - xval
            i = 0
            while x[i] <= 0:
                i += 1
            y0 = self['y'][i-1]
            y1 = self['y'][i]
            x0 = self['x'][i-1]
            x1 = self['x'][i]
            dy = y1 - y0
            dx = x1 - x0
            yval = dy / dx * (xval - x0) + y0
            return yval



