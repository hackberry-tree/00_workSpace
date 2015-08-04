#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
atat/checkrelaxの結果を受け取って
歪みが大きい構造(0.1以上)を除外する
歪みが大きい構造dirにerrorファイルを置く
"""
import os
import sys
from subprocess import Popen, PIPE

__date__ = "Sep 1 2014"

def main():
    """main"""
    if len(sys.argv) == 1:
        touch_errors()
    else:
        if sys.argv[1] == "--fitout":
            revise_fitout()
        else:
            touch_errors(float(sys.argv[1]))

def get_distorted(thresh=0.1):
    """
    構造歪みの大きい dir を return する
    """
    checkrelax = Popen('checkrelax',
                       stdout=PIPE).communicate()[0].split('\n')[:-1]

    vals = [x for x in checkrelax if x.split()[0] != 'nan']
    nons = [x for x in checkrelax if x.split()[0] == 'nan']
    checkrelax = vals + nons

    distortion = [float(x.split()[0]) for x in checkrelax]

    # checkrelaxの値が0.1以上になる行を判定
    threshold = 0
    while distortion[threshold] <= thresh:
        threshold += 1
    distorted_list = checkrelax[threshold:]
    distorted = '\n'.join(distorted_list)
    with open('listDistorted.txt', 'w') as wfile:
        wfile.write(distorted)

    print('{0} of {1} structures have too large distoion.'
          .format(len(distortion), len(distortion[threshold:])))
    print('{0} structures are O.K.'.format(threshold))
    distorted_path = [x.split()[1].split('/')[0]
                      for x in distorted_list]
    return distorted_path

def touch_errors(thresh=0.1):
    """
    構造歪みの大きい(0.1以上)のdirにtouch errorする
    """
    distorted_path = get_distorted(thresh)
    for path in distorted_path:
        out_path = os.path.join(path, 'error')
        with open(out_path, 'w') as wfile:
            wfile.write('')

def revise_fitout(thresh=0.1):
    """
    fit.out を修正する
    """
    distorted_path = get_distorted(thresh)
    with open("fit.out", 'r') as rfile:
        lines = rfile.readlines()
    out_line = ""
    for line in lines:
        if line.split()[5] in distorted_path:
            # out_line += "#" + line
            pass
        else:
            out_line += line
    with open("fit_wo_distorted.out", 'w') as wfile:
        wfile.write(out_line)


if __name__ == "__main__":
    main()
