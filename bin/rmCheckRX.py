#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
atat/checkrelaxの結果を受け取って
歪みが大きい構造(0.1以上)を除外する
歪みが大きい構造dirにerrorファイルを置く
"""
import os
from subprocess import Popen, PIPE

def main():
    """main"""
    touch_errors()

def touch_errors():
    """
    構造歪みの大きい(0.1以上)のdirにtouch errorする
    """
    checkrelax = Popen('checkrelax', stdout=PIPE).communicate()[0]
    distortion = [float(x.split()[0])
                  for x in checkrelax.split('\n')[:-1]]
    print(distortion)
    # checkrelaxの値が0.1以上になる行を判定
    threshold = 0
    while distortion[threshold] <= 0.1:
        threshold += 1
    distorted_list = checkrelax.split('\n')[threshold:-1]
    distorted = '\n'.join(distorted_list)
    with open('listDistorted.txt', 'w') as wfile:
        wfile.write(distorted)

    distorted_path = [x.split()[1].split('/')[0]
                      for x in distorted_list]
    for path in distorted_path:
        out_path = os.path.join(path, 'error')
        with open(out_path, 'w') as wfile:
            wfile.write('')

if __name__ == "__main__":
    main()
