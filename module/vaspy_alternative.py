#!/opt/anaconda/bin/python
# -*- coding: utf-8 -*-
"""alt parameter"""
import math
import vaspy

def main():
    """alt parameter"""
    x = 0.07*2*math.pi
    vaspy.MakeInputs.all('.', kp_rx=0.18, kp_soc=0.15)


if __name__ == '__main__':
    main()
