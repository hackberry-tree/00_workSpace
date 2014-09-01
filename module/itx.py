#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make .itx(igor graph) file
"""

class Wave(object):
    """
    Handling wave

    input parameters
    data : list
    label : strings
    """

    def __init__(self, label, data):
        self.label = label
        self.data = data

    @property
    def str_data(self):
        """
        self.dataの中身をstrに変換してreturn
        """
        return [str(x) for x in self.data]

    def __str__(self):
        lines = "WAVES/D {0}\n".format(self.label)
        lines += "BEGIN\n"
        lines += "\t"
        lines += "\n\t".join(self.str_data) + "\n"
        lines += "END\n"
        lines += "X SetScale/P x 0,1,\"\", {0}; ".format(self.label)
        lines += "SetScale y 0,0,\"\", {0}\n".format(self.label)
        return lines

class Preferences(object):
    """
    Preferencesを作成
    """
    def __init__(self):
        self.lines = "X Preferences 0\n"

    def _add(self, string):
        """
        行を追加
        共通部分を省略するためのmethod
        """
        self.lines += "X\t{0};DelayUpdate\n".format(string)

    def vertical3(self, waves):
        """
        縦三つの図をプロット
        """
        wvx, wv1, wv2, wv3 = waves
        self._add("Display /W=(400,88,857,720) "
                  "{0} vs {1}".format(wvx.label, wv1.label))
        self._add("AppendToGraph/L=y2/B=x2 "
                  "{0} vs {1}".format(wvx.label, wv2.label))
        self._add("AppendToGraph/L=y3/B=x3 "
                  "{0} vs {2}".format(wvx.label, wv3.label))
        self._add("AppendToGraph/L=y4/B=x4 "
                  "{0} vs {0}".format(wvx.label))
        self._add("AppendToGraph/L=y5/B=x5 "
                  "{0} vs {0}".format(wvx.label))

        # 不明
