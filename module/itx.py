#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make .itx(igor graph) file
"""

class Wave(object):
    """
    Handling wave

    input parameters
    data: list
    name: strings
    label: strings
    """

    def __init__(self, name, data, label=None):
        self.name = name
        self.data = data
        self.label = label

    @property
    def str_data(self):
        """
        self.dataの中身をstrに変換してreturn
        """
        return [str(x) for x in self.data]

    def to_itx(self):
        """
        itxのフォーマットで出力
        """
        lines = "WAVES/D {0}\n".format(self.name)
        lines += "BEGIN\n"
        lines += "\t"
        lines += "\n\t".join(self.str_data) + "\n"
        lines += "END\n"
        lines += "X SetScale/P x 0,1,\"\", {0}; ".format(self.name)
        lines += "SetScale y 0,0,\"\", {0}\n".format(self.name)
        return lines

class PlotCommandsMixin(object):
    """
    waves作成以外のcommand関連
    """
    def append_to_graph(self, xwave, ywave, axis_num=0):
        """
        waveをグラフに追加
        (Display & AppendToGraph コマンド)
        arguments:
        wave: 追加するwave
        axis_num: 追加する軸の番号
        """
        if axis_num == 0:
            size = [str(x) for x in self.size]
            cmd = "Display /W=({0})".format(",".join(size))
        else:
            cmd = "AppendToGraph/L=y{0}/B=x{0}".format(axis_num)
        data = "{0} vs {1}".format(xwave.name, ywave.name)
        self._add({cmd: data})

    def modify_graph(self, *lines):
        """
        ModifyGraph コマンド
        """
        modifies = [{"ModifyGraph": x} for x in lines]
        self._add(*modifies) #pylint: disable=W0142

    def label(self, axis, label):
        """
        Label コマンド
        """
        self._add({"Label": "{0} \"{1}\"".format(axis, label)})

    def set_axis(self, axis, axis_range):
        """
        SetAxis command
        arguments:
            axis: axis_name (str)
            axis_range: range of axis (list of float)
        """
        if type(axis_range) == Wave:
            dw = (max(axis_range.data) - min(axis_range.data)) / 5.
            axis_range = [min(axis_range.data)-dw, max(axis_range.data)+dw]
        self._add({"SetAxis": "{0} {1[0]},{1[1]}".format(axis, axis_range)})

    def text_box(self, text, pos):
        text_box = ("TextBox/C/N=text0/F=0/M/A=LT"
                    "/X={0[0]}/Y={0[1]}".format(pos))
        input_text = r"\\F'Times New Roman'\\Z20{0}".format(text)
        self._add({text_box: "\"{0}\"".format(input_text)})


    def _add(self, *add_cmds):
        """
        self.commandに追加する
        """
        for add_cmd in add_cmds:
            self.command_lines.append(add_cmd)


    def vertical3(self, waves):
        """
        縦三つの図をプロット
        """
        wvx, wv0, wv1, wv2 = waves
        self.append_to_graph(wv0, wvx, 0)
        self.append_to_graph(wv1, wvx, 1)
        self.append_to_graph(wv2, wvx, 2)
        # 反対側のx軸を描く為のダミー
        self.append_to_graph(wvx, wvx, 3)
        self.append_to_graph(wvx, wvx, 4)
        # ダミーのプロットは消去
        self.modify_graph("hideTrace({0})=1".format(wvx.name))
        self.modify_graph("hideTrace({0}#1)=1".format(wvx.name))

        # makerの種類 4: maker & lines
        self.modify_graph("mode=4")
        # makerの形状 19: 丸
        self.modify_graph("marker=19")
        # makerのsize
        self.modify_graph("msize=5")
        # makerのedgeを表示
        self.modify_graph("useMrkStrokeRGB=1")

        # plotのフォント
        self.modify_graph("font=\"Times New Roman\"")
        self.modify_graph("fSize=16")

        # plotのmarginやoffset
        self.modify_graph("lblMargin(left)=0", "lblMargin(bottom)=0")
        self.modify_graph("standoff(left)=0", "standoff(bottom)=0",
                          "standoff(y1)=0", "standoff(x1)=0",
                          "standoff(y2)=0", "standoff(x2)=0",
                          "standoff(x3)=0", "standoff(x4)=0")
        self.modify_graph("axOffset(left)=0", "axOffset(bottom)=0")

        # プロットのカラー
        # green
        cwv0 = "rgb({0})=(16385,65535,41303)".format(wv0.name)
        cwv0_edge = "mrkStrokeRGB({0})=(2,39321,1)".format(wv0.name)
        # cyan
        cwv1 = "rgb({0})=(32768,54615,65535)".format(wv1.name)
        cwv1_edge = "mrkStrokeRGB({0})=(0,0,65535)".format(wv1.name)
        # magenta
        cwv2 = "rgb({0})=(65535,49151,62258)".format(wv2.name)
        cwv2_edge = "mrkStrokeRGB({0})=(65535,0,52428)".format(wv2.name)
        self.modify_graph(cwv0, cwv0_edge, cwv1, cwv1_edge, cwv2, cwv2_edge)

        # 目盛りの種類 2: 内側のみ (デフォルトは外側)
        self.modify_graph("tick(left)=2", "tick(bottom)=2",
                          "tick(y1)=2", "tick(x1)=2",
                          "tick(y2)=2", "tick(x2)=2")
        # 反対の軸をon
        self.modify_graph("mirror(left)=1", "mirror(bottom)=1",
                          "mirror(y1)=1", "mirror(y2)=1")
        # 軸のラベルを消去
        self.modify_graph("noLabel(x1)=1", "noLabel(x2)=1", "noLabel(x3)=1",
                          "noLabel(x4)=1", "noLabel(y3)=1", "noLabel(y4)=1")
        # 軸を消去
        self.modify_graph("axThick(y3)=0", "axThick(y4)=0")

        # ラベルの位置
        self.modify_graph("lblPos(left)=70", "lblPos(bottom)=50",
                          "lblPos(y1)=70", "lblPos(y2)=70",
                          "lblLatPos(left)=0", "lblLatPos(bottom)=0")
        # 軸の描写および範囲
        self.modify_graph("freePos(y1)={0,bottom}",
                          "freePos(x1)={0.345,kwFraction}",
                          "freePos(y2)={0,bottom}",
                          "freePos(x2)={0.69,kwFraction}",
                          "freePos(x3)={0.31,kwFraction};",
                          "freePos(x4)={0.655,kwFraction}")
        self.modify_graph("axisEnab(left)={0,0.31}",
                          "axisEnab(y1)={0.345,0.655}",
                          "axisEnab(y2)={0.69,1}")

        self.label('left', wv0.label)
        self.label('bottom', wvx.label)
        self.label('y1', wv1.label)
        self.label('y2', wv2.label)

        self.set_axis('left', wv0)
        dx = (max(wvx.data) - min(wvx.data)) / 10.
        xrange = [min(wvx.data)-dx, max(wvx.data)+dx]
        self.set_axis('bottom', xrange)
        self.set_axis('y1', wv1)
        self.set_axis('y2', wv2)
        self.set_axis('x1', xrange)
        self.set_axis('x2', xrange)
        self.set_axis('x3', xrange)
        self.set_axis('x4', xrange)

        self.text_box(self.title, [5, 2.5])

class Produce(PlotCommandsMixin):
    """
    itxファイルを出力
    arguments:
        title: プロットの名称
        waves: list of waves objects
    """
    def __init__(self, title, waves):
        self.size = [400, 100, 850, 720]
        self.title = title
        self.waves = waves
        self.command_lines = []

    def to_itx(self):
        lines = "IGOR\n"
        lines += "X NewDataFolder/S/O {0}\n\n".format(self.title)
        for wave in self.waves:
            lines += wave.to_itx() + "\n"
        lines += "X Preferences 0\n"
        for command in self.command_lines:
            cmd = list(command.items())[0]
            lines += "X\t{0[0]} {0[1]};DelayUpdate\n".format(cmd)
        return lines


