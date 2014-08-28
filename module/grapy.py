#!/usr/bin/python
# -*- coding: utf-8 -*-
"""arrayデータをplotする"""
import solid
import pylab


class PlotStyle(object):
    """
    pylabのstyleを設定
    """
    figsize_dict = [{'figsize': (5, 5/1.618)},
                    {'figsize': (5, 10/1.618)},
                    {'figsize': (5, 5*1.618)}]

    axis_dict = {'energy': r"Energy ($eV$/atom)",
                 'mag': r"Mag. ($\mu_B$/atom)",
                 'c/a': r"$c/a$",
                 'volume': r"Volume ($\AA^3$)",
                 'enthalpy': r"Enthalpy ($kJ/mol$)",
                 'mae': r"MAE ($\mu eV$/atom)",
                 'order': ""}

    def __init__(self, figs_num):
        """
        style: dict形式
        """
        self.style = None
        self.fig = pylab.figure(**self.figsize_dict[figs_num-1])
        self.figs_num = figs_num
        default_font = {'font.family': 'Times New Roman'}
        pylab.rcParams.update(default_font)

    def set_style(self, color, style='default', label=None):
        """
        self.styleにdictをセット
        """
        style_dict = {'default': self.__plotstyle00,
                      'line': self.__plotstyle01}
        self.style = style_dict[style](color)
        if label:
            self.style.update({'label': label})

    def __plotstyle00(self, color):
        """
        デフォルトのプロットスタイル
        丸と線
        外からは呼ばない
        """
        style = self.make_colorstyle(color)
        style.update({'markersize': 6, 'linewidth': 1.25,
                      'markeredgewidth': 1.5,
                      'marker': 'o', 'linestyle': '-'})
        return style

    def __plotstyle01(self, color):
        """
        lineのみ
        """
        style = self.make_colorstyle(color)
        style.update({'markersize': 1, 'marker': 'o', 'linestyle': ''})
        return style

    @staticmethod
    def make_colorstyle(color):
        """
        プロットの色を記述するstyle_dictを作成
        """
        copa = solid.COLOR_PALETTE
        style = {'color': copa[color][0], 'markeredgecolor': copa[color][0],
                 'markerfacecolor': copa[color][1]}
        return style

    @staticmethod
    def set_fontsize(axis, size):
        """
        ラベルのフォントサイズを変更する
        """
        axis.tick_params(axis='both', which='major', labelsize=size)

    def set_title(self, title):
        """タイトルをセットする"""
        pos = {1: {'y': 1.00},
               2: {'y': 0.97},
               3: {'y': 0.95}}
        self.fig.suptitle(title, fontsize='xx-large', **pos[self.figs_num])


class Vertical(PlotStyle):
    """
    グラフを縦に並べてプロットする
    """
    def __init__(self, figs_num):
        """
        figs_num: num of figs
        fig: object of pylab.figure
        axis_dict: labels of dict
        ax1, ax2, ax3: axis objects
        """
        PlotStyle.__init__(self, figs_num)
        self.figs_num = figs_num
        self.ax1, self.ax2, self.ax3 = None, None, None
        self.axis_list = [self.ax1, self.ax2, self.ax3]
        self.set_axis(figs_num)

    def set_axis(self, figs_num):
        """
        subplotを作成
        最大3つまで
        """
        pos = 11 + self.figs_num * 100
        self.ax1 = self.fig.add_subplot(pos)
        if figs_num >= 2:
            pos = 12 + self.figs_num * 100
            pylab.setp(self.ax1.get_xticklabels(), visible=False)
            self.ax2 = self.fig.add_subplot(pos, sharex=self.ax1)
        if figs_num >= 3:
            pos = 13 + self.figs_num * 100
            self.ax3 = self.fig.add_subplot(pos, sharex=self.ax1)
            pylab.setp(self.ax2.get_xticklabels(), visible=False)
        self.fig.subplots_adjust(hspace=0.03)
        self.axis_list = [self.ax1, self.ax2, self.ax3]

    def set_single_data(self, axis, data, xaxis, yaxis):
        """
        データをセットする
        """
        axis = self.axis_list[axis-1]

        xdata = data[xaxis]
        ydata = data[yaxis]
        axis.plot(xdata, ydata, **self.style)

        xlabel = self.set_axis_label(xaxis)
        ylabel = self.set_axis_label(yaxis)
        axis.set_xlabel(xlabel, size='large')
        axis.set_ylabel(ylabel, size='large', multialignment='left')

    def set_axis_label(self, axis_name):
        """
        軸のラベルを設定
        axis_dictに入っていない場合は空白('')を返す
        """
        try:
            label = self.axis_dict[axis_name]
        except KeyError:
            label = ''
        return label

    @staticmethod
    def plot(dst):
        """
        出力する
        """
        if not dst:
            pass
        elif dst == 'show':
            pylab.show()
        else:
            pylab.savefig(dst, bbox_inches="tight", pad_inches=0.1)
            pylab.close()
            pylab.cla()
            pylab.clf()

    def adjust_yscale(self, axis):
        """
        y軸のスケールを10%広げる
        """
        axis = self.axis_list[axis-1]
        deltay = (axis.get_ylim()[1] - axis.get_ylim()[0]) / 10
        axis.set_ylim(axis.get_ylim()[0]-deltay, axis.get_ylim()[1]+deltay)

    def adjust_auto(self):
        """
        全てのy軸についてスケールを10%広げる
        """
        for i in range(0, self.figs_num):
            self.adjust_yscale(i+1)

    def set123(self, data, xaxis, *yaxis):
        """
        最大3つまでのデータをsetする
        """
        for i in range(0, len(yaxis)):
            self.set_single_data(i+1, data, xaxis, yaxis[i])


def test():
    """
    For test
    """
    xdata = [1, 2]
    ydata = [3, 1]
    zdata = [2, 1]

    data = {'c/a': xdata, 'energy': ydata, 'mag': zdata}
    plot = Vertical(2)
    plot.set_style('blue')
    plot.set123(data, 'c/a', 'energy', 'mag')
    plot.set_style('magenta')
    plot.set_single_data(2, data, 'c/a', 'energy')

    pylab.show()


if __name__ == '__main__':
    test()
