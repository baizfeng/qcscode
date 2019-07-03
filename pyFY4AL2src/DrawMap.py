#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/8/6 17:06
# @Author  : baizhaofeng
# @File    : DrawMap.py

import matplotlib

matplotlib.use("Pdf")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import numpy.ma as ma


class CMap:
    fig = plt.figure(figsize=(27.48, 32), dpi=50)  # 图像的长*高=2748*3000
    axes2 = fig.add_axes([0., 0., 1., 0.14125], facecolor='#e8e8e8')
    axes1 = fig.add_axes([0., 0.14125, 1., 0.85875])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高

    # 加气象局logo
    axicon = fig.add_axes([0.85, 0.01, 0.15, 0.05])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'), origin='upper')
    axicon.axis('off')

    def __init__(self, mainTitle, dataName, lat, lon):
        self.mainTitle = mainTitle
        self.dataName = dataName
        self.lat = lat
        self.lon = lon

    def CreateMap(self):
        m = Basemap(projection='nsper', lat_0=0, lon_0=104.7, resolution='l', ax=CMap.axes1)
        m.drawmapboundary(linewidth=0.5)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)

        m.drawparallels(range(-90, 90, 10), linewidth=1.)
        m.drawmeridians(range(0, 360, 10), linewidth=1.)
        if (self.dataName in ('LST', 'LSE')):
            m.drawlsmask(land_color='#bebebe', ocean_color='#01008a')  # 陆地，海洋的颜色

        CMap.axes2.spines['left'].set_visible(False)
        CMap.axes2.spines['top'].set_visible(False)

        # 标题
        CMap.axes1.set_xlabel(self.mainTitle,
                              family='Times New Roman',
                              fontsize=42,
                              fontweight='bold', labelpad=20)
        return m

    def AddDiscMap(self, m, data, values, colors, names):
        """
        为底图叠加离散型数据
        :param m: 底图
        :param data: 离散数据
        :param values: 数据分级值
        :param colors: 颜色
        :param names: 数据各个值代表的名称
        :return:
        """
        x, y = m(self.lon, self.lat)
        m.contourf(x, y, data, values, colors=colors)

        for i in range(len(colors)):
            if (i > 4):
                xx = 0.1 + 0.15 * (i - 5)
                yy = 0.1
            else:
                xx = 0.1 + 0.15 * i
                yy = 0.4
            CMap.axes2.add_patch(plt.Rectangle((xx, yy), 0.04, 0.2, color=colors[i]))
            CMap.axes2.text(xx + 0.04, yy + 0.05, names[i], fontsize=28)
        return True

    def AddContiMap(self, m, data, cmap, levels, dataUnit):
        """
        为底图添加连续型数据
        :param m: 底图
        :param data: 连续型数据
        :param cmap: 图例颜色
        :param levels: 图例分级
        :param fig:
        :param ax:
        :param dataUnit:
        :return:
        """
        cax = CMap.fig.add_axes([0.1, 0.08, 0.7, 0.02])
        x, y = m(self.lon, self.lat)
        cs = m.contourf(x, y, data, cmap=cmap, levels=levels)
        cb = plt.colorbar(cs, cax=cax, orientation='horizontal', format="%.2f")
        cb.ax.tick_params(labelsize=28)
        CMap.axes2.text(0.801, 0.6, dataUnit, fontsize=32)  # 添加单位
        return True

    def AddBarbsMap(self, m, data):
        """
        为底图添加连续型数据
        :param m: 底图
        :param data: 连续型数据
        :return:
        """
        wind_speed = data[0]
        wind_direction = data[1]
        pressure = data[2]
        wind_mdirection = ma.masked_values(wind_direction, -999)
        wind_mspeed = ma.masked_values(wind_speed, -999)
        mpressure = ma.masked_values(pressure, -999)

        x, y = m(self.lon, self.lat)

        U = -wind_mspeed * np.sin(np.deg2rad(wind_mdirection)) * 2.5
        V = -wind_mspeed * np.cos(np.deg2rad(wind_mdirection)) * 2.5

        m.barbs(x[mpressure < 400], y[mpressure < 400], U[mpressure < 400], V[mpressure < 400],
                pivot='middle', flagcolor='red', lw=0.5, length=5)
        m.barbs(x[(mpressure >= 400) & (mpressure <= 700)], y[(mpressure >= 400) & (mpressure <= 700)],
                U[(mpressure >= 400) & (mpressure <= 700)],
                V[(mpressure >= 400) & (mpressure <= 700)], pivot='middle', flagcolor='#00ff00', lw=0.5, length=5)
        m.barbs(x[mpressure > 700], y[mpressure > 700], U[mpressure > 700], V[mpressure > 700],
                pivot='middle', flagcolor='#0009fa', lw=0.5, length=5)
        colors = ["red", "#00ff00", "#0009fa"]
        names = ["<400hPa", "400~700hPa", ">700hPa"]
        for i in range(len(colors)):
            if (i > 5):
                xx = 0.1 + 0.15 * (i - 6)
                yy = 0.1
            else:
                xx = 0.1 + 0.15 * i
                yy = 0.4
            CMap.axes2.add_patch(plt.Rectangle((xx, yy), 0.04, 0.2, color=colors[i]))
            CMap.axes2.text(xx + 0.04, yy + 0.05, names[i], fontsize=58)
        return True

    def AddSpecial(self, m, data, varName):
        x, y = m(self.lon, self.lat)
        if (varName == "RSR"):
            m.contourf(x, y, data, [65531.1, 65532.1], colors=["#c8c8c8"])
            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color="#c8c8c8"))
            CMap.axes2.text(0.151, 0.21, "Night", fontsize=32)
        elif (varName in ('CPD_COT', 'CPD_CER', 'CPD_IWP', 'CPD_LWP', 'CPN_COT', 'CPN_CER', 'CPN_IWP', 'CPN_LWP')):
            dcolors = ["#c8c8c8", "#9f9f9f", "#7f7f7f"]
            m.contourf(x, y, data, [65526.1, 65527.1, 65531.1, 65532.1], colors=dcolors)
            m.contourf(x, y, data, [-999.1, -998.1], colors=['black'])

            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color=dcolors[0]))
            CMap.axes2.text(0.151, 0.21, "Sea Ice", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.22, 0.2), 0.05, 0.15, color=dcolors[1]))
            CMap.axes2.text(0.271, 0.21, "Clear", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.34, 0.2), 0.05, 0.15, color=dcolors[2]))
            CMap.axes2.text(0.391, 0.21, "Night", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.46, 0.2), 0.05, 0.15, color="black"))
            CMap.axes2.text(0.511, 0.21, "Invalid", fontsize=32)
        elif (varName == "SSI"):
            m.contourf(x, y, data, [65531.1, 65532.1], colors=["#c8c8c8"])
            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color="#c8c8c8"))
            CMap.axes2.text(0.151, 0.21, "Sunzenith Angle>70", fontsize=32)
        elif (varName in ("CTH", "CTP", "CTT")):
            m.contourf(x, y, data, [-999.1, -998.1], colors=["#c8c8c8"])
            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color="#c8c8c8"))
            CMap.axes2.text(0.151, 0.21, "Invalid", fontsize=32)
        elif (varName in ("NIR Black sky albedo", "NIR White sky albedo", "SW Black sky albedo",
                          "SW White sky albedo", "VIS Black sky albedo", "VIS White sky albedo")):
            m.contourf(x, y, data.data, [-999.1, -998.1], colors=["#c8c8c8"])
            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color="#c8c8c8"))
            CMap.axes2.text(0.151, 0.21, "Invalid", fontsize=32)
        elif (varName in ("LPW_HIGH", "LPW_LOW", "LPW_MID", "TPW")):
            m.contourf(x, y, data, [65533.1, 65534.1], colors=["#c8c8c8"])
            m.contourf(x, y, data, [-999.1, -998.1], colors=['black'])

            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color="#c8c8c8"))
            CMap.axes2.text(0.151, 0.21, "Cloud", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.22, 0.2), 0.05, 0.15, color='black'))
            CMap.axes2.text(0.271, 0.21, "Invalid", fontsize=32)
        elif (varName == 'DSD'):
            dst = ma.masked_values(data, 32766.)
            dst = ma.masked_values(dst, 32764.)
            m.contourf(x, y, dst, [10, 25], colors=['red'])
            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color="red"))
            CMap.axes2.text(0.151, 0.21, "DUST AREA(DST>=10)", fontsize=32)
        elif (varName == 'AOD'):
            dcolors = ["#000188", "#a65026", "black", "#ffffff", "#575757"]
            m.contourf(x, y, data, [65529.1, 65530.1, 65531.1, 65532.1, 65533.1, 65534.1], colors=dcolors)
            m.contourf(x, y, data, [-999.1, -998.1], colors=['#c8c8c8'])

            CMap.axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.05, 0.15, color=dcolors[0]))
            CMap.axes2.text(0.151, 0.21, "Sea", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.22, 0.2), 0.05, 0.15, color=dcolors[1]))
            CMap.axes2.text(0.271, 0.21, "Land", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.34, 0.2), 0.05, 0.15, color=dcolors[2]))
            CMap.axes2.text(0.391, 0.21, "Night", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.46, 0.2), 0.05, 0.15, color=dcolors[3]))
            CMap.axes2.text(0.511, 0.21, "Cloud", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.58, 0.2), 0.05, 0.15, color=dcolors[4]))
            CMap.axes2.text(0.631, 0.21, "SatZen>72", fontsize=32)

            CMap.axes2.add_patch(plt.Rectangle((0.74, 0.2), 0.05, 0.15, color="#c8c8c8"))
            CMap.axes2.text(0.791, 0.21, "Invalid", fontsize=32)
        elif (varName == 'FHS'):
            m.scatter(x[data == 10], y[data == 10], s=100, color="red")
        return True
