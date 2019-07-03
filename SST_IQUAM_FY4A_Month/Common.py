#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/8/17 11:26
# @Author  : baizhaofeng
# @File    : Common.py


import collections
from xml.dom import minidom
import os
import datetime
from netCDF4 import Dataset
import numpy.ma as ma
import numpy as np
import matplotlib

matplotlib.use("Pdf")
import matplotlib.pyplot as plt
from scipy import spatial, stats
import re
import matplotlib.mlab as mlab
# import seaborn as sns
from mpl_toolkits.basemap import Basemap


def QCSlog(info, xmlPath):
    logpath = xmlPath[:-4] + ".log"
    infos = "%s  %s \n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), info)
    print infos
    with open(logpath, 'a') as f:
        f.write(infos)


def Read_xml(xmlpath):
    """
解析XML文件
    :param xmlpath:xml路径
    :return:xml文件中的结点，字典类型
    """
    with open(xmlpath, 'r') as fh:
        dom = minidom.parse(fh)
        root = dom.documentElement

        xmlnodes = collections.OrderedDict()

        filenode = root.getElementsByTagName('inputFileName')[0]
        xmlnodes['filepath'] = filenode.childNodes[0].nodeValue

        txtnode = root.getElementsByTagName('txtPath')[0]
        xmlnodes['txtpath'] = txtnode.childNodes[0].nodeValue

        typenode = root.getElementsByTagName('type')[0]
        xmlnodes['type'] = typenode.childNodes[0].nodeValue

        outputnode = root.getElementsByTagName('outputPath')[0]
        xmlnodes['outputpath'] = outputnode.childNodes[0].nodeValue

        outputxmlnode = root.getElementsByTagName('outputXml')[0]
        xmlnodes['xmlPath'] = outputxmlnode.childNodes[0].nodeValue

        warnnode = root.getElementsByTagName('warningValue')[0]
        xmlnodes['warnvalue'] = warnnode.childNodes[0].nodeValue

        return xmlnodes


def GetBaseInfo(inputFileName):
    """
    得到数据基本信息：时间、数据名称、数据源名称、分辨率、数据单位
    :param filename:
    :return:
    """
    filename = os.path.basename(inputFileName)
    fy4_time = datetime.datetime.strptime(filename.split('_')[9][:-8], '%Y%m')  # 时间
    dataRes = filename[74:79]
    baseInfo = {"time": fy4_time,
                "dataRes": dataRes,
                "filename": filename[:-3]}
    return baseInfo


def Name(filename, fytime, outputpath):
    scattername = "%s%s_QCS_IQUAM_SCAT_ORIG.png" % (outputpath, filename)
    histname = "%s%s_QCS_IQUAM_HIST_BIAS.png" % (outputpath, filename)
    iquamname = "%s%s_QCS_IQUAM_MAPG_ORIG.png" % (outputpath, filename)
    fy4name = "%s%s_QCS_MAPN_ORIG.png" % (outputpath, filename)
    scattertitle = "FY4A_IQUAM_SST %s" % (fytime.__format__('%Y-%m'))
    histtitle = "FY4A_IQUAM_SST %s" % (fytime.__format__('%Y-%m'))
    iquamtitle = "IQUAM_SST_MAPG %s" % (fytime.__format__('%Y-%m'))  # 标题
    fy4title = "FY4A_SST_MAPN %s" % (fytime.__format__('%Y-%m'))  # 标题
    names = {'scatname': scattername,
             'histname': histname,
             'iquamname': iquamname,
             'fy4name': fy4name,
             'scattitle': scattertitle,
             'histtitle': histtitle,
             'iquamtitle': iquamtitle,
             'fy4title': fy4title}
    return names


def MatchTXT(fy4time, nodes):
    txtpath = nodes['txtpath']
    pattern = ['QCS_IQUAM_DATA_BIAS.txt', nodes['type'], fy4time.__format__('%Y%m')]
    result = []  # 所有的文件

    for filename in os.listdir(txtpath):
        if (pattern[0] in filename) and (pattern[1] in filename) and (pattern[2] in filename):
            result.append(filename)
    return txtpath + result[0]


def ReadTXTData(filepath):
    with open(filepath, 'r') as r:
        lines = r.readlines()
    with open(filepath, 'w') as w:
        for l in lines:
            lst = l.split(' ')
            if len(lst) == 7:
                w.write(l)

    arr = np.loadtxt(filepath, delimiter=" ", usecols=(1, 2, 3, 4, 5, 6), dtype=np.float)
    arr = np.unique(arr, axis=0)
    exSST = arr[:, 0]
    fySST = arr[:, 1]
    fyDQF = arr[:, 2]
    plon = arr[:, 3]
    plat = arr[:, 4]
    ptype = arr[:, 5]

    bias = fySST[np.where(fySST != -999)] - exSST[np.where(fySST != -999)]

    mdat = {"exSST": exSST, "fySST": fySST, "bias": bias, "fyDQF": fyDQF, "plon": plon, "plat": plat, "ptype": ptype}
    return mdat


def Read_Data(filepath):
    """
    读取数据，根据变量名
    :param filepath:
    :return: 数据，array类型
    """
    # 读取FY4数据
    with Dataset(filepath, 'r') as fid:
        fy4_data = fid.variables['SST Monthly Mean Product'][:]

    if (isinstance(fy4_data, ma.MaskedArray)):
        fy4_data = fy4_data.data
    fy4_data[(fy4_data >= 65530) | (fy4_data == -888)] = -999
    fy4_mdata = ma.masked_values(fy4_data, -999)
    return fy4_mdata


def Read_LUT(prjType, filepath, dataRes):
    """
    根据数据的投影类型，分辨率读取经纬度查找表
    :param prjType:
    :param filepath:
    :param dataRes:
    :return:
    """
    if ((prjType == 'GEO') | (prjType == 'NUL')):
        lutpath = filepath
        with Dataset(lutpath, 'r') as fid:
            lat = fid.variables['lat'][:]
            lon = fid.variables['lon'][:]
    elif (prjType == 'NOM'):
        if (dataRes == '4000M'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/FullMask_Grid_4000_1047E_BSQ.nc'
        elif (dataRes == '1000M'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/FullMask_Grid_1000_1047E_BSQ.nc'
        elif (dataRes == '012KM'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/FullMask_Grid_12000_1047E_BSQ.nc'

        with Dataset(lutpath, 'r') as fid:
            lat = fid.variables['LAT'][:]
            lon = fid.variables['LON'][:]
        lat = ma.masked_values(lat, -999.)
        lon = ma.masked_values(lon, -999.)
    return lat, lon


def Calcate(mdat):
    exSST = mdat['exSST']
    fySST = mdat['fySST']
    fyDQF = mdat['fyDQF']
    metrics = collections.OrderedDict()
    bias0 = fySST[np.where(fyDQF == 0)] - exSST[np.where(fyDQF == 0)]
    bias1 = fySST[np.where(fyDQF == 1)] - exSST[np.where(fyDQF == 1)]
    bias2 = fySST[np.where(fyDQF == 2)] - exSST[np.where(fyDQF == 2)]
    bias3 = fySST[np.where(fySST != -999)] - exSST[np.where(fySST != -999)]

    if (bias0.shape[0] != 0):
        slope0, intercept0, r_value0, p_value0, std_err0 = stats.linregress(fySST[np.where(fyDQF == 0)],
                                                                            exSST[np.where(fyDQF == 0)])
    else:
        slope0, intercept0, r_value0, p_value0, std_err0 = [np.nan, np.nan, np.nan, np.nan, np.nan]
    if (bias1.shape[0] != 0):
        slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(fySST[np.where(fyDQF == 1)],
                                                                            exSST[np.where(fyDQF == 1)])
    else:
        slope1, intercept1, r_value1, p_value1, std_err1 = [np.nan, np.nan, np.nan, np.nan, np.nan]
    if (bias2.shape[0] != 0):
        slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(fySST[np.where(fyDQF == 2)],
                                                                            exSST[np.where(fyDQF == 2)])
    else:
        slope2, intercept2, r_value2, p_value2, std_err2 = [np.nan, np.nan, np.nan, np.nan, np.nan]
    if (bias3.shape[0] != 0):
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(fySST[np.where(fySST != -999)],
                                                                            exSST[np.where(fySST != -999)])
    else:
        slope3, intercept3, r_value3, p_value3, std_err3 = [np.nan, np.nan, np.nan, np.nan, np.nan]

    num0 = bias0.shape[0]  # 样本数量
    mean0 = round(np.mean(bias0), 4)  # 平均偏差
    absmean0 = round(np.mean(np.abs(bias0)), 4)  # AE绝对平均偏差
    std0 = round(np.sqrt(np.square(bias0).sum() / (num0-1)), 4)  # 标准差
    max0 = round(np.max(bias0), 4) if num0 != 0 else np.nan  # 最大值
    min0 = round(np.min(bias0), 4) if num0 != 0 else np.nan  # 最小值
    median0 = round(np.median(bias0), 4)  # 中位数
    rmse0 = round(np.sqrt(np.square(bias0).sum() / num0), 4)  # 均方根误差
    skew0 = round(stats.skew(bias0), 4)  # 偏度系数
    kurt0 = round(stats.kurtosis(bias0), 4)  # 峰度系数

    num1 = bias1.shape[0]  # 样本数量
    mean1 = round(np.mean(bias1), 4)  # 平均偏差
    absmean1 = round(np.mean(np.abs(bias1)), 4)  # AE绝对平均偏差
    std1 = round(np.sqrt(np.square(bias1).sum() / (num1-1)), 4)  # 标准差
    max1 = round(np.max(bias1), 4) if num1 != 0 else np.nan  # 最大值
    min1 = round(np.min(bias1), 4) if num1 != 0 else np.nan  # 最小值
    median1 = round(np.median(bias1), 4)  # 中位数
    rmse1 = round(np.sqrt(np.square(bias1).sum() / num1), 4)  # 均方根误差
    skew1 = round(stats.skew(bias1), 4)  # 偏度系数
    kurt1 = round(stats.kurtosis(bias1), 4)  # 峰度系数

    num2 = bias2.shape[0]  # 样本数量
    mean2 = round(np.mean(bias2), 4)  # 平均偏差
    absmean2 = round(np.mean(np.abs(bias2)), 4)  # AE绝对平均偏差
    std2 = round(np.sqrt(np.square(bias2).sum() / (num2-1)), 4)  # 标准差
    max2 = round(np.max(bias2), 4) if num2 != 0 else np.nan  # 最大值
    min2 = round(np.min(bias2), 4) if num2 != 0 else np.nan  # 最小值
    median2 = round(np.median(bias2), 4)  # 中位数
    rmse2 = round(np.sqrt(np.square(bias2).sum() / num2), 4)  # 均方根误差
    skew2 = round(stats.skew(bias2), 4)  # 偏度系数
    kurt2 = round(stats.kurtosis(bias2), 4)  # 峰度系数

    num3 = bias3.shape[0]  # 样本数量
    mean3 = round(np.mean(bias3), 4)  # 平均偏差
    absmean3 = round(np.mean(np.abs(bias3)), 4)  # AE绝对平均偏差
    std3 = round(np.sqrt(np.square(bias3).sum() / (num3-1)), 4)  # 标准差
    max3 = round(np.max(bias3), 4) if num3 != 0 else np.nan  # 最大值
    min3 = round(np.min(bias3), 4) if num3 != 0 else np.nan  # 最小值
    median3 = round(np.median(bias3), 4)  # 中位数
    rmse3 = round(np.sqrt(np.square(bias3).sum() / num3), 4)  # 均方根误差
    skew3 = round(stats.skew(bias3), 4)  # 偏度系数
    kurt3 = round(stats.kurtosis(bias3), 4)  # 峰度系数

    metrics['QualID0_NUM'] = num0
    metrics['QualID0_MAX'] = max0
    metrics['QualID0_MIN'] = min0
    metrics['QualID0_MEDIAN'] = median0
    metrics['QualID0_MEAN'] = mean0
    metrics['QualID0_AE'] = absmean0
    metrics['QualID0_STD'] = std0
    metrics['QualID0_RMSE'] = rmse0
    metrics['QualID0_SKEW'] = skew0
    metrics['QualID0_KURT'] = kurt0
    metrics['QualID0_CORR'] = r_value0

    metrics['QualID1_NUM'] = num1
    metrics['QualID1_MAX'] = max1
    metrics['QualID1_MIN'] = min1
    metrics['QualID1_MEDIAN'] = median1
    metrics['QualID1_MEAN'] = mean1
    metrics['QualID1_AE'] = absmean1
    metrics['QualID1_STD'] = std1
    metrics['QualID1_RMSE'] = rmse1
    metrics['QualID1_SKEW'] = skew1
    metrics['QualID1_KURT'] = kurt1
    metrics['QualID1_CORR'] = r_value1

    metrics['QualID2_NUM'] = num2
    metrics['QualID2_MAX'] = max2
    metrics['QualID2_MIN'] = min2
    metrics['QualID2_MEDIAN'] = median2
    metrics['QualID2_MEAN'] = mean2
    metrics['QualID2_AE'] = absmean2
    metrics['QualID2_STD'] = std2
    metrics['QualID2_RMSE'] = rmse2
    metrics['QualID2_SKEW'] = skew2
    metrics['QualID2_KURT'] = kurt2
    metrics['QualID2_CORR'] = r_value2

    metrics['QualID3_NUM'] = num3
    metrics['QualID3_MAX'] = max3
    metrics['QualID3_MIN'] = min3
    metrics['QualID3_MEDIAN'] = median3
    metrics['QualID3_MEAN'] = mean3
    metrics['QualID3_AE'] = absmean3
    metrics['QualID3_STD'] = std3
    metrics['QualID3_RMSE'] = rmse3
    metrics['QualID3_SKEW'] = skew3
    metrics['QualID3_KURT'] = kurt3
    metrics['QualID3_CORR'] = r_value3

    param = {"slope0": slope0,
             "slope1": slope1,
             "slope2": slope2,
             "slope3": slope3,
             "intercept0": intercept0,
             "intercept1": intercept1,
             "intercept2": intercept2,
             "intercept3": intercept3,
             "QualID0_CORR": r_value0,
             "QualID1_CORR": r_value1,
             "QualID2_CORR": r_value2,
             "QualID3_CORR": r_value3}

    return metrics, param


def DrawScatter(mdat, param, totalNum, names):
    exSST = mdat['exSST']
    fySST = mdat['fySST']
    fyDQF = mdat['fyDQF']

    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)

    fig.suptitle(names['scattitle'], fontsize=14, fontweight='bold')
    H0 = ax.scatter(fySST[np.where(fyDQF == 0)], exSST[np.where(fyDQF == 0)], c='r', marker='o')
    H1 = ax.scatter(fySST[np.where(fyDQF == 1)], exSST[np.where(fyDQF == 1)], c='g', marker='D')
    H2 = ax.scatter(fySST[np.where(fyDQF == 2)], exSST[np.where(fyDQF == 2)], c='b', marker='x')

    ax.set_xlabel('FY4A SST', fontsize=12)
    ax.set_ylabel('IQUAM SST', fontsize=12)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')

    ax.plot([-5, 45], [-5, 45], color='black', zorder=3)  # 画参考线
    plt.xlim(-5, 45)
    plt.ylim(-5, 45)
    ax.legend([H0, H1, H2], ['Best', 'Good', 'Bad'], loc="lower right")

    txtstr = "N=%s\nCORR=%s\ny=%s+%sx" % (totalNum, round(param['QualID3_CORR'], 4), round(param['intercept3'], 4),
                                          round(param['slope3'], 4))
    ax.text(0.04, 1., txtstr,
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=10)

    plt.savefig(names['scatname'])
    plt.close
    return True


def DrawHist(bias, metric, names):
    # 创建直方图
    # 第一个参数为待绘制的定量数据
    # 第二个参数为划分的区间个数
    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)  # 输出图片的大小
    fig.suptitle(names['histtitle'], fontsize=14, fontweight='bold')
    n, bins, patches = ax.hist(bias, bins=20, density=True, edgecolor='k', facecolor='white', zorder=2)  # bins是分为20个柱子
    # add a 'best fit' line
    y = mlab.normpdf(bins, metric['QualID3_MEAN'], metric['QualID3_STD'])
    l = ax.plot(bins, y, '--', color='black', linewidth=1)

    # kde=True表示是否显示拟合曲线，如果为False则只出现直方图,kde_kws（拟合曲线的设置）、hist_kws（直方柱子的设置）
    # sns.distplot(bias, bins=20, kde=True,
    #              kde_kws={"color": "black", "lw": 1, 'linestyle': '--'},
    #              hist_kws={"histtype": 'bar', "edgecolor": 'black', "facecolor": 'white', "lw": 1})

    # plt.xlim(-5, 5)
    ax.set_xlabel('FY4A-iQUAM(' + u'\u2103' + ')', fontsize=12)
    ax.set_ylabel('Number Density', fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    # ax.grid(True, linestyle='dashed')

    txtstr = "N=%s\n\nMIN=%s\nMAX=%s\n\nMEAN=%s\nSTD=%s\n\nAE=%s\nMEDIAN=%s\n\nRMSE=%s\nSKEW=%s\nKURT=%s" % (
        metric['QualID3_NUM'], metric['QualID3_MIN'], metric['QualID3_MAX'], metric['QualID3_MEAN'],
        metric['QualID3_STD'], metric['QualID3_AE'], metric['QualID3_MEDIAN'], metric['QualID3_RMSE'],
        metric['QualID3_SKEW'], metric['QualID3_KURT'])

    ax.text(0.04, 1., txtstr,
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=10)

    plt.savefig(names['histname'])
    plt.close
    return True


def DrawIQUAMMap(mdat, names):
    data = mdat['exSST']
    lat = mdat['plat']
    lon = mdat['plon']
    ptype = mdat['ptype']
    fig = plt.figure(figsize=(24, 21), dpi=100)  # 图像的长*高=2748*3000
    axes2 = fig.add_axes([0., 0., 1., 0.08], facecolor='#e8e8e8')
    axes1 = fig.add_axes([0.04, 0.11, 0.96, 0.88])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高
    cax = fig.add_axes([0.3, 0.03, 0.5, 0.02])  # 图例

    m = Basemap(projection='cea', resolution='l', ax=axes1, llcrnrlat=-60, llcrnrlon=44, urcrnrlon=165, urcrnrlat=60)
    x, y = m(lon, lat)
    m.drawcoastlines()
    m.drawcountries()
    m.drawlsmask(land_color='#bebebe')  # 陆地的颜色
    m.drawparallels(range(-90, 90, 15), labels=[1, 0, 0, 0], fontsize=20)
    m.drawmeridians(range(0, 360, 20), labels=[0, 0, 0, 1], fontsize=20)

    cs1 = m.scatter(x[ptype == 1], y[ptype == 1], c=data[ptype == 1], s=200, marker='o', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)  # s指定散点图点的大小
    cs2 = m.scatter(x[ptype == 2], y[ptype == 2], c=data[ptype == 2], s=200, marker='D', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)
    cs3 = m.scatter(x[ptype == 3], y[ptype == 3], c=data[ptype == 3], s=200, marker='x', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)
    cs4 = m.scatter(x[ptype == 4], y[ptype == 4], c=data[ptype == 4], s=200, marker='s', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)

    cs5 = m.scatter(x[ptype == 5], y[ptype == 5], c=data[ptype == 5], s=200, marker='+', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)
    cs6 = m.scatter(x[ptype == 6], y[ptype == 6], c=data[ptype == 6], s=200, marker='*', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)
    cs7 = m.scatter(x[ptype == 7], y[ptype == 7], c=data[ptype == 7], s=200, marker='1', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)
    cs8 = m.scatter(x[ptype == 8], y[ptype == 8], c=data[ptype == 8], s=200, marker='^', cmap=plt.cm.jet, vmin=-2,
                    vmax=36)
    cb = plt.colorbar(cs8, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=20)
    axes2.text(0.8, 0.4, "Unit:(" + u'\u2103' + ")", fontsize=22)  # 添加单位

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)

    fig.legend((cs1, cs2, cs3, cs4, cs5, cs6, cs7, cs8), (
        'Ship', 'Drifting buoy', 'Tropical moored buoy', 'Coastal moored buoy', 'Argo float', 'High res drifter',
        'Imos', 'Crw buoy'), ncol=2, loc='lower left', facecolor='#e8e8e8', fontsize=18)

    cax.set_title(names['iquamtitle'],
                  family='Times New Roman',
                  fontsize=32,
                  fontweight='bold',pad=10)
    #添加Logo
    axicon=fig.add_axes([0.85,0.01,0.15,0.05])	
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'),origin='upper')
    axicon.axis('off')	

    fig.savefig(names['iquamname'])
    plt.close()
    return 0


def DrawMap_NOM(fy4_mdata, lat, lon, names):
    fig = plt.figure(figsize=(27.48, 30), dpi=100)  # 图像的长*高=2748*3000
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高
    m = Basemap(projection='ortho', lat_0=0, lon_0=104.7, resolution='l', ax=axes1)
    m.drawmapboundary(color='black', linewidth=3.0)
    m.drawcoastlines(color='black', linewidth=3.0)
    m.drawcountries(color='black', linewidth=3.0)

    m.drawparallels(range(-90, 90, 10), color='black', linewidth=3.0)
    m.drawmeridians(range(0, 360, 10), color='black', linewidth=3.0)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)

    cax = fig.add_axes([0.15, 0.03, 0.7, 0.02])
    x, y = m(lon, lat)
    cs = m.contourf(x, y, fy4_mdata, cmap=plt.cm.jet, levels=np.arange(-2, 36, 2))
    cb = plt.colorbar(cs, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=28)
    axes2.text(0.851, 0.4, "Unit:(" + u'\u2103' + ")", fontsize=32)  # 添加单位

    # 标题
    axes1.set_xlabel(names['fy4title'],
                     family='Times New Roman',
                     fontsize=42,
                     fontweight='bold', labelpad=20)
    fig.savefig(names['fy4name'])
    plt.close()
    return 0


def LogXML(logMsg, logType, nodes, fyTime=True, metric=True, names=True):
    xmlfilePath = nodes['xmlPath']
    inputFilename = nodes['filepath']
    warnValue = float(nodes['warnvalue'])

    # 在内存中创建一个空的文档
    doc = minidom.Document()
    # 创建一个根节点Managers对象
    root = doc.createElement('XML')
    # 设置根节点的属性
    root.setAttribute('identify', 'QC')
    # 将根节点添加到文档对象中
    doc.appendChild(root)

    nodeOutputFiles = doc.createElement('OutputFiles')
    nodeOutputFile = doc.createElement('OutputFile')
    nodeOutputFiles.appendChild(nodeOutputFile)

    nodeInputFilename = doc.createElement('InputFilename')

    # 给叶子节点name设置一个文本节点，用于显示文本内容
    nodeInputFilename.appendChild(doc.createTextNode(str(os.path.basename(inputFilename))))

    if (names != True):
        nodeOutputFilename0 = doc.createElement("OutputFilename")
        nodeOutputFilename0.appendChild(doc.createTextNode(str(os.path.basename(names['fy4name']))))
        nodeOutputFile.appendChild(nodeOutputFilename0)

        nodeOutputFilename1 = doc.createElement("OutputFilename")
        nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['iquamname']))))
        nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeOutputFilename2 = doc.createElement("OutputFilename")
        nodeOutputFilename2.appendChild(doc.createTextNode(str(os.path.basename(names['scatname']))))
        nodeOutputFile.appendChild(nodeOutputFilename2)

        nodeOutputFilename3 = doc.createElement("OutputFilename")
        nodeOutputFilename3.appendChild(doc.createTextNode(str(os.path.basename(names['histname']))))
        nodeOutputFile.appendChild(nodeOutputFilename3)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', 'FY4A')
        nodeDeviationstatic.setAttribute('checkVariety', 'IQUAM')
        nodeDeviationstatic.setAttribute('dataFormat', 'month')
        nodeDeviationstatic.setAttribute('dataType', 'NOM')
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', 'SST')
        metricStr = "{%s}" % re.sub('[()\'\[\]]', '',
                                    str(metric).replace("',", ":").replace("OrderedDict", "")).replace(' ', '')
        nodeDeviationstatic.setAttribute('values', metricStr)
        nodeOutputFile.appendChild(nodeDeviationstatic)

        if metric['QualID3_AE'] > warnValue:
            warntype = 1 		
            nodewarnType = doc.createElement("warningType")
            nodewarnType.appendChild(doc.createTextNode(str(warntype)))
            nodewarnMsg = doc.createElement("warningMsg")
            nodewarnMsg.appendChild(doc.createTextNode("info" if not warntype else "精度预警：FY4A SST产品绝对偏差超过%s ℃" % warnValue))

            nodewarning = doc.createElement("warning")
            nodewarning.appendChild(nodewarnType)
            nodewarning.appendChild(nodewarnMsg)
            root.appendChild(nodewarning)

    # 将各叶子节点添加到父节点OutputFilename中，
    # 最后将OutputFilename添加到根节点XML中
    nodeOutputFile.appendChild(nodeInputFilename)
    root.appendChild(nodeOutputFiles)

    nodeLog = doc.createElement('log')

    nodeLogtype = doc.createElement('logType')
    nodeLogtype.appendChild(doc.createTextNode(str(logType)))

    nodeLogmsg = doc.createElement('logMsg')
    nodeLogmsg.appendChild(doc.createTextNode(str(logMsg)))
    nodeLog.appendChild(nodeLogtype)
    nodeLog.appendChild(nodeLogmsg)
    root.appendChild(nodeLog)

    # 开始写xml文档
    fp = open(xmlfilePath, 'w')
    doc.writexml(fp, indent='\t', addindent='\t', newl='\n', encoding="utf-8")
    fp.close()