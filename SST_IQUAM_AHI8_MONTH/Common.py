#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/5/7 16:50
# @Author  : baizhaofeng
# @Email   : baizhaofeng@piesat.cn
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
import seaborn as sns
from mpl_toolkits.basemap import Basemap
import struct
import calendar


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

        timenode = root.getElementsByTagName('timeEnd')[0]
        xmlnodes['time'] = timenode.childNodes[0].nodeValue

        txtnode = root.getElementsByTagName('inputFilePath')[0]
        xmlnodes['txtpath'] = txtnode.childNodes[0].nodeValue

        typenode = root.getElementsByTagName('type')[0]
        xmlnodes['type'] = typenode.childNodes[0].nodeValue

        outputnode = root.getElementsByTagName('outputPath')[0]
        xmlnodes['outputpath'] = outputnode.childNodes[0].nodeValue

        outputxmlnode = root.getElementsByTagName('outputXml')[0]
        xmlnodes['xmlPath'] = outputxmlnode.childNodes[0].nodeValue

        warnnode = root.getElementsByTagName('warningValue')[0]
        xmlnodes['warnvalue'] = warnnode.childNodes[0].nodeValue

        nanode = root.getElementsByTagName('isna')[0]
        xmlnodes['isna'] = nanode.childNodes[0].nodeValue

        return xmlnodes


def Name(datatime, outputpath, isna):
    ahi8time = datetime.datetime.strptime(datatime[:6], '%Y%m')
    # 获取当月第一天的星期和当月的总天数
    _, monthRange = calendar.monthrange(ahi8time.year, ahi8time.month)
    last = datetime.date(ahi8time.year, ahi8time.month, monthRange)
    filename = "AHI8-_AGRI--_N_DISK_1407E_L2-_SST-_MULT_NOM_%s_%s_4000M_V0001" % (
    ahi8time.__format__('%Y%m%d%H%M%S'), last.__format__(
        '%Y%m%d%H%M%S'))
    if (isna == 'false'):
        scattername = "%s%s_QCS_IQUAM_SCATMONTH_ORIG.png" % (outputpath, filename)
        histname = "%s%s_QCS_IQUAM_HISTMONTH_BIAS.png" % (outputpath, filename)
        iquamname = "%s%s_QCS_IQUAM_MAPGMONTH_ORIG.png" % (outputpath, filename)
        scattertitle = "AHI8_IQUAM_SST %s" % (ahi8time.__format__('%Y-%m'))
        histtitle = "AHI8_IQUAM_SST %s" % (ahi8time.__format__('%Y-%m'))
        iquamtitle = "IQUAM_SST_MAPG %s" % (ahi8time.__format__('%Y-%m'))  # 标题
    elif (isna == 'true'):
        scattername = "%s%s_QCS_IQUAM_SCATMONTH_ORIG_NA.png" % (outputpath, filename)
        histname = "%s%s_QCS_IQUAM_HISTMONTH_BIAS_NA.png" % (outputpath, filename)
        iquamname = "%s%s_QCS_IQUAM_MAPGMONTH_ORIG_NA.png" % (outputpath, filename)
        scattertitle = "AHI8_IQUAM_SST %s" % (ahi8time.__format__('%Y-%m'))
        histtitle = "AHI8_IQUAM_SST %s" % (ahi8time.__format__('%Y-%m'))
        iquamtitle = "IQUAM_SST_MAPG %s" % (ahi8time.__format__('%Y-%m'))  # 标题
    names = {'scatname': scattername,
             'histname': histname,
             'iquamname': iquamname,
             'scattitle': scattertitle,
             'histtitle': histtitle,
             'iquamtitle': iquamtitle}
    return names


def MatchTXT(fy4timestr, nodes):
    fy4time = datetime.datetime.strptime(fy4timestr[:6], '%Y%m')
    txtpath = nodes['txtpath']
    isna = nodes['isna']
    if (isna == 'false'):
        pattern = ['QCS_IQUAM_DATA_BIAS.txt', nodes['type'], fy4time.__format__('%Y%m')]
    elif (isna == 'true'):
        pattern = ['QCS_IQUAM_DATA_BIAS_NA.txt', nodes['type'], fy4time.__format__('%Y%m')]
    result = []  # 所有的文件

    for filename in os.listdir(txtpath):
        if (pattern[0] in filename) and (pattern[1] in filename) and (pattern[2] in filename):
            result.append(filename)
    return txtpath + result[0]


def ReadTXTData(filepath):
    with open(filepath, 'r') as r:
        lines = r.readlines()

    lsts = []
    for l in lines:
        lst = l.split(' ')
        if (len(lst) == 7):
            lsts.append(lst[1:7])
    arr = np.array(lsts, dtype=np.float)
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


def Calcate(mdat):
    exSST = mdat['exSST']
    fySST = mdat['fySST']
    metrics = collections.OrderedDict()
    bias3 = mdat['bias']

    if (bias3.shape[0] != 0):
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(fySST.flatten(), exSST.flatten())
    else:
        slope3, intercept3, r_value3, p_value3, std_err3 = [np.nan, np.nan, np.nan, np.nan, np.nan]

    num3 = bias3.shape[0]  # 样本数量
    mean3 = round(np.mean(bias3), 4)  # 平均偏差
    absmean3 = round(np.mean(np.abs(bias3)), 4)  # AE绝对平均偏差
    std3 = round(np.sqrt(np.square(bias3).sum() / (num3 - 1)), 4)  # 标准差
    max3 = round(np.max(bias3), 4) if num3 != 0 else np.nan  # 最大值
    min3 = round(np.min(bias3), 4) if num3 != 0 else np.nan  # 最小值
    median3 = round(np.median(bias3), 4)  # 中位数
    rmse3 = round(np.sqrt(np.square(bias3).sum() / num3), 4)  # 均方根误差
    skew3 = round(stats.skew(bias3), 4)  # 偏度系数
    kurt3 = round(stats.kurtosis(bias3), 4)  # 峰度系数

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

    param = {"slope3": slope3,
             "intercept3": intercept3,
             "QualID3_CORR": r_value3}

    return metrics, param


def DrawScatter(mdat, param, totalNum, names):
    exSST = mdat['exSST']
    fySST = mdat['fySST']

    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)

    fig.suptitle(names['scattitle'], fontsize=14, fontweight='bold')
    H0 = ax.scatter(fySST, exSST, c='r', marker='o', s=5)

    ax.set_xlabel('AHI8 SST', fontsize=12)
    ax.set_ylabel('IQUAM SST', fontsize=12)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')

    ax.plot([-5, 45], [-5, 45], color='black', zorder=3)  # 画参考线
    plt.xlim(-5, 45)
    plt.ylim(-5, 45)

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
    # n, bins, patches = ax.hist(bias, bins=20, density=True, edgecolor='k', facecolor='white', zorder=2)  # bins是分为20个柱子
    # # add a 'best fit' line
    # y = mlab.normpdf(bins, metric['QualID3_MEAN'], metric['QualID3_STD'])
    # l = ax.plot(bins, y, '--', color='black', linewidth=1)

    # kde=True表示是否显示拟合曲线，如果为False则只出现直方图,kde_kws（拟合曲线的设置）、hist_kws（直方柱子的设置）
    sns.distplot(bias, bins=30, kde=True,
                 kde_kws={"color": "black", "lw": 1, 'linestyle': '--'},
                 hist_kws={"histtype": 'bar', "edgecolor": 'black', "facecolor": 'white', "lw": 1})

    plt.xlim(-10, 10)
    ax.set_xlabel('AHI8-IQUAM(' + u'\u2103' + ')', fontsize=12)
    ax.set_ylabel('Number Density', fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')

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
    cax = fig.add_axes([0.3, 0.03, 0.6, 0.02])  # 图例

    m = Basemap(projection='cea', resolution='l', ax=axes1, llcrnrlat=-60, llcrnrlon=80, urcrnrlon=200, urcrnrlat=60)
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
    axes2.text(0.9, 0.4, "Unit:(" + u'\u2103' + ")", fontsize=22)  # 添加单位

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)

    fig.legend((cs1, cs2, cs3, cs4, cs5, cs6, cs7, cs8), (
        'Ship', 'Drifting buoy', 'Tropical moored buoy', 'Coastal moored buoy', 'Argo float', 'High res drifter',
        'Imos', 'Crw buoy'), ncol=2, loc='lower left', facecolor='#e8e8e8', fontsize=18)

    cax.set_title(names['iquamtitle'],
                  family='Times New Roman',
                  fontsize=32,
                  fontweight='bold')

    fig.savefig(names['iquamname'])
    plt.close()
    return 0


def LogXML(logMsg, logType, nodes, metric=True, names=True):
    fyTimestr = nodes['time']
    fyTime = datetime.datetime.strptime(fyTimestr[:6], '%Y%m')
    xmlfilePath = nodes['xmlPath']
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

    if (names != True):
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
        nodeDeviationstatic.setAttribute('bechekVariety', 'AHI8')
        nodeDeviationstatic.setAttribute('checkVariety', 'IQUAM')
        nodeDeviationstatic.setAttribute('dataFormat', 'month')
        nodeDeviationstatic.setAttribute('dataType', 'NOM')
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', 'AHI8#SST')
        metricStr = "{%s}" % re.sub('[()\'\[\]]', '',
                                    str(metric).replace("',", ":").replace("OrderedDict", "")).replace(' ', '')
        nodeDeviationstatic.setAttribute('values', metricStr)
        nodeOutputFile.appendChild(nodeDeviationstatic)

        if metric['QualID3_AE'] > warnValue:
            warntype = 1 		
            nodewarnType = doc.createElement("warningType")
            nodewarnType.appendChild(doc.createTextNode(str(warntype)))
            nodewarnMsg = doc.createElement("warningMsg")
            nodewarnMsg.appendChild(doc.createTextNode("info" if not warntype else "精度预警：AHI8 SST产品绝对偏差超过%s ℃" % warnValue))

            nodewarning = doc.createElement("warning")
            nodewarning.appendChild(nodewarnType)
            nodewarning.appendChild(nodewarnMsg)
            root.appendChild(nodewarning)

    # 将各叶子节点添加到父节点OutputFilename中，
    # 最后将OutputFilename添加到根节点XML中
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
