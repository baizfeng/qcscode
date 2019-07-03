#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/7/11 14:24
# @Author  : baizhaofeng
# @File    : Func.py

from xml.dom import minidom
import os
import datetime
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib

matplotlib.use("Pdf")
import matplotlib.pyplot as plt
import collections
from mpl_toolkits.basemap import Basemap
from scipy import ndimage
import re
import glob


def QCSlog(info, xmlPath):
    logpath = xmlPath[:-4] + ".log"
    infos = "%s  %s \n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), info)
    print infos
    with open(logpath, 'a') as f:
        f.write(infos)


# 解析XML文件
def Read_xml(xmlpath):
    xmlnodes = collections.OrderedDict()
    with open(xmlpath, 'r') as fh:
        dom = minidom.parse(fh)
        root = dom.documentElement

        filenode = root.getElementsByTagName('inputFileName')[0]
        xmlnodes['filepath'] = filenode.childNodes[0].nodeValue

        auxnode = root.getElementsByTagName('auxiliaryDataPath')[0]
        xmlnodes['auxpath'] = auxnode.childNodes[0].nodeValue

        outputnode = root.getElementsByTagName('outputPath')[0]
        xmlnodes['outputpath'] = outputnode.childNodes[0].nodeValue

        outputxmlnode = root.getElementsByTagName('outputXml')[0]
        xmlnodes['xmlPath'] = outputxmlnode.childNodes[0].nodeValue

        prjTypenode = root.getElementsByTagName('type')[0]
        xmlnodes['prjType'] = prjTypenode.childNodes[0].nodeValue

        timenode = root.getElementsByTagName('timeValue')[0]
        xmlnodes['timevalue'] = timenode.childNodes[0].nodeValue

        spacenode = root.getElementsByTagName('spaceValue')[0]
        xmlnodes['spacevalue'] = spacenode.childNodes[0].nodeValue

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
    dirname = os.path.dirname(inputFileName)
    fy4_time = datetime.datetime.strptime(filename.split('_')[9], '%Y%m%d%H%M%S')  # 时间
    dataName = filename[30:33]  # 数据名称
    dataSour = filename[0:4]  # 数据源名称
    dataRes = filename[74:79]
    baseInfo = {"dataName": dataName,
                "dataSour": dataSour,
                "dataRes": dataRes,
                "time": fy4_time,
                "filename": filename[:-3],
                "dirname": dirname}
    return baseInfo


def Name(infos, nodes, filenames, times):
    filename = infos['filename']
    dirname = infos['dirname']
    time0 = infos['time']
    prjType = nodes['prjType']
    outputpath = nodes['outputpath']
    auxpath = nodes['auxpath']

    if (prjType == 'NOM'):
        prj = "MAPN"
    elif (prjType == "GLL"):
        prj = "MAPG"

    t1filepath = "%s/%s" % (dirname, filenames[0])
    t2filepath = "%s/%s" % (dirname, filenames[1])
    tbb_t0path = "%s/%s/%s/%s.NC" % (
        auxpath, time0.__format__('%Y'), time0.__format__('%Y%m%d'), filename.replace('CIX', 'TBB'))
    tbb_t1path = "%s/%s/%s/%s" % (
        auxpath, times[0].__format__('%Y'), times[0].__format__('%Y%m%d'), filenames[0].replace('CIX', 'TBB'))
    tbb_t2path = "%s/%s/%s/%s" % (
        auxpath, times[1].__format__('%Y'), times[1].__format__('%Y%m%d'), filenames[1].replace('CIX', 'TBB'))
    histname = "%s%s_QCS_FY4A_HIST.png" % (outputpath, filename)
    histtitle = "FY4A_FY4A_CIX %s" % (time0.__format__('%Y-%m-%d %H:%M:%S'))

    t0map = "%s%s_QCS_FY4A_%s_ORIG.png" % (outputpath, filename, prj)
    t0title = "FY4A_CIX_%s %s" % (prj, time0.__format__('%Y-%m-%d %H:%M:%S'))  # 标题
    t1map = "%s%s_QCS_%s_ORIG.png" % (outputpath, filenames[0][:-3], prj)
    t1title = "FY4A_CIX_%s %s" % (prj, times[0].__format__('%Y-%m-%d %H:%M:%S'))  # 标题
    t2map = "%s%s_QCS_%s_ORIG.png" % (outputpath, filenames[1][:-3], prj)
    t2title = "FY4A_CIX_%s %s" % (prj, times[1].__format__('%Y-%m-%d %H:%M:%S'))  # 标题

    names = {"t1filepath": t1filepath,
             "t2filepath": t2filepath,
             "tbb_t0path": tbb_t0path,
             "tbb_t1path": tbb_t1path,
             "tbb_t2path": tbb_t2path,
             "histname": histname,
             "histtitle": histtitle,
             "t0map": t0map,
             "t1map": t1map,
             "t2map": t2map,
             "t0title": t0title,
             "t1title": t1title,
             "t2title": t2title}
    return names


# 时间匹配：
# path-待匹配数据路径
# base_time-匹配的时间
# timeThreshold-时间阈值，FY4时间向前找雷达时间半个小时
def TimeMatch(filepath, base_time, timeThreshold):
    files = glob.glob(os.path.join(filepath, '*REGC*.NC'))
    dfiles = []
    dtimes = []
    for file in files:
        file = os.path.basename(file)
        timeStr = file.split('_')[9]
        time = datetime.datetime.strptime(timeStr, '%Y%m%d%H%M%S')
        if ((time - base_time).total_seconds() / 60 <= float(timeThreshold)) & (
                (time - base_time).total_seconds() / 60 > 0):
            dfiles.append(file)
            dtimes.append(time)
    dfiles.sort()
    dtimes.sort()
    return dfiles, dtimes


# 读取数据
def ReadDatas(fy4_path):
    # 读取FY4数据
    with Dataset(fy4_path, 'r') as fid_fy4:
        fy4_data = fid_fy4.variables['RDC_Rank'][:]  # Convective_Initiation,RDC_Rank
        fy4_data.dtype = np.int32
    fy4_data = fy4_data.byteswap()
    return fy4_data


def Read_LUT(prjType, filepath, dataRes, rows):
    """
    根据数据的投影类型，分辨率读取经纬度查找表
    :param prjType:
    :param filepath:
    :param dataRes:
    :return:
    """
    if (prjType == 'GLL'):
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
            lat = fid.variables['LAT'][:][:rows, :]
            lon = fid.variables['LON'][:][:rows, :]
        lon[np.where(lon > 180.)] = lon[np.where(lon > 180.)] - 360.
    return lat, lon


def ReadTBB(tbbpath):
    with Dataset(tbbpath, 'r') as fid:
        dataset = fid.variables['NOMChannel12'][:]
    mdataset = ma.masked_where((dataset < 0) | (dataset > 400), dataset)
    return mdataset


def LabelT0(t0data):
    t0data[t0data != -1] = 0
    t0data[t0data == -1] = 1
    t0data_dil = ndimage.binary_dilation(t0data, structure=np.ones((7, 7))).astype(np.int)

    # now identify the objects and remove those above a threshold
    # Generate a structuring element that will consider features connected even if they touch diagonally
    s = ndimage.generate_binary_structure(2, 2)
    Zlabeled, Nlabels = ndimage.measurements.label(t0data_dil, structure=s)
    Zlabeled_v2 = Zlabeled * t0data
    return t0data, Zlabeled_v2


def LabelT1T2(t1data):
    t1data[(t1data == 2) | (t1data == 3) | (t1data == 4) | (t1data == 5) | (
            t1data == 6)] = 1  # 将1:Exctinct;2:preserve;3:deepdevelop，归为1
    t1data[t1data != 1] = 0

    t1data_dil = ndimage.binary_dilation(t1data, structure=np.ones((7, 7))).astype(np.int)

    s = ndimage.generate_binary_structure(2, 2)

    labeled_array, Nlabels = ndimage.measurements.label(t1data_dil, structure=s)
    labeled_array_v2 = labeled_array * t1data

    # m, n = labeled_array.shape
    #
    # for x in range(m - 2):
    #     for y in range(n - 2):
    #         print x, y
    #         scanwindow = labeled_array[x:x + 3, y:y + 3]
    #         scanwindow_ma = ma.masked_values(scanwindow, 0)
    #         minvalue = np.min(scanwindow_ma)
    #         scanwindow[scanwindow > minvalue] = minvalue
    #         labeled_array[x:x + 3, y:y + 3] = scanwindow
    return t1data, labeled_array_v2


def writenetcdf(data, filepath):
    ncfile = Dataset(filepath, mode='w')
    lat = ncfile.createDimension('lat', 1108)  # latitude axis
    lon = ncfile.createDimension('lon', 2748)  # longitude axis
    QREF = ncfile.createVariable('COMP', np.int32, ('lat', 'lon'))
    QREF[:, :] = data
    ncfile.close()


def CheckTB(dataT0, t0Label, dataT1, t1Label, dataT2, t2Label, names):
    TBBt0path = names['tbb_t0path']
    TBBt1path = names['tbb_t1path']
    TBBt2path = names['tbb_t2path']
    deta1 = dataT0 + dataT1 + dataT2
    t0list = np.unique(t0Label[deta1 == 3])  # 匹配上的区域
    t1list = np.unique(t1Label[deta1 == 3])
    t2list = np.unique(t2Label[deta1 == 3])

    t0label_size = [(t0Label == labelt0).sum() for labelt0 in t0list]
    t1label_size = [(t1Label == labelt1).sum() for labelt1 in t1list]
    t2label_size = [(t2Label == labelt2).sum() for labelt2 in t2list]

    dt0 = []
    dt1 = []
    dt2 = []
    t0minlist = []
    t1minlist = []
    t2minlist = []

    TBBt0data = ReadTBB(TBBt0path)
    TBBt1data = ReadTBB(TBBt1path)
    TBBt2data = ReadTBB(TBBt2path)
    for t0 in t0list:
        t0min = np.min(TBBt0data[t0Label == t0])
        t0max = np.max(TBBt0data[t0Label == t0])
        t0minlist.append(t0min)
        dt0.append(np.abs(t0min - t0max))
    for t1 in t1list:
        t1min = np.min(TBBt1data[t1Label == t1])
        t1max = np.max(TBBt1data[t1Label == t1])
        t1minlist.append(t1min)
        dt1.append(np.abs(t1min - t1max))
    for t2 in t2list:
        t2min = np.min(TBBt2data[t2Label == t2])
        t2max = np.max(TBBt2data[t2Label == t2])
        t2minlist.append(t2min)
        dt2.append(np.abs(t2min - t2max))
    indicator = collections.OrderedDict()
    indicator['t0label_size'] = t0label_size
    indicator['t1label_size'] = t1label_size
    indicator['t2label_size'] = t2label_size
    indicator['t0minlist'] = t0minlist
    indicator['t1minlist'] = t1minlist
    indicator['t2minlist'] = t2minlist
    indicator['dt0'] = dt0
    indicator['dt1'] = dt1
    indicator['dt2'] = dt2

    return indicator


# 计算检验指标
def CalMetric(indicator, t0Label):
    total = np.sum(t0Label != 0)

    t0label_size = np.array(indicator['t0label_size'])
    t1label_size = np.array(indicator['t1label_size'])
    t2label_size = np.array(indicator['t2label_size'])
    t0minlist = np.array(indicator['t0minlist'])
    t1minlist = np.array(indicator['t1minlist'])
    t2minlist = np.array(indicator['t2minlist'])
    dt0 = np.array(indicator['dt0'])
    dt1 = np.array(indicator['dt1'])
    dt2 = np.array(indicator['dt2'])

    num = len(t0label_size)
    flag = np.ones(num)
    for i in range(num):
        if ((t2label_size[i] > t1label_size[i] > t0label_size[i]) &  # 条件1：区域尺寸变大
                (t2minlist[i] < t1minlist[i] < t0minlist[i]) &  # 条件2：TB最小值下降
                (dt2[i] > dt1[i] > dt0[i])):  # 条件2：△T=TBmin-TBmax 增大
            flag[i] = 1
        else:
            flag[i] = 0

    Tpixel = t0label_size[flag == 1].sum()
    POD = round(float(Tpixel) / total, 4)
    return POD


# 绘制柱状图
def DrawBar(pod, names):
    barname = names['histname']
    histtitle = names['histtitle']
    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)  # 输出图片的大小
    colors = ["#c72e29", "#016392", "#be9c2e", "#098154"]
    ax.bar(0.5, pod, tick_label="POD", color=colors, width=0.05, zorder=2)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)

    # 循环，为每个柱形添加文本标注
    # for xx, yy in zip(range(len(rates)), rates.values()):
    if pod == np.inf:
        ax.text(0.5, 0.00, "NAN", ha='center')
    elif pod >= 0:
        ax.text(0.5, pod + 0.01, str(pod), ha='center')
    else:
        ax.text(0.5, pod - 0.04, str(pod), ha='center')

    ax.set_ylabel('Rates', fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')
    fig.suptitle(histtitle, fontsize=14, fontweight='bold')

    plt.savefig(barname)
    plt.close
    return True


# 绘制产品分布图
def DrawMap(data, lat, lon, titlename, filename, tbbpath):
    mlat = ma.masked_values(lat, -999.)
    mlon = ma.masked_values(lon, -999.)
    mdata = ma.masked_values(data, 255)
    mmdata = ma.masked_values(mdata, 0)

    fig = plt.figure(figsize=(27.48, 30), dpi=100)  # 图像的长*高=2748*3000
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
    m = Basemap(projection='ortho', lat_0=0, lon_0=104.7, resolution='l', ax=axes1)

    x, y = m(mlon, mlat)
    m.drawmapboundary(color='white')
    m.drawcoastlines(color='white')
    m.drawcountries(color='white')

    parallels = np.arange(-90., 90., 10.)
    m.drawparallels(parallels, color='white')
    meridians = np.arange(-180., 180., 10.)
    m.drawmeridians(meridians, color='white')

    if os.path.exists(tbbpath):
        row = lat.shape[0]
        mtbb = ReadTBB(tbbpath)[:row, :]
        m.contourf(x, y, mtbb, cmap=plt.cm.gray_r)

    colors = ["#ff0000", "#04ac30", "#81b805", "#b7fe9a", "#f7f80f", "#ff7a00", "#448764"]
    # -1:CI;1:Exctinct;2:preserve;3:deepdevelop;4:rdc4;5:rdc5;6:rdc6
    cs = m.contourf(x, y, mmdata, [-1.9, -0.9, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1], colors=colors)

    # 加图例
    axes2.add_patch(plt.Rectangle((0.1, 0.5), 0.04, 0.3, color=colors[0]))  # 左下起点，长，宽，颜色
    axes2.add_patch(plt.Rectangle((0.25, 0.5), 0.04, 0.3, color=colors[1]))
    axes2.add_patch(plt.Rectangle((0.4, 0.5), 0.04, 0.3, color=colors[2]))
    axes2.add_patch(plt.Rectangle((0.55, 0.5), 0.04, 0.3, color=colors[3]))
    axes2.add_patch(plt.Rectangle((0.7, 0.5), 0.04, 0.3, color=colors[4]))
    axes2.add_patch(plt.Rectangle((0.85, 0.5), 0.04, 0.3, color=colors[5]))
    axes2.add_patch(plt.Rectangle((0.1, 0.1), 0.04, 0.3, color=colors[6]))
    axes2.add_patch(plt.Rectangle((0.25, 0.1), 0.04, 0.3, color='white'))

    axes2.text(0.14, 0.55, 'CI', fontsize=30)
    axes2.text(0.29, 0.55, 'Extinct', fontsize=30)
    axes2.text(0.44, 0.55, 'Preserve', fontsize=30)
    axes2.text(0.59, 0.55, 'Deep Develop', fontsize=30)
    axes2.text(0.74, 0.55, 'RDC4', fontsize=30)
    axes2.text(0.89, 0.55, 'RDC5', fontsize=30)
    axes2.text(0.14, 0.15, 'RDC6', fontsize=30)
    axes2.text(0.29, 0.15, 'Space', fontsize=30)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)

    axes1.text(0.01, 0.99, titlename,
               horizontalalignment='left',
               verticalalignment='top',
               transform=axes1.transAxes,
               family='Times New Roman',
               fontsize=42,
               fontweight='bold')

    fig.savefig(filename)
    plt.close()
    return True


# 绘制产品分布图
def DrawMap1(data, lat, lon, titlename, filename, tbbpath):
    mlat = ma.masked_values(lat, -999.)
    mlon = ma.masked_values(lon, -999.)
    mdata = ma.masked_values(data, 255)
    mmdata = ma.masked_values(mdata, 0)

    fig = plt.figure(figsize=(25, 21), dpi=100)
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')

    m = Basemap(projection='cyl', resolution='l', ax=axes1,llcrnrlat=0, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)

    x, y = m(mlon, mlat)
    m.drawmapboundary(color='white')
    m.drawcoastlines(color='white')
    m.drawcountries(color='white')
    parallels = np.arange(-90., 90., 10.)
    m.drawparallels(parallels, color='white',linewidth=1, labels=[1, 0, 0, 0], fontsize=16)
    meridians = np.arange(-180., 180., 10.)
    m.drawmeridians(meridians, color='white',linewidth=1, labels=[1, 0, 0, 0], fontsize=16)
    if os.path.exists(tbbpath):
        row = lat.shape[0]
        mtbb = ReadTBB(tbbpath)[:row, :]
        m.contourf(x, y, mtbb, cmap=plt.cm.gray_r)
    colors = ["#ff0000", "#04ac30", "#81b805", "#b7fe9a", "#f7f80f", "#ff7a00", "#448764"]
    # -1:CI;1:Exctinct;2:preserve;3:deepdevelop;4:rdc4;5:rdc5;6:rdc6
    cs = m.contourf(x, y, mmdata, [-1.9, -0.9, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1], colors=colors)
    # 加图例
    axes2.add_patch(plt.Rectangle((0.1, 0.4), 0.04, 0.2, color=colors[0]))  # 左下起点，长，宽，颜色
    axes2.add_patch(plt.Rectangle((0.25, 0.4), 0.04, 0.2, color=colors[1]))
    axes2.add_patch(plt.Rectangle((0.4, 0.4), 0.04, 0.2, color=colors[2]))
    axes2.add_patch(plt.Rectangle((0.55, 0.4), 0.04, 0.2, color=colors[3]))
    axes2.add_patch(plt.Rectangle((0.7, 0.4), 0.04, 0.2, color=colors[4]))
    axes2.add_patch(plt.Rectangle((0.85, 0.4), 0.04, 0.2, color=colors[5]))
    axes2.add_patch(plt.Rectangle((0.1, 0.1), 0.04, 0.2, color=colors[6]))
    axes2.add_patch(plt.Rectangle((0.25, 0.1), 0.04, 0.2, color='white'))

    axes2.text(0.14, 0.45, 'CI', fontsize=24)
    axes2.text(0.29, 0.45, 'Extinct', fontsize=24)
    axes2.text(0.44, 0.45, 'Preserve', fontsize=24)
    axes2.text(0.59, 0.45, 'Deep Develop', fontsize=24)
    axes2.text(0.74, 0.45, 'RDC4', fontsize=24)
    axes2.text(0.89, 0.45, 'RDC5', fontsize=24)
    axes2.text(0.14, 0.15, 'RDC6', fontsize=24)
    axes2.text(0.29, 0.15, 'Space', fontsize=24)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)

    axes2.text(0.35, 0.7, unicode(titlename, "utf8"), fontproperties='SimHei', fontsize=30)

    fig.savefig(filename)
    plt.close()
    return 0


def LogXML(logMsg, logType, xmlfilePath, inputFilename, fyTime=True, metric=True, names=True):
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
        nodeOutputFilename0.appendChild(doc.createTextNode(str(os.path.basename(names['t0map']))))
        nodeOutputFile.appendChild(nodeOutputFilename0)

        nodeOutputFilename1 = doc.createElement("OutputFilename")
        nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['t1map']))))
        nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeOutputFilename2 = doc.createElement("OutputFilename")
        nodeOutputFilename2.appendChild(doc.createTextNode(str(os.path.basename(names['t2map']))))
        nodeOutputFile.appendChild(nodeOutputFilename2)

        nodeOutputFilename3 = doc.createElement("OutputFilename")
        nodeOutputFilename3.appendChild(doc.createTextNode(str(os.path.basename(names['histname']))))
        nodeOutputFile.appendChild(nodeOutputFilename3)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', 'FY4A')
        nodeDeviationstatic.setAttribute('checkVariety', "FY4A")
        nodeDeviationstatic.setAttribute('dataFormat', '15mm')
        nodeDeviationstatic.setAttribute('dataType', 'NOM')
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', 'CIX')
        metricStr = "{POD:%s}" % re.sub('[()\'\[\]]', '',
                                        str(metric).replace("',", ":").replace("OrderedDict", "")).replace(' ', '')
        nodeDeviationstatic.setAttribute('values', metricStr)
        nodeOutputFile.appendChild(nodeDeviationstatic)

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