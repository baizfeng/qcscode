#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/6/21 15:06
# @Author  : Bai
# @File    : Func.py


from xml.dom import minidom
import os
import datetime
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from scipy import spatial
import matplotlib

matplotlib.use("Pdf")
import matplotlib.pyplot as plt
import collections
from mpl_toolkits.basemap import Basemap, interp
import struct
from matplotlib.colors import from_levels_and_colors
from scipy import ndimage
import re
import glob
from matplotlib.font_manager import FontProperties


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

        externode = root.getElementsByTagName('exterDataPath')[0]
        xmlnodes['exterpath'] = externode.childNodes[0].nodeValue

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
    fy4_time = datetime.datetime.strptime(filename.split('_')[9], '%Y%m%d%H%M%S')  # 时间
    dataName = filename[30:33]  # 数据名称
    dataSour = filename[0:4]  # 数据源名称
    dataRes = filename[74:79]
    baseInfo = {"dataName": dataName,
                "dataSour": dataSour,
                "dataRes": dataRes,
                "time": fy4_time,
                "filename": filename[:-3]}
    return baseInfo


def Name(infos, nodes, times):
    filename = infos['filename']
    time0 = infos['time']
    prjType = nodes['prjType']
    outputpath = nodes['outputpath']
    auxpath = nodes['auxpath']

    if (prjType == 'NOM'):
        prj = "MAPN"
    elif (prjType == "GLL"):
        prj = "MAPG"

    histname = "%s%s_QCS_RADA_HIST.png" % (outputpath, filename)  # 直方图名称
    histtitle = "FY4A_RADA_CIX %s" % (time0.__format__('%Y-%m-%d %H:%M:%S'))  # 直方图标题

    radamap = "%s%s_QCS_RADA_COMPARE_RADA.png" % (outputpath, filename)  # 插值后雷达拼图名称
    radatitle = "全国雷达拼图\n%s-%s UTC" % (times[0].__format__('%Y-%m-%d %H:%M:%S'),
                                       times[len(times) - 1].__format__('%H:%M:%S'))  # 插值后雷达拼图标题

    biasmap = "%s%s_QCS_RADA_COMPARE_CHECK.png" % (outputpath, filename)  # 偏差分布图名称
    biastitle = "FY4A CIX与雷达拼图对比分布图\nFY4A Time:%s UTC RADA Time:%s-%s UTC" % \
                (time0.__format__('%Y-%m-%d %H:%M:%S'), times[0].__format__('%Y-%m-%d %H:%M:%S'),
                 times[len(times) - 1].__format__('%H:%M:%S'))  # 偏差分布图标题

    fy4map = "%s%s_QCS_RADA_COMPARE_FY4A.png" % (outputpath, filename)
    fy4title = "FY4A_CIX_%s %s" % (prj, time0.__format__('%Y-%m-%d %H:%M:%S'))  # 标题

    tbbpath = "%s/%s/%s/%s.NC" % (
        auxpath, time0.__format__('%Y'), time0.__format__('%Y%m%d'), filename.replace('CIX', 'TBB'))

    names = {"histname": histname,
             "histtitle": histtitle,
             "radamap": radamap,
             "radatitle": radatitle,
             "biasmap": biasmap,
             "biastitle": biastitle,
             "fy4map": fy4map,
             "fy4title": fy4title,
             "tbbpath": tbbpath}
    return names


# 时间匹配：
# path-待匹配数据路径
# base_time-匹配的时间
# timeThreshold-时间阈值，FY4时间向前找雷达时间半个小时
def TimeMatch(exterpath, base_time, timeThreshold):
    # files = os.listdir(exterpath)
    files = glob.glob(os.path.join(exterpath, 'Z_RADA_C_BABJ*.latlon'))
    dfiles = []
    dtimes = []
    for file in files:
        file = os.path.basename(file)
        timeStr = file.split('_')[4]
        time = datetime.datetime.strptime(timeStr, '%Y%m%d%H%M%S')
        # utctime = time - datetime.timedelta(hours=8)
        if ((base_time - time).total_seconds() / 60 <= float(timeThreshold)) & (
                (base_time - time).total_seconds() / 60 >= 0):
            dfiles.append(file)
            dtimes.append(time)
    if (len(dfiles) == 0):
        raise IndexError("Cannot find rada data in the special time threshold!")
    return dfiles, dtimes


# 利用最大值法，处理匹配到的多条雷达数据
def MaxRada(filenames, exterpath):
    radas = []
    # 读取数据
    for filename in filenames:
        filename = "%s%s" % (exterpath, filename)
        rada = ReadRADA(filename)[0]
        radas.append(rada.flatten())
    radas = np.array(radas).T

    rada_file = "%s%s" % (exterpath, filenames[0])
    lat = ReadRADA(rada_file)[1]
    lon = ReadRADA(rada_file)[2]

    radaMax = np.max(radas, axis=1).reshape(rada.shape)
    return radaMax, lat, lon


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


def ReadRADA(radapath):
    f = open(radapath, 'rb').read()
    # 经纬度信息（南纬、西经，北纬，东经）
    slat = list(struct.unpack('f', f[180:184]))[0]
    wlon = list(struct.unpack('f', f[184:188]))[0]
    nlat = list(struct.unpack('f', f[188:192]))[0]
    elon = list(struct.unpack('f', f[192:196]))[0]
    # 行列数
    rows = list(struct.unpack('i', f[204:208]))[0]
    cols = list(struct.unpack('i', f[208:212]))[0]
    # 单层数据字节数
    levelbytes = list(struct.unpack('i', f[224:228]))[0]
    # 数值放大倍数
    amp = list(struct.unpack('h', f[230:232]))[0]
    # 根据经纬度生成初始矩阵，及经纬度栅格资料
    lon_g1km = np.linspace(wlon, elon, cols)
    lat_g1km = np.linspace(nlat, slat, rows)
    # lon2d, lat2d = np.meshgrid(lon_g1km, lat_g1km)
    dbzdata = np.zeros((len(lat_g1km), len(lon_g1km)), np.float32)

    # 数据区
    rowcount = 256
    bytenum = 0
    n = 0
    # 以每层所占最大字符数为界，如果行数&列数&连续字段数同时输出为-1，则结束循环。
    while (bytenum <= levelbytes):
        # 数据开始行，列，连续字段数据量信息
        datainfor = np.array(struct.unpack('3h', f[rowcount:rowcount + 6]))
        if (datainfor[0] == -1 and datainfor[1] == -1 and datainfor[2] == -1):
            break
        bytestr = str(datainfor[2]) + 'h'
        dbzdata[datainfor[0], datainfor[1]:datainfor[1] + datainfor[2]] \
            = np.array(struct.unpack(bytestr, f[rowcount + 6:rowcount + 6 + datainfor[2] * 2]))
        rowcount = rowcount + 6 + datainfor[2] * 2
        bytenum = bytenum + 6 + datainfor[2] * 2
        n = n + 1
    dbzdata = dbzdata / float(amp)
    dbzdata[(dbzdata < 0) | (dbzdata > 400)] = 0
    mdbzdata = ma.masked_values(dbzdata, 0)
    return mdbzdata, lat_g1km, lon_g1km


def ReadTBB(tbbpath):
    with Dataset(tbbpath, 'r') as fid:
        dataset = fid.variables['NOMChannel12'][:]
    mdataset = ma.masked_where((dataset < 0) | (dataset > 400), dataset)
    return mdataset


# 最邻近插值函数：
# 参数说明：tlat，tlon为目标数据经纬度
# slat,slon，sdata为待插值经纬度和数据
# 注意：经纬度范围要保持一致
def KDTreeSample(tlat, tlon, slat, slon, sdata, spaceThreshold):
    tree = spatial.cKDTree(zip(slat.flatten(), slon.flatten()))  # 先将待插值数据建树，latlon不能有NA值
    d, inds = tree.query(zip(tlat.flatten(), tlon.flatten()), k=1, distance_upper_bound=float(spaceThreshold))
    nearest = -999 * np.ones(tlat.size)  # 阈值之外的赋值为-999
    nearest[~np.isinf(d)] = sdata.flatten()[inds[~np.isinf(d)]]
    return nearest.reshape(tlat.shape)


def BasemapInterp(tlat, tlon, slat, slon, sdata):
    lat_new = np.flipud(slat)
    sdata_new = np.flipud(sdata)

    rada_nom = interp(sdata_new, slon, lat_new, tlon, tlat, checkbounds=False, masked=-999., order=1)
    return rada_nom


# 剔除雷达数据中的小斑点
def postprocess(data):
    data[data < 35.] = 0
    data[data >= 35.] = 1

    # now identify the objects and remove those above a threshold
    # Generate a structuring element that will consider features connected even if they touch diagonally
    s = ndimage.generate_binary_structure(2, 2)
    Zlabeled, Nlabels = ndimage.measurements.label(data, structure=s)

    label_size = [(Zlabeled == label).sum() for label in range(Nlabels + 1)]
    # for label, size in enumerate(label_size): print("label %s is %s pixels in size" % (label, size))

    # # now remove the labels
    for label, size in enumerate(label_size):
        if size < 9:
            data[Zlabeled == label] = 0
    return data


# 计算检验指标
def CalMetric(data1, data2):
    comp = -999 * np.ones(data1.shape)
    # metric = collections.OrderedDict()

    TP = np.where((data1 == -1) & (data2 == 1))
    FP = np.where((data1 != -1) & (data2 == 1))
    FN = np.where((data1 == -1) & (data2 == 0))
    TN = np.where((data1 != -1) & (data2 == 0))

    comp[TP] = 1
    comp[FP] = 2
    comp[FN] = 3
    comp[TN] = 4

    TP_L = len(TP[0])
    FP_L = len(FP[0])
    FN_L = len(FN[0])
    # TN_L = len(TN[0])

    # metric["POD"] = round(float(TP_L) / (TP_L + FN_L), 4) if (TP_L + FN_L) != 0 else np.inf
    # metric["FAR"] = round(float(FN_L) / (FN_L + TP_L), 4) if (FN_L + TP_L) != 0 else np.inf
    # metric["CSI"] = round(float(TP_L) / (TP_L + FP_L + FN_L), 4) if (TP_L + FP_L + FN_L) != 0 else np.inf
    pod = round(float(TP_L) / (TP_L + FN_L), 4) if (TP_L + FN_L) != 0 else np.inf
    return pod, comp


# 绘制柱状图
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


def DrawNOMMap(data, titlename, filename, fy4_lat, fy4_lon):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    mlat = ma.masked_values(fy4_lat, -999.)
    mlon = ma.masked_values(fy4_lon, -999.)
    mdata = ma.masked_values(data, -999.)

    fig = plt.figure(figsize=(25, 21), dpi=100)
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
    axes1 = fig.add_axes([0.03, 0.11, 0.95, 0.88])
    cax = fig.add_axes([0.15, 0.03, 0.7, 0.02])

    m = Basemap(projection='cyl', resolution='l', ax=axes1,llcrnrlat=0, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)

    x, y = m(mlon, mlat)
    m.drawcoastlines()
    m.drawcountries(linewidth=1.)

    m.drawparallels(np.arange(-90., 90., 10.),linewidth=1, labels=[1, 0, 0, 0], fontsize=16)
    m.drawmeridians(np.arange(-180., 180., 10.),linewidth=1, labels=[0, 0, 0, 1], fontsize=16)

    v = np.linspace(10., 75.0, 14, endpoint=True)
    colors = ["#01a0f6", "#00ecec", "#00d800", "#019000",
              "#ffff00", "#e7c000", "#ff9000", "#ff0000",
              "#d60000", "#c00000", "#ff00f0", "#9600b4", "#ad90f0"]
    cs = m.contourf(x, y, mdata, v, colors=colors)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)
    axes2.set_xticks([])
    axes2.set_yticks([])
    #axes1.spines['left'].set_visible(False)
    #axes1.spines['top'].set_visible(False)
    #axes1.set_xticks([])
    #axes1.set_yticks([])

    axes1.text(0.01, 0.99, unicode(titlename, "utf8"), fontproperties=font, fontsize=36,
               horizontalalignment='left',
               verticalalignment='top',
               transform=axes1.transAxes)

    # 加图例
    cb = plt.colorbar(cs, ticks=v, orientation='horizontal', cax=cax)
    cb.ax.set_title(unicode("组合反射率(单位：dBZ)", "utf8"), fontproperties=font, fontsize=24)
    cb.ax.tick_params(labelsize=20)

    fig.savefig(filename)
    plt.close()
    return True


def DrawRADAMap(data, lon, lat, titlename, filename):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    lon2d, lat2d = np.meshgrid(lon, lat)

    fig = plt.figure(figsize=(12, 13))
    axes1 = fig.add_axes([0.03, 0.03, 0.95, 0.95])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高
    axes2 = fig.add_axes([0.048, 0.93, 0.4, 0.05])
    axes3 = fig.add_axes([0.048, 0.03, 0.09, 0.3])
    cax = fig.add_axes([0.075, 0.04, 0.03, 0.25])

    m = Basemap(projection='poly',
                lat_0=35, lon_0=110, llcrnrlon=80, llcrnrlat=3.01, urcrnrlon=140, urcrnrlat=53.123,
                resolution='l', ax=axes1)

    x, y = m(lon2d, lat2d)
    m.drawcoastlines()
    m.drawcountries(linewidth=1.)

    parallels = np.arange(-90., 90., 10.)
    m.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=10)
    meridians = np.arange(-180., 180., 10.)
    m.drawmeridians(meridians, labels=[0, 0, 0, 1], fontsize=10)

    v = np.linspace(10., 75.0, 14, endpoint=True)
    colors = ["#01a0f6", "#00ecec", "#00d800", "#019000",
              "#ffff00", "#e7c000", "#ff9000", "#ff0000",
              "#d60000", "#c00000", "#ff00f0", "#9600b4", "#ad90f0"]
    cs = m.contourf(x, y, data, v, colors=colors)

    axes2.set_xticks([])
    axes2.set_yticks([])
    axes3.set_xticks([])
    axes3.set_yticks([])
    axes2.text(0.05, 0.1, unicode(titlename, "utf8"), fontproperties=font, fontsize=20)
    cb = plt.colorbar(cs, ticks=v, cax=cax)
    cb.ax.set_title(unicode("组合反射率\ndBZ", "utf8"), fontproperties=font, fontsize=14)

    fig.savefig(filename, dpi=200)
    plt.close()
    return True


def DrawCompMap(data, titlename, filename, fy4_lat, fy4_lon):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    fig = plt.figure(figsize=(25, 21), dpi=100)
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')

    m = Basemap(projection='cyl', resolution='l', ax=axes1,llcrnrlat=0, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)

    x, y = m(fy4_lon, fy4_lat)
    m.drawcoastlines()
    m.drawcountries()

    m.drawparallels(np.arange(-90., 90., 10.),linewidth=1, labels=[1, 0, 0, 0], fontsize=16)
    m.drawmeridians(np.arange(-180., 180., 10.),linewidth=1, labels=[0, 0, 0, 1], fontsize=16)

    colors = ["#ff7473", "#ffc952", "#47b8e0", "white"]
    # cmap, norm = from_levels_and_colors([0, 1.1, 2.1, 3.1, 4.1], colors)

    m.contourf(x, y, data, [0, 1.1, 2.1, 3.1, 4.1], colors=colors)

    # 加图例
    axes2.add_patch(plt.Rectangle((0.1, 0.2), 0.04, 0.2, color=colors[0]))  # 左下起点，长，宽，颜色
    axes2.add_patch(plt.Rectangle((0.25, 0.2), 0.04, 0.2, color=colors[1]))
    axes2.add_patch(plt.Rectangle((0.4, 0.2), 0.04, 0.2, color=colors[2]))
    axes2.add_patch(plt.Rectangle((0.55, 0.2), 0.04, 0.2, color=colors[3]))

    axes2.text(0.14, 0.25, 'Hit', fontsize=24)
    axes2.text(0.29, 0.25, 'Miss', fontsize=24)
    axes2.text(0.44, 0.25, 'False alarm', fontsize=24)
    axes2.text(0.59, 0.25, 'Correct negative', fontsize=24)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)
    axes1.spines['left'].set_visible(False)
    axes1.spines['top'].set_visible(False)
    axes1.set_xticks([])
    axes1.set_yticks([])

    axes2.text(0.1, 0.5, unicode(titlename, "utf8"), fontproperties=font, fontsize=30)

    fig.savefig(filename)
    plt.close()
    return True


# 绘制产品分布图
def DrawMap(data, lat, lon, titlename, filename, tbbpath):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    mlat = ma.masked_values(lat, -999.)
    mlon = ma.masked_values(lon, -999.)
    mdata = ma.masked_values(data, 255)
    mmdata = ma.masked_values(mdata, 0)

    fig = plt.figure(figsize=(25, 21), dpi=100)
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')

    m = Basemap(projection='ortho', lat_0=30, lon_0=104.7, resolution='l',
                llcrnrx=-3300 * 1000, llcrnry=-2000 * 1000,
                urcrnrx=+3200 * 1000, urcrnry=+3000 * 1000,
                ax=axes1)

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
    cs = m.contourf(x, y, mmdata, [-1.9, -0.9, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1], colors=colors, V=1)

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

    axes2.text(0.35, 0.7, unicode(titlename, "utf8"), fontproperties=font, fontsize=30)

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
        nodeOutputFilename0.appendChild(doc.createTextNode(str(os.path.basename(names['fy4map']))))
        nodeOutputFile.appendChild(nodeOutputFilename0)

        nodeOutputFilename1 = doc.createElement("OutputFilename")
        nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['radamap']))))
        nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeOutputFilename2 = doc.createElement("OutputFilename")
        nodeOutputFilename2.appendChild(doc.createTextNode(str(os.path.basename(names['biasmap']))))
        nodeOutputFile.appendChild(nodeOutputFilename2)

        nodeOutputFilename3 = doc.createElement("OutputFilename")
        nodeOutputFilename3.appendChild(doc.createTextNode(str(os.path.basename(names['histname']))))
        nodeOutputFile.appendChild(nodeOutputFilename3)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', 'FY4A')
        nodeDeviationstatic.setAttribute('checkVariety', "RADA")
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