#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/10/23 17:45
# @Author  : baizhaofeng
# @File    : Common_NA.py


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
from matplotlib.colors import LogNorm
from scipy import spatial, stats
import re
import matplotlib.mlab as mlab
import seaborn as sns
from mpl_toolkits.basemap import Basemap
import glob


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

        externode = root.getElementsByTagName('exterDataPath')[0]
        xmlnodes['exterpath'] = externode.childNodes[0].nodeValue

        outputnode = root.getElementsByTagName('outputPath')[0]
        xmlnodes['outputpath'] = outputnode.childNodes[0].nodeValue

        outputxmlnode = root.getElementsByTagName('outputXml')[0]
        xmlnodes['xmlPath'] = outputxmlnode.childNodes[0].nodeValue

        prjTypenode = root.getElementsByTagName('type')[0]
        xmlnodes['prjType'] = prjTypenode.childNodes[0].nodeValue

        chVarnode = root.getElementsByTagName('checkVariety')[0]
        xmlnodes['cVar'] = chVarnode.childNodes[0].nodeValue

        timenode = root.getElementsByTagName('timeValue')[0]
        xmlnodes['timevalue'] = timenode.childNodes[0].nodeValue

        spacenode = root.getElementsByTagName('spaceValue')[0]
        xmlnodes['spacevalue'] = spacenode.childNodes[0].nodeValue

        warnnode = root.getElementsByTagName('warningValue')[0]
        xmlnodes['warnvalue'] = warnnode.childNodes[0].nodeValue

        nanode = root.getElementsByTagName('abnormal')[0]
        xmlnodes['na'] = nanode.childNodes[0].nodeValue

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


def TimeMatch(path, base_time, timeThreshold):
    """
    时间匹配：根据FY4时间得到检验数据
    :param path:
    :param base_time:
    :return:
    """
    files = glob.glob(os.path.join(path, 'AHI8-_AGRI*N_DISK*SST*NOM*V0001.NC'))  # AHI8-_AGRI*N_DISK*FHS*NOM*V0001.NC
    dfiles = []
    dtimes = []
    for filepath in files:
        file = os.path.basename(filepath)
        timeStr = file.split('_')[9]
        time = datetime.datetime.strptime(timeStr, '%Y%m%d%H%M%S')
        if (np.abs((time - base_time).total_seconds() / 60) <= float(timeThreshold)):
            dfiles.append(filepath)
            dtimes.append(time)
    dfiles.sort()
    dtimes.sort()
    if (len(dfiles) == 0):
        raise NameError("未找到检验源数据")
    return dfiles, dtimes


def Name(filename, fytime, outputpath, prjType, checkVar):
    if (prjType == 'NOM'):
        prj = "MAPN"
    elif (prjType == "GLL"):
        prj = "MAPG"

    scattername = "%s%s_QCS_%s_SCAT_ORIG_NA.png" % (outputpath, filename, checkVar)
    histname = "%s%s_QCS_%s_HIST_BIAS_NA.png" % (outputpath, filename, checkVar)
    fy4name = "%s%s_QCS_%s_ORIG_NA.png" % (outputpath, filename, prj)
    checkname = "%s%s_QCS_%s_%s_ORIG_NA.png" % (outputpath, filename, checkVar, prj)
    biasname = "%s%s_QCS_%s_%s_BIAS_NA.png" % (outputpath, filename, checkVar, prj)
    scattertitle = "FY4A_%s_SST %s" % (checkVar, fytime.__format__('%Y-%m-%d %H:%M:%S'))
    histtitle = "FY4A-%s_SST %s" % (checkVar, fytime.__format__('%Y-%m-%d %H:%M:%S'))
    fy4title = "FY4A_SST_%s %s" % (prj, fytime.__format__('%Y-%m-%d %H:%M:%S'))
    biastitle = "FY4A-%s_SST_%s %s" % (checkVar, prj, fytime.__format__('%Y-%m-%d %H:%M:%S'))  # 标题
    checktitle = "%s_SST_%s %s" % (checkVar, prj, fytime.__format__('%Y-%m-%d %H:%M:%S'))  # 标题

    names = {'scatname': scattername,
             'histname': histname,
             'biasname': biasname,
             'fy4name': fy4name,
             'checkname': checkname,
             'scattitle': scattertitle,
             'histtitle': histtitle,
             'fy4title': fy4title,
             'checktitle': checktitle,
             'biastitle': biastitle}
    return names


def Read_Data(filepath):
    """
    读取数据，根据变量名
    :param filepath:
    :return: SST数据,质量标识[0:good,1:mid,2:bad]，array类型
    """
    # 读取FY4数据
    with Dataset(filepath, 'r') as fid:
        fy4_data = fid.variables["SST"][:]
        dqf = fid.variables["DQF"][:]
    if (isinstance(fy4_data, ma.MaskedArray)):
        fy4_data = fy4_data.data
    return fy4_data, dqf


def Read_exterData(filepath):
    """
    读取数据，根据变量名
    :param filepath:
    :return: SST数据,质量标识[0:good,1:mid,2:bad]，array类型
    """
    # 读取FY4数据
    with Dataset(filepath, 'r') as fid:
        fy4_data = fid.variables["SST"][:]
    if (isinstance(fy4_data, ma.MaskedArray)):
        fy4_data = fy4_data.data
    return fy4_data


def Read_LUT(prjType, filepath, dataRes):
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
        elif (dataRes == '2000M'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/FullMask_Grid_2000_1047E_BSQ.nc'
        elif (dataRes == '012KM'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/FullMask_Grid_12000_1047E_BSQ.nc'

        with Dataset(lutpath, 'r') as fid:
            lat = fid.variables['LAT'][:]
            lon = fid.variables['LON'][:]
        lon[np.where(lon > 180.)] = lon[np.where(lon > 180.)] - 360.
    return lat, lon


def Read_cLUT():
    """
    根据数据的投影类型，分辨率读取经纬度查找表
    :param prjType:
    :param filepath:
    :param dataRes:
    :return:
    """
    lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/AHI8_OBI_4000M_NOM_LATLON.HDF'

    with Dataset(lutpath, 'r') as fid:
        lat = fid.variables['Lat'][:]
        lon = fid.variables['Lon'][:]
    return lat, lon


def lon_lat_to_cartesian(lon, lat):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    R = 6367  # radius of the Earth in kilometers
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x = R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x, y, z


def Interp(tlat, tlon, slat, slon, sdata, spaceThreshold):
    # 第一步：标记NAN值,目标经纬度不需要标记为NA
    slat[np.where(slat == 65534.)] = np.nan  # 注意小数点
    slon[np.where(slon == 65534.)] = np.nan

    # 2.转换经纬度到笛卡尔坐标系
    xs, ys, zs = lon_lat_to_cartesian(slon[~np.isnan(slon)], slat[~np.isnan(slat)])
    xt, yt, zt = lon_lat_to_cartesian(tlon.flatten(), tlat.flatten())

    tree = spatial.cKDTree(zip(xs, ys, zs))  # 先将待插值数据建树，latlon不能有NA值
    d, inds = tree.query(zip(xt, yt, zt), k=1,distance_upper_bound=float(spaceThreshold))  # 阈值为0.04°
    nearest = -999 * np.ones(tlat.size)  # 阈值之外的赋值为-999
    nearest[~np.isinf(d)] = sdata[~np.isnan(slat)][inds[~np.isinf(d)]]
    return nearest.reshape(tlat.shape)


def CalcateDQF(fy4data, exdata, bias, dqf):
    fy4data0 = fy4data[np.where(dqf == 0)]
    exdata0 = exdata[np.where(dqf == 0)]
    bias0 = bias[np.where(dqf == 0)]
    fy4data0 = fy4data0[~fy4data0.mask]
    exdata0 = exdata0[~exdata0.mask]
    bias0 = bias0[~bias0.mask]

    fy4data1 = fy4data[np.where(dqf == 1)]
    exdata1 = exdata[np.where(dqf == 1)]
    bias1 = bias[np.where(dqf == 1)]
    fy4data1 = fy4data1[~fy4data1.mask]
    exdata1 = exdata1[~exdata1.mask]
    bias1 = bias1[~bias1.mask]

    fy4data2 = fy4data[np.where(dqf == 2)]
    exdata2 = exdata[np.where(dqf == 2)]
    bias2 = bias[np.where(dqf == 2)]
    fy4data2 = fy4data2[~fy4data2.mask]
    exdata2 = exdata2[~exdata2.mask]
    bias2 = bias2[~bias2.mask]

    fy4data3 = fy4data[~fy4data.mask]
    exdata3 = exdata[~exdata.mask]
    bias3 = bias[~bias.mask]
    # 指标计算
    slope0, intercept0, r_value0, p_value0, std_err0 = stats.linregress(fy4data0, exdata0)
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(fy4data1, exdata1)
    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(fy4data2, exdata2)
    slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(fy4data3, exdata3)

    metrics = collections.OrderedDict()
    metrics['QualID0_NUM'] = fy4data0.size  # 总数
    metrics['QualID0_MAX'] = round(bias0.max(), 4)  # 最大值
    metrics['QualID0_MIN'] = round(bias0.min(), 4)  # 最小值
    metrics['QualID0_MEDIAN'] = round(np.median(bias0), 4)  # 中位数
    metrics['QualID0_MEAN'] = round(bias0.mean(), 4)  # 平均值
    metrics['QualID0_AE'] = round(np.abs(bias0).mean(), 4)  # 绝对值平均数
    metrics['QualID0_STD'] = round(np.sqrt(np.square(bias0).sum() / (fy4data0.size - 1)), 4)
    metrics['QualID0_RMSE'] = round(np.sqrt(np.square(bias0).sum() / fy4data0.size), 4)  # 均方根误差
    metrics['QualID0_SKEW'] = round(stats.skew(bias0), 4)  # 偏度系数
    metrics['QualID0_KURT'] = round(stats.kurtosis(bias0), 4)  # 峰度系数
    metrics['QualID0_CORR'] = round(r_value0, 4)

    metrics['QualID0_slope'] = round(slope0, 4)
    metrics['QualID0_intercept'] = round(intercept0, 4)

    metrics['QualID1_NUM'] = fy4data1.size  # 总数
    metrics['QualID1_MAX'] = round(bias1.max(), 4)  # 最大值
    metrics['QualID1_MIN'] = round(bias1.min(), 4)  # 最小值
    metrics['QualID1_MEDIAN'] = round(np.median(bias1), 4)  # 中位数
    metrics['QualID1_MEAN'] = round(bias1.mean(), 4)  # 平均值
    metrics['QualID1_AE'] = round(np.abs(bias1).mean(), 4)  # 绝对值平均数
    metrics['QualID1_STD'] = round(np.sqrt(np.square(bias1).sum() / (fy4data1.size - 1)), 4)
    metrics['QualID1_RMSE'] = round(np.sqrt(np.square(bias1).sum() / fy4data1.size), 4)  # 均方根误差
    metrics['QualID1_SKEW'] = round(stats.skew(bias1), 4)  # 偏度系数
    metrics['QualID1_KURT'] = round(stats.kurtosis(bias1), 4)  # 峰度系数
    metrics['QualID1_CORR'] = round(r_value1, 4)

    metrics['QualID1_slope'] = round(slope1, 4)
    metrics['QualID1_intercept'] = round(intercept1, 4)

    metrics['QualID2_NUM'] = fy4data2.size  # 总数
    metrics['QualID2_MAX'] = round(bias2.max(), 4)  # 最大值
    metrics['QualID2_MIN'] = round(bias2.min(), 4)  # 最小值
    metrics['QualID2_MEDIAN'] = round(np.median(bias2), 4)  # 中位数
    metrics['QualID2_MEAN'] = round(bias2.mean(), 4)  # 平均值
    metrics['QualID2_AE'] = round(np.abs(bias2).mean(), 4)  # 绝对值平均数
    metrics['QualID2_STD'] = round(np.sqrt(np.square(bias2).sum() / (fy4data2.size - 1)), 4)
    metrics['QualID2_RMSE'] = round(np.sqrt(np.square(bias2).sum() / fy4data2.size), 4)  # 均方根误差
    metrics['QualID2_SKEW'] = round(stats.skew(bias2), 4)  # 偏度系数
    metrics['QualID2_KURT'] = round(stats.kurtosis(bias2), 4)  # 峰度系数
    metrics['QualID2_CORR'] = round(r_value2, 4)

    metrics['QualID2_slope'] = round(slope2, 4)
    metrics['QualID2_intercept'] = round(intercept2, 4)

    metrics['QualID3_NUM'] = fy4data3.size  # 总数
    metrics['QualID3_MAX'] = round(bias3.max(), 4)  # 最大值
    metrics['QualID3_MIN'] = round(bias3.min(), 4)  # 最小值
    metrics['QualID3_MEDIAN'] = round(np.median(bias3), 4)  # 中位数
    metrics['QualID3_MEAN'] = round(bias3.mean(), 4)  # 平均值
    metrics['QualID3_AE'] = round(np.abs(bias3).mean(), 4)  # 绝对值平均数
    metrics['QualID3_STD'] = round(np.sqrt(np.square(bias3).sum() / (fy4data3.size - 1)), 4)
    metrics['QualID3_RMSE'] = round(np.sqrt(np.square(bias3).sum() / fy4data3.size), 4)  # 均方根误差
    metrics['QualID3_SKEW'] = round(stats.skew(bias3), 4)  # 偏度系数
    metrics['QualID3_KURT'] = round(stats.kurtosis(bias3), 4)  # 峰度系数
    metrics['QualID3_CORR'] = round(r_value3, 4)

    metrics['QualID3_slope'] = round(slope3, 4)
    metrics['QualID3_intercept'] = round(intercept3, 4)
    return metrics


# 绘制散点图
def DrawScatter(data1, data2, metric, filename, scattitle, checkVar):
    data1 = data1[~data1.mask]
    data2 = data2[~data2.mask]

    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)
    fig.suptitle(scattitle, fontsize=18, fontweight='bold')
    H = ax.hist2d(data1, data2, bins=100, norm=LogNorm(), cmap=plt.get_cmap('jet'), zorder=2)
    ax.set_xlabel('FY4A SST'+'(' + u'\u2103' + ')', fontsize=14)
    ax.set_ylabel(checkVar+" SST"+'(' + u'\u2103' + ')' , fontsize=14)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')

    ax.plot([-5, 45], [-5, 45], color='black', zorder=3)  # 画参考线
    plt.xlim(-5, 45)
    plt.ylim(-5, 45)
    fig.colorbar(H[3], ax=ax)

    txtstr = "%s%s\n%s%s\n%s%s%s%s%s" % (
        'N=', metric['QualID3_NUM'], 'CORR=', metric['QualID3_CORR'], 'y=', metric['QualID3_intercept'], '+',
        metric['QualID3_slope'], 'x')
    ax.text(0.04, 1., txtstr,
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=10)

    x = np.linspace(np.min([data1, data2]), np.max([data1, data2]))
    y = metric['QualID3_intercept'] + metric['QualID3_slope'] * x
    ax.plot(x, y, color='r', zorder=4)  # 画参考线
    plt.savefig(filename)
    plt.close
    return 0


# 绘制偏差直方图
# 直方图类似于柱状图，是用柱的高度来指代频数，不同的是其将定量数据划分为若干连续的区间，在这些连续的区间上绘制柱。
def DrawHist(bias, metric, filename, histtitle, checkVar):
    # 创建直方图
    # 第一个参数为待绘制的定量数据
    # 第二个参数为划分的区间个数
    bias = bias[~bias.mask]
    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)  # 输出图片的大小
    fig.suptitle(histtitle, fontsize=14, fontweight='bold')
    # n, bins, patches = ax.hist(bias, bins=80, density=True, edgecolor='k', facecolor='white', zorder=2)  # bins是分为20个柱子
    # # add a 'best fit' line
    # y = mlab.normpdf(bins, metric['QualID3_MEAN'], metric['QualID3_STD'])
    # l = ax.plot(bins, y, '--', color='black', linewidth=1)

    # kde=True表示是否显示拟合曲线，如果为False则只出现直方图,kde_kws（拟合曲线的设置）、hist_kws（直方柱子的设置）
    sns.distplot(bias, bins=80, hist=True, kde=True,
                 kde_kws={"color": "black", "lw": 1, 'linestyle': '--'},
                 hist_kws={"histtype": 'bar', "edgecolor": 'black', "facecolor": 'white', "lw": 1})

    plt.xlim(-8, 8)
    ax.set_xlabel("FY4A-" + checkVar + " SST" + '(' + u'\u2103' + ')', fontsize=14)
    ax.set_ylabel('Number Density', fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')

    txtstr = "%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s\n%s%s" % (
        'N=', metric['QualID3_NUM'], 'MIN=', metric['QualID3_MIN'], 'MAX=', metric['QualID3_MAX'], 'MEAN=',
        metric['QualID3_MEAN'], 'AE=', metric['QualID3_AE'],
        'MEDIAN=', metric['QualID3_MEDIAN'], 'RMSE=', metric['QualID3_RMSE'],
        'SKEW=', metric['QualID3_SKEW'], 'KURT=', metric['QualID3_KURT'])

    ax.text(0.04, 1., txtstr,
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=10)

    plt.savefig(filename)
    plt.close
    return 0


def DrawNOMMap(data, lat, lon, titlename, filename, level=np.arange(-5, 36, 1)):
    fig = plt.figure(figsize=(27.48, 30), dpi=100)  # 图像的长*高=2748*3000
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
    cax = fig.add_axes([0.1, 0.03, 0.65, 0.02])

    mlat = ma.masked_values(lat, -999.)
    mlon = ma.masked_values(lon, -999.)
    m = Basemap(projection='nsper', lat_0=0, lon_0=104.7, resolution='l',
                ax=axes1)  # resoluiton:c (crude), l (low), i (intermediate), h (high), f (full)
    x, y = m(mlon, mlat)
    m.drawcoastlines()
    m.drawcountries()
    # m.drawlsmask(land_color='#bebebe', ocean_color='#01008a')  # 陆地，海洋的颜色
    m.drawparallels(range(-90, 90, 10))
    m.drawmeridians(range(0, 360, 10))

    cs = m.contourf(x, y, data, cmap=plt.cm.jet, levels=level)  # 定义分级：levels=np.arange(250, 355, 5)
    cb = plt.colorbar(cs, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=28)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)
    axes2.text(0.76, 0.4, "Unit:(" + u'\u2103' + ")", fontsize=32)

    cax.set_title(titlename,
                  family='Times New Roman',
                  fontsize=42,
                  fontweight='bold',pad=20)
    #添加Logo
    axicon=fig.add_axes([0.85,0.01,0.15,0.05])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'),origin='upper')
    axicon.axis('off')	

    fig.savefig(filename)
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
        # nodeOutputFilename1 = doc.createElement("OutputFilename")
        # nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['fy4name']))))
        # nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeOutputFilename2 = doc.createElement("OutputFilename")
        nodeOutputFilename2.appendChild(doc.createTextNode(str(os.path.basename(names['scatname']))))
        nodeOutputFile.appendChild(nodeOutputFilename2)

        nodeOutputFilename3 = doc.createElement("OutputFilename")
        nodeOutputFilename3.appendChild(doc.createTextNode(str(os.path.basename(names['histname']))))
        nodeOutputFile.appendChild(nodeOutputFilename3)

        nodeOutputFilename4 = doc.createElement("OutputFilename")
        nodeOutputFilename4.appendChild(doc.createTextNode(str(os.path.basename(names['biasname']))))
        nodeOutputFile.appendChild(nodeOutputFilename4)

        # nodeOutputFilename5 = doc.createElement("OutputFilename")
        # nodeOutputFilename5.appendChild(doc.createTextNode(str(os.path.basename(names['checkname']))))
        # nodeOutputFile.appendChild(nodeOutputFilename5)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', 'FY4A')
        nodeDeviationstatic.setAttribute('checkVariety', nodes['cVar'])
        nodeDeviationstatic.setAttribute('dataFormat', '15mm')
        nodeDeviationstatic.setAttribute('dataType', nodes['prjType'])
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', 'SST-NA')
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