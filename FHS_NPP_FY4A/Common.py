#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/9 14:51
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
from matplotlib.patches import Ellipse
from mpl_toolkits.basemap import Basemap
import glob
from matplotlib.font_manager import FontProperties


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

        timenode = root.getElementsByTagName('timeValue')[0]
        xmlnodes['timevalue'] = timenode.childNodes[0].nodeValue

        spacenode = root.getElementsByTagName('spaceValue')[0]
        xmlnodes['spacevalue'] = spacenode.childNodes[0].nodeValue

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
    时间匹配：根据FY4时间得到检验数据,并进一步判断在时间范围内的数据是否在FY4圆盘内
    :param path:
    :param base_time:
    :return:
    """
    path = path + base_time.__format__('%Y') + "/" + base_time.__format__('%j')
    files = glob.glob(os.path.join(path, '*.nc'))  # AHI8-_AGRI*N_DISK*FHS*NOM*V0001.NC
    dfiles = []
    dtimes = []
    for filepath in files:
        file = os.path.basename(filepath)
        timeStr = file[7:19]
        time = datetime.datetime.strptime(timeStr, '%Y%j.%H%M')
        if (np.abs((time - base_time).total_seconds() / 60) <= float(timeThreshold)):
            dfiles.append(filepath)
            dtimes.append(time)
    dfiles.sort()
    dtimes.sort()
    if (len(dfiles) == 0):
        raise NameError("未找到检验源数据")
    return dfiles, dtimes


def Name(filename, fy4time, ahi8time, outputpath):
    """
    定义输出的文件名或文件标题
    :param baseInfo:
    :param nodes:
    :param fy3times:
    :return:
    """
    histname = "%s%s_QCS_NPP_HIST_BIAS.png" % (outputpath, filename)
    histtitle = u"FY4A时间:%s UTC NPP时间:%s--%s UTC" % (
        fy4time.__format__('%Y-%m-%d %H:%M:%S'), ahi8time[0].__format__('%Y-%m-%d %H:%M:%S'),
        ahi8time[len(ahi8time) - 1].__format__('%Y-%m-%d %H:%M:%S'))
    biastitle = u"FY4A-NPP火点像元位置匹配示意图\n FY4A时间:%s UTC  NPP时间:%s--%s UTC" % (
        fy4time.__format__('%Y-%m-%d %H:%M:%S'), ahi8time[0].__format__('%Y-%m-%d %H:%M:%S'),
        ahi8time[len(ahi8time) - 1].__format__('%H:%M:%S'))
    BIASMap = "%s%s_QCS_NPP_MAPN_BIAS.png" % (outputpath, filename)
    names = {'histname': histname,
             'histtitle': histtitle,
             'biastitle': biastitle,
             'BIASMap': BIASMap}
    return names


def Read_LUT():
    """
    根据数据的投影类型，分辨率读取经纬度查找表
    :param prjType:
    :param filepath:
    :param dataRes:
    :return:
    """
    lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/FullMask_Grid_2000_1047E_BSQ.nc'

    with Dataset(lutpath, 'r') as fid:
        lat = fid.variables['LAT'][:]
        lon = fid.variables['LON'][:]
    lon[np.where(lon > 180.)] = lon[np.where(lon > 180.)] - 360.

    return lat, lon


def Read_Data(filepath):
    """
    读取数据，根据变量名
    :param filepath:
    :return: FHS数据,10:火点，array类型
    """
    # 读取FY4数据
    with Dataset(filepath, 'r') as fid:
        fy4_data = fid.variables["FHS"][:]
        fpt = fid.variables['FPT'][5:]
    if (isinstance(fy4_data, ma.MaskedArray)):
        fy4_data = fy4_data.data
    if (fy4_data.dtype == 'uint16'):
        fy4_data.dtype = np.int16
        fy4_data = fy4_data.byteswap()
    fy4_mdata = ma.masked_values(fy4_data, 255)
    lst = [str(fpt[i]).split()[2] for i in range(fpt.shape[0])]
    reli = np.array(lst, dtype='int')

    return fy4_mdata, reli


def ReadExterData(filespath):
    """
    读取检验源数据
    :param mod_path:
    :return:
    """
    lats = None
    lons = None
    for filepath in filespath:
        with Dataset(filepath, 'r') as fid:
            lat = fid.variables["FP_latitude"]
            lon = fid.variables["FP_longitude"]
            if lat.shape[0]:
                if lats is None:
                    lats = lat
                    lons = lon
                else:
                    lats = np.hstack((lats, lat))
                    lons = np.hstack((lons, lon))
    return lats, lons


def Interp(tlat, tlon, slat, slon, spaceThreshold):
    tree = spatial.cKDTree(zip(slat, slon))  # 先将待插值数据建树，latlon不能有NA值
    d, inds = tree.query(zip(tlat.flatten(), tlon.flatten()), k=1, distance_upper_bound=float(spaceThreshold))
    nearest = -999 * np.ones(tlat.size)  # 阈值之外的赋值为-999
    nearest[~np.isinf(d)] = 1
    nearest = ma.masked_values(nearest, -999)
    return nearest.reshape(tlat.shape)


def MatchFY4_AHI8(fy4_mdata, mod_mdata, Reli):
    """
    FY4与MOD空间匹配，若FY4火点7*7范围内有MOD的火点，则表示该火点匹配成功，匹配上的点为3，返回二者匹配的行列号和匹配结果
    :param fy4_mdata:
    :param mod_mdata:
    :return:
    """
    FY4 = np.array(np.where(fy4_mdata == 10)).T
    AHI8 = np.array(np.where(mod_mdata == 1)).T
    ahi8_Line = AHI8[:, 0]
    ahi8_Column = AHI8[:, 1]
    m1 = FY4.shape[0]
    m2 = AHI8.shape[0]
    m_a = np.ones(m1).reshape(m1, 1)  # 匹配结果初始化A中火点默认为1
    m_b = np.ones(m2).reshape(m2, 1) * 2  # 匹配结果初始化B中火点默认为2
    # 遍历FY4的行，哪些行的7*7范围内有AHI8的值的索引
    ind = [AHI8[(ahi8_Line >= lineColumn[0] - 3) &
                (ahi8_Line <= lineColumn[0] + 3) &
                (ahi8_Column >= lineColumn[1] - 3) &
                (ahi8_Column <= lineColumn[1] + 3)].shape[0] != 0 for lineColumn in FY4]
    m_a[ind] = 3  # 匹配上的火点为3

    Am = np.hstack((FY4, m_a))
    Bm = np.hstack((AHI8, m_b))
    M = np.vstack((Am, Bm))
    # -----------统计个数------------
    totalNum = np.array([sum(Reli == 1), sum(Reli == 2), sum(Reli == 3), sum(Reli == 4)])
    trueNum = np.array([sum((m_a.flatten() == 3) & (Reli == 1)),  # 轻度火点像元
                        sum((m_a.flatten() == 3) & (Reli == 2)),  # 中度火点像元
                        sum((m_a.flatten() == 3) & (Reli == 3)),  # 重度火点像元
                        sum((m_a.flatten() == 3) & (Reli == 4))])  # 云区火点像元
    statisRes = np.array([totalNum, trueNum])
    return M, statisRes


def DrawHist(statisRes, histname, histtitle):
    trueNum = statisRes[1]
    totalNum = statisRes[0]
    acc = np.flipud(trueNum * 1.0 / totalNum)

    # plt.rcParams['font.family'] = ['SimHei']  # 用来正常显示中文标签
    font = FontProperties(fname=r"/FY4APGSQCS/QCS/src/FHS/FHS_NPP_FY4A/Logo/simhei.ttf")

    matplotlib.rcParams['font.sans-serif'] = ['SimHei']
    matplotlib.rcParams['font.family'] = 'sans-serif'
    fig = plt.figure(figsize=(10.24, 5), dpi=100)  # 图像的长*高=2748*3000
    axes1 = fig.add_axes([0., 0., 1., 0.7], facecolor='#F5F2E1')
    axes2 = fig.add_axes([0.75, 0.1, 0.2, 0.5], facecolor='#F5F2E1')
    fig.text(0.05, 0.9, u"FY4A-NPP火点匹配准确率", fontsize=18, fontproperties=font)
    fig.text(0.05, 0.8, histtitle, fontsize=16, fontproperties=font)

    axes1.text(0.05, 0.95, u'      匹配类型                 总个数                 准确个数                 准确率', fontsize=14,
               fontweight='bold', fontproperties=font)

    # axes1.text(0.125, 0.8, u'火区', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.72, u'轻度火点像元', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.56, u'中度火点像元', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.4, u'重度火点像元', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.24, u'云区火点像元', fontsize=14, fontproperties=font)

    # axes1.text(0.36, 0.8, totalNum[0], fontsize=14)
    axes1.text(0.36, 0.72, totalNum[0], fontsize=14)
    axes1.text(0.36, 0.56, totalNum[1], fontsize=14)
    axes1.text(0.36, 0.4, totalNum[2], fontsize=14)
    axes1.text(0.36, 0.24, totalNum[3], fontsize=14)

    # axes1.text(0.6, 0.8, trueNum[0], fontsize=14)
    axes1.text(0.6, 0.72, trueNum[0], fontsize=14)
    axes1.text(0.6, 0.56, trueNum[1], fontsize=14)
    axes1.text(0.6, 0.4, trueNum[2], fontsize=14)
    axes1.text(0.6, 0.24, trueNum[3], fontsize=14)

    axes1.axhline(y=.92, xmin=0.05, xmax=0.95, color='k')
    axes1.axhline(y=.1, xmin=0.05, xmax=0.95, color='k')
    axes1.text(0.05, 0.05, u'数据来源:国家卫星气象中心', fontsize=10, fontproperties=font)
    axes1.spines['left'].set_visible(False)
    axes1.spines['top'].set_visible(False)

    axes2.barh([0.5, 1.5, 2.5, 3.5], acc, color=["#fcb300", "#007989", "#ffaa91", "#7a0251"],
               height=0.4)
    axes2.set_ylim(0, 4)
    axes2.set_xlim(0, 1)
    # axes2.text(acc[0] - 0.2, 0.5, round(acc[0], 2), fontsize=14, color='white')
    axes2.text(acc[0] - 0.2, 0.4, round(acc[0], 2), fontsize=14)
    axes2.text(acc[1] - 0.2, 1.4, round(acc[1], 2), fontsize=14)
    axes2.text(acc[2] - 0.2, 2.4, round(acc[2], 2), fontsize=14)
    axes2.text(acc[3] - 0.2, 3.4, round(acc[3], 2), fontsize=14)
    axes2.axis('off')
    fig.savefig(histname, facecolor="#F5F2E1")
    plt.close()
    return acc


def DrawBIASMap(data, lat, lon, titlename, filename):
    font = FontProperties(fname=r"/FY4APGSQCS/QCS/src/FHS/FHS_NPP_FY4A/Logo/simhei.ttf")

    matplotlib.rcParams['font.sans-serif'] = ['SimHei']
    matplotlib.rcParams['font.family'] = 'sans-serif'

    bias = np.zeros(lat.shape)
    data = data.astype(int)
    bias[data[:, 0], data[:, 1]] = data[:, 2]
    mbias = ma.masked_values(bias, 0)

    fig = plt.figure(figsize=(27.48, 30), dpi=100)  # 图像的长*高=2748*3000
    axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
    axes1 = fig.add_axes([0., 0.084, 1., 0.916])
    mlat = ma.masked_values(lat, -999.)
    mlon = ma.masked_values(lon, -999.)
    m = Basemap(projection='nsper', lat_0=0, lon_0=104.7, resolution='l', ax=axes1)
    x, y = m(mlon, mlat)
    m.drawcoastlines()
    m.drawcountries()
    m.drawlsmask(land_color='#bebebe', ocean_color='#01008a')  # 陆地，海洋的颜色
    m.drawparallels(np.arange(-90., 90., 10.))
    m.drawmeridians(np.arange(0., 360., 10.))

    m.scatter(x[mbias == 1], y[mbias == 1], s=50, color="red")
    m.scatter(x[mbias == 2], y[mbias == 2], s=50, color="blue")
    m.scatter(x[mbias == 3], y[mbias == 3], s=50, color="yellow")

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)
    axes1.spines['bottom'].set_linewidth(2.)
    axes1.spines['left'].set_linewidth(2.)
    axes1.spines['top'].set_linewidth(2.)
    axes1.spines['right'].set_linewidth(2.)

    axes1.set_xlabel(titlename,
                     family='Times New Roman',
                     fontsize=40,
                     fontweight='bold', labelpad=20, fontproperties=font)

    colors = ["red", "blue", "yellow"]
    names = [u"FY4A未匹配像元", u"NPP未匹配像元", u"匹配像元"]
    for i in range(len(colors)):
        xx = 0.25 + 0.15 * i
        yy = 0.25
        axes2.add_patch(Ellipse(xy=(xx, yy), width=0.015, height=0.2, color=colors[i]))
        axes2.text(xx + 0.01, yy - 0.04, names[i], fontsize=30, fontproperties=font)
    #添加Logo
    axicon=fig.add_axes([0.86,0.01,0.15,0.05])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'),origin='upper')
    axicon.axis('off')		

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
        nodeOutputFilename0.appendChild(doc.createTextNode(str(os.path.basename(names['histname']))))
        nodeOutputFile.appendChild(nodeOutputFilename0)

        nodeOutputFilename1 = doc.createElement("OutputFilename")
        nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['BIASMap']))))
        nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', 'FY4A')
        nodeDeviationstatic.setAttribute('checkVariety', "NPP")
        nodeDeviationstatic.setAttribute('dataFormat', '15mm')
        nodeDeviationstatic.setAttribute('dataType', 'NOM')
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', 'FHS')
        metricStr = "{MildPixel:%s, ModeratePixel:%s, SeverePixel:%s,CloudPixel:%s}" % (
            metric[0], metric[1], metric[2], metric[3])
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