#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/9/27 10:01
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
from scipy import spatial
import glob
from shapely.geometry import Point, Polygon
import codecs
from matplotlib.font_manager import FontProperties
from matplotlib.collections import PatchCollection
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
    xmlnodes = collections.OrderedDict()
    with open(xmlpath, 'r') as fh:
        dom = minidom.parse(fh)
        root = dom.documentElement
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

        checkVarnode = root.getElementsByTagName('checkVariety')[0]
        xmlnodes['checkVar'] = checkVarnode.childNodes[0].nodeValue

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


def TimeMatch(path, base_time, timeThreshold, checkVar):
    """
    时间匹配：根据FY4时间得到检验数据
    :param path:
    :param base_time:
    :return:
    """
    path = path + "/*"
    files = glob.glob(os.path.join(path, '*PLSF_*_%s*.txt' % checkVar))
    dfiles = []
    dtimes = []

    for filepath in files:
        file = os.path.basename(filepath)	
        if (checkVar == 'FY3B')and(len(file)==48):
            timeStr = file[30:44]
        elif(checkVar == 'FY3D')and(len(file)==49):
            timeStr = file[31:45]

        time = datetime.datetime.strptime(timeStr, '%Y%m%d%H%M%S')
        if (np.floor(np.abs((time - base_time).total_seconds() / 60)) <= float(timeThreshold)):
            dfiles.append(filepath)
            dtimes.append(time)
    dfiles.sort()
    dtimes.sort()
    if (len(dfiles) == 0):
        raise NameError("未找到检验源数据")
    return dfiles, dtimes


def Name(baseInfo, nodes, fy3times, checkVar):
    """
    定义输出的文件名或文件标题
    :param baseInfo:
    :param nodes:
    :param fy3times:
    :return:
    """
    outputpath = nodes['outputpath']
    filename = baseInfo['filename']
    fy4time = baseInfo['time']

    txtpath = "%s%s_QCS_%s_DATA_MATCH.TXT" % (outputpath, filename, checkVar)
    ncpath = "%s%s_QCS_%s_DATA_BIAS.NC" % (outputpath, filename, checkVar)
    histname = "%s%s_QCS_%s_HIST_BIAS.png" % (outputpath, filename, checkVar)
    histtitle = u"FY4A时间:%s UTC %s时间:%s UTC" % (
        fy4time.__format__('%Y-%m-%d %H:%M:%S'), checkVar, fy3times[0].__format__('%Y-%m-%d %H:%M:%S'))
    pixelCondMap = "%s%s_QCS_%s_MAPN_PIXEL_CONDITION.png" % (outputpath, filename, checkVar)
    areaSketchMap = "%s%s_QCS_%s_MAPN_AREA_SKETCH.png" % (outputpath, filename, checkVar)
    pixelSketchMap = "%s%s_QCS_%s_MAPN_PIXEL_SKETCH.png" % (outputpath, filename, checkVar)

    names = {'txtpath': txtpath,
             'histname': histname,
             'ncpath': ncpath,
             'histtitle': histtitle,
             'pixelCondMap': pixelCondMap,
             'areaSketchMap': areaSketchMap,
             'pixelSketchMap': pixelSketchMap}
    return names


def ConvertData(fy4_path, txtpath):
    """
    读取NC文件，并将其输出为TXT
    :param fy4_path: nc文件路径
    :param txtpath: txt文件路径
    :return:
    """
    # 读取FY4数据
    with Dataset(fy4_path, 'r') as fid_fy4:
        fhs = fid_fy4.variables['FHS'][:]
        fpt = fid_fy4.variables['FPT'][4:]
        BNL = fid_fy4.variables['FPT'][2:3]
    beginLine = int(str(BNL).split()[1])
    endLine = int(str(BNL).split()[3].replace(r"\r\n']", ""))

    if (fhs.dtype == 'uint16'):
        fhs.dtype = np.int16
        fhs = fhs.byteswap()
    lst = []
    for i in range(fpt.shape[0]):
        tmp = str(fpt[i]).split()
        lst.append(tmp)
    arr = np.array(lst)
    row, col = np.where(fhs == 10)
    # row = row + 1
    # col = col + 1
    col = col.astype(str)
    row = row.astype(str)
    rowstr = np.insert(col, 0, 'Line.')
    colstr = np.insert(row, 0, 'Column.')
    arr = np.insert(arr, 1, rowstr, axis=1)
    arr = np.insert(arr, 2, colstr, axis=1)
    # 写出到txt文件
    if (os.path.exists(txtpath)):
        os.remove(txtpath)
    with open(txtpath, "a") as f:
        for i in range(arr.shape[0]):
            write_str = '%s %s %s %s %s %s %s %s %s %s %s\n' % (
                arr[i, 0], arr[i, 1], arr[i, 2], arr[i, 3], arr[i, 4], arr[i, 5], arr[i, 6], arr[i, 7], arr[i, 8],
                arr[i, 9], arr[i, 10])
            f.write(write_str)
    return txtpath, beginLine, endLine


def Read_Data(txtpath):
    """
    读取FY4txt数据
    :param txtpath:
    :return:
    """
    arr = np.loadtxt(txtpath, skiprows=1)
    return arr


def Read_exterData(files):
    """
    读取FY3数据，若匹配到多个文件，则合并为一个
    :param files:
    :return:
    """
    arr = np.empty(shape=(0, 11))
    for filepath in files:
        filecp = codecs.open(filepath, encoding='cp1252')
        tmpArr = np.loadtxt(filecp, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12))
        arr = np.vstack((arr, tmpArr))
    return arr


def Read_LUT(prjType, filepath, dataRes, beginLine, endLine):
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
            lat = fid.variables['LAT'][beginLine:endLine + 1, :]
            lon = fid.variables['LON'][beginLine:endLine + 1, :]
        lon[np.where(lon > 180.)] = lon[np.where(lon > 180.)] - 360.
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


def PixelNearFY3(fy4dat, fy3dat, spaceThreshold):
    """
    输入FY4,FY3数据，得到在FY3范围内的FY4火点，将FY3插值到FY4上面后，进一步判断在FY3范围内的FY4火点是否在FY3一定空间阈值范围内
    :param fy4dat:
    :param fy3dat:
    :return:fy4dat_:完全在FY3范围内的点;nearOrNot:FY4火点邻近区域是否有FY3点，有1;无0
    """
    m = fy4dat.shape[0]
    tlat = fy4dat[:, 5]
    tlon = fy4dat[:, 6]

    slat = fy3dat[:, 1]
    slon = fy3dat[:, 2]

    fy3minlat = slat.min() - 1
    fy3maxlat = slat.max() + 1
    fy3minlon = slon.min() - 1
    fy3maxlon = slon.max() + 1

    # 2.转换经纬度到笛卡尔坐标系
    xs, ys, zs = lon_lat_to_cartesian(slon, slat)
    xt, yt, zt = lon_lat_to_cartesian(tlon, tlat)

    tree = spatial.cKDTree(zip(xs, ys, zs))
    d, inds = tree.query(zip(xt, yt, zt), k=1, distance_upper_bound=float(spaceThreshold))
    nearOrNot = np.zeros(m)
    nearOrNot[d != np.inf] = 1
    fy4dat_ = fy4dat[(fy3maxlat >= tlat) & (tlat >= fy3minlat) & (fy3maxlon >= tlon) & (tlon >= fy3minlon)]
    return fy4dat_, nearOrNot


def PixelRange(fy4dat, LAT, LON):
    """
    输入火点行列号，得到其经纬度范围，左上角，右下角，左下角，右上角
    :param fy4dat: FY4数据（包括行列号信息）
    :param LAT: 纬度
    :param LON: 经度
    :return:
    """
    lines = fy4dat[:, 2]  # 行号
    lines = lines.astype("int")
    columns = fy4dat[:, 1]  # 列号
    columns = columns.astype("int")

    lons = np.empty(shape=(0, 5))
    lats = np.empty(shape=(0, 5))
    for line, column in zip(lines, columns):
        # 火点经纬度
        clat = LAT[line, column]
        clon = LON[line, column]

        # 火点左上角经纬度
        ltlat = LAT[line - 1, column - 1]
        ltlon = LON[line - 1, column - 1]
        # 火点左下角经纬度
        lblat = LAT[line + 1, column - 1]
        lblon = LON[line + 1, column - 1]
        # 火点右上角经纬度
        rtlat = LAT[line - 1, column + 1]
        rtlon = LON[line - 1, column + 1]
        # 火点右下角经纬度
        rblat = LAT[line + 1, column + 1]
        rblon = LON[line + 1, column + 1]
        tmpLon = np.array(
            [(ltlon + clon) / 2., (lblon + clon) / 2., (rblon + clon) / 2., (rtlon + clon) / 2., clon]).reshape(1, 5)
        tmpLat = np.array(
            [(ltlat + clat) / 2., (lblat + clat) / 2., (rblat + clat) / 2., (rtlat + clat) / 2., clat]).reshape(1, 5)
        lons = np.vstack((lons, tmpLon))
        lats = np.vstack((lats, tmpLat))
    ranges = {"LON": lons, "LAT": lats}
    return ranges


def InterpFY4(ranges, spaceThreshold):
    """
    插值:输入火点像元的经纬度范围，进行插值，有可能出现某个火点的一个角无效插值情况，因此对无效插值的火点进行剔除，
    返回插值后等经纬格网最邻近的行列号
    :param ranges:
    :param spaceThreshold:
    :return:
    """
    m, n = ranges["LON"].shape
    tlon2d, tlat2d = np.meshgrid(np.arange(70, 140, 0.01), np.arange(55, 15, -0.01))
    x, y = tlon2d.shape
    slons = ranges["LON"].flatten(order='F')
    slats = ranges["LAT"].flatten(order='F')

    # 2.转换经纬度到笛卡尔坐标系
    xs, ys, zs = lon_lat_to_cartesian(slons, slats)
    xt, yt, zt = lon_lat_to_cartesian(tlon2d.flatten(), tlat2d.flatten())

    tree = spatial.cKDTree(zip(xt, yt, zt))
    d, inds = tree.query(zip(xs, ys, zs), k=1, distance_upper_bound=float(spaceThreshold))
    inds_arr = np.array(inds).reshape(m, n, order='F')
    d_arr = np.array(d).reshape(m, n, order='F')
    rowValid = np.where(np.sum(d_arr, axis=1) != np.inf)[0]  # 四个角点均插值为有效值的行
    inds_arr = inds_arr[rowValid, :]  # 剔除含有inf值的行
    # 将索引转换为行列号
    row_arr = (inds_arr / y)
    col_arr = (inds_arr % y)
    return row_arr, col_arr, rowValid


def PixelGLL(rows, cols):
    """
    输入四个角点的行列号，得到在其四边形内部的点的行列号
    :param rows:array类型，shape=(m,5) 5列分别对应左上，左下，右下，右上，中心点的行号
    :param cols:array类型，shape=(m,5) 5列分别对应左上，左下，右下，右上，中心点的列号
    :return:
    """
    rowCol = np.empty(shape=(0, 2))
    for row, col in zip(rows[:, 0:4], cols[:, 0:4]):
        rowMin = row.min()
        rowMax = row.max()

        colMin = col.min()
        colMax = col.max()
        axisX = np.repeat(np.arange(rowMin, rowMax + 1), colMax - colMin + 1)
        axisY = np.tile(np.arange(colMin, colMax + 1), rowMax - rowMin + 1)
        polygon = Polygon([(row[0], col[0]), (row[1], col[1]), (row[2], col[2]), (row[3], col[3])])

        for ax, ay in zip(axisX, axisY):
            point = Point(ax, ay)
            if polygon.intersects(point):
                rowCol = np.vstack((rowCol, point))
    return rowCol


def InterpFY3(fy3dat, spaceThreshold):
    slat = fy3dat[:, 1]
    slon = fy3dat[:, 2]
    tlon2d, tlat2d = np.meshgrid(np.arange(70, 140, 0.01), np.arange(55, 15, -0.01))
    x, y = tlon2d.shape
    # 2.转换经纬度到笛卡尔坐标系
    xs, ys, zs = lon_lat_to_cartesian(slon, slat)
    xt, yt, zt = lon_lat_to_cartesian(tlon2d.flatten(), tlat2d.flatten())

    tree = spatial.cKDTree(zip(xt, yt, zt))
    d, inds = tree.query(zip(xs, ys, zs), k=1, distance_upper_bound=float(spaceThreshold))
    # 将索引转换为行列号
    row = inds[d != np.inf] / y
    col = inds[d != np.inf] % y
    rowCol = np.vstack((row, col)).T
    return rowCol


def intersect2d(X, Y):
    """
    Function to find intersection of two 2D arrays.
    Returns index of rows in X that are common to Y.
    """
    X = np.tile(X[:, :, None], (1, 1, Y.shape[0]))
    Y = np.swapaxes(Y[:, :, None], 0, 2)
    Y = np.tile(Y, (X.shape[0], 1, 1))
    eq = np.all(np.equal(X, Y), axis=1)
    eq = np.any(eq, axis=1)
    return np.nonzero(eq)[0]


def MatchFY(A, B):
    """
    :param A: FY4火点的行列号,array(m,2)
    :param B: FY3火点的行列号,array(n,2)
    :return: 重复的个数为x,则M,array(m+n-x,3),第1列为行号，第2列为列号，第3列为匹配结果，FY4为1，FY3为2，匹配上的为3
    """
    m1 = A.shape[0]
    m2 = B.shape[0]

    ind1 = intersect2d(A, B)  # A中的行在B中的索引
    ind2 = intersect2d(B, A)  # B中的行在A中的索引

    m_a = np.ones(m1).reshape(m1, 1)  # 匹配结果初始化A中火点默认为1
    m_b = np.ones(m2).reshape(m2, 1) * 2  # 匹配结果初始化B中火点默认为2

    m_a[ind1] = 3  # 重复的火点为3
    Am = np.hstack((A, m_a))
    Bm = np.delete(np.hstack((B, m_b)), ind2, 0)

    M = np.vstack((Am, Bm))
    return M


def ReJectFY3(MatchPixel):
    """
    剔除FY3 7*7范围内无FY4点的火点
    :param MatchPixel:
    :return:
    """
    FY3 = MatchPixel[MatchPixel[:, 2] == 2, :]
    FY4 = MatchPixel[MatchPixel[:, 2] == 1, :]

    FY4_Line = FY4[:, 0]
    FY4_Column = FY4[:, 1]

    # 遍历FY3的行，哪些行的7*7范围内有FY4的值的索引
    ind = [FY4[(FY4_Line >= lineColumn[0] - 3) &
               (FY4_Line <= lineColumn[0] + 3) &
               (FY4_Column >= lineColumn[1] - 3) &
               (FY4_Column <= lineColumn[1] + 3)].shape[0] != 0 for lineColumn in FY3]

    MatchPixel_ = np.delete(MatchPixel, np.where(MatchPixel[:, 2] == 2)[0], axis=0)
    MatchPixel_reject = np.vstack((MatchPixel_, FY3[ind]))
    return MatchPixel_reject


def Match_Dot_Area(fy4dat, row_arr, col_arr, nearOrNot, rowValid):
    centerPixel_row = 4000-row_arr[:, 4]
    centerPixel_col = col_arr[:, 4]
    MatchDot = np.vstack((centerPixel_row, centerPixel_col, nearOrNot[rowValid])).T
    Reli = fy4dat[rowValid, 4]  # 火点像元可信度

    spot_no = fy4dat[rowValid, 3]  # 火区编号
    pixel_Size = fy4dat[rowValid, 7]  # 火点像元面积
    sno_unique = np.unique(spot_no)  # 火区编号唯一值

    nearOrNot_ = np.zeros(nearOrNot[rowValid].shape[0])  # 一个火区的像元匹配上，则认为该火区匹配上
    for sno in sno_unique:
        nearOrNot_[np.where(spot_no == sno)[0]] = nearOrNot[rowValid][np.where(spot_no == sno)[0]].max()

    inds = [np.where(spot_no == sno)[0][0] for sno in sno_unique]  # Pixel_Hot_Spot_NO第一次出现的索引
    spot_area = [pixel_Size[np.where(spot_no == sno)[0]].sum() for sno in sno_unique]
    MatchArea = np.vstack((centerPixel_row[inds], centerPixel_col[inds], nearOrNot_[inds], spot_area)).T
    TrueNum = np.array([sum(nearOrNot_[inds] == 1),  # 匹配上的火区的个数
                        sum((nearOrNot[rowValid] == 1) & (Reli == 1)),  # 轻度火点像元
                        sum((nearOrNot[rowValid] == 1) & (Reli == 2)),  # 中度火点像元
                        sum((nearOrNot[rowValid] == 1) & (Reli == 3)),  # 重度火点像元
                        sum((nearOrNot[rowValid] == 1) & (Reli == 4))])  # 云区火点像元
    return MatchDot, MatchArea, TrueNum


def DrawHist(fy4dat_check, trueNum, histname, histtitle, checkVar):
    Reli = fy4dat_check[:, 4]  # 火点像元可信度
    # FY3观测区内FY4火区个数和轻度火点像元，中度火点像元,重度火点像元,云区火点像元个数统计
    totalNum = np.array(
        [np.unique(fy4dat_check[:, 3]).size, sum(Reli == 1), sum(Reli == 2), sum(Reli == 3), sum(Reli == 4)])
    acc = np.flipud(trueNum * 1.0 / totalNum)

    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")

    matplotlib.rcParams['font.sans-serif'] = ['SimHei']
    matplotlib.rcParams['font.family'] = 'sans-serif'
    fig = plt.figure(figsize=(10.24, 5), dpi=100)  # 图像的长*高=2748*3000
    axes1 = fig.add_axes([0., 0., 1., 0.7], facecolor='#F5F2E1')
    axes2 = fig.add_axes([0.75, 0.1, 0.2, 0.5], facecolor='#F5F2E1')
    fig.text(0.05, 0.9, u"FY4A-%s火点匹配准确率" % checkVar, fontsize=18, fontproperties=font)
    fig.text(0.05, 0.8, histtitle, fontsize=16, fontproperties=font)

    axes1.text(0.05, 0.95, u'      匹配类型                 总个数                 准确个数                 准确率', fontsize=14,
               fontweight='bold', fontproperties=font)

    axes1.text(0.125, 0.8, u'火区', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.65, u'轻度火点像元', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.5, u'中度火点像元', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.35, u'重度火点像元', fontsize=14, fontproperties=font)
    axes1.text(0.09, 0.2, u'云区火点像元', fontsize=14, fontproperties=font)

    axes1.text(0.36, 0.8, totalNum[0], fontsize=14)
    axes1.text(0.36, 0.65, totalNum[1], fontsize=14)
    axes1.text(0.36, 0.5, totalNum[2], fontsize=14)
    axes1.text(0.36, 0.35, totalNum[3], fontsize=14)
    axes1.text(0.36, 0.2, totalNum[4], fontsize=14)

    axes1.text(0.6, 0.8, trueNum[0], fontsize=14)
    axes1.text(0.6, 0.65, trueNum[1], fontsize=14)
    axes1.text(0.6, 0.5, trueNum[2], fontsize=14)
    axes1.text(0.6, 0.35, trueNum[3], fontsize=14)
    axes1.text(0.6, 0.2, trueNum[4], fontsize=14)

    axes1.axhline(y=.92, xmin=0.05, xmax=0.95, color='k')
    axes1.axhline(y=.1, xmin=0.05, xmax=0.95, color='k')
    axes1.text(0.05, 0.05, u'数据来源:国家卫星气象中心', fontsize=10, fontproperties=font)
    axes1.spines['left'].set_visible(False)
    axes1.spines['top'].set_visible(False)

    axes2.barh([0.6, 1.6, 2.7, 3.7, 4.7], acc, color=["#ff5961", "#fcb300", "#007989", "#ffaa91", "#7a0251"],
               height=0.4)
    axes2.set_ylim(0, 5)
    axes2.set_xlim(0, 1)

    axes2.text(acc[0] - 0.2, 0.5, round(acc[0], 2), fontsize=14)
    axes2.text(acc[1] - 0.2, 1.5, round(acc[1], 2), fontsize=14)
    axes2.text(acc[2] - 0.2, 2.6, round(acc[2], 2), fontsize=14)
    axes2.text(acc[3] - 0.2, 3.6, round(acc[3], 2), fontsize=14)
    axes2.text(acc[4] - 0.2, 4.6, round(acc[4], 2), fontsize=14)
    axes2.axis('off')
    fig.savefig(histname, facecolor="#F5F2E1")
    plt.close()
    return acc


def PixelConditionMap(fire, fy4time, fy3time, output, checkVar):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    data = np.ones(shape=(4000, 7000)) * -999
    if fire.shape[0] > 0:
        fire = fire.astype("int")
        data[fire[:, 0], fire[:, 1]] = fire[:, 2]
        data = data[::-1]
    mdata = ma.masked_values(data, -999)
    mdata = np.flipud(mdata)	
    lon = np.arange(70, 140, 0.01)
    lat = np.arange(55, 15, -0.01)
    lon2d, lat2d = np.meshgrid(lon, lat)

    minLon = np.min(lon2d[(mdata == 2) | (mdata == 3)])
    maxLon = np.max(lon2d[(mdata == 2) | (mdata == 3)])
    minLat = np.min(lat2d[(mdata == 2) | (mdata == 3)])
    maxLat = np.max(lat2d[(mdata == 2) | (mdata == 3)])

    # 定义底图
    fig = plt.figure(figsize=(100, 100), dpi=100)
    axes1 = fig.add_axes([0., 0., 0.7, 0.4])
    axes2 = fig.add_axes([0.6, 0.005, 0.1, 0.1])

    m = Basemap(projection='cyl', resolution='l', ax=axes1,
                llcrnrlat=15, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)
    m2 = Basemap(projection='cyl', resolution='l', ax=axes2,
                 llcrnrlat=2, llcrnrlon=105, urcrnrlon=123, urcrnrlat=23)
    m.drawparallels(range(-90, 90, 5), labels=[1, 0, 0, 0], fontsize=38, color='#33b0ef')
    m.drawmeridians(range(0, 360, 5), labels=[0, 0, 0, 1], fontsize=38, color='#33b0ef')
    m2.drawparallels(range(-90, 90, 5), color='#33b0ef')
    m2.drawmeridians(range(0, 360, 5), color='#33b0ef')
    m.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/country', 'country')
    m.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/provice', 'provice')
    m2.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/country', 'country')
    m2.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/provice', 'provice')
    patches = []
    patches2 = []
    for info, shape, info2, shape2 in zip(m.country_info, m.country, m2.country_info, m2.country):
        if info['SOC'] == 'CHN':
            patches.append(matplotlib.patches.Polygon(np.array(shape), True))
        if info2['SOC'] == 'CHN':
            patches2.append(matplotlib.patches.Polygon(np.array(shape2), True))
    axes1.add_collection(PatchCollection(patches, facecolor="#bebebe"))
    axes2.add_collection(PatchCollection(patches2, facecolor="#bebebe"))
    axes1.set_title(u"FY4A-%s火点像元位置匹配实况图" % checkVar, fontproperties=font, fontsize=120)

    # 添加本地数据
    x, y = m(lon2d, lat2d)
    x2, y2 = m2(lon2d, lat2d)
    cm_light = matplotlib.colors.ListedColormap(['red', 'blue', 'yellow'])
    m.pcolormesh(x, y, mdata, cmap=cm_light)
    m2.pcolormesh(x2, y2, mdata, cmap=cm_light)
    # 绘制FY3B范围的矩形框
    xmin, ymin = m(minLon - 0.1, minLat - 0.1)
    xmax, ymax = m(maxLon + 0.1, maxLat + 0.1)
    po = matplotlib.patches.Polygon([[xmax, ymin], [xmin, ymin], [xmin, ymax], [xmax, ymax]], fill=False,
                                    edgecolor="#00ffff", linewidth=4)
    axes1.add_patch(po)
    # 加气象局logo
    axicon = fig.add_axes([0.605, 0.355, 0.09, 0.04])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/%s/logo.png' % checkVar), origin='upper')
    axicon.axis('off')
    # 添加时间
    axtime = fig.add_axes([0.007, 0.368, 0.135, 0.023])
    axtime.text(0.04, 0.1, "FY4A Time:%s\n%s Time:%s" % (
        fy4time.__format__('%Y-%m-%d %H:%M'), checkVar, fy3time.__format__('%Y-%m-%d %H:%M')), fontsize=64)
    axtime.set_xticks([])
    axtime.set_yticks([])
    # 添加图例
    axlegend = fig.add_axes([0.005, 0.005, 0.07, 0.05])
    axlegend.imshow(plt.imread('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/%s/Lengend-PixelCondition.png' % checkVar), origin='upper')
    axlegend.axis('off')
    # 保存图片
    fig.savefig(output, bbox_inches='tight')
    plt.close()
    return 0, [xmin, xmax, ymin, ymax]


def AreaSketchMap(firezone, fy4time, fy3time, output, rect, checkVar):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    xmin, xmax, ymin, ymax = rect

    lon = np.arange(70, 140, 0.01)
    lat = np.arange(55, 15, -0.01)
    lon2d, lat2d = np.meshgrid(lon, lat)

    lat2d = lat2d[::-1]
    lon2d = lon2d[::-1]

    slat = lat2d[firezone[:, 0].astype("int"), firezone[:, 1].astype("int")]
    slon = lon2d[firezone[:, 0].astype("int"), firezone[:, 1].astype("int")]

    # 定义底图
    fig = plt.figure(figsize=(100, 100), dpi=100)
    axes1 = fig.add_axes([0., 0., 0.7, 0.4])
    axes2 = fig.add_axes([0.6, 0.005, 0.1, 0.1])

    m = Basemap(projection='cyl', resolution='l', ax=axes1,
                llcrnrlat=15, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)
    m2 = Basemap(projection='cyl', resolution='l', ax=axes2,
                 llcrnrlat=2, llcrnrlon=105, urcrnrlon=123, urcrnrlat=23)
    m.drawparallels(range(-90, 90, 5), labels=[1, 0, 0, 0], fontsize=38, color='#33b0ef')
    m.drawmeridians(range(0, 360, 5), labels=[0, 0, 0, 1], fontsize=38, color='#33b0ef')
    m2.drawparallels(range(-90, 90, 5), color='#33b0ef')
    m2.drawmeridians(range(0, 360, 5), color='#33b0ef')
    m.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/country', 'country')
    m.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/provice', 'provice')
    m2.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/country', 'country')
    m2.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/provice', 'provice')
    patches = []
    patches2 = []
    for info, shape, info2, shape2 in zip(m.country_info, m.country, m2.country_info, m2.country):
        if info['SOC'] == 'CHN':
            patches.append(matplotlib.patches.Polygon(np.array(shape), True))
        if info2['SOC'] == 'CHN':
            patches2.append(matplotlib.patches.Polygon(np.array(shape2), True))
    axes1.add_collection(PatchCollection(patches, facecolor="#bebebe"))
    axes2.add_collection(PatchCollection(patches2, facecolor="#bebebe"))
    axes1.set_title(u"FY4A-%s火区位置匹配示意图" % checkVar, fontproperties=font, fontsize=120)

    # 添加本地数据
    x, y = m(slon, slat)
    x2, y2 = m2(slon, slat)
    m.scatter(x[firezone[:, 2] == 0], y[firezone[:, 2] == 0], s=firezone[:, 3][firezone[:, 2] == 0], c='red')
    m.scatter(x[firezone[:, 2] == 1], y[firezone[:, 2] == 1], s=firezone[:, 3][firezone[:, 2] == 1], c='yellow')
    m2.scatter(x2[firezone[:, 2] == 0], y2[firezone[:, 2] == 0], s=firezone[:, 3][firezone[:, 2] == 0], c='red')
    m2.scatter(x2[firezone[:, 2] == 1], y2[firezone[:, 2] == 1], s=firezone[:, 3][firezone[:, 2] == 1], c='yellow')
    # 绘制FY3B范围的矩形框
    po = matplotlib.patches.Polygon([[xmax, ymin], [xmin, ymin], [xmin, ymax], [xmax, ymax]], fill=False,
                                    edgecolor="#00ffff", linewidth=4)
    axes1.add_patch(po)
    # 加气象局logo
    axicon = fig.add_axes([0.605, 0.355, 0.09, 0.04])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/%s/logo.png' % checkVar), origin='upper')
    axicon.axis('off')
    # 添加时间
    axtime = fig.add_axes([0.007, 0.368, 0.135, 0.023])
    axtime.text(0.04, 0.1, "FY4A Time:%s\n%s Time:%s" % (
        fy4time.__format__('%Y-%m-%d %H:%M'), checkVar, fy3time.__format__('%Y-%m-%d %H:%M')), fontsize=64)
    axtime.set_xticks([])
    axtime.set_yticks([])
    # 添加图例
    axlegend = fig.add_axes([0.005, 0.005, 0.07, 0.05])
    axlegend.imshow(plt.imread('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/%s/Lengend-AreaSketch.png' % checkVar), origin='upper')
    axlegend.axis('off')
    # 保存图片
    fig.savefig(output, bbox_inches='tight')
    plt.close()
    return 0


def PixelSketchMap(firepixel, fy4time, fy3time, output, rect, checkVar):
    font = FontProperties(fname="/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/simhei.ttf")
    xmin, xmax, ymin, ymax = rect

    lon = np.arange(70, 140, 0.01)
    lat = np.arange(55, 15, -0.01)
    lon2d, lat2d = np.meshgrid(lon, lat)

    lat2d = lat2d[::-1]
    lon2d = lon2d[::-1]

    slat = lat2d[firepixel[:, 0].astype("int"), firepixel[:, 1].astype("int")]
    slon = lon2d[firepixel[:, 0].astype("int"), firepixel[:, 1].astype("int")]

    # 定义底图
    fig = plt.figure(figsize=(100, 100), dpi=100)
    axes1 = fig.add_axes([0., 0., 0.7, 0.4])
    axes2 = fig.add_axes([0.6, 0.005, 0.1, 0.1])

    m = Basemap(projection='cyl', resolution='l', ax=axes1,
                llcrnrlat=15, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)
    m2 = Basemap(projection='cyl', resolution='l', ax=axes2,
                 llcrnrlat=2, llcrnrlon=105, urcrnrlon=123, urcrnrlat=23)
    m.drawparallels(range(-90, 90, 5), labels=[1, 0, 0, 0], fontsize=38, color='#33b0ef')
    m.drawmeridians(range(0, 360, 5), labels=[0, 0, 0, 1], fontsize=38, color='#33b0ef')
    m2.drawparallels(range(-90, 90, 5), color='#33b0ef')
    m2.drawmeridians(range(0, 360, 5), color='#33b0ef')
    m.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/country', 'country')
    m.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/provice', 'provice')
    m2.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/country', 'country')
    m2.readshapefile('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/provice', 'provice')
    patches = []
    patches2 = []
    for info, shape, info2, shape2 in zip(m.country_info, m.country, m2.country_info, m2.country):
        if info['SOC'] == 'CHN':
            patches.append(matplotlib.patches.Polygon(np.array(shape), True))
        if info2['SOC'] == 'CHN':
            patches2.append(matplotlib.patches.Polygon(np.array(shape2), True))
    axes1.add_collection(PatchCollection(patches, facecolor="#bebebe"))
    axes2.add_collection(PatchCollection(patches2, facecolor="#bebebe"))
    axes1.set_title(u"FY4A-%s火点像元位置匹配示意图" % checkVar, fontproperties=font, fontsize=120)

    # 添加本地数据
    x, y = m(slon, slat)
    x2, y2 = m2(slon, slat)
    m.scatter(x[firepixel[:, 2] == 0], y[firepixel[:, 2] == 0], c='red')
    m.scatter(x[firepixel[:, 2] == 1], y[firepixel[:, 2] == 1], c='yellow')
    m2.scatter(x2[firepixel[:, 2] == 0], y2[firepixel[:, 2] == 0], c='red')
    m2.scatter(x2[firepixel[:, 2] == 1], y2[firepixel[:, 2] == 1], c='yellow')
    # 绘制FY3B范围的矩形框
    po = matplotlib.patches.Polygon([[xmax, ymin], [xmin, ymin], [xmin, ymax], [xmax, ymax]], fill=False,
                                    edgecolor="#00ffff", linewidth=4)
    axes1.add_patch(po)
    # 加气象局logo
    axicon = fig.add_axes([0.605, 0.355, 0.09, 0.04])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/%s/logo.png' % checkVar), origin='upper')
    axicon.axis('off')
    # 添加时间
    axtime = fig.add_axes([0.007, 0.368, 0.135, 0.023])
    axtime.text(0.04, 0.1, "FY4A Time:%s\n%s Time:%s" % (
        fy4time.__format__('%Y-%m-%d %H:%M'), checkVar, fy3time.__format__('%Y-%m-%d %H:%M')), fontsize=64)
    axtime.set_xticks([])
    axtime.set_yticks([])
    # 添加图例
    axlegend = fig.add_axes([0.005, 0.005, 0.07, 0.05])
    axlegend.imshow(plt.imread('/FY4APGSQCS/QCS/src/FHS/FHS_FY3B_FY4A/Logo/%s/Lengend-PixelSketch.png' % checkVar), origin='upper')
    axlegend.axis('off')
    # 保存图片
    fig.savefig(output, bbox_inches='tight')
    plt.close()
    return 0


def LogXML(logMsg, logType, xmlfilePath, inputFilename, checkVar, fyTime=True, metric=True, names=True):
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
        nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['pixelCondMap']))))
        nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeOutputFilename2 = doc.createElement("OutputFilename")
        nodeOutputFilename2.appendChild(doc.createTextNode(str(os.path.basename(names['areaSketchMap']))))
        nodeOutputFile.appendChild(nodeOutputFilename2)

        nodeOutputFilename3 = doc.createElement("OutputFilename")
        nodeOutputFilename3.appendChild(doc.createTextNode(str(os.path.basename(names['pixelSketchMap']))))
        nodeOutputFile.appendChild(nodeOutputFilename3)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', 'FY4A')
        nodeDeviationstatic.setAttribute('checkVariety', checkVar)
        nodeDeviationstatic.setAttribute('dataFormat', '15mm')
        nodeDeviationstatic.setAttribute('dataType', 'NOM')
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', 'FHS')
        metricStr = "{Region:%s, MildPixel:%s, ModeratePixel:%s, SeverePixel:%s,CloudPixel:%s}" % (
            metric[0], metric[1], metric[2], metric[3], metric[4])
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