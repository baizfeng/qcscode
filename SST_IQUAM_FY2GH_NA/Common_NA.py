#!/usr/bin/env python
# _*_ coding:UTF-8 _*_
# @Time:2018/12/2012:43
# @Author:baizhaofeng
# @File:Common_NA.py

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

        outTXTnode = root.getElementsByTagName('outputTXTPath')[0]
        xmlnodes['txtPath'] = outTXTnode.childNodes[0].nodeValue

        platTypenode = root.getElementsByTagName('platformtype')[0]
        xmlnodes['platformtype'] = platTypenode.childNodes[0].nodeValue

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
    fy4_time = datetime.datetime.strptime(filename[17:30], '%Y%m%d_%H%M')  # 时间
    dataName = filename[5:8]  # 数据名称
    dataSour = filename[0:4]  # 数据源名称
    baseInfo = {"dataName": dataName,
                "dataSour": dataSour,
                "time": fy4_time,
                "filename": filename[:-4]}
    return baseInfo


def Name(filename, fytime, plattype, outputpath, outTXTpath, prjType,dataSour):
    if (prjType == 'NOM'):
        prj = "MAPN"
    elif (prjType == "GLL"):
        prj = "MAPG"

    shortName = Plattypeprocess(plattype)
    if (len(shortName) == 8):
        scattername = "%s%s_QCS_IQUAM_SCAT_ORIG_NA.png" % (outputpath, filename)
        histname = "%s%s_QCS_IQUAM_HIST_BIAS_NA.png" % (outputpath, filename)
        iquamname = "%s%s_QCS_IQUAM_%s_NA.png" % (outputpath, filename, prj)
        scattertitle = "%s_IQUAM_SST %s" % (dataSour,fytime.__format__('%Y-%m-%d %H:%M'))
        histtitle = "%s_IQUAM_SST %s" % (dataSour,fytime.__format__('%Y-%m-%d %H:%M'))
        iquamtitle = "IQUAM_SST_%s %s-01 0000--%s" % (prj,
                                                      fytime.__format__('%Y-%m'),
                                                      fytime.__format__('%Y-%m-%d %H%M'))  # 标题
        # 带标签中间文件
        txtnamelabel = "%s%s_QCS_IQUAM(FULL)_DATA_BIAS_NA.txt" % (outTXTpath, filename[:-7])
        txtname = "%s%s_QCS_IQUAM_DATA_BIAS_NA.txt" % (outTXTpath, filename[:-7])  # 全类型中间文件
    else:
        platshortStr = '-'.join(shortName)
        scattername = "%s%s_QCS_IQUAM-%s_SCAT_ORIG_NA.png" % (outputpath, filename, platshortStr)
        histname = "%s%s_QCS_IQUAM-%s_HIST_BIAS_NA.png" % (outputpath, filename, platshortStr)
        iquamname = "%s%s_QCS_IQUAM-%s_%s_NA.png" % (outputpath, filename, platshortStr, prj)
        scattertitle = "%s_IQUAM(%s)_SST %s" % (dataSour,platshortStr, fytime.__format__('%Y-%m-%d %H:%M'))
        histtitle = "%s_IQUAM(%s)_SST %s" % (dataSour,platshortStr, fytime.__format__('%Y-%m-%d %H:%M'))
        iquamtitle = "IQUAM(%s)_SST_%s %s-01 0000--%s" % (
            platshortStr, prj, fytime.__format__('%Y-%m'), fytime.__format__('%Y-%m-%d %H%M'))  # 标题
        # 带标签中间文件
        txtnamelabel = "%s%s_QCS_IQUAM-%s_DATA_BIAS_NA.txt" % (outTXTpath, filename[:-7], platshortStr)
        txtname = "%s%s_QCS_IQUAM_DATA_BIAS_NA.txt" % (outTXTpath, filename[:-7])  # 全类型中间文件
    names = {'scatname': scattername,
             'histname': histname,
             'iquamname': iquamname,
             'scattitle': scattertitle,
             'histtitle': histtitle,
             'iquamtitle': iquamtitle,
             'txtnamelabel': txtnamelabel,
             'txtname': txtname}
    return names


def CoarseTimeMatch(path, base_time):
    """
    时间匹配：根据FY4时间得到检验数据
    :param path:
    :param base_time:
    :return:
    """
    exterPath = "%s%s/" % (path, base_time.__format__('%Y'))
    files = os.listdir(exterPath)
    base_time2str = base_time.__format__('%Y%m')
    for file in files:
        timeStr = file.split('-')[0]
        if (timeStr == base_time2str):
            mfile = file
    try:
        mfile
    except NameError:
        var_exists = False
    else:
        var_exists = True
    if (not var_exists):
        raise NameError("未找到检验源数据")
    return exterPath + mfile


def Read_HData(filepath, dataSour):
    """
    读取数据，根据变量名
    :param filepath:
    :return: SST数据,质量标识[0:good,1:mid,2:bad]，array类型
    """
    # 读取FY4数据
    with Dataset(filepath, 'r') as fid:
        fy4_data = fid.variables["%s SST 3-Hour Product" % dataSour][:]
        dqf = fid.variables["SST_Quality"][:]
    if (isinstance(fy4_data, ma.MaskedArray)):
        fy4_data = fy4_data.data
    if (isinstance(dqf, ma.MaskedArray)):
        dqf = dqf.data
    fy4_data[(fy4_data > 40) | (fy4_data < -2)] = -999
    dqf[(dqf == 254) | (dqf == 255)] = 127.
    return fy4_data, dqf


def Read_GData(filepath, dataSour):
    """
    读取数据，根据变量名
    :param filepath:
    :return: SST数据,质量标识[0:good,1:mid,2:bad]，array类型
    """
    # 读取FY4数据
    with Dataset(filepath, 'r') as fid:
        fy4_data = fid.variables["%s SST 3-Hour Product" % dataSour][:]
    if (isinstance(fy4_data, ma.MaskedArray)):
        fy4_data = fy4_data.data
    fy4_data[(fy4_data > 40) | (fy4_data < -2)] = -999
    return fy4_data


def Read_LUT(prjType, filepath, dataSour):
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
        if (dataSour == 'FY2F'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/fy2f2288.dat'
        elif (dataSour == 'FY2G'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/fy2g2288.dat'
        elif (dataSour == 'FY2H'):
            lutpath = r'/FY4APGSQCS/QCS/pyFY4AL2src/LUT/fy2h2288.dat'

        pic = np.zeros((2, 2288, 2288))
        with open(lutpath, "rb") as f:
            for b in range(2):
                for i in range(2288):
                    for j in range(2288):
                        data = f.read(4)
                        elem = struct.unpack("f", data)[0]
                        pic[b][i][j] = elem
        if (dataSour == 'FY2H'):
            lon = pic[0] + 79
        else:
            lon = pic[0]
        lat = pic[1]
    return lat, lon


def ReadTXTDataH(filepath):
    with open(filepath, 'r') as r:
        lines = r.readlines()

    lsts = []
    for l in lines:
        lst = l.split(' ')
        if len(lst) == 7:
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


def ReadTXTDataG(filepath):
    with open(filepath, 'r') as r:
        lines = r.readlines()

    lsts = []
    for l in lines:
        lst = l.split(' ')
        if len(lst) == 6:
            lsts.append(lst[1:6])
    arr = np.array(lsts, dtype=np.float)
    arr = np.unique(arr, axis=0)
    exSST = arr[:, 0]
    fySST = arr[:, 1]
    plon = arr[:, 2]
    plat = arr[:, 3]
    ptype = arr[:, 4]
    bias = fySST[np.where(fySST != -999)] - exSST[np.where(fySST != -999)]

    mdat = {"exSST": exSST, "fySST": fySST, "bias": bias, "plon": plon, "plat": plat, "ptype": ptype}
    return mdat


def ArraytoTime(yearmonth, day, hour, minute):
    """
    输入日、时、分数组，转换为时间类型list
    :param yearmonth:
    :param day:
    :param hour:
    :param minute:
    :return:
    """
    times = np.c_[day, hour, minute]
    times = times.astype(np.str)
    extertimes = []
    for time in times:
        timestr = "%s-%s" % (str(yearmonth), '-'.join(time))
        extertimes.append(datetime.datetime.strptime(timestr, '%Y%m-%d-%H-%M'))
    return np.array(extertimes)


def Read_exterData(exterfilepath):
    """
    读取检验源数据信息
    :param exterfilepath: 检验源数据路径
    :return: 检验源数据信息，字典类型
    """
    filename = os.path.basename(exterfilepath)

    yearMonth = filename.split('-')[0]
    with Dataset(exterfilepath, 'r') as fid:
        lat = fid.variables["lat"][:]
        lon = fid.variables["lon"][:]

        day = fid.variables["day"][:]
        hour = fid.variables["hour"][:]
        minute = fid.variables["minute"][:]

        platform = fid.variables["platform_type"][:]  # 获取数据平台类型
        level = fid.variables["quality_level"][:]  # 质量等级 best:5
        sst = fid.variables["sst"][:] - 273.15  # 将K转换为℃

    exterTime = ArraytoTime(yearMonth, day, hour, minute)

    exterInfo = {"lat": lat,
                 "lon": lon,
                 "times": exterTime,
                 "platform": platform,
                 "level": level,
                 "sst": sst}
    return exterInfo


def FineTimeMatch(base_time, exterTime, threshold, exterinfo):
    """
    时间精匹配：找出时间阈值范围内的检验源
    :param base_time:FY4时间
    :param exterTime:检验源时间，array
    :param threshold:时间阈值
    :param exterinfo:检验源数据
    :return:
    """
    index = np.where(np.abs(exterTime - base_time) <= datetime.timedelta(minutes=int(threshold)))
    indexlevel = np.where(exterinfo['level'][index] == 5)
    extermatched = {"lat": exterinfo['lat'][index][indexlevel],
                    "lon": exterinfo['lon'][index][indexlevel],
                    "platform": exterinfo['platform'][index][indexlevel],
                    "sst": exterinfo['sst'][index][indexlevel]}
    return extermatched


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


def InterpH(tlat, tlon, slat, slon, sdata, dqf, spaceThreshold):
    """
    空间插值
    :param tlat: 目标纬度
    :param tlon: 目标精度
    :param slat: 待插值纬度
    :param slon: 待插值精度
    :param sdata: 待插值数据
    :param spaceThreshold: 空间阈值
    :return: 插值后的数据，维度和目标纬度一致
    """
    if (np.ndim(slat) == 1):
        slon, slat = np.meshgrid(slon, slat)

    slat[np.where(slat == 99999.)] = np.nan
    slon[np.where(slon == 99999.)] = np.nan

    # 2.转换经纬度到笛卡尔坐标系
    # xs, ys, zs = lon_lat_to_cartesian(slon[~np.isnan(slon)], slat[~np.isnan(slat)])
    # xt, yt, zt = lon_lat_to_cartesian(tlon, tlat)

    tree = spatial.cKDTree(zip(slon[~np.isnan(slon)], slat[~np.isnan(slat)]))  # 先将待插值数据建树，latlon不能有NA值
    d, inds = tree.query(zip(tlon, tlat), k=1, distance_upper_bound=float(spaceThreshold))  # 4km

    nearest = -999 * np.ones(tlat.size)  # 阈值之外的赋值为-999
    nearestdqf = 127 * np.ones(tlat.size)
    nearest[~np.isinf(d)] = sdata[~np.isnan(slat)][inds[~np.isinf(d)]]
    nearestdqf[~np.isinf(d)] = dqf[~np.isnan(slat)][inds[~np.isinf(d)]]
    return nearest, nearestdqf


def InterpG(tlat, tlon, slat, slon, sdata, spaceThreshold):
    """
    空间插值
    :param tlat: 目标纬度
    :param tlon: 目标精度
    :param slat: 待插值纬度
    :param slon: 待插值精度
    :param sdata: 待插值数据
    :param spaceThreshold: 空间阈值
    :return: 插值后的数据，维度和目标纬度一致
    """
    if (np.ndim(slat) == 1):
        slon, slat = np.meshgrid(slon, slat)

    slat[np.where(slat == 99999.)] = np.nan
    slon[np.where(slon == 99999.)] = np.nan

    # 2.转换经纬度到笛卡尔坐标系
    # xs, ys, zs = lon_lat_to_cartesian(slon[~np.isnan(slon)], slat[~np.isnan(slat)])
    # xt, yt, zt = lon_lat_to_cartesian(tlon, tlat)

    tree = spatial.cKDTree(zip(slon[~np.isnan(slon)], slat[~np.isnan(slat)]))  # 先将待插值数据建树，latlon不能有NA值
    d, inds = tree.query(zip(tlon, tlat), k=1, distance_upper_bound=float(spaceThreshold))  # 4km

    nearest = -999 * np.ones(tlat.size)  # 阈值之外的赋值为-999
    nearest[~np.isinf(d)] = sdata[~np.isnan(slat)][inds[~np.isinf(d)]]
    return nearest


def Strprocess(mystr):
    """
    字符串处理，将输入的数据获取途径，替换为对应的ID
    :param mystr:
    :return:list类型
    """
    s = str(mystr).strip().lower().split(',')
    pattern = {'ship': 1, 'drifting buoy': 2, 'tropical moored buoy': 3, 'coastal moored buoy': 4, 'argo float': 5,
               'high resolution drifter': 6, 'imos': 7, 'crw buoy': 8}
    rep = [pattern[x] if x in pattern else x for x in s]
    return rep


def Plattypeprocess(mystr):
    """
    字符串处理，将输入的数据获取途径，替换为对应的简写
    :param mystr:
    :return:list类型
    """
    s = str(mystr).strip().upper().split(',')
    pattern1 = {'SHIP': 'SHIP', 'DRIFTING BUOY': 'DRIFTBUOY', 'TROPICAL MOORED BUOY': 'TROPBUOY',
                'COASTAL MOORED BUOY': 'COASTBUOY', 'ARGO FLOAT': 'ARGO',
                'HIGH RESOLUTION DRIFTER': 'HIGHGRESDRIFT', 'IMOS': 'IMOS', 'CRW BUOY': 'CRWBUOY'}
    rep = [pattern1[x] if x in pattern1 else x for x in s]
    return rep


def Condi(platID, plattype):
    if (len(platID) == 1):
        platindex = np.where(plattype == platID[0])
    elif (len(platID) == 2):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]))
    elif (len(platID) == 3):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]) | (plattype == platID[2]))
    elif (len(platID) == 4):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]) | (plattype == platID[2]) |
                             (plattype == platID[3]))
    elif (len(platID) == 5):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]) | (plattype == platID[2]) |
                             (plattype == platID[3]) | (plattype == platID[4]))
    elif (len(platID) == 6):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]) | (plattype == platID[2]) |
                             (plattype == platID[3]) | (plattype == platID[4]) | (plattype == platID[5]))
    elif (len(platID) == 7):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]) | (plattype == platID[2]) |
                             (plattype == platID[3]) | (plattype == platID[4]) | (plattype == platID[5]) |
                             (plattype == platID[6]))
    elif (len(platID) == 8):
        platindex = np.where((plattype == platID[0]) | (plattype == platID[1]) | (plattype == platID[2]) |
                             (plattype == platID[3]) | (plattype == platID[4]) | (plattype == platID[5]) |
                             (plattype == platID[6]) | (plattype == platID[7]))
    return platindex


def SpaceMatchH(fy4, fy4_dqf, exterMatched, platformtype):
    plattype = exterMatched['platform']
    iquam = exterMatched['sst']
    lon = exterMatched['lon']
    lat = exterMatched['lat']

    platID = Strprocess(platformtype)
    platindex = Condi(platID, plattype)

    exSST = iquam[platindex]
    fySST = fy4[platindex]
    fyDQF = fy4_dqf[platindex]
    plon = lon[platindex]
    plat = lat[platindex]
    ptype = plattype[platindex]

    validindex = np.where((fySST != -999.0) & (fyDQF != 127.0))

    leng = len(validindex[0])
    if (leng == 0):
        raise IndexError("Though find iquam data in the special time threshold,but can't find valid region")

    exSST = exSST[validindex].reshape(leng, 1)
    fySST = fySST[validindex].reshape(leng, 1)
    fyDQF = fyDQF[validindex].reshape(leng, 1)
    plon = plon[validindex].reshape(leng, 1)
    plat = plat[validindex].reshape(leng, 1)
    ptype = ptype[validindex].reshape(leng, 1)
    bias = fySST - exSST

    mdat = {'exSST': exSST,
            'fySST': fySST,
            'fyDQF': fyDQF,
            'bias': bias,
            'plon': plon,
            'plat': plat,
            'ptype': ptype}

    return mdat


def WriteToTXTH(mdat, fytime, filepath):
    """
    输出为TXT文件
    :param exSST: 经过和fy[时间匹配]指定[platform Type]检验源数据
    :param fySST: 经过和指定[platform Type]ex数据进行[空间匹配]的fy数据,去除了无效值
    :param fyDQF: 经过和指定[platform Type]ex数据进行[空间匹配]的fy dqf数据,去除了无效值
    :param plon: 经度
    :param plat: 纬度
    :return:
    """
    np.set_printoptions(suppress=True)
    timestr = fytime.__format__('%Y%m%d%H%M%S')
    leng = mdat['fySST'].shape[0]
    vtime = np.array([timestr] * leng).reshape(leng, 1)
    newarr = np.hstack((vtime, mdat['exSST'], mdat['fySST'], mdat['fyDQF'], mdat['plon'], mdat['plat'], mdat['ptype']))
    newarr = np.unique(newarr, axis=0)
    if (newarr.shape[1] == 7):
        # 判断文件是否存在，若不存在
        if not os.path.exists(filepath):
            with open(filepath, "w") as f:
                for i in range(newarr.shape[0]):
                    write_str = '%s %s %s %s %s %s %s\n' % (
                        newarr[i, 0], newarr[i, 1], newarr[i, 2], newarr[i, 3], newarr[i, 4], newarr[i, 5],
                        newarr[i, 6])
                    f.write(write_str)
        else:
            with open(filepath, "a") as f:
                for i in range(newarr.shape[0]):
                    write_str = '%s %s %s %s %s %s %s\n' % (
                        newarr[i, 0], newarr[i, 1], newarr[i, 2], newarr[i, 3], newarr[i, 4], newarr[i, 5],
                        newarr[i, 6])
                    f.write(write_str)


def CalcateH(mdat):
    exSST = mdat['exSST']
    fySST = mdat['fySST']
    fyDQF = mdat['fyDQF']
    metrics = collections.OrderedDict()
    bias0 = fySST[np.where(fyDQF == 0)] - exSST[np.where(fyDQF == 0)]
    bias1 = fySST[np.where(fyDQF == 1)] - exSST[np.where(fyDQF == 1)]
    bias2 = fySST[np.where(fyDQF == 2)] - exSST[np.where(fyDQF == 2)]
    bias3 = mdat['bias']

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
        slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(fySST.flatten(), exSST.flatten())
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
    metrics['QualID0_CORR'] = round(r_value0,4)

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
    metrics['QualID1_CORR'] = round(r_value1,4)

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
    metrics['QualID2_CORR'] = round(r_value2,4)

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
    metrics['QualID3_CORR'] = round(r_value3,4)

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


def SpaceMatchG(fy4, exterMatched, platformtype):
    plattype = exterMatched['platform']
    iquam = exterMatched['sst']
    lon = exterMatched['lon']
    lat = exterMatched['lat']

    platID = Strprocess(platformtype)
    platindex = Condi(platID, plattype)

    exSST = iquam[platindex]
    fySST = fy4[platindex]
    plon = lon[platindex]
    plat = lat[platindex]
    ptype = plattype[platindex]

    validindex = np.where(fySST != -999.0)

    leng = len(validindex[0])
    if (leng == 0):
        raise IndexError("Though find iquam data in the special time threshold,but can't find valid region")

    exSST = exSST[validindex].reshape(leng, 1)
    fySST = fySST[validindex].reshape(leng, 1)
    plon = plon[validindex].reshape(leng, 1)
    plat = plat[validindex].reshape(leng, 1)
    ptype = ptype[validindex].reshape(leng, 1)
    bias = fySST - exSST

    mdat = {'exSST': exSST,
            'fySST': fySST,
            'bias': bias,
            'plon': plon,
            'plat': plat,
            'ptype': ptype}

    return mdat


def WriteToTXTG(mdat, fytime, filepath):
    """
    输出为TXT文件
    :param exSST: 经过和fy[时间匹配]指定[platform Type]检验源数据
    :param fySST: 经过和指定[platform Type]ex数据进行[空间匹配]的fy数据,去除了无效值
    :param fyDQF: 经过和指定[platform Type]ex数据进行[空间匹配]的fy dqf数据,去除了无效值
    :param plon: 经度
    :param plat: 纬度
    :return:
    """
    np.set_printoptions(suppress=True)
    timestr = fytime.__format__('%Y%m%d%H%M')
    leng = mdat['fySST'].shape[0]
    vtime = np.array([timestr] * leng).reshape(leng, 1)
    newarr = np.hstack((vtime, mdat['exSST'], mdat['fySST'], mdat['plon'], mdat['plat'], mdat['ptype']))
    newarr = np.unique(newarr, axis=0)
    if (newarr.shape[1] == 6):
        # 判断文件是否存在，若不存在
        if not os.path.exists(filepath):
            with open(filepath, "w") as f:
                for i in range(newarr.shape[0]):
                    write_str = '%s %s %s %s %s %s\n' % (
                        newarr[i, 0], newarr[i, 1], newarr[i, 2], newarr[i, 3], newarr[i, 4], newarr[i, 5])
                    f.write(write_str)
        else:
            with open(filepath, "a") as f:
                for i in range(newarr.shape[0]):
                    write_str = '%s %s %s %s %s %s\n' % (
                        newarr[i, 0], newarr[i, 1], newarr[i, 2], newarr[i, 3], newarr[i, 4], newarr[i, 5])
                    f.write(write_str)


def CalcateG(mdat):
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
    std3 = round(np.sqrt(np.square(bias3).sum() / (num3-1)), 4)  # 标准差
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
    metrics['QualID3_CORR'] = round(r_value3,4)

    param = {"slope3": slope3,
             "intercept3": intercept3,
             "QualID3_CORR": r_value3}

    return metrics, param


def DrawScatterH(mdat, param, totalNum, names):
    exSST = mdat['exSST']
    fySST = mdat['fySST']
    fyDQF = mdat['fyDQF']

    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)

    fig.suptitle(names['scattitle'], fontsize=14, fontweight='bold')
    H0 = ax.scatter(fySST[np.where(fyDQF == 0)], exSST[np.where(fyDQF == 0)], c='r', marker='o')
    H1 = ax.scatter(fySST[np.where(fyDQF == 1)], exSST[np.where(fyDQF == 1)], c='g', marker='D')
    H2 = ax.scatter(fySST[np.where(fyDQF == 2)], exSST[np.where(fyDQF == 2)], c='b', marker='x')

    ax.set_xlabel('FY2H SST'+'(' + u'\u2103' + ')', fontsize=12)
    ax.set_ylabel('IQUAM SST'+'(' + u'\u2103' + ')', fontsize=12)

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


def DrawScatterG(mdat, param, totalNum, names):
    exSST = mdat['exSST']
    fySST = mdat['fySST']

    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)

    fig.suptitle(names['scattitle'], fontsize=14, fontweight='bold')
    H0 = ax.scatter(fySST, exSST, c='r', marker='o')

    ax.set_xlabel('FY2G SST'+'(' + u'\u2103' + ')', fontsize=12)
    ax.set_ylabel('IQUAM SST'+'(' + u'\u2103' + ')', fontsize=12)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.grid(True, linestyle='dashed')

    ax.plot([-5, 45], [-5, 45], color='black', zorder=3)  # 画参考线
    plt.xlim(-5, 45)
    plt.ylim(-5, 45)
    # ax.legend([H0], ['Best', 'Good', 'Bad'], loc="lower right")

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


def DrawHist(bias, metric, names, dataSour):
    # 创建直方图
    # 第一个参数为待绘制的定量数据
    # 第二个参数为划分的区间个数
    fig, ax = plt.subplots(figsize=(8, 5), dpi=100)  # 输出图片的大小
    fig.suptitle(names['histtitle'], fontsize=14, fontweight='bold')
    # n, bins, patches = ax.hist(bias, bins=20, density=True, edgecolor='k', facecolor='white', zorder=2)  # bins是分为20个柱子
    # if (bias.shape[0] > 5):
    #     # add a 'best fit' line
    #     y = mlab.normpdf(bins, metric['QualID3_MEAN'], metric['QualID3_STD'])
    #     l = ax.plot(bins, y, '--', color='black', linewidth=1)

    # kde=True表示是否显示拟合曲线，如果为False则只出现直方图,kde_kws（拟合曲线的设置）、hist_kws（直方柱子的设置）
    sns.distplot(bias, bins=40, hist=True, kde=True,
                 kde_kws={"color": "black", "lw": 1, 'linestyle': '--'},
                 hist_kws={"histtype": 'bar', "edgecolor": 'black', "facecolor": 'white', "lw": 1})

    plt.xlim(-6, 4)
    ax.set_xlabel(dataSour + '-IQUAM(' + u'\u2103' + ')', fontsize=12)
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
                  fontsize=20,
                  fontweight='bold',pad=10)
    #添加Logo
    axicon=fig.add_axes([0.85,0.01,0.15,0.05])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'),origin='upper')
    axicon.axis('off')	

    fig.savefig(names['iquamname'])
    plt.close()
    return 0


def LogXML(logMsg, logType, nodes, dataSour, fyTime=True, metric=True, names=True):
    xmlfilePath = nodes['xmlPath']
    inputFilename = nodes['filepath']
    plattype = nodes['platformtype']
    warnValue = float(nodes['warnvalue'])

    shortName = Plattypeprocess(plattype)
    if (len(shortName) == 8):
        checkVar = 'IQUAM'
    else:
        checkVar = "IQUAM-%s" % '-'.join(shortName)
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
        nodeOutputFilename1 = doc.createElement("OutputFilename")
        nodeOutputFilename1.appendChild(doc.createTextNode(str(os.path.basename(names['iquamname']))))
        nodeOutputFile.appendChild(nodeOutputFilename1)

        nodeOutputFilename2 = doc.createElement("OutputFilename")
        nodeOutputFilename2.appendChild(doc.createTextNode(str(os.path.basename(names['scatname']))))
        nodeOutputFile.appendChild(nodeOutputFilename2)

        nodeOutputFilename3 = doc.createElement("OutputFilename")
        nodeOutputFilename3.appendChild(doc.createTextNode(str(os.path.basename(names['histname']))))
        nodeOutputFile.appendChild(nodeOutputFilename3)

        nodeOutputFilename4 = doc.createElement("OutputFilename")
        nodeOutputFilename4.appendChild(doc.createTextNode(str(os.path.basename(names['txtnamelabel']))))
        nodeOutputFile.appendChild(nodeOutputFilename4)

        nodeOutputFilename5 = doc.createElement("OutputFilename")
        nodeOutputFilename5.appendChild(doc.createTextNode(str(os.path.basename(names['txtname']))))
        nodeOutputFile.appendChild(nodeOutputFilename5)

        nodeDeviationstatic = doc.createElement('Deviationstatic')
        nodeDeviationstatic.setAttribute('bechekVariety', dataSour)
        nodeDeviationstatic.setAttribute('checkVariety', checkVar)
        nodeDeviationstatic.setAttribute('dataFormat', '15mm')
        nodeDeviationstatic.setAttribute('dataType', 'NOM')
        nodeDeviationstatic.setAttribute('productDate', fyTime.__format__('%Y%m%d%H%M%S'))
        nodeDeviationstatic.setAttribute('productVariety', dataSour+'#SST-NA')
        metricStr = "{%s}" % re.sub('[()\'\[\]]', '',
                                    str(metric).replace("',", ":").replace("OrderedDict", "")).replace(' ', '')
        nodeDeviationstatic.setAttribute('values', metricStr)
        nodeOutputFile.appendChild(nodeDeviationstatic)

        if metric['QualID3_AE'] > warnValue:
            warntype = 1 		
            nodewarnType = doc.createElement("warningType")
            nodewarnType.appendChild(doc.createTextNode(str(warntype)))
            nodewarnMsg = doc.createElement("warningMsg")
            nodewarnMsg.appendChild(doc.createTextNode("info" if not warntype else "精度预警：FY2 SST产品绝对偏差超过%s ℃" % warnValue))

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