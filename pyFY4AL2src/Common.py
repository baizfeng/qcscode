#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/8/3 14:59
# @Author  : baizhaofeng
# @File    : Common.py

import collections
from xml.dom import minidom
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import datetime
import os
from DrawMap import *
from mpl_toolkits.basemap import Basemap
from PIL import Image


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
        xmlnodes['inputFileName'] = filenode.childNodes[0].nodeValue

        dsNamenode = root.getElementsByTagName('identify')[0]
        xmlnodes['varName'] = dsNamenode.childNodes[0].nodeValue

        outputnode = root.getElementsByTagName('outputPath')[0]
        xmlnodes['pngPath'] = outputnode.childNodes[0].nodeValue

        prjTypenode = root.getElementsByTagName('type')[0]
        xmlnodes['prjType'] = prjTypenode.childNodes[0].nodeValue

        outputxmlnode = root.getElementsByTagName('outputXml')[0]
        xmlnodes['xmlPath'] = outputxmlnode.childNodes[0].nodeValue

        endfilenode = root.getElementsByTagName('end_fileName')[0]
        xmlnodes['endFileName'] = endfilenode.childNodes[0].nodeValue

    return xmlnodes


def Var(varName, fy4time, fy4Endtime):
    dhour = ((fy4Endtime - fy4time).seconds + 1) / 3600.
    if (varName == 'Precipitation'):
        if (dhour == 0.25):
            varName = 'Precipitation_0.25h'
        elif (dhour == 1):
            varName = 'Precipitation_1h'
        elif (dhour == 3):
            varName = 'Precipitation_3h'
        elif (dhour == 6):
            varName = 'Precipitation_6h'
        elif (dhour == 24):
            varName = 'Precipitation_24h'
    else:
        varName = varName
    return varName


def GetBaseInfo(inputFileName, varName):
    """
    得到数据基本信息：时间、数据名称、数据源名称、分辨率、数据单位
    :param filename:
    :return:
    """
    filename = os.path.basename(inputFileName)
    fy4_time = datetime.datetime.strptime(filename.split('_')[9], '%Y%m%d%H%M%S')  # 时间
    fy4_endtime = datetime.datetime.strptime(filename.split('_')[10], '%Y%m%d%H%M%S')  # 时间
    dataName = filename[30:33]  # 数据名称
    dataSour = filename[0:4]  # 数据源名称
    dataRes = filename[74:79]

    varName = Var(varName, fy4_time, fy4_endtime)

    baseInfo = {"dataName": dataName,
                "dataSour": dataSour,
                "dataRes": dataRes,
                "dataTime": fy4_time,
                "dataEndTime": fy4_endtime,
                "dataUnit": GetUnit(varName),
                "filename": filename[:-3]}
    return baseInfo, varName


def GetUnit(varName):
    """
根据产品的变量名，得到数据的单位
    :param varName: 产品中的变量名称
    :return: 数据单位，str类型
    """
    if (varName in ("SST", "SST Monthly Mean Product")):
        unit = "Unit:(" + u'\u2103' + ")"
    elif (varName in ("ACI", "NIR Black sky albedo", "NIR White sky albedo", "SW Black sky albedo",
                      "SW White sky albedo", "VIS Black sky albedo", "VIS White sky albedo", "RDC_Rank",
                      "SNC", "CLM", "CLP", "CLT", "FOG", "CPD_COT", "CPN_COT", "DSD", "AOD", "AMV", "FHS")):
        unit = "Unit:(" + "NULL" + ")"
    elif (varName == 'OLR'):
        unit = "Unit:(" + r'$W/M^2$' + ")"
    elif (varName in ("LSE_8.5um", "LSE_10.8um", "LSE_12.0um")):
        unit = ""
    elif (varName in ('LST', 'NOMChannel07', 'NOMChannel08', 'NOMChannel09', 'NOMChannel10',
                      'NOMChannel11', 'NOMChannel12', 'NOMChannel13', 'NOMChannel14')):
        unit = "Unit:(" + r'$K$' + ")"
    elif (varName in ("LPW_HIGH", "LPW_LOW", "LPW_MID", "TPW")):  # LPW
        unit = "Unit:(" + r'$g/kg$' + ")"
    elif (varName == 'CTH'):
        unit = "Unit:(" + r'$km$' + ")"
    elif (varName == 'TFTP_Z_depth'):
        unit = "Unit:(" + r'$km$' + ")"
    elif (varName == 'SSI'):
        unit = "Unit:(" + r'$W/m^2$' + ")"
    elif (varName == 'CTT'):
        unit = "Unit:(" + r'$k$' + ")"
    elif (varName in ('DLR', 'ULR')):
        unit = "Unit:(" + r'$W/M^2$' + ")"
    elif (varName in ('CPD_CER', 'CPN_CER')):
        unit = "Unit:(" + r'$um$' + ")"
    elif (varName in ('CPD_LWP', 'CPN_LWP')):
        unit = "Unit:(" + r'$g/m^2$' + ")"
    elif (varName in ('CPD_IWP', 'CPN_IWP')):
        unit = "Unit:(" + r'$g/m^2$' + ")"
    elif (varName == 'CTP'):
        unit = "Unit:(" + r'$hPa$' + ")"
    elif (varName == 'RSR'):
        unit = "Unit:(" + r'$W/M^2$' + ")"
    elif (varName in ('Precipitation_0.25h', 'Precipitation_1h', 'Precipitation_3h', 'Precipitation_6h',
                      'Precipitation_24h')):
        unit = "Unit:(" + r'$mm$' + ")"
    return unit


def myfunc(a):
    if len(a) < 4:
        return np.inf
    else:
        return a[3]


def Read_Data(filepath, varName):
    """
    读取数据，根据变量名
    :param filepath:
    :param varName:
    :return: 数据，array类型
    """
    if (varName in ("LSE_8.5um", "LSE_10.8um", "LSE_12.0um")):
        band = varName.split('_')[1]
        varName = varName.split('_')[0]
    if (varName in ('CPD_COT', 'CPD_CER', 'CPD_IWP', 'CPD_LWP', 'CPN_COT', 'CPN_CER', 'CPN_IWP', 'CPN_LWP')):
        varName = varName.split('_')[1]
    if (varName in ('Precipitation_0.25h', 'Precipitation_1h', 'Precipitation_3h', 'Precipitation_6h',
                    'Precipitation_24h')):
        varName = varName.split('_')[0]

    if (varName == "ACI"):
        with Dataset(filepath, 'r') as fid:
            ACI_R = fid.variables["Channel0161"][:]
            ACI_G = fid.variables["Channel0083"][:]
            ACI_B = fid.variables["Channel0065"][:]
        fy4_data = np.dstack((ACI_R, ACI_G, ACI_B))
    elif (varName == "AMV"):
        with Dataset(filepath, 'r') as fid:
            wind_direction = fid.variables["wind_direction"][:]
            wind_speed = fid.variables["wind_speed"][:]
            pressure = fid.variables["pressure"][:]
        fy4_data = np.array([wind_speed, wind_direction, pressure])
    elif (varName == "DSD"):
        with Dataset(filepath, 'r') as fid:
            dst = fid.variables["DST"][:]
            iddi = fid.variables["IDDI_BK"][:]
        if (isinstance(dst, ma.MaskedArray)):
            dst = dst.data
        if (isinstance(iddi, ma.MaskedArray)):
            iddi = iddi.data
        if (dst.dtype == 'uint16'):
            dst.dtype = np.int16
            dst = dst.byteswap()
        fy4_data = np.dstack([dst, iddi])
    elif (varName == "CLM"):
        with Dataset(filepath, 'r') as fid:
            fy4_data = fid.variables["CLM"][:]
			
        with Dataset("/FY4APGSQCS/QCS/pyFY4AL2src/LUT/Land1_Ocean0.NC", 'r') as fid:
            qcbin4 = fid.variables["mask"][:]			
        if (isinstance(fy4_data, ma.MaskedArray)):
            fy4_data = fy4_data.data
			
        fy4_data[(fy4_data == 2) & (qcbin4 == 1)] = 4
        fy4_data[(fy4_data == 2) & (qcbin4 == 0)] = 5
        fy4_data[(fy4_data == 3) & (qcbin4 == 1)] = 6
        fy4_data[(fy4_data == 3) & (qcbin4 == 0)] = 7	
    else:
        # 读取FY4数据
        with Dataset(filepath, 'r') as fid:
            fy4_data = fid.variables[varName][:]
        if ((np.ndim(fy4_data) == 3) & (varName == 'AOD')):
            fy4_data = fy4_data[1]
        elif ((np.ndim(fy4_data) == 3) & (varName == 'LSE')):
            if (band == '8.5um'):
                fy4_data = np.round(fy4_data[0] * 10000, 0)
            elif (band == '10.8um'):
                fy4_data = np.round(fy4_data[1] * 10000, 0)
            elif (band == '12.0um'):
                fy4_data = np.round(fy4_data[2] * 10000, 0)
        if (isinstance(fy4_data, ma.MaskedArray)):
            fy4_data = fy4_data.data
        if (fy4_data.dtype == 'uint16'):
            fy4_data.dtype = np.int16
            fy4_data = fy4_data.byteswap()
        elif (fy4_data.dtype == 'uint32'):
            fy4_data.dtype = np.int32
            fy4_data = fy4_data.byteswap()
        elif ((fy4_data.dtype == 'float32') & (varName != 'LSE')):
            fy4_data = fy4_data.astype(np.int16)  # 数据类型转换
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
    if (prjType in ('GEO', 'NUL')):
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
            lat = fid.variables['LAT'][:][:rows, :]
            lon = fid.variables['LON'][:][:rows, :]
        lat = ma.masked_values(lat, -999.)
        lon = ma.masked_values(lon, -999.)
    return lat, lon


def MaskInValid(fy4_data, varName):
    """
    无效值掩膜
    :param fy4_data:
    :param varName:
    :return: 掩膜数组
    """
    if (varName in ("SST", "SST Monthly Mean Product")):
        fy4_data[(fy4_data >= 65530) | (fy4_data == -888)] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "ACI"):
        fy4_mdata = ma.masked_values(fy4_data, 65535.)
        fy4_mdata = (fy4_mdata * 255.999).astype(np.uint8)
    elif (varName in ("NIR Black sky albedo", "NIR White sky albedo", "SW Black sky albedo",
                      "SW White sky albedo", "VIS Black sky albedo", "VIS White sky albedo")):
        fy4_data[fy4_data > 1000.] = -999.
        fy4_mdata_ = ma.masked_values(fy4_data, -999.)
        fy4_mdata = fy4_mdata_ * 0.001
    elif (varName == "RDC_Rank"):
        fy4_data[fy4_data == 0] = 7
        fy4_mdata = ma.masked_values(fy4_data, 255)
    elif (varName == "OLR"):
        fy4_data[(fy4_data > 450) | (fy4_data < 40)] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "LST"):
        fy4_data[fy4_data > 500] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName in ("LPW_HIGH", "LPW_LOW", "LPW_MID", "TPW")):
        fy4_data[(fy4_data < 65534) & (fy4_data > 10)] = -999
        fy4_data[fy4_data < 0] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName == "SNC"):
        fy4_mdata = ma.masked_values(fy4_data, 243)
    elif (varName in ('NOMChannel07', 'NOMChannel08', 'NOMChannel09', 'NOMChannel10',
                      'NOMChannel11', 'NOMChannel12', 'NOMChannel13', 'NOMChannel14')):
        fy4_data[(fy4_data > 400) | (fy4_data < 0)] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "CLM"):
        fy4_data[fy4_data > 7] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "CLP"):
        fy4_data[fy4_data > 5] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "CLT"):
        fy4_data[fy4_data > 9] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "CTH"):
        fy4_mdata = ma.masked_values(fy4_data, 65535)
        fy4_mdata = fy4_mdata / 1000.0
        fy4_mdata[(fy4_mdata > 27) | (fy4_mdata < 0)] = -999
    elif (varName == "TFTP_Z_depth"):
        fy4_data[(fy4_data > 1) | (fy4_data < 0)] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "SSI"):
        fy4_data[fy4_data == 65535] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "CTT"):
        fy4_data[(fy4_data < 160) | (fy4_data > 320)] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName in ("DLR", "ULR")):
        fy4_data[(fy4_data < 32761) & (fy4_data > 500)] = -999
        fy4_data[fy4_data < 50] = -999
        fy4_mdata = ma.masked_values(fy4_data, 32766)
    elif (varName == "FOG"):
        fy4_data[fy4_data == 65535] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName in ("CPD_CER", "CPN_CER")):
        fy4_data[(fy4_data < 65527) & (fy4_data > 100)] = -999
        fy4_data[fy4_data < 0] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName in ("CPD_COT", "CPN_COT")):
        fy4_data[(fy4_data < 65527) & (fy4_data > 150)] = -999
        fy4_data[fy4_data < 0] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName in ("CPD_LWP", "CPN_LWP")):
        fy4_data[(fy4_data < 65527) & (fy4_data > 1000)] = -999
        fy4_data[fy4_data < 0] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName in ("CPD_IWP", "CPN_IWP")):
        fy4_data[(fy4_data < 65527) & (fy4_data > 1500)] = -999
        fy4_data[fy4_data < 0] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName == "CTP"):
        fy4_data[(fy4_data < 0) | (fy4_data > 1100)] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName == "RSR"):
        fy4_data[fy4_data == 65535] = -999
        fy4_mdata = ma.masked_values(fy4_data, -999)
    elif (varName == "DSD"):
        dst = fy4_data[:, :, 0]
        iddi = fy4_data[:, :, 1]
        mdst = ma.masked_values(dst, 32766)
        middi = ma.masked_values(iddi, 65535)
        middi = ma.masked_values(middi, 65534)
        fy4_mdata = np.array([mdst, middi])
    elif (varName == "FHS"):
        fy4_mdata = ma.masked_values(fy4_data, 255)
    elif (varName == "Precipitation_24h"):
        fy4_mdata = ma.masked_values(fy4_data, 65535.)
        fy4_mdata = ma.masked_values(fy4_mdata, 0.)
        fy4_mdata[fy4_mdata > 100] = 100
    elif (varName == "Precipitation_6h"):
        fy4_mdata = ma.masked_values(fy4_data, 65535.)
        fy4_mdata = ma.masked_values(fy4_mdata, 0.)
        fy4_mdata[fy4_mdata > 50] = 50
    elif (varName in ("Precipitation_3h", "Precipitation_1h", "Precipitation_0.25h")):
        fy4_mdata = ma.masked_values(fy4_data, 65535.)
        fy4_mdata = ma.masked_values(fy4_mdata, 0.)
        fy4_mdata[fy4_mdata > 20] = 20
    elif (varName == "AOD"):
        fy4_data[(fy4_data < 65530) & (fy4_data > 2)] = -999
        fy4_data[fy4_data < 0] = -999
        fy4_mdata = ma.masked_values(fy4_data, 65535)
    elif (varName in ("LSE_8.5um", "LSE_10.8um", "LSE_12.0um")):
        fy4_data[(fy4_data > 10000) | (fy4_data < 0)] = -999
        fy4_mdata_ = ma.masked_values(fy4_data, -999)
        fy4_mdata = fy4_mdata_ * 0.0001
    elif (varName == "AMV"):
        wind_speed = fy4_data[0]
        wind_direction = fy4_data[1]
        pressure = fy4_data[2]
        wind_mdirection = ma.masked_values(wind_direction, -999)
        wind_mspeed = ma.masked_values(wind_speed, -999)
        mpressure = ma.masked_values(pressure, -999)
        fy4_mdata = np.array([wind_mspeed, wind_mdirection, mpressure])
    return fy4_mdata


def MapTitle(baseInfo, varName, prjType):
    """
    根据文件基本信息，变量，投影，输出产品图标题
    :param baseInfo:
    :param varName:
    :param prjType:
    :return: title,str类型
    """
    if (varName == baseInfo['dataName']):
        if (varName == 'AMV'):
            maintitle = "%s_%s_%s_%s_%s_%s" % (
                baseInfo['dataSour'], baseInfo['dataName'], baseInfo['dataTime'].__format__('%Y%m%d_%H%M%S'),
                baseInfo['dataRes'], baseInfo['filename'][35:39], prjType)
        else:
            maintitle = "%s_%s_%s_%s_%s" % (
                baseInfo['dataSour'], baseInfo['dataName'], baseInfo['dataTime'].__format__('%Y%m%d_%H%M%S'),
                baseInfo['dataRes'], prjType)
    else:
        if (varName in ('CPD_COT', 'CPD_CER', 'CPD_IWP', 'CPD_LWP', 'CPN_COT', 'CPN_CER', 'CPN_IWP', 'CPN_LWP')):
            varName = varName.split('_')[1]
        if (varName in ("LSE_8.5um", "LSE_10.8um", "LSE_12.0um")):
            varName = varName.split('_')[1]
        if (varName in ('NOMChannel07', 'NOMChannel08', 'NOMChannel09', 'NOMChannel10',
                        'NOMChannel11', 'NOMChannel12', 'NOMChannel13', 'NOMChannel14')):
            dict1 = {'NOMChannel07': '3.72um(high)', 'NOMChannel08': '3.72um(low)', 'NOMChannel09': '6.25um',
                     'NOMChannel10': '7.10um', 'NOMChannel11': '8.50um', 'NOMChannel12': '10.08um',
                     'NOMChannel13': '12um', 'NOMChannel14': '13.5um'}
            varName = dict1[varName]

        maintitle = "%s_%s_%s_%s_%s_%s" % (
            baseInfo['dataSour'], baseInfo['dataName'], baseInfo['dataTime'].__format__('%Y%m%d_%H%M%S'),
            baseInfo['dataRes'], varName, prjType)
    return maintitle


def MapName(pngPath, baseInfo, varName, prjType):
    """
    输出图片路径，数据投影类型，输出文件名
    :param pngPath:
    :param baseInfo:
    :param varName:
    :param prjType:
    :return: 文件名全路径，str类型
    """
    if (prjType in ('NOM', 'NUL')):
        if (varName == baseInfo['dataName']):
            outfileName = "%s%s_QCS_MAPN_ORIG.png" % (pngPath, baseInfo['filename'])
        else:
            if (varName in ('CPD_COT', 'CPD_CER', 'CPD_IWP', 'CPD_LWP', 'CPN_COT', 'CPN_CER', 'CPN_IWP', 'CPN_LWP')):
                varName = varName.split('_')[1]

            outfileName = "%s%s_QCS_MAPN_ORIG_%s.png" % (pngPath, baseInfo['filename'], varName)
    elif (prjType == 'GEO'):
        if (varName == baseInfo['dataName']):
            outfileName = "%s%s_QCS_MAPG_ORIG.png" % (pngPath, baseInfo['filename'])
        else:
            outfileName = "%s%s_QCS_MAPG_ORIG_%s.png" % (pngPath, baseInfo['filename'], varName)
    return outfileName


def GetType(dataName):
    """
    将数据分为两大类型：连续型(True)和离散型（Flase）
    :param dataName:
    :return:
    """
    if (dataName in ('SST', 'LSA', 'ACI', 'OLR', 'LST', 'LPW', 'TBB', 'CTH', 'TFP', 'SSI', 'CTT', 'DLR',
                     'CPD', 'CTP', 'RSR', 'DSD', 'LDA', 'OCA', 'LSE', 'ULR', 'CPN', 'QPE')):
        return True
    elif (dataName in ("SNC", 'CIX', 'CLM', 'CLP', 'CLT', 'FOG', 'FHS')):
        return False


# 自定义colormap
def colormap():
    return matplotlib.colors.LinearSegmentedColormap.from_list('cmap', ["#0d1691", "#140a9f", "#1805a7", "#1b02ae",
                                                                        "#1f00b4", "#1f01ba", "#1c0ac0", "#1918c5",
                                                                        "#1525cb", "#1333d0", "#1040d6", "#0d4ddb",
                                                                        "#0a5ae1", "#0666e6", "#0573ec", "#057eee",
                                                                        "#0d96ee", "#12a2ed", "#17afed", "#1ab9ed",
                                                                        "#1fc5ed", "#23d1ec", "#28ddec", "#2ce9ec",
                                                                        "#35eee7", "#2ce9ec", "#5feec9", "#6ceebe",
                                                                        "#7aeeb5", "#87eeac", "#94eea1", "#a1ee96",
                                                                        "#aeee8d", "#b7ed86", "#bde982", "#c3e57f",
                                                                        "#c9e17b", "#cfdd78", "#d4d974", "#dad571",
                                                                        "#e0d16c", "#e6cd69", "#ecc965", "#e9c261",
                                                                        "#dfb85b", "#d5af54", "#caa44e", "#c09a47",
                                                                        "#b58f42", "#ab853c", "#a17b35", "#96712e",
                                                                        "#8c6628"], 256)


def DrawPara_con(varName):
    if (varName in ("SST", "SST Monthly Mean Product")):
        cmap = plt.cm.jet
        levels = np.arange(-2, 36, 2)
    elif (varName in ("NIR Black sky albedo", "NIR White sky albedo", "SW Black sky albedo",
                      "SW White sky albedo", "VIS Black sky albedo", "VIS White sky albedo")):
        cmap = plt.cm.jet
        levels = np.arange(0, 1.1, 0.1)
    elif (varName == 'OLR'):
        cmap = plt.cm.jet
        levels = np.arange(40, 460, 10)
    elif (varName == 'LST'):
        cmap = plt.cm.jet
        levels = np.arange(250, 360, 10)
    elif (varName in ("LPW_HIGH", "LPW_LOW", "LPW_MID", "TPW")):
        cmap = plt.cm.jet
        levels = np.arange(0, 10.5, 0.5)
    elif (varName in ('NOMChannel07', 'NOMChannel08', 'NOMChannel09', 'NOMChannel10',
                      'NOMChannel11', 'NOMChannel12', 'NOMChannel13', 'NOMChannel14')):
        cmap = plt.cm.jet
        levels = np.arange(250, 360, 10)
    elif (varName == 'CTH'):
        cmap = plt.cm.jet
        levels = np.arange(0, 28., 1)
    elif (varName == "TFTP_Z_depth"):
        cmap = plt.cm.jet
        levels = np.arange(0, 1.1, 0.1)
    elif (varName == 'SSI'):
        cmap = plt.cm.jet
        levels = np.arange(0, 1510, 10)
    elif (varName == 'CTT'):
        cmap = plt.cm.jet
        levels = np.arange(160, 330, 10)
    elif (varName == 'ULR'):
        cmap = plt.cm.jet
        levels = np.arange(50, 510, 10)
    elif (varName == 'DLR'):
        cmap = plt.cm.jet
        levels = np.arange(50, 760, 10)
    elif (varName == 'CPN_CER'):
        cmap = plt.cm.jet
        levels = np.arange(2, 51, 1)
    elif (varName == 'CPD_CER'):
        cmap = plt.cm.jet
        levels = np.arange(0, 51, 1)
    elif (varName == 'CPN_COT'):
        cmap = plt.cm.jet
        levels = np.arange(1, 5.1, 0.1)
    elif (varName == 'CPD_COT'):
        cmap = plt.cm.jet
        levels = np.arange(0, 5.1, 0.1)
    elif (varName == 'CPN_LWP'):
        cmap = plt.cm.jet
        levels = np.arange(25, 105, 5)
    elif (varName == 'CPD_LWP'):
        cmap = plt.cm.jet
        levels = np.arange(0, 510, 10)
    elif (varName == 'CPN_IWP'):
        cmap = plt.cm.jet
        levels = np.arange(25, 180, 5)
    elif (varName == 'CPD_IWP'):
        cmap = plt.cm.jet
        levels = np.arange(0, 510, 10)
    elif (varName == 'CTP'):
        cmap = plt.cm.jet_r
        levels = np.arange(0, 1150, 50)
    elif (varName == 'RSR'):
        cmap = plt.cm.jet
        levels = np.arange(0, 1050, 50)
    elif (varName == 'DSD'):
        cmap = colormap()
        levels = np.arange(0, 41, 1)
    elif (varName == 'AOD'):
        cmap = plt.cm.jet
        levels = np.arange(0, 1.1, 0.1)
    elif (varName in ("Precipitation_3h", "Precipitation_1h", "Precipitation_0.25h")):
        cmap = plt.cm.jet
        levels = np.arange(0, 21, 1)
    elif (varName == "Precipitation_6h"):
        cmap = plt.cm.jet
        levels = np.arange(0, 51, 1)
    elif (varName == "Precipitation_24h"):
        cmap = plt.cm.jet
        levels = np.arange(0, 101, 1)
    elif (varName in ("LSE_8.5um", "LSE_10.8um", "LSE_12.0um")):
        cmap = plt.cm.jet
        levels = np.arange(0.8, 1.005, 0.005)
    return cmap, levels


def DrawPara_dis(varName):
    if (varName == "RDC_Rank"):
        colors = ["#ff0000", "#04ac30", "#81b805", "#b7fe9a", "#f7f80f", "#ff7a00", "#448764", "#c8c8c8"]
        names = ['CI', 'Extinct', 'Preserve', 'Deep Develop', 'RDC4', 'RDC5', 'RDC6', 'Invalid']
        # -1:CI;1:Exctinct;2:preserve;3:deepdevelop;4:rdc4;5:rdc5;6:rdc6
        values = [-1.9, -0.9, 1.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1]
    elif (varName == 'SNC'):
        colors = ["#000000", "#faff01", "#d8bfd6", "#a65026", "#000188", "#ffffff", "#008b8b", "#00ffff", "#a60100"]
        names = ['Bad Data', 'Uncertain', 'Night', 'Clear Land', 'Water', 'Cloud', 'Ice', 'Snow', 'Saturated']
        # 0:NoData;1:Unknown;11:Night;25:Land;39:Water;50:Cloud;100:IceCover;200:SnowCover;254:Saturated
        values = [-0.9, 0.1, 1.1, 11.1, 25.1, 39.1, 50.1, 100.1, 200.1, 254.1]
    elif (varName == 'CLM'):
        colors = ["#f7f7f7", "#a9a9a9", "#00ff00", "#4169e1", "#008000", "#191970"]
        names = ['Cloudy', 'Prob Cloudy', 'Prob LandClear', 'Prob WaterClear', 'LandClear', 'WaterClear']
        # 0:Cloud;1:Probably cloud;4:Probably land clear;5:Probably ocean clear;6:Land Clear;7:Ocean Clear
        values = [-1, 0, 1, 4, 5, 6, 7]
    elif (varName == 'CLP'):
        colors = ["#d0fefe", "#0000f4", "#20a5ff", "#21ffaa", "#ff0000", "#ffff00"]
        names = ['Clear', 'Water', 'Super Cooled', 'Mixed', 'Ice', 'Uncertain']
        # 0:Clear;1:Water;2:Super Cooled;3:Mixed;4:Ice;5:Uncertain;
        values = [-1.1, 0.1, 1.1, 2.1, 3.1, 4.1, 5.1]
    elif (varName == 'CLT'):
        colors = ["#d0fefe", "#0000f4", "#20a5ff", "#21ffaa", "#ff0000", "#b414ff", "#69ff00", "#ffff00"]
        names = ['Clear', 'Water', 'Super Cooled', 'Mixed', 'Ice', 'Cirrus', 'Overlap', 'Uncertain']
        # 0:Clear;2:Water;3:Super Cooled;4:Mixed;5:Ice;6:Cirrus;7:Overlap;9:Uncertain;
        values = [-1.1, 0.1, 2.1, 3.1, 4.1, 5.1, 6.1, 7.1, 9.1]
    elif (varName == 'FOG'):
        colors = ["#c8c8c8", "#f9f87d", "#01ddfc", "#0329f8"]
        names = ['Fill Value', 'FOG', 'Ice Cloud', 'Clear Sky']
        # 0:Fill Value;100:FOG;65519:Ice Cloud;65520:Clear Sky;
        values = [-1.1, 0.1, 100.1, 65519.1, 65520.1]
    elif (varName == 'FHS'):
        colors = ["#fe0000", "#646464", "#505050", "#a0a0a0", "#00ff01", "#c83c3b", "#0aa09f", "#969634", "#01008a",
                  "#c8c8c8"]
        names = ['Fire point', 'Fillvalue', 'SZA>80', 'Flare angle<30', 'Land', 'BT3.9<200K',
                 'BT10.8<200K', 'Desert', 'Water', 'Cloud']
        # 10:fire point,40:fillvalue,50:satallite zenithangle>80,60:flare angle<30,100:land,126:BT3.9um<200K,
        # 127:BT10.8um<200,150:Desert,153:Water,200:cloud01,205:cloud02,210:cloud03,215:cloud04,220:cloud05,
        values = [9.1, 10.1, 40.1, 50.1, 60.1, 100.1, 126.1, 127.1, 150.1, 153.1, 220.1]
    return colors, names, values


def DrawMap_GLL(fy4_mdata, lat, lon, mainTitle, dataUnit, outfileName):
    fig = plt.figure(figsize=(13, 15), dpi=100)  # 图像的长*高=2748*3000
    axes1 = fig.add_axes([0.05, 0.1, 0.94, 0.9])  # 加两个panel，地图区和图例区域，四参数分别为：左边、下边框离边缘距离百分比，绘图区的宽和高
    axes2 = fig.add_axes([0., 0., 1., 0.1], facecolor='#e8e8e8')
    cax = fig.add_axes([0.1, 0.03, 0.65, 0.02])  # 图例

    lon2d, lat2d = np.meshgrid(lon, lat)

    m = Basemap(lon_0=104.7, resolution='l', ax=axes1, llcrnrlat=-90, llcrnrlon=0., urcrnrlat=90, urcrnrlon=180.)
    x, y = m(lon2d, lat2d)
    m.drawcoastlines()
    m.drawcountries()
    m.drawparallels(range(-90, 90, 30), labels=[1, 0, 0, 0], fontsize=14)
    m.drawmeridians(range(0, 360, 45), labels=[0, 0, 0, 1], fontsize=14)

    cs = m.contourf(x, y, fy4_mdata, cmap=plt.cm.jet, levels=np.arange(-2, 36, 2))  # 定义分级
    cb = plt.colorbar(cs, cax=cax, orientation='horizontal')

    cb.ax.tick_params(labelsize=18)

    axes2.spines['left'].set_visible(False)
    axes2.spines['top'].set_visible(False)

    axes2.text(0.21, 0.8, mainTitle,
               horizontalalignment='left',
               verticalalignment='top',
               family='Times New Roman',
               fontsize=24,
               fontweight='bold')
    axes2.text(0.751, 0.35, dataUnit, fontsize=20)

    # 添加logo
    axicon = fig.add_axes([0.85, 0.01, 0.15, 0.05])
    axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'), origin='upper')
    axicon.axis('off')

    fig.savefig(outfileName)
    plt.close()
    return 0


def DrawMap_NOM(fy4_mdata, lat, lon, mainTitle, dataName, varName, dataUnit, outfileName):
    if (dataName == 'QPE'):
        fig = plt.figure(figsize=(25, 21), dpi=50)
        axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
        axes1 = fig.add_axes([0.03, 0.11, 0.95, 0.88])
        cax = fig.add_axes([0.15, 0.03, 0.6, 0.02])

        m = Basemap(projection='cyl', resolution='l', ax=axes1, llcrnrlat=0, llcrnrlon=70, urcrnrlon=140, urcrnrlat=55)
        x, y = m(lon, lat)
        # m.drawmapboundary()
        m.drawcoastlines()
        m.drawcountries()

        m.drawparallels(np.arange(-90., 90., 10.), linewidth=1, labels=[1, 0, 0, 0], fontsize=16)
        m.drawmeridians(np.arange(-180., 180., 10.), linewidth=1, labels=[0, 0, 0, 1], fontsize=16)
        cmap, levels = DrawPara_con(varName)
        cs = m.contourf(x, y, fy4_mdata, cmap=cmap, levels=levels)
        cb = plt.colorbar(cs, cax=cax, orientation='horizontal')
        cb.ax.tick_params(labelsize=18)
        # 标题
        cax.set_title(mainTitle,
                      family='Times New Roman',
                      fontsize=32,
                      fontweight='bold')
        axes2.text(0.751, 0.4, dataUnit, fontsize=24)  # 添加单位
        axes2.spines['left'].set_visible(False)
        axes2.spines['top'].set_visible(False)
        # axes1.spines['left'].set_visible(False)
        # axes1.spines['top'].set_visible(False)
        # 添加logo
        axicon = fig.add_axes([0.85, 0.01, 0.15, 0.05])
        axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'), origin='upper')
        axicon.axis('off')
        fig.savefig(outfileName)
    elif (dataName == 'ACI'):
        fig = plt.figure(figsize=(27.48, 30), dpi=50)  # 图像的长*高=2748*3000
        axes2 = fig.add_axes([0., 0., 1., 0.084], facecolor='#e8e8e8')
        axes1 = fig.add_axes([0., 0.084, 1., 0.916])
        img = Image.fromarray(fy4_mdata)
        axes1.imshow(img)
        axes1.set_xlabel(mainTitle,
                         family='Times New Roman',
                         fontsize=42,
                         fontweight='bold', labelpad=20)
        axes1.spines['left'].set_visible(False)
        axes1.spines['bottom'].set_visible(False)
        axes1.spines['top'].set_visible(False)
        axes2.spines['left'].set_visible(False)
        axes2.spines['top'].set_visible(False)
        axes1.set_xticks([])
        axes1.set_yticks([])
        # 加气象局logo
        axicon = fig.add_axes([0.85, 0.01, 0.15, 0.05])
        axicon.imshow(plt.imread('/FY4APGSQCS/QCS/pyFY4AL2src/Logo/logo.jpg'), origin='upper')
        axicon.axis('off')
        fig.savefig(outfileName)
    elif (dataName == 'DSD'):
        # 初始化底图
        initmap = CMap(mainTitle, dataName, lat, lon)
        m = initmap.CreateMap()
        cmap, levels = DrawPara_con(varName)
        data = fy4_mdata[1]
        mdata = ma.masked_values(data, 65534.)
        mdata[(mdata > 40) & (mdata < np.max(mdata))] = 40

        initmap.AddContiMap(m, mdata, cmap, levels, dataUnit)
        initmap.AddSpecial(m, fy4_mdata[0], varName)
        initmap.fig.savefig(outfileName)
    else:
        # 初始化底图
        initmap = CMap(mainTitle, dataName, lat, lon)
        m = initmap.CreateMap()
        # 叠加数据,如果是连续型数据用AddContiMap方法，如果是离散型数据用AddDiscMap方法
        if (GetType(dataName)):
            cmap, levels = DrawPara_con(varName)
            initmap.AddContiMap(m, fy4_mdata, cmap, levels, dataUnit)
        else:
            colors, names, values = DrawPara_dis(varName)
            initmap.AddDiscMap(m, fy4_mdata, values, colors, names)

        # 特殊图例标注
        initmap.AddSpecial(m, fy4_mdata, varName)
        initmap.fig.savefig(outfileName)
    return 0


def DrawMap_NUL(fy4_mdata, lat, lon, mainTitle, dataName, outfileName):
    # 初始化底图
    initmap = CMap(mainTitle, dataName, lat, lon)
    initmap.fig.set_size_inches(54.96, 60.)
    initmap.fig.set_dpi(100)
    m = initmap.CreateMap()
    initmap.axes1.set_xlabel(mainTitle,
                             family='Times New Roman',
                             fontsize=82,
                             fontweight='bold', labelpad=20)
    # 叠加风类型数据
    initmap.AddBarbsMap(m, fy4_mdata)

    initmap.fig.savefig(outfileName)
    return 0


def LogXML(logMsg, logType, xmlfilePath, inputFilename, outputFilename=True):
    # 在内存中创建一个空的文档
    doc = minidom.Document()
    # 创建一个根节点Managers对象
    root = doc.createElement('XML')
    # 设置根节点的属性
    root.setAttribute('identify', 'L2')
    # 将根节点添加到文档对象中
    doc.appendChild(root)

    nodeOutputFiles = doc.createElement('OutputFiles')
    nodeOutputFile = doc.createElement('OutputFile')
    nodeOutputFiles.appendChild(nodeOutputFile)
    nodeInputFilename = doc.createElement('InputFilename')
    # 给叶子节点name设置一个文本节点，用于显示文本内容
    nodeInputFilename.appendChild(doc.createTextNode(str(inputFilename)))
    if (outputFilename != True):
        nodeOutputFilename = doc.createElement("OutputFilename")
        nodeOutputFilename.appendChild(doc.createTextNode(str(outputFilename)))
        nodeOutputFile.appendChild(nodeOutputFilename)

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
