#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/12/27 14:06
# @Author  : baizhaofeng
# @File    : FY4_IQUAM_SST_HIST.py

from Common import *
import time
import sys
import os
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        # xmlpath = r"D:\DATA\QCS\SST\SST\FY4A-_AGRI--_N_DISK_1047E_L3-_SST-_MULT_NOM_20181214000000_20181214235959_4000M_ADV01_HIST.xml"
        xmlpath=sys.argv[1]
        nodes = Read_xml(xmlpath)
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])
        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'],nodes['platformtype'], nodes['outputpath'],baseInfo['dataSour'],nodes['isna'])
        txtfile = MatchTXT(baseInfo['time'], nodes)
        QCSlog("日累积文件数据读取完成："+str(txtfile), nodes['xmlPath'])

        # 读取中匹配的中间文件
        mdat = ReadTXTData(txtfile,baseInfo['time'])
        # 3.读取数据
        fy4_mdata = Read_Data(nodes['filepath'])

        # 4.读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['type'], nodes['filepath'])
        QCSlog("数据读取完成！", nodes['xmlPath'])
        # 7.计算检验指标
        metric, param = Calcate(mdat)
        QCSlog("检验指标计算完成！", nodes['xmlPath'])
        # 8.绘制散点图
        scatter_flag = DrawScatter(mdat, param, metric['QualID3_NUM'], names)
        QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 9.绘制直方图
        hist_flag = DrawHist(mdat['bias'], metric, names)
        QCSlog("直方图绘图完成:"+str(names['histname']), nodes['xmlPath'])

        map_flag = DrawIQUAMMap(mdat, names)
        QCSlog("空间分布图绘图完成:" + str(names['iquamname']), nodes['xmlPath'])

        flag = DrawMap_NOM(fy4_mdata, lat, lon, names)
        QCSlog("FY4 日空间分布图绘图完成:"+str(names['fy4name']), nodes['xmlPath'])

        # 输出日志xml
        LogXML(flag, "info", nodes,baseInfo['dataSour'], baseInfo['time'], metric, names)

    except Exception, ex:
        LogXML(1, ex, nodes,baseInfo['dataSour'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
