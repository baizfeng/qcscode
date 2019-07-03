#!/usr/bin/env python
# _*_ coding:UTF-8 _*_
# @Time:2018/12/2014:48
# @Author:baizhaofeng
# @File:FY2GH_IQUAM_SST_HIST.py

from Common import *
import time
import sys
import os
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        # xmlpath = r"D:\DATA\QCS\SST\FY2H\FY2H_SST_MLT_NOM_20181126_AOAD_DHIST.xml"
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
        if (baseInfo['dataSour'] == 'FY2H'):
            mdat = ReadTXTDataH(txtfile,baseInfo['time'])
        elif (baseInfo['dataSour'] == 'FY2G'):
            mdat = ReadTXTDataG(txtfile,baseInfo['time'])
        # 3.读取数据
        fy4_mdata = Read_Data(nodes['filepath'],baseInfo['dataSour'])

        # 4.读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['type'], nodes['filepath'], baseInfo['dataSour'])
        QCSlog("数据读取完成！", nodes['xmlPath'])
        # 7.计算检验指标
        if (baseInfo['dataSour'] == 'FY2H'):
            metric, param = CalcateH(mdat)
            QCSlog("检验指标计算完成！", nodes['xmlPath'])
            # 8.绘制散点图
            scatter_flag = DrawScatterH(mdat, param, metric['QualID3_NUM'], names)
            QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])
        elif(baseInfo['dataSour'] == 'FY2G'):
            metric, param = CalcateG(mdat)
            QCSlog("检验指标计算完成！", nodes['xmlPath'])
            # 8.绘制散点图
            scatter_flag = DrawScatterG(mdat, param, metric['QualID3_NUM'], names)
            QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 9.绘制直方图
        hist_flag = DrawHist(mdat['bias'], metric, names,baseInfo['dataSour'])
        QCSlog("直方图绘图完成:"+str(names['histname']), nodes['xmlPath'])

        map_flag = DrawIQUAMMap(mdat, names)
        QCSlog("空间分布图绘图完成:" + str(names['iquamname']), nodes['xmlPath'])

        flag = DrawMap_NOM(fy4_mdata, lat, lon, names,baseInfo['dataSour'])
        QCSlog("FY2 空间分布图绘图完成:"+str(names['fy4name']), nodes['xmlPath'])

        # 输出日志xml
        LogXML(flag, "info", nodes,baseInfo['dataSour'], baseInfo['time'], metric, names)

    except Exception, ex:
        LogXML(1, ex, nodes,baseInfo['dataSour'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
