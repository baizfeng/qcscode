#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/5/7 16:50
# @Author  : baizhaofeng
# @Email   : baizhaofeng@piesat.cn
# @File    : AHI8_IQUAM_SST_MONTH.py

from Common import *
import time
import sys
import os
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = 'D:\DATA\QCS\SST\AHI8\AHI8-_AGRI--_N_DISK_1407E_L3-_SST-_MULT_NOM_20181210020000_20181210021000_4000M_V0001_IQUAM.xml'
        nodes = Read_xml(xmlpath)
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(nodes['time'], nodes['outputpath'], nodes['isna'])
        txtfile = MatchTXT(nodes['time'], nodes)

        # 读取中匹配的中间文件
        mdat = ReadTXTData(txtfile)
        # 7.计算检验指标
        metric, param = Calcate(mdat)
        QCSlog("检验指标计算完成！", nodes['xmlPath'])

        # 8.绘制散点图
        scatter_flag = DrawScatter(mdat, param, metric['QualID3_NUM'], names)
        QCSlog("散点图绘图完成：" + str(names['scatname']), nodes['xmlPath'])

        # 9.绘制直方图
        hist_flag = DrawHist(mdat['bias'], metric, names)
        QCSlog("直方图绘图完成:" + str(names['histname']), nodes['xmlPath'])

        iquammap_flag = DrawIQUAMMap(mdat, names)
        QCSlog("偏差空间分布图绘图完成:" + str(names['iquamname']), nodes['xmlPath'])

        # 输出日志xml
        LogXML(iquammap_flag, "info", nodes, metric, names)

    except Exception, ex:
        LogXML(1, ex, nodes)
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
