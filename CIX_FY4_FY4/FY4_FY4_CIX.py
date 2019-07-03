#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/7/11 14:20
# @Author  : baizhaofeng
# @File    : FY4_FY4_CIX.py

import warnings
from Func import *
import sys

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 读取xml文件
        xmlpath = sys.argv[1]
        nodes = Read_xml(xmlpath)
        QCSlog("XML文件解析完成", nodes['xmlPath'])
        # =========================Part1:时间匹配====================================#
        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])
        fy4_filenames, fy4_times = TimeMatch(baseInfo['dirname'], baseInfo['time'], nodes['timevalue'])
        QCSlog("时间匹配完成：" + str(fy4_filenames), nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(baseInfo, nodes, fy4_filenames, fy4_times)
        # =========================Part2:读取数据======================================#
        # 读取t0时刻数据
        t0data = ReadDatas(nodes['filepath'])
        dataT0, t0Label = LabelT0(t0data.copy())
        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataRes'], t0data.shape[0])

        # 读取t1,t2时刻数据
        t1data = ReadDatas(names['t1filepath'])
        t2data = ReadDatas(names['t2filepath'])
        dataT1, t1Label = LabelT1T2(t1data.copy())
        dataT2, t2Label = LabelT1T2(t2data.copy())
        QCSlog("空间匹配完成", nodes['xmlPath'])

        indicator = CheckTB(dataT0, t0Label, dataT1, t1Label, dataT2, t2Label, names)
        # =========================Part3:计算检验指标===================================#
        metric = CalMetric(indicator, t0Label)
        QCSlog("检验指标计算完成", nodes['xmlPath'])

        # =========================Part4:绘柱状图========================================#
        bar_flag = DrawBar(metric, names)
        QCSlog("结果直方图绘图已完成：" + str(names['histname']), nodes['xmlPath'])

        # =========================Part5:绘图======================================#
        # 绘制FY4A T0时刻空间分布图
        fy4_map_flag = DrawMap1(t0data, lat, lon, names['t0title'], names['t0map'],names['tbb_t0path'])
        QCSlog("FY4A T0时刻空间分布图输出完成：" + str(names['t0map']), nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['t0map'] + " -o " + names['t0map']
        os.system(cmd)
        # 绘制FY4A T1时刻空间分布图
        fy4_map1_flag = DrawMap1(t1data, lat, lon, names['t1title'], names['t1map'],names['tbb_t1path'])
        QCSlog("FY4A T1时刻空间分布图输出完成：" + str(names['t1map']), nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['t1map'] + " -o " + names['t1map']
        os.system(cmd)

        # 绘制FY4A T2时刻空间分布图
        fy4_map2_flag = DrawMap1(t2data, lat, lon, names['t2title'], names['t2map'],names['tbb_t2path'])
        QCSlog("FY4A T2时刻空间分布图输出完成：" + str(names['t2map']), nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['t2map'] + " -o " + names['t2map']
        os.system(cmd)
        # 输出日志xml
        LogXML(fy4_map2_flag, "info", nodes['xmlPath'], nodes['filepath'], baseInfo['time'], metric, names)
    except Exception as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
