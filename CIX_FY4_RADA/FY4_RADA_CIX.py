#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/6/21 15:01
# @Author  : Bai
# @File    : FY4_RADA_CIX.py


import warnings
from Func import *
import sys

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 读取xml文件
        # xmlpath = r'D:\QCSDATA\CIX\FY4A-_AGRI--_N_DISK_1047E_L2-_CIX-_MULT_NOM_20180606030000_20180606031459_4000M_V0001.xml'
        xmlpath = sys.argv[1]
        nodes = Read_xml(xmlpath)
        QCSlog("XML文件解析完成", nodes['xmlPath'])
        # =========================Part1:时间匹配====================================#
        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])
        rada_filenames, rada_times = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['timevalue'])
        QCSlog("时间匹配完成：" + str(rada_filenames), nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(baseInfo, nodes, rada_times)
        # ===============Part2:匹配到多个雷达数据，进行合并==============================#
        rada, rlat, rlon = MaxRada(rada_filenames,nodes['exterpath'])
        # =========================Part3:读取数据======================================#
        # 读取数据
        tdata = ReadDatas(nodes['filepath'])
        QCSlog("数据读取完成", nodes['xmlPath'])

        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataRes'], tdata.shape[0])

        # # =========================Part4:空间匹配===========================================#
        rada_nom = BasemapInterp(lat, lon, rlat, rlon, rada)
        QCSlog("RADA数据插值完成", nodes['xmlPath'])

        # 雷达数据剔除小斑点
        rada_nom_pro = postprocess(rada_nom.copy())
        rada_nom_pro[rada_nom == -999] = -999
        QCSlog("雷达插值数据后处理完成", nodes['xmlPath'])

        # =========================Part5:计算检验指标===================================#
        metric, comp = CalMetric(tdata, rada_nom_pro.copy())
        QCSlog("检验指标计算完成", nodes['xmlPath'])

        # =========================Part7:绘柱状图========================================#
        bar_flag = DrawBar(metric, names)
        QCSlog("直方图绘图已完成：" + str(names['histname']), nodes['xmlPath'])

        # =========================Part7:绘空间分布图=======================================#
        # 绘制插值后的rada空间分布图
        rada_NOMmap_flag = DrawNOMMap(rada_nom, names['radatitle'], names['radamap'], lat, lon)
        QCSlog("插值后RADA分布图输出完成：" + str(names['radamap']), nodes['xmlPath'])

        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['radamap'] + " -o " + names['radamap']
        os.system(cmd)

        # 绘制对比分布图
        biasmap = DrawCompMap(comp, names['biastitle'], names['biasmap'], lat.copy(), lon.copy())
        QCSlog("对比分布图输出完成：" + str(names['biasmap']), nodes['xmlPath'])

        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['biasmap'] + " -o " + names['biasmap']
        os.system(cmd)

        # # 绘制原雷达拼图
        # radamap_filename = "%s%s.png" % (nodes['outputpath'], rada_filenames[0])
        # radatile = "全国雷达拼图\n%s-%s UTC" % (rada_times[0].__format__('%Y-%m-%d %H:%M:%S'),
        #                                   rada_times[len(rada_times) - 1].__format__('%H:%M:%S'))  # 标题
        #
        # rada_map_flag = DrawRADAMap(rada, rada_lon, rada_lat, radatile, radamap_filename)
        # print '雷达拼图产品分布图输出完成：', radamap_filename

        # 绘制FY4A原空间分布图
        fy4map = DrawMap(tdata, lat, lon, names['fy4title'], names['fy4map'], names['tbbpath'])
        QCSlog("FY4A空间分布图输出完成：" + str(names['fy4map']), nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['fy4map'] + " -o " + names['fy4map']
        os.system(cmd)
        # 输出日志xml
        LogXML(fy4map, "info", nodes['xmlPath'], nodes['filepath'], baseInfo['time'], metric, names)
    except IndexError as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
