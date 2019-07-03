#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/1/9 14:51
# @Author  : baizhaofeng
# @File    : FY4A_NPP_FHS.py

from Common import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        # xmlpath = r"D:\DATA\QCS\FHS\FY4\FY4A-_AGRI--_N_DISK_1047E_L2-_FHS-_MULT_NOM_20180920020000_20180920021459_2000M_V0001_NPP.xml"
        xmlpath=sys.argv[1]
        nodes = Read_xml(xmlpath)
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        modfiles, modtimes = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['timevalue'])
        QCSlog("时间匹配完成,共匹配到%s块npp数据" % len(modfiles), nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], modtimes, nodes['outputpath'])

        # 3.读取FY4数据
        fy4_mdata, Reliability = Read_Data(nodes['filepath'])
        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        FY4_LAT, FY4_LON = Read_LUT()
        QCSlog("FY4A 火点数据读取完成", nodes['xmlPath'])

        MOD_LAT, MOD_LON = ReadExterData(modfiles)
        QCSlog("NPP 火点数据读取完成", nodes['xmlPath'])

        mod_mnom = Interp(FY4_LAT, FY4_LON, MOD_LAT, MOD_LON, nodes['spacevalue'])
        QCSlog("空间插值完成", nodes['xmlPath'])

        MatchPixel, statisRes = MatchFY4_AHI8(fy4_mdata, mod_mnom, Reliability)
        QCSlog("FY4 NPP 火点像元匹配完成", nodes['xmlPath'])

        # 输出统计结果直方图
        metric = DrawHist(statisRes, names['histname'], names['histtitle'])
        QCSlog("统计结果直方图输出完成" + str(names['histname']), nodes['xmlPath'])

        bias_flag = DrawBIASMap(MatchPixel, FY4_LAT, FY4_LON, names['biastitle'], names['BIASMap'])
        QCSlog("匹配结果空间分布图输出完成" + str(names['BIASMap']), nodes['xmlPath'])

        # 输出XML文件
        LogXML(bias_flag, "info", nodes['xmlPath'], nodes['filepath'], baseInfo['time'], metric, names)

    except IndexError as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
    except Exception as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
