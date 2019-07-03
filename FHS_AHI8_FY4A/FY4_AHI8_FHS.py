#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/10/16 9:17
# @Author  : baizhaofeng
# @File    : FY4_AHI8_FHS.py


from Common import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = r"D:\QCSDATA\FHS\FY4\FY4A-_AGRI--_N_DISK_1047E_L2-_FHS-_MULT_NOM_20180920040000_20180920041459_2000M_V0001.xml"
        nodes = Read_xml(xmlpath)
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        ahi8files, ahi8times = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['timevalue'])
        QCSlog("时间匹配完成:" + str(ahi8files), nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], ahi8times, nodes['outputpath'])
        # 3.读取数据
        fy4_mdata, Reliability = Read_Data(nodes['filepath'])

        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        FY4_LAT, FY4_LON = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataRes'])
        QCSlog("FY4A 火点数据读取完成", nodes['xmlPath'])

        ahi8_mdata = Read_exterData(ahi8files[0])
        AHI8_LAT, AHI8_LON = Read_cLUT()
        QCSlog("AHI8 火点数据读取完成", nodes['xmlPath'])

        ahi8_mnom = Interp(FY4_LAT, FY4_LON, AHI8_LAT, AHI8_LON, ahi8_mdata, nodes['spacevalue'])
        QCSlog("空间插值完成", nodes['xmlPath'])

        MatchPixel, statisRes = MatchFY4_AHI8(fy4_mdata, ahi8_mnom, Reliability)
        QCSlog("FY4 AHI8 火点像元匹配完成", nodes['xmlPath'])

        # 输出统计结果直方图
        metric = DrawHist(statisRes, names['histname'], names['histtitle'])
        QCSlog("统计结果直方图输出完成" + str(names['histname']), nodes['xmlPath'])

        bias_flag = DrawBIASMap(MatchPixel, FY4_LAT, FY4_LON, names['biastitle'], names['BIASMap'])
        QCSlog("匹配结果空间分布图输出完成" + str(names['BIASMap']), nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['BIASMap'] + " -o " + names['BIASMap']
        os.system(cmd)

        # 输出日志xml
        LogXML(bias_flag, "info", nodes['xmlPath'], nodes['filepath'], baseInfo['time'], metric, names)

    except IndexError as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
