#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/8/3 14:55
# @Author  : baizhaofeng
# @File    : Draw_FY4A_L2_MAIN.py

from Common import *
import time
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = r"D:\DATA\QCS\L2\CLM\FY4A-_AGRI--_N_DISK_1047E_L2-_CLM-_MULT_NOM_20190524020000_20190524021459_4000M_V0001.xml"
        nodes = Read_xml(xmlpath)
        # 2.数据的基本信息
        baseInfo, nodes['varName'] = GetBaseInfo(nodes['inputFileName'], nodes['varName'])

        # 3.读取数据
        fy4_data = Read_Data(nodes['inputFileName'], nodes['varName'])
        # 4.读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['inputFileName'], baseInfo['dataRes'], fy4_data.shape[0])

        # 5.掩膜无效值
        fy4_mdata = MaskInValid(fy4_data, nodes['varName'])

        # 6.产品图标题
        mainTitle = MapTitle(baseInfo, nodes['varName'], nodes['prjType'])

        # 7.输出产品图名称
        outfileName = MapName(nodes['pngPath'], baseInfo, nodes['varName'], nodes['prjType'])

        # 8.绘图
        if (nodes['prjType'] == 'GEO'):
            flag = DrawMap_GLL(fy4_mdata, lat, lon, mainTitle, baseInfo['dataUnit'], outfileName)
        elif (nodes['prjType'] == 'NOM'):
            flag = DrawMap_NOM(fy4_mdata, lat, lon, mainTitle, baseInfo['dataName'], nodes['varName'],
                               baseInfo['dataUnit'], outfileName)
        elif (nodes['prjType'] == 'NUL'):
            flag = DrawMap_NUL(fy4_mdata, lat, lon, mainTitle, baseInfo['dataName'], outfileName)
        # 9.压缩图片
        cmd = "pngquant" + " --force " + outfileName + " -o " + outfileName
        os.system(cmd)
        # 输出日志xml
        LogXML(flag, "info", nodes['xmlPath'], os.path.basename(nodes['inputFileName']), os.path.basename(outfileName))

    except Exception, ex:
        print ex
        LogXML(1, ex, nodes['xmlPath'], os.path.basename(nodes['inputFileName']))
    else:
        print "程序处理成功！"
