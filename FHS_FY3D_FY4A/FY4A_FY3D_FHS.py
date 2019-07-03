#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/9/27 10:01
# @Author  : baizhaofeng
# @File    : FY4A_FY3B_FHS.py

from Common import *
import sys
import warnings

warnings.filterwarnings("ignore")
np.set_printoptions(suppress=True)
if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = r"D:\QCSDATA\FHS\FY4\FY4A-_AGRI--_N_REGC_1047E_L2-_FHS-_MULT_NOM_20180901071500_20180901071916_2000M_V0001.xml"
        nodes = Read_xml(xmlpath)

        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        fy3files, fy3times = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['timevalue'])
        QCSlog("时间匹配完成:" + str(fy3files), nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(baseInfo, nodes, fy3times)

        # 读取NC数据并输出txt
        txtpath, beginline, endline = ConvertData(nodes['filepath'], names['txtpath'])
        # txtpath, beginline, endline = [names['txtpath'], 0, 5495]
        QCSlog("FY4A 火点数据转换完成:" + str(txtpath), nodes['xmlPath'])
        fy4dat = Read_Data(txtpath)
        QCSlog("FY4A 火点数据读取完成", nodes['xmlPath'])

        fy3dat = Read_exterData(fy3files)
        QCSlog("FY3D 火点数据读取完成", nodes['xmlPath'])
        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataRes'], beginline, endline)

        # 判断FY4火点一定空间阈值范围内是否有FY3火点，fy4dat_check为在FY3范围内的FY4所有火点, nearOrNot为判断结果，1为邻近区域有FY3火点，0为无火点
        fy4dat_check, nearOrNot = PixelNearFY3(fy4dat, fy3dat, nodes['spacevalue'])

        # 得到经纬度范围
        ranges = PixelRange(fy4dat, lat, lon)

        # 插值得到最邻近经纬度
        row_arr, col_arr, rowValid = InterpFY4(ranges, nodes['spacevalue'])  # FY4插值为等经纬度后四个角点和中心点的行列号

        QCSlog("FY4A 火点数据插值完成", nodes['xmlPath'])
        fy4PixelGLL = PixelGLL(row_arr, col_arr)  # FY4插值为等经纬度后的行列号，四个角点组成的多边形内部的点
        QCSlog("FY4A 内部点查找完成", nodes['xmlPath'])
        fy3PixelGLL = InterpFY3(fy3dat, nodes['spacevalue'])
        QCSlog("FY3 火点数据插值完成", nodes['xmlPath'])
        # 插值后的FY4,FY3火点匹配
        MatchPixel = MatchFY(fy4PixelGLL, fy3PixelGLL)
        QCSlog("FY4 FY3 火点像元匹配完成", nodes['xmlPath'])
        # 进一步对MatchPixel进行处理，如果插值后FY3周围7*7范围内有FY4的点则显示该点，如果没有，则不显示该点
        MatchPixel_reject = ReJectFY3(MatchPixel)

        MatchDot, MatchArea, trueNum = Match_Dot_Area(fy4dat, row_arr, col_arr, nearOrNot, rowValid)
        QCSlog("FY4 FY3 火点火区匹配完成", nodes['xmlPath'])
        # 输出统计结果直方图
        metric = DrawHist(fy4dat_check, trueNum, names['histname'], names['histtitle'])
        QCSlog("统计结果直方图输出完成" + str(names['histname']), nodes['xmlPath'])

        pixelCondition, rect = PixelConditionMap(MatchPixel_reject, baseInfo['time'], fy3times[0],
                                                 names['pixelCondMap'])
        QCSlog("火点像元位置匹配实况图绘图成功", nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['pixelCondMap'] + " -o " + names['pixelCondMap']
        os.system(cmd)

        areaSketch = AreaSketchMap(MatchArea, baseInfo['time'], fy3times[0], names['areaSketchMap'], rect)
        QCSlog("火区位置匹配示意图绘图成功", nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['areaSketchMap'] + " -o " + names['areaSketchMap']
        os.system(cmd)

        pixelSketch = PixelSketchMap(MatchDot, baseInfo['time'], fy3times[0], names['pixelSketchMap'], rect)
        QCSlog("火点像元位置匹配示意图绘图成功", nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['pixelSketchMap'] + " -o " + names['pixelSketchMap']
        os.system(cmd)

        # 输出日志xml
        LogXML(pixelSketch, "info", nodes['xmlPath'], nodes['filepath'], baseInfo['time'], metric, names)

    except IndexError as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes['xmlPath'], nodes['filepath'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
