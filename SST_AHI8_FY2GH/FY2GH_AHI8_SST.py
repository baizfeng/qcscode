#!/usr/bin/env python
# _*_ coding:UTF-8 _*_
# @Time:2018/12/1915:15
# @Author:baizhaofeng
# @File:FY2H_AHI8_SST_NA.py


from Common import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = r"D:\DATA\QCS\SST\FY2H\FY2H_SST_003_NOM_20181210_0200_AHI8.xml"
        nodes = Read_xml(xmlpath)
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        ahi8files, ahi8times = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['timevalue'])
        QCSlog("时间匹配完成:" + str(ahi8files), nodes['xmlPath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], baseInfo['dataSour'], nodes['outputpath'],
                     nodes['prjType'], nodes['cVar'])
        # 3.读取数据
        if (baseInfo['dataSour'] == 'FY2H'):
            fy4_data, dqf = Read_HData(nodes['filepath'], baseInfo['dataSour'])
        elif (baseInfo['dataSour'] == 'FY2G'):
            fy4_data = Read_GData(nodes['filepath'], baseInfo['dataSour'])

        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        FY4_LAT, FY4_LON = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataSour'])
        QCSlog("FY2 SST数据读取完成", nodes['xmlPath'])

        ahi8_data = Read_exterData(ahi8files[0])
        AHI8_LAT, AHI8_LON = Read_cLUT()
        QCSlog("AHI8 SST数据读取完成", nodes['xmlPath'])

        ahi8_nom = Interp(FY4_LAT, FY4_LON, AHI8_LAT, AHI8_LON, ahi8_data, nodes['spacevalue'])
        QCSlog("空间插值完成", nodes['xmlPath'])

        # # 绘图前数据预处理，将无效值掩摸
        mask = (fy4_data < -5) | (fy4_data > 45) | (ahi8_nom < -5) | (ahi8_nom > 45)  # 有效值为[-5,45]
        exter_mdata = ma.masked_where(mask, ahi8_nom)
        fy4_mdata = ma.masked_where(mask, fy4_data)
        if (baseInfo['dataSour'] == 'FY2H'):
            dqf_m = ma.masked_where(mask, dqf)
        bias_ma = fy4_mdata - exter_mdata

        # 7.计算检验指标
        if (baseInfo['dataSour'] == 'FY2H'):
            metric = CalcateDQF(fy4_mdata, exter_mdata, bias_ma, dqf_m)
        elif (baseInfo['dataSour'] == 'FY2G'):
            metric = Calcate(fy4_mdata, exter_mdata, bias_ma)
        QCSlog("检验指标计算完成", nodes['xmlPath'])

        # 绘制散点图
        scatter_flag = DrawScatter(fy4_mdata, exter_mdata, metric, names['scatname'], names['scattitle'], nodes['cVar'],
                                   baseInfo['dataSour'])
        QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 绘制偏差统计直方图
        hist_flag = DrawHist(bias_ma, metric, names['histname'], names['histtitle'], nodes['cVar'],
                             baseInfo['dataSour'])
        QCSlog("直方图绘图完成:" + str(names['histname']), nodes['xmlPath'])

        bias_flag = DrawNOMMap(bias_ma, FY4_LAT, FY4_LON, baseInfo['dataSour'], names['biastitle'], names['biasname'],
                               np.arange(-5, 6, 1))
        QCSlog("偏差空间分布图输出完成：" + str(names['biasname']), nodes['xmlPath'])
        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['biasname'] + " -o " + names['biasname']
        os.system(cmd)
        # 输出日志xml
        LogXML(bias_flag, "info", nodes, baseInfo['dataSour'], baseInfo['time'], metric, names)
    except IndexError as ex:
        LogXML(1, ex, nodes, baseInfo['dataSour'])
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes, baseInfo['dataSour'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
