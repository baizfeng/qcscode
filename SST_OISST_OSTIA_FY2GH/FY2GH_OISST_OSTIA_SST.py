#!/usr/bin/env python
# _*_ coding:UTF-8 _*_
# @Time:2018/12/1916:48
# @Author:baizhaofeng
# @File:FY2GH_OISST_OSTIA_SST.py

from Common import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        # xmlpath = r"D:\DATA\QCS\SST\FY2G\FY2G_SST_003_NOM_20181130_0200_OISST.xml"
        xmlpath=sys.argv[1]
        nodes = Read_xml(xmlpath)
        print nodes

        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], baseInfo['dataSour'], nodes['outputpath'],
                     nodes['prjType'], nodes['cVar'])

        check_filePath = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['cVar'])
        QCSlog("时间匹配完成:" + str(check_filePath), nodes['xmlPath'])

        # 3.读取数据
        if (baseInfo['dataSour'] == 'FY2H'):
            fy4_data, dqf = Read_HData(nodes['filepath'], baseInfo['dataSour'])
        elif (baseInfo['dataSour'] == 'FY2G'):
            fy4_data = Read_GData(nodes['filepath'], baseInfo['dataSour'])
        print '666666666'
        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataSour'])
        # 读取ExterData数据
        exter_data = Read_exterData(check_filePath, nodes['prjType'], nodes['cVar'])
        QCSlog("数据读取完成", nodes['xmlPath'])

        # # 绘图前数据预处理，将无效值掩摸
        mask = (fy4_data < -5) | (fy4_data > 45) | (exter_data < -5) | (exter_data > 45)  # 有效值为[-5,45]
        exter_mdata = ma.masked_where(mask, exter_data)
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
        scatter_flag = DrawScatter(fy4_mdata, exter_mdata, metric, names['scatname'], names['scattitle'], nodes['cVar'],baseInfo['dataSour'])
        QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 绘制偏差统计直方图
        hist_flag = DrawHist(bias_ma, metric, names['histname'], names['histtitle'], nodes['cVar'],baseInfo['dataSour'])
        QCSlog("直方图绘图完成:" + str(names['histname']), nodes['xmlPath'])

        bias_flag = DrawNOMMap(bias_ma, lat, lon, baseInfo['dataSour'], names['biastitle'], names['biasname'],
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
