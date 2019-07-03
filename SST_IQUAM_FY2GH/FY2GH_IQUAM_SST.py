#!/usr/bin/env python
# _*_ coding:UTF-8 _*_
# @Time:2018/12/209:20
# @Author:baizhaofeng
# @File:FY2GH_IQUAM_SST.py


from Common import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = r"D:\DATA\QCS\SST\FY2G\FY2G_SST_003_NOM_20181120_0200_IQUAM.xml"
        nodes = Read_xml(xmlpath)

        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], nodes['platformtype'], nodes['outputpath'],
                     nodes['txtPath'], nodes['prjType'], baseInfo['dataSour'])

        iQUAM_filePath = CoarseTimeMatch(nodes['exterpath'], baseInfo['time'])
        QCSlog("时间粗匹配完成:" + str(iQUAM_filePath), nodes['xmlPath'])

        # 3.读取数据
        if (baseInfo['dataSour'] == 'FY2H'):
            fy4_data, dqf = Read_HData(nodes['filepath'], baseInfo['dataSour'])
        elif (baseInfo['dataSour'] == 'FY2G'):
            fy4_data = Read_GData(nodes['filepath'], baseInfo['dataSour'])
        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        FY4_LAT, FY4_LON = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataSour'])
        QCSlog("FY2 SST数据读取完成", nodes['xmlPath'])
        # 读取ExterData数据
        exterInfo = Read_exterData(iQUAM_filePath)
        QCSlog("IQUAM数据读取完成", nodes['xmlPath'])

        # 4.时间精匹配
        exterMatched = FineTimeMatch(baseInfo['time'], exterInfo['times'], nodes['timevalue'], exterInfo)
        QCSlog("时间精匹配完成", nodes['xmlPath'])

        if (exterMatched['lat'].shape[0] == 0):
            raise IndexError("Can't find iquam %s data in the special time threshold" % nodes['platformtype'])
        platform = "ship,drifting buoy,tropical moored buoy,coastal moored buoy,argo float,high resolution drifter,imos,crw buoy"

        # 5.空间插值：根据海温的经纬度，在FY4圆盘内找出里离检验源检测点最近的点的经纬度，个数应该和exterMatched中变量长度一致
        if (baseInfo['dataSour'] == 'FY2H'):
            fy4near, ndqf = InterpH(exterMatched['lat'], exterMatched['lon'], FY4_LAT, FY4_LON, fy4_data, dqf,
                                    nodes['spacevalue'])
            QCSlog("空间插值完成", nodes['xmlPath'])
            # 6.空间匹配：根据检验源的数据获取途径不同，分别和FY4A的不同质量数据进行比较，FY4质量等级包括：[0:good,1:mid,2:bad],
            mdat = SpaceMatchH(fy4near, ndqf, exterMatched, nodes['platformtype'])
            mdat_all = SpaceMatchH(fy4near, ndqf, exterMatched, platform)
            # 将带标记的匹配结果输出
            WriteToTXTH(mdat, baseInfo['time'], names['txtnamelabel'])
            # 将全样本的匹配结果输出
            WriteToTXTH(mdat_all, baseInfo['time'], names['txtname'])
            # 7.计算检验指标
            metric, param = CalcateH(mdat)
            QCSlog("检验指标计算完成", nodes['xmlPath'])

            # 8.绘制散点图
            scatter_flag = DrawScatterH(mdat, param, metric['QualID3_NUM'], names)
            QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])
        elif (baseInfo['dataSour'] == 'FY2G'):
            fy4near = InterpG(exterMatched['lat'], exterMatched['lon'], FY4_LAT, FY4_LON, fy4_data, nodes['spacevalue'])
            QCSlog("空间插值完成", nodes['xmlPath'])
            # 6.空间匹配：根据检验源的数据获取途径不同，分别和FY4A的不同质量数据进行比较，FY4质量等级包括：[0:good,1:mid,2:bad],
            mdat = SpaceMatchG(fy4near, exterMatched, nodes['platformtype'])
            mdat_all = SpaceMatchG(fy4near, exterMatched, platform)
            # 将带标记的匹配结果输出
            WriteToTXTG(mdat, baseInfo['time'], names['txtnamelabel'])
            # 将全样本的匹配结果输出
            WriteToTXTG(mdat_all, baseInfo['time'], names['txtname'])
            # 7.计算检验指标
            metric, param = CalcateG(mdat)
            QCSlog("检验指标计算完成", nodes['xmlPath'])
            # 8.绘制散点图
            scatter_flag = DrawScatterG(mdat, param, metric['QualID3_NUM'], names, baseInfo['dataSour'])
            QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 9.绘制直方图
        hist_flag = DrawHist(mdat['bias'], metric, names, baseInfo['dataSour'])
        QCSlog("直方图绘图完成" + str(names['histname']), nodes['xmlPath'])

        shortName = Plattypeprocess(nodes['platformtype'])
        if (baseInfo['dataSour'] == 'FY2H'):
            if (len(shortName) == 8):
                # 读取中匹配的中间文件
                mdat_ = ReadTXTDataH(names['txtname'])
            else:
                mdat_ = ReadTXTDataH(names['txtnamelabel'])
        elif (baseInfo['dataSour'] == 'FY2G'):
            if (len(shortName) == 8):
                # 读取中匹配的中间文件
                mdat_ = ReadTXTDataG(names['txtname'])
            else:
                mdat_ = ReadTXTDataG(names['txtnamelabel'])
        QCSlog("TXT数据读取完成", nodes['xmlPath'])

        map_flag = DrawIQUAMMap(mdat_, names)
        QCSlog("空间分布图绘图完成:" + str(names['iquamname']), nodes['xmlPath'])

        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['iquamname'] + " -o " + names['iquamname']
        os.system(cmd)

        # 输出日志xml
        LogXML(map_flag, "info", nodes, baseInfo['dataSour'], baseInfo['time'], metric, names)
    except IndexError as ex:
        LogXML(1, ex, nodes, baseInfo['dataSour'])
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes, baseInfo['dataSour'])
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
