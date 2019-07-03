#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2018/9/5 13:31
# @Author  : baizhaofeng
# @File    : FY4_IQUAM_SST_NA.py


from Common_NA import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        xmlpath = sys.argv[1]
        # xmlpath = r"D:\DATA\QCS\SST\SST\FY4A-_AGRI--_N_DISK_1047E_L2-_SST-_MULT_NOM_20181120020000_20181120021459_4000M_V0001_IQUAM_NA.xml"
        nodes = Read_xml(xmlpath)
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], nodes['platformtype'], nodes['outputpath'],
                     nodes['txtPath'], nodes['prjType'])

        iQUAM_filePath = CoarseTimeMatch(nodes['exterpath'], baseInfo['time'])
        QCSlog("时间粗匹配完成:" + str(iQUAM_filePath), nodes['xmlPath'])

        # 3.读取数据
        fy4_mdata, dqf = Read_Data(nodes['filepath'])

        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataRes'])
        # 读取ExterData数据
        exterInfo = Read_exterData(iQUAM_filePath)
        QCSlog("数据读取完成", nodes['xmlPath'])

        # 4.时间精匹配
        exterMatched = FineTimeMatch(baseInfo['time'], exterInfo['times'], nodes['timevalue'], exterInfo)
        QCSlog("时间精匹配完成", nodes['xmlPath'])

        if (exterMatched['lat'].shape[0] == 0):
            raise IndexError("Can't find iquam %s data in the special time threshold" % nodes['platformtype'])

        # 5.空间插值：根据海温的经纬度，在FY4圆盘内找出里离检验源检测点最近的点的经纬度，个数应该和exterMatched中变量长度一致
        fy4near, ndqf = Interp(exterMatched['lat'], exterMatched['lon'], lat, lon, fy4_mdata, dqf, nodes['spacevalue'])
        QCSlog("空间插值完成", nodes['xmlPath'])

        # 6.空间匹配：根据检验源的数据获取途径不同，分别和FY4A的不同质量数据进行比较，FY4质量等级包括：[0:good,1:mid,2:bad],
        mdat = SpaceMatch(fy4near, ndqf, exterMatched, nodes['platformtype'])
        platform = "ship,drifting buoy,tropical moored buoy,coastal moored buoy,argo float,high resolution drifter,imos,crw buoy"
        mdat_all = SpaceMatch(fy4near, ndqf, exterMatched, platform)
        # 剔除异常值-----------------------------------------------------
        mdat['exSST'] = mdat['exSST'][np.abs(mdat['bias']) < float(nodes['na'])]
        mdat['fySST'] = mdat['fySST'][np.abs(mdat['bias']) < float(nodes['na'])]
        mdat['fyDQF'] = mdat['fyDQF'][np.abs(mdat['bias']) < float(nodes['na'])]
        mdat['plon'] = mdat['plon'][np.abs(mdat['bias']) < float(nodes['na'])]
        mdat['plat'] = mdat['plat'][np.abs(mdat['bias']) < float(nodes['na'])]
        mdat['ptype'] = mdat['ptype'][np.abs(mdat['bias']) < float(nodes['na'])]
        mdat['bias'] = mdat['bias'][np.abs(mdat['bias']) < float(nodes['na'])]
        leng = mdat['bias'].shape[0]
        mdat['exSST'] = mdat['exSST'].reshape(leng, 1)
        mdat['fySST'] = mdat['fySST'].reshape(leng, 1)
        mdat['fyDQF'] = mdat['fyDQF'].reshape(leng, 1)
        mdat['plon'] = mdat['plon'].reshape(leng, 1)
        mdat['plat'] = mdat['plat'].reshape(leng, 1)
        mdat['ptype'] = mdat['ptype'].reshape(leng, 1)
        mdat['bias'] = mdat['bias'].reshape(leng, 1)

        mdat_all['exSST'] = mdat_all['exSST'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        mdat_all['fySST'] = mdat_all['fySST'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        mdat_all['fyDQF'] = mdat_all['fyDQF'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        mdat_all['plon'] = mdat_all['plon'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        mdat_all['plat'] = mdat_all['plat'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        mdat_all['ptype'] = mdat_all['ptype'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        mdat_all['bias'] = mdat_all['bias'][np.abs(mdat_all['bias']) < float(nodes['na'])]
        leng = mdat_all['bias'].shape[0]
        mdat_all['exSST'] = mdat_all['exSST'].reshape(leng, 1)
        mdat_all['fySST'] = mdat_all['fySST'].reshape(leng, 1)
        mdat_all['fyDQF'] = mdat_all['fyDQF'].reshape(leng, 1)
        mdat_all['plon'] = mdat_all['plon'].reshape(leng, 1)
        mdat_all['plat'] = mdat_all['plat'].reshape(leng, 1)
        mdat_all['ptype'] = mdat_all['ptype'].reshape(leng, 1)
        mdat_all['bias'] = mdat_all['bias'].reshape(leng, 1)
        # --------------------------------------------------------------
        # 将带标记的匹配结果输出
        WriteToTXT(mdat, baseInfo['time'], names['txtnamelabel'])
        # 将全样本的匹配结果输出
        WriteToTXT(mdat_all, baseInfo['time'], names['txtname'])

        # 7.剔除异常值,计算检验指标
        metric, param = CalcateH(mdat)
        QCSlog("检验指标计算完成", nodes['xmlPath'])

        # 8.绘制散点图
        scatter_flag = DrawScatterH(mdat, param, metric['QualID3_NUM'], names)
        QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 9.绘制直方图
        hist_flag = DrawHist(mdat['bias'], metric, names)
        QCSlog("直方图绘图完成" + str(names['histname']), nodes['xmlPath'])

        shortName = Plattypeprocess(nodes['platformtype'])
        if (len(shortName) == 8):
            # 读取中匹配的中间文件
            mdat_ = ReadTXTData(names['txtname'])
        else:
            mdat_ = ReadTXTData(names['txtnamelabel'])
        QCSlog("TXT数据读取完成", nodes['xmlPath'])

        map_flag = DrawIQUAMMap(mdat_, names)
        QCSlog("空间分布图绘图完成:" + str(names['iquamname']), nodes['xmlPath'])

        # 9.压缩图片
        cmd = "pngquant" + " --force " + names['iquamname'] + " -o " + names['iquamname']
        os.system(cmd)

        # 输出日志xml
        LogXML(map_flag, "info", nodes, baseInfo['time'], metric, names)

    except IndexError as ex:
        LogXML(1, ex, nodes)
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes)
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])
