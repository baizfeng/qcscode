#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/2/18 15:23
# @Author  : baizhaofeng
# @File    : AHI8_OISST_OSTIA_SST_NA.py


from Common_NA import *
import sys
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":
    try:
        # 1.解析XML文件
        # xmlpath = r"D:\DATA\QCS\SST\AHI8\AHI8-_AGRI--_N_DISK_1407E_L2-_SST-_MULT_NOM_20181210020000_20181210021000_4000M_V0001_OISST_NA.xml"
        xmlpath=sys.argv[1]
        nodes = Read_xml(xmlpath)
        nodes['cVar']=nodes['cVar'].split('_')[0]
        QCSlog("xml文件解析完成！", nodes['xmlPath'])

        # 2.数据的基本信息
        baseInfo = GetBaseInfo(nodes['filepath'])

        # 输出的文件名和标题
        names = Name(baseInfo['filename'], baseInfo['time'], nodes['outputpath'], nodes['prjType'], nodes['cVar'])

        check_filePath = TimeMatch(nodes['exterpath'], baseInfo['time'], nodes['cVar'])
        QCSlog("时间匹配完成:" + str(check_filePath), nodes['xmlPath'])

        # 3.读取数据
        fy4_data, dqf = Read_Data(nodes['filepath'])
        # 读取经纬度数据,如果是NOM，则经纬度查找表为绝对路径
        lat, lon = Read_LUT(nodes['prjType'], nodes['filepath'], baseInfo['dataRes'])
        # 读取ExterData数据
        exter_data = Read_exterData(check_filePath, nodes['prjType'],nodes['cVar'])
        QCSlog("数据读取完成", nodes['xmlPath'])

        # # 绘图前数据预处理，将无效值掩摸
        mask = (fy4_data < -5) | (fy4_data > 45) | (exter_data < -5) | (exter_data > 45)  # 有效值为[-5,45]
        exter_mdata = ma.masked_where(mask, exter_data)
        fy4_mdata = ma.masked_where(mask, fy4_data)
        dqf_m = ma.masked_where(mask, dqf)
        bias_ma = fy4_mdata - exter_mdata

        # 处理异常值-------------------------------------------------------------------------------------------------
        mask_na = np.abs(bias_ma) >= float(nodes['na'])
        bias_ma = ma.masked_where(mask_na, bias_ma)
        exter_mdata = ma.masked_where(mask_na, exter_mdata)
        fy4_mdata = ma.masked_where(mask_na, fy4_mdata)
        dqf_m = ma.masked_where(mask_na, dqf_m)
        # -----------------------------------------------------------------------------------------------------------

        # 7.计算检验指标
        metric = Calcate(fy4_mdata, exter_mdata, bias_ma, dqf_m)
        QCSlog("检验指标计算完成", nodes['xmlPath'])

        # 绘制散点图
        scatter_flag = DrawScatter(fy4_mdata, exter_mdata, metric, names['scatname'], names['scattitle'], nodes['cVar'])
        QCSlog("散点图绘图完成:" + str(names['scatname']), nodes['xmlPath'])

        # 绘制偏差统计直方图
        hist_flag = DrawHist(bias_ma, metric, names['histname'], names['histtitle'], nodes['cVar'])
        QCSlog("直方图绘图完成:" + str(names['histname']), nodes['xmlPath'])

        if (nodes['prjType'] == 'GLL'):
            bias_flag = DrawGLLMap(bias_ma, lat, lon, names['biastitle'], names['biasname'], np.arange(-5, 6, 1))
            QCSlog("偏差空间分布图输出完成：" + str(names['biasname']), nodes['xmlPath'])
            # 9.压缩图片
            cmd = "pngquant" + " --force " + names['biasname'] + " -o " + names['biasname']
            os.system(cmd)
        elif (nodes['prjType'] == 'NOM'):
            bias_flag = DrawNOMMap(bias_ma, lat, lon, names['biastitle'], names['biasname'], np.arange(-5, 6, 1))
            QCSlog("偏差空间分布图输出完成：" + str(names['biasname']), nodes['xmlPath'])
            # 9.压缩图片
            cmd = "pngquant" + " --force " + names['biasname'] + " -o " + names['biasname']
            os.system(cmd)
        # 输出日志xml
        LogXML(bias_flag, "info", nodes, baseInfo['time'], metric, names)

    except IndexError as ex:
        LogXML(1, ex, nodes)
        QCSlog(str(ex), nodes['xmlPath'])
    except Exception as ex:
        LogXML(1, ex, nodes)
        QCSlog(str(ex), nodes['xmlPath'])
    else:
        QCSlog("程序处理成功！", nodes['xmlPath'])