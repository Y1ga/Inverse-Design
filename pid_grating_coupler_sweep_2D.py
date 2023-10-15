"""
    Copyright (c) 2020 Ansys Inc. """
'''
    根据初始（4，1）numpy生成一个耦合器base '''

######## IMPORTS ########
# General purpose imports
import os,sys
print(sys.path);
sys.path.append("W:\\Software\\Lumerical\\v202\\api\\python\\")
import scipy as sp
import numpy as np
import json
from lumjson import LumEncoder, LumDecoder

import lumapi
######## 参数初始化 单位：m ########
######## OPTIMIZABLE GEOMETRY ########
lambda_c = 1.55e-6 
bandwidth_in_nm = 0     #< Only optimize for center frequency of 1550nm
F0 = 0.95
# 初始高度
height = 220e-9
etch_depth = 80e-9
y0 = 0
# 初polygon的坐标x（固定值）
x_begin = -5.1e-6
x_end = 22e-6
# 槽()/尺的数目
n_grates = 25

indexSi = 3.47668
indexSiO2 = 1.44401

data_file = "pid_grating_coupler_initial_params.json"
base_file = "pid_grating_coupler_2D_TE_base.fsp"

# input: （4, 1)numpy： F, R, a, b & (int)n_grates: 周期数
# output: (104, 2): polygon的104个点坐标
def grating_params_pos(params, n_grates):
    y0      = 0

    y3 = y0+height
    # 厚度 = max厚度-刻蚀高度
    y1 = y3-etch_depth
    # 1e-6转化为微米单位
    x_start = params[0]*1e-6  #< First parameter is the starting position
    R  = params[1]*1e6        #< second parameter (unit is 1/um)
    a  = params[2]            #< Third parameter (dim-less)
    b  = params[3]            #< Fourth parameter (dim-less)

    x0 = x_start
    # verts: polygon第一个多边形的4个点坐标
    verts = np.array( [[x_begin,y0],[x_begin,y3],[x0,y3],[x0,y1]] )
    
    ## Iterate over all but the last tooth
    for i in range(n_grates-1):
        # 公式： 填充因子F随波导位置而变化
        F = F0-R*(x0-x_start)
        # Lambda：光栅的空间变化周期数
        Lambda = lambda_c / (a+F*b)
        # 公式：1个槽/尺4个坐标(x,y)的迭代
        x1 = x0 + (1-F)*Lambda    #< Width of the etched region
        x2 = x0 + Lambda          #< Rest of cell
        # 更新verts，1次更新1个槽/尺=4个点(x,y)
        verts = np.concatenate((verts,np.array([[x1,y1],[x1,y3],[x2,y3],[x2,y1]])),axis=0)
        x0 = x2

    ## Last tooth is special
    F = F0-R*(x0-x_start)
    Lambda = lambda_c / (a+F*b)
    x1 = x0 + (1-F)*Lambda        #< Width of the etched region
    verts = np.concatenate((verts,np.array([[x1,y1],[x1,y3],[x_end,y3],[x_end,y0]])),axis=0) 

    return verts


if __name__ == "__main__":
    # 打开写了x_start;F;a;b的json文件并读取
    with open(data_file) as fh:
        # json文件读取，为一个4×1np数组
        initial_params = json.load(fh, cls=LumDecoder)["initial_params"][0]
                
    # Alternate starting point
    # 2.5：初始位置坐标；0.03：切趾因子（曲率半径？）；
    # ?a = neff_thin - index_SiO2*sin(theta);
    # ?b = neff_thick - neff_thin;
    # initial_params = [ -2.5, 0.03, 2.4, 0.5369]
    # 初始化1次光栅耦合器=25个槽/尺
    # vtx(104, 2)：存储25个矩形+1个初始坐标点+1个结束坐标点=100+2+2=104， 2表示二维
    ############## 生成一个polygon base ############
    vtx = grating_params_pos(initial_params, 25)
    # fdtd，启动！
    if os.path.exists(base_file):
        with lumapi.FDTD(filename=base_file) as fdtd:
            # 设置fdtd文件中的基本参数
            fdtd.addpoly()
            fdtd.set("vertices", vtx)
            fdtd.set("x", 0)
            fdtd.set("y", 0)
            fdtd.set("index", indexSi)
            fdtd.setglobalsource("center wavelength", lambda_c)
            fdtd.setglobalsource("wavelength span", 0)
            fdtd.save()
            # runsweep: 优化，启动！
            fdtd.runsweep("sweep source position")
            # 打印优化后的透射率——波长曲线： sweep
            # sweep_pos: (21, 1): 表示波长
            # getsweepdata()：局部；getsweetresult:全部
            sweep_pos = fdtd.getsweepdata("sweep source position", "x")
            sweep_T = fdtd.getsweepdata("sweep source position", "T")
            # 提取最优值
            Tmax = np.amax(sweep_T)
            # 这个提取为什么看着这么麻烦？
            Tpos = sweep_pos[np.where(sweep_T == np.amax(sweep_T))[0][0]][0]
            #
            print("Max transmission:", Tmax*100, "%")
            print("Position", Tpos*1e6, "um")
            
            fdtd.setnamed("source", "x", Tpos)
            # 选取所有名为“polygon”的structure
            # 然后删除
            fdtd.select("polygon")
            fdtd.delete()
            fdtd.save()            
    else:
        print("base file doesn't exist...")
