#!/usr/bin/env python3
# coding: utf-8
"""
@author: Yiran Wu  wuyiran55@outlook.com
@last modified by: Yiran Wu
@file:ts_nts.py
@time:12/20/22 9:06 AM
"""

import matplotlib.pyplot as plt
import copy

'''
colormaps['stereo_30'] = ["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#FFFF33",
                                  "#F781BF", "#999999", "#E5D8BD", "#B3CDE3", "#CCEBC5", "#FED9A6", "#FBB4AE",
                                  "#8DD3C7", "#BEBADA", "#80B1D3", "#B3DE69", "#FCCDE5", "#BC80BD", "#FFED6F",
                                  "#8DA0CB", "#E78AC3", "#E5C494", "#CCCCCC", "#FB9A99", "#E31A1C", "#CAB2D6",
                                  "#6A3D9A", "#B15928"]
'''
def plot_ts_nts(ts_file,nts_file, output):
    ts_color="#E41A1C"
    nts_color='#377EB8'

    # reading TS/NTS RPKM
    ts = np.loadtxt(ts_file)
    nts = np.loadtxt(nts_file)

    plt.plot(ts[:,:1],ts[:,1:2],color=ts_color,linewidth=1,label='TS' )
    plt.plot(nts[:,:1],nts[:,1:2],color=nts_color,linewidth=1,label="NTS")

    plt.legend(frameon=False,loc="upper right",fontsize='large') #设置图例无边框，将图例放在左上角
    plt.rcParams['figure.figsize']=(6.0,4.0) #图形大小
    plt.rcParams['savefig.dpi'] = 200 #图片像素
    plt.rcParams['figure.dpi'] = 200 #分辨率


    #plt.legend(loc='upper right', frameon=False)
    plt.ylabel('RPKM', fontsize=12)
    plt.xlabel('Bin', fontsize=12)
    plt.title('Effect of Transcription on AFB1')

    ori_xticks = np.arange(0,151,25)
    replaced_xticks = copy.deepcopy(ori_xticks.tolist())
    replaced_xticks[1] = 'TSS'
    replaced_xticks[-3] = 'TES'

    plt.xticks(ori_xticks,replaced_xticks)

    #设置图框线粗细
    bwith = 1.5 #边框宽度设置为2
    TK = plt.gca()#获取边框
    TK.spines['bottom'].set_linewidth(bwith)
    TK.spines['left'].set_linewidth(bwith)
    TK.spines['top'].set_linewidth(bwith)
    TK.spines['right'].set_linewidth(bwith)
    plt.savefig(output, dpi=600)
