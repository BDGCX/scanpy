# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#Tab:自动补齐
#！：调用系统命令，！shell命令
#？：查看对象属性
#？？：查看源代码
#output=!cmd args:执行cmd并赋值
#ipython-pylab启动可集成绘图窗口
'''
StackedVlnPlot
'''
#Scanpy读取loom文件转换为能够操作的anndata对象，以下是python代码
import scanpy as sc
import loompy as lp 
import scanyuan as scy
adata = sc.read_loom("C:/Users/bdgcx/OneDrive/single cell/GSE127235_DKD1_mice/sce.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
marker_genes = ["Cd300lg","Egfl7","Rasgrp3","Ptprb",
                "Cyyr1", "Gpihbp1", "Ramp3","Tspan7",
                "Adgrl4","Fabp4","Sfrp2","Agtr1a","Mgp",
                "Ptn", "Hopx","S1pr3","Cd248","P2rx1",
                "Lhfp","Fhl2","Clic3","Cdkn1c","Enpep",
                "Dpp4","Plce1","Magi2","Rhpn1","Srgap1",
                "Ildr2","Robo2","Mgat5"]
order=['TC','IC','POD','MC','EC']
orders=['EC','MC','POD','IC','TC']
palettes=[ "#F8766D","#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
ax = sc.pl.stacked_violin(adata, marker_genes, 
                          groupby='new_clusters', 
                          rotation=90,
                          figsize=[12,4],
                          categories_order=order,#设置分组顺序
                          #swap_axes=True,#转换X轴和Y轴标签：
                          row_palette=palettes
                          )

#旋转图片角度
from PIL import Image
import matplotlib.pyplot as plt
img = Image.open('./figures/stacked_violin.png')
img.transpose(Image.ROTATE_270)
#保存图片
img.save('stackedviolin.png')
plt.savefig("gene_Vio_celltype.pdf")

#scanyuan包
import scanpy as sc
import scanyuan as scy
ax = scy.stacked_violin_t(adata, marker_genes, 
                          figsize=[8,4], 
                          groupby='new_clusters',
                          order=order,
                          row_palette=palettes)    

