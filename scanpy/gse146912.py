# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:49 2022

@author: bdgcx
"""

'''
scanpy
'''
#========================================================
'''
Preprocessing and clustering
'''
#加载模块
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, facecolor='white',color_map = 'viridis_r')

# create the file that will store the analysis results
results_file = 'C:/Users/bdgcx/OneDrive/python/scanpy/sce_GSE146912.h5ad'

#Scanpy读取loom文件(Seurat处理的对象转化为loom文件)转换为能够操作的anndata对象
import loompy as lp 
adata = sc.read_loom("C:/Users/bdgcx/OneDrive/python/scanpy/adata_GSE146912.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32') 
adata
#adata是一个AnnData 对象,核心就是一个cell×gene的二维表
#.X的对象cell相关的信息记录在.obs中，obs是存储注释信息的dataframe,属性gene的信息记录在.var中，其他的信息在.uns中
adata.obs['treatment']
adata.to_df()
type(adata.obs['nCount_RNA'])

#cell QC

#Show those genes that yield the highest fraction of counts in each single cell, across all cells
sc.pl.highest_expr_genes(adata, n_top=20, )

#Basic filtering:
sc.pp.filter_cells(adata, min_genes=200)
#filtered out 1 cells that have less than 200 genes expressed
sc.pp.filter_genes(adata, min_cells=10)

#A violin plot of some of the computed quality measures:

# 'n_genes_by_counts'=='n_genes'=='nFeature_RNA':the number of genes expressed in the count matrix
# 'total_counts'=='n_counts'=='nCount_RNA':the total counts per cell 
# 'percent_mito':the percentage of counts in mitochondrial genes
sc.pl.violin(adata, ['nFeature_RNA','nCount_RNA'],
             jitter=0.4, multi_panel=True)
sc.pl.violin(adata, [ 'percent.mt','percent.ribo'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, x='nCount_RNA', y='percent.mt')
sc.pl.scatter(adata, x='nCount_RNA', y='nFeature_RNA')

#根据统计的n_genes 和线粒体基因含量percent_mito 进一步过滤
adata = adata[adata.obs.nCount_RNA < 2500, :]
adata.obs['percent_mito']=adata.obs['percent.mt']
adata = adata[adata.obs.percent_mito < 5, :]

#特征提取
#在特征提取之前要保证细胞之间是有可比性的，一般用的是归一化的方法，
#得到高变基因之后，为了使同一个基因在不同细胞之间具有可比性采用标准化

#Total-count normalize
sc.pp.normalize_total(adata, target_sum=1e4)
#Logarithmize the data
sc.pp.log1p(adata)

#Identify highly-variable genes
sc.pp.highly_variable_genes(adata, 
                            min_mean=0.1, 
                            max_mean=3, 
                            min_disp=0.5)

sc.pl.highly_variable_genes(adata)

#Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
adata.raw = adata

adata = adata[:, adata.var.highly_variable]

#回归每个细胞的总计数和线粒体基因表达百分比的影响
sc.pp.regress_out(adata, ['nCount_RNA', 'percent_mito', 'percent.ribo'])
#Scale each gene to unit variance. Clip values exceeding standard deviation 10
sc.pp.scale(adata, max_value=10)

#Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')

sc.pl.pca(adata, color='Kdr')
sc.pl.pca_variance_ratio(adata, log=True)
adata.write(results_file)
adata

#Computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
#Clustering the neighborhood graph
sc.tl.leiden(adata)
#Embedding the neighborhood graph
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')

sc.tl.umap(adata)
sc.pl.umap(adata, color=['Kdr', 'Cd74', 'Lyz2'])

sc.tl.tsne(adata)
sc.pl.tsne(adata, color='cell type',legend_loc='on data',)
#plot the scaled and corrected gene expression by explicitly stating that you don’t want to use
sc.pl.umap(adata, color=['Cd74', 'Lyz2'], use_raw=False)

sc.pl.umap(adata, color=['Phase', 'Kdr', 'Cd74'])

#embedding的散点图其它画法
#Visualization of gene expression and other variables
#For the scatter plots, the value to plot is given as the color argument. This can be any gene or any column in .obs

# rc_context is used for the figure size, in this case 4x4
with rc_context({'figure.figsize': (4, 4)}):
    sc.pl.umap(adata, color='Cd79a')

# compute clusters using the leiden method and store the results with the name `clusters`
sc.tl.leiden(adata, key_added='clusters', resolution=0.1)

with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color='clusters', add_outline=True, legend_loc='on data',
               legend_fontsize=12, legend_fontoutline=2,frameon=False,
               title='clustering of cells', palette='Set1')    

with rc_context({'figure.figsize': (3, 3)}):
     sc.pl.umap(adata, color=['Pdgfrb', 'Nphs1','Kdr',
                              'Cd74','Atp1b1', 
                              'n_counts', 'clusters'], 
                s=50, frameon=False, ncols=4, 
                vmax='p99')
stress = ["Fos",#"Scos3",
          "Jun","Junb","Atf3","Egr1",
          "Cebpd","Hspa1b","Hsp90aa1"]
 
with rc_context({'figure.figsize': (3, 3)}):
     sc.pl.umap(adata, color=["Fos",
               "Jun","Junb","Atf3","Egr1",
               "Cebpd","Hspa1b","Hsp90aa1", 'clusters'], 
                s=50, frameon=False, ncols=3, 
                vmax='p99')
with rc_context({'figure.figsize': (3, 3)}):
     sc.pl.umap(adata, color=['treatment', 'clusters'], 
                s=50, frameon=False, ncols=3, 
                vmax='p99',use_raw=False)     
#Finding marker genes
sc.tl.rank_genes_groups(adata, 'clusters', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(results_file)

adata = sc.read(results_file)

#查看差异分析的结果
sc.get.rank_genes_groups_df(adata, group=['0'])
#指定哪两组进行差异分析
sc.tl.rank_genes_groups(adata, 'clusters', groups=['0'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)

#差异基因的可视化
marker_genes = ["Flt1", "Pecam1",#EC
                  "Gata3", "Pdgfrb",#MC
                  "Nphs1", "Pdpn",#POD
                  "Acta2", "Myh11",#SMC
                  "Lyz2", "Ptprc",#IC
                  "Atp1b1","Cdh16",#TEC
                  "Cldn1","Pax8"]#PEC
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
#Get a table with the scores and groups
pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 
                                    'logfoldchanges',
                                    'pvals',
                                    'pvals_adj',
                                    'scores']}).head(5)

sc.pl.violin(adata, marker_genes, groupby='clusters')
sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
#Visualize marker genes using dotplot   
sc.pl.rank_genes_groups_dotplot(adata, n_genes=4)

sc.pl.rank_genes_groups_dotplot(adata, 
                                n_genes=4, 
                                values_to_plot='logfoldchanges', 
                                min_logfoldchange=3, 
                                vmax=7, vmin=-7, 
                                cmap='bwr') 
#Focusing on particular groups
sc.pl.rank_genes_groups_dotplot(adata, 
                                n_genes=30, 
                                values_to_plot='logfoldchanges', 
                                min_logfoldchange=4, 
                                vmax=7, vmin=-7, 
                                cmap='bwr', 
                                groups=['0', '4']) 
# scale and store results in layer
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X  
adata.write(results_file)
#Visualize marker genes using matrixplot
sc.pl.rank_genes_groups_matrixplot(adata, n_genes=10,
                                   #use_raw=True,
                                   vmin=-3, 
                                   vmax=3, 
                                   #layer='scaled',
                                   cmap='bwr')
#Visualize marker genes using stacked violin plots
sc.pl.rank_genes_groups_stacked_violin(adata, 
                                       n_genes=3, 
                                       cmap='viridis_r')  

import loompy as lp 
import scanyuan as scy
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
ax = scy.stacked_violin_t(adata, marker_genes, 
                          figsize=[8,4], 
                          groupby='new_clusters',
                          order=order,
                          row_palette=palettes)    

#Visualize marker genes using heatmap
sc.pl.rank_genes_groups_heatmap(adata, n_genes=3, 
                                #use_raw=False, 
                                swap_axes=True, 
                                vmin=-2, vmax=2, 
                                cmap='bwr', 
                                #layer='scaled', 
                                figsize=(10,7), 
                                show=False)
sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, 
                                #use_raw=False, 
                                swap_axes=True, 
                                show_gene_labels=False,
                                vmin=-3, vmax=3, 
                                cmap='bwr')
#Visualize marker genes using tracksplot
sc.pl.rank_genes_groups_tracksplot(adata, n_genes=3)
adata = sc.read(results_file)

sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)

#定义细胞类型
#方法1
new_cluster_names = [
    'MC', 'IC',
    'EC1', 'EC2',
    'EC3', 'POD1',
    'TC', 'POD2']
adata.rename_categories('leiden', new_cluster_names)

sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False)
sc.pl.dotplot(adata, marker_genes, groupby='leiden')
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90)

#方法2
#Identification of clusters based on known marker genes    
marker_genes_dict = {
    "EC" :["Flt1","Tie1","Pecam1", "Kdr","Emcn","Cdh5"],
    "MC" :["Pdgfrb", "Gata3", "Des", "Itga8"],
    "POD":["Nphs1", "Nphs2", "Pdpn", "Wt1", 
           "Mafb", "Synpo", "Cdkn1c", "Ptpro"],
    "IC" :["Ptprc", "Lyz1", "Csf1r", "Itgam", "Cd3d", 
          "Ms4a1","Lyz2","Cd74","H2-Aa"],
    "PEC":["Cldn1", "Pax8"],
    "SMC" :["Acta2", "Myh11", "Tagln"],
    "TEC" :["Fxyd2","Slc12a1", "Slc12a3", "Slc14a2",
            "Aqp1","Aqp2","Umod","Atp1b1","Cdh16"]
}

sc.pl.dotplot(adata, marker_genes_dict, 'clusters', 
              dendrogram=True)    

# create a dictionary to map cluster to annotation label
cluster2annotation = {
     '0': 'EC',
     '1': 'MC',
     '2': 'EC',
     '3': 'IC',
     '4': 'IC',
     '5': 'IC',
     '6': 'POD',
     '7': 'TC'
}

# add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function
adata.obs['cell type'] = adata.obs['clusters'].map(cluster2annotation).astype('category')    
sc.pl.dotplot(adata, marker_genes_dict, 
              groupby='cell type', 
              dendrogram=True)   
#umap图
sc.pl.umap(adata, color='cell type',
           legend_loc='on data',
           frameon=False, 
           legend_fontsize=10, 
           legend_fontoutline=2) 

adata1=adata[adata.obs['treatment']].isin(['Control'])
#violin plot
with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata, ['Cd74', 'Pdgfrb'], groupby='cell type' )    
    
with rc_context({'figure.figsize': (4.5, 3)}):
    sc.pl.violin(adata, ['n_genes', 'percent_mito'], 
                 groupby='cell type', 
                 stripplot=False, inner='box')
   

# use stripplot=False to remove the internal dots, inner='box' adds a boxplot inside violins    

#stacked-violin plot  
sc.pl.stacked_violin(adata, marker_genes_dict, 
                          groupby='cell type', 
                          #dendrogram=True,
                          swap_axes=False
                          ) 
#matrixplot
sc.pl.matrixplot(adata, marker_genes_dict, 
                 'cell type', 
                 #dendrogram=True, 
                 cmap='Blues', 
                 standard_scale='var', 
                 colorbar_title='column scaled\nexpression') 

# scale and store results in layer
adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X    
sc.pl.matrixplot(adata, marker_genes_dict, 
                 groupby='cell type', 
                 dendrogram=True,
                 colorbar_title='mean z-score', 
                 #layer=('scaled'), 
                 vmin=-2, vmax=2, cmap='RdBu_r')    
#Combining plots in subplots
import matplotlib.pyplot as plt

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20,4), gridspec_kw={'wspace':0.9})

ax1_dict = sc.pl.dotplot(adata, marker_genes_dict, groupby='cell type', ax=ax1, show=False)
ax2_dict = sc.pl.stacked_violin(adata, marker_genes_dict, groupby='cell type', ax=ax2, show=False)
ax3_dict = sc.pl.matrixplot(adata, marker_genes_dict, groupby='cell type', ax=ax3, show=False, cmap='viridis')    
    
#Heatmaps
sc.pl.heatmap(adata, marker_genes_dict, 
                   groupby='cell type', 
                   cmap='viridis', 
                   dendrogram=True)    

sc.pl.heatmap(adata, marker_genes_dict, 
                   groupby='cell type', 
                   layer='scaled', 
                   vmin=-2, vmax=2, 
                   cmap='RdBu_r', 
                   dendrogram=True, 
                   swap_axes=True, 
                   figsize=(11,4))    
    
#Tracksplot
sc.pl.tracksplot(adata, marker_genes_dict, 
                      groupby='cell type', 
                      dendrogram=True)  



#Comparison of marker genes using split violin plots
with rc_context({'figure.figsize': (9, 1.5)}):
    sc.pl.rank_genes_groups_violin(adata, n_genes=20, 
                                   jitter=False)  
#Plot correlation
sc.pl.correlation_matrix(adata, 'cell type', 
                              figsize=(5,3.5))   
adata.write(results_file) 
#保存分析结果以供下一次调用
adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading
#If you want to share this file with people who merely want to use it for visualization, a simple way to reduce the file size is by removing the dense scaled and corrected data matrix. The file still contains the raw data used in the visualizations in adata.raw.
adata.raw.to_adata().write('./adata_withoutX.h5ad')
#If you want to export to “csv”, you have the following options
# Export single fields of the annotation of observations
# adata.obs[['n_counts', 'louvain_groups']].to_csv(
#     './write/pbmc3k_corrected_louvain_groups.csv')

# Export single columns of the multidimensional annotation
# adata.obsm.to_df()[['X_pca1', 'X_pca2']].to_csv(
#     './write/pbmc3k_corrected_X_pca.csv')

# Or export everything except the data using `.write_csvs`.
# Set `skip_data=False` if you also want to export the data.
# adata.write_csvs(results_file[:-5], )


#Load dataset
results_file = 'C:/Users/bdgcx/OneDrive/python/scanpy/sce.h5ad'
adata =sc.read(results_file)
adata
