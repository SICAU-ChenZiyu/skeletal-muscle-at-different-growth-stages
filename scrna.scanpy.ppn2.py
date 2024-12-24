import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import scanpy.external as sce
import os
import operator
import time
import re
from collections import Counter
from collections import OrderedDict

sc.settings.verbosity = 3 ##输出日志的详细程度
sc.logging.print_header() ##日志标题
sc.settings.set_figure_params(dpi=300, facecolor='white') ##图片参数
##########
dir_name = 'seruat_v3.method.nb10.pc30.hg1500.leiden'
os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'
os.chdir('/') ##修改路径
if not os.path.exists('../' + dir_name):
	os.makedirs('../' + dir_name)	
results_file = '../' + dir_name + '/scrna.h5ad' ##设置结果文件保存路径
nonormal_results_file = '../' + dir_name + '/scrna.nonormal.h5ad'
samples_name = os.listdir('.')
if 'cache' in samples_name:
    samples_name.remove('cache')
marker_file = 'marker.list.gid_3.0.txt'

scRNAlist = []
mito_genes = ["ND1", "ND2", "COX1","COX2","ATP8","ATP6","COX3","ND3", "ND4L","ND4", "ND5", "ND6", "CYTB"]
for sample in samples_name:
	sample_file = sample + '/05.Doublet' + '/' + sample + '_doublet_feature_bc_matrix'
	locals()['sc_' + sample] = sc.read_10x_mtx(sample_file, cache=True)
	locals()['sc_' + sample].obs["sample"] = sample
	locals()['sc_' + sample].obs["tissue"] = 'LDM'
	locals()['sc_' + sample].obs["time"] = sample.split('_')[0]
	locals()['sc_' + sample].obs.index = [sample + '_' + x for x in locals()['sc_' + sample].obs.index]
	locals()['sc_' + sample].obs_names_make_unique()
	scRNAlist.append(locals()['sc_' + sample])
scrna = scRNAlist[0].concatenate(scRNAlist[1:len(scRNAlist)], index_unique=None)
mt = []
for i in scrna.var.index:
    if i in mito_genes:
        mt.append(True)
    else:
        mt.append(False)
scrna.var["mt"] = mt
scrna.var["rp"] = scrna.var_names.str.match(r'^RP[SL]')
sc.pp.calculate_qc_metrics(scrna, qc_vars=["mt","rp"], inplace=True)

sc.pp.filter_cells(scrna, min_genes=200)
sc.pp.filter_genes(scrna, min_cells=3)
### # mito_genes = scrna.var_names.str.startswith('MT-')
sc.pl.violin(scrna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, stripplot=False, show=False)
plt.savefig('../' + dir_name + "/QC_violin.pdf", dpi=300, bbox_inches='tight')

scrna = scrna[scrna.obs.n_genes_by_counts > 200, :]
scrna = scrna[scrna.obs.n_genes_by_counts < 4000, :]
sc.pp.filter_cells(scrna, min_counts=200)
sc.pp.filter_cells(scrna, max_counts=10000)
scrna = scrna[scrna.obs.pct_counts_mt < 5, :]
scrna.write(nonormal_results_file)
sc.pl.violin(scrna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True, stripplot=False, show=False)
plt.savefig('../' + dir_name + "/QC_violin.filtered.pdf", dpi=300, bbox_inches='tight')
sc.pp.normalize_total(scrna, target_sum=1e4)
sc.pp.log1p(scrna)
scrna.raw = scrna
sc.pp.highly_variable_genes(scrna, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.highly_variable_genes(scrna, flavor="seurat_v3", n_top_genes=3000)
sc.pl.highly_variable_genes(scrna, show=False)
plt.savefig('../' + dir_name + "/highly_variable_genes.pdf", dpi=300, bbox_inches='tight')
scrna = scrna[:, scrna.var.highly_variable]
sc.pp.regress_out(scrna, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(scrna, max_value=10)

###PCA
sc.tl.pca(scrna, svd_solver='arpack')# svd_solver 指定奇异值分解 SVD 的方法
sc.pl.pca_variance_ratio(scrna, log=True, n_pcs=50, show=False)
plt.savefig('../' + dir_name + "/PCA.pdf", dpi=300, bbox_inches='tight')
sce.pp.harmony_integrate(scrna, 'batch')
scrna.write(results_file)
###UMAP
scrna = sc.read_h5ad(results_file)
sc.pp.neighbors(scrna, n_neighbors=10, n_pcs=30, use_rep='X_pca')
scrna.write(results_file)

scrna = sc.read_h5ad(results_file)
my_colormap = LinearSegmentedColormap.from_list("", ["#BAB8BC", "#0000FF"])
#my_colormap = LinearSegmentedColormap.from_list("", ["#ACC7E7", "#CF383D"])
for louvain_num in ['0.5']:
	if not os.path.exists('../' + dir_name + '/louvain' + louvain_num):
		os.makedirs('../' + dir_name + '/louvain' + louvain_num)
	sc.tl.leiden(scrna, resolution=float(louvain_num))
	cell_num = pd.DataFrame([Counter(scrna.obs["leiden"].tolist())])
	cell_num.T.to_csv('../' + dir_name + "/louvain" + louvain_num + "/Cell_Num.csv",header=False)

	sc.tl.umap(scrna)
	
	sc.pl.umap(scrna, color=['leiden'], legend_loc='on data', legend_fontsize=10, legend_fontoutline=1, show=False)
	plt.savefig('../' + dir_name + "/louvain" + louvain_num + "/Umap.pdf", dpi=300, bbox_inches='tight')
	sc.pl.umap(scrna, color=['pct_counts_mt'], color_map=my_colormap, vmax=5, show=False)
	plt.savefig('../' + dir_name + "/louvain" + louvain_num + "/MTpct_Umap.pdf", dpi=300, bbox_inches='tight')
	sc.pl.umap(scrna, color=['pct_counts_rp'], color_map=my_colormap, vmax=5, show=False)
	plt.savefig('../' + dir_name + "/louvain" + louvain_num + "/RPpct_Umap.pdf", dpi=300, bbox_inches='tight')
	
	os.makedirs('../' + dir_name + '/louvain' + louvain_num + '/Each_Time')
	samples_list = scrna.obs['time'].tolist()
	samples_list_T = list(set(scrna.obs['time'].tolist()) & set(samples_list))
	samples_list_T.sort(key=samples_list.index)
	for batch in samples_list_T:
	    scrna_each = scrna[scrna.obs['time'] == batch, :]
	    sc.pl.umap(scrna_each, color=['leiden'], title=[batch], show=False)
	    plt.savefig('../' + dir_name + '/louvain' + louvain_num + '/Each_Time/' + batch + "_Each_Umap.png", dpi=300, bbox_inches='tight')
	
	os.makedirs('../' + dir_name + '/louvain' + louvain_num + '/Each_Sample')
	samples_list = scrna.obs['sample'].tolist()
	samples_list_T = list(set(scrna.obs['sample'].tolist()) & set(samples_list))
	samples_list_T.sort(key=samples_list.index)
	for batch in samples_list_T:
	    scrna_each = scrna[scrna.obs['sample'] == batch, :]
	    sc.pl.umap(scrna_each, color=['leiden'], title=[batch], show=False)
	    plt.savefig('../' + dir_name + '/louvain' + louvain_num + '/Each_Sample/' + batch + "_Each_Umap.png", dpi=300, bbox_inches='tight')

	sc.tl.dendrogram(scrna, groupby='leiden', n_pcs=30, use_rep='X_pca')

	marker_file = "/Lustre02/ChenZiyu/data/scrna-seq/Sus_scrofa/CQ.muscle/02.result/marker_file.txt"
	marker_data = pd.read_table(marker_file, header=None)
	marker_data.columns = ['cell','gene']
	marker_genes = marker_data['gene'].tolist()
	marker_genes_T = list(set(scrna.raw.var_names.tolist()) & set(marker_genes))
	marker_genes_T.sort(key=marker_genes.index)
	marker_data_T = marker_data[marker_data['gene'].isin(marker_genes_T)]
	marker_data_dic = marker_data_T.groupby('cell')['gene'].apply(list).to_dict()
	sc.pl.violin(scrna, marker_genes_T[0:5], groupby='cell_class',rotation=90, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Marker_Violin1.pdf", dpi=300, bbox_inches='tight')
	sc.pl.violin(scrna, marker_genes_T[5:10], groupby='cell_class',rotation=90, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Marker_Violin2.pdf", dpi=300, bbox_inches='tight')
	sc.pl.violin(scrna, marker_genes_T[10:15], groupby='cell_class',rotation=90, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Marker_Violin3.pdf", dpi=300, bbox_inches='tight')
	sc.pl.violin(scrna, marker_genes_T[15:20], groupby='cell_class',rotation=90, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Marker_Violin4.pdf", dpi=300, bbox_inches='tight')
	sc.pl.violin(scrna, marker_genes_T[20:25], groupby='cell_class',rotation=90, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Marker_Violin5.pdf", dpi=300, bbox_inches='tight')

	marker_data_dic_sort = OrderedDict()
	cell_list = ["Myo_I","Myo_IIB","Myo_IIA","MUSC","Myoblast","Adipocyte","FAPs","SMC","Blood_Endo","Lymphatic_Endo","BT_Cell","Macrophages"]
	for key_name in cell_list:
		marker_data_dic_sort[key_name] = marker_data_dic.get(key_name)
	sc.pl.dotplot(scrna, marker_data_dic_sort, 'leiden', dendrogram=False, show=False,categories_order=["3","0","1","2","14","5","8","6","12","4","11","9","13","7","10"])
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Marker_Dot.pdf", dpi=300, bbox_inches='tight')
	if not os.path.exists('../' + dir_name + '/louvain' + louvain_num + '/Each_Marker'):
		os.makedirs('../' + dir_name + '/louvain' + louvain_num + '/Each_Marker')
	for cell_name in marker_data_dic.keys():
		each_marker_genes = marker_data_dic[cell_name]
		sc.pl.umap(scrna, color=each_marker_genes, color_map=my_colormap, vmax=5, show=False)
		plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Each_Marker/" + cell_name + "_Marker_Umap.pdf", dpi=300, bbox_inches='tight')
	scrna.write(results_file)

	###找差异基因
	#scrna = sc.read_h5ad(results_file)
	scrna.uns['log1p']["base"] = None
	sc.tl.rank_genes_groups(scrna, 'leiden', method='wilcoxon')
	celltype = scrna.obs['leiden'].unique().tolist()
	df = sc.get.rank_genes_groups_df(scrna,group=celltype)
	df.to_csv('../' + dir_name + '/louvain' + louvain_num + "/Dif_Gene_cluster.csv", index=False, header=True, sep='\t')
	sc.pl.rank_genes_groups(scrna, n_genes=30, sharey=False, fontsize=5, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Dif_Gene_cluster.pdf", dpi=300, bbox_inches='tight')
	scrna.write(results_file)

	scrna.obs.to_csv('../' + dir_name + "/louvain" + louvain_num + "/cell_class_obs_matrix", index=True, header=True, sep='\t')

	###注释
	#scrna = sc.read_h5ad(results_file)
	cell_class = ['Myo_IIB1', 'Myo_IIB2', 'Myo_IIB3','Myo_I1','FAPs1','Myo_IIA1','Myoblast1','BT_Cell1','MUSC1','Blood_Endo1','Macrophages1','SMC1','Adipocyte1','Lymphatic_Endo1','Myo_IIB4']
	scrna.obs['cell_class'] = scrna.obs['leiden']
	scrna.rename_categories('cell_class', cell_class)
	scrna.obs['cell_class'] = scrna.obs['cell_class'].str[:-1]
	#scrna.uns['cell_class_colors'] = ['#FF851B','#DA3D3E','#E78DCC','#FFA09F','#237AB5','#C4CA80','#9B6C63','#C57BF8','#4AAD81','#1CBFD0','#B1C6D8']
	sc.pl.umap(scrna, color=['cell_class'], legend_loc='on data', legend_fontsize=7, legend_fontoutline=1, show=False)
	plt.savefig('../' + dir_name + "/louvain" + louvain_num + "/Cell_Class_Umap.pdf", dpi=300, bbox_inches='tight')
	
	scrna.write(results_file)
	
	scrna.obs.to_csv('../' + dir_name + "/louvain" + louvain_num + "/cell_class_obs_matrix", index=True, header=True, sep='\t')

	###找差异基因
	scrna.uns['log1p']["base"] = None
	sc.tl.rank_genes_groups(scrna, 'cell_class', method='wilcoxon')
	celltype = scrna.obs['cell_class'].unique().tolist()
	df = sc.get.rank_genes_groups_df(scrna,group=celltype)
	df.to_csv('../' + dir_name + '/louvain' + louvain_num + "/Dif_Gene_class.csv", index=False, header=True, sep=',')
	sc.pl.rank_genes_groups(scrna, n_genes=30, sharey=False, fontsize=5, show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Dif_Gene_class.pdf", dpi=300, bbox_inches='tight')
	scrna.write(results_file)

	sc.tl.dendrogram(scrna, groupby='cell_class', n_pcs=30, use_rep='X_pca_harmony')
	sc.pl.dendrogram(scrna, 'cell_class',orientation='left')
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Dendorgram_cell_class.pdf", dpi=300, bbox_inches='tight')

	celltype=scrna.obs['cell_class'].unique().tolist() #把所有细胞簇种类拿出来
	deg=sc.get.rank_genes_groups_df(scrna,group=celltype) #把所有细胞簇对应的deg拿出来
	top=deg.groupby('group')
	top5=[] #同样以top5举例
	for i in ["Myo_I","Myo_IIA", "Myo_IIB","MUSC","Myoblast","BT_Cell","Macrophages","Blood_Endo","Lymphatic_Endo","Adipocyte","FAPs","SMC"]: #分群提取top5
##'Myo_IIB1', 'Myo_IIB2', 'Myo_IIB3','Myo_I1','FAPs','Myo_IIA1','Myoblast1','BT_Cell1','MUSC1','Blood_Endo.1','Macrophages1','SMC1','Adipocyte1','Lymphatic_Endo.1'
		tmp=top.get_group(str(i))
		tmp=tmp.sort_values('scores',ascending=False) #按scores排序
		tmp=tmp[tmp['logfoldchanges']>2]
		top5.append(tmp['names'].head(10).tolist())
       #array list 转为 list
	top5=np.array(top5)
	top5=np.reshape(top5,10*len(celltype),'C').tolist()
	print(top5)
	sc.pl.heatmap(scrna, top5, groupby='cell_class', swap_axes=False,vmax=3, dendrogram=True, use_raw=True,show_gene_labels=False,cmap=my_colormap,figsize=[15,18],show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Dif_Gene_heatmap.pdf", dpi=300, bbox_inches='tight')
	sc.pl.heatmap(scrna, top5, groupby='cell_class', swap_axes=False,vmax=3, categories_order=["Adipocyte","FAPs","Myo_I", "Myo_IIA", "Myo_IIB","Myoblast","MUSC", "SMC","Blood_Endo","Lymphatic_Endo","BT_Cell","Macrophages"], use_raw=True,show_gene_labels=True,cmap=my_colormap,figsize=[15,18],show=False)
	sc.pl.rank_genes_groups_heatmap(scrna, n_genes=20, use_raw=True, swap_axes=False,show_gene_labels=True,cmap='bwr',show=False)
	plt.savefig('../' + dir_name + '/louvain' + louvain_num + "/Dif_Gene_heatmap.pdf", dpi=300, bbox_inches='tight')
	
	#scrna = sc.read_h5ad(results_file)

	df = pd.DataFrame(data=scrna.obsm['X_umap'],columns=['umap1','umap2'])
	scrna.obs.to_csv('../' + dir_name + "/louvain" + louvain_num + "/obs_matrix", index=True, header=True, sep='\t')
	df.to_csv('../' + dir_name + "/louvain" + louvain_num + "/umap_matrix", index=False, header=True, sep='\t')
	