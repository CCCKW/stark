from .analysis import MultiMcoolReader
from .plot import DrawCellInternerCount, DrawChromatainInternerCount, DrawDistanceCellInternerCount, DrawDistanceChromatainInternerCount
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import anndata as ad

def pp(opt):
    opt['output'] = opt['output'] + '/software_analysis'

    if not os.path.exists(opt['output']):
        os.makedirs(opt['output'])

    reader  = MultiMcoolReader(directory=opt['dir'], resolution=500000)
    print('extracting mcool files')
    reader.extract_mcool_files()
    print('reading mcool files')
    reader.read_mcool_files()
    print('getting cooler')
    coolers = reader.get_coolers()  
    print('calculating cell counts')
    countlist, ChromCountList = reader.count_values()
    print('calculating distance cell counts')
    distanceCountList, distanceChromCountList, chromLens = reader.count_distance_interaction()
    print('calculating chrom interaction')
    chromInteraction = reader.calculate_Chorm_interactions()
    print('calculating cell similarities')
    cellSimilarities, CellsEmbedding = reader.calculateAjcMatEmbadding()
    print('Drawing fig1-6')
    DrawCellInternerCount(CellCounts = countlist, ChromLenList = chromLens, save=opt['output']) 
    DrawChromatainInternerCount(ChromIntraction = chromInteraction, ChromLenList = chromLens, save=opt['output'])
    DrawDistanceCellInternerCount(CellCounts = distanceCountList, ChromLenList = chromLens, save=opt['output'])
    DrawDistanceChromatainInternerCount(ChromatainCounts=distanceChromCountList, save=opt['output'])

    print('Clustering')
    data = np.array(reader.CellFeature)
    adata = ad.AnnData(data)
    adata.obs['SRR'] = reader.mcool_names
    adata.obs['experiment'] = reader.GSE
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, svd_solver='auto')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=7)
    sc.tl.umap(adata)

    sc.tl.leiden(adata, resolution=opt['resolution'])
    for cr in ['experiment','leiden']:
        dot = adata.obsm['X_umap']
        if cr == 'leiden':
            label = list(np.unique(adata.obs[cr]).astype(int))
        else:
            label = list(np.unique(adata.obs[cr]))
        label = sorted(label)
        ax = plt.subplot(111)
        for cls in label:
            dot = adata.obsm['X_umap'][adata.obs[cr]==str(cls)]
            plt.scatter(dot[:,0],dot[:,1])
        plt.legend(label,loc='best')
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('UMAP1')
        plt.ylabel('UMAP2')
        plt.savefig(opt['output'] + '/UMAP_{}_{}.png'.format(cr,opt['resolution']),dpi=300)
        plt.close()