import os
import json
import pandas as pd
import pickle
import sys
from code.ckw_generatematrix_cell import generatematrix_cell
from code.ckw_embed_concatcell_chr import embed_concatcell_chr
from code.ckw_impute_cell import impute_cell
from code.ckw_embed_mergechr import embed_mergechr
from joblib import Parallel, delayed
import h5py
import umap
import numpy as np
import matplotlib.pyplot as plt
from umap import UMAP
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt
def workflow(config):
    config_path = sys.argv[1]
    config = json.load(open(config_path))
    print(config)
    filelist = config['filelist']
    n_jobs = config['n_jobs']

    pair_list = []
    pair_name = []
    filelist = pd.read_csv(filelist, header=None)
    for i in range(len(filelist)):
        pair = filelist.iloc[i, 0]
        pair_list.append(pair)
        name = os.path.basename(pair)
        name = name.split('.')[0]  # Remove file extension
        pair_name.append(name)

    print(f'Total {len(pair_list)} pairs found ')
    if os.path.exists(config['temp_dir'] + "/embedding.png"):
        print("Embedding image already exists. Exiting to avoid overwriting.")
        sys.exit("Embedding image already exists. Exiting to avoid overwriting.")
    else:
        logfile = open(config['temp_dir'] + "/task_log.txt", "w", buffering=1)
        logfile.write("Starting scHiCluster pipeline...\n")
        for key, value in config.items():
            logfile.write(f"{key}: {value}\n")
            print(f"{key}: {value}")    
        
        # delete the dir:cell_matrix,emb_matrix,imputed_matrix,merged_embed_matrix
        for dir_name in ['cell_matrix', 'emb_matrix', 'imputed_matrix', 'merged_embed_matrix']:
            dir_path = os.path.join(config['temp_dir'], dir_name)
            if os.path.exists(dir_path):
                logfile.write(f"Deleting directory: {dir_path}\n")
                for file in os.listdir(dir_path):
                    file_path = os.path.join(dir_path, file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)
                os.rmdir(dir_path)
                logfile.write(f"Deleted directory: {dir_path}\n")
            else:
                logfile.write(f"Directory does not exist: {dir_path}\n")

    ## 1.generatematrix_cell
    ## 2.impute_cell
    ## 3.embed_concatcell_chr
    ## 4. embed_mergechr


    ## Parallel processing generatematrix_cell
    logfile.write("Generating matrix for each cell pair...\n")
    Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(generatematrix_cell)(config_path,  pair) for pair in pair_list
    )
    print('                                             ')
    print('                                             ')
    print('            generatematrix_cell done         ')
    print('                                             ')
    print('                                             ')

    # Parallel processing impute_cell
    logfile.write("Imputing cell matrices...\n")
    chrom_list = config['chrom_list']
    for c in chrom_list:
        Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(impute_cell)(config_path, pair_name[i], c) for i in range(len(pair_name))
        )
    print('                                             ')
    print('                                             ')
    print('            impute_cell done                 ')
    print('                                             ')
    print('                                             ')

    # Parallel processing embed_concatcell_chr
    # first get list
    logfile.write("Embedding concatenated cell matrices for each chromosome...\n")
    imputed_path = config['temp_dir'] + '/imputed_matrix'
    embed_path = config['temp_dir'] + '/embed_matrix'
    os.makedirs(embed_path, exist_ok=True)
    chrom_dict = {}
    chrom_list = config['chrom_list']
    for c in chrom_list:
        chrom_matrix_list = []
        for pair in pair_name:
            tmp = imputed_path + f'/{pair}/{pair}_{c}_pad1_std1_rp0.5_sqrtvc.hdf5'
            chrom_matrix_list.append(tmp)
        chrom_dict[c] = chrom_matrix_list
    # embed_concatcell_chr(config,
    #     chrom_dict[chrom_list[0]],  # Using the first chromosome's list for embedding
    #     outprefix= os.path.join(embed_path, f'embed_{chrom_list[0]}_pad1_std1_rp0.5_sqrtvc'))
    logfile.write("Starting embedding concatenation for each chromosome...\n")
    Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(embed_concatcell_chr)(config, 
                                    chrom_dict[c], 
                                    outprefix=os.path.join(embed_path, f'embed_{c}_pad1_std1_rp0.5_sqrtvc'))
        for c in chrom_list
    )
    print('                                             ')
    print('                                             ')
    print('          embed_concatcell_chr done          ')
    print('                                             ')
    print('                                             ')

    logfile.write("Merging embeddings across chromosomes...\n")
    embed_npy_list = []
    for npy in os.listdir(embed_path):
        if npy.endswith('.npy'):
            embed_npy_list.append(os.path.join(embed_path, npy))
    print(f'Total {len(embed_npy_list)} embed files found in {embed_path}')
    merged_embed_path = config['temp_dir'] + '/merged_embed_matrix'
    os.makedirs(merged_embed_path, exist_ok=True)
    outprefix = os.path.join(merged_embed_path, 'embed_merged_pad1_std1_rp0.5_sqrtvc')
    embed_mergechr(embed_npy_list, outprefix, dim=20, norm_sig=True)


    # Load the merged embedding
    logfile.write("Loading merged embeddings...\n")
    merged_embed = config['temp_dir'] + '/merged_embed_matrix/embed_merged_pad1_std1_rp0.5_sqrtvc.svd20.hdf5'
    with h5py.File(merged_embed , 'r') as f:
        cell_embeddings = f['data'][()]
        
    with open(config['temp_dir'] + '/cell_embeddings.pickle', 'wb') as f:
        pickle.dump(cell_embeddings, f)

    filelist = pd.read_csv(config['filelist'], header=None)
    cell_type = [x.split('/')[-1].split('.')[0].split('_')[-1] for x in filelist[0].tolist()]

    logfile.write("Cell embeddings loaded successfully.\n")
    logfile.write("Plotting embeddings...\n")
    fig = plt.figure(figsize=(14, 5))
    ax = plt.subplot(1, 2, 1)
    vec = PCA(n_components=2).fit_transform(cell_embeddings)
    sns.scatterplot(x=vec[:, 0], y=vec[:, 1], hue=cell_type, ax=ax, s=6, linewidth=0)
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
    ax = plt.subplot(1, 2, 2)
    vec = UMAP(n_components=2).fit_transform(cell_embeddings)
    sns.scatterplot(x=vec[:, 0], y=vec[:, 1], hue=cell_type, ax=ax, s=6, linewidth=0)
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., ncol=1)
    plt.tight_layout()
    plt.savefig(config['temp_dir'] + "/embedding.png")
    plt.show()

    logfile.write("Embedding plot saved successfully.\n")
    logfile.close()
