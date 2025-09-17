import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import numpy as np
from sc3dg.clustering.re_DeepNanoHiC.Wrapper import *

import json



def workflow(config):


    tmp_config = json.load(open(config))
    os.makedirs(tmp_config['temp_dir'], exist_ok=True)



    if os.path.exists(tmp_config['temp_dir'] + "/embedding.png"):
        print("Embedding image already exists. Exiting to avoid overwriting.")
        sys.exit("Embedding image already exists. Exiting to avoid overwriting.")
    else:
        logfile = open(tmp_config['temp_dir'] + "/task_log.txt", "w", buffering=1)
        logfile.write("Starting Higashi model training...\n")
        for key, value in tmp_config.items():
            logfile.write(f"{key}: {value}\n")
            print(f"{key}: {value}")

    # generage data
    if not os.path.exists(tmp_config['temp_dir'] + "/data.txt") and tmp_config["input_format"] == "higashi_v1":
        logfile.write("Generating data...\n")
        logfile.write("Processing data...\n")
        logfile.write("pair to cooler...\n")
        cload_script = '/cluster2/home/Kangwen/hic/0_temp/method/DeepNanoHi-C-main/cload.py'
        filelist = pd.read_csv(tmp_config['data_dir'] + '/filelist.txt', header=None)
        tmp = filelist.iloc[:, 0].values
        pair_path = os.path.dirname(tmp[0]) 
        print(pair_path)
        coolpath = tmp_config['temp_dir'] + '/cool'
        os.makedirs(coolpath, exist_ok=True)
        genome_ref_path = tmp_config['genome_reference_path']
        cmd = f'python {cload_script} {pair_path} {coolpath} {genome_ref_path} {10000}'

        os.system(cmd)
        logfile.write("pair to cooler done.\n")

        logfile.write("cool to mcool...\n")
        zoom_script = '/cluster2/home/Kangwen/hic/0_temp/method/DeepNanoHi-C-main/re_zoom.py'
        cmd = f'python {zoom_script} {coolpath} {tmp_config["temp_dir"]}/mcool'
        os.makedirs(tmp_config['temp_dir'] + '/mcool', exist_ok=True)
        os.system(cmd)
        
        mcool_path = tmp_config['temp_dir'] + '/mcool'
        datalist = []
        cell = []
        ct = []
        logfile.write("mcool to data...\n")
        count = 0
        for mcool in os.listdir(mcool_path):
            if mcool.endswith('mcool'):
                ct.append(mcool.split('.')[0].split('_')[-1])
                cell.append(mcool.split('.')[0])
                clr = cooler.Cooler(os.path.join(mcool_path, mcool) + '::resolutions/' + str(tmp_config['resolution']))
                pixels = clr.pixels()[:]
                bins = clr.bins()[:]
                pixels['chr1'] = list(bins.iloc[list(pixels['bin1_id']), 0])
                pixels['chr2'] = list(bins.iloc[list(pixels['bin2_id']), 0])
                pixels['start1'] = list(bins.iloc[list(pixels['bin1_id']), 1])
                pixels['start2'] = list(bins.iloc[list(pixels['bin2_id']), 1])
                pixels['cell_id'] = count
                pixels = pixels.loc[:,['cell_id', 'chr1', 'start1', 'chr2', 'start2', 'count']]
                pixels.columns = ['cell_id', 'chrom1', 'pos1', 'chrom2', 'pos2', 'count']
                datalist.append(pixels)
                count += 1
        data = pd.concat(datalist, ignore_index=True)
        data.to_csv(tmp_config['temp_dir'] + "/data.txt", index=False, sep='\t')
        ct = {'cell type': ct}
        pickle.dump(ct, open(tmp_config['temp_dir'] + '/label_info.pickle', 'wb'))
    else:
        filelist = pd.read_csv(tmp_config['data_dir'] + '/filelist.txt', header=None)
        filelist = filelist.iloc[:, 0].values
        ct = [os.path.basename(f).split('.')[0].split('_')[-1] for f in filelist]
        ct = {'cell type': ct}
        pickle.dump(ct, open(tmp_config['temp_dir'] + '/label_info.pickle', 'wb'))
        logfile.write("Data already exists. Skipping data generation.\n")

    model = DeepNanoHiC(config)
    model.process_data()
    logfile.write("Data processing completed.\n")
    model.prep_model()
    logfile.write("Model preparation completed.\n")

    # Training
    model.train_for_embeddings()
    logfile.write("Model training completed.\n")

    cell_embeddings = model.fetch_cell_embeddings()
    logfile.write("Saving cell embeddings...\n")
    # save embedding
    with open(tmp_config['temp_dir'] + '/cell_embeddings.pickle', 'wb') as f:
        pickle.dump(cell_embeddings, f)
    logfile.write("Cell embeddings saved successfully.\n")

    cell_type = model.label_info['cell type']
    from umap import UMAP
    from sklearn.decomposition import PCA
    import seaborn as sns
    import matplotlib.pyplot as plt

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
    plt.savefig(model.config['temp_dir'] + "/embedding.png")
    plt.show()

    logfile.write("UMAP and PCA plots saved successfully.\n")
    logfile.write("Higashi model training completed successfully.\n")
    logfile.close()