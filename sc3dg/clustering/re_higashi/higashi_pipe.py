import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import sys
import pickle
import numpy as np

import json
from sc3dg.clustering.re_higashi.higashi2.Higashi_wrapper import *

def workflow(config):
    higashi_model = Higashi(config)

    tmp_config = json.load(open(config))

    if os.path.exists(tmp_config['data_dir'] + "/embedding.png"):
        print("Embedding image already exists. Exiting to avoid overwriting.")
        sys.exit("Embedding image already exists. Exiting to avoid overwriting.")
    else:
        logfile = open(tmp_config['temp_dir'] + "/task_log.txt", "w", buffering=1)
        logfile.write("Starting Higashi model training...\n")
        for key, value in tmp_config.items():
            logfile.write(f"{key}: {value}\n")
            print(f"{key}: {value}")


    logfile.write("Processing data...\n")
    higashi_model.process_data()


    filelist = pd.read_csv(higashi_model.config['data_dir'] + '/filelist.txt',header=None)
    nums = filelist.shape[0]

    # 加载label
    label_info = {'name': np.arange(nums ),
                'lb': np.array([x.split('/')[-1].split('.')[0].split('_')[-1] for x in filelist[0]]),
                }

    pickle.dump(label_info, open(higashi_model.config['data_dir'] + "/label_info.pickle", "wb"))




    logfile.write("Processing data...\n")
    higashi_model.prep_model()

    logfile.write("Training Higashi model...\n")
    higashi_model.train_for_embeddings()


    logfile.write("Fetching cell embeddings...\n")
    cell_embeddings = higashi_model.fetch_cell_embeddings()

    logfile.write("Saving cell embeddings...\n")
    # save embedding
    with open(higashi_model.config['data_dir'] + '/cell_embeddings.pickle', 'wb') as f:
        pickle.dump(cell_embeddings, f)

    logfile.write("Cell embeddings saved successfully.\n")
    from umap import UMAP
    from sklearn.decomposition import PCA
    import seaborn as sns
    import matplotlib.pyplot as plt
    cell_type = higashi_model.label_info['lb']
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
    plt.savefig(higashi_model.config['data_dir'] + "/embedding.png")
    plt.show()

    logfile.write("UMAP and PCA plots saved successfully.\n")
    logfile.write("Higashi model training completed successfully.\n")
    logfile.close()


