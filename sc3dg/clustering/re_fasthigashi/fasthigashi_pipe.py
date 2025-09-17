import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from fasthigashi.FastHigashi_Wrapper import *
import sys
import pickle
import numpy as np
from sc3dg.clustering.re_fasthigashi.higashi2.Higashi_wrapper import *

def workflow(config):
    higashi_model = Higashi(config)

    tmp_config = json.load(open(config))
    tmp_config = json.load(open(config))
    if os.path.exists(tmp_config['data_dir'] + "/embedding.png"):
        print("Embedding image already exists. Exiting to avoid overwriting.")
        sys.exit("Embedding image already exists. Exiting to avoid overwriting.")
    else:
        logfile = open(tmp_config['temp_dir'] + "/task_log.txt", "w", buffering=1)
        logfile.write("Starting FastHigashi model training...\n")
        # print tmp_config
        for key, value in tmp_config.items():
            logfile.write(f"{key}: {value}\n")
            print(f"{key}: {value}")
        
        

    logfile.write("Processing data...\n")
    higashi_model.process_data()


    filelist = pd.read_csv(higashi_model.config['data_dir'] + '/filelist.txt',header=None)
    nums = filelist.shape[0]
    label_info = {'name': np.arange(nums ),
                'lb': np.array([x.split('/')[-1].split('.')[0].split('_')[-1] for x in filelist[0]]),
                }

    pickle.dump(label_info, open(higashi_model.config['data_dir'] + "/label_info.pickle", "wb"))






    logfile.write("initializing FastHigashi model...\n")
    fh_model = FastHigashi(config_path=config,
                        path2input_cache=tmp_config['temp_dir'],
                        path2result_dir=tmp_config['temp_dir'],
                        off_diag=100, # 0-100th diag of the contact maps would be used.
                        filter=False, # fit the model on high quality cells, transform the rest
                        do_conv=False,# linear convolution imputation
                        do_rwr=False, # partial random walk with restart imputation
                        do_col=False, # sqrt_vc normalization
                        no_col=False) # force to not do sqrt_vc normalization

    # Pack from sparse mtx to tensors
    # Initialize FastHigashi model
    logfile.write("Packing data into tensors...\n")
    fh_model.prep_dataset()


    logfile.write("Running FastHigashi model...\n")
    filelist = pd.read_csv(fh_model.config['data_dir'] + '/filelist.txt', header=None)
    # 默认分解的中间dim是256，但如果细胞小于256个，则取细胞数//2
    if 256 > filelist.shape[0]:
        rank = filelist.shape[0] // 2
    else:
        rank = 256
    fh_model.run_model(dim1=.6,
                    rank=rank,
                    n_iter_parafac=1,
                    extra="")

    logfile.write("Fetching cell embeddings...\n")
    cell_embeddings = fh_model.fetch_cell_embedding(final_dim=256)# dict
    # save embedding
    logfile.write("Saving cell embeddings...\n")
    with open(fh_model.config['data_dir'] + '/cell_embeddings.pickle', 'wb') as f:
        pickle.dump(cell_embeddings, f)

    cell_embeddings = cell_embeddings['embed_all']
    logfile.write("Cell embeddings saved successfully.\n")



    logfile.write("Plotting embeddings...\n")
    from umap import UMAP
    from sklearn.decomposition import PCA
    import seaborn as sns
    import matplotlib.pyplot as plt
    cell_type = fh_model.label_info['lb']
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
    plt.savefig(fh_model.config['data_dir'] + "/embedding.png")
    plt.show()

    logfile.write("UMAP and PCA plots saved successfully.\n")
    logfile.write("FastHigashi model training completed successfully.\n")
    logfile.close()

