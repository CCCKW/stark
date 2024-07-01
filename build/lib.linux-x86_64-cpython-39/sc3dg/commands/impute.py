import numpy as np
from scipy.sparse import csr_matrix, diags, eye
from scipy.sparse.linalg import norm
from scipy.ndimage import gaussian_filter
from functools import reduce
import cooler
import logging
import pandas as pd
import click
import os
from joblib import Parallel, delayed

def calculate_sparsity(sparse_matrix):
    """Calculate the sparsity of a sparse matrix."""
    rows, cols = sparse_matrix.shape
    sparsity = sparse_matrix.nnz / (rows * cols)
    return sparsity

def perform_random_walk(matrix_P, restart_prob, tolerance):
    """Perform Random Walk with Restart (RWR) on a sparse matrix."""
    if restart_prob == 1:
        return matrix_P

    num_genes = matrix_P.shape[0]
    identity_matrix = eye(num_genes, dtype=np.float32)
    transition_matrix = matrix_P.copy()

    for i in range(100):
        new_transition_matrix = matrix_P.dot(transition_matrix * (1 - restart_prob) + restart_prob * identity_matrix)
        delta = norm(transition_matrix - new_transition_matrix)
        transition_matrix = new_transition_matrix.copy()
        if delta < tolerance:
            break

    return transition_matrix

def make_cooler_bins_chr(chr_name, length, resolution):
    total_bin = []
    for i in range(0, length, resolution):
        total_bin.append([chr_name, i, min(length, i+resolution)])
    df = pd.DataFrame(total_bin, columns=['chrom', 'start', 'end'])
    return df

# 其他杂染色体可以自行添加，字典里有就会读取没有就跳过
chr_list = ['chr' + str(i) for i in [k for k in range(1, 23)] + ['X', 'Y']]

def merge_segement(hic_dict,
                   save_path,
                   resolution, #分辨率要已知
                   chr_list = chr_list,
                   scale_size = 10000, # 概率矩阵上整体放大倍数
                  ):
    total_matrix_list = []
    total_bin_list = []
    stack = 0
    for chr_name in chr_list:
        if chr_name in hic_dict:
            bin_df = make_cooler_bins_chr(chr_name,
                                          length=hic_dict[chr_name].shape[0] * resolution,
                                          resolution=resolution)
            total_bin_list.append(bin_df)
            # 先算单独的坐标，再加stack
            chr_matrix = hic_dict[chr_name]
            rows, cols = np.nonzero(chr_matrix)
            #print(rows, cols)
            counts = np.array(hic_dict[chr_name][rows, cols]) * scale_size
            rows += stack; cols += stack
            matrix_sp = pd.DataFrame(np.column_stack((rows, cols, counts[0])))
            matrix_sp = matrix_sp[matrix_sp.iloc[:, 0] <= matrix_sp.iloc[:, 1]]
            total_matrix_list.append(matrix_sp)
            stack += bin_df.shape[0]

    bin_df_total = reduce(lambda x, y: pd.concat([x,y], axis = 0), total_bin_list)
    bin_df_total = bin_df_total.reset_index(drop=True)

    sp_df_total = reduce(lambda x, y: pd.concat([x,y], axis = 0), total_matrix_list)
    sp_df_total = sp_df_total.reset_index(drop=True)
    sp_df_total.columns = ['bin1_id', 'bin2_id', 'count']
    cooler.create_cooler(cool_uri=save_path,
                         bins=bin_df_total, pixels=sp_df_total)



def deal_imoute( mcool, resolution, output, use_logscale, padding, std_dev, restart_prob, tolerance):
    hic_data = cooler.Cooler(f"{mcool}::/resolutions/{resolution}")
    rwr_results = {}

    for chromosome in hic_data.chromnames:
        contact_matrix = hic_data.matrix(balance=False, sparse=True).fetch(chromosome)

        # Log transform
        if use_logscale:
            contact_matrix.data = np.log2(contact_matrix.data + 1)

        # Remove diagonal before convolution
        contact_matrix = contact_matrix - diags(contact_matrix.diagonal())

        # Apply Gaussian filter
        if padding > 0:
            contact_matrix = gaussian_filter(contact_matrix.astype(np.float32).toarray(), std_dev, order=0, mode='mirror', truncate=padding)
            contact_matrix = csr_matrix(contact_matrix)

        # Remove diagonal before RWR
        contact_matrix = contact_matrix - diags(contact_matrix.diagonal())

        # Prepare matrix for RWR
        matrix_B = contact_matrix + diags((contact_matrix.sum(axis=0).A.ravel() == 0).astype(int))
        degree_matrix = diags(1 / matrix_B.sum(axis=0).A.ravel())
        transition_probability_matrix = degree_matrix.dot(matrix_B).astype(np.float32)
        imputed_matrix = perform_random_walk(transition_probability_matrix, restart_prob, tolerance)

        # Normalize
        imputed_matrix += imputed_matrix.T
        degree = imputed_matrix.sum(axis=0).A.ravel()
        degree[degree == 0] = 1
        normalization_matrix = diags(1 / np.sqrt(degree))
        imputed_matrix = normalization_matrix.dot(imputed_matrix).dot(normalization_matrix)

        # Mask the lower triangle of the imputed matrix
        upper_triangle_indices = np.triu_indices(imputed_matrix.shape[0], 0)
        upper_triangle_mask = csr_matrix((np.ones(len(upper_triangle_indices[0])), (upper_triangle_indices[0], upper_triangle_indices[1])),
                                        imputed_matrix.shape, dtype=np.float32)
        imputed_matrix = imputed_matrix.tocsr().multiply(upper_triangle_mask)

        # Store RWR data
        rwr_results[chromosome] = imputed_matrix

    # Save the RWR results to cool output path
    merge_segement(rwr_results, save_path= output, resolution=resolution)



@click.command('impute', short_help='impute the Hi-C data')
@click.option('--mcool-path',  help='cool file path',required=True, default=None)
@click.option('--resolution', help='resolution',required=True, default=None, type=int)
@click.option('--output', help='output file path',required=True, default=None)
@click.option('--use_logscale', help='use logscale',required=False, default=False, type=bool)
@click.option('--padding', help='padding',required=False, default=1, type=int)
@click.option('--std_dev', help='std_dev',required=False, default=1, type=float)
@click.option('--restart_prob', help='restart_prob',required=False, default=0.5, type=float)
@click.option('--tolerance', help='tolerance',required=False, default=1e-6, type=float)
@click.option('--nthread', help='number of threads',required=False, default=4, type=int)
def impute(mcool_path, 
            resolution, 
            output, 
            use_logscale=False, 
            padding=1, 
            std_dev=1, 
            restart_prob=0.5,
        tolerance=1e-6,
        nthread=4):
    """
    Impute missing values in Hi-C data using Random Walk with Restart (RWR).

    Args:
        cool_file_url (str): URL or path to the .cool file.
        resolution (int): Resolution of the Hi-C data.
        output_file_path (str): Path to save the imputed matrices.
        use_logscale (bool): Whether to log-transform the data.
        padding (int): Padding for Gaussian filter.
        std_dev (float): Standard deviation for Gaussian filter.
        restart_prob (float): Restart probability for RWR.
        tolerance (float): Convergence tolerance for RWR.
        window_size_bp (int): Size of the window for local smoothing in base pairs.
        step_size_bp (int): Step size for sliding window in base pairs.
        output_distance_bp (int): Output distance for storing matrices in base pairs.
        min_cutoff (int): Minimum cutoff for values in the matrix.
    """
    task = []
    for mc in os.listdir(mcool_path):
        if mc.endswith('mcool'):
            task.append( mc)
            
    print(task)
    Parallel(n_jobs=nthread)(delayed(deal_imoute)(mcool=mcool_path + '/' + mc, 
                                                resolution=resolution, 
                                                output=output + '/' + mc.split('.')[0] + '_impute.mcool',
                                                use_logscale=use_logscale, 
                                                padding=padding, 
                                                std_dev=std_dev, 
                                                restart_prob=restart_prob,
                                                tolerance=tolerance) for mc in task)
    




