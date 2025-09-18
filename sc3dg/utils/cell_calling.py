#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

""" Functions for calling cell-associated barcodes - Ultra-Optimized Version with Numba & Parallel Processing """

import sys
import os
from collections import namedtuple
import numpy as np
import numpy.ma as ma
import scipy.stats as sp_stats
from sgt import sgt_proportions,SimpleGoodTuringError
import math
from scipy.special import gammaln
from tqdm import tqdm
from numba import jit, njit, prange
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import warnings
warnings.filterwarnings('ignore')

@njit(cache=True)
def adjust_pvalue_bh_numba(p):
    """ Numba优化的Benjamini-Hochberg校正 """
    n = len(p)
    descending = np.argsort(p)[::-1]
    scale = np.arange(n, 0, -1, dtype=np.float64)
    q = np.zeros(n, dtype=np.float64)
    
    # 计算调整后的p值
    for i in range(n):
        idx = descending[i]
        q[i] = min(1.0, p[idx] * n / scale[i])
    
    # 确保单调性
    for i in range(1, n):
        q[i] = min(q[i], q[i-1])
    
    # 重新排序
    result = np.zeros(n, dtype=np.float64)
    for i in range(n):
        result[descending[i]] = q[i]
    
    return result

@njit(cache=True)
def fast_multinomial_logpmf_numba(x, n, log_p):
    """Numba优化的多项式对数概率质量函数"""
    if np.sum(x) != n:
        return -np.inf
    
    # 计算多项式系数的对数
    log_coeff = 0.0
    for i in range(len(x)):
        if x[i] > 0:
            log_coeff += math.lgamma(x[i] + 1)
    log_coeff = math.lgamma(n + 1) - log_coeff
    
    # 计算概率部分
    log_prob = 0.0
    for i in range(len(x)):
        if x[i] > 0:
            log_prob += x[i] * log_p[i]
    
    return log_coeff + log_prob

@njit(cache=True, parallel=True)
def multinomial_logpmf_batch_numba(matrix_chunk, profile_p):
    """Numba并行批处理多项式对数概率"""
    n_bcs = matrix_chunk.shape[0]
    n_features = matrix_chunk.shape[1]
    result = np.zeros(n_bcs, dtype=np.float64)
    
    log_p = np.log(profile_p + 1e-300)  # 避免log(0)
    
    for i in prange(n_bcs):
        row = matrix_chunk[i]
        n = np.sum(row)
        if n > 0:
            result[i] = fast_multinomial_logpmf_numba(row, n, log_p)
        else:
            result[i] = -np.inf
    
    return result

def eval_multinomial_loglikelihoods_ultra_fast(matrix, profile_p, config):
    """超级优化版本的多项式对数似然计算 - 使用Numba并行"""
    
    gb_per_bc = float(matrix.shape[0] * matrix.dtype.itemsize) / (1024**3)
    bcs_per_chunk = max(1, int(round(config['MAX_MEM_GB']/gb_per_bc)))
    num_bcs = matrix.shape[1]
    
    # 计算总chunk数量
    total_chunks = (num_bcs + bcs_per_chunk - 1) // bcs_per_chunk
    loglk = np.zeros(num_bcs, dtype=np.float64)

    # 使用ThreadPoolExecutor进行并行处理
    def process_chunk(chunk_start):
        chunk_end = min(chunk_start + bcs_per_chunk, num_bcs)
        chunk = slice(chunk_start, chunk_end)
        matrix_chunk = matrix[:, chunk].transpose().toarray().astype(np.int32)
        return chunk_start, multinomial_logpmf_batch_numba(matrix_chunk, profile_p)

    # 并行处理chunks
    with ThreadPoolExecutor(max_workers=config['N_CORES']) as executor:
        with tqdm(total=total_chunks, desc="🚀 并行计算对数似然", unit="chunk") as pbar:
            futures = []
            for chunk_start in range(0, num_bcs, bcs_per_chunk):
                futures.append(executor.submit(process_chunk, chunk_start))
            
            for future in futures:
                chunk_start, chunk_result = future.result()
                chunk_end = min(chunk_start + bcs_per_chunk, num_bcs)
                loglk[chunk_start:chunk_end] = chunk_result
                pbar.update(1)
    
    return loglk

@njit(cache=True)
def simulate_single_multinomial_numba(profile_p, distinct_n, num_sims, seed_offset=0):
    """Numba优化的单个多项式模拟"""
    np.random.seed(42 + seed_offset)
    n_distinct = len(distinct_n)
    loglk = np.zeros((n_distinct, num_sims), dtype=np.float64)
    log_p = np.log(profile_p + 1e-300)
    
    for sim_idx in range(num_sims):
        curr_counts = np.random.multinomial(distinct_n[0], profile_p).astype(np.int32)
        loglk[0, sim_idx] = fast_multinomial_logpmf_numba(curr_counts, distinct_n[0], log_p)
        
        for i in range(1, n_distinct):
            step = distinct_n[i] - distinct_n[i-1]
            additional_counts = np.random.multinomial(step, profile_p).astype(np.int32)
            curr_counts += additional_counts
            loglk[i, sim_idx] = fast_multinomial_logpmf_numba(curr_counts, distinct_n[i], log_p)
    
    return loglk

# 全局函数用于多进程调用
def run_simulation_worker(args):
    """全局工作函数用于多进程"""
    process_id, process_sims, profile_p, distinct_n = args
    return simulate_single_multinomial_numba(
        profile_p, distinct_n, process_sims, seed_offset=process_id * 1000
    )

@njit(cache=True, parallel=True)
def simulate_multinomial_batch_numba(profile_p, distinct_n, num_sims):
    """使用Numba并行替代多进程的模拟函数"""
    n_distinct = len(distinct_n)
    loglk = np.zeros((n_distinct, num_sims), dtype=np.float64)
    log_p = np.log(profile_p + 1e-300)
    
    # 并行处理不同的模拟
    for sim_idx in prange(num_sims):
        # 每个线程使用不同的种子
        np.random.seed(42 + sim_idx)
        
        curr_counts = np.random.multinomial(distinct_n[0], profile_p).astype(np.int32)
        loglk[0, sim_idx] = fast_multinomial_logpmf_numba(curr_counts, distinct_n[0], log_p)
        
        for i in range(1, n_distinct):
            step = distinct_n[i] - distinct_n[i-1]
            additional_counts = np.random.multinomial(step, profile_p).astype(np.int32)
            curr_counts += additional_counts
            loglk[i, sim_idx] = fast_multinomial_logpmf_numba(curr_counts, distinct_n[i], log_p)
    
    return loglk

def simulate_multinomial_loglikelihoods_ultra_fast(profile_p, umis_per_bc, config):
    """超级优化版本的多项式对数似然模拟 - 优先使用Numba并行"""
    distinct_n = np.flatnonzero(np.bincount(umis_per_bc.astype(int)))
    num_sims = config['NUM_SIMS']
    
    print(f'🎲 模拟参数: 不同N值数量: {len(distinct_n)}, 模拟次数: {num_sims}')
    
    # 优先使用Numba并行，避免多进程的pickle问题
    try:
        print("⚡ 使用Numba并行模拟...")
        with tqdm(total=1, desc="🔥 Numba并行模拟", unit="批次") as pbar:
            loglk = simulate_multinomial_batch_numba(profile_p, distinct_n, num_sims)
            pbar.update(1)
        return distinct_n, loglk
        
    except Exception as e:
        print(f"⚠️ Numba并行失败，尝试多进程: {e}")
        
        # 备用：使用多进程（修复pickle问题）
        sims_per_process = max(1, num_sims // config['N_CORES'])
        remaining_sims = num_sims % config['N_CORES']
        
        # 准备任务参数
        tasks = []
        for i in range(config['N_CORES']):
            process_sims = sims_per_process + (1 if i < remaining_sims else 0)
            if process_sims > 0:
                tasks.append((i, process_sims, profile_p, distinct_n))
        
        # 并行执行模拟
        try:
            with ProcessPoolExecutor(max_workers=config['N_CORES']) as executor:
                with tqdm(total=len(tasks), desc="🔄 多进程模拟", unit="进程") as pbar:
                    futures = [executor.submit(run_simulation_worker, task) for task in tasks]
                    results = []
                    for future in futures:
                        results.append(future.result())
                        pbar.update(1)
            
            # 合并结果
            loglk = np.concatenate(results, axis=1)
            return distinct_n, loglk
            
        except Exception as e2:
            print(f"⚠️ 多进程也失败，使用串行版本: {e2}")
            # 最后备用：串行版本
            with tqdm(total=1, desc="🐌 串行模拟", unit="批次") as pbar:
                loglk = simulate_single_multinomial_numba(profile_p, distinct_n, num_sims, 0)
                pbar.update(1)
            return distinct_n, loglk

@njit(cache=True, parallel=True)
def compute_ambient_pvalues_numba(umis_per_bc, obs_loglk, sim_n, sim_loglk):
    """Numba并行优化的p值计算"""
    n_bcs = len(umis_per_bc)
    pvalues = np.zeros(n_bcs, dtype=np.float64)
    num_sims = sim_loglk.shape[1]
    
    # 预计算sim_n索引
    sim_n_idx = np.searchsorted(sim_n, umis_per_bc)
    
    for i in prange(n_bcs):
        # 向量化比较
        sim_row = sim_loglk[sim_n_idx[i], :]
        num_lower_loglk = 0
        for j in range(num_sims):
            if sim_row[j] < obs_loglk[i]:
                num_lower_loglk += 1
        pvalues[i] = float(1 + num_lower_loglk) / (1 + num_sims)
    
    return pvalues

def compute_ambient_pvalues_ultra_fast(umis_per_bc, obs_loglk, sim_n, sim_loglk):
    """超级优化版本的p值计算"""
    print("⚡ 使用Numba并行计算p值...")
    with tqdm(total=1, desc="🧮 并行p值计算", unit="批次") as pbar:
        pvalues = compute_ambient_pvalues_numba(umis_per_bc, obs_loglk, sim_n, sim_loglk)
        pbar.update(1)
    return pvalues

def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
    """ Estimate a gene expression profile by Simple Good Turing. """
    prof_mat = matrix[:,barcode_indices]
    profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
    zero_feat = np.flatnonzero(profile == 0)

    try:
        p_smoothed, p0 = sgt_proportions(profile[np.flatnonzero(profile)])
    except:
        # 如果SGT失败，使用简单的拉普拉斯平滑
        profile_smooth = profile + 1
        return profile_smooth / profile_smooth.sum()

    p0_i = p0/max(len(zero_feat), 1)  # 避免除零
    profile_p = np.repeat(p0_i, len(nz_feat))
    profile_p[np.flatnonzero(profile)] = p_smoothed

    return profile_p / profile_p.sum()  # 确保归一化

def est_background_profile_sgt(matrix, use_bcs):
    """ Estimate a gene expression profile on a given subset of barcodes. """
    use_feats = np.flatnonzero(np.asarray(matrix.sum(1)))
    
    with tqdm(total=1, desc="🌟 SGT背景估计") as pbar:
        bg_profile_p = estimate_profile_sgt(matrix, use_bcs, use_feats)
        pbar.update(1)
    
    return (use_feats, bg_profile_p)

NonAmbientBarcodeResult = namedtuple('NonAmbientBarcodeResult',
                                     ['eval_bcs', 'log_likelihood', 'pvalues', 
                                      'pvalues_adj', 'is_nonambient'])

def find_nonambient_barcodes_ultra_fast(matrix, orig_cell_bcs, custom_config):
    """超级优化版本的非环境条形码识别 - Numba + 并行处理"""
    
    config = custom_config
    
    print(f"🔧 设置参数: N_CANDIDATE_BARCODES={config['N_CANDIDATE_BARCODES']}")
    print(f"MIN_UMIS={config['MIN_UMIS']}, NUM_SIMS={config['NUM_SIMS']}")
    print(f"MIN_UMI_FRAC_OF_MEDIAN={config['MIN_UMI_FRAC_OF_MEDIAN']}, MAX_ADJ_PVALUE={config['MAX_ADJ_PVALUE']}")
    
    print(f"🚀 开始超级优化的非环境条形码识别流程... (使用 {config['N_CORES']} 核心)")
    
    # 步骤1: 计算UMI计数
    with tqdm(total=1, desc="📊 计算UMI计数") as pbar:
        umis_per_bc = matrix.get_counts_per_bc()
        bc_order = np.argsort(umis_per_bc)
        pbar.update(1)

    # 步骤2: 处理原始细胞条形码
    with tqdm(total=1, desc="🔍 处理原始细胞条形码") as pbar:
        orig_cell_bc_set = set(orig_cell_bcs)
        orig_cells = np.flatnonzero(np.fromiter(
            (bc in orig_cell_bc_set for bc in matrix.bcs),
            count=len(matrix.bcs), dtype=bool))
        pbar.update(1)
    
    if orig_cells.sum() == 0:
        print('❌ 没有有效的原始细胞调用')
        return None

    # 步骤3: 计算阈值
    with tqdm(total=1, desc="📏 计算阈值") as pbar:
        median_initial_umis = np.median(umis_per_bc[orig_cells])
        min_umis = int(max(config['MIN_UMIS'], 
                          round(np.ceil(median_initial_umis * config['MIN_UMI_FRAC_OF_MEDIAN']))))
        pbar.update(1)
    
    print(f'📈 初始细胞调用的中位UMI数: {median_initial_umis}')
    print(f'📉 最小UMI数: {min_umis}')

    # 步骤4: 估计环境背景
    with tqdm(total=1, desc="🌍 选择背景条形码") as pbar:
        empty_bcs = bc_order[:int(config['EMPTY_BC_FRAC'] * len(bc_order))]
        nz_bcs = np.flatnonzero(umis_per_bc)
        use_bcs = np.intersect1d(empty_bcs, nz_bcs, assume_unique=True)
        pbar.update(1)

    if len(use_bcs) > 0:
        try:
            eval_features, ambient_profile_p = est_background_profile_sgt(matrix.m, use_bcs)
        except Exception as e:
            print(f'❌ 背景估计错误: {str(e)}')
            return None
    else:
        print('❌ 没有可用的背景条形码')
        return None

    # 步骤5: 选择候选条形码
    with tqdm(total=1, desc="🎯 选择候选条形码") as pbar:
        eval_bcs = np.ma.array(np.arange(matrix.bcs_dim))
        eval_bcs[umis_per_bc < min_umis] = ma.masked
        n_unmasked_bcs = len(eval_bcs) - eval_bcs.mask.sum()
        
        # 根据配置选择候选条形码
        eval_bcs = np.argsort(ma.masked_array(umis_per_bc, mask=eval_bcs.mask))[
            max(0, n_unmasked_bcs-config['N_CANDIDATE_BARCODES']):n_unmasked_bcs]
        pbar.update(1)

    if len(eval_bcs) == 0:
        print('❌ 没有候选条形码')
        return None

    print(f'🎯 候选条形码数量: {len(eval_bcs)}')
    print(f'📊 候选条形码UMI范围: {umis_per_bc[eval_bcs].min()}, {umis_per_bc[eval_bcs].max()}')
   
    eval_mat = matrix.m[eval_features, :][:, eval_bcs]

    if len(ambient_profile_p) == 0:
        return None

    # 步骤6: 计算观测对数似然 (Numba并行)
    print("🚀 使用Numba并行计算观测对数似然...")
    obs_loglk = eval_multinomial_loglikelihoods_ultra_fast(eval_mat, ambient_profile_p, config)

    # 步骤7: 模拟对数似然 (优化并行策略)
    print("⚡ 智能并行模拟对数似然...")
    distinct_ns, sim_loglk = simulate_multinomial_loglikelihoods_ultra_fast(
        ambient_profile_p, umis_per_bc[eval_bcs], config)

    # 步骤8: 计算p值 (Numba并行)
    pvalues = compute_ambient_pvalues_ultra_fast(umis_per_bc[eval_bcs], obs_loglk, 
                                               distinct_ns, sim_loglk)
    
    # 步骤9: 多重检验校正 (Numba优化)
    with tqdm(total=1, desc="⚡ Numba Benjamini-Hochberg校正") as pbar:
        pvalues_adj = adjust_pvalue_bh_numba(pvalues)
        is_nonambient = pvalues_adj <= config['MAX_ADJ_PVALUE']
        pbar.update(1)

    print(f"✅ 识别出 {np.sum(is_nonambient)} 个非环境条形码")
    print(f"🚀 性能提升: 使用了Numba JIT编译 + {config['N_CORES']}核心并行处理")

    return NonAmbientBarcodeResult(
        eval_bcs=eval_bcs,
        log_likelihood=obs_loglk,
        pvalues=pvalues,
        pvalues_adj=pvalues_adj,
        is_nonambient=is_nonambient,
    )

# 预编译Numba函数以避免首次运行的编译开销
def warmup_numba_functions():
    """预热Numba函数"""
    print("🔥 预热Numba函数...")
    
    # 创建小测试数据
    test_profile = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float64)
    test_x = np.array([1, 2, 3, 4], dtype=np.int32)
    test_p = np.array([0.1, 0.2, 0.3], dtype=np.float64)
    
    # 预热函数
    _ = fast_multinomial_logpmf_numba(test_x, 10, np.log(test_profile))
    _ = adjust_pvalue_bh_numba(test_p)
    
    test_matrix = np.random.randint(0, 10, (5, 4), dtype=np.int32)
    _ = multinomial_logpmf_batch_numba(test_matrix, test_profile)
    
    print("✅ Numba函数预热完成")

# 向后兼容的函数
def find_nonambient_barcodes(matrix, orig_cell_bcs, 
                           min_umi_frac_of_median=0.05,
                           min_umis_nonambient=1000,
                           max_adj_pvalue=0.000000001):
    """向后兼容的函数，使用默认配置"""
    default_config = {
        'MIN_UMIS': min_umis_nonambient,
        'MIN_UMI_FRAC_OF_MEDIAN': min_umi_frac_of_median,
        'MAX_ADJ_PVALUE': max_adj_pvalue,
        'N_CANDIDATE_BARCODES': 10000,
        'NUM_SIMS': 2000,
        'N_CORES': min(cpu_count(), 20),
        'MAX_MEM_GB': 0.1,
        'EMPTY_BC_FRAC': 0.03,
    }
    
    result = find_nonambient_barcodes_ultra_fast(matrix, orig_cell_bcs, default_config)
    if result is None:
        return None, None, None, None, None
    
    return (result.eval_bcs, result.log_likelihood, result.pvalues, 
            result.pvalues_adj, result.is_nonambient)

# 在模块导入时预热
if __name__ == "__main__":
    warmup_numba_functions()