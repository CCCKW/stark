# -*- coding: utf-8 -*-
import umap
import numpy as np
from sklearn.manifold import MDS, TSNE, LocallyLinearEmbedding, Isomap
from scipy.stats import norm
def gaussian_distance(self, matrix1, matrix2, sigma):
        # 计算高斯核矩阵
        kernel1 = np.exp(-np.sum((matrix1[:, np.newaxis] - matrix1) ** 2, axis=2) / (2 * sigma ** 2))
        kernel2 = np.exp(-np.sum((matrix2[:, np.newaxis] - matrix2) ** 2, axis=2) / (2 * sigma ** 2))

        # 计算高斯距离
        distance = np.sum((kernel1 - kernel2) ** 2)

        return distance


def cosine_similarity(self, matrix_A, matrix_B):
        ####计算两个矩阵的余弦相似性
        vector_A = matrix_A.ravel()  # 使用ravel()代替flatten()，无需创建副本
        vector_B = matrix_B.ravel()

        dot_product = np.dot(vector_A, vector_B)
        norm_A = np.linalg.norm(vector_A)
        norm_B = np.linalg.norm(vector_B)

        similarity = dot_product / (norm_A * norm_B)
        return similarity
# PCA :---------------------------------------------------
def scale_matrix(self, matrix, target_size):
        n = matrix.shape[0]  # 原始矩阵的大小
        scale_factor = n / target_size  # 计算缩放比例

        # 构建目标矩阵的行和列索引
        row_indices = np.arange(target_size) * scale_factor
        col_indices = np.arange(target_size) * scale_factor

        # 使用索引进行插值计算
        row_floor = row_indices.astype(int)
        row_ceil = np.ceil(row_indices).astype(int)
        col_floor = col_indices.astype(int)
        col_ceil = np.ceil(col_indices).astype(int)

        # 计算四个邻近像素值
        Q11 = matrix[row_floor[:, None], col_floor]
        Q12 = matrix[row_floor[:, None], col_ceil]
        Q21 = matrix[row_ceil[:, None], col_floor]
        Q22 = matrix[row_ceil[:, None], col_ceil]

        # 双线性插值计算目标矩阵
        target_matrix = (Q11 * (row_ceil[:, None] - row_indices[:, None]) * (col_ceil - col_indices) +
                         Q21 * (row_indices[:, None] - row_floor[:, None]) * (col_ceil - col_indices) +
                         Q12 * (row_ceil[:, None] - row_indices[:, None]) * (col_indices - col_floor) +
                         Q22 * (row_indices[:, None] - row_floor[:, None]) * (col_indices - col_floor))

        return target_matrix
        
def PCAEmbedding(X, dim=2):
    """
    
    Principal components analysis，PCA.
    
    
    Parameters
    ----------
    X : numpy.ndarray
        Coordinates of input data points.
        
    dim : int, optional
        The dimension of the projected points. 
        The default is 2.


    Returns
    -------
    Y : numpy.ndarray
        Coordinates of the projected points.

    """
    
    X=X-np.mean(X,axis=0)
    U, S, V = np.linalg.svd(X, full_matrices=False, compute_uv=True)
    Y = np.dot(V[:dim,:],X.T).T
    return(Y)
    

def MDSEmbedding(X, dim=2):
    """
    Multidimensional Scaling, MDS.

    Parameters
    ----------
    X : numpy.ndarray
        Coordinates of input data points.

    dim : int, optional
        The dimension of the projected points.
        The default is 2.

    Returns
    -------
    Y : numpy.ndarray
        Coordinates of the projected points.

    """
    mds = MDS(n_components=dim, normalized_stress='auto')
    Y = mds.fit_transform(X)

    return Y


def tSNEEmbedding(X, dim=2, perplexity=30.0, learning_rate=200.0, n_iter=1000):
    """
    t-Distributed Stochastic Neighbor Embedding，t-SNE.

    Parameters
    ----------
    X : numpy.ndarray
        Coordinates of input data points.

    dim : int, optional
        The dimension of the projected points.
        The default is 2.

    perplexity : float, optional
        The perplexity parameter for t-SNE.
        The default is 30.0.

    learning_rate : float, optional
        The learning rate for t-SNE optimization.
        The default is 200.0.

    n_iter : int, optional
        The number of iterations for t-SNE optimization.
        The default is 1000.

    Returns
    -------
    Y : numpy.ndarray
        Coordinates of the projected points.
    """
    tsne = TSNE(n_components=dim, perplexity=perplexity, learning_rate=learning_rate, n_iter=n_iter, method='exact')
    Y = tsne.fit_transform(X)
    return Y


def UMAPEmbedding(X, dim=2):
    """
    Uniform Manifold Approximation and Projection (UMAP).

    Parameters
    ----------
    X : numpy.ndarray
        Coordinates of input data points.

    dim : int, optional
        The dimension of the projected points.
        The default is 2.

    Returns
    -------
    Y : numpy.ndarray
        Coordinates of the projected points.
    """

    reducer = umap.UMAP(n_components=dim)
    Y = reducer.fit_transform(X)

    return Y


def LLEEmbedding(X, dim=2):
    """
    Locally Linear Embedding (LLE).

    Parameters
    ----------
    X : numpy.ndarray
        Coordinates of input data points.

    dim : int, optional
        The dimension of the projected points.
        The default is 2.

    Returns
    -------
    Y : numpy.ndarray
        Coordinates of the projected points.
    """

    lle = LocallyLinearEmbedding(n_components=dim)
    Y = lle.fit_transform(X)

    return Y


def IsomapEmbedding(X, dim=2, n_neighbors=5):
    """
    Isometric Mapping (Isomap).

    Parameters
    ----------
    X : numpy.ndarray
        Coordinates of input data points.

    dim : int, optional
        The dimension of the projected points.
        The default is 2.

    n_neighbors : int, optional
        The number of neighbors to consider for each data point.
        The default is 5.

    Returns
    -------
    Y : numpy.ndarray
        Coordinates of the projected points.
    """

    isomap = Isomap(n_components=dim, n_neighbors=n_neighbors)
    Y = isomap.fit_transform(X)

    return Y










