# -*- coding: utf-8 -*-
import umap
import numpy as np
from sklearn.manifold import MDS, TSNE, LocallyLinearEmbedding, Isomap


# PCA :---------------------------------------------------

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










