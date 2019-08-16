
import numpy as np
from sklearn.preprocessing import LabelEncoder
from anndata import AnnData
from typing import Optional
from time import time

def calinski_harabasz_score(X, labels):
    """Compute the Calinski and Harabasz score.
    It is also known as the Variance Ratio Criterion.
    The score is defined as ratio between the within-cluster dispersion and
    the between-cluster dispersion.
    Read more in the :ref:`User Guide <calinski_harabasz_index>`.
    Parameters
    ----------
    X : array-like, shape (``n_samples``, ``n_features``)
        List of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (``n_samples``,)
        Predicted labels for each sample.
    Returns
    -------
    score : float
        The resulting Calinski-Harabasz score.
    References
    ----------
    .. [1] `T. Calinski and J. Harabasz, 1974. "A dendrite method for cluster
       analysis". Communications in Statistics
       <https://www.tandfonline.com/doi/abs/10.1080/03610927408827101>`_
    """
    #X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)

    n_samples, _ = X.shape
    n_labels = len(le.classes_)

    #check_number_of_labels(n_labels, n_samples)

    extra_disp, intra_disp = 0., 0.
    mean = np.mean(X, axis=0)
    for k in range(n_labels):
        cluster_k = X[labels == k]
        mean_k = np.mean(cluster_k, axis=0)
        extra_disp += len(cluster_k) * np.sum((mean_k - mean) ** 2)
        intra_disp += np.sum((cluster_k - mean_k) ** 2)

    return (1. if intra_disp == 0. else
            extra_disp * (n_samples - n_labels) /
            (intra_disp * (n_labels - 1.)))
set_n = 5
def parameter_suggestion(
    adata: AnnData,
    set_n: int = set_n,
    Bin_n: Optional[list] = None,
    Eps: Optional[list] = None,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    t0=time()
    adata = adata.copy() if copy else adata
    Bin_n = Bin_n if Bin_n else [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    Eps = Eps if Eps else [1.2, 1.6, 1.9, 2.1, 2.3, 2.7]
    feature_data = adata.obsm['X_pca'].astype('Float64')
    CHixNobs_values = {}
    for bin_n in Bin_n:
        for eps in Eps:
            #print("0 runing time: "+ str(round(time()-t1,3)))
            sc_FlowGrid(adata,bin_n,eps)
            #print("1 runing time: "+ str(round(time()-t1,3)))
            label_data = adata.obs['binN_'+str(bin_n)+'_eps_'+ str(eps)+'_FlowGrid'].tolist()
            #print("2 runing time: "+ str(round(time()-t1,3)))
            #return label_data
            CHixNobs_values['binN_'+str(bin_n)+'_eps_'+str(eps)+'_FlowGrid'] = calinski_harabasz_score(feature_data, label_data)\
                                                                        * len(set(label_data))
            
            
    maxn_CHixNobs_values = sorted(CHixNobs_values, key=CHixNobs_values.get, reverse=True)[:set_n]
    print(str(set_n)+ " sets of parameters are recommended.\n" +"Suggestion completed in : "+ str(round(time()-t0,3)) + " seconds.")

    return maxn_CHixNobs_values
