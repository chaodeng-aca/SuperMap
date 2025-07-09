# -*- coding: utf-8 -*-
import argparse
from sklearn import metrics
import numpy as np
import pandas as pd
import sklearn.neighbors
from sklearn.metrics import f1_score
from scipy.spatial.distance import cdist

with open("tmp_file/true_label.txt", "r") as file:
    vector1 = file.read().split(",")
with open("tmp_file/pred_label.txt", "r") as file:
    vector2 = file.read().split(",")

vector1 = np.array(vector1)
vector2 = np.array(vector2)

acc = sum(vector1 == vector2)/len(vector1)
ARI = metrics.adjusted_rand_score(vector1,vector2)
NMI =  metrics.normalized_mutual_info_score(vector1,vector2)
f1 = f1_score(vector1, vector2, average='macro')

matrix = pd.read_csv("tmp_file/tmp_cord.csv")
matrix = matrix.to_numpy()



def knn_purity_score(knn: np.ndarray, labels: np.ndarray) -> float:
    """Compute the kNN purity score, averaged over all observations.
    For one observation, the purity score is the percentage of
    nearest neighbors that share its label.

    Args:
        knn (np.ndarray):
            The knn, shaped (n_obs, k). The i-th row should contain integers
            representing the indices of the k nearest neighbors.
        labels (np.ndarray):
            The labels, shaped (n_obs)

    Returns:
        float: The purity score.
    """
    # Check the dimensions of the input.
    assert knn.shape[0] == labels.shape[0]

    # Initialize a list of purity scores.
    score = 0

    # Iterate over the observations.
    for i, neighbors in enumerate(knn):

        # Do the neighbors have the same label as the observation?
        matches = labels[neighbors] == labels[i]

        # Add the purity rate to the scores.
        score += np.mean(matches) / knn.shape[0]

    # Return the average purity.
    return score


def embedding_to_knn(
    embedding: np.ndarray, k: float = 0.01, metric: str = "euclidean"
) -> np.ndarray:
    """Convert embedding to knn

    Args:
        embedding (np.ndarray): The embedding (n_obs, n_latent)
        k (int, optional): The number of nearest neighbors. Defaults to 15.
        metric (str, optional): The metric to compute neighbors with. Defaults to "euclidean".

    Returns:
        np.ndarray: The knn (n_obs, k)
    """
    k = max(round(embedding.shape[0] * k), 1)
    # Initialize the knn graph.
    knn = np.zeros((embedding.shape[0], k), dtype=int)

    # Compute pariwise distances between observations.
    distances = cdist(embedding, embedding, metric=metric)

    # Iterate over observations.
    for i in range(distances.shape[0]):

        # Get the `max_neighbors` nearest neighbors.
        knn[i] = distances[i].argsort()[1 : k + 1]

    # Return the knn graph.
    return knn

knn_graph = embedding_to_knn(embedding=matrix)
knn_score = knn_purity_score(knn=knn_graph,labels = np.concatenate((vector1, vector1), axis = 0))

def mean_average_precision(
        x: np.ndarray, y: np.ndarray, neighbor_frac: float = 0.01, **kwargs
) -> float:
    r"""
    Mean average precision

    Parameters
    ----------
    x
        Coordinates
    y
        Cell type labels
    neighbor_frac
        Nearest neighbor fraction
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`

    Returns
    -------
    map
        Mean average precision
    """
    k = max(round(y.shape[0] * neighbor_frac), 1)
    nn = sklearn.neighbors.NearestNeighbors(
        n_neighbors=min(y.shape[0], k + 1), **kwargs
    ).fit(x)
    nni = nn.kneighbors(x, return_distance=False)
    match = np.equal(y[nni[:, 1:]], np.expand_dims(y, 1))
    return np.apply_along_axis(_average_precision, 1, match).mean().item()


def _average_precision(match: np.ndarray) -> float:
    if np.any(match):
        cummean = np.cumsum(match) / (np.arange(match.size) + 1)
        return cummean[match].mean().item()
    return 0.0

MAP = mean_average_precision(matrix,np.concatenate((vector1, vector1), axis = 0))


print('acc',acc,"ARI",ARI,'NMI',NMI,'F1',f1,'KNN_score',knn_score,'MAP',MAP)
