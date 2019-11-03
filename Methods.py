import random
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, KernelPCA, FastICA, FactorAnalysis, TruncatedSVD, SparsePCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA


def get_average_cell(cells):
    return np.mean (cells, axis=0)


def different(cell1, cell2):
    return np.linalg.norm (cell1 - cell2)


def gene_MatchLimit_check(cell1, cell2, limit):
    total_gene_num = len (cell1)
    judge1 = (cell1 == 0)
    judge2 = (cell2 == 0)
    judge = judge1 + judge2
    effective_gene_num = len (judge) - np.sum (judge)
    if effective_gene_num >= total_gene_num * limit:
        return True
    else:
        return False


def disMatrix_to_bfMatrix(disMatrix, clusters):
    n = len (clusters)
    to_clusters = ["to_" + c for c in clusters]
    bfMatrix = pd.DataFrame (np.zeros (shape=(n, n)), index=to_clusters, columns=clusters)
    for cluster1 in clusters:
        for cluster2 in [c2 for c2 in clusters if c2 != cluster1]:
            med = np.median (disMatrix[cluster1][cluster2])
            probability_obj = sum (disMatrix[cluster1][cluster2] <= med) / len (disMatrix[cluster1][cluster2])
            probability_refs = probability_obj
            print ("{} to {}".format (cluster1, cluster2))
            print ("probability_obj: {}".format (probability_obj))
            for cluster3 in [c3 for c3 in clusters if c3 != cluster1 and c3 != cluster2]:
                probability_refs += sum (disMatrix[cluster1][cluster3] <= med) / len (disMatrix[cluster1][cluster3])
                print ("{}:".format (cluster3), end=" ")
                print (sum (disMatrix[cluster1][cluster3] <= med) / len (disMatrix[cluster1][cluster3]))
            bfMatrix[cluster1]["to_" + cluster2] = probability_obj / probability_refs
            print (probability_obj / probability_refs)
            print ("")
        print (bfMatrix[cluster1])
        print ("")
    return bfMatrix


def disMatrix_to_meanMatrix(disMatrix, clusters):
    n = len (clusters)
    to_clusters = ["to_" + c for c in clusters]
    bfMatrix = pd.DataFrame (np.zeros (shape=(n, n)), index=to_clusters, columns=clusters)
    for i in range (len (clusters) - 1):
        cluster1 = clusters[i]
        for j in range (i + 1, len (clusters)):
            cluster2 = clusters[j]
            med = np.median (disMatrix[cluster1][cluster2])
            bfMatrix[cluster1]["to_" + cluster2] = bfMatrix[cluster2]["to_" + cluster1] = med
    return bfMatrix


def corrMatrix_to_bfMatrix(corrMatrix, clusters):
    n = len (clusters)
    to_clusters = ["to_" + c for c in clusters]
    bfMatrix = pd.DataFrame (np.zeros (shape=(n, n)), index=to_clusters, columns=clusters)
    for cluster1 in clusters:
        for cluster2 in [c2 for c2 in clusters if c2 != cluster1]:
            med = np.median (corrMatrix[cluster1][cluster2])
            probability_obj = sum (corrMatrix[cluster1][cluster2] >= med) / len (corrMatrix[cluster1][cluster2])
            probability_refs = probability_obj
            print ("{} to {}".format (cluster1, cluster2))
            print ("probability_obj: {}".format (probability_obj))
            for cluster3 in [c3 for c3 in clusters if c3 != cluster1 and c3 != cluster2]:
                probability_refs += sum (corrMatrix[cluster1][cluster3] >= med) / len (corrMatrix[cluster1][cluster3])
                print ("{}:".format (cluster3), end=" ")
                print (sum (corrMatrix[cluster1][cluster3] <= med) / len (corrMatrix[cluster1][cluster3]))
            bfMatrix[cluster1]["to_" + cluster2] = probability_obj / probability_refs
            print (probability_obj / probability_refs)
            print ("")
        print (bfMatrix[cluster1])
        print ("")
    return bfMatrix


def corrMatrix_to_meanMatrix(corrMatrix, clusters):
    n = len (clusters)
    to_clusters = ["to_" + c for c in clusters]
    bfMatrix = pd.DataFrame (np.zeros (shape=(n, n)), index=to_clusters, columns=clusters)
    for i in range (len (clusters)):
        cluster1 = clusters[i]
        bfMatrix[cluster1]["to_" + cluster1] = 1
        if i < len (clusters) - 1:
            for j in range (i + 1, len (clusters)):
                cluster2 = clusters[j]
                med = np.median (corrMatrix[cluster1][cluster2])
                bfMatrix[cluster1]["to_" + cluster2] = bfMatrix[cluster2]["to_" + cluster1] = med
    return bfMatrix


def get_distance(cluster1, cluster2, average_number, caculation_number, match_limit=0):
    if match_limit == 0:
        num1 = len (cluster1)
        num2 = len (cluster2)
        if average_number == 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = cluster1[random.randint (0, num1 - 1)]
                random_cell2 = cluster2[random.randint (0, num2 - 1)]
                diff = different (random_cell1, random_cell2) / len (random_cell1)
                result.append (diff)
                turn = turn + 1
            return np.array (result)

        elif average_number > 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = get_average_cell (cluster1[np.array (random.sample (range (0, num1), average_number))])
                random_cell2 = get_average_cell (cluster2[np.array (random.sample (range (0, num2), average_number))])
                diff = different (random_cell1, random_cell2) / len (random_cell1)
                result.append (diff)
                turn = turn + 1
            return np.array (result)

    elif match_limit > 0:
        num1 = len (cluster1)
        num2 = len (cluster2)
        if average_number == 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = cluster1[random.randint (0, num1 - 1)]
                random_cell2 = cluster2[random.randint (0, num2 - 1)]
                if gene_MatchLimit_check (random_cell1, random_cell2, match_limit):
                    diff = different (random_cell1, random_cell2) / len (random_cell1)
                    result.append (diff)
                    turn = turn + 1
            return np.array (result)

        elif average_number > 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = get_average_cell (cluster1[np.array (random.sample (range (0, num1), average_number))])
                random_cell2 = get_average_cell (cluster2[np.array (random.sample (range (0, num2), average_number))])
                if gene_MatchLimit_check (random_cell1, random_cell2, match_limit):
                    diff = different (random_cell1, random_cell2) / len (random_cell1)
                    result.append (diff)
                    turn = turn + 1
            return np.array (result)


def get_corr(cluster1, cluster2, average_number, caculation_number, match_limit):
    if match_limit == 0:
        num1 = len (cluster1)
        num2 = len (cluster2)
        if average_number == 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = cluster1[random.randint (0, num1 - 1)]
                random_cell2 = cluster2[random.randint (0, num2 - 1)]
                corr = np.corrcoef (random_cell1, random_cell2)[0][1]
                result.append (corr)
                turn = turn + 1
            return np.array (result)

        elif average_number > 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = get_average_cell (cluster1[np.array (random.sample (range (0, num1), average_number))])
                random_cell2 = get_average_cell (cluster2[np.array (random.sample (range (0, num2), average_number))])
                corr = np.corrcoef (random_cell1, random_cell2)[0][1]
                result.append (corr)
                turn = turn + 1
            return np.array (result)

    elif match_limit > 0:
        num1 = len (cluster1)
        num2 = len (cluster2)
        if average_number == 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = cluster1[random.randint (0, num1 - 1)]
                random_cell2 = cluster2[random.randint (0, num2 - 1)]
                if gene_MatchLimit_check (random_cell1, random_cell2, match_limit):
                    corr = np.corrcoef (random_cell1, random_cell2)[0][1]
                    result.append (corr)
                    turn = turn + 1
            return np.array (result)

        elif average_number > 1:
            result = []
            turn = 1
            while turn <= caculation_number:
                random_cell1 = get_average_cell (cluster1[np.array (random.sample (range (0, num1), average_number))])
                random_cell2 = get_average_cell (cluster2[np.array (random.sample (range (0, num2), average_number))])
                if gene_MatchLimit_check (random_cell1, random_cell2, match_limit):
                    corr = np.corrcoef (random_cell1, random_cell2)[0][1]
                    result.append (corr)
                    turn = turn + 1
            return np.array (result)


def get_matrix_dist(data, lab, clusters, average_number, caculation_number, match_limit=0):
    n = len (clusters)
    matrix = {clusters[p]: {clusters[q]: [] for q in range (p + 1, n)} for p in range (n)}
    for i in range (n):
        for j in range (i + 1, n):
            cluster1 = data[lab[clusters[i]]]
            cluster2 = data[lab[clusters[j]]]
            distances = get_distance (cluster1, cluster2, average_number, caculation_number, match_limit)
            matrix[clusters[i]][clusters[j]] = matrix[clusters[j]][clusters[i]] = distances
    return matrix


def get_matrix_corr(data, lab, clusters, average_number, caculation_number, match_limit=0):
    n = len (clusters)
    matrix = {clusters[p]: {clusters[q]: [] for q in range (p + 1, n)} for p in range (n)}
    for i in range (n):
        for j in range (i + 1, n):
            cluster1 = data[lab[clusters[i]]]
            cluster2 = data[lab[clusters[j]]]
            corr = get_corr (cluster1, cluster2, average_number, caculation_number, match_limit)
            matrix[clusters[i]][clusters[j]] = matrix[clusters[j]][clusters[i]] = corr
    return matrix
