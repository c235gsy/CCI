import comment_lines
import globalvar as gl
import Methods
import time
import numpy as np
import re
import write_output
from sklearn.decomposition import PCA, KernelPCA, FastICA, FactorAnalysis, TruncatedSVD, SparsePCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA


def sort_key(s):
    if s:
        try:
            c = tuple (map (int, re.findall ("\d+", s)))
        except:
            c = (-1)
        return c


# def auto_KernelPCA(X, n_components):
#     from sklearn.metrics import mean_squared_error
#     best_parameters = ""
#     best_score = 0.0
#     for kernel in ["poly", "rbf", "sigmoid", "cosine"]:
#         kpca = KernelPCA (n_components=n_components, kernel=kernel, fit_inverse_transform=True)
#         X_reduced = kpca.fit_transform (X)
#         X_preimage = kpca.inverse_transform (X_reduced)
#         score = mean_squared_error (X, X_preimage)
#         if score > best_score:
#             best_score = score
#             best_parameters = kernel
#     print (best_parameters, best_score)
#     rbf_pca = KernelPCA (n_components=n_components, kernel=best_parameters, fit_inverse_transform=True)
#     X_reduced = rbf_pca.fit_transform (X)
#     # X_preimage = rbf_pca.inverse_transform(X_reduced)
#     # mean_squared_error(X, X_preimage)
#     return X_reduced


def matrix2dic(matrix):
    Graph = {}
    for index in matrix.index:
        result = matrix.loc[index]
        result = result.sort_values (ascending=False)
        # third = 0
        third = result[2]
        first = result[0]
        next_points = [point[3:] for point in result.index[(result >= third) & (result > 0.5 * first)]]
        # next_points = [point[3:] for point in result.index[(result >= third)]]
        Graph[index] = next_points


def show_dics(dics):
    for key, value in dics.items ():
        print (key, ":", value.keys ())


'''
Corr_Flag, Corr_CalNum, Corr_AvgNum, Corr_MatchLimit = False, 1000, 10, 0.05
Diff_Flag, Diff_CalNum, Diff_AvgNum, Diff_MatchLimit = False, 1000, 10, 0.05
PCA_Flag, PCA_CalNum, PCA_AvgNum, PCA_Set = False, 1000, 30
Gene_Set_Flag, Gene_Set_CalNum, geneListFile, geneSetFile, Gene_Set_Set = False, 1000, "", "", "all"
objName = ""
exprFile = ""
identFile = ""
outputFile = ""
'''

start = time.time ()
comment_lines.get_sets ()
expre_data = np.loadtxt (gl.get_value ("exprFile"), dtype=str)[:, 2:].astype (np.float64)
# expre_data = np.loadtxt (gl.get_value ("exprFile"), dtype=np.float64)
# ident_data = np.loadtxt (gl.get_value ("identFile"), dtype=str, delimiter="$$$$$$$$$$$$$$$")
ident_data_open = open (gl.get_value ("identFile"))
ident_data = np.array ([l.strip () for l in ident_data_open.readlines ()], dtype=str)
ident_data_open.close ()

clusters = sorted (list (set (ident_data)), key=sort_key)

lab = {c: ident_data == c for c in clusters}

out_put_file = open (gl.get_value ("outputFile") + "_infor.txt", "w")
write_output.comment_all (out_put_file, gl.get_value ("comment"), lab, expre_data)
out_put_file.close ()

print (time.time () - start)
bigNum = max ([sum (value) for value in lab.values ()])
lab = {key: np.random.choice (np.where (value)[0], bigNum) for key, value in lab.items ()}
wholeData = np.concatenate ([expre_data[value] for value in lab.values ()], axis=0)
wholeIdent = np.concatenate ([ident_data[value] for value in lab.values ()], axis=0)
print (time.time () - start)

global matrix_Diff, matrix_Corr, matrix_PCA, matrix_KPCA, matrix_ICA, matrix_FA, matrix_TSVD, matrix_LDA, matrix_SPCA
global expre_PCA, expre_KPCA, expre_ICA, expre_FA, expre_TSVD, expre_LDA, expre_SPCA

if gl.get_value ("objName") == "ALL":

    if gl.get_value ("Corr_Flag"):
        matrix_Corr = Methods.get_matrix_corr (data=expre_data, lab=lab, clusters=clusters,
                                               average_number=gl.get_value ("Corr_AvgNum"),
                                               caculation_number=gl.get_value ("Corr_CalNum"),
                                               match_limit=gl.get_value ("Corr_MatchLimit"))
        # print ("get Corr matrix")
        # show_dics(matrix_Corr)
        matrix_Corr_BF = Methods.corrMatrix_to_bfMatrix (matrix_Corr, clusters)
        matrix_Corr_mean = Methods.corrMatrix_to_meanMatrix (matrix_Corr, clusters)
        # print ("get Corr BF matrix")
        # print (matrix_Corr_BF)
        # print ("get Corr mean matrix")
        # print (matrix_Corr_mean)
        matrix_Corr_BF.to_csv (gl.get_value ("outputFile") + "_Corr_BF.txt", sep='\t', header=True, index=True)
        matrix_Corr_mean.to_csv (gl.get_value ("outputFile") + "_Corr_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("Diff_Flag"):
        matrix_Diff = Methods.get_matrix_dist (data=expre_data, lab=lab, clusters=clusters,
                                               average_number=gl.get_value ("Diff_AvgNum"),
                                               caculation_number=gl.get_value ("Diff_CalNum"),
                                               match_limit=gl.get_value ("Diff_MatchLimit"))
        # print ("get Diff matrix")
        # show_dics(matrix_Diff)
        matrix_Diff_BF = Methods.disMatrix_to_bfMatrix (matrix_Diff, clusters)
        matrix_Diff_mean = Methods.disMatrix_to_meanMatrix (matrix_Diff, clusters)
        # print ("get Diff BF matrix")
        # print (matrix_Diff_BF)
        # print ("get Diff mean matrix")
        # print (matrix_Diff_mean)
        matrix_Diff_BF.to_csv (gl.get_value ("outputFile") + "_Diff_BF.txt", sep='\t', header=True, index=True)
        matrix_Diff_mean.to_csv (gl.get_value ("outputFile") + "_Diff_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("PCA_Flag"):
        pca = PCA (n_components=gl.get_value ("PCA_n_components"), svd_solver="full")
        pca.fit (wholeData)
        expre_PCA = pca.transform (expre_data)
        # print ("get PCA data")
        # np.savetxt("./all/for_guosiyuan/expr_ref_all_PCA.txt", expre_PCA)
        matrix_PCA = Methods.get_matrix_dist (data=expre_PCA, lab=lab, clusters=clusters,
                                              average_number=gl.get_value ("PCA_AvgNum"),
                                              caculation_number=gl.get_value ("PCA_CalNum"))
        # print ("get PCA matrix")
        # show_dics(matrix_PCA)
        matrix_PCA_BF = Methods.disMatrix_to_bfMatrix (matrix_PCA, clusters)
        matrix_PCA_mean = Methods.disMatrix_to_meanMatrix (matrix_PCA, clusters)
        # print ("get PCA BF matrix")
        # print (matrix_PCA_BF)
        # print ("get PCA mean matrix")
        # print (matrix_PCA_mean)
        matrix_PCA_BF.to_csv (gl.get_value ("outputFile") + "_PCA_BF.txt", sep='\t', header=True, index=True)
        matrix_PCA_mean.to_csv (gl.get_value ("outputFile") + "_PCA_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("KPCA_Flag"):
        # if gl.get_value ("KPCA_Kernel") == "auto":
        #     expre_KPCA = auto_KernelPCA (expre_data, n_components=gl.get_value ("KPCA_n_components"))
        # else:
        kpca = KernelPCA (n_components=gl.get_value ("KPCA_n_components"), kernel=gl.get_value ("KPCA_Kernel"))
        kpca.fit (wholeData)
        expre_KPCA = kpca.transform (expre_data)
        # print ("get KPCA data")
        matrix_KPCA = Methods.get_matrix_dist (data=expre_KPCA, lab=lab, clusters=clusters,
                                               average_number=gl.get_value ("KPCA_AvgNum"),
                                               caculation_number=gl.get_value ("KPCA_CalNum"))
        # print ("get KPCA matrix")
        # show_dics(matrix_KPCA)
        matrix_KPCA_BF = Methods.disMatrix_to_bfMatrix (matrix_KPCA, clusters)
        matrix_KPCA_mean = Methods.disMatrix_to_meanMatrix (matrix_KPCA, clusters)
        # print ("get KPCA BF matrix")
        # print (matrix_KPCA_BF)
        # print ("get KPCA mean matrix")
        # print (matrix_KPCA_mean)
        matrix_KPCA_BF.to_csv (gl.get_value ("outputFile") + "_KPCA_BF.txt", sep='\t', header=True, index=True)
        matrix_KPCA_mean.to_csv (gl.get_value ("outputFile") + "_KPCA_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("SPCA_Flag"):
        spca = SparsePCA (n_components=gl.get_value ("SPCA_n_components"))
        spca.fit (wholeData)
        expre_SPCA = spca.transform (expre_data)
        # print ("get SPCA data")
        matrix_SPCA = Methods.get_matrix_dist (data=expre_SPCA, lab=lab, clusters=clusters,
                                               average_number=gl.get_value ("SPCA_AvgNum"),
                                               caculation_number=gl.get_value ("SPCA_CalNum"))
        # print ("get SPCA matrix")
        matrix_SPCA_BF = Methods.disMatrix_to_bfMatrix (matrix_SPCA, clusters)
        matrix_SPCA_mean = Methods.disMatrix_to_meanMatrix (matrix_SPCA, clusters)
        # print ("get SPCA BF matrix")
        # print (matrix_SPCA_BF)
        matrix_SPCA_BF.to_csv (gl.get_value ("outputFile") + "_SPCA_BF.txt", sep='\t', header=True, index=True)
        matrix_SPCA_mean.to_csv (gl.get_value ("outputFile") + "_SPCA_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("ICA_Flag"):
        ica = FastICA (n_components=gl.get_value ("ICA_n_components"))
        ica.fit (wholeData)
        expre_ICA = ica.transform (expre_data)
        # print ("get ICA data")
        matrix_ICA = Methods.get_matrix_dist (data=expre_ICA, lab=lab, clusters=clusters,
                                              average_number=gl.get_value ("ICA_AvgNum"),
                                              caculation_number=gl.get_value ("ICA_CalNum"))
        # print ("get ICA matrix")
        matrix_ICA_BF = Methods.disMatrix_to_bfMatrix (matrix_ICA, clusters)
        matrix_ICA_mean = Methods.disMatrix_to_meanMatrix (matrix_ICA, clusters)
        # print ("get ICA BF matrix")
        # print (matrix_ICA_BF)
        matrix_ICA_BF.to_csv (gl.get_value ("outputFile") + "_ICA_BF.txt", sep='\t', header=True, index=True)
        matrix_ICA_mean.to_csv (gl.get_value ("outputFile") + "_ICA_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("FA_Flag"):
        fa = FactorAnalysis (n_components=gl.get_value ("FA_n_components"))
        fa.fit (wholeData)
        expre_FA = fa.transform (expre_data)
        # print ("get FA data")
        matrix_FA = Methods.get_matrix_dist (data=expre_FA, lab=lab, clusters=clusters,
                                             average_number=gl.get_value ("FA_AvgNum"),
                                             caculation_number=gl.get_value ("FA_CalNum"))
        # print ("get FA matrix")
        matrix_FA_BF = Methods.disMatrix_to_bfMatrix (matrix_FA, clusters)
        matrix_FA_mean = Methods.disMatrix_to_meanMatrix (matrix_FA, clusters)
        # print ("get FA BF matrix")
        # print (matrix_FA_BF)
        matrix_FA_BF.to_csv (gl.get_value ("outputFile") + "_FA_BF.txt", sep='\t', header=True, index=True)
        matrix_FA_mean.to_csv (gl.get_value ("outputFile") + "_FA_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("TSVD_Flag"):
        tsvd = TruncatedSVD (n_components=gl.get_value ("TSVD_n_components"))
        tsvd.fit (wholeData)
        expre_TSVD = tsvd.transform (expre_data)
        # print ("get TSVD data")
        matrix_TSVD = Methods.get_matrix_dist (data=expre_TSVD, lab=lab, clusters=clusters,
                                               average_number=gl.get_value ("TSVD_AvgNum"),
                                               caculation_number=gl.get_value ("TSVD_CalNum"))
        # print ("get TSVD matrix")
        matrix_TSVD_BF = Methods.disMatrix_to_bfMatrix (matrix_TSVD, clusters)
        matrix_TSVD_mean = Methods.disMatrix_to_meanMatrix (matrix_TSVD, clusters)
        # print ("get TSVD BF matrix")
        # print (matrix_TSVD_BF)
        matrix_TSVD_BF.to_csv (gl.get_value ("outputFile") + "_TSVD_BF.txt", sep='\t', header=True, index=True)
        matrix_TSVD_mean.to_csv (gl.get_value ("outputFile") + "_TSVD_mean.txt", sep='\t', header=True, index=True)

    if gl.get_value ("LDA_Flag"):
        if gl.get_value ("LDA_n_components") == "max":
            LDA_n_components = len (clusters) - 1
        else:
            LDA_n_components = gl.get_value ("LDA_n_components")

        lda = LDA (n_components=LDA_n_components)
        lda.fit (wholeData, wholeIdent)
        expre_LDA = lda.transform (expre_data)
        # print ("get LDA data")
        matrix_LDA = Methods.get_matrix_dist (data=expre_LDA, lab=lab, clusters=clusters,
                                              average_number=gl.get_value ("LDA_AvgNum"),
                                              caculation_number=gl.get_value ("LDA_CalNum"))
        # print ("get LDA matrix")
        matrix_LDA_BF = Methods.disMatrix_to_bfMatrix (matrix_LDA, clusters)
        matrix_LDA_mean = Methods.disMatrix_to_meanMatrix (matrix_LDA, clusters)
        # print ("get LDA BF matrix")
        # print (matrix_LDA_BF)
        matrix_LDA_BF.to_csv (gl.get_value ("outputFile") + "_LDA_BF.txt", sep='\t', header=True, index=True)
        matrix_LDA_mean.to_csv (gl.get_value ("outputFile") + "_LDA_mean.txt", sep='\t', header=True, index=True)

# theMatrix.to_csv(gl.get_value("outputFile"), sep='\t', header=True, index=True)
# out_put_file = open(gl.get_value("outputFile"), "w")
# write_output.comment_and_title(out_put_file, gl.get_value("comment"), gl.get_value("objName"), lab, expre_data)
# end = time.time()
# out_put_file.write('Running time: %s Seconds' %(end - start))
# out_put_file.close()
