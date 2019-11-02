import sys
import getopt
import globalvar as gl


# global Corr_Flag, Corr_CalNum, Corr_AvgNum, Corr_MatchLimit
# global Diff_Flag, Diff_CalNum, Diff_AvgNum, Diff_MatchLimit
# global PCA_Flag, PCA_CalNum, PCA_Set
# global Gene_Set_Flag, Gene_Set_CalNum, geneListFile, geneSetFile, Gene_Set_Set
# global objName
# global exprFile
# global identFile
# global outputFile

gl._init()


def get_sets():
    if sys.argv[1] != "-h" and sys.argv[1] != "--help":
    #if True:
        comment_file = sys.argv[1]
        #comment_file = "./comment.txt"
        comment = open(comment_file, "r").readlines()
        gl.set_value("comment", comment)

        if len(comment) > 0:

            for com in comment:
                if com[0] != "#":
                    com = com.strip()
                    com = com.split(" ")
                    #print(com)

                    while com[-1] == " ":
                        del com[-1]

                    if com[0] == "Corr":
                        #Corr_Flag = True
                        gl.set_value('Corr_Flag', True)
                        opts, args = getopt.getopt(com[1:], "c:a:m:", ["Corr_CalNum=", "Corr_AvgNum=", "Corr_MatchLimit="])
                        for opt, arg in opts:
                            if opt in ["-c", "--Corr_CalNum"]:
                                #Corr_CalNum = int(arg)
                                gl.set_value("Corr_CalNum", int(arg))
                            if opt in ["-a", "--Corr_AvgNum"]:
                                #Corr_AvgNum = int(arg)
                                gl.set_value("Corr_AvgNum", int(arg))
                            if opt in ["-m", "--Corr_MatchLimit"]:
                                #Corr_MatchLimit = float(arg)
                                gl.set_value ("Corr_MatchLimit", float(arg))

                    if com[0] == "Diff":
                        #Diff_Flag = True
                        gl.set_value("Diff_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:a:m:", ["Diff_CalNum=", "Diff_AvgNum=", "Diff_MatchLimit="])
                        for opt, arg in opts:
                            if opt in ["-c", "--Diff_CalNum"]:
                                #Diff_CalNum = int(arg)
                                gl.set_value("Diff_CalNum", int(arg))
                            if opt in ["-a", "--Diff_AvgNum"]:
                                #Diff_AvgNum = int(arg)
                                gl.set_value("Diff_AvgNum", int(arg))
                            if opt in ["-m", "--Diff_MatchLimit"]:
                                #Diff_MatchLimit = float(arg)
                                gl.set_value("Diff_MatchLimit", float(arg))

                    if com[0] == "PCA":
                        #PCA_Flag = True
                        gl.set_value("PCA_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:", ["PCA_CalNum=","PCA_n_components=","PCA_AvgNum="])
                        for opt, arg in opts:
                            if opt in ["-c", "--PCA_CalNum"]:
                                #PCA_CalNum = int(arg)
                                gl.set_value("PCA_CalNum", int(arg))
                            if opt in ["-a", "--PCA_AvgNum"]:
                                #PCA_CalNum = int(arg)
                                gl.set_value("PCA_AvgNum", int(arg))
                            if opt in ["-n", "--PCA_n_components"]:
                                if float(arg) > 1:
                                    #PCA_Set = int(arg)
                                    gl.set_value("PCA_n_components", int(arg))
                                else:
                                    #PCA_Set = float(arg)
                                    gl.set_value("PCA_n_components", float(arg))

                    if com[0] == "SPCA":
                        #PCA_Flag = True
                        gl.set_value("SPCA_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:", ["SPCA_CalNum=",
                                                                       "SPCA_n_components=", "SPCA_AvgNum="])
                        for opt, arg in opts:
                            if opt in ["-c", "--SPCA_CalNum"]:
                                #PCA_CalNum = int(arg)
                                gl.set_value("SPCA_CalNum", int(arg))
                            if opt in ["-a", "--SPCA_AvgNum"]:
                                #PCA_CalNum = int(arg)
                                gl.set_value("SPCA_AvgNum", int(arg))
                            if opt in ["-n", "--SPCA_n_components"]:
                                gl.set_value("SPCA_n_components", int(arg))

                    if com[0] == "KPCA":
                        #KPCA_Flag = True
                        gl.set_value("KPCA_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:k:",
                                                   ["KPCA_CalNum=", "KPCA_n_components=", "KPCA_AvgNum=", "KPCA_Kernel="])
                        for opt, arg in opts:
                            if opt in ["-c", "--KPCA_CalNum"]:
                                #KPCA_CalNum = int(arg)
                                gl.set_value("KPCA_CalNum", int(arg))
                            if opt in ["-a", "--KPCA_AvgNum"]:
                                #KPCA_CalNum = int(arg)
                                gl.set_value("KPCA_AvgNum", int(arg))
                            if opt in ["-n", "--KPCA_n_components"]:
                                gl.set_value("KPCA_n_components", int(arg))
                            #if opt in ["-g", "--KPCA_Gamma"]:
                                #gl.set_value("KPCA_Gamma", float(arg))
                            if opt in ["-k", "--KPCA_Kernel"]:
                                gl.set_value("KPCA_Kernel", str(arg))

                    if com[0] == "ICA":
                        #ICA_Flag = True
                        gl.set_value("ICA_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:", ["ICA_CalNum=","ICA_n_components=","ICA_AvgNum="])
                        for opt, arg in opts:
                            if opt in ["-c", "--ICA_CalNum"]:
                                #ICA_CalNum = int(arg)
                                gl.set_value("ICA_CalNum", int(arg))
                            if opt in ["-a", "--ICA_AvgNum"]:
                                #ICA_CalNum = int(arg)
                                gl.set_value("ICA_AvgNum", int(arg))
                            if opt in ["-n", "--ICA_n_components"]:
                                    gl.set_value("ICA_n_components", int(arg))

                    if com[0] == "FA":
                        #FA_Flag = True
                        gl.set_value("FA_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:", ["FA_CalNum=","FA_n_components=","FA_AvgNum="])
                        for opt, arg in opts:
                            if opt in ["-c", "--FA_CalNum"]:
                                #FA_CalNum = int(arg)
                                gl.set_value("FA_CalNum", int(arg))
                            if opt in ["-a", "--FA_AvgNum"]:
                                #FA_CalNum = int(arg)
                                gl.set_value("FA_AvgNum", int(arg))
                            if opt in ["-n", "--FA_n_components"]:
                                    gl.set_value("FA_n_components", int(arg))

                    if com[0] == "LDA":
                        #LDA_Flag = True
                        gl.set_value("LDA_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:", ["LDA_CalNum=","LDA_n_components=","LDA_AvgNum="])
                        for opt, arg in opts:
                            if opt in ["-c", "--LDA_CalNum"]:
                                #LDA_CalNum = int(arg)
                                gl.set_value("LDA_CalNum", int(arg))
                            if opt in ["-a", "--LDA_AvgNum"]:
                                #LDA_CalNum = int(arg)
                                gl.set_value("LDA_AvgNum", int(arg))
                            if opt in ["-n", "--n_components"]:
                                if arg != "max":
                                    gl.set_value("LDA_n_components", int(arg))
                                else:
                                    gl.set_value ("LDA_n_components", "max")

                    if com[0] == "TSVD":
                        #TSVD_Flag = True
                        gl.set_value("TSVD_Flag", True)
                        opts, args = getopt.getopt(com[1:], "c:n:a:", ["TSVD_CalNum=","TSVD_Set=","TSVD_AvgNum="])
                        for opt, arg in opts:
                            if opt in ["-c", "--TSVD_CalNum"]:
                                #TSVD_CalNum = int(arg)
                                gl.set_value("TSVD_CalNum", int(arg))
                            if opt in ["-a", "--TSVD_AvgNum"]:
                                #TSVD_CalNum = int(arg)
                                gl.set_value("TSVD_AvgNum", int(arg))
                            if opt in ["-n", "--TSVD_n_components"]:
                                    gl.set_value("TSVD_n_components", int(arg))

                    if com[0] == "Input":
                        opts, args = getopt.getopt(com[1:], "o:i:e:", ["objName=","identFile=","exprFile="])
                        for opt,arg in opts:
                            if opt in ["-o","--objName"]:
                                #objName = arg
                                gl.set_value("objName", arg)
                            if opt in ["-i","--identFile"]:
                                #identFile = arg
                                gl.set_value ("identFile", arg)
                            if opt in ["-e","--exprFile"]:
                                #exprFile = arg
                                gl.set_value ("exprFile", arg)

                    if com[0] == "Output":
                        opts, args = getopt.getopt(com[1:], "o:", ["outputFile="])
                        for opt, arg in opts:
                            if opt in ["-o", "--outputFile"]:
                                #outputFile = arg
                                gl.set_value("outputFile", arg)
