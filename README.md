# CCI

We developed a new approach, counterpart composite index (CCI), to search healthy counterpart of each leukemia cell subpopulation by integrating multiple statistics to project leukemia cells onto healthy hematopoietic cells. CCI combines multiple measures, i.e. expression difference, networks, embedding space in dimension reduction, GSEA et al., to improve statistical power. Facilitated by the comprehensive cell atlas of BMMCs and HSPCs, CCI project leukemia subclones to normal cells to characterize the cell types or progressive status of leukemias.

Run this line of code to start the calculation :
```
python CCI_main.py comment.txt
```
The **test_data** could be used directly.

The file comment.txt is something like :
```
Input -o ALL
Input -e Path_to_expression_file
Input -i Path_to_identity_file
Output -o Path_to_write_the_results
Corr -c 10000 -a 5 -m 0.05
Diff -c 10000 -a 5 -m 0.05
PCA -c 10000 -a 1 -n 30
KPCA -c 10000 -a 1 -n 30 -k cosine
ICA -c 10000 -a 1 -n 30
FA -c 10000 -a 1 -n 30
LDA -c 10000 -a 1 -n max
TSVD -c 10000 -a 1 -n 30
```
The output contains three kinds of files : 

*Mean file* : <br>..._**Method**_mean.txt : The _median matrix_ of multiple calculations 

*BF file* : <br>..._**Method**_BF.txt : The _BF matrix_ of multiple calculations 

*Information file* : <br>..._**infor**.txt : The information of calculations and data 

|  comment | meaning  |
|:---:|---|
| Input -o |  **Not useful at the moment** , but need to equal to **ALL**  |
| Input -e | Path to the **expression file** |
| Input -i  | Path to the **identity file** |
| Output -o | Path to write the results |
| Corr | Calculate the **Correlation Coefficient** of all the expression data of two cellular genes |
| Diff  | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes |
| PCA | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes after **PCA** dimensionality reduction |
| KPCA | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes after **Kernel PCA** dimensionality reduction |
| ICA | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes after **ICA** dimensionality reduction |
| FA | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes after **Factor Analysis** dimensionality reduction |
| LDA | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes after **Linear Discriminant Analysis** dimensionality reduction |
| TSVD | Calculate the average **Euclidean Distance** of all the expression data of two cellular genes after **Truncated SVD** dimensionality reduction |
| -c  | Number of calculations |
| -a | Number of cells to get an average cell |
| -n | Number of data dimensions left after dimension reduction. If **-n max** in LDA model, it means that **n = min(num_features, mun_groups)-1** |
| -m | Limit of the rate of genes where both tow average cells have meaningful expression value|
| -k | The **kernel function** used in **Kernel PCA** |
| # | **Make a line of commands invalid** |




