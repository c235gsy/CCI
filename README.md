# CCI

We developed a new approach, counterpart composite index (CCI), to search healthy counterpart of each leukemia cell subpopulation by integrating multiple statistics to project leukemia cells onto healthy hematopoietic cells. CCI combines multiple measures, i.e. expression difference, networks, embedding space in dimension reduction, GSEA et al., to improve statistical power. Facilitated by the comprehensive cell atlas of BMMCs and HSPCs, CCI project leukemia subclones to normal cells to characterize the cell types or progressive status of leukemias.

Run this line of code to start the calculation :
```
python CCI_main.py comment.txt
```

The file comment.txt is something like :
```
Input -o ALL
Input -e Path_to_expression_file
Input -i Path_to_identity_file
Output -o Path_to_write_the_results
Corr -c 10000 -a 5 -m 0.05
Diff -c 10000 -a 5 -m 0.05
PCA -c 10000 -a 1 -n 30
SPCA -c 10000 -a 1 -n 30
KPCA -c 10000 -a 1 -n 30 -k cosine
ICA -c 10000 -a 1 -n 30
FA -c 10000 -a 1 -n 30
LDA -c 10000 -a 1 -n max
TSVD -c 10000 -a 1 -n 30
# -c Number of calculations
# -a Number of cells to get an average cell
# -n Number of data dimensions left after dimension reduction
# -m Limit of the rate of genes where both tow average cells have meaningful expression data


```
