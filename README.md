# ClustOfVar
This R package is dedicated to the clustering  of variables. Variables can be quantitative, qualitative or a mixture of both. It provides hierarchichal and k-means clustering of a set of variables. The center of a cluster of variables is a synthetic variable but is not a ’mean’ as for classical k-means or Ward clustering. This synthetic variable is the first principal component calculated by PCAmix. The homogeneity of a cluster of variables is defined as the sum of the correlation ratio (for qualitative variables) and the squared correlation (for quantitative variables) between the variables and the center of the cluster, which is in all cases a numerical variable. This package deals with datasets of thousands of variables like gene expression data and can also be used for dimension reduction purpose.

## Install

To install the current development version from github, use :

```{r eval=FALSE}
devtools::install_github("chavent/ClustOfVar")
# This needs the devtools package to be installed :
# install.packages("devtools")
```
