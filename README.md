# GGEE


GGEE implements the group lasso Gene-Gene Eigen Epistasis (G-GEE) method for detecting epistasis at the gene level. The proposed functions allow to compute interaction variables for each pair and to fit a general model with a Group Lasso penalty. The package allows to generate gene structured genotype data and continuous phenotype.


## Installation

The package can be installed from github:

```{r}
# install.packages("devtools")
devtools::install_github("vstanislas/GGEE")
```

Some functions of GGEE use C++ code. On Linux compilation can be interupt depending on the g++ version. The solution is to create a .R/Makevars file in your home directory that contains:
```{r}
CXXFLAGS= -std=c++11
```




## Examples and descriptions
A detailled example can be find in GGEEvignette.pdf


