# PieBW-GL
MP-lasso chart: Multi-level polar chart for visualizing group lasso analysis of high dimensional data

# Notice
All source codes were listed in file "PieBW_gglasso.R", "PieBW_SGL.R" for implemeting

# Citation
- Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for Computing Group-Lasso Penalized Learning Problems,'' Statistics and Computing. 25(6), 1129-1141.
- Simon, N., Friedman, J., Hastie, T., and Tibshirani, R. (2011) A Sparse-Group Lasso
- Xiaoxuan Liang, Aaron Cohen, Anibal Solón Heinsfeld, Franco Pestilli, Daniel J. McDonald (2022) sparsegl: An R Package for Estimating Sparse Group Lasso

# Example
Try run_example.R to get a quick start

example.csv contains data used for simulation study in the paper and generated as follow
```
set.seed(1010)
n <- 100
p <- 200
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta <- c(
  rep(3, 4), #absolute_group_mean=3, absolute_max_beta=3
  c(8, -1, 1, 0, 0), #absolute_group_mean=2, absolute_max_beta=8
  rep(-1, 6), 
  c(1, -3, 3, 0, 0),
  rep(0, (p - 20)) #unsignificant group
)
groups <- c(rep(1,4), rep(2,5), rep(3,6), rep(4,5) ,rep(5:(p / 5), each = 5))
groups_name <- paste0("Group : ", groups)

eps <- rnorm(n, mean = 0, sd = 1)
y <- X %*% beta + eps # continuous response variable
pr <- 1 / (1 + exp(-X %*% beta))
y0 <- rbinom(n, 1, pr) # binary response variable
```
# Usage
In R:
- Load package for use
```
library(gglasso)
library(ggplot2)
library(dplyr)
library(forcats)
library(gridExtra)
library(SGL)
library(ggiraph)
library(patchwork)
```
- Import example data
```
example <- read.csv("example.csv")
X <- example$X
y <- example$y
y0 <- example$y0
groups <- example$groups
groups_name <- example$groups_name
```
- Implement cross-validation
```
cv_gl <- cv.gglasso(X, y, group=groups, loss="ls", nfolds=3)
cv_gl2 <- cv.gglasso(X, y0, group=groups, loss="logit", nfolds=3)
cv_sgl <- cvSGL(data=list(x=X,y=y), index=groups, type="linear", nfold=3, alpha=0.5, standardize=FALSE)
cv_sgl2 <- cvSGL(data=list(x=X,y=y0), index=groups, type="logit", nfold=3, alpha=0.5, standardize=FALSE)
```
- Visualizing group lasso analysis
```
source("PieBW_gglasso")
source("PieBW_SGL")

PieBW_gglasso(cv_object = cv_gl, group = groups_name, lambda.type = "min", sort.type = "max")
PieBW_gglasso(cv_object = cv_gl, group = groups_name, lambda.type = "min", sort.type = "mean")
PieBW_gglasso(cv_object = cv_gl2, group = groups_name, lambda.type = "1se", sort.type = "max")
PieBW_gglasso(cv_object = cv_gl2, group = groups_name, lambda.type = "1se", sort.type = "mean")

PieBW_SGL(cv_object = cv_sgl, group = groups_name, lambda.type = "min", sort.type = "max")
PieBW_SGL(cv_object = cv_sgl, group = groups_name, lambda.type = "min", sort.type = "mean")
PieBW_SGL(cv_object = cv_sgl2, group = groups_name, lambda.type = "1se", sort.type = "max")
PieBW_SGL(cv_object = cv_sgl2, group = groups_name, lambda.type = "1se", sort.type = "mean")
```
- Input
  + cv_object : fitted cv.gglasso() object or cvSGL() object
  + group : group name, character or integer.
  + lambda.type : selection of optimal lambda value. if lambda.type=”min”, then select 
                 lambda vaule which gives minimum loss for model. if lambda.type=”1se”, then select
                 lambda vaule which gives minimum loss in one stadard error range.
  + sort.type : selection of sort type for group importance. if sort.type=”min”, then group
               sorted by absolute average of coefficients in each group. if sort.type=”max”, then
               group sorted by absolute average of coefficients in each group.
