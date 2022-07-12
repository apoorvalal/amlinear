rm(list = ls())

# %%
fns = list.files("R/", ".R", full.names = T)
sapply(fns, source, .GlobalEnv)

# %%
library(LalRUtils); library(CVXR); library(glmnet)

# %%
data(lalonde.exp)
yn = 're78'; wn = "treat"; Xn = setdiff(colnames(lalonde.exp), c(yn, wn))
Y = lalonde.exp[[yn]]; W = lalonde.exp[[wn]]
X = lalonde.exp[, Xn]

# %%
average_partial_effect(X, Y, W)
# %%
