# RobustDC
This R package implements a Huber-type robust distance correlation measure (called "Robust DC" or "RDC") for sure screening in ultrahigh-dimensional heavy-tailed data.

# Installation
To install the `RobustDC` package, you will first need to install `devtools` package and then execute the following code: 
```
devtools::install_github('kehongjie/RobustDC')
```

# Main Functions
Two main functions in this package are `RDC` and `RDC.lasso`. `RDC` is for screening with robust distance correlation, and `RDC.lasso` is for 
fitting the final LASSO regularized regression model after RDC. You can always use the following command to see more details:
```
library(RobustDC)
?RDC
?RDC.lasso
```
