robust_cpp
==========

The functions in `robustmethods_CPP.R` are adapted from a subset of robust statistical methods written by Dr. Rand Wilcox ([R-Forge link](https://r-forge.r-project.org/projects/wrs/); a full collection of his functions can also be found in the file `Rallfun-v22.R` in this github). Different from their originals, the functions in `robustmethod_CPP.R` are partially implemented in `C++`. Due to the iterative nature of these methods, the original R functions can be relatively slow especially when they are used in bootstrapping. The partial C++ implememtations provide substantial speed increases.

The raw C++ code is in `robustmethods_CPP.cpp`; `robustmethods_CPP.so` is the compiled shared object, for 64-bit R on Mac. To use them, follow the steps below:

* (1). Install the R package `RcppArmadillo` (`Rcpp` will be automatically installed too).
* (2). Download and `source()` `Rallfun-v22.R` and `robustmethods_CPP.R`. The first file is needed because of the ancillary functions it stores that are required for function executions.
* (3). Download and `dyn.load()` `robustmethods_CPP.so`
* (4). Run the functions in `robustmethods_CPP.R`. These functions will call the C++ code that has been dynamically loaded.


A list of the R originals and their C++ counterparts as well as their benchmark results are below:


###Benchmark

(Not all functions' benchmark results are shown here. Complete results will be added soon).

    #install.packages("RcppArmadillo")   #Install the package, if it hasn't been installed. 
                                         #Rcpp will automatically be installed too.
    
    require("RcppArmadillo")
    dyn.load(file.choose())    #choose robustmethods_CPP.so
    source(file.choose())      #choose Rallfun-v22.R
    source(file.choose())      #choose robustmethods_CPP.R
    
We will use the function `regci()` from `Rallfun-v22.R` to test the performances of the original functions and their C++ counterparts. `regci()` computes a .95 confidence interval for each of the parameters of a linear regression equation using the percentile bootstrap method. The user can specify the regression estimator to be used (e.g., Theil-Sen (default), least squares, etc) as well as the number of bootstrap samples (nboot=599 by default).

    #Prepare data:
    set.seed(1)
    dv = rnorm(10)
    set.seed(2)
    iv = matrix(rnorm(30), ncol=3)
  
    require("rbenchmark")   #Load the package `rbenchmark

####1. `stsreg` and `stsreg_C`
`stsreg`  computes Theil-Sen regression estimator. It uses Gauss-Seidel algorithm when there is more than one predictor. We will call this function and its C++ counterpart using `regci()`.


    benchmark(replications=1, 
              regci(iv, dv, regfun=stsreg),
              regci(iv, dv, regfun=stsreg_C)
              )
                                  test replications elapsed relative user.self sys.self user.child
    2 regci(iv, dv, regfun = stsreg_C)            1   1.919    1.000     1.918    0.009          0
    1   regci(iv, dv, regfun = stsreg)            1 176.860   92.163   173.182    1.125          0
  




####2. `tshdreg` and `tshdreg_C`
`tshdreg` compute Theil-Sen regression estimator. It uses back-fitting when there is more than one predictor and estimates intercept using Harrel-Davis estimator.


    benchmark(replications=1, 
              regci(iv, dv, regfun=tshdreg),
              regci(iv, dv, regfun=tshdreg_C)
              )

                                   test replications elapsed relative user.self sys.self user.child sys.child
    2 regci(iv, dv, regfun = tshdreg_C)            1   2.185    1.000     2.176    0.017          0         0
    1   regci(iv, dv, regfun = tshdreg)            1  11.885    5.439    11.781    0.110          0         0



####3. `fdepthv2` and `fdepthv2_C`
`fdepthv2` determines the depth of points in `pts` relative to points in `m`. It draws a line between each pair of distinct points and determines depth of the projected points. The final depth of a point is its minimum depth among all projections.

This function allows data to have a singular covariance matrix and it provides a more accurate approximation of halfspace depth. 

`plotit=TRUE` creates a scatterplot when working with bivariate data and pts=NA.

    set.seed(1)
    m = matrix(rnorm(200), ncol=2)
    
    benchmark(replications=1, 
              fdepthv2(m),
              fdepthv2_C(m)
              )

               test replications elapsed relative user.self sys.self user.child sys.child
    2 fdepthv2_C(m)            1   2.012    1.000     1.799    0.131          0         0
    1   fdepthv2(m)            1  14.090    7.003    13.147    0.960          0         0


![plot](http://img844.imageshack.us/img844/7411/spp.png)

