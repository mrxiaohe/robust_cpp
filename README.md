robust_cpp
==========

The functions in `robustmethods_CPP.R` are adapted from a subset of robust statistical methods written by Dr. Rand Wilcox ([R-Forge link](https://r-forge.r-project.org/projects/wrs/); a full collection of his functions can also be found in the file `Rallfun-v22.R` in this github). Different from their originals, the functions in `robustmethod_CPP.R` are partially implemented in `C++`. Due to the iterative nature of these methods, the original R functions can be relatively slow especially when they are used in bootstrapping. The partial C++ implememtations provide substantial speed increases.

##A Mac binary named `WRScpp` is now available to be installed following the instructions via the link below:

##https://github.com/mrxiaohe/WRScpp
