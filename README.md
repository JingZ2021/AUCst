
# AUCst

<!-- badges: start -->
<!-- badges: end -->

Simulated example dataset and sample code for calculating the two-dimensional AUC for one longitudinal biomarker, under the regular visit scenario.

## Installation

You can install the development version of AUCst like so:

``` r
library(devtools)
install_github("JingZ2021/AUCst",force = TRUE)
library(AUCst)
```

## Example

R file contians functions used to implement the proposed method.

``` r
library(AUCst)
## basic example code
da.sim.short <- AUCst::da.sim.short
da.sim.long <- AUCst::da.sim.long
sort(unique(da.sim.long$vtime))  
obj.i=AUCst::main1.sub.func(da=da.sim.short, 
                            da.long=da.sim.long,
                            sk=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2), 
                            tseq.eval=c(0.4,1.11,1.50,1.90), resample =0, nsap=3)
obj.i$AUC          
#>             [,1]      [,2]      [,3]      [,4]
#> AUC.st 0.8298518 0.8108635 0.7480704 0.6767092
#> AUC.st 0.8370427 0.8107604 0.7437303 0.6677119
#> AUC.st 0.8446459 0.8123010 0.7420099 0.6624297
#> AUC.st 0.8522235 0.8150176 0.7423900 0.6603039
#> AUC.st        NA 0.8184398 0.7443130 0.6607144
#> AUC.st        NA 0.8221120 0.7472090 0.6630092
#> AUC.st        NA 0.8256029 0.7505132 0.6665251
#> AUC.st        NA 0.8285070 0.7536762 0.6706026
#> AUC.st        NA 0.8304380 0.7561665 0.6745941
#> AUC.st        NA 0.8310145 0.7574653 0.6778676
#> AUC.st        NA 0.8298386 0.7570525 0.6798030
#> AUC.st        NA        NA 0.7543852 0.6797820
obj.i$theta
#> [1]  0.9220660  0.8021574  2.3696226  0.3791758 -2.2648142 -1.0925795  0.5177396 -0.4832387  0.1538552
#> [10]  0.4035965
obj.i$SE
#>              [,1]        [,2]       [,3]       [,4]
#>  [1,] 0.008601802 0.004602132 0.01445407 0.03311861
#>  [2,] 0.003474240 0.006420346 0.01939001 0.03336046
#>  [3,] 0.006042407 0.006324712 0.02171633 0.03542930
#>  [4,] 0.012034498 0.004582757 0.02146406 0.03526624
#>  [5,]          NA 0.003079514 0.02016032 0.03292638
#>  [6,]          NA 0.004550213 0.01949487 0.02962096
#>  [7,]          NA 0.006941513 0.02048720 0.02682730
#>  [8,]          NA 0.008144392 0.02271480 0.02557692
#>  [9,]          NA 0.007029369 0.02493444 0.02564940
#> [10,]          NA 0.002787379 0.02633627 0.02581909
#> [11,]          NA 0.005716703 0.02768279 0.02539897
#> [12,]          NA          NA 0.03236712 0.02654454
obj.i$conv
#> [1] 0
```

## Reference

Zhag, J., Ning, J., Huang X, & Li R. (2021). On the Time-varying Predictive Performance of Longitudinal Biomarkers: Measure and Estimation. Statistics in Medicine,
40(23):5065-5077. <https://doi.org/10.1002/sim.9111>
