# Equalden.HD
Testing the Equality of a High Dimensional Set of Densities

[SiDOR Group](http://sidor.uvigo.es/en/). [University of Vigo.](http://uvigo.gal/)

## Install the development version from GitHub
```r
devtools::install_github('sidoruvigo/Equalden.HD')
```
## Description (package)
This package implements three different k-sample tests in the setting of low sample size and large dimension. The first one, proposed in Zhan and Hart (2012), test the null hypothesis that all the k-samples come from a single distribution when the number of samples is large, the sample size is small and the samples are independent of each other. The other tests proposed in Cousido-Rocha et. al (2018) are an adaptation of the test of Zhan and Hart (2012) to the setting of dependent sample sizes. The type of dependence considered in such tests is absolutely regular and strongly mixing (see Doukhan, 1995). The standarized version of the test statistic and the p-values are computed among other things.

## Details
+ Package: Equalden.HD
+ Version: 1.0
+ Maintainer: Marta Cousido Rocha martacousido@uvigo.es
+ License: GPL-2

## Equalden.test.HD

### Description (Equalden.test.HD)
The function includes the k-sample test proposed in Zhan and Hart (2012), “indep” and its adaptations for dependent data proposed in Cousido-Rocha et.al (2018), “dep.boot” and “dep.spect”. The k-sample test of Zhan and Hart (2012) test the null hypothesis that all the k-samples come from a single distribution when the number of samples is large, the sample size is small and the samples are independent of each other. The statistic of Zhan and Hart (2012) is based on a comparison of k sample-specific kernel density estimates with a kernel density estimate computed from the pooled sample. An alternative expression of this statistic shows that it can be interpreted as a difference between the intra-samples variability and the inter-samples variability. This statistic is standarized using a variance estimator suitable when the k samples are independent of each other. The asymptotic normality (when k tends to infinity) of the standardized version of the statistic is used to compute the corresponding p-value. Cousido-Rocha et. al (2018) proposed two adaptions of the test of Zhan and Hart (2012) for the setting of dependent samples. These tests consider the statistic proposed in Zhan and Hart (2012) but standarize it using variance estimators suitables when the samples are weak dependent (mixing conditions see Doukhan, 1995). One of tests, “dep.boot”, standarize the statistic using a variance estimator based on the dependent multiplier bootstrap (B ̈uhlmann, 1993, Section 3.3), whereas the other test, “dep.spect”, uses a variance estimator based on the spectral analysis theory. Both tests performs similar, however the “dep.spect” test is computationally more efficient than the “dep.boot” test. Cousido-Rocha et. al (2018) concluded based on their simulation study that for independent samples the tests “dep.boot” and “dep.spect” are equal or more powerful than the test of Zhan and Hart (2012) althought they are protected against possible dependency.

### Parameters
+ **X**: A matrix where each row is one of the k-samples.
+ **method**: the k-sample test. See details.

### Return
 A list with class "htest" containing the following components:
+ **standarized statistic**: the value of the standarized statistic.
+ **p.value**: the p-value for the test.
+ **statistic**: the value of the statistic.
+ **variance**: the value of the variance estimator.
+ **m**: number of significant lags for the variance estimator if the method is “dep.spect” or “dep.boot”. Null if the method is “indep” since no measure of the dependence between the samples is considered in this method.
+ **method**: a character string indicating what k-test was performed.
+ **data.name**: a character string giving the name of the data.


### Usage
Equalden.HD(X, method = c("indep", "dep.boot", "dep.spect"))

### Example
```r
library("Equalden.HD")
n <- 2
k <- 100
X <- matrix(rnorm(n * k), ncol = 2)
res <- Equalden.test.HD(X,  method = "indep")
```

## References
+ Cousido-Rocha, M., de Uña-Álvarez, J. (2018). Testing equality of a large number of densities under mixing conditions. Test (preprinted)
+ Doukhan, P. (1995) Mixing: Properties and Examples. Springer-Verlag, New York.
+ Zhan, D., Hart, J. (2012) Testing equality of a large number of densities. Biometrika, 99, 1-17.


## Authors
+ Cousido-Rocha, Marta.
+ de Uña-Álvarez, Jacobo.
+ Soage González, José Carlos.
