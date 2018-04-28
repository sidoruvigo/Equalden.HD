# Equalden.HD
Testing the Equality of a High Dimensional Set of Densities

[SiDOR Group](http://sidor.uvigo.es/en/). [University of Vigo.](http://uvigo.gal/)

## Install the development version from GitHub
```r
devtools::install_github('sidoruvigo/Equalden.HD')
```
## Description (package)
This package implements three different methods to test the null hypothesis that a large number k of samples have a common density. The sample size can be as small as 2. These methods are particularly well suited to the low sample size, high dimensional setting (n << k). The first method, proposed by Zhan and Hart (2012), was developed to test the null hypothesis when the samples are independent of each other. The other tests, proposed by Cousido-Rocha et al. (2018), are adaptations of the test in Zhan and Hart (2012) for the setting in which the samples are weakly dependent. The standarized version of each test statistic and its p-value are computed among other things.

## Details (package)
+ Package: Equalden.HD
+ Version: 1.0
+ Maintainer: Marta Cousido Rocha martacousido@uvigo.es
+ License: GPL-2

## Equalden.test.HD

### Description (Equalden.test.HD)

Performs the k-sample test proposed by Zhan and Hart (2012) for the low sample size, high dimensional setting with independent samples, and its extensions for dependent samples proposed by Cousido-Rocha et al. (2018).

### Parameters
+ **X**: A matrix where each row is one of the k-samples.
+ **method**: the k-sample test. By default the “dep.spect” method is computed. See details.

### Details (Equalden.test.HD)
The function implements the k-sample test proposed by Zhan and Hart (2012), method=“indep”, and its extensions for dependent data proposed by Cousido-Rocha et al. (2018), method=“dep.boot” and “dep.spect”. The method proposed by Zhan and Hart (2012) serves to test the null hypothesis that the k-samples have a common distribution. It is suitable when the k samples are independent and the number of samples k is large, and it works for sample sizes as small as 2. The statistic in Zhan and Hart (2012) is based on a comparison between the k sample-specific kernel density estimates and the kernel density estimate computed from the pooled sample. An alternative expression of this statistic shows that it can be interpreted as a difference between the intra-samples variability and the inter samples variability. This statistic is standarized using a variance estimator which is valid for independent samples. The asymptotic normality (when k tends to infinity) of the standardized version of the statistic is used to compute the corresponding p-value. Cousido-Rocha et al. (2018) proposed two corrections of the test of Zhan and Hart (2012) for dependent samples. These tests standarize the statistic proposed in Zhan and Hart (2012) by using variance estimators which are suitable when the samples are weakly dependent. The method “dep.boot” implements the dependent multiplier bootstrap to estimate the variance, whereas the method “dep.spect” uses a variance estimator based on the spectral analysis theory. Both tests perform similarly, but the “dep.spect” test tends to be computationally more efficient than the “dep.boot” test. Cousido-Rocha et al. (2018) showed through simulations that, for independent samples, the tests “dep.boot” and “dep.spect” may be more powerful than the test in Zhan and Hart (2012) despite of being protected against possible dependences.

### Return
 A list containing the following components:
+ **standarized statistic**: the value of the standarized statistic.
+ **p.value**: the p-value for the test.
+ **statistic**: the value of the statistic.
+ **variance**: the value of the variance estimator.
+ **m**: number of significant lags for the variance estimator if the method is “dep.spect” or “dep.boot”. Null if the method is “indep” since no measure of the dependence between the samples is considered in this method.
+ **k**: number of samples or populations.
+ **n**: sample size.
+ **method**: a character string indicating what k-test was performed.
+ **data.name**: a character string giving the name of the data.


### Usage
```r
Equalden.test.HD(X, method = c("indep", "dep.boot", "dep.spect"))
```

### Example
```r
library("Equalden.HD")
n <- 2
k <- 100
set.seed(1234)
X <- matrix(rnorm(n * k), ncol = 2)
res <- Equalden.test.HD(X,  method = "indep")

### The statistic and the variance estimator
res$statistic
res$variance
### The number of samples and sample size
res$k
res$n

### Real data analysis. We test the null hypothesis that 1000 randomly selected genes
### measured in patients with BRCA2 mutations have a common distribution. We use the test
### proposed in Cousido-Rocha et al. (2018) since correlation among expression levels of
### different genes on the same individual is expected.
data(Hedenfalk)
X <- Hedenfalk
k <- dim(X)[1]
### We eliminate the additive patients effects by substracting to each column its sample mean.
BRCA2 <- sweep(X[, 8:15], 2, apply(X[, 8:15], 2, mean))
set.seed (1234)
k <- sample(1:k, 1000)
res1 <- Equalden.test.HD(BRCA2[k, ], method = "dep.boot")
res1
res2 <- Equalden.test.HD(BRCA2[k, ], method = "dep.spect")
res2
```

## References
+ Cousido-Rocha, M., de Uña-Álvarez, J., and Hart, J.(2018). Testing equality of a large number of densities under mixing conditions. Preprint.
+ Zhan, D., Hart, J. (2012) Testing equality of a large number of densities. Biometrika, 99, 1-17.


## Authors
+ Marta Cousido-Rocha.
+ José Carlos Soage González.
+ Jacobo de Uña-Álvarez.
+ Jeffrey D. Hart.
