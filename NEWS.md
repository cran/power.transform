# Version 1.0.1

## Changes

- The empirical central normality test is now dependent on the number of samples.
  The test can now be called using `ecn.test`. `cn.test` is a stricter central
  normality test whose test statistics are determined from strictly normal
  distributions, instead of normal distributions with up to 10% outliers.
  
- The robust location- and shift-invariant transformations now use weights
  optimised for achieving central normality.


# Version 1.0.0

This is the initial public release of the `power.transform` package.
