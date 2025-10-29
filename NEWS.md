# Version 1.0.2

## Changes

- `ecn.test` was renamed to `ecn_test`. It now additionally takes `tau`, `n`,
  and `kappa` as arguments (see function documentation).

- `cn.test` was deprecated. Use `ecn_test` with `kappa = 1.0` instead.

- When checking if required packages are installed, `power.transform` now caches
  results. This prevents unnecessary look-up.

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
