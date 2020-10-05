# betafunctions v. 1.3.1

- Added `true.model` argument allowing for greater control over the true-score estimation procedure in `Beta.tp.fit()`, `LL.CA`, and `LL.ROC`.

	- The `true.model` allows for specifying whether to fit four- or two-parameter Beta distribution to the estimated moments of the true-score distribution.

	- Currently, the options for the argument are "4P" and "2P". Further options might be added in the future.

- Added the possibility of specifying a Beta error model for the `LL.ROC()` function by way of the `true.model` argument.

- Additional correction of typographical errors in documentation.

---

# betafunctions v. 1.3.0

- Fixed some typographical errors in documentation.

- Minor modification to `LL.CA()` where the function now won't automatically terminate if true-score distribution fitting procedure produces `NA` or `NaN` estimates.

	- True-score distribution fitting now performed using the new `Beta.tp.fit()` function.

	- The `grainsize` argument is now fully removed from `LL.CA()`.

	- Added additional checks and diagnostics in `LL.CA()`, issuing warnings for aberrant events.

	- If a list of parameter values are supplied in place of a score-vector, it is now necessary to supply an `etl` (effective test length) parameter as well. See documentation for the `ETL()` function for more information.

- Fixed an error in the fitting procedure of `Beta.4P.fit()` which occurred during positive skew.

- Added the `Beta.tp.fit()` function for estimating four-parameter Beta distribution parameters for an underlying true-score distribution, assuming that the observations are generated from a Beta-Binomial model.

---

# betafunctions v. 1.2.1

- Fixed the `LL.ROC()` function which stopped working after the previous update.

	- Included possibility of specifying the number of points for which the ROC curve is drawn and the AUC statistic is calculated by way of the new `grainsize`argument.

---

# betafunctions v. 1.2.0

- Added functionality for calculating consistency statistics by way of the new `ccStats()` function.

- Substantial developments concerning primarily the `LL.CA()` function, adding functionality and improving performance.

	- Added calculation of consistency statistics by calling the new `ccStats()` function.

	- Added `output` argument indicating which statistics to calculate. Default is to compute both accuracy and consistency.

	- Calculation of distribution-based output now utilizes the `integrate()` function.

		- Subsequently rendered the `grainsize` argument inert, as it was used for the previous method.

	- Added possibility of supplying list of custom parameter values for the four-parameter Beta true-score distribution, forgoing the need to estimate one.

	- Added `override` argument providing the possibility of overriding the default fail-safe reverting to a two-parameter Beta true-score distribution, should the fitting procedure produce impermissible parameter estimates.

---

# betafunctions v. 1.1.0

- Fixed small typographical errors in documentation.

- Added the `cba()` function for computing Cronbach's alpha.