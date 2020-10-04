###### `betafunctions` v. 1.3.0

- Fixed some typographical errors in documentation.

- Minor modification to `LL.CA()` where the function now won't automatically terminate if `Beta.4P.fit()` produces `NaN` estimates.

	- If `Beta.4P.fit()` produces `NaN` estimates, the default failsafe of reverting to a two-parameter solution is employed.

	- The `grainsize` argument is now fully removed from `LL.CA()`.

	- Added additional checks and diagnostics in `LL.CA()`, issuing warnings for aberrant events.

- Fixed an error in the fitting procedure of `Beta.4P.fit()` which occurred during positive skew.

- Added the `RF.Beta.4P.fit()` function for estimating four-parameter Beta distribution parameters from data generated from a Beta compound Binomial model.

---

###### `betafunctions` v. 1.2.1

- Fixed the `LL.ROC()` function which stopped working after the previous update.

	- Included possibility of specifying the number of points for which the ROC curve is drawn and the AUC statistic is calculated by way of the new `grainsize`argument.

---

###### `betafunctions` v. 1.2.0

- Added functionality for calculating consistency statistics by way of the new `ccStats()` function.

- Substantial developments concerning primarily the `LL.CA()` function, adding functionality and improving performance.

	- Added calculation of consistency statistics by calling the new `ccStats()` function.

	- Added `output` argument indicating which statistics to calculate. Default is to compute both accuracy and consistency.

	- Calculation of distribution-based output now utilizes the `integrate()` function.

		- Subsequently rendered the `grainsize` argument inert, as it was used for the previous method.

	- Added possibility of supplying list of custom parameter values for the four-parameter Beta true-score distribution, forgoing the need to estimate one.

	- Added `override` argument providing the possibility of overriding the default fail-safe reverting to a two-parameter Beta true-score distribution, should the fitting procedure produce impermissible parameter estimates.

---

###### `betafunctions` v. 1.1.0

- Fixed small typographical errors in documentation.

- Added the `cba()` function for computing Cronbach's alpha.