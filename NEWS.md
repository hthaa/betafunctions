# betafunctions v. 1.4.1

- The `true.model` argument of the `Beta.tp.fit()` function now includes a `"3P"` option, allowing for the specification of one location-parameter (l or u) and estimating the remaining location-parameter and the shape-parameters (alpha and beta) so as to make the resulting distribution have the same skewness and kurtosis as the estimated true-score ditribution.

- The `AMS()` and `BMS()` functions now issue warnings if there was not enough information supplied to calculate the target parameter.

- The `ETL()` function: For the sake of argument consistency, the `l` and `u` arguments are not renamed `min` and `max`, respectively.

# betafunctions v. 1.4.0

- Added a number of new functions for working with Beta distributions.

	- The new `LABMSU()` function allows for finding the lower-bound parameter for the four-parameter Beta distribution by supplying the shape-parameters and moments of the resulting distribution, and (optionally) the upper-bound location parameter.
	
	- The new `UABMSU()` function allows for finding the upper-bound parameter for the four-parameter Beta distribution by supplying the shape-parameters and moments of the resulting distribution, and (optionally) the lower-bound location parameter.

- Added additional functionality to some existing functions, allowing for specifying the lower- and upper-bounds, and/or moments of the resulting distribution.

	- The `AMS()` function now includes `l` and `u` arguments, finding the Alpha shape parameter necessary to produce a Beta distribution with target moments and specified lower and upper bounds of the resulting distribution.

	- The `BMS()` function now includes `l` and `u` arguments, finding the Beta shape parameter necessary to produce a Beta distribution with target moments and specified lower and upper bounds of the resulting distribution.

	- The `Beta.2p.fit()` function, essentially a wrapper-function for `AMS()` and `BMS()`, now includes `l` and `u` arguments for specifying the bounds for the resulting distribution.

- Functions relating to estimating classification accuracy and consistency:
	
	- `Beta.tp.fit()`, `LL.CA()`, and `LL.ROC()`: 
	
		- Added `true.model` argument allowing for greater control over the true-score estimation procedure in.

		- The `true.model` allows for specifying whether to fit four- or two-parameter Beta distribution to the estimated moments of the true-score distribution.

		- Currently, the options for the argument are "4P" and "2P". Further options might be added in the future.

	- The `Beta.tp.fit()` function: 

		- Now allows for specifying the lower- and/or upper bounds of the two-parameter solutions.

		- Now allows for specifying that the fitted distribution should estimate two parameters rather than the default of four. Estimates then the parameters necessary to produce a distribution with the same mean and variance as the estimated true-score distribution, given specified moments, shape- and location parameters (default is `l = 0` and `u = 1`).

		- Added the option of specifying the true-score distribution moments (sans a functional form) as output rather than the estimated parameters of the Beta true-score distribution in the new `output` argument. The default is to retrieve the estimated parameters of the Beta true-score distribution.

		- When `parameters` is specified in the `output` argument, the effective test length is included as part of the output. As such, the default output of the `Beta.tp.fit()` function is now complete for the purposes of being used as input for the `LL.CA()` function.

	- The `LL.CA()` function:

		- The `override` and `failsafe` arguments:

			- The `override` argument is rendered inert, to be removed in a later version.

			- The `failsafe` argument is introduced to replace the `override` argument. It is essentially an inversion of the `override` argument.

			- The above change is to streamline and make consistent arguments structures across functions.


- General code cleaning:

	- Changed some argument names so as to have the names be consistent across functions. Argument locations in the function calls remain unchanged.

		- Function arguments previously named `lt` now named `lower.tail`.

		- Function arguments previously named `a` or `b` now named `alpha` and `beta` (respectively).

		- Function arguments previously named `var` now named `variance`.


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