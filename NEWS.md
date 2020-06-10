betafunctions v. 1.2.0
- Added functionality for calculating consistency statistics by way of the new ccStats() function.
- Substantial developments concerning primarily the LL.CA() function, adding functionality and improving performance.
-- Added calculation of consistency statistics.
-- Added "output" argument indicating which statistics to calculate. Default is to compute both accuracy and consistency.
-- Calculation of distribution-based output now utilizes the integrate() function.
--- Subsequently removed the "grainsize" argument, as it was used for the previous method.
-- Added possibility of supplying list of custom parameter values for the four-parameter Beta true-score distribution, forgoing the need to estimate one.
-- Added "override" argument providing the possibility of overriding the default fail-safe reverting to a two-parameter Beta true-score distribution, should the fitting procedure produce impermissible parameter estimates.



betafunctions v. 1.1.0
- Fixed small typographical errors in documentation.
- Added function for computing Cronbach's alpha.