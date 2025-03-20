# WIP xsHospital Empirical Bayes

This program estimates hospital level empirical bayes models. It adapts teacher value-added methods for hospital settings with several key modifications:

1. Allows users to control for hospital volume effects
2. Handles hospital-specific structure, with one "classroom" per hospital-year 
3. Provides additional leave-out estimators and intermediate outputs
4. Includes hospital-specific adjustments

```stata
Syntax:
vamhosp depvar, hospitalid(varname) year(varname) [options]

Required Arguments:
- depvar: Dependent variable (outcome measure)
- hospitalid: Hospital identifier 
- year: Year identifier

Optional Arguments:
- class: Ward/unit identifier (defaults to 1 class per hospital-year)
- by: Estimate separately by groups
- controls: Additional control variables
- shrinkage_target: Variables to control for before shrinkage estimation
- absorb: Fixed effects to absorb
- tfx_resid: Hospital fixed effects residuals
- data: Data handling options ("preserve", "tv", "merge tv", etc.)
- output: Output file path prefix
- driftlimit: Maximum number of lags (-1 for all)
- leaveout_years: Year ranges to leave out
- leaveout_vars: Variable mappings
```


## Example:

```
// Load sample hospital data
use hospital_data.dta, clear

// Run basic model with controls
vamhosp mortality, hospitalid(providerid) year(year) ///
    controls(age female comorbid) 

// Run model with volume adjustment
vamhosp mortality, hospitalid(providerid) year(year) ///
    controls(age female comorbid) ///
    shrinkage_target(log_volume)

// Run model with leave-out years
vamhosp mortality, hospitalid(providerid) year(year) ///
    controls(age female comorbid) ///
    leaveout_years("2010,2015") ///
    leaveout_vars(tv_early tv_late)
```


