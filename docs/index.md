# Hospital-Level Empirical Bayes Estimation

hospital_ebayes estimates hospital-level empirical Bayes models using value-added methods adapted from education research. The command implements several key modifications for healthcare settings:

1. Controls for hospital volume effects
2. Allows for the addition of shrinkage target 
3. Provides additional leave-out estimators and intermediate outputs 
4. Includes hospital-specific adjustments

## Syntax

hospital_ebayes depvar, hospitalid(varname) year(varname) [options]

### Required Options
- hospitalid(varname): hospital identifier
- year(varname): year identifier

### Optional Options
- controls(varlist): additional control variables
- shrinkage_target(varlist): variables to control for before shrinkage
- absorb(varname): fixed effects to absorb
- tfx_resid: hospital fixed effects residuals
- by(varname): estimate separately by groups
- data(string): data handling ("preserve", "tv", "merge tv")
- output(string): output file path prefix
- driftlimit(#): maximum number of lags (-1 for all)
- debug: display additional diagnostic information during estimation
- leaveout_years(string): relative year ranges to exclude from estimation
- leaveout_vars(string): names for variables that will store leave-out estimates

## Examples

Basic model with controls:

```stata
hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid)
```

This basic model estimates hospital effects on mortality while controlling for patient characteristics (age, gender, and comorbidities).
Model with Volume Adjustment

*Model with Volume shrinkage target*

```stata
hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) shrinkage_target(log_volume)
```

This model adds a volume shrinkage target by including hospital log volume.
Model with Leave-out Periods

*Model with Leave-out Periods*
```stata
hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) ///
    leaveout_years("-3, 1,3") leaveout_vars(tv_pre tv_post)
```

- For each year t, this model creates:
- tv_pre: estimates excluding data from years t-3 to infitnity
- tv_post: estimates excluding data from years t+1 to t+3 (following 3 years)


## Controls versus shrinkage target 

## Testing >
