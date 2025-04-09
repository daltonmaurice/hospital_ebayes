*! version 0.0.1  February 2019 Maurice Dalton, daltonm
/* Based on original code written by Michael Stepner, forked */

/*******************************************************************************
Hospital Value-Added Model with Leave-Out Estimators
--------------------------------------------------------------------------------

This program estimates hospital value-added models using various leave-out 
estimators. It adapts teacher value-added methods for hospital settings with 
several key modifications:

1. Allows users to controls for hospital volume effects
2. Handles hospital-specific structure, note we impose one "classroom" per hospital-year
3. Provides additional leave-out estimators and intermediate outputs
4. Includes hospital-specific adjustments

Required Arguments:
- depvar:       Dependent variable (outcome measure)
- hospitalid:   Hospital identifier
- year:         Year identifier

Optional Arguments:
- class:         Ward/unit identifier - leftover from original code but not used. We set this 
                 to 1 for all observations, so there is one class per hospital-year
- by:          Estimate separately by groups
- controls:    Additional control variables
- shrinkage_target: Variables to control for before shrinkage estimation
- absorb:      Fixed effects to absorb
- tfx_resid:   Hospital fixed effects residuals
- data:        Data handling options ("preserve", "tv", "merge tv", etc.)
- output:      Output file path prefix
- driftlimit:  Maximum number of lags (-1 for all)
- leaveout_years: New parameter for year ranges to leave out
- leaveout_vars: New parameter for variable mappings
- debug:       New debug parameter

Usage Example:
    vamhclose score, hospitalid(hospital) year(year) ///
        controls(age female) shrinkage_target(volume)

    // Example with leave-out estimation
    vamhclose score, hospitalid(hospital) year(year) ///
        controls(age female) shrinkage_target(volume) ///
        leaveout_years("-2,2 -1,1") leaveout_vars("tv_2yr tv_1yr")
        
    /* The leaveout example above will:
    1. Create tv_2yr using data excluding 2 years before/after current year
    2. Create tv_1yr using data excluding 1 year before/after current year
    Format is "before,after" where negative numbers are years before */

Notes:
- Requires Stata 10.2+
- Missing values in key variables are automatically dropped
- Hospital IDs should be consistent across years

Authors:
Maurice Dalton 
Doug Staiger 
---
Based on vam.ado written by Michael Stepner version 2.0.1  27jul2013.
*******************************************************************************/

// Include Mata functions
include "hospital_ebayes_mata.mata"

cap program drop  hospital_ebayes
program define hospital_ebayes
version 10.2

set more off
syntax varname(ts fv), hospitalid(varname) year(varname) [class(varname) ///
    by(varlist) ///
    shrinkage_target(varlist) ///
    controls(varlist ts fv) absorb(varname) tfx_resid(varname) ///
    data(string) output(string) output_addvars(varlist) ///
    driftlimit(integer -1) ///
    leaveout_years(string) /// New parameter for year ranges to leave out
    leaveout_vars(string) /// New parameter for variable mappings
    debug] /// New debug parameter

* By default we use 1 class or ward per hospital. We didnt feel there was 
* a direct comparable unit to classrooms within a hospital.
if "`class'" == "" {
    tempvar class_var 
    egen `class_var'=group(`hospitalid' `year')
    local class `class_var' 
}

* Error checks
local depvar `varlist'

capture confirm variable score_r, exact
if (_rc==0) {
    di as error "The dataset loaded in memory when vam is run cannot have a variable named score_r."
    exit 110
}

capture confirm variable tv, exact
if (_rc==0) {
    di as error "The dataset loaded in memory when vam is run cannot have a variable named tv."
    exit 110
}


if ("`leaveout_years'"!="") {
    // Parse the leaveout rules
    local n_rules = 0
    foreach rule in `leaveout_years' {
        local ++n_rules
        tokenize "`rule'", parse(",")
        local rule_`n_rules'_before "`1'"
        local rule_`n_rules'_after "`3'"
    }
    
    // Parse variable names and validate count matches
    tokenize `leaveout_vars'
    local n_vars = 0
    foreach v in `leaveout_vars' {
        local ++n_vars
    }
    
    if (`n_rules' != `n_vars') {
        di as error "Error: Number of leaveout rules (`n_rules') does not match number of leaveout variables (`n_vars')"
        di as error "leaveout_years: `leaveout_years'"
        di as error "leaveout_vars: `leaveout_vars'"
        exit 198
    }
    
    // Parse variable names
    forvalues i = 1/`n_rules' {
        local var_`i' "``i''"
        capture confirm variable ``i'', exact
        if (_rc==0) {
            di as error "The dataset loaded in memory cannot have a variable named ``i''."
            exit 110
        }
    }
}

local merge_tv=0
local merge_resid=0
if ("`data'"=="") local data="preserve"
else {
    if !inlist("`data'","preserve","tv","merge tv","merge score_r","merge tv score_r","merge score_r tv","variance") {
        di as error "Not a valid argument for data. Choose either 'preserve', 'tv', 'merge [tv AND/OR score_r]', or 'variance'."
        exit 198
    }
    else {
        tokenize "`data'"
        if ("`1'")=="merge" {
            if ("`2'"=="tv") | ("`3'"=="tv") local merge_tv=1
            if ("`2'"=="score_r") | ("`3'"=="score_r") local merge_resid=1
        }
    }
}

if "`tfx_resid'"!="" & "`absorb'"!="" {
    di as error "Cannot specify an absorb variable and a tfx_resid variable simultaneously."
    exit 198
}

* If output was left blank, set a tempfile for the tv output
if `"`output'"'=="" {
    tempfile output
    local nooutput=1
}
else local nooutput=0

* Start log
if (`nooutput'!=1) log using `"`output'_log"', replace name(t) text

* Process by variables
if ("`by'"!="") {
    tempvar byvar
    egen `byvar'=group(`by'), label
    sum `byvar', meanonly
    local by_vals=`r(max)'
}
else local by_vals=1

****************

preserve

*** Run through separately for each by-value.
local firstloop=1
forvalues l=1/`by_vals' {

    if (`firstloop'!=1) restore, preserve

    *** Print heading (with by-variable identifier if applciable)
    di "{txt}{hline}"
    if ("`by'"!="") {
        local bylabel : label `byvar' `l', strict
        di "{bf:-> by variables:} `by' = `bylabel'"
    }

    *** Drop invalid observations ***
    qui drop if missing(`hospitalid',`year',`class')

    *** Keep only the correct by-value
    if ("`by'"!="") qui keep if `byvar'==`l'

    *** Run regression
    di "run regressions residualizing dependent variable for controls"
    * If absorb or tfx_resid is not empty (only one is non-empty, otherwise an error was thrown), use areg
    if "`absorb'"!="" | "`tfx_resid'"!="" {
        areg `depvar' `controls' , absorb(`absorb'`tfx_resid')
    }
    * If absorb and tfx_resid are both empty, run regular regression
    else {
        reg `depvar' `controls'
    }

    *** Predict residuals
    sort `hospitalid' `year' `class'
    * If tfx_resid is empty, predict residuals
    if "`tfx_resid'"=="" {
        predict score_r1 if e(sample),r
    }
    * If tfx_resid was specified, predict residuals + absorbed teacher fixed effects
    else {
        qui predict score_r1 if e(sample), dresiduals
    }
    ** Adjust for shrinkage target if specified
    if "`shrinkage_target'" != "" {
        reg score_r1 `shrinkage_target'
        qui predict score_r if e(sample), res
        qui predict y_shrinktarget if e(sample), xb
            
        // Check if y_shrinktarget was created successfully
        capture confirm variable y_shrinktarget
        if _rc {
            di as error "Error: Failed to create y_shrinktarget variable"
            exit 111
        }
    }
    else {
        gen score_r = score_r1
    }
    *sum score_r, detail

    *** Save residuals to a dataset if merging them later
    if `merge_resid'==1 {
        tempfile resid_data_`l'
        qui save `"`resid_data_`l''"', replace
    }

    *** Save number of parameters

    tempname num_obs num_par

    scalar `num_obs' = e(N)

    * If absorb is not empty (and tfx_resid is), save (number of slopes + number of clusters + 1)
    if "`absorb'"!="" {
        scalar `num_par' = e(df_m) + e(df_a) + 1
    }
    * Otherwise, save (number of slopes + 1)
    else {
        scalar `num_par' = e(df_m) + 1
    }

    *** Create var for number of students in class
    tempvar n_tested
    qui bys `hospitalid' `year' `class': egen `n_tested' = count(score_r)
    *** Compute total variance ***
    tempvar class_mean index mshrinktarget
    qui by `hospitalid' `year' `class': egen `class_mean' = mean(score_r)
    qui by `hospitalid' `year' `class': g `index' = _n
    if "`shrinkage_target'" != "" {
        qui by `hospitalid' `year' `class': egen `mshrinktarget' = mean(y_shrinktarget)
    }

    tempname var_total
    qui sum score_r
    /// from looking I think this might var(A_it)
    scalar `var_total' = r(Var)*((`num_obs' - 1)/(`num_obs' - `num_par'))

    *** Compute individual variance (i.e. within class variance)
    *--> note that we use rmse instead of direct variance of residuals here to deal with fact that class effects have not been shrunk
    tempname num_class var_ind var_class

    tempvar individual_dev_from_class
    qui gen `individual_dev_from_class' = score_r - `class_mean'

    qui count if `index'==1 & `n_tested'!=0
    scalar `num_class' = r(N)

    qui sum `individual_dev_from_class'
    ///\hat{sigma_{epsilon}}^2
    scalar `var_ind' = r(Var)*((`num_obs' - 1)/(`num_obs' - `num_class' - `num_par' + 1))


    ********** Collapse to class-level data **********

    qui by `hospitalid' `year' `class': keep if _n==1


    *** Estimate covariance of two classes for same hospital in the same year
    set seed 9827496
    tempvar rand classnum
    g `rand'=uniform()
    bys `hospitalid' `year' (`rand'): gen `classnum'=_n

    * If there are multiple classes per hospital-year cell, compute the covariance.
    * Otherwise set to 0. Will display as missing in output, but internally set to 0 because it will never appear in the VCV, but the way things are coded requires that it be non-missing.
    tempname cov_sameyear corr_sameyear obs_sameyear
    qui sum `classnum'
    if (r(max)==1) {
        local missing_sameyear=1
        scalar `cov_sameyear'=0
    }
    else {
        local missing_sameyear=0
        tempvar identifier
        egen `identifier'=group(`hospitalid' `year')
        qui tsset `identifier' `classnum' /*, noquery*/
        qui corr `class_mean' f.`class_mean' [aw=`n_tested'+f.`n_tested'], cov
        scalar `cov_sameyear'=r(cov_12)
        scalar `corr_sameyear'=r(cov_12) / ( sqrt(r(Var_1)) * sqrt(r(Var_2)) )
        scalar `obs_sameyear'=r(N)
    }

    *** Compute the variance of the class-level shock.  Hits al lkids in the class in the same way, but is unrelated across classes even taught by the same teacher in the same year.
    /// this is variance_theta
    scalar `var_class' = `var_total' - `var_ind' - `cov_sameyear'
    if (`var_class'<0) {
        di as error "Note: var_class has been computed as being less than 0."
        di "var_class is defined as = var_total - var_ind - cov_sameyear."
        di "Computed variances: var_total, var_ind, cov_sameyear, var_class"
        di `var_total',`var_class',`var_ind',`cov_sameyear'
        di "This negative variance can occur because cov_sameyear is calculated using only the subsample of observations that teach multiple classes per year (in the same by-group)."
    }

    /* 2019-02-07 D.Staiger : change to code to allow us to make the M invertable using an eigen value trick.  */
    /* This will only work when we have only one classroom per teacher (e.g. no classrooms within hospital). */
    /* If you have multiple classrooms per teacher the code should work fine. This resets the diagonal of */
    /* the M matrix to be our estimate of the hospital-level variance (the original code set this to 0, */
    /* and put the hospital level variance into the class level variance (so it was part of weight). */
    /* Now, M will be what we want (mumu), i.e. it will have a the hospital variance along the diagnol */
    if (`missing_sameyear'==1) {
	scalar `cov_sameyear' = `var_class'
	scalar `var_class' = 0
    }
    // <END>
    tempvar weight
    qui g `weight'=1/(`var_class' + `var_ind'/`n_tested')

    *** Keep teacher-years which have no weight

    tempvar excess_weight
    qui gen `excess_weight'=(missing(`weight'))

    qui replace `weight'=1 if missing(`weight')
    * note: adding this weight doesn't affect the class_mean, because missing observations are not included
    * in the mean computation.  it only affects the rawsum of weight, and so we remove it afterward.


********** Collapse to teacher-year level data using precision weights **********
if "`shrinkage_target'" != "" {
    collapse (mean) `class_mean' `mshrinktarget'  (rawsum) `weight' `n_tested' `excess_weight' [aw=`weight'], by(`hospitalid' `year' `by') fast
}
else {
    collapse (mean) `class_mean' (rawsum) `weight' `n_tested' `excess_weight' [aw=`weight'], by(`hospitalid' `year' `by') fast
}

* Remove the excess weight used to keep missing scores
qui replace `weight'=`weight'-`excess_weight'

///DRIFT
*** Estimate the covariance of years t and t+i for every i, and store in vector m
qui tsset `hospitalid' `year'/*, noquery*/

tempvar minyear maxyear diff validyear minvalidyear maxvalidyear diffvalid

qui bys `hospitalid': egen `minyear'=min(`year')
qui by `hospitalid': egen `maxyear'=max(`year')
qui g `diff'=`maxyear'-`minyear'
qui sum `diff'
local maxspan=`r(max)'

qui gen `validyear'=`year' if !missing(`class_mean')
qui by `hospitalid': egen `minvalidyear'=min(`validyear')
qui by `hospitalid': egen `maxvalidyear'=max(`validyear')
qui g `diffvalid'=`maxvalidyear'-`minvalidyear'
qui sum `diffvalid'
local maxscorespan=`r(max)'

if (`maxscorespan'<`maxspan') & (`driftlimit'<=0) {
    di as error _n	"error: The maximum lags of teacher data is `maxspan', but the maximum lags of teacher data with class scores is `maxscorespan'."
    di as error		"       You must either set driftlimit() <= `maxscorespan', or drop observations so that the spans are no longer mismatched."
    exit 499
}
if (`driftlimit'>`maxscorespan') {
    di as error "error: driftlimit(`driftlimit') was specified, which is greater than the number of lags (`maxscorespan') in the data."
    exit 499
}

mata:CC=compute_cov_corr("`class_mean'","`n_tested'",`maxscorespan',"`hospitalid'")

if (`driftlimit'>0)	mata:m=create_m(CC[.,1],st_numscalar("`cov_sameyear'"),`maxspan',`driftlimit')
else				mata:m=create_m(CC[.,1],st_numscalar("`cov_sameyear'"))

/* Code addition by D.Staiger 2019-02-07 - to match changes made to other code.  */
di "Standard deviations: total, classes, students, Hospital same year"
if (`missing_sameyear'==0) di sqrt(`var_total'),sqrt(`var_class'),sqrt(`var_ind'),sqrt(`cov_sameyear')
else di sqrt(`var_total'),sqrt(`var_class'),sqrt(`var_ind'),sqrt(`cov_sameyear')

/* OLD CODE */
/* *** Print estimated variances and covariances */
/* if (`missing_sameyear'==0) di sqrt(`var_total'),sqrt(`var_class'),sqrt(`var_ind'),sqrt(`cov_sameyear') */
/* else di sqrt(`var_total'),sqrt(`var_class'),sqrt(`var_ind'),. */


di "Covariances (left), correlations (middle), observations (right).  Row i indicates the relation between year t and t+i:"
mata:CC[.,1..3]

di "Covariances used for VA computations:"
mata: m[2..length(m)]'

if (`driftlimit'>0) {
    di "Drift limit specified:"
    di `driftlimit'

    di "Covariances used for VA computations:"
    mata: m[2..length(m)]'
}

mata:check_m_nomissing(m)

*** Accumulate the estimated variances/covariances/correlations across by-vals
if (`firstloop'==1) {
    mata:cov_lag_accum= CC[.,1]
    mata:corr_lag_accum= CC[.,2]
    mata:obs_lag_accum= CC[.,3]
    mata:cov_se_lag_accum= CC[.,4]
    mata:var_total_accum=	st_numscalar("`var_total'")
    mata:var_class_accum=	st_numscalar("`var_class'")
    mata:var_ind_accum=	st_numscalar("`var_ind'")

    if (`missing_sameyear'==1) {
        mata:cov_sameyear_accum=.
        mata:corr_sameyear_accum=.
        mata:obs_sameyear_accum=0
    }
    else {
        mata:cov_sameyear_accum=st_numscalar("`cov_sameyear'")
        mata:corr_sameyear_accum=st_numscalar("`corr_sameyear'")
        mata:obs_sameyear_accum=st_numscalar("`obs_sameyear'")
    }
}
else {
    mata:cov_lag_accum=		rightAppendMatrices(cov_lag_accum,CC[.,1])
    mata:corr_lag_accum=	rightAppendMatrices(corr_lag_accum,CC[.,2])
    mata:obs_lag_accum=		rightAppendMatrices(obs_lag_accum,CC[.,3])
    mata:cov_se_lag_accum=	rightAppendMatrices(cov_se_lag_accum,CC[.,4])
    mata:var_total_accum=	var_total_accum,st_numscalar("`var_total'")
    mata:var_class_accum=	var_class_accum,st_numscalar("`var_class'")
    mata:var_ind_accum=		var_ind_accum,st_numscalar("`var_ind'")

    if (`missing_sameyear'==1) {
        mata:cov_sameyear_accum= cov_sameyear_accum,.
        mata:corr_sameyear_accum= corr_sameyear_accum,.
        mata:obs_sameyear_accum= obs_sameyear_accum,.
    }
    else {
        mata:cov_sameyear_accum=cov_sameyear_accum,st_numscalar("`cov_sameyear'")
        mata:corr_sameyear_accum=corr_sameyear_accum,st_numscalar("`corr_sameyear'")
        mata:obs_sameyear_accum=obs_sameyear_accum,st_numscalar("`obs_sameyear'")
    }
}

*********

* Count the number of obs for each hospital
sort `hospitalid' `year'
tempvar obs_hosp
by `hospitalid': egen `obs_hosp'=count(`hospitalid')

* Compute teacher VA
qui gen float tv=.

if ("`leaveout_years'"!="") {
    foreach v in `leaveout_vars' {
        qui gen float `v'=.
    }
    // Call mata function with leaveout parameters and debug flag
    mata: driftcalclist(vectorToStripeDiag(m), "`hospitalid'", "`year'", "`class_mean'", "`weight'", "`obs_hosp'", "tv", "`leaveout_years'", "`leaveout_vars'", "`debug'"!="")
}
else {
    // Call mata function without leaveout parameters but with debug flag
    mata: driftcalclist(vectorToStripeDiag(m), "`hospitalid'", "`year'", "`class_mean'", "`weight'", "`obs_hosp'", "tv", "", "", "`debug'"!="")
}

* Save the VA estimates to a dataset
local shrinkage_vars_to_keep 
if "`shrinkage_target'" != "" {
    local shrinkage_vars_to_keep   `mshrinktarget'
}
local leaveout_vars_to_keep 
if "`leaveout_years'" != "" {
    local leaveout_vars_to_keep  `leaveout_vars'
}
des *,full 
keep `hospitalid' `year' `by' tv `shrinkage_vars_to_keep' `leaveout_vars_to_keep'



///need to add back the hospital charactericis portion
if "`shrinkage_target'" != "" {
    if "`leaveout_vars'" != "" {
        foreach v in `leaveout_vars' {
            gen `v'_shrinktgt =  `v' + `mshrinktarget'        
            replace `v'_shrinktgt=`mshrinktarget' if `v'==.
        }
    }
    gen shrinktarget_base=`mshrinktarget'
}

if (`firstloop'!=1) {
    append using `"`output'"', nolabel
}
qui save `"`output'"', replace

* Turn firstloop counter off
local firstloop=0


di "{txt}{hline}"

* Save VA estimates
if "`output_addvars'"!="" quietly {
    restore, preserve
    keep `hospitalid' `year' `by' `output_addvars'
    bys `hospitalid' `year' `by' `output_addvars': keep if _n==1
    merge m:1 `hospitalid' `year' `by' using `"`output'"', nogen nolabel
}
sort `hospitalid' `year' `by'
qui save `"`output'"', replace

* Save "variances / covariances / correlations" dataset to csv
if ("`by'"!="") {
    local bylabels=""
    forvalues i=1/`by_vals' {
        local bylabel : label `byvar' `i', strict
        local bylabel=subinstr("`bylabel'"," ","_",.)
        local bylabels `bylabels' _`bylabel'
    }
    mata:saveVariancesToDataset(cov_lag_accum, corr_lag_accum, obs_lag_accum, cov_se_lag_accum, var_total_accum, var_class_accum, var_ind_accum, cov_sameyear_accum, corr_sameyear_accum, obs_sameyear_accum, tokens(st_local("bylabels")))
}
else mata:saveVariancesToDataset(cov_lag_accum, corr_lag_accum, obs_lag_accum, cov_se_lag_accum, var_total_accum, var_class_accum, var_ind_accum, cov_sameyear_accum, corr_sameyear_accum, obs_sameyear_accum, "")
if (`nooutput'!=1) qui outsheet using `"`output'_variance.csv"', comma replace


* Load the correct output dataset
tokenize "`data'"
if inlist("`1'","preserve","merge") {
    restore

    if (`merge_resid'==1) {
        if ("`byvar'"!="") qui keep if missing(`hospitalid',`year',`class',`byvar')
        else qui keep if missing(`hospitalid',`year',`class')
        forvalues l=1/`by_vals' {
            append using `"`resid_data_`l''"', nolabel
        }
    }
    if (`merge_tv'==1) qui merge m:1 `hospitalid' `year' `by' `output_addvars' using `"`output'"', nogen nolabel
    /* else "`data'"=="preserve", and that is already loaded. */
}
else {
    restore, not

    if ("`data'"=="tv") use `"`output'"', clear
    /* else "`data'"=="variance", and that is already loaded. */
}

* Close log
if (`nooutput'!=1) log close t
}
end
