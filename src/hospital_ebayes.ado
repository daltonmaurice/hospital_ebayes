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
    leaveout_vars(string)]  /// New parameter for variable mappings

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
    
    // Parse variable names
    tokenize `leaveout_vars'
    forvalues i = 1/`n_rules' {
        local var_`i' "``i''"
        capture confirm variable ``i'', exact
        if (_rc==0) {
            di as error "The dataset loaded in memory cannot have a variable named ``i''."
            exit 110
        }
        qui gen float ``i'' = .
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
    qui sum score_r, detail

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

    // Call mata function with leaveout parameters
    mata: driftcalclist(vectorToStripeDiag(m), "`hospitalid'", "`year'", "`class_mean'", "`weight'", "`obs_hosp'", "tv", "`leaveout_years'", "`leaveout_vars'")
}
else {
    // Call mata function without leaveout parameters
    mata: driftcalclist(vectorToStripeDiag(m), "`hospitalid'", "`year'", "`class_mean'", "`weight'", "`obs_hosp'", "tv")
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


version 11
set matastrict on

mata:
    real rowvector computeweights(real matrix M, real scalar i, real colvector c, | real colvector weights) {
        // Add safety checks
        if (rows(M) != cols(M)) {
            printf("Error: Non-square matrix M (%f x %f)\n", rows(M), cols(M))
            _error(3205, "Matrix must be square")
        }
        
        if (missing(M)) {
            printf("Warning: Matrix M contains missing values\n")
        }
        
        real matrix X
        real matrix L
        real matrix vcv
        real matrix Mpos

	// construct matrix A which is used to select the relevant elements of M in constructing the VCV matrix
	real matrix temp
	real matrix A
	temp=designmatrix(c)

        /* ************************************************************************  */
        /* *** Make M matrix which is off diagnol */
        /* ************************************************************************  */
        /* Base of code adapted from Doug Staiger, added 8/30/2019 */
        /* NOW fix vcv so that it is pos semi def (with block/n will always */
        /* be invertable see higham, NJ, 1988 "computing a nearest symetric */
        /* pos sem def matrix I do this by maintianing the estimates of sd */
        /* of each signal, and fixing the corr matrix so take pos semi def */
        /* part of vcv, use it to estimate corr(vcv), then */
        /* vcvpos = corr(vcv):*(sd*sd') */
        X=.
        L=.
        symeigensystem(M,X,L)
        Mpos = X*diag(L:*(L:>=0))*X'
        /* The original code just used M everywhere, which is a matrix that is fed into this */
        A = temp, J(rows(c),cols(Mpos)-cols(temp),0)
        /* use A to select elements of M and build the VCV.  The second term adjusts the diagonal */
        /* elements of the VCV matrix to account for the class-level and individual-level shocks */
        /* We want to make the underlying signal matrix */
        if (args()==4) vcv=A*Mpos*A' + diag(1:/weights)
        else vcv=A*Mpos*A'
        // phi is the vector of autocovariances, selected correctly using the matrix A.
        real rowvector phi
        phi=Mpos[i,.]*A'

        /* return the vector of weights, choose the VCV that D.Staiger */
        /* coded  to always be pos semi def */
        return    (phi*cholinv(vcv))
}



real matrix compute_cov_corr(string scalar scores_var, string scalar weight_var, real scalar dim, string scalar hospitalid_var) {

    // pre-allocate matrix
    real matrix CC
    CC = J(dim,4,.)

    // Fill cov's and corr's: between time t and t+i
    real scalar i
    real scalar tstat
    for (i=1; i<=dim; i++) {
        // check that there are >=2 obs, in order to compute covariance
        stata(invtokens(("quietly count if !missing(",scores_var,",f",strofreal(i),".",scores_var,")"),""))
        if (st_numscalar("r(N)")>1) {
            stata(invtokens(("quietly corr ",scores_var," f",strofreal(i),".",scores_var," [aw=",weight_var,"+f",strofreal(i),".",weight_var,"], cov"),""))
            CC[i,1]=st_numscalar("r(cov_12)")
            CC[i,2]=CC[i,1] / ( sqrt(st_numscalar("r(Var_1)")) * sqrt(st_numscalar("r(Var_2)")) )
        }
        CC[i,3]=st_numscalar("r(N)")

        // Compute SE for covariance estimate
        if (st_numscalar("r(N)")>1) {
            stata(invtokens(("quietly reg ",scores_var," f",strofreal(i),".",scores_var," [aw=",weight_var,"+f",strofreal(i),".",weight_var,"], cluster(",hospitalid_var,")"),""))
            tstat=st_matrix("e(b)")[1,1] / sqrt( st_matrix("e(V)")[1,1] )
            CC[i,4]=abs(CC[i,1]/tstat)
        }
    }

    return (CC)
}

real rowvector create_m(real colvector lag_covariances, real scalar cov_sameyear, | real scalar lagdim, real scalar driftlimit) {
    // Add debugging
    printf("lag_covariances dimensions: %f x %f\n", rows(lag_covariances), cols(lag_covariances))
    printf("cov_sameyear: %f\n", cov_sameyear)
    
    real rowvector m

    if (args()==2)	m=cov_sameyear,lag_covariances'
else {
    if (length(lag_covariances)<driftlimit) _error("driftlimit specified is higher than the number of lags in the dataset")
    m=cov_sameyear,lag_covariances'[1..driftlimit],J(1,lagdim-driftlimit,lag_covariances[driftlimit])
}

return (m)
}

void check_m_nomissing(real rowvector m) {
    if (missing(m)>0) _error("covariance vector contains missing values")
}

real matrix vectorToStripeDiag(real vector m) {
    // Add debugging
    printf("Input vector m dimensions: %f x %f\n", rows(m), cols(m))
    
    real scalar dim
    dim = length(m)

    // pre-allocate matrix M
    real matrix M
    M=J(dim,dim,.)

    // fill lower triangle of M
    real scalar i
    real scalar j
    for (i=1; i<=dim; i++) {
        for (j=i; j<=dim; j++) {
            M[j,i]=m[j-i+1]
        }
    }

    _makesymmetric(M)
    return (M)
}

real matrix rightAppendMatrices(real matrix A, real matrix B) {
    real scalar rA
    real scalar rB
    rA=rows(A)
    rB=rows(B)

    if (rA==rB)		return (A,B)
    else if (rA<rB)	return ( ( A \ J(rB-rA,cols(A),.) ) , B )
    else			return ( A , ( B \ J(rA-rB,cols(B),.) ) )
}

void saveVariancesToDataset(real matrix cov_lag_accum, real matrix corr_lag_accum, real matrix obs_lag_accum, real matrix cov_se_lag_accum, real rowvector var_total_accum, real rowvector var_class_accum, real rowvector var_ind_accum, real rowvector cov_sameyear_accum, real rowvector corr_sameyear_accum, real rowvector obs_sameyear_accum, string rowvector suffixes) {

    stata("clear")

    // count number of lags, create correct number of obs, generate variable for number of lags
    real scalar n_lags
    n_lags=rows(cov_lag_accum)

    real scalar null
    null=st_addvar("int","lag")

    st_addobs(n_lags)
    stata("qui replace lag=_n")
    st_addobs(1)

    // generate output variables
    st_store(1::n_lags, st_addvar("float", "cov_lag":+suffixes), cov_lag_accum)
    st_store(1::n_lags, st_addvar("float", "corr_lag":+suffixes), corr_lag_accum)
    st_store(1::n_lags, st_addvar("float", "obs_lag":+suffixes), obs_lag_accum)
    st_store(1::n_lags, st_addvar("float", "cov_se_lag":+suffixes), cov_se_lag_accum)
    st_store(n_lags+1, st_addvar("float", "var_total":+suffixes), var_total_accum)
    st_store(n_lags+1, st_addvar("float", "var_class":+suffixes), var_class_accum)
    st_store(n_lags+1, st_addvar("float", "var_ind":+suffixes), var_ind_accum)
    st_store(n_lags+1, st_addvar("float", "cov_sameyear":+suffixes), cov_sameyear_accum)
    st_store(n_lags+1, st_addvar("float", "corr_sameyear":+suffixes), corr_sameyear_accum)
    st_store(n_lags+1, st_addvar("float", "obs_sameyear":+suffixes), obs_sameyear_accum)
}

real scalar driftcalc(real matrix M, real scalar i, real colvector c, real colvector weights, real colvector scores) {

    // b is the vector of weights
    real rowvector b
    b=computeweights(M, i, c, weights)
    // return the computed tv estimate -- where it basically is summing up all the
    // scores * weight - by matrix mulitplication of row and column vector
    return (b*scores)
}


void driftcalclist(real matrix M, string scalar hospitalid_var, string scalar time_var, 
    string scalar scores_var, string scalar weights_var, string scalar hospobs_var, 
    string scalar va_var, | string scalar leaveout_years, string scalar leaveout_vars) {
    
    // Declare all variables upfront
    real scalar nobs, obs, hospitalid, obs_hosp, time, new_hospitalid, new_time, year_index, i
    real matrix Z, Z_hosp, Z_obs, Z_quasi
    
    nobs = st_nobs()
    
    // Get variable indices
    real scalar hospitalid_var_ind, time_var_ind, hospobs_var_ind, va_var_ind
    hospitalid_var_ind = st_varindex(hospitalid_var)
    time_var_ind = st_varindex(time_var)
    hospobs_var_ind = st_varindex(hospobs_var)
    va_var_ind = st_varindex(va_var)
    
    // Create view of variables
    st_view(Z=., ., (hospitalid_var, time_var, weights_var, scores_var))
    
    // Initialize
    hospitalid = .
    time = .
    
    // Loop over observations
    for (obs=1; obs<=nobs; obs++) {
        new_hospitalid = _st_data(obs, hospitalid_var_ind)
        new_time = _st_data(obs, time_var_ind)
        
        // Only perform calculations for new hospital-year
        if (new_time != time | new_hospitalid != hospitalid) {
            time = new_time
            
            if (new_hospitalid != hospitalid) {
                hospitalid = new_hospitalid
                obs_hosp = _st_data(obs, hospobs_var_ind)
                st_subview(Z_hosp=., Z, (obs, obs+obs_hosp-1), .)
                year_index = min(Z_hosp[.,2])-1
            }
            
            // Get observations excluding current year
            Z_obs = select(Z_hosp, Z_hosp[.,2]:!=time)
            Z_obs = select(Z_obs, Z_obs[.,4]:!=.)
            
            // Compute standard VA
            if (rows(Z_obs) > 0) {
                st_store(obs, va_var_ind, 
                    driftcalc(M, time-year_index, Z_obs[.,2]:-year_index, Z_obs[.,3], Z_obs[.,4]))
            }
            string vector lyears, lvars
            lyears = tokens(leaveout_years)
            lvars = tokens(leaveout_vars)
            // Compute leaveout estimates if specified
            if (args()>7) {
                for (i=1; i<=length(leaveout_years); i++) {
                    string scalar before, after
                    _parse_rule(lyears[i], before, after)
                    
                    // Get base observations
        
                    Z_quasi = Z_obs
                    
                    // Apply filters if valid
                    if (before != "" & before != " ") {
                        Z_quasi = select(Z_quasi, Z_quasi[.,2] :< (time + strtoreal(before)))
                    }
                    if (after != "" & after != " ") {
                        Z_quasi = select(Z_quasi, Z_quasi[.,2] :> (time + strtoreal(after)))
                    }
                    
                    if (rows(Z_quasi) > 0) {
                        st_store(obs, st_varindex(lvars[i]), 
                            driftcalc(M, time-year_index, Z_quasi[.,2]:-year_index, Z_quasi[.,3], Z_quasi[.,4]))
                    }
                }
            }
        }
    }
}

// Helper function to parse leaveout rules
void _parse_rule(string scalar rule, string scalar before, string scalar after) {
    string vector parts
    parts = tokens(rule, ",")
    before = parts[1]
    after = parts[3]
}
end

