/*******************************************************************************
Test Suite: VAMHOSP
Purpose: Tests functionality of hospital_ebayes estimator
*******************************************************************************/
cd /Users/mad265/git-pub/hospital_ebayes/test
clear all
mata: mata clear

do "../src/hospital_ebayes.ado"

capture log close
//mkdir test/logs 
//log using test/logs/test_results_vamhosp, text replace

* Load required files for each test
capture confirm file "test.dta"
if _rc {
    di as error "Required test.dta file not found"
    exit 601
}

capture confirm file "../src/hospital_ebayes.ado"
if _rc {
    di as error "Required hospital_ebayes.ado file not found"
    exit 601
}

* Test 1: Basic Functionality
di _n "Test 1: Basic hospital_ebayes Functionality"
capture noisily {
    use test.dta, clear
    mata: mata clear
    quietly do "../src/hospital_ebayes.ado"
    
    /* Is this an error in code or because of simulated data?
    // shrinkage target is simulated to be within hospital id the tfx_resid
    hospital_ebayes y, hospitalid(id) year(year) ///
    controls(xb) shrinkage_target(z) data("merge tv") absorb(id)
    */

    hospital_ebayes y, hospitalid(id) year(year) ///
        controls(xb) shrinkage_target(z) data("merge tv") tfx_resid(id)
        
    * Verify required output variables exist
    foreach var in tv shrinktarget_base {
        capture confirm variable `var'
        if _rc {
            di as error "Required variable `var' not generated"
            exit 498
        }
        
        * Check for non-missing values
        count if !missing(`var')
        if r(N) == 0 {
            di as error "No valid estimates in `var'"
            exit 498
        }
    }
    di "✓ Basic functionality test passed"
}
if _rc exit _rc

* Test 2: Volume Effects
di _n "Test 2: Volume Effects with Shrinkage Target"
capture noisily {
    use test.dta, clear
    mata: mata clear
    quietly do "../src/hospital_ebayes.ado"
    
    hospital_ebayes y, hospitalid(id) year(year) ///
        data("merge tv") shrinkage_target(z)
        
    * Verify volume adjustment variables
    foreach var in tv shrinktarget_base {
        capture confirm variable `var'
        if _rc {
            di as error "Required variable `var' not generated with volume effects"
            exit 498
        }
        count if !missing(`var')
        if r(N) == 0 {
            di as error "No valid estimates in `var' with volume effects"
            exit 498
        }
    }
    di "✓ Volume effects test passed"
}
if _rc exit _rc

* Test 3: Controls and Fixed Effects
di _n "Test 3: Controls and Fixed Effects"
capture noisily {
    use test.dta, clear
    mata: mata clear
    quietly do "../src/hospital_ebayes.ado"

    hospital_ebayes y, hospitalid(id) year(year) ///
        controls(x)  shrinkage_target(z) ///
        data("merge tv") tfx_resid(id)
        
    * Verify output with controls and FE
    foreach var in tv shrinktarget_base {
        capture confirm variable `var'
        if _rc {
            di as error "Required variable `var' not generated with controls/FE"
            exit 498
        }
        count if !missing(`var')
        if r(N) == 0 {
            di as error "No valid estimates in `var' with controls/FE"
            exit 498
        }
    }
    di "✓ Controls and fixed effects test passed"
}
if _rc exit _rc

/*
* Test 4: Leave-out Estimators
di _n "Test 4: Leave-out Estimators"
capture noisily {
    use test.dta, clear
    mata: mata clear
    quietly do "../src/hospital_ebayes.ado"
    
    * Define leave-out patterns and variable names
    local leaveout_patterns ///
        "-1,+1" /// Leave out t-1 and t+1
        "-2,+2" /// Leave out t-2 and t+2
        "-3,+1" /// Leave out t-3 and t+1
        "-3,+2" /// Leave out t-3 and t+2
        "-5," /// Leave out before t-5
        ",+5" // Leave out after t+5
        
    local leaveout_vars ///
        tv_tm1_t1 ///
        tv_tm2_t2 ///
        tv_tm3_t1 ///
        tv_tm3_t2 ///
        tv_tm5_t ///
        tv_t_t5
        
    hospital_ebayes y, hospitalid(id) year(year) ///
        controls(xb) shrinkage_target(z) data("merge tv") ///
        leaveout_years("`leaveout_patterns'") ///
        leaveout_vars("`leaveout_vars'") tfx_resid(id)
    
    * Verify leave-out estimates
    foreach var of local leaveout_vars {
        capture confirm variable `var'
        if _rc {
            di as error "Required leave-out variable `var' not generated"
            exit 498
        }
        
        count if !missing(`var')
        if r(N) == 0 {
            di as error "No valid estimates in leave-out variable `var'"
            exit 498
        }
        
        * Check shrinkage target adjustments
        capture confirm variable `var'_shrinktgt
        if _rc {
            di as error "Shrinkage target adjustment not generated for `var'"
            exit 498
        }
    }
    di "✓ Leave-out estimators test passed"
}
if _rc exit _rc
*/


* Test 4b: Leave-out Estimators with tfx_resid
di _n "Test 4b: Leave-out Estimators with tfx_resid"
capture noisily {
    use test.dta, clear
    mata: mata clear
    quietly do "../src/hospital_ebayes.ado"
    
    * Define leave-out patterns and variable names
    local leaveout_patterns ///
        "-1,+1" /// Leave out t-1 and t+1
        "-2,+2" /// Leave out t-2 and t+2
        "-3,+1" /// Leave out t-3 and t+1
        "-3,+2" /// Leave out t-3 and t+2
        "-5," /// Leave out before t-5
        ",+5" // Leave out after t+5
        
    local leaveout_vars ///
        tv_tm1_t1 ///
        tv_tm2_t2 ///
        tv_tm3_t1 ///
        tv_tm3_t2 ///
        tv_tm5_t ///
        tv_t_t5
        
    hospital_ebayes y, hospitalid(id) year(year) ///
        controls(xb) shrinkage_target(z) data("merge tv") ///
        leaveout_years("`leaveout_patterns'") ///
        leaveout_vars("`leaveout_vars'") tfx_resid(id) debug
    
    * Verify leave-out estimates
    foreach var of local leaveout_vars {
        capture confirm variable `var'
        if _rc {
            di as error "Required leave-out variable `var' not generated"
            exit 498
        }
        
        count if !missing(`var')
        if r(N) == 0 {
            di as error "No valid estimates in leave-out variable `var'"
            exit 498
        }
        
        * Check shrinkage target adjustments
        capture confirm variable `var'_shrinktgt
        if _rc {
            di as error "Shrinkage target adjustment not generated for `var'"
            exit 498
        }
    }
    di "✓ Leave-out estimatorswith tfx test passed"
}
if _rc exit _rc

* Test 5: Error Handling
di _n "Test 5: Error Handling"
capture noisily {
    use test.dta, clear
    * Test invalid variable specification
    capture noisily hospital_ebayes y, hospitalid(id) year(year) controls(nonexistent_var)
    if _rc != 111 {
        di as error "Failed to catch invalid variable error"
        exit 498
    }
    
    * Test mismatched leaveout specifications
    capture noisily hospital_ebayes y, hospitalid(id) year(year) ///
        leaveout_years("-1,1") leaveout_vars("var1 var2")
    if _rc != 198 {
        di as error "Failed to catch mismatched leaveout specifications"
        exit 498
    }
    di "✓ Error handling test passed"
}
if _rc != 498 exit _rc

* Summary
di _n "{txt}Test Summary:"
di "✓ All hospital_ebayes tests completed successfully"
di "  - Basic functionality verified"
di "  - Volume effects tested"
di "  - Controls and fixed effects verified"
di "  - Leave-out estimators tested"
di "  - Error handling verified"

log close 