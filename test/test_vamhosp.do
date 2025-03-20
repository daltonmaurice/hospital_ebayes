/*******************************************************************************
Test Suite: VAMHOSP
Purpose: Tests functionality of vamhosp estimator
*******************************************************************************/

clear all
mata: mata clear
capture log close
//mkdir test/logs 
//log using test/logs/test_results_vamhosp, text replace
cd /Users/mad265/git-pub/hospital_ebayes
cd test
use test.dta, clear

* Test 1: Basic Functionality
di "Test 1: Basic vamhosp Functionality"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp y , hospitalid(id)  year(year) ///
         data("merge tv") shrinkage_target(z) controls(xb) ///
         tfx_resid(id) 
        if _rc != 0 {
            di as error "Basic vamhosp failed with error code: " _rc
            exit _rc
        }
        
        count if !missing(tv)
        if r(N) == 0 {
            di as error "No value-added estimates generated"
            exit 498
        }
        di "✓ Basic functionality test passed"
    restore
}
if _rc != 0 exit _rc

* Test 2: Volume Effects
di "Test 2: Volume Effects with Shrinkage Target"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp dum_y, hospitalid(hospid) year(year) data("merge tv") shrinkage_target(lnvol)
        if _rc != 0 {
            di as error "Volume effects test failed with error code: " _rc
            exit _rc
        }
        
        count if !missing(tv)
        if r(N) == 0 {
            di as error "No value-added estimates generated with volume effects"
            exit 498
        }
        
        * Check if shrinkage target adjustment was applied
        count if !missing(hchar)
        if r(N) == 0 {
            di as error "No shrinkage target adjustments generated"
            exit 498
        }
        di "✓ Volume effects test passed"
    restore
}
if _rc != 0 exit _rc

* Test 3: Controls and Fixed Effects
di "Test 3: Controls and Fixed Effects"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado

        vamhosp dum_y, hospitalid(hospid) year(year) data("merge tv") ///
            controls(x) absorb(hospid) shrinkage_target(lnvol)
        if _rc != 0 {
            di as error "Controls and FE test failed with error code: " _rc
            exit _rc
        }
        
        count if !missing(tv)
        if r(N) == 0 {
            di as error "No value-added estimates generated with controls and FE"
            exit 498
        }
        di "✓ Controls and fixed effects test passed"
    restore
}
if _rc != 0 exit _rc

* Test 4: Drift Calculations
di "Test 4: Drift Calculations"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp dum_y, hospitalid(hospid) year(year) data("merge tv") char(lnvol)
        if _rc != 0 {
            di as error "Drift calculations failed with error code: " _rc
            exit _rc
        }
        
        foreach var in tv_2yr_l tv_2yr_f {
            count if !missing(`var')
            if r(N) == 0 {
                di as error "No estimates generated for `var'"
                exit 498
            }
        }
        di "✓ Drift calculations test passed"
    restore
}
if _rc != 0 exit _rc

* Test 5: Leave-out Estimators
di "Test 5: Leave-out Estimators"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        
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
            
        vamhosp y, hospitalid(id) year(year) ///
            controls(xb) char(z) data("merge tv") ///
            leaveout_years("`leaveout_patterns'") ///
            leaveout_vars("`leaveout_vars'")
        
        * Check if leave-out estimates were generated
        foreach var of local leaveout_vars {
            count if !missing(`var')
            if r(N) == 0 {
                di as error "No estimates generated for `var'"
                exit 498
            }
            
            * Check correlation with true mu
            corr `var' mu
            if abs(r(rho)) < 0.2 {
                di as error "Low correlation with true hospital effects for `var'"
                exit 498
            }
        }
        
        * Check that estimates differ appropriately
        * For example, tm1_t1 should differ from tm2_t2 due to different leave-out patterns
        corr tv_tm1_t1 tv_tm2_t2
        if abs(r(rho)) > 0.99 {
            di as error "Leave-out estimates too similar despite different patterns"
            exit 498
        }
        
        * Check that estimates follow expected drift pattern
        * Earlier years should have lower correlation with current estimates
        local prev_corr = 1
        foreach var in tv_tm1_t1 tv_tm2_t2 tv_tm3_t1 {
            corr tv `var'
            if r(rho) > `prev_corr' {
                di as error "Leave-out estimates don't follow expected drift pattern"
                exit 498
            }
            local prev_corr = r(rho)
        }
        
        di "✓ Leave-out estimators test passed"
    restore
}
if _rc != 0 exit _rc

* Test 6: Error Handling
di "Test 6: Error Handling"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        capture noisily vamhosp dum_y, hospitalid(nonexistent) year(year) data("merge tv") char(lnvol)
        if _rc == 0 {
            di as error "Failed to catch invalid variable error"
            exit 498
        }
        di "✓ Error handling test passed"
    restore
}
if _rc != 0 exit _rc

* Summary
di _n "Test Summary:"
di "✓ All vamhosp tests completed"
di "  - Basic functionality verified with continuous outcome"
di "  - Binary outcome estimation checked"
di "  - Hospital characteristics incorporation verified"
di "  - Drift calculations validated"
di "  - Leave-out estimators verified"
di "  - Error handling verified"

log close 