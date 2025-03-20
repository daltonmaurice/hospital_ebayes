/*******************************************************************************
Test Suite: VAMHOSP
Purpose: Tests functionality of vamhosp estimator
*******************************************************************************/

clear all
mata: mata clear
capture log close
log using test_results_vamhosp, text replace

use test.dta, clear

* Test 1: Basic Functionality
di "Test 1: Basic vamhosp Functionality"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp y , hospitalid(id)  year(year) ///
         data("merge tv") char(z) controls(xb) ///
         absorb(id)  quasi before(1)
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
exit
* Test 2: Volume Effects
di "Test 2: Volume Effects"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp dum_y, hospitalid(hospid) year(year) data("merge tv") hospchar(lnvol)
        if _rc != 0 {
            di as error "Volume effects test failed with error code: " _rc
            exit _rc
        }
        
        count if !missing(tv)
        if r(N) == 0 {
            di as error "No value-added estimates generated with volume effects"
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
        egen fe_group = group(class)
        vamhosp dum_y, hospitalid(hospid) year(year) data("merge tv") ///
            controls(x) absorb(fe_group) char(lnvol)
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

* Test 5: By-Group Analysis
di "Test 5: By-Group Analysis"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        gen group = mod(_n, 2)
        vamhosp dum_y, hospitalid(hospid) year(year) data("merge tv") by(group) char(lnvol)
        if _rc != 0 {
            di as error "By-group analysis failed with error code: " _rc
            exit _rc
        }
        
        forvalues g = 0/1 {
            count if !missing(tv) & group == `g'
            if r(N) == 0 {
                di as error "No estimates generated for group `g'"
                exit 498
            }
        }
        di "✓ By-group analysis test passed"
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
di "  - Basic functionality verified"
di "  - Volume effects checked"
di "  - Controls and fixed effects tested"
di "  - Drift calculations validated"
di "  - By-group analysis confirmed"
di "  - Error handling verified"

log close 