/*******************************************************************************
Test Suite: VAMHOSP
Purpose: Tests functionality of hospital_ebayes estimator
*******************************************************************************/
cd /Users/mad265/git-pub/hospital_ebayes/src 
mata: mata clear

* Ensure logs directory exists
capture mkdir ../test/logs

* Close any existing logs and start new one
capture log close
capture log using ../test/logs/test_results_vamhosp, text replace
if _rc {
    di as error "Failed to create log file. Error code: " _rc
    exit 498
}

* Test counter and failure tracking
local test_number = 0
local failed_tests = 0

* Load required files for each test
capture confirm file "../test/test.dta"
if _rc {
    di as error "Required test.dta file not found"
    exit 601
}

capture confirm file "hospital_ebayes.ado"
if _rc {
    di as error "Required hospital_ebayes.ado file not found"
    exit 601
}

* Test 1: Basic Functionality
local ++test_number
di _n "Test `test_number': Basic hospital_ebayes Functionality"
capture noisily {
    use ../test/test.dta, clear
    mata: mata clear
    quietly do "hospital_ebayes.ado"
    
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
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 2: Volume Effects
local ++test_number
di _n "Test `test_number': Volume Effects with Shrinkage Target"
capture noisily {
    use ../test/test.dta, clear
    mata: mata clear
    quietly do "hospital_ebayes.ado"
    
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
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 3: Controls and Fixed Effects
local ++test_number
di _n "Test `test_number': Controls and Fixed Effects"
capture noisily {
    use ../test/test.dta, clear
    mata: mata clear
    quietly do "hospital_ebayes.ado"

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
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 4: Leave-out Estimators with tfx_resid
local ++test_number
di _n "Test `test_number': Leave-out Estimators with tfx_resid"
capture noisily {
    use ../test/test.dta, clear
    mata: mata clear
    quietly do "hospital_ebayes.ado"
    
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
    di "✓ Leave-out estimators with tfx test passed"
}
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 5: Error Handling
local ++test_number
di _n "Test `test_number': Error Handling"
capture noisily {
    use ../test/test.dta, clear
    * Test invalid variable specification
    capture noisily hospital_ebayes y, hospitalid(id) year(year) controls(nonexistent_var)
    if _rc != 111 {
        di as error "Failed to catch invalid variable error, got error code: " _rc
        exit 498
    }
    
    * Test mismatched leaveout specifications
    capture noisily {
        use ../test/test.dta, clear
        hospital_ebayes y, hospitalid(id) year(year) ///
            leaveout_years("-1,1") leaveout_vars("var1 var2") data("merge tv")
        
        * If we get here, the command succeeded when it should have failed
        di as error "Expected command to fail with mismatched leaveout specifications"
        exit 498
    }
    if _rc == 0 {
        di as error "Command succeeded when it should have failed"
        exit 498
    }
    di "✓ Error handling test passed"
}
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Final test results
di _n "Test Summary"
di "=============="
di "Total tests run: `test_number'"
di "Tests failed: `failed_tests'"

if `failed_tests' > 0 {
    di as error "Some tests failed. Check log for details."
    log close
    exit 9
}
else {
    di as txt "All tests passed successfully."
    log close
} 