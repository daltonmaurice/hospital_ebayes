/*******************************************************************************
Test Suite: hospital_ebayes 
Purpose: Comprehensive testing of hospital_ebayes functionality
Reference: Stata 18 Programming Manual [P] error
*******************************************************************************/

capture log close
log using test_results_hospital_ebayes, text replace

* Load ado file once at the start
mata: mata clear
do ../src/hospital_ebayes.ado

* Helper program to check correlation
capture program drop check_correlation
program define check_correlation
    args var1 var2 threshold
    corr `var1' `var2'
    if abs(r(rho) - 1) >= `threshold' {
        di as error "Correlation test failed. Expected correlation near 1, got: " r(rho)
        exit 9
    }
    di "✓ Correlation between `var1' and `var2': " r(rho)
end

* Test counter
local test_number = 0
local failed_tests = 0

* Test 1: Basic Functionality
local ++test_number
di _n "Test `test_number': Basic Functionality"
capture noisily {
    use test.dta, clear
    hospital_ebayes y, hospitalid(id) year(year) data("merge tv")
}
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 2: Controls
local ++test_number
di _n "Test `test_number': Controls"
capture noisily {
    use test.dta, clear
    hospital_ebayes y, hospitalid(id) year(year) controls(xb) data("merge tv")
}
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 3: Leave-out Estimators
local ++test_number
di _n "Test `test_number': Leave-out Estimators"
capture noisily {
    mata: mata clear
 do ../src/hospital_ebayes.ado

    use test.dta, clear
    hospital_ebayes y, hospitalid(id) year(year) ///
        leaveout_years("-2,2 -1,1") leaveout_vars("tv_2yr tv_1yr") data("merge tv")
}
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 4: Shrinkage Targets
local ++test_number
di _n "Test `test_number': Shrinkage Targets"
capture noisily {
    use test.dta, clear
    hospital_ebayes y, hospitalid(id) year(year) ///
        shrinkage_target(z) data("merge tv")
}
if _rc {
    local ++failed_tests
    di as error "Test `test_number' failed with error code: " _rc
}
else di "✓ Test `test_number' passed"

* Test 5: Data Handling Options
foreach opt in "preserve" "tv" "merge tv" "variance" {
    local ++test_number
    di _n "Test `test_number': Data Handling Option - `opt'"
    capture noisily {
        use test.dta, clear
        hospital_ebayes y, hospitalid(id) year(year) data("`opt'")
    }
    if _rc {
        local ++failed_tests
        di as error "Test `test_number' failed with error code: " _rc
    }
    else di "✓ Test `test_number' passed"
}

* Test 6: Error Cases - Invalid driftlimit
local ++test_number
di _n "Test `test_number': Error Cases - Invalid driftlimit"
capture noisily {
    use test.dta, clear
    hospital_ebayes y, hospitalid(id) year(year) driftlimit(999)
}
if _rc != 499 {  // Expecting specific error code 499
    local ++failed_tests
    di as error "Test `test_number' failed: Expected error 499, got " _rc
}
else di "✓ Test `test_number' passed"

* Test 7: Error Cases - Invalid variable names
local ++test_number
di _n "Test `test_number': Error Cases - Invalid variable names"
capture noisily {
    use test.dta, clear
    gen tv = 1
    hospital_ebayes y, hospitalid(id) year(year)
}
if _rc != 110 {  // Expecting specific error code 110 for variable already exists
    local ++failed_tests
    di as error "Test `test_number' failed: Expected error 110, got " _rc
}
else di "✓ Test `test_number' passed"

* Test 8: Combined Features
local ++test_number
di _n "Test `test_number': Combined Features"
capture noisily {
    use test.dta, clear
    hospital_ebayes y, hospitalid(id) year(year) ///
        controls(xb) ///
        leaveout_years("-1,1") leaveout_vars("tv_1yr") ///
        shrinkage_target(z) ///
        data("merge tv")
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
    exit 9
}
else {
    di as txt "All tests passed successfully."
}

log close 