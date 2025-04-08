*! Test comparing hospital_ebayes.ado and vam.ado outputs
*! This test ensures that hospital_ebayes.ado with leave-out parameters
*! produces equivalent results to vam.ado with quasi-experiment parameters

* Set up test environment
clear all
set more off
capture program drop hospital_ebayes
capture program drop vam

* Load the programs


* Load existing test data
use "test.dta", clear

* Sort data
sort id year

* Run hospital_ebayes with leave-out parameters
* Configure leave-out to match vam's quasi-experiment approach:
* - tv_2yr_l: exclude previous year (t-1)
* - tv_2yr_f: exclude next year (t+1)
* - tv_ss: exclude {t-2,t-1,t+1,t+2}
mata: mata clear
include "../src/hospital_ebayes.ado"

hospital_ebayes y, hospitalid(id) year(year) ///
    controls(xb) tfx_resid(id) ///
    leaveout_years("-1,0 0,1 -2,2") ///
    leaveout_vars("tv_2yr_l tv_2yr_f tv_ss") ///
    output("test_hospital_ebayes")

* Run vam with quasi-experiment
mata: mata clear
include "../test/vam.ado"
gen class=1
vam y, teacher(id) year(year) class(class) ///
    controls(xb) tfx_resid(id) ///
    quasiexperiment ///
    output("test_vam")

* Merge results for comparison
use "test_hospital_ebayes.dta", clear
rename tv tv_hospital_ebayes
rename tv_2yr_l tv_2yr_l_hospital_ebayes
rename tv_2yr_f tv_2yr_f_hospital_ebayes
rename tv_ss tv_ss_hospital_ebayes
tempfile hospital_ebayes_results
save `hospital_ebayes_results'

use "test_vam.dta", clear
rename tv tv_vam
rename tv_2yr_l tv_2yr_l_vam
rename tv_2yr_f tv_2yr_f_vam
rename tv_ss tv_ss_vam
merge 1:1 id year using `hospital_ebayes_results'

* Compare results
* Check if differences are within numerical precision
gen diff_tv = abs(tv_hospital_ebayes - tv_vam)
gen diff_2yr_l = abs(tv_2yr_l_hospital_ebayes - tv_2yr_l_vam)
gen diff_2yr_f = abs(tv_2yr_f_hospital_ebayes - tv_2yr_f_vam)
gen diff_ss = abs(tv_ss_hospital_ebayes - tv_ss_vam)

* Set tolerance for numerical differences
local tolerance = 1e-10

* Check if any differences exceed tolerance
quietly sum diff_*
local max_diff = `r(max)'

if `max_diff' > `tolerance' {
    display as error "Test failed: Differences between hospital_ebayes and vam exceed tolerance"
    display "Maximum difference: `max_diff'"
    display "Tolerance: `tolerance'"
    exit 498
} else {
    display as text "Test passed: hospital_ebayes and vam produce equivalent results"
    display "Maximum difference: `max_diff'"
    display "Tolerance: `tolerance'"
}

* Clean up temporary files
erase "test_hospital_ebayes.dta"
erase "test_vam.dta" 