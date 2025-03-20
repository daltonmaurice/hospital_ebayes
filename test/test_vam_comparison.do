/*******************************************************************************
Test Suite: VAM Comparison
Purpose: Compares results between vamhosp  and vam implementations
*******************************************************************************/

capture log close
log using test_results_vam, text replace

use test.dta, clear
cap restore 
* Test 1: Basic Functionality Comparison
di "Test 1: Basic Functionality Comparison"
capture noisily {
    preserve
        * Run vamhclose
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp price , hospitalid(hospid) year(year) data("merge tv") char(lnvol) 
        if _rc != 0 {
            di as error "vamhosp failed with error code: " _rc
            exit _rc
        }
        rename tv tv_hosp
        
        * Run vam
        mata: mata clear
        do ../src/vam_hclose.ado
        gen tt=1
        egen cgroup=group(hospid year tt)

        vam dum_y, teacher(hospid) year(year) data("merge tv") hospchar(lnvol) class(cgroup)
        if _rc != 0 {
            di as error "vam failed with error code: " _rc
            exit _rc
        }
        rename tv tv_vam
        
        * Compare results
        capture noisily corr tv_hclose tv_vam
        if _rc != 0 {
            di as error "Correlation calculation failed with error code: " _rc
            exit _rc
        }
        local corr = r(rho)
        
        if abs(`corr' - 1) >= 0.01 {
            di as error "Correlation test failed. Expected correlation near 1, got: `corr'"
            exit 9
        }
        di "✓ Basic functionality correlation: `corr'"
    restore
}
if _rc != 0 exit _rc

* Test 2: Volume Effects
di "Test 2: Volume Effects Comparison"
capture noisily {
    preserve
        mata: mata clear
        do ../src/hospital_ebayes.ado
        vamhosp dum_y, hospitalid(hospid) year(year) char(lnvol) charflag(1) data('merge tv')
        rename tv_hosp_bss_t_t1_t2 tv_hosp_vol
        
        mata: mata clear
        do ../src/vam_hclose.ado
        gen cgroup=1
        vam dum_y, teacher(hospid) year(year) hospchar(lnvol) hospvolumeflag(1) data('merge tv') class(cgroup)
        rename tv tv_vam_vol
        
        corr tv_hclose_vol tv_vam_vol
        assert abs(r(rho) - 1) < 0.01
        di "✓ Volume effects correlation: " r(rho)
    restore
}


log close 