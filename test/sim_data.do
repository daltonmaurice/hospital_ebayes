* sim_data.do */
/*  Doug Staiger, 3/17/2025 */
/*  simulates data to test our empirical Bayes code */
/*   for estimating hospital quality with drift */
/*  saves data in sim_data */
/*   and then checks to be sure we can recover the parameters */
/*  (note I do this for continuous dep var y) */
/*  (I could replace y with dummy dy, but not sure what params would recover) */

cd "/Users/mad265/git-pub/ebayes-hospital/src"
cap log close
log using log/sim_data, text replace

clear all
set seed 7

/* set globals here */
global delta=.95 /* ar(1) correlation */
global nhosp=4000 /* number of hospitals */
global nyear=20 /* number of years per hospital */
global npat=50 /* number of patients per hospital */
/* var mu is normalized to 1 */
global sd_z = .1 /* SD of the hospital-level shrinkage target */
global sd_xb = 1 /* SD of the patient-level risk adjuster index xb */
global sd_resid = 5 /* SD of patient-level residual */

/* finished with globals */

set obs $nhosp
gen mu1=rnormal()
forval i=2/$nyear {
local j = `i'-1
gen mu`i' = $delta*mu`j' + sqrt(1-$delta^2)*rnormal()
}
summ mu*
corr mu*

gen id=_n
reshape long mu, i(id) j(year)
summ
gen z=$sd_z*rnormal()
gen muz=mu+z
expand $npat
gen xb=$sd_xb*rnormal()
gen y=muz+xb+$sd_resid*rnormal()
gen dy=(y>$sd_resid)
summ
compress
save test, replace

/* can we correctly recover the moments? */

* regress y on xb and hosp-year fixed effects and save fe+e=y-xb
qui {
egen id_yr=group(id year)
areg y xb, absorb(id_yr)
predict fe_e, dresid
global sd_e=e(rmse)
* regress fe_e on z, save residual = mu+e and get the rmse = sd(e)
reg fe_e z
predict mu_e, resid
* also get sd of mu_e for later 
summ mu_e
global sd_mu_e=r(sd)
* here is estimate of sd of mu 
noi disp "The estimated signal variance of mu is = " (($sd_mu_e)^2 - ($sd_e)^2) "; it should be = 1"


* collapse data to hospital-year level
collapse mu_e, by(id year)
tsset id year
global loopyr=$nyear-1
qui forval i=0/$loopyr {
	corr mu_e L`i'.mu_e, cov
	noi disp "the estimated covar at " `i' " lags is = " r(cov_12) "; it should be = " (`i'!=0)*($delta)^(`i') + (`i'==0)*((1+(($sd_resid)^2)/$npat))
}

}
clear

/* repeat replacing y with dy */
/* not sure what estimates should be exactly, but assume that covariance is (delta^lag)*signal variance of mu. */
/*  and used sd_e from fe regression in place of sd_resid */

use test, clear
replace y=dy
* regress y on xb and hosp-year fixed effects and save fe+e=y-xb
qui {
egen id_yr=group(id year)
areg y xb, absorb(id_yr)
predict fe_e, dresid
global sd_e=e(rmse)
* regress fe_e on z, save residual = mu+e and get the rmse = sd(e)
reg fe_e z
predict mu_e, resid
* also get sd of mu_e for later 
summ mu_e
global sd_mu_e=r(sd)
* here is estimate of sd of mu 
noi disp "The estimated signal variance of mu is = " (($sd_mu_e)^2 - ($sd_e)^2) "; Use this estimate to calculate expected below"
global sig_var=(($sd_mu_e)^2 - ($sd_e)^2)

* collapse data to hospital-year level
collapse mu_e, by(id year)
tsset id year
global loopyr=$nyear-1
qui forval i=0/$loopyr {
	corr mu_e L`i'.mu_e, cov
	noi disp "the estimated covar at " `i' " lags is = " r(cov_12) "; it should be = " (`i'!=0)*($sig_var)*($delta)^(`i') + (`i'==0)*(($sig_var+(($sd_e)^2)/$npat))
}


}
