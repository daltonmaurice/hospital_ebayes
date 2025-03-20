
/* now try estimating the rolling filter */
/*  This data is especially easy, because everyone has the same volume. */
/*  But more generally we'll have to incorporate volume */
/*   (maybe by just running this code on quintiles of volume or something) */

/* set globals here */
global delta=.95 /* ar(1) correlation */
global nhosp=4000 /* number of hospitals */
global nyear=20 /* number of years per hospital */
global npat=50 /* number of patients per hospital */
/* var mu is normalized to 1 */
global sd_z = .1 /* SD of the hospital-level shrinkage target */
global sd_xb = 1 /* SD of the patient-level risk adjuster index xb */
global sd_resid = 5 /* SD of patient-level residual */

set seed 7
clear
use test,clear 
* first create hospital resids by year
*  initialize the variables
gen fe_e=.
gen mu_e=.
gen zbhat=.
gen nobs=.

qui forval y=1/$nyear {
  areg y xb if year==`y', absorb(id)
  predict fe_e`y' if e(sample), dresid
  replace fe_e=fe_e`y' if year==`y'
  replace nobs=e(sample) if year==`y'
  global sd_e`y'=e(rmse)
  global var_e`y'=${sd_e`y'}^2
* regress fe_e on z, save residual = mu+e and get the rmse = sd(e)
  reg fe_e z if year==`y'
  predict mu_e`y' if e(sample), resid
  replace mu_e=mu_e`y' if year==`y'
  predict zbhat`y' if e(sample)
  replace zbhat=zbhat`y' if year==`y'
* also get sd of mu_e for later 
  summ mu_e`y' if year==`y'
  global sd_mu_e`y'=r(sd)
* here is estimate of sd of mu 
  global sig_var`y'=(($sd_mu_e`y')^2 - ($sd_e`y')^2)
  noi disp "The estimated signal variance of mu in year `y' is = " $sig_var`y'
  noi disp "the estimated residual variance of y in year `y' is = " (${sd_e`y'})^2
}

* collapse the data to hospital-year level.
collapse (mean) mu_e zbhat (count) nobs, by(id year)
tsset id year
* initiate the prediction at 0
gen muhat=.
gen muhat_prior=0 if year==1
global var_muhat=0

* now loop over years and use shrinkage target
qui forval t=1/$nyear {
	replace muhat_prior=L.muhat if `t'!=1
	reg mu_e muhat_prior [aw=nobs] if year==`t'
	predict muhat0 if e(sample)
	predict resid0 if e(sample), resid
	summ resid0 if year==`t' [aw=nobs]
	global totvar=(r(sd))^2
	gen residvar0=${var_e`t'}/nobs if year==`t'
	summ residvar0 if year==`t' [aw=nobs]
	global avgresidvar=r(mean)
	global sigvar=$totvar - $avgresidvar
	global var_muhat1=($delta^2)*$var_muhat + ([1-($delta^2)*$var_muhat]^2)/[1-($delta^2)*$var_muhat + ($sd_resid^2)/$npat]
	noi disp "for year `t': sigvar = " %3.2f $sigvar "; (= " %3.2f 1-$var_muhat " in theory)"
	disp "avg residvar = " $avgresidvar
	replace muhat = muhat0 + resid0*($sigvar)/($sigvar + residvar0) if year==`t'
	global var_muhat = $var_muhat1
	drop muhat0 resid0 residvar0
}  

* here is another slower (but easily implemented) way using mixed

qui forval t=1/$nyear {
	clear
	local t0 = `t'-1
	use test,clear 
	keep if year==`t'
	gen muhat_prior=0
	if `t'>1 {
	  merge m:1 id using muhat_`t0'
	  keep if _merge==1 | _merge==3
 	  replace muhat_prior=muhat if muhat!=.
	  drop muhat
	}
	mixed y z x muhat_prior || id:
	predict muhat, reff
	replace muhat = muhat+_b[muhat_prior]*muhat_prior if year>1
	keep muhat id
	collapse muhat, by(id)
	save muhat_`t', replace
	noi disp "for year `t': sigvar =" %3.2f exp(_b[lns1_1_1:_cons])^2
}

* here is the slower way using melogit

qui forval t=1/$nyear {
	clear
	local t0 = `t'-1
	use test,clear 
	keep if year==`t'
	gen muhat_prior=0
	if `t'>1 {
	  merge m:1 id using muhat_`t0'
	  keep if _merge==1 | _merge==3
 	  replace muhat_prior=muhat if muhat!=.
	  drop muhat
	}
	melogit dy z x muhat_prior || id:
	predict muhat, reff
	replace muhat = muhat+_b[muhat_prior]*muhat_prior if year>1
	keep muhat id
	collapse muhat, by(id)
	save muhat_`t', replace
	noi disp "for year `t': sigvar =" %3.2f _b[/:var(_cons[id])]
}


* here is another slower (but easily implemented) way using mixed
* for comparison, estimate how cms-style 3-year measure predicts the next year

* using mixed

qui forval t=1/$nyear {
	clear
	use test,clear
	keep if year>=`t'-2 & year<=`t'
	mixed y z x || id:
	predict muhat_cms, reff
	keep muhat_cms id
	collapse muhat_cms, by(id)
	save muhat_cms_`t', replace
	noi disp `t'
}

qui forval t=1/$nyear {
	clear
	local t0 = `t'-1
	use test,clear 
	keep if year==`t'
	gen muhat_prior=0
	if `t'==1 {
		gen muhat_cms=0
	}
	if `t'>1 {
	  merge m:1 id using muhat_`t0'
	  keep if _merge==1 | _merge==3
 	  replace muhat_prior=muhat if muhat!=.
	  drop muhat
	  drop _merge
	  merge m:1 id using muhat_cms_`t0'
	  keep if _merge==1 | _merge==3
	}
	mixed y z x muhat_prior || id:
	predict muhat, reff
	replace muhat = muhat+_b[muhat_prior]*muhat_prior if year>1
	noi disp "for year `t': sigvar =" %3.2f exp(_b[lns1_1_1:_cons])^2
	mixed y z x muhat_cms || id:
	noi disp "CMS year `t': sigvar =" %3.2f exp(_b[lns1_1_1:_cons])^2
	keep muhat id
	collapse muhat, by(id)
	save muhat_`t', replace
}

* here is the slower way using melogit
qui forval t=1/$nyear {
	clear
	use test,clear 
	keep if year>=`t'-2 & year<=`t'
	melogit dy z x || id:
	predict muhat_cms, reff
	keep muhat_cms id
	collapse muhat_cms, by(id)
	save muhat_cms_`t', replace
}

qui forval t=1/$nyear {
	clear
	local t0 = `t'-1
	use test,clear 
	keep if year==`t'
	gen muhat_prior=0
	if `t'==1 {
		gen muhat_cms=0
	}
	if `t'>1 {
	  merge m:1 id using muhat_`t0'
	  keep if _merge==1 | _merge==3
 	  replace muhat_prior=muhat if muhat!=.
	  drop muhat
	  drop _merge
	  merge m:1 id using muhat_cms_`t0'
	  keep if _merge==1 | _merge==3
	}
	
	melogit dy z x muhat_prior || id:
	predict muhat, reff
	replace muhat = muhat+_b[muhat_prior]*muhat_prior if year>1
	noi disp "for year `t': sigvar =" %3.2f _b[/:var(_cons[id])]
	melogit dy z x muhat_cms || id:
	noi disp "CMS year `t': sigvar =" %3.2f _b[/:var(_cons[id])]
	keep muhat id
	collapse muhat, by(id)
	save muhat_`t', replace
}
