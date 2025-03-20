/*******************************************************************************
Program: create dummy datas to test 
Author: Unknown
Date: Unknown
Purpose: Creates simulated hospital data to test hospital_ebayes.ado 

Description:
This program generates synthetic hospital data with the following structure:
- 2000 hospitals with varying patient volumes (10-210 patients per hospital)
- Hospital-specific random effects (normal distribution: mean=8000, sd=4000)
- Patient-level covariate with sd=0.5
- Binary outcome variable generated through a logistic process

Key Variables:
- hospid: Unique hospital identifier
- nobs: Number of observations (patients) per hospital
- b0: Hospital-specific random effect
- x: Patient-level covariate
- lnvol: Natural log of hospital volume
- y_linear: Linear predictor for outcome
- prob: Probability of outcome
- dum_y: Binary outcome variable (0/1)

Model Specifications:
- Logit model with hospital random effects
- Volume-outcome relationship included (coefficient: -0.1)
- Base probability of outcome ≈ 0.08 (intercept: -2)

Output:
- Saves dataset as '../test/test.dta'
*******************************************************************************/

version 18
clear all
cap log close
! mkdir log 
log using log/explore_ebayes, text replace

set seed 12345

* Parameters
global N_hospitals = 100     // Reasonable number of hospitals
global N_years = 5          // 5-year panel (2015-2019)
global N_patients_mean = 200 // Average patients per hospital-year
global base_year = 2015     // Starting year

* Create hospital-level dataset first
clear
set obs $N_hospitals

* Generate hospital IDs and fixed characteristics
gen long hospitalid = _n
gen double baseline_quality = rnormal(0, 0.5)  // Underlying hospital quality

* Expand to create panel
expand $N_years
bysort hospitalid: gen year = $base_year + _n - 1

* Generate time-varying hospital characteristics
gen double volume = round(exp(rnormal(5, 0.5)))  // Log-normal distribution for volume
replace volume = volume * (1 + 0.2*rnormal()) // Add some year-to-year variation
replace volume = max(50, volume) // Ensure minimum volume

* Generate time-varying quality (AR(1) process)
gen double quality_shock = rnormal(0, 0.3)
bysort hospitalid (year): replace quality_shock = ///
    0.7 * quality_shock[_n-1] + 0.3 * rnormal(0, 0.3) if _n > 1

* Combine fixed and time-varying quality
gen double true_quality = baseline_quality + quality_shock

* Generate varying number of patients per hospital-year
gen long n_patients = round($N_patients_mean * (1 + 0.3*rnormal()))
replace n_patients = max(30, n_patients) // Ensure minimum patients

* Expand to patient level
expand n_patients
bysort hospitalid year: gen patient_id = _n

* Generate patient characteristics
gen double age = 65 + 15*rnormal()
replace age = max(40, min(90, age))  // Bound between 40-90
gen byte female = runiform() < 0.5
gen double severity = max(0, rnormal(1, 0.3))
gen byte emergency = runiform() < 0.2

* Generate class/ward assignment (3 wards per hospital)
gen int ward = ceil(runiform()*3)
egen class = group(hospitalid ward)

* Generate outcome with proper error structure
* Component weights based on typical healthcare quality measures
gen double patient_effect = 0.7*rnormal()  // Patient-specific variation
gen double ward_effect = 0.3*rnormal()     // Ward-specific variation
bysort class year: replace ward_effect = ward_effect[1] // Same effect within ward-year

* Construct final outcome score
* Including:
* - True hospital quality (persistent + time-varying)
* - Volume effect (diminishing returns)
* - Patient characteristics
* - Ward effects
* - Random noise
gen double score = ///
    true_quality + ///                    // Hospital quality
    -0.1 * log(volume) + ///             // Volume effect (negative = better with volume)
    0.02 * (age-65) + ///                // Age effect
    0.1 * female + ///                    // Gender effect
    0.3 * severity + ///                  // Severity effect
    0.2 * emergency + ///                 // Emergency effect
    ward_effect + ///                     // Ward-specific effect
    patient_effect                        // Patient-specific noise

* Add some missing values (randomly) to test handling
replace score = . if runiform() < 0.01
replace volume = . if runiform() < 0.005

* Clean and organize final dataset
keep hospitalid year class ward volume age female severity emergency score true_quality
order hospitalid year class ward volume age female severity emergency score true_quality
sort hospitalid year class

* Save labels
label var hospitalid "Hospital ID"
label var year "Year"
label var class "Ward/Unit ID"
label var ward "Ward Number"
label var volume "Hospital Volume"
label var age "Patient Age"
label var female "Female Patient"
label var severity "Case Severity"
label var emergency "Emergency Admission"
label var score "Outcome Score"
label var true_quality "True Hospital Quality"

* Basic data checks
assert hospitalid != .
assert year != .
assert class != .
assert ward != .
misstable summarize

* Save dataset
save "test2.dta", replace

* Display summary statistics
di _n "Summary Statistics:"
sum score volume true_quality, detail

di _n "Correlation Structure:"
preserve
    collapse (mean) score true_quality (first) volume, by(hospitalid year)
    xtset hospitalid year
    di _n "Year-to-year correlation of hospital scores:"
    corr L.score score
    di _n "Correlation between true quality and observed scores:"
    corr true_quality score
restore

log close