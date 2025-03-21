{smcl}
{* *! version 0.0.1 February 2024}{...}
{vieweralsosee "vam" "help vam"}{...}
{viewerjumpto "Syntax" "hospital_ebayes##syntax"}{...}
{viewerjumpto "Description" "hospital_ebayes##description"}{...}
{viewerjumpto "Options" "hospital_ebayes##options"}{...}
{viewerjumpto "Examples" "hospital_ebayes##examples"}{...}
{title:Title}

{phang}
{bf:hospital_ebayes} {hline 2} Hospital-Level Empirical Bayes Estimation

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:hospital_ebayes} {depvar}, {opt hosp:italid(varname)} {opt y:ear(varname)} [{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt hosp:italid(varname)}}hospital identifier{p_end}
{synopt:{opt y:ear(varname)}}year identifier{p_end}

{syntab:Optional}
{synopt:{opt c:ontrols(varlist)}}additional control variables{p_end}
{synopt:{opt shrink:age_target(varlist)}}variables to control for before shrinkage{p_end}
{synopt:{opt abs:orb(varname)}}fixed effects to absorb{p_end}
{synopt:{opt tfx_resid}}hospital fixed effects residuals{p_end}
{synopt:{opt by(varname)}}estimate separately by groups{p_end}
{synopt:{opt data(string)}}data handling ("preserve", "tv", "merge tv"){p_end}
{synopt:{opt out:put(string)}}output file path prefix{p_end}
{synopt:{opt drift:limit(#)}}maximum number of lags (-1 for all){p_end}
{synopt:{opt leave:out_years(string)}}year ranges to leave out{p_end}
{synopt:{opt leave:out_vars(string)}}variable mappings for leave-out estimates{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:hospital_ebayes} estimates hospital-level empirical Bayes models using value-added methods adapted from education research. 
The command implements several key modifications for healthcare settings:

{pstd}
1. Controls for hospital volume effects{break}
2. Handles hospital-specific structure (one "classroom" per hospital-year){break}
3. Provides additional leave-out estimators and intermediate outputs{break}
4. Includes hospital-specific adjustments

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt hospitalid(varname)} specifies the hospital identifier variable.

{phang}
{opt year(varname)} specifies the year identifier variable.

{dlgtab:Optional}

{phang}
{opt controls(varlist)} specifies additional control variables to include in the model.

{phang}
{opt shrinkage_target(varlist)} specifies variables to control for before shrinkage estimation, 
typically used for volume adjustment.

{phang}
{opt absorb(varname)} specifies fixed effects to be absorbed.

{phang}
{opt data(string)} specifies data handling options:
{break}    "preserve" - preserves original data
{break}    "tv" - keeps only value-added estimates
{break}    "merge tv" - merges value-added estimates back to original data

{phang}
{opt leaveout_years(string)} specifies year ranges to exclude from estimation.

{phang}
{opt leaveout_vars(string)} specifies variable names for storing leave-out estimates.

{marker examples}{...}
{title:Examples}

{pstd}Basic model with controls{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid)}{p_end}

{pstd}Model with volume adjustment{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) shrinkage_target(log_volume)}{p_end}

{pstd}Model with leave-out years{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) leaveout_years("2010,2015") leaveout_vars(tv_early tv_late)}{p_end}

{marker authors}{...}
{title:Authors}

{pstd}Maurice Dalton{break}
Doug Staiger{break}

{pstd}
Based on {cmd:vam} by Michael Stepner.

{marker references}{...}
{title:References}

{phang}
Working paper available at: {browse "https://www.nber.org/system/files/working_papers/w31789/w31789.pdf"}

{title:Also see}

{psee}
Online: {help vam}
{p_end} 