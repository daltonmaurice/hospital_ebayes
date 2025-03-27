{smcl}
{* *! version 0.0.1 March 2025}{...}
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
{synopt:{opt debug}}display additional diagnostic information during estimation{p_end}
{synopt:{opt leave:out_years(string)}}relative year ranges to exclude from estimation{p_end}
{synopt:{opt leave:out_vars(string)}}names for variables that will store leave-out estimates{p_end}
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
{opt leaveout_years(string)} specifies relative year ranges to exclude from estimation. Years are specified 
relative to the target year t (current year). For example, "-2,2" excludes observations 2 years before and 
after the target year, while "-3,-1" excludes the three years prior to the target year.

{phang}
{opt leaveout_vars(string)} specifies names for variables that will store the leave-out estimates. Multiple 
variable names should be provided in order, corresponding to each leave-out period specified in leaveout_years.

{phang}
{opt debug} displays additional diagnostic information during the estimation process, including:
{break}    - Intermediate calculation steps
{break}    - Memory usage statistics
{break}    - Timing information for each estimation stage
{break}    - Variable construction details

{marker examples}{...}
{title:Examples}

{pstd}Basic model with controls{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid)}{p_end}

{pstd}Model with volume adjustment{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) shrinkage_target(log_volume)}{p_end}

{pstd}Model with leave-out periods{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) ///}{break}
{phang2}{cmd:    leaveout_years("-3,-1 1,3") leaveout_vars(tv_pre tv_post)}{p_end}

{pstd}In the above example, for each year t:{p_end}
{phang2}- tv_pre will contain estimates excluding data from t-3 to t-1 (prior 3 years){p_end}
{phang2}- tv_post will contain estimates excluding data from t+1 to t+3 (following 3 years){p_end}

{pstd}Model with multiple non-consecutive leave-out periods{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female comorbid) ///}{break}
{phang2}{cmd:    leaveout_years("-2,-1 1,1 2,3") leaveout_vars(tv_pre tv_adjacent tv_future)}{p_end}

{pstd}In this example, for each year t:{p_end}
{phang2}- tv_pre excludes years t-2 to t-1 (prior 2 years){p_end}
{phang2}- tv_adjacent excludes year t+1 (adjacent year){p_end}
{phang2}- tv_future excludes years t+2 to t+3 (2-3 years ahead){p_end}

{pstd}Model with debug output{p_end}
{phang2}{cmd:. hospital_ebayes mortality, hospitalid(providerid) year(year) controls(age female) debug}{p_end}

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