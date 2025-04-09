

version 11
set matastrict on

mata:
    /**
     * Computes optimal weights for empirical Bayes estimation
     * 
     * @param M Square matrix of covariances
     * @param i Index for selecting row of covariances
     * @param c Vector of class/time indicators
     * @param weights Optional vector of weights for each observation
     * @returns Vector of optimal weights for shrinkage
     */
    real rowvector computeweights(real matrix M, real scalar i, real colvector c, | real colvector weights) {
        // Add safety checks
        if (rows(M) != cols(M)) {
            printf("Error: Non-square matrix M (%f x %f)\n", rows(M), cols(M))
            _error(3205, "Matrix must be square")
        }
        
        if (missing(M)) {
            printf("Warning: Matrix M contains missing values\n")
        }
        
        real matrix X
        real matrix L
        real matrix vcv
        real matrix Mpos

	// construct matrix A which is used to select the relevant elements of M in constructing the VCV matrix
	real matrix temp
	real matrix A
	temp=designmatrix(c)

        /* ************************************************************************  */
        /* *** Make M matrix which is off diagnol */
        /* ************************************************************************  */
        /* Base of code adapted from Doug Staiger, added 8/30/2019 */
        /* NOW fix vcv so that it is pos semi def (with block/n will always */
        /* be invertable see higham, NJ, 1988 "computing a nearest symetric */
        /* pos sem def matrix I do this by maintianing the estimates of sd */
        /* of each signal, and fixing the corr matrix so take pos semi def */
        /* part of vcv, use it to estimate corr(vcv), then */
        /* vcvpos = corr(vcv):*(sd*sd') */
        X=.
        L=.
        symeigensystem(M,X,L)
        Mpos = X*diag(L:*(L:>=0))*X'
        /* The original code just used M everywhere, which is a matrix that is fed into this */
        A = temp, J(rows(c),cols(Mpos)-cols(temp),0)
        /* use A to select elements of M and build the VCV.  The second term adjusts the diagonal */
        /* elements of the VCV matrix to account for the class-level and individual-level shocks */
        /* We want to make the underlying signal matrix */
        if (args()==4) vcv=A*Mpos*A' + diag(1:/weights)
        else vcv=A*Mpos*A'
        // phi is the vector of autocovariances, selected correctly using the matrix A.
        real rowvector phi
        phi=Mpos[i,.]*A'

        /* return the vector of weights, choose the VCV that D.Staiger */

        /* coded  to always be pos semi def */
        return    (phi*cholinv(vcv))
}


/**
 * Computes covariances and correlations between scores across time periods
 * 
 * @param scores_var Name of scores variable
 * @param weight_var Name of weights variable  
 * @param dim Number of time periods to compute
 * @param hospitalid_var Name of hospital ID variable
 * @returns Matrix of covariances, correlations, observations and standard errors
 */
real matrix compute_cov_corr(string scalar scores_var, string scalar weight_var, real scalar dim, string scalar hospitalid_var) {

    // pre-allocate matrix
    real matrix CC
    CC = J(dim,4,.)

    // Fill cov's and corr's: between time t and t+i
    real scalar i
    real scalar tstat
    for (i=1; i<=dim; i++) {
        // check that there are >=2 obs, in order to compute covariance
        stata(invtokens(("quietly count if !missing(",scores_var,",f",strofreal(i),".",scores_var,")"),""))
        if (st_numscalar("r(N)")>1) {
            stata(invtokens(("quietly corr ",scores_var," f",strofreal(i),".",scores_var," [aw=",weight_var,"+f",strofreal(i),".",weight_var,"], cov"),""))
            CC[i,1]=st_numscalar("r(cov_12)")
            CC[i,2]=CC[i,1] / ( sqrt(st_numscalar("r(Var_1)")) * sqrt(st_numscalar("r(Var_2)")) )
        }
        CC[i,3]=st_numscalar("r(N)")

        // Compute SE for covariance estimate
        if (st_numscalar("r(N)")>1) {
            stata(invtokens(("quietly reg ",scores_var," f",strofreal(i),".",scores_var," [aw=",weight_var,"+f",strofreal(i),".",weight_var,"], cluster(",hospitalid_var,")"),""))
            tstat=st_matrix("e(b)")[1,1] / sqrt( st_matrix("e(V)")[1,1] )
            CC[i,4]=abs(CC[i,1]/tstat)
        }
    }

    return (CC)
}

/**
 * Creates covariance vector from lag covariances and same-year covariance
 * 
 * @param lag_covariances Vector of covariances at different lags
 * @param cov_sameyear Covariance within same year
 * @param lagdim Optional dimension of output vector
 * @param driftlimit Optional maximum number of lags to use
 * @returns Vector of covariances
 */
real rowvector create_m(real colvector lag_covariances, real scalar cov_sameyear, | real scalar lagdim, real scalar driftlimit) {
    // Add debugging
    printf("lag_covariances dimensions: %f x %f\n", rows(lag_covariances), cols(lag_covariances))
    printf("cov_sameyear: %f\n", cov_sameyear)
    
    real rowvector m

    if (args()==2)	m=cov_sameyear,lag_covariances'
else {
    if (length(lag_covariances)<driftlimit) _error("driftlimit specified is higher than the number of lags in the dataset")
    m=cov_sameyear,lag_covariances'[1..driftlimit],J(1,lagdim-driftlimit,lag_covariances[driftlimit])
}

return (m)
}

/**
 * Checks that covariance vector contains no missing values
 * 
 * @param m Vector of covariances to check
 */
void check_m_nomissing(real rowvector m) {
    if (missing(m)>0) _error("covariance vector contains missing values")
}

/**
 * Converts vector to symmetric matrix with stripe diagonal pattern
 * 
 * @param m Input vector
 * @returns Symmetric matrix with stripe diagonal pattern
 */
real matrix vectorToStripeDiag(real vector m) {
    // Add debugging
    printf("Input vector m dimensions: %f x %f\n", rows(m), cols(m))
    
    real scalar dim
    dim = length(m)

    // pre-allocate matrix M
    real matrix M
    M=J(dim,dim,.)

    // fill lower triangle of M
    real scalar i
    real scalar j
    for (i=1; i<=dim; i++) {
        for (j=i; j<=dim; j++) {
            M[j,i]=m[j-i+1]
        }
    }

    _makesymmetric(M)
    return (M)
}

/**
 * Appends two matrices horizontally, padding with missing values if needed
 * 
 * @param A First matrix
 * @param B Second matrix to append
 * @returns Combined matrix
 */
real matrix rightAppendMatrices(real matrix A, real matrix B) {
    real scalar rA
    real scalar rB
    rA=rows(A)
    rB=rows(B)

    if (rA==rB)		return (A,B)
    else if (rA<rB)	return ( ( A \ J(rB-rA,cols(A),.) ) , B )
    else			return ( A , ( B \ J(rA-rB,cols(B),.) ) )
}

/**
 * Saves variance components to Stata dataset
 * 
 * @param cov_lag_accum Matrix of lag covariances
 * @param corr_lag_accum Matrix of lag correlations
 * @param obs_lag_accum Matrix of lag observations
 * @param cov_se_lag_accum Matrix of lag covariance standard errors
 * @param var_total_accum Vector of total variances
 * @param var_class_accum Vector of class-level variances
 * @param var_ind_accum Vector of individual-level variances
 * @param cov_sameyear_accum Vector of same-year covariances
 * @param corr_sameyear_accum Vector of same-year correlations
 * @param obs_sameyear_accum Vector of same-year observations
 * @param suffixes Vector of variable name suffixes
 */
void saveVariancesToDataset(real matrix cov_lag_accum, real matrix corr_lag_accum, real matrix obs_lag_accum, real matrix cov_se_lag_accum, real rowvector var_total_accum, real rowvector var_class_accum, real rowvector var_ind_accum, real rowvector cov_sameyear_accum, real rowvector corr_sameyear_accum, real rowvector obs_sameyear_accum, string rowvector suffixes) {

    stata("clear")

    // count number of lags, create correct number of obs, generate variable for number of lags
    real scalar n_lags
    n_lags=rows(cov_lag_accum)

    real scalar null
    null=st_addvar("int","lag")

    st_addobs(n_lags)
    stata("qui replace lag=_n")
    st_addobs(1)

    // generate output variables
    st_store(1::n_lags, st_addvar("float", "cov_lag":+suffixes), cov_lag_accum)
    st_store(1::n_lags, st_addvar("float", "corr_lag":+suffixes), corr_lag_accum)
    st_store(1::n_lags, st_addvar("float", "obs_lag":+suffixes), obs_lag_accum)
    st_store(1::n_lags, st_addvar("float", "cov_se_lag":+suffixes), cov_se_lag_accum)
    st_store(n_lags+1, st_addvar("float", "var_total":+suffixes), var_total_accum)
    st_store(n_lags+1, st_addvar("float", "var_class":+suffixes), var_class_accum)
    st_store(n_lags+1, st_addvar("float", "var_ind":+suffixes), var_ind_accum)
    st_store(n_lags+1, st_addvar("float", "cov_sameyear":+suffixes), cov_sameyear_accum)
    st_store(n_lags+1, st_addvar("float", "corr_sameyear":+suffixes), corr_sameyear_accum)
    st_store(n_lags+1, st_addvar("float", "obs_sameyear":+suffixes), obs_sameyear_accum)
}

/**
 * Calculates drift-adjusted value-added estimate
 * 
 * @param M Matrix of covariances
 * @param i Time period index
 * @param c Vector of class/time indicators
 * @param weights Vector of observation weights
 * @param scores Vector of scores
 * @returns Drift-adjusted value-added estimate
 */
real scalar driftcalc(real matrix M, real scalar i, real colvector c, real colvector weights, real colvector scores) {

    // b is the vector of weights
    real rowvector b
    b=computeweights(M, i, c, weights)
    // return the computed tv estimate -- where it basically is summing up all the
    // scores * weight - by matrix mulitplication of row and column vector
    return (b*scores)
}

/**
 * Calculates drift-adjusted value-added estimates for a list of hospitals
 * 
 * @param M Matrix of covariances
 * @param hospitalid_var Hospital ID variable name
 * @param time_var Time variable name
 * @param scores_var Scores variable name
 * @param weights_var Weights variable name
 * @param hospobs_var Hospital observations variable name
 * @param va_var Value-added output variable name
 * @param leaveout_years Optional string specifying years to leave out
 * @param leaveout_vars Optional string specifying variables for leave-out estimates
 * @param debug Optional debug flag
 */
void driftcalclist(real matrix M, string scalar hospitalid_var, string scalar time_var, 
    string scalar scores_var, string scalar weights_var, string scalar hospobs_var, 
    string scalar va_var, string scalar leaveout_years, string scalar leaveout_vars,
    | real scalar debug) {
    
    // Declare all variables upfront
    real scalar nobs, obs, hospitalid, obs_hosp, time, new_hospitalid, new_time, year_index, i
    real matrix Z, Z_hosp, Z_obs, z_quasi
    
    nobs = st_nobs()
    
    // Get variable indices
    real scalar hospitalid_var_ind, time_var_ind, hospobs_var_ind, va_var_ind
    hospitalid_var_ind = st_varindex(hospitalid_var)
    time_var_ind = st_varindex(time_var)
    hospobs_var_ind = st_varindex(hospobs_var)
    va_var_ind = st_varindex(va_var)
    
    // Create view of variables
    st_view(Z=., ., (hospitalid_var, time_var, weights_var, scores_var))
    
    // Initialize
    hospitalid = .
    time = .
    
    // Loop over observations
    for (obs=1; obs<=nobs; obs++) {
        new_hospitalid = _st_data(obs, hospitalid_var_ind)
        new_time = _st_data(obs, time_var_ind)
        
        // Only perform calculations for new hospital-year
        if (new_time != time | new_hospitalid != hospitalid) {
            time = new_time
            
            if (new_hospitalid != hospitalid) {
                hospitalid = new_hospitalid
                obs_hosp = _st_data(obs, hospobs_var_ind)
                st_subview(Z_hosp=., Z, (obs, obs+obs_hosp-1), .)
                year_index = min(Z_hosp[.,2])-1
            }
            
            // Get observations excluding current year
            Z_obs = select(Z_hosp, Z_hosp[.,2]:!=time)
            Z_obs = select(Z_obs, Z_obs[.,4]:!=.)
            
            // Compute standard VA
            if (rows(Z_obs) > 0) {
                st_store(obs, va_var_ind, 
                    driftcalc(M, time-year_index, Z_obs[.,2]:-year_index, Z_obs[.,3], Z_obs[.,4]))
            }

            // Compute leaveout estimates if specified
            if (args()>7) {
                string vector lyears, lvars
                lyears = tokens(leaveout_years)
                lvars = tokens(leaveout_vars)
                
                // Modify debug prints to check debug flag
                if (args()>9 & debug) {
                    printf("leaveout_years: %s\n", leaveout_years)
                    printf("leaveout_vars: %s\n", leaveout_vars)
                    printf("Number of rules: %f\n", length(lyears))
                }
                
                for (i=1; i<=length(lyears); i++) {
                    string scalar before, after
                    _parse_rule(lyears[i], before, after)
                    
                    // Wrap debug prints in debug flag check
                    if (args()>9 & debug) {
                        printf("Rule %f:  before=%s, after=%s\n", i, before, after)
                        printf("Initial Z_quasi rows: %f\n", rows(Z_obs))
                        
                        if (strlen(before)) {
                            printf("Filtering out years >= %f (current year %f + offset %f)\n", 
                                   time + strtoreal(before), time, strtoreal(before))
                        }
                        
                        if (strlen(after)) {
                            printf("Filtering out years <= %f (current year %f + offset %f)\n", 
                                   time + strtoreal(after), time, strtoreal(after))
                        }
                                      
                        printf("Creating variable: %s (index: %f)\n", lvars[i], st_varindex(lvars[i]))
                    }
                    
                    // Get base observations
                    z_quasi = Z_obs
                    
                    // Apply filters if valid
                    real scalar before_val, after_val
                    
                    if (strlen(before)) {
                        before_val = strtoreal(before)
                                   }
                    
                    if (strlen(after)) {
                        after_val = strtoreal(after)

                    }
                    
                    // apply the appropriate filter based on which values are present
                    if (strlen(before)==0 & strlen(after)) {
                        z_quasi = select(z_quasi, z_quasi[.,2] :> (time + after_val))
                    }
                    else if (strlen(before) & strlen(after)==0) {
                        z_quasi = select(z_quasi, z_quasi[.,2] :< (time + before_val))
                    }
                    else if (strlen(before) & strlen(after)) {
                        z_quasi = select(z_quasi, (z_quasi[.,2] :< (time + before_val)) + (z_quasi[.,2] :> (time + after_val)))
                    }
                    if (args()>9 & debug) {
                        printf("After filter: %f observations\n", rows(z_quasi))
                    }
                    
                    if (rows(z_quasi) > 0) {
                        st_store(obs, st_varindex(lvars[i]), 
                            driftcalc(M, time-year_index, z_quasi[.,2]:-year_index, z_quasi[.,3], z_quasi[.,4]))
                    } else {
                        printf("Warning: No observations left after filtering for rule %f\n", i)
                    }
                }
            }
        }
    }
}

/**
 * Helper function to parse leave-out rules into before/after components
 * 
 * @param rule String containing leave-out rule
 * @param before String to store before component
 * @param after String to store after component
 */
void _parse_rule(string scalar rule, string scalar before, string scalar after) {
    string vector parts
    parts = tokens(rule, ",")
    
    // Handle empty before value
    if (length(parts) == 3) {
        before = parts[1]
        after = parts[3]   
    }
    // Handle after value
    if (length(parts) == 2) {
        if (parts[1] == ",") {
            before = ""  
            after = parts[2]
        } 
        if (parts[2]==",") {
            before = parts[1]
            after = ""
        }
    }
}
end 

