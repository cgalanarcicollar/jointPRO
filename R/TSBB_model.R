TSBB <- function(fixed.formula = NULL, random.formula = NULL,  
                 m, data_long, time_col = "time", outcome_col = "outcome", 
                 id_col ="id", survival.formula = NULL, Time = NULL, 
                 status = NULL, data_surv = list(), Time_col = "Time", 
                 event_col = "status", nDim = NULL, type = NULL) {  
  
  if( is.null(type) == TRUE){
    type <- "p"
  }
  
  if(!type %in% c("p","mp")){
    stop("relationship between survival and longitudinal data should be p or mp")
  }
  
  # Multivariate case
  if (!is.null(nDim)) {
    if (nDim != as.integer(nDim)) {
      stop("The number of dimensions (nDim) must be an integer.")
    }
    if (nDim <= 0) {
      stop("The number of dimensions (nDim) must be positive.")
    }
  } else {
    nDim <- 1
  }
  
  # Handle case where a single fixed or random formula is applied to all dimensions
  if (nDim > 1 && length(fixed.formula) == 1) {
    fixed.formula <- rep(fixed.formula, nDim)  # Duplicate the formula for all dimensions
  }
  if (nDim > 1 && length(random.formula) == 1) {
    random.formula <- rep(random.formula, nDim)  # Duplicate the formula for all dimensions
  }
  
  # Validate input dimensions
  if (length(fixed.formula) != nDim) {
    stop("Each longitudinal outcome must have a corresponding fixed formula.")
  }
  if (length(random.formula) != nDim) {
    stop("Each longitudinal outcome must have a corresponding random formula.")
  }
  if (length(m) != 1 && length(m) != nDim) {
    stop("m must be a single value or a vector with a length equal to the number of longitudinal outcomes.")
  }
  
  # Ensure m is correctly aligned with the number of outcomes
  if (length(m) == 1) {
    m <- rep(m, nDim)  # Expand single value of m to match nDim
  }
  
  # Validate `data_long` structure
  if (!all(outcome_col %in% names(data_long))) {
    missing_cols <- setdiff(outcome_col, names(data_long))
    stop(paste("The longitudinal data must contain all specified outcome columns. Missing columns:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Validate survival data
  if (is.null(Time) || is.null(status)) {
    surv_response <- model.response(model.frame(survival.formula, data = data_surv))
    Time <- surv_response[, 1]
    status <- surv_response[, 2]
  }
  nRand <- length(Time) #number of subjects
  
  # Initialize storage for longitudinal models and predictions
  longitudinal_models <- list()
  predicted_means <- list()
  
  # Fit longitudinal models for each outcome
  for (i in seq_len(nDim)) {
    
    # Get the name of the current outcome column
    current_outcome <- outcome_col[i]
    
    # Subset data for this outcome
    outcome_data <- data_long[c(id_col, time_col, current_outcome)]  # Select relevant columns
    
    # Rename the outcome column to a standard name for consistent processing
    names(outcome_data)[names(outcome_data) == current_outcome] <- "y"
    
    
    # Generate longitudinal model matrices
    fixed.mf <- model.frame(formula = fixed.formula[[i]], data = outcome_data)
    X <- model.matrix(attr(fixed.mf, "terms"), data = fixed.mf)
    
    random.mf <- model.frame(formula = update(random.formula[[i]], ~ . - 1), data = outcome_data)
    nComp <- dim(random.mf)[2] #Number of random components, number of u,v...
    Z <- model.matrix(random.formula[[i]], data = outcome_data)
    
    # Validate longitudinal response
    y <- model.response(fixed.mf)
  #  if (length(y) %% nDim != 0) {
   #   stop("The dependent variable must have the same number of observations in all dimensions.")
  #  }
    if (!all(as.integer(y) == y) || max(y - m[i]) > 0 || min(y) < 0) {
      stop(paste("y for outcome", i, "must be integers bounded between 0 and m."))
    }
    
    # Fit longitudinal model
    BB <- try(PROreg::BBmm(fixed.formula = fixed.formula[[i]], Z = Z, nRandComp = c(rep(nRand,nComp)),
                           m = m[i],  # Use specific m for this outcome
                           data = outcome_data, maxiter = 100, show = FALSE))
    if (inherits(BB, "try-error")) {
      stop(paste("Error in BBmm model fit for longitudinal outcome", i))
    }
    
    # Store model summary
    longitudinal_models[[i]] <- try(summary(BB), silent = TRUE)
    if (inherits(longitudinal_models[[i]], "try-error")) {
      stop(paste("Error summarizing BBmm model for longitudinal outcome", i))
    }
    
    # STEP 2: Compute estimated values
    
    upt_fixed.formula <- as.formula(gsub(time_col,Time_col, deparse(fixed.formula[[i]])))
    upt_fixed.formula <- sub("^[^~]*", "", deparse(upt_fixed.formula))
    
    upt_random.formula <- as.formula(gsub(time_col, Time_col, deparse(random.formula[[i]])))
    
    X_upt= NULL
    Z_upt= NULL
    
    upt_fixed.mf <- model.frame(formula=upt_fixed.formula, data=data_surv)
    X_upt <- model.matrix(attr(upt_fixed.mf, "terms"), data=upt_fixed.mf)
    
    upt_random.mf <- model.frame(formula=update(upt_random.formula,~.-1), data=data_surv)
    Z_upt <-  model.matrix(upt_random.mf , data = data_surv )
    
    
    
    eta_y <- try(X_upt %*% BB$fixed.coef + Z_upt %*% BB$random.coef, silent = TRUE)
    
    if (inherits(eta_y, "try-error")) {
      stop(paste("Error computing eta_y for longitudinal outcome", i))
    }
    
    if(type=="mp"){
    predicted_means[[i]] <- m[i] * exp(eta_y) / (1 + exp(eta_y))  # Use specific m here
    }else{
      predicted_means[[i]] <-  exp(eta_y) / (1 + exp(eta_y))  # No use of m here
      
    }
    # Add predicted means to survival data as exY_i
    data_surv[[paste0("exY_", i)]] <- predicted_means[[i]]
  }
  

  
  # Construct survival formula dynamically
  survival_covariates <- paste0("exY_", seq_len(nDim), collapse = " + ")
  survival_formula_updated <- as.formula(paste("survival::Surv(",Time_col,",",event_col,") ~", survival_covariates))
  
  # Fit survival model
  TS <- try(survival::coxph(survival_formula_updated, data = data_surv), silent = TRUE)
  if (inherits(TS, "try-error")) {
    stop("Error in survival model fit.")
  }
  
  sum.TS <- summary(TS)
  
  # Return Results
  
  out <- list(LongitudinalModels = longitudinal_models, SurvivalModel = sum.TS)
  class(out) <- "TSBB"
  return(out)
}
