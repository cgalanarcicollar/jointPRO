JMBB <- function(data_long =list(), id_col= "id", outcome_col = "outcome", m, meas_col ="meas",
                 Time_col = "Time", status_col = "status", data_surv = list(), nDim = NULL,
                 chains = 3, iter = 2000, warmup = 1000, type = NULL, reg = FALSE) {  
  
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
  
  # Regularization is not allowed for `nDim = 1`
  if (nDim == 1 && reg != FALSE) {
    stop("Regularization is not available for models with nDim = 1.")
  }
  
  # `mp` is only allowed for `nDim = 1`
  if (type == "mp" && nDim != 1) {
    stop("'mp' type is only supported for nDim = 1.")
  }
  
  # Validate the `reg` argument
  valid_regs <- c(FALSE, "LASSO", "RIDGE", "HORSESHOE", "ELASTICN")
  if (!reg %in% valid_regs) {
    stop("Invalid value for 'reg'. Choose from FALSE, 'LASSO', 'RIDGE', 'HORSESHOE', or 'ELASTICN'.")
  }
  
  # Dynamically create `file_map` based on `type` and `reg`
  if (type == "p") {
    file_map <- list(
      `1` = "stan/JMBB_lp.stan",
      `3` = if (reg == FALSE) "stan/JMBB_SGRQ.stan" else paste0("stan/SGRQ_", reg, ".stan"),
      `8` = if (reg == FALSE) "stan/JMBB_SF36.stan" else paste0("stan/SF36_", reg, ".stan")
    )
  } else if (type == "mp") {
    # For `mp`, we only support `nDim = 1`
    file_map <- list(
      `1` = "stan/JMBB.stan"
    )
  }

  
  # Check if `nDim` exists in the map and assign the corresponding file
  if (as.character(nDim) %in% names(file_map)) {
    file <- file_map[[as.character(nDim)]]
  } else {
    stop("Dimension should be 1, 3, or 8")
  }
  
  
  # Validate input dimensions
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
  
  for (i in seq(length(outcome_col))) {
    # Force the column to become a plain atomic vector
    assign(paste0("y", i), unlist(data_long[[outcome_col[i]]], use.names = FALSE))
    assign(paste0("m", i), m[i])
  }
  
  

  start <-c()
  stop <-c()
  
  for(i in seq(length(unique(data_long$id)))){
    start[i] <- min(which(data_long[[id_col]] == i))
    stop[i] <- max(which(data_long[[id_col]] == i))
  }
  
  
  ###### JM #########
  glq <- statmod::gauss.quad(15, kind = "legendre")
  xk <- glq$nodes   # nodes
  wk <- glq$weights # weights
  K <- length(xk)   # K-points
  
  
  
  # Initialize data_stan with fixed elements
  data_stan <- list(
    N = nrow(data_long),
    n = nrow(data_surv),
    times = unlist(data_long[[meas_col]], use.names = FALSE), 
    ID = unlist(data_long[[id_col]], use.names = FALSE),
    start=start,
    stop=stop,
    Time = as.vector(data_surv[[Time_col]]),
    status = as.vector(data_surv[[status_col]]),
    K= K, 
    xk= xk,
    wk= wk)


  for (i in seq(length(outcome_col))) {
    data_stan[[paste0("m", i)]] <- get(paste0("m", i))
    data_stan[[paste0("y", i)]] <- get(paste0("y", i))
  }  
 

  pars <- c()
  
  for (i in seq(length(outcome_col))) {
    pars[i] <- paste0("betas", i)
  }
  for (i in seq(length(outcome_col))) {
    pars <- c(pars,paste0("phi", i))
  }
  pars <- c(pars,"Alpha","Sigma","nu","gamma")
  
  
 
  
  out <- rstan::stan(file = file,
                     data = data_stan, pars=pars,
                     #init = 0, 
                     chains = chains,
                     cores = getOption("mc.cores",chains),
                     iter = iter,
                     warmup = warmup)
  # Return Results
  return(out)
}
