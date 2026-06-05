# ============================================================
# Simulation-Based Predicted Probabilities — Observed Value Approach
# For logit GLM models (glm with family = binomial(link = "logit"))
#
#   The simulation steps are:
#   1. Draw S ~ MVN(coef, vcov) to capture estimation uncertainty
#   2. For each draw, compute predicted probabilities over all observations
#   3. Average across observations (OVA) -> one value per simulation draw
#   4. Summarise across draws -> mean + 95% CI
#
# Functions:
#   sim_ova_pp()       — baseline OVA: overall average Pr(Y=1)
#   sim_ova_range()    — OVA over a sequence of values for one focal variable
#   sim_ova_fd()       — first difference between two values of a focal variable
# ============================================================

library(MASS)   # for mvrnorm


# ------------------------------------------------------------
# HELPER: draw simulated coefficient vectors
# ------------------------------------------------------------
.draw_S <- function(model, nsim, seed, vcov_matrix = NULL) {
  set.seed(seed)
  vc    <- if (is.null(vcov_matrix)) vcov(model) else vcov_matrix #extract a variance-covariance matrix
  coefs <- if (inherits(model, "merMod")) fixef(model) else coef(model) #extract coefs or fixed effects if glmer()
  mvrnorm(nsim, coefs, vc)
}

# logit inverse-link
.inv_logit <- function(x) 1 / (1 + exp(-x))


# ============================================================
# 1. sim_ova_pp()
#    Baseline OVA predicted probability
#
# Arguments:
#   model   — fitted glm object (binomial logit)
#   nsim    — number of simulation draws (default 1000)
#   seed    — random seed (default 1234)
#   ci_lvl  — confidence level for interval (default 0.95)
#
# Returns a list:
#   $mean   — mean predicted probability
#   $lower  — lower CI bound
#   $upper  — upper CI bound
#   $sims   — full vector of nsim average predicted probabilities
# ============================================================

sim_ova_pp <- function(model, nsim = 1000, seed = 1234, ci_lvl = 0.95) {
  
  S <- .draw_S(model, nsim, seed) #draw from mvrnorm
  X <- model.matrix(model) # get the dataset from the model
  
  alpha <- (1 - ci_lvl) / 2
  
  # For each simulation draw s: compute Pr(Y=1|X) for every obs, then average
  pp_sims <- apply(S, 1, function(s) mean(.inv_logit(X %*% s)))
  
  list(
    mean  = mean(pp_sims),
    lower = quantile(pp_sims, alpha,       names = FALSE),
    upper = quantile(pp_sims, 1 - alpha,   names = FALSE),
    sims  = pp_sims
  )
}


# ============================================================
# 2. sim_ova_range()
#    OVA predicted probabilities over a sequence of values for one
#    focal variable, holding all other variables at their observed values.
#
# Arguments:
#   model       — fitted model
#   focal_var   — character: column name as it appears in model.matrix(model)
#                 Use colnames(model.matrix(model)) to check exact names.
#                 For log-transformed vars already in the data, use that name
#                 (e.g. "log_gdp_pcap_l"). For squared terms use "I(var^2)".
#   focal_range — numeric vector of values to set the focal variable to.
#                 If NULL, uses seq(min, max, length.out = n_vals).
#   n_vals      — number of values in auto-generated range (default 50)
#   nsim        — number of simulation draws (default 1000)
#   seed        — random seed (default 1234)
#   ci_lvl      — confidence level for interval (default 0.95)
#
# Returns a named list:
#   $result  — data.frame with columns: focal_value, mean, lower, upper
#   $sims    — matrix of dim (nsim x length(focal_range)) of average PPs
# ============================================================

sim_ova_range <- function(model,
                          focal_var,
                          focal_range = NULL,
                          n_vals      = 50,
                          nsim        = 1000,
                          seed        = 1234,
                          ci_lvl      = 0.95,
                          vcov_matrix  = NULL) {
  
  S <- .draw_S(model, nsim, seed, vcov_matrix)           # draw from mvrnorm
  X <- model.matrix(model)                               # get the dataset from the model
  re        <- ranef(model)$cow                          # random intercepts per country
  group_ids <- as.character(model@flist$cow)             # which country each row belongs to
  obs_re    <- re[group_ids, "(Intercept)"]              # vector, length = nrow(X)
  
  # Check focal variable exists in model matrix
  col_idx <- which(colnames(X) == focal_var)
  if (length(col_idx) == 0) {
    stop(
      paste0(
        "'", focal_var, "' not found in model matrix.\n",
        "Available columns:\n  ",
        paste(colnames(X), collapse = "\n  ")
      )
    )
  }
  
  # Build focal range
  if (is.null(focal_range)) {
    focal_range <- seq(
      min(X[, col_idx], na.rm = TRUE),
      max(X[, col_idx], na.rm = TRUE),
      length.out = n_vals
    )
  }
  
  n_scenarios <- length(focal_range)
  alpha       <- (1 - ci_lvl) / 2
  
  # Work with arrays for OVA
  cases <- array(rep(X, n_scenarios), dim = c(dim(X), n_scenarios))
  for (i in seq_len(n_scenarios)) {
    cases[, col_idx, i] <- focal_range[i]
  }
  
  # For each scenario and each simulation draw: average Pr(Y=1) across obs
  val <- matrix(NA, nrow = nsim, ncol = n_scenarios)
  for (i in seq_len(n_scenarios)) {
    val[, i] <- apply(S, 1, function(s) mean(.inv_logit(cases[, , i] %*% s + obs_re))) #add a random effect
  }
  
  result <- data.frame(
    focal_value = focal_range,
    mean        = apply(val, 2, mean),
    lower       = apply(val, 2, quantile, alpha,     names = FALSE),
    upper       = apply(val, 2, quantile, 1 - alpha, names = FALSE)
  )
  
  list(result = result, sims = val)
}


# ============================================================
# 3. sim_ova_fd()
#    OVA first difference: Pr(Y=1 | focal_var = val2) - Pr(Y=1 | focal_var = val1)
#    All other variables remain at their observed values.
#
# Arguments:
#   model      — fitted glm object (binomial logit)
#   focal_var  — character: column name in model.matrix(model)
#   val1       — value 1
#   val2       — value 2
#   nsim       — number of simulation draws (default 1000)
#   seed       — random seed (default 1234)
#   ci_lvl     — confidence level (default 0.95)
#
# Returns a named list:
#   $mean  — mean first difference
#   $lower — lower CI bound
#   $upper — upper CI bound
#   $sims  — vector of nsim first differences (for histograms etc.)
#   $pp1   — nsim-length vector of avg PPs at val1
#   $pp2   — nsim-length vector of avg PPs at val2
# ============================================================

sim_ova_fd <- function(model,
                       focal_var,
                       val1,
                       val2,
                       nsim   = 1000,
                       seed   = 1234,
                       ci_lvl = 0.95,
                       vcov_matrix  = NULL) {
  
  out <- sim_ova_range(
    model       = model,
    focal_var   = focal_var,
    focal_range = c(val1, val2), #put the focal values
    nsim        = nsim,
    seed        = seed,
    ci_lvl      = ci_lvl,
    vcov_matrix = vcov_matrix 
  )
  
  pp1 <- out$sims[, 1] #value 1
  pp2 <- out$sims[, 2] #value 2
  fd  <- pp2 - pp1 #calculate first differences
  
  alpha <- (1 - ci_lvl) / 2
  
  list(
    mean  = mean(fd),
    lower = quantile(fd, alpha,     names = FALSE),
    upper = quantile(fd, 1 - alpha, names = FALSE),
    sims  = fd,
    pp1   = pp1,
    pp2   = pp2
  )
}