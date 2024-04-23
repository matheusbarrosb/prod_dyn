#' Population dynamics model to estimate gross fish biomass production
#'
#' Converts population demographic rates into cohort accumulated biomass
#' @param GR Fish growth rate during the juvenile phase (mm/day)
#' @param N0 Initial number of fish in the cohort 
#' @param M_ref Reference mortality rate (yr-1)
#' @param L_ref Reference size at which M = M_ref (in mm)
#' @param Lm Size at maturity (mm)
#' @param K von Bertalanffy growth rate (yr-1)
#' @param Linf von Bertalanffy asymptotic size (mm)
#' @param t0 von Bertalanffy theoretical age at length zero (unitless)
#' @param a Intercept of length-weight relationship in the form W = aL^b
#' @param b Slope of length-weight relationship in the form W = aL^b
#' @param L0 Initial size / size at hatching (mm)
#' @param t.steps Number of time steps to simulate population dynamics (in days)
#' @param plot Logical ('TRUE' or 'FALSE')
#' @param juvGrowthType how the juveniles grow. Either 'linear' for linear growth or 'VB' for von Bertalanffy
#' @return Simulated biomass trajectories
#' @examples 
#' pred_B(GR = 0.5, M_ref = 1.6, L_ref = 131.5, Lm = 131.5, K = 0.33, Linf = 219,
#' t0 = -1.1, t.steps = 365*10, L0 = 20, N0 = 0.4, a = 0.0064, b = 2.6,
#' plot = TRUE)
#' @export

pred_B <- function(parList,
                   GR = NULL, N0 = NULL, M_ref = NULL, L_ref = NULL,
                   Lm = NULL, K = NULL,  Linf = NULL, t0 = NULL,
                   alpha = NULL, beta = NULL, L0 = NULL, t.steps = NULL,
                   plot = FALSE, juvGrowthType = "linear") {
  

  if (juvGrowthType == "linear") {
  L <- rep(NA, t.steps)
  L_prior <- rep(NA, t.steps)
  for(t in 2:length(L)) {
    
    L_prior[1] <- L0
    L_prior[t] <- L_prior[t-1]+GR
    
  }
  
  for (t in 1:length(L)) {
    
    L[t] <- ifelse(L_prior[t] < Lm,
                   L_prior[t]*1,
                   (Linf*(1-exp(- (K/365) *(t- (t0/365))))))
    
  }
  
  L_VB <- rep(NA, t.steps)
  for(t in 1:length(L_VB)) {L_VB[t] <- Linf*(1-exp(-(K/365)*(t - (t0/365))))}
  
  breakpoint <- c(which.min(abs(Lm-L_prior))) # find value closest to breakpoint Lm
  corr.point <- c(which.min(abs(L_prior[breakpoint]-L_VB)))
  
  L_corr <- rep(NA, t.steps)
  for(t in 1:length(L_corr)) {
    #apply a correction factor to make VB growth continue from the last point
    L_corr[t] <- ifelse(L_prior[t] < Lm,
                        L_prior[t],
                        L_VB[t] + (Lm - L_VB[breakpoint+1] ))
  }
  
  # estimation of mortality
  # mortality is size-dependent and is calculated based on known values for a given length
  M <- rep(NA, t.steps)
  for(t in 1:length(M)) {M[t] <- (M_ref/365)*(L_ref/L_corr[t])}
  
  # NUMBERS-AT-AGE
  N <- rep(NA, t.steps)
  for(t in 2:length(N)) {
    N[1] <- N0
    N[t] <- N[t-1]*exp(-M[t])
  }
  
  # INDIVIDUAL BIOMASS
  B <- rep(NA, t.steps)
  for(t in 1:length(B)) {B[t] <- (alpha*L_corr[t])^beta}
  
  # COHORT BIOMASS
  CB <-rep(NA, t.steps)
  for(t in 1:length(CB)) {CB[t] <- N[t] * B[t]}
  
  P_acc <- cumsum(CB)
  GP_max <- max(P_acc)
  
  # PLOTTING
  if (plot == TRUE) {
    
    plot(P_acc, xlab = "Time (days)",
         ylab = "Gross production", type = "l")
    
  }
  
  # ORGANIZING OUTPUT
  vars <- data.frame(L_VB, L_corr, N, B, CB, P_acc)
  names(vars) <- c("LVB", "Pred_size", "Nt", "B", "CB", "GP")
  output <- list(vars, GP_max)
  names(output) <- c("sim_values", "final_GP")
  
  } else if (juvGrowthType == "VB") {
    
    # SIZE USING VON BERTALANFFY
    L <- rep(NA, t.steps)
    for (t in 1:length(L)) {L[t] = Linf*(1-exp(-(K/365)*(t - (t0/365))))}
    
    # MORTALITY
    M <- rep(NA, t.steps)
    for(t in 1:length(M)) {M[t] <- (M_ref/365)*(L_ref/L[t])}
    
    # NUMBERS-AT-AGE
    N <- rep(NA, t.steps)
    for(t in 2:length(N)) {
      N[1] <- N0
      N[t] <- N[t-1]*exp(-M[t])
    }
    
    # INDIVIDUAL BIOMASS
    B <- rep(NA, t.steps)
    for(t in 1:length(B)) {B[t] <- (alpha*L[t])^beta}
    
    # COHORT BIOMASS
    CB <-rep(NA, t.steps)
    for(t in 1:length(CB)) {CB[t] <- N[t] * B[t]}
    
    P_acc <- cumsum(CB) # derived quantities
    GP_max <- max(P_acc)
    
    # PLOTTING
    if (plot == TRUE) {
      
      plot(P_acc, xlab = "Time (days)",
           ylab = "Gross production", type = "l")
      
    }
    
    # ORGANIZING OUTPUT
    vars <- data.frame(L, N, B, CB, P_acc)
    names(vars) <- c("Pred_size", "Nt", "B", "CB", "GP")
    output <- list(vars, GP_max)
    names(output) <- c("sim_values", "final_GP")
    
  } else {stop("Argument 'juvGrowthType' must be either 'linear' for assuming linear growth during the juvenile phase, or 'VB' for von-Bertalanffy growth during all life ")}
  
  return(output)
  
} # end of function
