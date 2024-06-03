################################################################################
#                       SIMULATED MODEL FUNCTIONS                              #
################################################################################

###Functions calculate simulated data value for each given parameter set, for decline, stable, and expansion models

###Parameter list
# p0 = (1000, 500, 100, 50), initial population size
# cp = (0.5, 0.1, 0.05), crash proportion
# tl = (200, 100, 50, 30), trajectory length (from ts)
# ts = (200, 100, 50, 30), trajectory start (ybp)
# ss = (200, 100, 50, 20), sample size


library(tidyr)

#==============================================================================#
#                         Expansion model function
#==============================================================================#

expan_xy <- function(p0, cp, tl, ts, yr) {
  xn <- p0;
  x1 <- p0 * cp;
  if (yr > (ts - 1)) {
    x <- x1;
  }
  tend <- ts - tl;
  r <- log(xn/x1)/tl;
  if ((yr < ts) & (yr >= tend)) {
    x  <- x1*exp(r*(ts - yr));
  }
  if (yr < tend) {
    x <- xn;
  }
  return(x);
}


#==============================================================================#
#                         Decline model function
#==============================================================================#
decline_xy <- function(p0, cp, tl, ts, yr) {
  x1 <- p0;
  xn <- p0 * cp;
  if (yr > (ts - 1)) {
    x <- x1;
  }
  tend <- ts - tl;
  r <- log(x1 - xn)/tl;
  if ((yr < ts) & (yr >= tend)) {
    x <- (x1 - xn)*exp((-r)*(ts - yr)) + xn;
  }
  if (yr < tend) {
    x <- xn;
  }
  return(x);
}



#==============================================================================#
#                Calculate simulated values across models
#==============================================================================#

###Parameters:
#mod = model
#c.prop = crash proportion
#all other parameters as above

model_plot <- function (mod, p0, c.prop, tl, ts, yr) {
  
  
  if (mod == 1) {
    x  <- decline_xy(p0, c.prop, tl, ts, yr)
  }
  else {
    if (mod == 2) {
      x  <- expan_xy(p0, c.prop, tl, ts, yr)
    }
    
    else {
      x  <- p0;
    }
  }
  return(x);
}




