fit_bd <-
 function (phylo, tot_time, f.lamb, f.mu, lamb_par, mu_par, f=1,
           meth = "Nelder-Mead", cst.lamb=FALSE, cst.mu=FALSE,
           expo.lamb=FALSE, expo.mu=FALSE, fix.mu=FALSE,
           cond="crown")
{
  if (!inherits(phylo, "phylo"))
      stop("object \"phylo\" is not of class \"phylo\"")

  nobs <- Ntip(phylo)

  if (fix.mu==FALSE)
  {
    init <- c(lamb_par,mu_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1+length(lamb_par)):length(init)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
# #       f.lamb.abs.scal <- function(t){abs(f.lamb(t,lamb_par))}
# #       f.lamb.par <- function(t){mapply(f.lamb.abs.scal, t)}
# #       f.mu.abs.scal <- function(t){abs(f.mu(t,mu_par))}
# #       f.mu.par <- function(t){mapply(f.mu.abs.scal,t)}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1) , lamb_par=temp$par[1:length(lamb_par)], mu_par=temp$par[(1+length(lamb_par)):length(init)])
  }

  else
  {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init)
    {
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t){abs(f.lamb(t,lamb_par))}
      f.mu.par <- function(t){abs(f.mu(t,mu_par))}
#       f.lamb.abs.scal <- function(t){abs(f.lamb(t,lamb_par))}
#       f.lamb.par <- function(t){mapply(f.lamb.abs.scal, t)}
#       f.mu.abs.scal <- function(t){abs(f.mu(t,mu_par))}
#       f.mu.par <- function(t){mapply(f.mu.abs.scal,t)}
      LH <- likelihood_bd(phylo,tot_time,f.lamb.par,f.mu.par,f,cst.lamb=cst.lamb,cst.mu=TRUE,expo.lamb=expo.lamb,cond=cond)
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    res <- list(model = "birth death", LH = -temp$value, aicc=2*temp$value+2*p+(2*p*(p+1))/(nobs-p-1),lamb_par=temp$par[1:length(lamb_par)])
  }

  return(res)
}