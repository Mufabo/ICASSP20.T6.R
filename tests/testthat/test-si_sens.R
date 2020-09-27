test_that("Why does it get stuck at ii_eps=9 ?", {
  library(ICASSP20.T6.R)
  library(zeallot)
  library(parallel)
  library(doParallel)
  library(foreach)

  # User Data
  MC <- 5 # number of Monte Carlo iterations
  epsilon <- 0.04 # percantage of replacement outliers
  N_k <- 50 # Number of samples per cluster

  em_bic <- matrix(c(1,1, 2,2, 2,4, 3,3, 3,4),5, 2, byrow = TRUE)
  embic_iter = nrow(em_bic)
  nu <- 3 # t
  qH <- 0.8 # Huber
  cT <- 4.685 # Tukey

  out_range <- matrix(c(-20, -20, 20, 20), 2, 2) # range of outliers
  step_eps <- 20 #10 # steps between outliers

  # Data Generation
  x <- seq(out_range[1, 1], out_range[1, 2], step_eps)
  y <- seq(out_range[2, 1], out_range[2, 2], step_eps)
  c(X, Y) %<-% pracma::meshgrid(x, y)
  eps_iter <- length(x)^2

  data <- array(0, c(N_k*3,2,eps_iter,MC))
  for(ii_eps in 1:eps_iter){
    for(ii_mc in 1:MC){

      c(data[,,ii_eps, ii_mc], labels, r, N, K_true, mu_true, S_true) %<-% data_31(N_k, 0)
      N_repl <- 1
      index_repl <- pracma::randperm(N, N_repl)
      data[index_repl, , ii_eps, ii_mc] <- c(X[ii_eps], Y[ii_eps])
    }
  }

  L_max <- 2*K_true # search range

  # Huber Parameters
  cH <- sqrt(stats::qchisq(qH, r))
  bH <- stats::pchisq(cH^2, r+2) + cH^2/r*(1-stats::pchisq(cH^2, r))
  aH <- gamma(r/2)/pi^(r/2) / ( (2*bH)^(r/2)*(gamma(r/2) - pracma::incgam(r/2, cH^2/(2*bH))) + (2*bH*cH^2*exp(-cH^2/(2*bH)))/(cH^2 - bH * r))

  g <- list(gaus = function(t) g_gaus(t, r),
            t = function(t) g_t(t, r, nu),
            huber = function(t) g_huber(t, r, list(cH, bH, aH)))

  rho <- list(gaus = function(t) rho_gaus(t, r),
              t = function(t) rho_t(t, r, nu),
              huber = function(t) rho_huber(t, r, list(cH, bH, aH)),
              tukey = function(t) rho_tukey(t, r, cT)
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, r, nu),
              huber = function(t) psi_huber(t, r, list(cH, bH)),
              tukey = function(t) psi_tukey(t, cT)
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, r, nu),
              huber = function(t) eta_huber(t, r, list(cH, bH)),
              tukey = function(t) eta_tukey(t, cT)
  )


  bic <- array(0, c(MC, eps_iter, L_max, 3, embic_iter))
  like <- array(0, c(MC, eps_iter, L_max, 3, embic_iter))
  pen <- array(0, c(MC, eps_iter, L_max, 3, embic_iter))

  for(ii_eps in 9:eps_iter){
    for(iMC in 1:MC){
      for(iEmBic in 1:embic_iter){
        for(ll in 1:L_max){
          c(mu_est, S_est, t, R) %<-% EM_RES(data[,,ii_eps, iMC], ll, g[[em_bic[iEmBic, 1]]], psi[[em_bic[iEmBic, 1]]])

          mem <- (R == apply(R, 1, max))

          c(bic[iMC, ii_eps, ll, 1, iEmBic], pen[iMC, ii_eps, ll, 1, iEmBic], like[iMC, ii_eps, ll, 1, iEmBic]) %<-% BIC_F(data[,,ii_eps,iMC], S_est, mu_est, t, mem, rho[[em_bic[iEmBic, 2]]], psi[[em_bic[iEmBic, 2]]], eta[[em_bic[iEmBic, 2]]])
          c(bic[iMC, ii_eps, ll, 2, iEmBic], pen[iMC, ii_eps, ll, 2, iEmBic], like[iMC, ii_eps, ll, 2, iEmBic]) %<-% BIC_A(S_est, t, mem, rho[[em_bic[iEmBic, 2]]], psi[[em_bic[iEmBic, 2]]], eta[[em_bic[iEmBic, 2]]])
          c(bic[iMC, ii_eps, ll, 3, iEmBic], pen[iMC, ii_eps, ll, 3, iEmBic], like[iMC, ii_eps, ll, 3, iEmBic]) %<-% BIC_S(S_est, t, mem , rho[[em_bic[iEmBic, 2]]])

        }
      }
    }

  }
})
