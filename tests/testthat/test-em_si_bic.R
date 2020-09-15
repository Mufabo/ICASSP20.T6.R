test_that("em si bic", {
  #### ####
  MC <- 5 # number of Monte Carlo iterations
  epsilon <- 0.04 # percantage of replacement outliers
  N_k <- 100 # Number of samples per cluster

  em_bic <- matrix(c(1,1, 2,2, 2,4, 3,3, 3,4),5, 2, byrow = TRUE)
  embic_iter = nrow(em_bic)
  nu <- 3 # t
  qH <- 0.8 # Huber
  cT <- 4.685 # Tukey

  tmp <- data_31(N_k, epsilon)
  data <- tmp$data
  labels_true <- tmp$labels
  r <- tmp$r
  N <- tmp$N
  K_true <- tmp$K_true
  mu_true <- tmp$mu_true
  S_true <- tmp$S_true
  L_max <- 2 * K_true # search range
  #### ####
  cH <- sqrt(stats::qchisq(qH, r))
  bH <- stats::pchisq(cH^2, r+2) + cH^2/r*(1-stats::pchisq(cH^2, r))
  aH <- gamma(r/2)/pi^(r/2) / ( (2*bH)^(r/2)*(gamma(r/2) - pracma::incgam(r/2, cH^2/(2*bH))) + (2*bH*cH^2*exp(-cH^2/(2*bH)))/(cH^2 - bH * r))

  #### ####

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

    #### ####
  iMC <- MC
  ii_embic <- embic_iter
  ll <- L_max

  #### ####
  data("em_si_bic_args")
  data("em_si_bic_outs")
  data("em_si_bic_data")
  #                 , test_args = list(mu.hat = t(em_si_bic_args$mu.Kmeans), clu.memb.kmeans = em_si_bic_args$clu.memb.kmeans)
  tmp <- EM_RES(em_si_bic_data[[1]], ll, g[[em_bic[ii_embic, 1]]], psi[[em_bic[ii_embic, 1]]]
                )

  #####
  mu_est <- tmp$mu_hat
  S_est <- tmp$S_hat
  t <- tmp$t
  R <- tmp$R
  mem <- (R == apply(R, 1, max))


  #### BIC ####
  bicf <- BIC_F(em_si_bic_data[[1]], S_est, mu_est, t, mem, rho[[em_bic[ii_embic, 2]]], psi[[em_bic[ii_embic, 2]]], eta[[em_bic[ii_embic, 2]]])
  bica <- BIC_A(S_est, t, mem, rho[[em_bic[ii_embic, 2]]], psi[[em_bic[ii_embic, 2]]], eta[[em_bic[ii_embic, 2]]])
  bics <- BIC_S(S_est, t, mem, rho[[em_bic[ii_embic, 2]]])

  bic[iMC, ll, , ii_embic] <- c(bicf$bic, bica$bic, bics$bic)
  like[iMC, ll, , ii_embic] <- c(bicf$like, bica$like, bics$like)
  pen[iMC, ll, , ii_embic] <- c(bicf$pen, bica$pen, bics$pen)
})
