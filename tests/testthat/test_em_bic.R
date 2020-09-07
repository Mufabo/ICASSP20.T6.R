test_that('gaus gaus: 1 1', {
  library(ICASSP20.T6.R)
  data("em_bic_test")
  tmp <- em_bic_test

  mat <- t(tmp$mu.Kmeans)
  mat[, c(1,2)] <- mat[, c(2,1)]

  clubs <- tmp$clu.memb.kmeans

  ones <- clubs == 1
  twos <- clubs  == 2

  clubs[ones] <- 2
  clubs[twos] <- 1


  #### EM test ####

  res <- EM_RES(tmp$data,
                tmp$ll,
                function(t) g_gaus(t, tmp$r),
                function(t) psi_gaus(t),
                test_args = list(mu.hat = mat, clu.memb.kmeans = clubs)
                )

  mem <- (res$R == apply(res$R, 1, max))

  expect_true(all(mem == tmp$mem.gg))
  expect_equal(tmp$mu.gg, res$mu_hat)
  expect_true(all(round(res$S_hat[,,3], 6) == round(tmp$S.gg[,,3], 6)))
  expect_true(all(round(res$S_hat[,,2], 6) == round(tmp$S.gg[,,2], 6)))
  expect_true(all(round(res$S_hat[,,1], 6) == round(tmp$S.gg[,,1], 6)))
  expect_equal(tmp$t.gg, res$t)
  expect_equal(round(tmp$R.gg,6), round(res$R,6))

  #### BIC ####
  g <- list(gaus = function(t) g_gaus(t, tmp$r),
            t = function(t) g_t(t, tmp$r, tmp$nu),
            huber = function(t) g_huber(t, tmp$r, list(cH, bH, aH)))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r),
              t = function(t) rho_t(t, tmp$r, tmp$nu),
              huber = function(t) rho_huber(t, tmp$r, list(cH, bH, aH)),
              tukey = function(t) rho_tukey(t, tmp$r, cT)
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r, tmp$nu),
              huber = function(t) psi_huber(t, tmp$r, list(tmp$cH, tmp$bH, tmp$aH)),
              tukey = function(t) psi_tukey(t, tmp$cT)
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r, tmp$nu),
              huber = function(t) eta_huber(t, tmp$r, list(tmp$cH, tmp$bH, tmp$aH)),
              tukey = function(t) eta_tukey(t, tmp$cT)
  )

  bica <- BIC_A(res$S_hat, res$t, mem, rho[[1]], psi[[1]], eta[[1]])
  bicf <- BIC_F(tmp$data, res$S_hat, res$mu_hat, res$t, mem, rho[[1]], psi[[1]], eta[[1]])
  bics <- BIC_S(res$S_hat, res$t, mem, rho[[1]])

  expect_true(round(tmp$bic.gg.A, 6) == round(bica$bic, 6))
  expect_true(round(tmp$pen.gg.A, 6) == round(bica$pen, 6))
  expect_true(round(tmp$like.gg.A, 6) == round(bica$like, 6))

  expect_true(round(tmp$bic.gg.F, 6) == round(bicf$bic, 6))
  expect_true(round(tmp$pen.gg.F, 6) == round(bicf$pen, 6))
  expect_true(round(tmp$like.gg.F, 6) == round(bicf$like, 6))

  expect_true(round(tmp$bic.gg.S, 6) == round(bics$bic, 6))
  expect_true(round(tmp$pen.gg.S, 6) == round(bics$pen, 6))
  expect_true(round(tmp$like.gg.S, 6) == round(bics$like, 6))


})

test_that('t t: 2 2', {
  library(ICASSP20.T6.R)
  data("em_bic_test")
  tmp <- em_bic_test

  mat <- t(tmp$mu.Kmeans)
  mat[, c(1,2)] <- mat[, c(2,1)]

  clubs <- tmp$clu.memb.kmeans

  ones <- clubs == 1
  twos <- clubs  == 2

  clubs[ones] <- 2
  clubs[twos] <- 1

  g <- list(gaus = function(t) g_gaus(t, tmp$r),
            t = function(t) g_t(t, tmp$r, tmp$nu),
            huber = function(t) g_huber(t, tmp$r, list(cH, bH, aH)))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r),
              t = function(t) rho_t(t, tmp$r, tmp$nu),
              huber = function(t) rho_huber(t, tmp$r, list(cH, bH, aH)),
              tukey = function(t) rho_tukey(t, tmp$r, cT)
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r, tmp$nu),
              huber = function(t) psi_huber(t, tmp$r, list(tmp$cH, tmp$bH, tmp$aH)),
              tukey = function(t) psi_tukey(t, tmp$cT)
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r, tmp$nu),
              huber = function(t) eta_huber(t, tmp$r, list(tmp$cH, tmp$bH, tmp$aH)),
              tukey = function(t) eta_tukey(t, tmp$cT)
  )
  #### EM test ####

  res <- EM_RES(tmp$data,
                tmp$ll,
                g[[2]],
                psi[[2]],
                test_args = list(mu.hat = mat, clu.memb.kmeans = clubs)
  )

  mem <- (res$R == apply(res$R, 1, max))

  expect_true(all(mem == tmp$mem.tt))
  expect_equal(tmp$mu.tt, res$mu_hat)
  expect_true(all(round(res$S_hat[,,3], 6) == round(tmp$S.tt[,,3], 6)))
  expect_true(all(round(res$S_hat[,,2], 6) == round(tmp$S.tt[,,2], 6)))
  expect_true(all(round(res$S_hat[,,1], 6) == round(tmp$S.tt[,,1], 6)))
  expect_equal(tmp$t.tt, res$t)
  expect_equal(round(tmp$R.tt,6), round(res$R,6))

  #### BIC ####


  bica <- BIC_A(res$S_hat, res$t, mem, rho[[2]], psi[[2]], eta[[2]])
  bicf <- BIC_F(tmp$data, res$S_hat, res$mu_hat, res$t, mem, rho[[2]], psi[[2]], eta[[2]])
  bics <- BIC_S(res$S_hat, res$t, mem, rho[[2]])

  expect_true(round(tmp$bic.tt.A, 6) == round(bica$bic, 6))
  expect_true(round(tmp$pen.tt.A, 6) == round(bica$pen, 6))
  expect_true(round(tmp$like.tt.A, 6) == round(bica$like, 6))

  expect_true(round(tmp$bic.tt.F, 6) == round(bicf$bic, 6))
  expect_true(round(tmp$pen.tt.F, 6) == round(bicf$pen, 6))
  expect_true(round(tmp$like.tt.F, 6) == round(bicf$like, 6))

  expect_true(round(tmp$bic.tt.S, 6) == round(bics$bic, 6))
  expect_true(round(tmp$pen.tt.S, 6) == round(bics$pen, 6))
  expect_true(round(tmp$like.tt.S, 6) == round(bics$like, 6))

})
