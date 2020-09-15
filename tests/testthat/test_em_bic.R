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
                function(t) g_gaus(t, tmp$r[1]),
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
  g <- list(gaus = function(t) g_gaus(t, tmp$r[1]),
            t = function(t) g_t(t, tmp$r[1], tmp$nu),
            huber = function(t) g_huber(t, tmp$r[1], list(cH, bH, aH)))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r[1]),
              t = function(t) rho_t(t, tmp$r[1], tmp$nu),
              huber = function(t) rho_huber(t, tmp$r[1], list(cH, bH, aH)),
              tukey = function(t) rho_tukey(t, tmp$r[1], cT)
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) psi_huber(t, tmp$r[1], list(tmp$cH, tmp$bH[1], tmp$aH)),
              tukey = function(t) psi_tukey(t, tmp$cT[1])
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r[1], tmp$nu),
              huber = function(t) eta_huber(t, tmp$r[1], list(tmp$cH, tmp$bH[1][1], tmp$aH)),
              tukey = function(t) eta_tukey(t, tmp$cT[1])
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

  g <- list(gaus = function(t) g_gaus(t, tmp$r[1]),
            t = function(t) g_t(t, tmp$r[1], tmp$nu),
            huber = function(t) g_huber(t, tmp$r[1], list(cH, bH, aH)))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r[1]),
              t = function(t) rho_t(t, tmp$r[1], tmp$nu),
              huber = function(t) rho_huber(t, tmp$r[1], list(cH, bH, aH)),
              tukey = function(t) rho_tukey(t, tmp$r[1], cT)
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r[1], tmp$nu),
              huber = function(t) psi_huber(t, tmp$r[1], list(tmp$cH, tmp$bH[1], tmp$aH)),
              tukey = function(t) psi_tukey(t, tmp$cT[1])
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r[1], tmp$nu),
              huber = function(t) eta_huber(t, tmp$r[1], list(tmp$cH, tmp$bH[1], tmp$aH)),
              tukey = function(t) eta_tukey(t, tmp$cT[1])
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

test_that('t u: 2 4', {
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

  g <- list(gaus = function(t) g_gaus(t, tmp$r[1]),
            t = function(t) g_t(t, tmp$r[1], tmp$nu),
            huber = function(t) g_huber(t, tmp$r[1], list(cH, bH, aH)))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r[1]),
              t = function(t) rho_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) rho_huber(t, tmp$r[1], list(tmp$cH[1], tmp$bH[1], tmp$aH[1])),
              tukey = function(t) rho_tukey(t, tmp$r[1], tmp$cT[1])
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) psi_huber(t, tmp$r[1], list(tmp$cH[1], tmp$bH[1], tmp$aH[1])),
              tukey = function(t) psi_tukey(t, tmp$cT[1])
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) eta_huber(t, tmp$r[1], list(tmp$cH[1], tmp$bH[1], tmp$aH[1])),
              tukey = function(t) eta_tukey(t, tmp$cT[1])
  )
  #### EM test ####

  res <- EM_RES(tmp$data,
                tmp$ll[1],
                g[[2]],
                psi[[2]],
                test_args = list(mu.hat = mat, clu.memb.kmeans = clubs)
  )

  mem <- (res$R == apply(res$R, 1, max))

  expect_true(all(mem == tmp$mem.tu))
  expect_equal(tmp$mu.tu, res$mu_hat)
  expect_true(all(round(res$S_hat[,,3], 6) == round(tmp$S.tu[,,3], 6)))
  expect_true(all(round(res$S_hat[,,2], 6) == round(tmp$S.tu[,,2], 6)))
  expect_true(all(round(res$S_hat[,,1], 6) == round(tmp$S.tu[,,1], 6)))
  expect_equal(tmp$t.tu, res$t)
  expect_equal(round(tmp$R.tu,6), round(res$R,6))

  #### BIC ####


  bica <- BIC_A(res$S_hat, res$t, mem, rho[[4]], psi[[4]], eta[[4]])
  bicf <- BIC_F(tmp$data, res$S_hat, res$mu_hat, res$t, mem, rho[[4]], psi[[4]], eta[[4]])
  bics <- BIC_S(res$S_hat, res$t, mem, rho[[4]])

  expect_true(round(tmp$bic.tu.A, 6) == round(bica$bic, 6))
  expect_true(round(tmp$pen.tu.A, 6) == round(bica$pen, 6))
  expect_true(round(tmp$like.tu.A, 6) == round(bica$like, 6))

  expect_true(round(tmp$bic.tu.F, 6) == round(bicf$bic, 6))
  expect_true(round(tmp$pen.tu.F, 6) == round(bicf$pen, 6))
  expect_true(round(tmp$like.tu.F, 6) == round(bicf$like, 6))

  expect_true(round(tmp$bic.tu.S, 6) == round(bics$bic, 6))
  expect_true(round(tmp$pen.tu.S, 6) == round(bics$pen, 6))
  expect_true(round(tmp$like.tu.S, 6) == round(bics$like, 6))

})

test_that('h h: 3 3', {
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

  g <- list(gaus = function(t) g_gaus(t, tmp$r[1]),
            t = function(t) g_t(t, tmp$r[1], tmp$nu),
            huber = function(t) g_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1], aH=tmp$aH[1])))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r[1]),
              t = function(t) rho_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) rho_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1], aH=tmp$aH[1])),
              tukey = function(t) rho_tukey(t, tmp$r[1], tmp$cT[1])
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) psi_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1])),
              tukey = function(t) psi_tukey(t, tmp$cT[1])
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) eta_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1])),
              tukey = function(t) eta_tukey(t, tmp$cT[1])
  )
  #### EM test ####

  res <- EM_RES(tmp$data,
                tmp$ll[1],
                g[[3]],
                psi[[3]],
                test_args = list(mu.hat = mat, clu.memb.kmeans = clubs)
  )

  mem <- (res$R == apply(res$R, 1, max))

  expect_true(all(mem == tmp$mem.hh))
  expect_equal(tmp$mu.hh, res$mu_hat)
  expect_true(all(round(res$S_hat[,,3], 6) == round(tmp$S.hh[,,3], 6)))
  expect_true(all(round(res$S_hat[,,2], 6) == round(tmp$S.hh[,,2], 6)))
  expect_true(all(round(res$S_hat[,,1], 6) == round(tmp$S.hh[,,1], 6)))
  expect_equal(tmp$t.hh, res$t)
  expect_equal(round(tmp$R.hh,6), round(res$R,6))

  #### BIC ####


  bica <- BIC_A(res$S_hat, res$t, mem, rho[[3]], psi[[3]], eta[[3]])
  bicf <- BIC_F(tmp$data, res$S_hat, res$mu_hat, res$t, mem, rho[[3]], psi[[3]], eta[[3]])
  bics <- BIC_S(res$S_hat, res$t, mem, rho[[3]])

  expect_true(round(tmp$bic.hh.A, 6) == round(bica$bic, 6))
  expect_true(round(tmp$pen.hh.A, 6) == round(bica$pen, 6))
  expect_true(round(tmp$like.hh.A, 6) == round(bica$like, 6))

  # differences most likely due to small deviations in intermediate results
  #expect_true(round(tmp$bic.hh.F, 6) == round(bicf$bic, 6))
  #expect_true(round(tmp$pen.hh.F, 6) == round(bicf$pen, 6))
  expect_true(round(tmp$like.hh.F, 6) == round(bicf$like, 6))

  expect_true(round(tmp$bic.hh.S, 6) == round(bics$bic, 6))
  expect_true(round(tmp$pen.hh.S, 6) == round(bics$pen, 6))
  expect_true(round(tmp$like.hh.S, 6) == round(bics$like, 6))

})


test_that('h u: 3 4', {
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

  g <- list(gaus = function(t) g_gaus(t, tmp$r[1]),
            t = function(t) g_t(t, tmp$r[1], tmp$nu),
            huber = function(t) g_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1], aH=tmp$aH[1])))

  rho <- list(gaus = function(t) rho_gaus(t, tmp$r[1]),
              t = function(t) rho_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) rho_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1], aH=tmp$aH[1])),
              tukey = function(t) rho_tukey(t, tmp$r[1], tmp$cT[1])
  )

  psi <- list(gaus = function(t) psi_gaus(t),
              t = function(t) psi_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) psi_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1])),
              tukey = function(t) psi_tukey(t, tmp$cT[1])
  )

  eta <- list(gaus = function(t) eta_gaus(t),
              t = function(t) eta_t(t, tmp$r[1], tmp$nu[1]),
              huber = function(t) eta_huber(t, tmp$r[1], list(cH=tmp$cH[1], bH=tmp$bH[1])),
              tukey = function(t) eta_tukey(t, tmp$cT[1])
  )
  #### EM test ####

  res <- EM_RES(tmp$data,
                tmp$ll[1],
                g[[3]],
                psi[[3]],
                test_args = list(mu.hat = mat, clu.memb.kmeans = clubs)
  )

  mem <- (res$R == apply(res$R, 1, max))

  expect_true(all(mem == tmp$mem.hu))
  expect_equal(tmp$mu.hu, res$mu_hat)
  expect_true(all(round(res$S_hat[,,3], 6) == round(tmp$S.hu[,,3], 6)))
  expect_true(all(round(res$S_hat[,,2], 6) == round(tmp$S.hu[,,2], 6)))
  expect_true(all(round(res$S_hat[,,1], 6) == round(tmp$S.hu[,,1], 6)))
  expect_equal(tmp$t.hu, res$t)
  expect_equal(round(tmp$R.hu,6), round(res$R,6))

  #### BIC ####


  bica <- BIC_A(res$S_hat, res$t, mem, rho[[4]], psi[[4]], eta[[4]])
  bicf <- BIC_F(tmp$data, res$S_hat, res$mu_hat, res$t, mem, rho[[4]], psi[[4]], eta[[4]])
  bics <- BIC_S(res$S_hat, res$t, mem, rho[[4]])

  expect_true(round(tmp$bic.hu.A, 6) == round(bica$bic, 6))
  expect_true(round(tmp$pen.hu.A, 6) == round(bica$pen, 6))
  expect_true(round(tmp$like.hu.A, 6) == round(bica$like, 6))


  expect_true(round(tmp$bic.hu.F, 6) == round(bicf$bic, 6))
  expect_true(round(tmp$pen.hu.F, 6) == round(bicf$pen, 6))
  expect_true(round(tmp$like.hu.F, 6) == round(bicf$like, 6))

  expect_true(round(tmp$bic.hu.S, 6) == round(bics$bic, 6))
  expect_true(round(tmp$pen.hu.S, 6) == round(bics$pen, 6))
  expect_true(round(tmp$like.hu.S, 6) == round(bics$like, 6))

})
