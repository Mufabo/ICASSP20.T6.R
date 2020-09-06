test_that("Em simple",{
  library(ClusterR)
  library(ICASSP20.T6.R)
  # t
  nu <- 3
  # Huber
  qH <- .8
  # Tukey
  cT <- 4.685
  r = 2

  data("em_test_1_args")
  test1 <- EM_RES(em_test_1_args$data, 2, function(t) g_t(t,r=r,nu=nu), function(t) psi_tukey(t, cT=cT), test_args = em_test_1_args)
  data("em_test_1_out")
  expect_equal(test1$mu_hat, em_test_1_out$mu)
  expect_true(all(test1$S_hat[,,2] == em_test_1_out$S[,,2]))
  expect_true(all(test1$S_hat[,,1] == em_test_1_out$S[,,1]))
  expect_equal(test1$t, em_test_1_out$t)
  expect_equal(test1$R, em_test_1_out$R)

})
