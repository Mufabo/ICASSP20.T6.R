test_that('gaus gaus', {
  data("em_bic_test")
  tmp <- em_bic_test

  mat <- t(tmp$mu.Kmeans)
  mat[, c(1,2)] <- mat[, c(2,1)]

  clubs <- tmp$clu.memb.kmeans

  ones <- clubs == 1
  twos <- clubs  == 2

  clubs[ones] <- 2
  clubs[twos] <- 1

  res <- EM_RES(tmp$data,
                tmp$ll,
                function(t) g_gaus(t, tmp$r),
                function(t) psi_gaus(t),
                test_args = list(mu.hat = mat, clu.memb.kmeans = clubs)
                )

  mem <- tmp$mem.gg


})
