test_that("Evaluation", {
  library(zeallot)
  data("si_sens")
  bic_final <- si_sens$bic.final
  eps_iter <- si_sens$eps.iter
  embic_iter <- si_sens$embic.iter
  K_true <- si_sens$K.true
  L_max <- si_sens$L.max
  MC <- si_sens$MC
  N_k <- si_sens$N.k

  #####
  p_under <- array(0, c(dim(bic_final)[4], eps_iter, embic_iter))
  p_det <- array(0, c(dim(bic_final)[4], eps_iter, embic_iter))
  p_over <- array(0, c(dim(bic_final)[4], eps_iter, embic_iter))

  K_true_det <- pracma::repmat(c(rep(K_true, K_true) == 1:K_true, numeric(L_max-K_true)), MC, 1) == 1

  K_true_under <- pracma::repmat(c(!(rep(K_true, K_true-1) == 1:(K_true-1)), numeric(L_max-K_true+1) ), MC, 1) == 1

  for(iEmBic in 1:embic_iter){
    for(iEpsilon in 1:eps_iter){
      for(k in 1:dim(bic_final)[4]){
        BICmax <- aperm(bic_final[,iEpsilon,,k,iEmBic, drop = FALSE]
                        , c(1, 3, 4, 2, 5)) == apply(aperm(bic_final[,iEpsilon,,k ,iEmBic ,drop = FALSE], c(1,3,4,2,5)), 1, max )

        p_under[k, iEpsilon, iEmBic] <- sum(BICmax[K_true_under])/MC
        p_det[k, iEpsilon, iEmBic] <- sum(BICmax[K_true_det])/MC
        p_over[k, iEpsilon, iEmBic] <- 1 - p_det[k, iEpsilon, iEmBic] - p_under[k, iEpsilon, iEmBic]
      }
    }
  }
#####
  #expect_equal(si_sens$BICmax, BICmax)
  expect_equal(si_sens$p.under, p_under)
  expect_equal(si_sens$p.det, p_det)
  expect_equal(si_sens$p.over, p_over)


# names = c("Finite", "Asymptotic", "Schwarz")
# g_names = c("Gaus", "t", "Huber", "Tukey")
# p_det_2 <- aperm(p_det, c(2,1,3))
# c(data2, labels, r, N, K_true, mu_true, S_true) %<-% data_31(N_k, 0)
# for(iEmBic in 1: embic_iter){
#   for(k_bic in 1:dim(bic_final)[4]){
#     layout(t(1:2))
#
#     plot_scatter(cbind(labels, data2), K_true, r)
#     Z <- pracma::Reshape(p_det_2[, k_bic, iEmBic], dim(X)[1], dim(X)[2])
#     graphics::contour(x, y, t(Z), add=TRUE, col = grDevices::rainbow(10))
#     graphics::title(paste("EM-", g_names[[em_bic[iEmBic, 1]]], ", BIC-", g_names[[em_bic[iEmBic, 2]]], names[k_bic]))
#
#     image(z=t(seq(0,1,.1)), col=grDevices::rainbow(10), axes=FALSE, main="Slope", cex.main=.8)
#     axis(4,cex.axis=0.8,mgp=c(0,.5,0))
#   }
# }

})
