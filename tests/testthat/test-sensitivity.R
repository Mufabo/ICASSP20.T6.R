test_that("Evaluation", {
  data("si_sens")
  bic_final <- si_sens$bic.final
  eps_iter <- si_sens$eps.iter
  embic_iter <- si_sens$embic.iter
  K_true <- si_sens$K.true
  L_max <- si_sens$L.max
  MC <- si_sens$MC

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

  expect_equal(si_sens$p.under, p_under)
  expect_equal(si_sens$p.det, p_det)
  expect_equal(si_sens$p.over, p_over)
})
