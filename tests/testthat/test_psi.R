#### psi_tukey ####

test_that("psi tukey", {
  a <- c(-0.4121757 ,-2.2634588 , 0.2893038,  0.1831577,  0.4861016)
  expect_equal(round(psi_tukey(a, .3), 4), c(15.5667, 341.8993, 0, 0, 0))
})
