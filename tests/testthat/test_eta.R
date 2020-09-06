library(ICASSP20.T6.R)

#### eta Huber ####

test_that("eta_huber raises an excpetion for invalid arguments" ,{
  expect_error(eta_huber(rnorm(4), 3, list(1,2,3,4,5)))
})

test_that("eta_huber computes correctly for arguments t and r", {
  eta <- eta_huber(c(0.8479573, -1.6692123, -0.8894032,  0.1656047), 3)
  expect_equal(eta, numeric(4))
})

test_that("eta_huber computes correctly for arguments t, r and qH", {
  eta <- eta_huber(4 + c(0.8479573, -1.6692123, -0.8894032,  0.1656047), 2, list(.3))
  res <- c(-0.0506  , -0.2188,   -0.1229,   -0.0685)
  expect_equal(round(eta, digits = 4), res)
})

test_that("eta_huber computes correctly for arguments t, r, cH and bH", {
  eta <- eta_huber(c(0.8479573, -1.6692123, -0.8894032,  0.1656047), 3, list(.4, .8))
  res <- c(-0.1390759, 0, 0, -3.6463192)
  expect_equal(eta, res)
})

#### eta t ####

test_that("eta_t computes correctly for arguments t, r, nu", {
  eta <- eta_t(c(-1.4456554,  0.9074437,  0.2274346, -0.8734823), .1, .3)
  res <- c(-0.1524, -0.1372, -0.7189, -0.6081)
  expect_equal(round(eta, digits = 4), res)
})
