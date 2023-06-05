test_that("rFisherBingham warns when mtop too low", {
  expect_warning(rFisherBingham(10,c(3,2,1),c(1,2,-3), mtop = 1))
})

test_that("small sphere distribution has some efficiency", {
  #suppose density is proportional to exp(-k0(m0.x-nu)^2 + k1m1.x)
  toFBparam <- function(k0, k1, m0, m1, nu){
    list(
      A = k0 * m0 %*% t(m0),
      k = 1,
      m0 = 2*nu*k0 * m0 - k1*m1
    )
  }
  FBparams <- toFBparam(0.5, 0.5, m0 = c(1,0,0), m1 = c(0, sqrt(2), sqrt(2)), nu = 0.2)
  expect_gt(attr(rFisherBingham(1000, mu = FBparams$m0, Aplus = FBparams$A), "summary")["efficiency"], 0.5)
  
  # circular mode
  FBparams <- toFBparam(0.5, 0, m0 = c(1,0,0), m1 = c(0, sqrt(2), sqrt(2)), nu = 0.2)
  expect_gt(attr(rFisherBingham(1000, mu = FBparams$m0, Aplus = FBparams$A), "summary")["efficiency"], 0.5)
  
  # high concentration, circular mode, efficiency goes to zero (consistent with S8.4 of Kent et al 2018, second last item)
  FBparams <- toFBparam(60, 0, m0 = c(1,0,0), m1 = c(0, sqrt(2), sqrt(2)), nu = 0.2)
  expect_lt(attr(rFisherBingham(1000, mu = FBparams$m0, Aplus = FBparams$A), "summary")["efficiency"], 0.5)
})
