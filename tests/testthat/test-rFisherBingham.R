test_that("rFisherBingham warns when mtop too low", {
  expect_warning(rFisherBingham(10,c(3,2,1),c(1,2,-3), mtop = 1))
})
