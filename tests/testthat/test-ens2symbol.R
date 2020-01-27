test_that("Ens2symbol works", {
  symbol <- ens2symbol("ENSG00000000003")
  expect_equal(symbol$hgnc_symbol, "TSPAN6")
})
