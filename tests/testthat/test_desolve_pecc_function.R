test_that("prototype deSolve pecc", {


  model <- "


  ##### Peccary Auto-Creation ###
  ##Pk.1comp.ka.Cl.Vd.difEq
  dtest <- 0
  ke <- Cl/Vd + dtest
  dGut_0 <- 50
  dGut <- - ka * Gut
  dCentral <- -ke * Central + ka * Gut
  Conc_output <- Central/ Vd
  "

totest <- deSolve_pecc(model)

finalvalue <- sum(!(totest$state == c("Gut" , "Central")))

finalvalue <- finalvalue + sum(!(totest$parameter == c("Cl", "Vd","ka")))

testthat::expect_equal( finalvalue , 0)

})



