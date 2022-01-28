# model <- "ke <- cl  /  v1
# ddepot  <-   -  ka * depot
# dcentral  <-  ka * depot  -  ke * central
# test_output <- central / v1"
# # model %>% cat
# #
# # OD_input <-   tibble(Output = c("test"),
# # Group = c("1","2"), Proto = "1", TimeSample = "c(1,5,10)", add = "0.3F", prop = "0.2F", nidgroup = 30, cov = c("cl = 3, v1 = 2", "cl = 2, v1 = 4"))
# #
# # OD_input <-   tibble(Output = c("test"),
# #                      Group = c("1","2"), Proto = "1", TimeSample = "c(1,5,10)", add = "0.3F", prop = "0.2F", nidgroup = 30, cov = c("", ""))
#
# # OD_input <-   tibble(Output = c("test"),
# #                      Group = c("1","2"), Proto = "1",
# #                      TimeSample = "c(1,5,10)",
# #                      add = "0.3F", prop = "0.2F",
# #                      nidgroup = 30,
# #                      cov = c("", ""))
#
# # #
# states <- tibble(Cmt = c("depot", "central"), t0 = c(0,0))
# # #
# events <- tibble(var = "depot", time = "c(0,40)", use = T, value = 50, method = "add", Proto = "1:2")
# # #
# parameters <- tibble(Param = c("ka", "cl", "v1"), Value = c(0.1,3,5)) %>%
#   mutate(E = c("Esti", "Fix","Esti"), Distrib = c("logN", "logN", "logN"))
# # #
# diagOmega <- tibble( ka = c("0.1", "0", "0"), cl = c("","0.1F", "0"), v1 = c("", "", "0") )
#
# output <- tibble(output = "test", YTYPE = NA, err_add = 0.1, err_prop = 0.3, export = T, rm = F)
#'@export
deSolve_to_nlmixr <- function(model, states, events, parameters, diagOmega, output){

  output$err_add <- as.double(gsub("F", "", output$err_add))
  output$err_prop <- as.double(gsub("F", "", output$err_prop))
## start reproduction
paramIniValues <- paste0("l", parameters$Param, " <- log(", parameters$Value, ")")

parameters %>%
  mutate(line = case_when(Distrib == "logN" ~ paste0(Param, "<- exp(l",Param," + eta.", Param,")"),
                          Distrib == "Norm" ~ paste0(Param, "<- exp(l",Param,") + eta.", Param),
                          Distrib == "NoVar" ~ paste0(Param, "<- exp(l",Param,")"))) %>%
  pull(line) -> lineDistrib
# parametersLines <- case_when(parameters$Distrib == "logN", "LogN", "a")
errorlines <- ""
if(output$err_prop > 0){

  errorlines <-  paste0(errorlines, "prop.err <- c(0,", output$err_prop,", 1)")
}

if(output$err_add > 0){
  errorlines <-  paste0(errorlines, "\nadd.err <- c(0,", output$err_add,", 1)")
}

namesdiagOmega <- names(diagOmega)
Omegas<- gsub("F", "", diag(as.matrix(diagOmega))) %>% as.double
OmegasIniValues <- paste0("eta.",namesdiagOmega, " ~ ", Omegas )

modelAnalysed <- deSolve_pecc(model)

modeltemp <- model
for(a in modelAnalysed$state){
  modeltemp <- gsub(paste0("d",a," *<- *"), paste0("d/dt(",a,") <- "), modeltemp)
}
modeltemp <- gsub("(_output)|(_plot)", "", modeltemp)


output %>%
  mutate(line = case_when(err_add == 0 & err_prop >0 ~ paste0(output, " ~ prop(prop.err)"),
                          err_add > 0 & err_prop == 0 ~ paste0(output, " ~ add(add.err)"),
                          err_add > 0 & err_prop >0 ~ paste0(output, " ~ prop(prop.err) + add(add.err)"),
                          T ~ "neederrormodel" )) %>%
  pull(line) -> errorLine

paste0("f <- function(){
       ini({\n",
       paste0(paramIniValues, collapse = "\n"),
       "\n",errorlines,"\n",
       paste0(OmegasIniValues, collapse = "\n"),
       "\n})
  model({\n",
       paste0(lineDistrib, collapse = "\n"),
       "\n\n",modeltemp,
       "\n\n",
       paste0(errorLine, collapse = "\n"),
       "\n})
}\n") -> fcreation

paste0(fcreation,
       "d <- read.table(\"D:/Peccary/Exemple_demo/DATA/Theoph.txt\", header = T, sep = \";\")",
       "\n\nfit.s <- nlmixr(f,d,est=\"saem\",control=saemControl(n.burn=500,n.em=500,print=50))\n
saveRDS(fit.s, file = \"D:/Peccary/Exemple_demo/nlmixr_test/rmodelShiny.nlmixr\")")
}

#
# sample()
#
# Theoph %>%
#   left_join(tibble(Subject = unique(Theoph$Subject), cov = sample(letters[1:3],12, replace = T),
#             bool = sample(0:1,12, replace = T))) %>%
#   mutate(doseCAT = case_when(Dose  < 3.5 ~ "<3.5",
#                              Dose  < 4.8 ~ ">3.5 & <4.8",
#                              T ~ ">4.8")) %>%
#   mutate(AMT = 0) %>%
#   rename(ID = Subject, TIME = Time, DV = conc) %>%
#   as_tibble-> temp
#
# temp %>%
#   mutate(MDV  = 0) %>%
#   bind_rows(
#     temp %>%
#       group_by(ID) %>%
#       slice(1) %>%
#       mutate(TIME = 0) %>%
#       mutate(AMT = Dose * 1000) %>%
#       mutate(MDV  = 1)
#
#   ) %>%
#   mutate(EVID = MDV) %>%
#   # mutate(LNDV = log(DV)) %>%
#   arrange(ID, desc(AMT), TIME) -> temp
#   write.table(temp, file = "D:/Peccary/Exemple_demo/DATA/Theoph.txt", row.names = F, quote = F, sep = ";")
#
# read.table("D:/Peccary/Exemple_demo/DATA/Theoph.txt", header = T, sep = ";")
#
# names(Theoph);names(Oral_1CPT)
# head(Oral_1CPT)
# d <- Oral_1CPT %>%
#   as.tibble
# d <- d[,names(d) != "SS"]
# d %>%
#   filter(ID %in% 1:10)
# # d <- nmDataConvert(d);
#
#  f <- function(){
#   ini({
#     lCl <- 3      #log Cl (L/hr)
#     lV <- log(90)   #log Vc (L)
#     lka <- 0.1      #log Ka (1/hr)
#     adderr <- c(0, 0.2, 1)
#     eta.Cl ~ 0.1 ## BSV Cl
#     eta.V ~ 0.1 ## BSV Vc
#     eta.ka ~ 0.1 ## BSV Ka
#   })
#   model({
#     ## First parameters are defined in terms of the initial estimates
#     ## parameter names.
#     Cl <- exp(lCl + eta.Cl)
#     V = exp(lV + eta.V)
#     ka <- exp(lka + eta.ka)
#     ## After the differential equations are defined
#     kel <- Cl / V;
#     d/dt(depot)    = -ka*depot;
#     d/dt(centr)  =  ka*depot-kel*centr;
#     ## And the concentration is then calculated
#     cp = centr / V;
#     ## Last, nlmixr is told that the plasma concentration follows
#     ## a proportional error (estimated by the parameter prop.err)
#     cp ~ add(adderr)
#   })
# }
#  path[[3]]
# path <- Sys.getenv("PATH")
# path <- c("C:/Rtools/bin", "C:/Rtools/mingw_64/bin", path)
# path <- paste(path,collapse=";")
# Sys.setenv(PATH=path)
# Sys.getenv("PATH")
# # fit <- nlmixr(model.function, rxode.dataset, est="est",control=estControl(options))
# setwd("D:/Peccary/Exemple_demo/nlmixr_test")
# fit.s <- nlmixr(f,temp,est="saem",control=saemControl(n.burn=500,n.em=500,print=50));
#
# saveRDS(fit.s, file = "model6.nlmixr")
#
# # test <- readRDS("D:/Peccary/Exemple_demo/nlmixr_test/model3.nlmixr")
#
# test <- createRun("D:/Peccary/Exemple_demo/nlmixr_test/model4.nlmixr")
# plot_pred(test)
# test$logLik
#
# folder <- createFolder("file:///D:/Peccary/Exemple_demo/nlmixr_test")
# plot_pred(folder(c(1:2, 5:6)), lowerlimit = 0.1)
# results(folder(5:6))
#
# plot_GOF(folder(c(1:2, 5:6)), ylog = T, xlog = T)%>%
#   arrange(name ) %>%
#   plot_invoke
