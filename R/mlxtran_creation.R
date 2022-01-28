# df_path <- "C:/Users/titi7/lixoft/monolix/monolix2019R2/demos/1.creating_and_using_models/1.1.libraries_of_models/data/theophylline_data.txt"
# header <- c(ID = "identifier", AMT = "amount", TIME = "time", CONC = "observation", WEIGHT = "covariate", SEX = "covariate")
# delimiter <- "tab"
# model_path <- "rezre"
# ytype = tibble(output = 4L, YTYPE = NA_integer_, err_add = "0.1",
#                err_prop = "0.3", export = TRUE, rm = FALSE) %>% filter(export == T)
# add_param = ""
# cmt <- tibble(Cmt = c("A1", "A2", "A3"),
#                  t0 = c("0", "0", "0"))
# times <- c(12)
# func <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {dA1 <- -KA*A1
#   dA2 <- KA*A1 + A3* Q/V2 -A2*(CL/V1+Q/V1)
#   dA3 <- A2* Q/V1-A3* Q/V2
#
#   test <- A2/(V1)
#   list(c(dA1, dA2, dA3),c(test_output = test_output))
#   })
# }
# outputcat = F
# parameter  <- tibble(Param = c("CL", "KA", "Q",
#  "V1", "V2"), Value = c("0.1", "0.1",
# "5", "3", "5"), Distrib = c("logN", "logN", "Norm", "logN",
# "Norm"), E = c("Esti", "Fix", "Input", "Esti", "Esti"))
#
# omega <- diag(parameter$Distrib) %>%
#   as.data.frame()
# rownames(omega) <- parameter$Param
# colnames(omega) <- parameter$Param
# omegas <- parameter$Distrib
# names(omegas) <- parameter$Param
#
# event <- tibble(Proto = "1", var = "A1", time = "0",
#                       value = 50, method = "add", use = TRUE, delete = FALSE, ADM = "",
#                       F = FALSE, tlag = FALSE, Perf = "None")
# output <- tibble(output = "test", YTYPE = NA_integer_, err_add = "0.1",
#                     err_prop = "0.3", export = TRUE, rm = FALSE)
#
# mlxtran_creation(df_path = df_path ,omega = omega, model_path, header = header, parameter = parameter, output = output, exportpath = "exportpath") %>%
#   cat
#' @export
mlxtran_creation <- function(df_path,model_path,omega, header, parameter, cmt,  mb_output, delimiter = "tab", exportpath = "exporpath", event){

  output <- mb_output
  output$YTYPE[is.na(output$YTYPE)] <- ""
print("omega gestion")
  ## initial gestion of omega
  if(names(omega)[[1]] == "Par") omega <- omega[-1]

  omegas <- diag(as.matrix(omega)) %>% as.double

  names(omegas) <- names(omega)

    ## initialisation bloc
  blocs <- list()


  ### observation columns

  obs <- paste0("y_", output$YTYPE)
  if(length(obs) == 1) obs <- "DV"

  replacementstring <-  paste0("observation, name = {", obs, "}, type =  {", paste0(rep("continuous", length(obs)), collapse = ", "),"}}")

  # if  ytype



  blocs[["datafile"]] <- paste0("<DATAFILE>

         [FILEINFO]
         file = '", df_path,"'
         delimiter = ", delimiter, "
         header = {", names(header) %>% paste(collapse =", "), "}

         [CONTENT]
         ",  imap(header[header != "drop"],~ paste0(.y, " = {use=", .x,"}")) %>%
           paste(collapse = "\n") %>%
             gsub(pattern = "observation}", replacement = replacementstring)
           )




  ### bloc model
  print("bloc model")

  ytypetemp <- paste0("y", mb_output$YTYPE) %>%
    gsub(pattern = "NA$", replacement = "")

  if(obs == "DV") ytypetemp <- "DV"
# a for additive, b for prop
 ytypetemp2 <-  mb_output %>%
   mutate(add = map2_chr(err_add, YTYPE,~ paste0(if_else(.x == 0, NA_character_, paste0("a",if_else(is.na(.y), "", as.character(.y)) ))))) %>%
   mutate(prop = map2_chr(err_prop, YTYPE,~ paste0(if_else(.x == 0, NA_character_, paste0("b",if_else(is.na(.y), "", as.character(.y)) ))))) %>%
   mutate(out = paste(add, prop, sep = ", ")) %>%
   pull(out) %>%
   paste0(collapse = ", ")


print("here")

definition <-  mb_output %>%
  mutate(YTYPE = if_else(is.na(YTYPE), "", as.character(YTYPE))) %>%
  # to add additive and prop !
  mutate(error_model = case_when(err_add > 0 & err_prop >0 ~ paste0("combined1(a", YTYPE,", b",YTYPE,")"))) %>%
  mutate(last = paste0("y", YTYPE, " ={distribution=normal, prediction = ", output, ", errorModel = ", error_model,"}")) %>%
  pull(last)

if(obs == "DV") definition <- gsub("^y", "DV", definition)

definitionparamater <-
  parameter %>%
  filter(E != "Input") %>%
  left_join(tibble(Param = names(omegas), omegas = omegas)) %>%
  distinct() %>%
  transmute(line = pmap(list(Param, Distrib,omegas), function(x,y,z){
    paste0(x, " = {distribution=", if_else(y == "logN", "logNormal", "normal"), ", typical=", x,"_pop, ",if_else(z == 0 | is.na(z), "no-variability",paste0("sd=omega_",x)),"}")
  }  )) %>%
  distinct() %>%
  pull(line) %>% paste(collapse = "\n")

  blocs[["model"]] <-  paste0("<MODEL>

         [INDIVIDUAL]
         input = {", c( paste0(parameter$Param[parameter$Param != "input"], "_pop"),  paste0("omega_",names(omegas)[omegas != 0]) ) %>% paste(collapse = ", "), "}

         DEFINITION:
         ", definitionparamater  ,"

         [LONGITUDINAL]
         input = {",ytypetemp2 ,"}

         file ='", model_path ,"'

         DEFINITION:
         ",definition ,"

         <FIT>
data = ", obs, "
model = ", ytypetemp )


# Bloc parameter ----------------------------------------------------------
## rajouter error et autre
  blocs[["parameter"]] <-  paste0("

         <PARAMETER>
         ", parameter %>%
            filter(E != "Input") %>%
           transmute(out = pmap_chr(list(Param, Value,E), function(x,y,z){
             paste0(x, "_pop = {value = ", y, ", method = ",if_else(z== "Esti", "MLE", "FIXED") , "}")
           })) %>%
           pull(out) %>% paste(collapse = "\n"),
           omegas[omegas != 0] %>%
             imap(~paste0("omega_", .y, " = {value = ", .x, ", method = MLE}")) %>%
              paste(collapse = "\n"),"\n",
           paste0("a",unique(output$YTYPE)," = {value=0.01, method=MLE}") %>% paste(collapse = "\n"), "\n",
           paste0("b",unique(output$YTYPE)," = {value=0.3, method=MLE}")  %>% paste(collapse = "\n")
           )


  blocs[["task"]] <-  paste0("<MONOLIX>

  [TASKS]
  populationParameters()
  individualParameters(method = {conditionalMean, conditionalMode })
  fim(run = false,method = StochasticApproximation)
  logLikelihood(run = false,method = ImportanceSampling)
  plotResult(method = {outputplot, indfits, obspred, vpc, npc, residualsscatter, residualsdistribution, parameterdistribution, covariatemodeldiagnosis, randomeffects, covariancemodeldiagnosis, blq, predictiondistribution, likelihoodcontribution, categorizedoutput, saemresults, condmeanresults, fisher, likelihoodresults })

  [SETTINGS]
  GLOBAL:
    exportpath = ", exportpath)


  blocs %>%
    map(~ gsub(x = .x,"\n *", "\n")) %>%
    paste(collapse = "\n\n")

}
# shinyApp(pecc_ui_reunif, pecc_server_reunif)
# shinyApp(pecc_ui, pecc_server)
# mlxtran_creation(df_path = df_path ,omega = omega, model_path = model_path, header = header, parameter = parameter, output = output) %>%
  # cat

