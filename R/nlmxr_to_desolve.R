# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_03/model_base.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_03/model_base_ExpansionCov_lessVari.mlxtran"

# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_03/modelTwoNonUCART.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_03/model_base_fix_growthElim.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/lymphodepletion2/Frieb1_Elim_full_sonstraint3.mlxtran"
### Find model file*
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_03/modelTwoNonUCART.mlxtran"
# path_mxtran <- "D:/these/Pecc_test/3_Models/1_Models//alemtuzumab/alemtuzumab_covBlaste.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/lymphodepletion/kinkout_FCkout_elimkout_RC.mlxtran"
# path_mxtran <-"D:/these/Pecc_test/3_Models/1_Models//000_20_04_03/kinkout.mlxtran"
# path_mxtran <- "D:/these/Pecc_test/3_Models/1_Models//lymphodepletion2/endogeneousProlif_no_transfert.mlxtran"
# InfoModelPecc <- monolix_to_desolve(path_mxtran =  path_mxtran)
# path_mxtran = "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_24/Cmax_1stTry.mlxtran"
# path_mxtran = "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_05_18/CmaxXfoldCombined_reworked7.mlxtran"
# path_mxtran = "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_05_18/CmaxXfoldCombined_reworked6_Egress3.mlxtran"
# teest <- monolix_to_desolve(path_mxtran)
# teest$df
# teest <- monolix_to_desolve(path_mxtran)
# teest$res
# path_nlmixr <- "D:/Peccary/Exemple_demo/nlmixr_test/model4.nlmixr"
path_nlmixr <- "file:///D:/Peccary/Exemple_demo/nlmixr_test/rmodelShiny.nlmixr"
# temp <- nlmixr_to_desolve(path_nlmixr)
#' @export
nlmixr_to_desolve <- function(path_nlmixr){

  output <- list()
  ### Find model file from mxtran

  path <-  path_nlmixr %>%
    gsub(pattern = "file:///", replacement = "")


  root <- gsub( gsub(".+/", "", path), "", path)


  model <- readRDS(path)


# read dataset ------------------------------------------------------------
  # names(teest$df)
  # model$origData
  # output[["df"]] <- model$dataSav %>%
  #   # slice(1) %>%
  #   rename(identifier = id,
  #          eventidentifier = evid,
  #          amount = amt,
  #          observation = dv,
  #          observationtype = cmt) %>%
  #   mutate(administration = observationtype) %>%
  #   mutate(eventidentifier = if_else(eventidentifier == 101 , 1L, eventidentifier))
  output[["df"]] <-  model$origData %>%
    rename(time = TIME, identifier = ID,
                   eventidentifier = EVID,
                   amount = AMT,
                   observation = DV) %>%
    mutate(administration = 1, observationtype = 1)
# head(output[["df"]])
  # resultats ---------------------------------------------------------------

  thetana <-  model$theta.names[grep("^l", model$theta.names)] %>%
    gsub(pattern = "^l", replacement =   "")

  model[ ,c("ID", thetana)] %>%
    distinct() %>%
    mutate(ID = as.double(as.character(ID))) %>%
    arrange(ID)  -> temp



  t(temp) %>%
    as.data.frame() %>%
    rownames_to_column()-> temp

  names(temp) <- c(paste0("id",unname(as.vector(temp[1,]))))
  names(temp)[[1]] <- "Param"
  temp <- temp %>%
    slice(-1)

  tibble(Param =   model$saem.theta.name, Pop =  model$saem.init$theta) %>%
    mutate(Param = gsub("^l", "", Param)) %>%
    right_join(temp)  -> output[["res"]]



    # retrieve input ----------------------------------------------------------

  output[["df"]] %>%
    filter(!is.na(amount)) -> inputbase


  namescompar <- model$model$pred.only$state

  inputbase %>%
    as_tibble %>%
    filter(identifier == output[["df"]]$identifier[[1]], amount > 0) %>%
    select(-eventidentifier, - observation, -identifier) %>%
    mutate(Proto = "1", method = factor("add", levels = c("add", "mult", "rep"))) %>%
    mutate(use = T, delete = F) %>%
    rename(ADM = observationtype, value = amount) %>%
    mutate(ADM = as.character(ADM)) %>%
    left_join(tibble(Proto = as.character(1:length(namescompar)), var = namescompar)) %>%
    select(Proto, var, time,value,  method, use, delete, ADM) %>%
    mutate(time = as.character(time)) %>%
    mutate(var = factor(var, levels = namescompar)) -> output[["input"]]






  # Model analysis ----------------------------------------------------------
  # model$model
  # model$rest
  # model$rxode



  # initial values ----------------------------------------------------------
  # teest$initial_values
tibble(Cmt = namescompar, t0 = "0") -> output[["initial_values"]]
  # Model ----------------------------------------------------------

  # handle equations
  # model$rxode


model$rxode %>%
  gsub(pattern = "d/dt\\(", replacement = "d") %>%
  gsub(pattern = "\\) *((=)|(<-)) *", replacement = " <- ") %>%
  gsub(pattern = ";.+", replacement = "") %>%
  gsub(pattern = " = ", replacement = " <- ") -> output[["model"]]


  gsub(".+cmt\\(","", model$rxode) %>%
    gsub(pattern = "\\).+", replacement = "") -> outputname

  output[["model"]] <-  output[["model"]] %>%
    gsub(pattern = paste0(outputname, " *<-"), replacement = paste0(outputname, "_output <-"))

  # initial values ----------------------------------------------------------

  output[["res"]] %>%
    select(Param, Pop) %>%
    rename(Value = Pop) %>%
    mutate(Value = as.character(Value)) %>%
    mutate(Distrib = factor("logN", levels = c("logN", "Norm", "NoVar"))) %>%
    mutate(E = factor("Esti", levels = c("Esti", "Fix", "Input"))) -> output[["values"]]

  # res (all values) --------------------------------------------------------

  ### Omega
  model$omegaR %>%
    as.data.frame() -> omega

  for(a in 0:(ncol(omega)-1)){

    new <- round(omega[[a + 1]],3)
    new[0:a] <- ""
    omega[[a + 1]] <- new
  }

  omega -> output[["matrix"]]
  ## Pop and individual values

  # teest$Todisplay

  # OUTPUTS -----------------------------------------------------------------
  # get values from table
  # Isoalte output
  output[["Todisplay"]] <- tibble(Plot = 1, Todisplay = unique(c(namescompar, outputname)),
         Check = T, Point = F, Filter_of_dataset = "", YTYPE = NA_character_)


  return(output)
}

