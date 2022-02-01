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
# path_mxtran = "file:///D:/Peccary/Exemple_demo/Simeoni/closeIV.mlxtran"
# path_mxtran = "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_10_01/firstTryNewdataset2covDOSE82.mlxtran"
# path_mxtran  = "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_10_01/copypreviousvalue_NOCOMPET.mlxtran"

# path_mxtran  = "D:/these/Pecc_test/3_Models/1_Models/000_20_10_06/Ref_previous_days_timeendexpan2limiteeffalem.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_10_08/lymphocyte_with_peak/lymphodepletionalwithpeakfreeCmaxc.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_10_08/All/experimental3_noytype992.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_12_15/NKTcellsElimAllo/NKTcellstest2diffpeaks_EllimAllo.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_21_01_11_5ytype/ref.mlxtran"
# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_21_01_11_5ytype/model_rewrite_test/test5lowHostTPeak3RealEndbaselineHostT4.mlxtran"
# demo <- monolix_to_desolve(path_mxtran)
# path_model <- "file:///D:/these/presentations/code_sing.mlxtran";
# path_mxtran <- "file:///D:/Peccary_Annexe/Exemple_demo/Brent_Monolix2021_model/DOX_SDM_L.mlxtran"
#' @export
monolix_to_desolve <- function(path_mxtran){

  output <- list()

  ### Find model file from mxtran

  path <-  path_mxtran %>%
    gsub(pattern = "file:///", replacement = "")


  root <- gsub( gsub(".+/", "", path), "", path)
  lines <- readLines(path)

  files <- lines[grep("file *= *", lines)]; file <- files[2] %>%
    gsub(pattern = "(file =)|(')|( )",replacement =  "")

  backwardpath <- sum(str_split(file, pattern = "/")[[1]] == "..")
  bloc_root <- str_split(root, pattern = "/")[[1]]
  bloc_root <- bloc_root[bloc_root != ""]
  root2 <- bloc_root[1:(length(bloc_root) - backwardpath)] %>%
    paste(collapse = "/")

  path_model <-  paste(root2, gsub("\\.\\./","", file), sep  = "/")


  model <- readLines(path_model)

  ### Find dataset file from mxtran


   dataset <- files[1] %>%
    gsub(pattern = "(file *=)|(')|( )",replacement =  "")

  backwardpath <- sum(str_split(dataset, pattern = "/")[[1]] == "..")
  bloc_root <- str_split(root, pattern = "/")[[1]]
  bloc_root <- bloc_root[bloc_root != ""]
  root <- bloc_root[1:(length(bloc_root) - backwardpath)] %>%
    paste(collapse = "/")

  #patwhay
  path_dataset <-  paste(root, gsub("\\.\\./","", dataset), sep  = "/")

  # need separator
  delimiter <- gsub("(.+=)| ","",lines[grep("^delimiter *=", lines)])
  if(delimiter == "semicolon")  delimiter = ";"
  if(delimiter == "comma")  delimiter = ","
  if(delimiter == "space")  delimiter = " "
  if(delimiter == "tab")  delimiter = "\t"
  #

  df <- read.table(path_dataset, header = T, sep = delimiter, na = ".") %>%
    as_tibble

  names(df) <- gsub("ï..", "", names(df))

  lines[(grep( "\\[CONTENT\\]" , lines)+1):(grep("<MODEL>", lines) -1)] %>%
    as_tibble %>%
    mutate(name = map_chr(value, ~ gsub("(=.+)| ", "", .x))) %>%
    mutate(use = map_chr(value, ~ gsub("(.+use *=)| |(,.+)|}", "", .x))) %>%
    select(-value) %>%
    filter(name != "")-> headerTable



  ## Regressor control
  lines[grep("use *= *regressor", lines)] %>%
    enframe(name = NULL) %>%
    mutate(dataset = gsub("(=.+)| ", "", value)) %>%
    select(-value) %>%
    mutate( model = model[grep("use *= *regressor", model)] %>%
              enframe(name = NULL) %>%
              mutate(model = gsub("(=.+)| ", "", value)) %>%
              pull(model)) -> regressors


  if(nrow(regressors)>0){
    for(a in 1:nrow(regressors)){
      # print(a)
      names(df)[which(names(df) == regressors[[1]][[a]])] <- regressors[[2]][[a]]
      headerTable[[1]][[which(headerTable$name == regressors[[1]][[a]])]] <- regressors[[2]][[a]]
    }
  }

  # headerTable$name <- gsub("_reg$", "", headerTable$name)

  df <- df[ , headerTable[[1]]]
  for(a in 1:nrow(headerTable)){
    if(! headerTable[[2]][[a]] %in% c("regressor", "covariate")) names(df)[names(df) == headerTable[[1]][[a]]] <-  headerTable[[2]][[a]]
  }

  if(! "administration" %in% names(df)) df$administration <- "1"


  output[["df"]] <- df
  ### Reading daaset

  # resultats ---------------------------------------------------------------


  ### Find result folder from mmxtran to get parameter values

  pathtemp <- str_split(path, "/")[[1]]
  pathtemp <- paste0(pathtemp[-length(pathtemp)], collapse = "/")
  lines[grep("exportpath *=", lines)] %>%
    gsub(pattern = "(.+=)|'| ", replacement = "") -> folder

  try({

    pathResult <- paste0(pathtemp, "/", folder)

    readLines(paste0(pathResult, "/summary.txt" )) -> linesRes

    if(length(grep("correlation *matrix", tolower(linesRes)))>0) linesRes <-
      linesRes[1:grep("correlation *matrix", tolower(linesRes))]

    linesRes[grep("_pop", linesRes)] %>%
      as_tibble %>%
      mutate(Param = map_chr(value,~ gsub("(_pop.+)| ", "", .x))) %>%
      mutate(value = map_chr(value,~ gsub("(.+ :)| ", " ", .x))) %>%
      mutate(value = map_chr(value,~ gsub("^ *", "", .x))) %>%
      mutate(value = map_chr(value,~ gsub(" .+", "", .x)))

  }) -> resultsPop


  outputRes <- resultsPop %>%
    select(Param, value) %>%
    rename(Pop = value)


  try({

    pathResult <- paste0(pathtemp, "/", folder)

    read.table(paste0(pathResult, "/IndividualParameters/estimatedIndividualParameters.txt" ),header = T, sep =  ",") %>%
      as_tibble-> linesRes

    linesRes %>%
      select(id, ends_with("_mean")) %>%
      gather(-id, key = "key", value = "value") %>%
      spread(key = "id", value = "value", sep = "") %>%
      rename(Param = key) %>%
      mutate(Param = gsub("_mean", "", Param))-> meanRes

    names(meanRes)[-1] <- paste0(  names(meanRes)[-1], "_mean")

    outputRes <- outputRes %>%
      left_join(meanRes)
  })

  try({ linesRes %>%
      select(id, ends_with("_mode")) %>%
      gather(-id, key = "key", value = "value") %>%
      spread(key = "id", value = "value", sep = "") %>%
      rename(Param = key) %>%
      mutate(Param = gsub("_mode", "", Param))-> modeRes

    names(modeRes)[-1] <- paste0(  names(modeRes)[-1], "_mode")

    outputRes <- outputRes %>%
      left_join(modeRes)
  })





  # for each regressos
  for(a in regressors$model){
    # print(a)
    if(!a %in% outputRes$Param){
      temp <- df[, c("identifier", a)] %>% distinct
      names(temp) <- c("identifier", "cov")


      outputRes %>%
        imap_dfr(function(x,y){

          if(y == "Param"){
            new <- a
          }else  if(y == "Pop"){

            new <- temp %>% arrange(cov) %>% slice(floor(nrow(temp)/2)) %>% pull(cov)

          }else{

            id <- gsub("(^id)|(_mean)|(_mode)", "", y)

            temp %>%
              mutate(identifier = as.character(identifier)) %>%
              filter(identifier == id) %>% pull(cov) -> new

            if(length(new) != 1) new <- NA
          }
          # print(c(x, new))
          # print(y)
          # print(length(c(x, new)))
          return(c(x, new))

        })  -> outputRes
    }
  }

# handle bioav (replace parametr name by BioAv_X, as used in Peccary)
# Create a table to be used throughout the script

  linewithBioAv <- model[grep("depot\\(.+, *p *=", model)]
  BioAvtable <- tibble( cmt = character(), previous = character())
    for(a in linewithBioAv){


      nameCmt <- gsub(".+target *=", "", a)
      nameCmt <- gsub(" *,.+","", nameCmt)

        previous <- gsub(".+, *p *=", "", a)
        previous <- gsub("(\\)| |,).+", "", previous)

        BioAvtable <- BioAvtable %>%
        add_row(cmt = nameCmt, previous = previous)

    }

  if(nrow(BioAvtable) >0){

    # update the table with new name
    BioAvtable <- BioAvtable %>%
      rowid_to_column() %>%
      mutate(new = paste0("BioAv_", rowid) ) %>%
      select(-rowid)

    #update param table
    for(a in 1:nrow(BioAvtable)){
    outputRes$Param[outputRes$Param == BioAvtable$previous[a]] <- BioAvtable$new[a]
    }

  }


  outputRes -> output[["res"]]
  #### Isolate working bloc


  balisePK <- grep("PK:", model)
  baliseEquation <- grep("EQUATION:", model)
  baliseOutput <- grep("OUTPUT:", model)

  blocPK <- model[(balisePK+1):(baliseEquation-1)]
  blocEq <- model[(baliseEquation+1):(baliseOutput-1)]
  bloqOutput <- model[(baliseOutput +1):length(model)]



  ## Retrive input
  # x <- "depot(target = UCART, adm = 4, p = Fadm4)"

  # retrieve input ----------------------------------------------------------


  blocPK[grep("depot", blocPK)] %>%
    map(function(x){

      str_split(pattern = ",",  gsub("(depot)|\\(|\\)", "", x) )[[1]] %>%
        as_tibble %>%
        mutate(type = gsub("(=.+)| ", "", value)) %>%
        mutate(value = gsub("(.+=)| ", "", value)) %>%
        spread(key = type, value = value)
    }) %>%
    bind_rows() -> inputbase#

  # find average input

  # Rule: take only patient with the highest number of different input, then compute standardised sum
  # of every input in order to take a patient with medium dose
  # then use this patient as the reference
  if( ! "administration" %in% names(df)){
    df$administration <- "1"
  }

  if( ! "adm" %in% names(inputbase)){
    inputbase$adm <- "1"
  }


  df %>%
    filter(eventidentifier == 1) %>%
    distinct(identifier, administration) %>%
    group_by(identifier) %>%
    tally %>%
    left_join(

      df %>%
        filter(eventidentifier == 1) %>%
        group_by(identifier, administration) %>%
        summarise(mean = mean(amount)) %>%
        group_by(administration) %>%
        mutate(meanall = mean(mean)) %>%
        mutate(ratio = mean/meanall) %>%
        group_by(identifier) %>%
        summarise(sum = sum(ratio))

    ) -> temp

  temp %>%
    filter(n == max(temp$n)) -> temp

  temp %>%
    arrange(sum) %>%
    slice(floor(nrow(temp)/2)) %>%
    pull(identifier) -> reference

  inputbase %>%
    left_join(

      df %>%
        filter(eventidentifier == 1, identifier == reference) %>%
        select(time, administration, amount) %>%
        rename(adm = administration) %>%
        mutate(adm  = as.character(adm))
    ) %>%
    table_input -> temp

  # Handle Bioav

  for(a in BioAvtable$cmt){

    temp$F[temp$var == a][[1]] <- TRUE

  }

  temp -> output[["input"]]
  ## Model analysis

  # Model analysis ----------------------------------------------------------


  # Find names of compartiment

  blocEq[grep("ddt_", blocEq)] %>%
    map(~ gsub("(=.+)| |(ddt_)","", .x)) %>%
    reduce(c) -> compartment





  # Model ----------------------------------------------------------
   modeldesolved <- monolix_to_desolve_str_model(model )
  # handle equations

  modeldesolved$model -> output[["model"]]
  modeldesolved$initial_values -> output[["initial_values"]]


  # Parameter values ----------------------------------------------------------


  ## Values initial inpupt
  baliseParamer <- grep( "<PARAMETER>", lines)
  baliseMonolix <- grep("<MONOLIX>"     , lines)
  blocInitialValues <- lines[(baliseParamer + 1 ):(baliseMonolix - 1)]

  # handle bioav (replace parametr name by BioAv_X, as used in Peccary)
  for(a in BioAvtable$previous){

    blocInitialValues <- gsub(paste0(a, "_pop"), paste0(BioAvtable$new[BioAvtable$previous == a], "_pop"), blocInitialValues)
    blocInitialValues <- gsub(paste0("omega_", a), paste0("omega_", BioAvtable$new[BioAvtable$previous == a]), blocInitialValues)

  }

  ### theta
  blocInitialValues[grep("_pop *=", blocInitialValues)] %>%
    as_tibble() %>%
    mutate(name = map_chr(value, ~ gsub("(=.*)| |(_pop)", "", .x))) %>%
    mutate(ini = map_chr(value, ~ gsub("(.+value *=)|(,.+)", "", .x))) %>%
    mutate(method = map_chr(value, ~ gsub("(.+method *=)|(})", "", .x))) %>%
    select(-value) %>%
    mutate(E = if_else(method == "FIXED", F, T)) %>%
    rename(Param = name, Value = ini, Distrib = method) %>%
    mutate(Distrib = case_when(Distrib == "FIXED" ~ "fix",
                               T ~ "logN")) %>%  # imperfect but enough for now....
    mutate(Distrib = factor(Distrib, levels = c("logN", "Norm", "fix"))) %>%
    table_param -> initialvalues


  # Regressor handling
  initialvalues %>%
    full_join(output[["res"]][c(1:2)]) %>%
    mutate(Value = if_else(is.na(Value), Pop, Value)) %>%
    select(-Pop) %>% distinct() -> initialvalues

  initialvalues -> output[["values"]]

  # res (all values) --------------------------------------------------------


  output[["res"]] %>%
    left_join(output[["values"]] %>% select(Param ,Value)) %>%
    rename(Initial = Value) %>%
    select(Param, Initial, Pop, everything()) -> output[["res"]]

  ### Omega
  blocInitialValues[grep("omega_", blocInitialValues)] %>%
    as_tibble() %>%
    mutate(name = map_chr(value, ~ gsub("(=.*)| |(omega_)", "", .x))) %>%
    mutate(ini = map_chr(value, ~ gsub("(.+value *=)|(,.+)", "", .x))) %>%
    select(-value) -> omega

  initialvalues %>%
    select(-Value, -Distrib, -E) %>%
    left_join(omega %>% rename(Param = name)) %>%
    mutate(ini = if_else(is.na(ini),"0",ini)) -> omega
  # try update with final value
  try({

    readLines(paste0(pathResult, "/summary.txt" )) -> omegasearch

    for(a in omega$Param){

      try({
        omegasearch[grep(paste0("omega_", a), omegasearch)][[1]] %>%
          gsub(pattern = ".+: *", replacement = "") %>%
          gsub(pattern = " .+", replacement = "")  -> restemp

        omega$ini[omega$Param == a] <- restemp

      }, silent = T)
    }
  })
  diag(omega$ini) %>%
    as.data.frame() -> omegamtrix
  rownames(omegamtrix) <- omega$Param
  colnames(omegamtrix) <- omega$Param

  for(a in 2:ncol(omegamtrix)) omegamtrix[,a ][1:(a-1)] <- ""
  omegamtrix[,1] <- as.character(omegamtrix[,1])
  omegamtrix -> output[["matrix"]]
  ## Pop and individual values



  # OUTPUTS -----------------------------------------------------------------
  # get values from table
  # Isoalte output
  if(length(grep("table *=", bloqOutput)) >0){
    bloqOutput2 <- bloqOutput[-(grep("table *=", bloqOutput):length(bloqOutput))]
  }else{
    bloqOutput2 <- bloqOutput

  }

  # and get all outputs
  str_split(gsub("(.+\\{)|\\}", "", bloqOutput2) %>%
              paste(collapse = ", "), ",")[[1]] %>%
    gsub(pattern = " ", replacement = "") -> realoutput


  realoutput <- realoutput[realoutput !=""]
  # realoutput[order(realoutput)]
  str_split(lines[grep("use=observation,", lines)] %>%
              gsub(pattern = ".+name=\\{", replacement = "") %>%
              gsub(pattern = "(\\}.+)|'| ", replacement = ""), ",")[[1]] %>%
    enframe -> temp

  temp %>%
    mutate(test = c(realoutput, rep("", nrow(temp) - length(realoutput)))) %>%
    rename(YTYPE = value, Todisplay = test) %>%
    select(-name) %>%
    filter(Todisplay != "")-> temp

  # if no several YTYPE
  if(nrow(temp %>% filter(Todisplay == temp[[2]][[1]]))>1){
    temp %>%
      mutate(YTYPE = "") %>%
      distinct -> temp

  }


  # add "output" to the model

  model_temp <- str_split(output$model, pattern = "\n")[[1]]

  for(a in realoutput[!realoutput %in% compartment]){

      model_temp[grep(a, model_temp)[[1]]] <-  gsub(a, paste0(a, "_output"), model_temp[grep(a, model_temp)[[1]]] )
  }

  output$model <- model_temp %>% paste0(collapse = "\n")

  # get all outputs (table + output)

  tibble(realoutput) %>%
    distinct %>%
    rename(Todisplay = realoutput) %>%
    mutate(Check = if_else(Todisplay %in% c(output[["initial_values"]]$Cmt, realoutput), T, F)) %>%
    mutate(Point = F, Filter_of_dataset = "") %>%
    left_join(temp) %>%
    mutate(Plot = 1) %>%
    select(Plot, everything()) ->  output[["Todisplay"]]

  return(output)
}

# text <- "file:///D:/these/TMM_models.txt"
# monolix_to_desolve_str_model(text)

#' Monolix to desolve
#' @export
monolix_to_desolve_str_model <- function(text){




# Store the ouwtput
output <- list()

if(length(text) > 1)  text <- paste0(text, collapse = "\n")
# Get line per line
  model <- str_split(text,"\n")[[1]]

# If the input is the path to the model file
if(length(model) == 1 & grepl(":/", model[[1]])){

  model <- readLines(text)

}

  # Model analysis ----------------------------------------------------------



  balisePK <- grep("PK(:|>)", model)
  baliseEquation <- grep("EQUATION(:|>)", model)
  baliseOutput <- grep("OUTPUT(:|>)", model)

  blocPK <- model[(balisePK+1):(baliseEquation-1)]
  blocEq <- model[(baliseEquation+1):(baliseOutput-1)]
  bloqOutput <- model[(baliseOutput +1):length(model)]

  ### STEP 1: Find names of compartment with initial values ----------------------------------------------------------

  blocEq[grep("ddt_", blocEq)] %>%
    map(~ gsub("(=.+)| |(ddt_)","", .x)) %>%
    reduce(c) -> compartment



  # initial values ----------------------------------------------------------


  as.tibble(compartment) %>%
    mutate(initial = map_chr(value, function(x){

      temp <- blocEq[grep(paste0(" ", x, "_0"),paste0(" ", blocEq))]

      if(length(temp) == 0) return("0")

      return(gsub("(.+=)| ", "", temp))

    })) %>%
    rename(Cmt = value, t0 = initial) -> output[["initial_values"]]

  ### STEP2: convert Model with desolve syntax ----------------------------------------------------------

  # handle equations
  bloqEq2 <- blocEq
  ## replace initialisation lines
  for(a in compartment){

    bloqEq2 <- gsub(paste0("^ *", a, "_0"), paste0("d", a, "_0"), bloqEq2)

  }

  # bloqEq2 <- blocEq[-grep(paste0("(", compartment,"_0)", collapse = "|") %>% paste0("|(t0 *=)"), blocEq)]

  paste0(bloqEq2, " ") %>%
    gsub(pattern = " *, *", replacement = " , ") %>%
    gsub(pattern = "=", replacement = " = ") %>%
    gsub(pattern = " =  = ", replacement = " == ") %>%
    gsub(pattern = " ! = ", replacement = " != ") %>%
    gsub(pattern = "ddt_", replacement = "d") %>%
    gsub(pattern = " = ", replacement = " <- ") %>%
    gsub(pattern = "\\^", replacement = "**") %>%
    gsub(pattern = "^ *else ",replacement = "}else{") %>%
    gsub(pattern = "^ *elseif",replacement = "}else if(") %>%
    gsub(pattern = "^ *end ",replacement = "}") %>%
    gsub(pattern = ";",replacement = "#") -> bloqEq2


  # "quiescence= 0"

  baliseIf <- grep("^ *if( |\\()", bloqEq2)

  for( n in baliseIf)  bloqEq2[n] <-paste0(gsub("^ *if", "if(",gsub("#.+", "",bloqEq2[n])),"){",if_else(grepl("#",bloqEq2[n]), gsub(".+#", "\n#",bloqEq2[n]),""))

  baliseelsfeIf <- grep(" *else if", bloqEq2)

  for( n in baliseelsfeIf)  bloqEq2[n] <- paste0(gsub("if", "if",bloqEq2[n]),"){")


  ## Additional treatment
  bloqEq2 <- bloqEq2[!grepl("odeType *<-", bloqEq2)]
  if(grepl("<.+>", bloqEq2) %>% sum > 1) bloqEq2 <- bloqEq2[1:(grep("<.+>", bloqEq2)[[1]] - 1)]



  bloqEq2 %>% paste0(collapse = "\n")  -> output[["model"]]


  ### Step3: take if possible some values for parameter ----------------------------------------------------------

  # initial values ----------------------------------------------------------
  testparam <- grep("PARAMETER(:|>)", blocEq)
  if(length(testparam) > 0){

  temp <-   blocEq[(testparam + 1):length(blocEq)]

  if(length(grep("<.+>", temp)) > 0) temp <- temp[1:(grep("<.+>", temp)[[1]] - 1)]

  parameter <- deSolve_pecc(output[["model"]])$parameter

  table_param( left_join(
    tibble(Param = parameter),
    tibble(Param = temp) %>%
      filter(Param != "") %>%
      mutate(Value = gsub(".+= *", "", Param)) %>%
      mutate(Param = gsub(" *=.*", "", Param))) ) -> output[["values"]]



  }




return(output)


}
