# model <- "\nke <- Cl/V\nConc_output <- ka/V * exp(- t * ke)"
# model <- "
#
# ## Absorption one compartment
# dinhib_plot <- exp(- t)
#   ke <- Cl/Vd * dinhib
#
# dGut_0 <- 50<
# dGut <- - ka * Gut + dummyparam
# dCentral <- -ke * Central + ka * Gut
# Conc_output <- Central/ Vd
# "
# #
# # # # #
# deSolve_pecc(model)


# model <- "dCentral <- Central * -ke
# Conc <- Central/Vd
# "
# deSModel <- deSolve_pecc(model)
#' deSolve_pecc
#'
#' @description Analyze simplified deSolve syntax to extract all informations
#'
#' @param model The model written by the user
#' @export
#'

deSolve_pecc <- function(model){

  # replace = by <-
  model <- gsub(" = ",  "<-", model)

  # split the line
  model_split <- str_split(model, "\n")[[1]] %>%
    gsub(pattern = "#.+", replacement = "") %>% # remove commentary
    gsub(pattern = "\\\t",replacement = "") %>% # remove \t
    gsub(pattern = "^ *", replacement = "") # remove initial spaces




# Find every word ---------------------------------------------------------

  paste(model_split, collapse = "\n") %>% # recompute as a single character
    gsub(pattern = "(_output)|(_plot)", replacement = "") %>% # remove the "_output" and "_plot"
    gsub(pattern = "(\\*)|\\\n|(<-)|>|<|\\+|\\(|\\)|-|\\/|,", replacement =   "  ") %>% # replace sign by space (note:
    # not all of them cauz we need to keep "_", hence we can not really use [:punct:].
    str_split(pattern = " ")  -> charac

  # remove all the key words ! Can be completed for sure
  charac <- charac[[1]][!(charac[[1]] %in% c("", "&", "|", "log", "exp", "if_else", "pnorm", "dnorm", "rnorm", ",", "sqrt", "max", "min",  "t", "if", "<", ">", "=", "==", "<=", ">=", "{", "}", "else", "}else{", "}else"))]
  charac <- charac[ is.na(as.double(charac))] # remove the absolute values
  charac <- unique(charac)


# Find computed values ---------------------------------------------------------

  computed <- model_split[grep("<-", model_split)] %>%
    gsub(pattern = "<-.+| ", replacement = "") %>%
    gsub(pattern = "(_output)|(_plot)", replacement = "") # remove the "_output" and "_plot"

  computed <- unique(computed)

# Find compartment ---------------------------------------------------------

  ## name of compartment (or declared variable starting with d)
  comp <- model_split[grep("^ *d.+<-", model_split)] %>% # keep all line starting with "d[Comp_name] <-"
    gsub(pattern = " *<-.+", replacement = "" ) # and keep what's before the "<-"


  ## new 04/04: possibility to explicityly give initial condition

  # find initial value declaration
  stateIV <- comp[grep("_0$", comp)] # among comp, does some of them finished by "_0"?
  if(sum(grep("_0$", comp))>0) comp <- comp[-grep("_0$", comp)] # if yes, remove them from comp

 # get the tibble with initial values
  stateIV %>%
    enframe(value = "Cmt") %>% # put cmt with initial value in a tibble
    mutate(value = map_chr(Cmt, function(x){ # for each of them

      model_split[grep(x, model_split)] %>% # refind the line where initial values is decleared
        gsub(pattern = "(.+<-)| ", replacement = "") # what is after the "<-" is the value of it (character as it can be parameter)

    })) %>%
    select(-name) -> initialCond


  # remove paramter startint with d from compartment
  textafterdeclaration <-  gsub("(.+<-)|(.+=)", "", model_split) %>% paste0(collapse = "") # merge of all text after "<-"
  wordafterdeclaration <- str_split(pattern = " ",gsub(textafterdeclaration, pattern = "\\+|\\*|-|\\+|\\(|\\)|\\{|\\}|>|<|=|/", replacement = " "))[[1]]
  wordafterdeclaration <- unique(wordafterdeclaration)


  comp <- unique(comp)
  comp <- gsub("(_plot)|(_output)", "", comp)


  for(a in unique(comp)){

    #if we use the same name (with "first d") AND
    #we do not use the name without "first d" ()
    if(sum(a %in% wordafterdeclaration) > 0  & sum(gsub("^d", "", a) %in% wordafterdeclaration) == 0 ){
      comp <- comp[-which(comp == a)]
    }
  }


# Find parameter ----------------------------------------------------------

  parameter <-  charac[- which(charac %in% c(comp, computed, gsub("^d", "", comp)))] %>%
    unique()

#  Reconstitute deSolve code -----------------------------------------


  ## 1 ) output of compartment

  ### initial compartment line
  state <- gsub("^d", "", comp)


  ## line of cmt for output (at the end of the deSolve bloc, eg  "list(c(##dGut, dCentral##), c(tumor = tumor, Conc = Conc))")
  if(length(comp) >=1){
    output_model <-  comp %>%
      reduce(paste) %>%
      gsub(pattern = " ", replacement = ", ")
  }else{

    output_model <- character()

    }


  ## 2 ) output manual
nontouchedcomputed <-  model_split[grep("<-", model_split)] %>%
  gsub(pattern = "<-.+| ", replacement = "")

  output_manual <-     nontouchedcomputed[grep("(output_)|(_output)", nontouchedcomputed)]



  # output <- output %>% reduce(paste)
  #
  if(length(output_manual) == 0){

    manueloutput <- ""
    manueloutput2 <- ""
  }else{


    temp <- output_manual %>%
      map(~ paste0(.x, " = ", .x)) %>%
      reduce(paste, sep = ", ")

    manueloutput <- paste0(",c(", temp, ")")

    manueloutput2 <- output_manual

  }

  ## 3 ) to plot manual
  toplot <-  nontouchedcomputed[grep("(plot_)|(_plot)", nontouchedcomputed)]

  if(length(toplot) == 0){

    manuelplot <- ""
    manuelplot2 <- ""
  }else{


    temp <- toplot %>%
      map(~ paste0(.x, " = ", .x)) %>%
      reduce(paste, sep = ", ")

    manuelplot <- paste0(",c(", temp, ")")

    manuelplot2 <- toplot

  }

  ## 4 ) Extra manipulation (removing useless part)

  model <- gsub("(plot_)|(_plot)|(output_)|(_output)", "", model)
  output_model <- gsub("(plot_)|(_plot)|(output_)|(_output)", "", output_model)
  manuelplot <- gsub("(plot_)|(_plot)|(output_)|(_output)", "", manuelplot)
  toplot <- gsub("(plot_)|(_plot)|(output_)|(_output)", "", toplot)
  manueloutput2 <- gsub("(plot_)|(_plot)|(output_)|(_output)", "", manueloutput2)
  manueloutput <- gsub("(plot_)|(_plot)|(output_)|(_output)", "", manueloutput)


## addition 04/04/20 of initial condition, to remove
  if(nrow(initialCond) >0){
  modeltemp<-  str_split(model, pattern = "\n")[[1]]

  modeltemp[-grep(paste0("(", initialCond$Cmt,")", collapse = "|"), modeltemp)] %>%
    paste0(collapse = "\n") -> model
  }
## end addition 04/04/20


  ## 5 ) deSolve output !
  if(length(state) > 0){

    modell <- paste0("
                         function(t, state, parameters) {
                         with(as.list(c(state, parameters)), {", model,

                     "
                         list(c(", output_model ,")", manueloutput,manuelplot, ")
                         })
                         } "
    )
  }else{


    modell <- paste0("function(t = times, parameter){\nwith(data = parameter,{\n", model,
                     "\n return(tibble(time = t,", paste0(manueloutput2, collapse = ","),"))\n})\n}")


  }

#  RxODE -----------------------------------------
  ## 23/08/2021 addition of rxode syntax
  modelrxode <- model

  for(a in state){

    modelrxode <-
      gsub(paste0("\n *d",a, " *<-"), paste0("\nd/dt(",a,") <-"), modelrxode) %>%
      gsub(pattern = paste0("^ *d",a, " *<-"), replacement = paste0("d/dt(",a,") <-"))
  }



  modelrxode <- expr(RxODE({
    !!!parse_exprs(modelrxode)
  }))
  ## end addition rxode

#  Final output -----------------------------------------
  list(model =  parse_expr(modell),
       modelrxode = modelrxode,

  state = state, parameter = parameter,

  output_manual = manueloutput2,
  toplot = toplot,
  initialCond = initialCond
  )

}

#' @export
model_pre_filled <- function(model){

  mode <- suppressWarnings(deSolve_pecc(model))

  listoutput <- list()

# Parameters
  listoutput$parameters <-  paste0("#Param: name parameter\n#Value: value of the parameters (to complete)
#Distrib: either \"logN\", \"Norm\" or \"Fix\"
#Estim: either \"Esti\", \"Fix\" or \"Input\"\n
parameters <- tribble(~Param, ~Value, ~ Distrib, ~E,\n",

         paste0("\"", mode$parameter, "\",     ,\"logN\", \"Esti\"", collapse = ",\n"),
         ")"
         )
  # States
listoutput$states <-
paste0("\n#state: name of compartment\n#t0: can be a string\n
states <- tribble(~Cmt, ~t0,\n",

          paste0("\"", mode$state, "\",  \"0\"", collapse = ",\n"),
                                 ")"
  )

npar <- length(mode$parameter)
 # Matrix_eta
listoutput$matrixEta <-paste0("matrix_eta <- tribble(", paste("~", mode$parameter, collapse = "," ), ",\n",
       map(1:npar, ~ paste0(paste0(c(rep("\"0\"", .x -1), "\"0.3\"", rep("NA", npar-.x)),collapse = ","), ", #", mode$parameter[[.x]] )) %>%
         paste0(collapse = "\n"),"\n)") %>%
  gsub(pattern = "\"0.3\", #", replacement = "\"0.3\" #")

# Administration
listoutput$adm <-paste0("
events <- tribble(~Proto, ~var, ~time, ~value, ~method,~ADM, ~tlag, ~F, ~Perf, ~use,
                  \"1\", \"", mode$state[[1]]  ,"\", \"0\",   , \"add\",\"\" , F, F, \"none\", T
                              )" )


# plottable
toplot <- c(mode$stat, mode$output_manual, mode$toplot)
toplot <- toplot[toplot !=""]

listoutput$plot_table <-paste0("
plot_table  <- tribble(~Plot, ~Todisplay, ~Point, ~Filter_of_dataset, ~YTYPE, ~Check,\n",
                 paste0(" 1, \"", toplot  ,"\", F, \"\", NA, T", collapse = ",\n"),
                              "\n)" )

# info-per-plot
listoutput$plot_table_cov <- "plot_table_cov  <- tibble(Plot = 1, wrap = \"Output\", ylog = T, xlog = F)"



# Les outputs
listoutput$paramInd <- paste0("\nnsimul =  0

paramInd <- random_etas(n = nsimul,parameters , matrix_eta, returnExpr = F)

simult <- make_simulations(paramInd, model, states, events, times = c(seq(0,100,1)), Progress = F, returnExpr = F)
plot_simulations(simult, plot_table = plot_table, plot_table_cov = plot_table_cov, returnExpr = F)                              "
)


# Les outputs
toplot2 <- c(mode$toplot, "")
listoutput$output <-
paste0("\nmb_output <- tibble(output = \"",toplot2[[1]], "\", YTYPE = \"\", err_add = 0.1, err_prop = 0.3, export = T, rm = F)")

# ode_monolix
listoutput$monolix <- paste0("\node_monolix(func = deSolve_pecc(model)$model,mb_output =mb_output , depot = events,parms = parameters, y = states, outputcat = F  ) %>% cat")

return( paste0(listoutput, collapse = "\n") %>% cat)
}


#Param: name parameter
#Value: value of the parameters (to complete)
#Distrib: either "logN", "Norm" or "Fix"
#Estim: either "Esti", "Fix" or "Input"


#
# model <- "V1_wt <- V1 * (WT/17.3) ** alloWT
# CL_wt <- CL *  (WT/17.3)**((alloCL1* WT)**alloCL2)
#
#
# K12 <- Q/V1_wt
# K21 <- Q/V2
# KE <- CL_wt/V1_wt
# V2 <-  0.7 * V1_wt
#
#
# dCentral <-  K21*Periph - K12* Central - VMAX* Central/(KM+ Central) - KE* Central
# dPeriph <- - K21*Periph  + K12*Central
#
# output_Conc <- Central / V1_w"
# model <-
# "PSTAR   = PT + QT + QP
# dC   = -KDE*C
# dPT  = lambdaP*PT*(1-PSTAR/K) + kQpP*QP - kPQ*PT - gamma*KDE*PT*C
# dQT  = kPQ*PT - gamma*KDE*QT*C
# dQP  = gamma*KDE*QT*C - kQpP*QP - deltaQP*QP
# output <- PSTAR"


# model <- "
#
# d<-a+b
#
# dX <-  a * X + Y * Z
# dY <-  b * (Y - Z)
# dZ <- -X * Y + c * Y - Z
#
# Conc <- X + Y
# Conc2 <- X + X
# "
# #
# model <- "dX <- -  ke * X"
# #
# model <- "
#
#
# dX <- -  ke * X * t
# dY <- -ke * 50 * Y
#
# plot_test <- X
# plot_test2 <- X
# output_conc <- X
# output_put <- Y
#
#
# # "

#
# deSolve_pecc(model)$model %>% cat
# #
# # model <- "
# # ke <- Cl / Vd
# #
# dX <- -ke * X + k21 * Y - k12 * X
# dY <-  - k21 * Y + k12 * X
#
# output_conc <- X /Vd
# "
# output <- "Conc Conc2"


# paste0(gsub(",", "= ,", output_model), "= ")
#
# parameters <- c(a = -8/3, b = -10, c = 28)
# state      <- c(X = 1, Y = 1, Z = 1)
#
#
# parameters <- c(ke = 0.1)
# state      <- c(X = 80)
# times      <- seq(0, 100, by = 0.01)
#
# outputlol <- "Conc Conc2"
# outputlol <- c("Conc", "Conc2")
#
#
#
# out <- ode(y = state, times = times, func = eval(parse_expr(deSolve_pecc(model)[[1]])), parms = parameters)
# plot(out)
#
#
# desolvepcc <- deSolve_pecc(model)
#
#
#
# deSolve_pecc(model)["output_manual1"]
#
# ode_nonmem(func = eval(parse_expr(deSolve_pecc(model)[[1]])), parms = parameters, y = state, outputcat = F)
# ode_monolix(func = eval(parse_expr(deSolve_pecc(model)[[1]])), parms = parameters, y = state, outputcat = F) %>% cat
#
