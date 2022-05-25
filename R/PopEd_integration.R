# model <- "ke <- cl  /  v1
# ddepot  <-   -  ka * depot
# dcentral  <-  ka * depot  -  ke * central
# test_plot <- central / v1"
# model %>% cat
#
# OD_input <-   tibble(Output = c("test"),
# Group = c("1","2"), Proto = "1", TimeSample = "c(1,5,10)", add = "0.3F", prop = "0.2F", nidgroup = 30, cov = c("cl = 3, v1 = 2", "cl = 2, v1 = 4"))
#
# OD_input <-   tibble(Output = c("test"),
#                      Group = c("1","2"), Proto = "1", TimeSample = "c(1,5,10)", add = "0.3F", prop = "0.2F", nidgroup = 30, cov = c("", ""), tlag = F)
#
# OD_input <-   tibble(Output = c("test"),
#                      TimeSample = "c(1,5,10)",
#                      add = "0.3F", prop = "0.2F",
#                      nidgroup = 30,
#                      cov = c("cl = 3, v1 = 2", "cl = 2, v1 = 4"))
#
# OD_input <-
# tibble(Output = c("test", "depot"),
#                    Group = c("1","1"), Proto = "1",
#                    TimeSample = "c(1,5,10)",
#                    add = "0.1", prop = "0.2",
#                    nidgroup = 30,
#                    cov = c("cl = 3, v1 = 2", "cl = 2, v1 = 4"))
#
# # #
# states <- tibble(Cmt = c("depot", "central"), t0 = c(0,0))
# #
# events <- tibble(tlag = F, ADM = 1, var = "depot", time = "c(0,40)", use = T, value = 50, method = "add", Proto = "1:2", F = FALSE, Perf = "rate", Perf_num = 5)
# events2 <- events
# parameters <- tibble(Param = c("ka", "cl", "v1"), Value = c(0.1,3,5)) %>%
#   mutate(E = c("Esti", "fix","Esti"), Distrib = c("logN", "logN", "logN"))
# # #
# diagOmega <- tibble( ka = c("0.1", "0", "0"), cl = c("","0.1F", "0"), v1 = c("", "", "0") )

# test <-
# rm(list = c( "feps_ODE_compiled")) # , "popedResult"
  # pecc_PopEd(model, OD_input, states, events, parameters, diagOmega, outputExpr = F)
 # test <-  pecc_PopEd(model, OD_input, states, events2, parameters, diagOmega, outputExpr = T);test
# test <- eval(test,envir = new.env(parent = global_env()))

#' @export
pecc_PopEd <- function(model, OD_input, states, events, parameters, diagOmega, outputExpr = F){

# Replace factorr by character for simplicity
OD_input <- OD_input %>%
  mutate(Proto = as.character(Proto),
         Group = as.character(Group),
         Output = as.character(Output))

# Create a list of expression
ex <- list()

# First, we need to load the package
ex$loadPackage <- expr(library(PopED))



# perfusion handling ------------------------------------------------------
evaltime <- F

if("Perf_num" %in% names(events)){


  events$Perf_num <- as.double( events$Perf_num)
  events$value <- as.double( events$value)
  evaltime <- T
  ## model modification

  toadd <- unique(events$var[events$Perf != "None"])

  model <- str_split(model, "\n")[[1]]

  for(a in toadd){

    model <- c(model, paste0("dPerf_", a, " <- 0"))
    states <- add_row(states, Cmt = paste0("Perf_", a), t0 = 0)

    linetemp <- grep(paste0("d", a, " *<-"), model)
    model[linetemp] <- paste0(model[linetemp], " + Perf_", a)

  }



  model <- paste0(model, collapse = "\n")



  ## perf on
  perfon <- events %>%
    filter(Perf != "None") %>%
    mutate(value = if_else(Perf == "rate", Perf_num, value/Perf_num)) %>%
    mutate(var = paste0("Perf_", var))


  ## perf off
  perfoff <- events %>%
    filter(Perf != "None") %>%
    mutate(time = if_else(Perf == "time", paste0(time, " +  ", Perf_num),
                          paste0(time, " + ", value/Perf_num))) %>%
    mutate(value = 0, method = "rep", var = paste0("Perf_",var))

  events <- bind_rows(perfon, perfoff, events %>% filter(Perf == "None"))


}
# end perfusion handling




# Desolve formula ---------------------------------------------------------


# Then, we need the deSolve formula with PopEd syntaxe
 modell <-
  deSolve_pecc(model)$model %>%
  deparse %>%
  gsub(pattern = " *function\\(t, state, parameters\\)", replacement = " function(Time, State, Pars)") %>%
  gsub(pattern = "state, parameters", replacement = "State, Pars") %>%
   paste0(collapse = "\n")

# Which is ofc an expression
ex$model <-  expr(modelpopEd <- !!parse_expr(modell))

#ex$dumb<- expr(print(modelpopEd))
# If there is a colonne "use" for events (in shinyapp), use it
if("use" %in% names(events)) events <- events %>%
  filter(use == T)

# evaluating and unnesting both time and protocol to have complete events df
events <- events %>%
  mutate(time = map(time, ~eval(parse_expr(.x)))) %>%
  unnest() %>%
  mutate(Proto = map(Proto, ~eval(parse_expr(.x)))) %>%
  unnest() %>%
  select(var, time, value, method, Proto) %>%
  mutate(var = as.character(var), method = as.character(method),
         time = as.double(time), value = as.double(value), Proto = as.character(Proto))

# Joining ODinput to protocol
events_ODinpu <- events %>%
  left_join(by = "Proto", OD_input %>% select(Proto, Group) %>%
              group_by(Group) %>%
              slice(1)) %>%
  mutate(Proto = as.integer(Proto)) %>%
  mutate(Group = as.integer(Group)) %>%
  arrange(Proto, time) %>%
  filter(!is.na(Group))

# And this is an expression to keep
ex$table <-
  expr(events_input <- data.frame(!!!map(events_ODinpu, ~ expr(!!.x))))

#ex$dumb2<- expr(print(events_input))
# handle covariate --------------------------------------------------------


## We need to know cov that will not be estimated (even if E == "Esti")
str_split(OD_input$cov[[1]], ",")[[1]] %>%
  gsub(pattern = "( *=.+)| ", replacement = "") -> covariates

## Only keep paramaters that we need to keep
names_par_est <- parameters %>%
  filter(E == "Esti", !Param %in% covariates) %>%
  pull(Param) %>%
  paste0("_pop")

# if(length(covariates) > 0 ){
# parameters <- parameters %>%
#   mutate(E = if_else(Param %in% covariates,"fix", E ))
# }
# output_name <- unique(OD_input$Output)[[1]]

# handling the biggest part -----------------------------------------------


# line containing initial values of compartiments
initial_values <- paste("c(", paste0(states %>%
                                       mutate(lines = paste0(Cmt, " = ", t0)) %>%
                                       pull(lines), collapse = ", "), ")")

## lines containing the output with proper switch
OD_input %>%
  distinct(Output) %>%
  rowid_to_column() %>%
  mutate(line = paste0("y[model_switch == ", rowid, "] <- out[ ,\"",Output,"\"][model_switch == ", rowid, "]")) %>%
  pull(line)  -> linestest

# If only one output, we can remove the modelèswitch part
if(length(linestest) == 1) linestest <- gsub("\\[model_switch == 1\\]", "", linestest)

# Finally, we can have the proper second model function!
ex$model2 <-
  expr(model2 <- function(model_switch, xt, parameters, poped.db){
    with(as.list(parameters),{
      A_ini <- !!parse_expr(initial_values)
      # eventdat <-  data.frame(!!!map(events, ~ expr(!!.x)))

      times_xt <- drop(xt)

      eventdat <- events_input %>%
        filter(Group == group)

      times <- sort(c(times_xt,eventdat$time))


      out <- ode(A_ini, times, modelpopEd, parameters, events = list(data = eventdat ))#atol=1e-13,rtol=1e-13)

      out <- out[match(times_xt,out[,"time"]),]
      y <- xt*0

      !!!parse_exprs(linestest)
      # y = out[ ,!!output_name]

      y=cbind(y)
      return(list(y=y,poped.db=poped.db))
    })
  }
  )

#ex$dumb3<- expr(print(model2))

#If a covariate is forced, we need to remove it from var-cov matrix
for(a in covariates[covariates %in% names(diagOmega) ]){

 ntorem <- which(names(diagOmega) == a)
 diagOmega <- diagOmega[-ntorem, - ntorem]
}

namesOmega <- names(diagOmega)

# For now, I only use diagonal so we can't have covariance
# Of course this is to be improved
diagOmega <- diagOmega %>%as.matrix() %>%  diag

# we estimte only if we can as.double a value, such as
# "0.01F" means we don't want to estimate this paramater
estimomega <- as.double(!is.na(as.double(diagOmega)))
# and then we remove the F part and transform to double
diagOmega <- as.double(gsub("F$","", diagOmega))
# andapply the propernames
names(diagOmega) <- namesOmega
names(estimomega) <- namesOmega

# If a omega is equal to 0 then it need to be removed
if(sum(diagOmega == 0) > 0){

  # first remove from estimatation
  estimomega <- estimomega[-which(diagOmega == 0)]
  # then to names
  namesOmega <- namesOmega[-which(diagOmega == 0)]
  # finally to values
  diagOmega <- diagOmega[-which(diagOmega == 0)]
}
# print("here")

# Take all the parameters
parameters %>%
  #remove covariates
  filter(! Param %in% covariates) %>%
  # attribute a number to each one of them
  rowid_to_column("nTheta") %>%
  # add omega value with a number for each of them
  left_join(tibble(Param = names(diagOmega), omega = diagOmega) %>%
              rowid_to_column("nOmega")) %>%
  ## compute the eta relationship from theta
  mutate(eta = case_when(Distrib == "logN" & !is.na(nOmega) ~ paste0(" * exp(b[",nOmega,"])"),
                         Distrib == "Norm" & !is.na(nOmega) ~ paste0(" + b[",nOmega,"]"),
                         T ~ "")) %>%
  # and write the complete line/formula of the parameter
  mutate(test = paste0(Param, " = ", "bpop[", nTheta,"]", eta))-> parameters2

# initiate a compting system
ncov <- compteur()
# compute the line for covariate
map_chr(covariates, ~ paste0(.x, " = a[",ncov()+ 1, "]")) %>%
  paste(collapse = ", ") %>%
  {paste0(", ", .)} -> covdecla

#remove the line if in fact there is no covariate
if(covariates[[1]] == "") covdecla <- ""

# parameters %>%
#   filter(! Param %in% covariates) %>%
#   rowid_to_column() %>%
#   mutate(test = paste0(Param, " = ", "bpop[", rowid,"] * exp(b[", rowid, "])")) %>%
#   pull(test) %>%
#   paste(collapse = ", ") %>%
#   {paste0("c(",.,", group = a[1]", covdecla,")")} %>%
#   parse_expr -> parshort
## carfeull !

#finally, we can create the all desired vector
parshort <- parse_expr(paste0("c(",paste0(parameters2$test, collapse = ","),", group = a[1]", covdecla,")"))
# parameters2
# print("parshort")
# print(parshort)


# print("vectorPar")
# print(vectorPar)

# and here is the fg function !
ex$fg <-
  expr(
    fg <- function(x,a,bpop,b,bocc){
      parameters= !!parshort

      return( parameters )
    }
  )

# ex$discrete <- expr(discrete_a <- cell(2,1))

# poped.database
# pk_sampling <- c(0.5,1,2,4)
# pk_sampling <- isolate(eval(parse_expr(input$samplingTime)))


# We need to compute initial value for theta
parameters %>%
  filter(! Param %in% covariates) %>%
  rowid_to_column() %>%
  mutate(test = paste0(Param , " = ", Value)) %>%
  pull(test) %>%
  paste(collapse = ", ") %>%
  {paste0("c(",.,")")} %>%
  parse_expr -> vectorPar
# print("poed")

# OD_input <- bind_rows(OD_input, OD_input)

# and handle error system...
OD_input_error  <- OD_input %>%
  group_by(Output) %>%
  slice(1)
## error model handling
# if there is only 1 output
if(length(unique(OD_input$Output)) == 1){
  OD_input_add_v <- gsub("F", "", unique(OD_input_error$add)) %>% as.double
  OD_input_add_fixed <- grep("F", unique(OD_input_error$add))
  if(length(OD_input_add_fixed) == 0){
    OD_input_add_fixed <- 1
  }else{
    OD_input_add_fixed <- 0
  }

  OD_input_prop_v <- gsub("F", "", unique(OD_input_error$prop)) %>% as.double
  OD_input_prop_fixed <- grep("F", unique(OD_input_error$prop))
  if(length(OD_input_prop_fixed) == 0){
    OD_input_prop_fixed <- 1
  }else{
    OD_input_prop_fixed <- 0
  }

  if(OD_input_add_v > 0 &  OD_input_prop_v > 0){
    fError_fun="feps.add.prop"
    sigmas_names <- c("add", "prop")[as.logical(c(OD_input_add_fixed, OD_input_prop_fixed))]
    sigmas <- expr(c(add = !!OD_input_add_v, prop = !!OD_input_prop_v))
    sigmasbool <-  expr(c(!!OD_input_add_fixed, !!OD_input_prop_fixed))
  }else if(OD_input_add_v == 0 &  OD_input_prop_v > 0){
    fError_fun="feps.prop"
    sigmas_names <- c("prop")[as.logical(OD_input_prop_fixed)]
    sigmas <- expr(c(prop = !!OD_input_prop_v))
    sigmasbool <-  expr(c(!!OD_input_prop_fixed))
  }else if(OD_input_add_v > 0 &  OD_input_prop_v == 0){
    fError_fun="feps.add"
    sigmas_names <- c("add")[as.logical(OD_input_add_fixed)]
    sigmas <- expr(c(add = !!OD_input_prop_v))
    sigmasbool <-  expr(c(!!OD_input_add_fixed))
  }
  # if there is more than 1 output now
}else{
   espcompt <- compteur()

  OD_input_error %>%
    rowid_to_column() %>%
    mutate(add = gsub("F", "", add)) %>%
    mutate(prop = gsub("F", "", prop)) %>%
    mutate(declar = pmap_chr(list(Output, add, prop,rowid), function(Output, add, prop,rowid){

      if(prop == 0) linesProp <- "" else linesProp <- paste0(" * (1 + epsi[, ",espcompt(),"])" )
      if(add == 0) linesAdd <- "" else linesAdd <-  paste0(" + epsi[, ",espcompt(),"]")

      paste0("Output", rowid, " <- y",linesProp, linesAdd)

    })) %>%
    mutate(modif = map_chr(rowid, function(rowid){

      paste0("y[MS == ", rowid, "] <- Output", rowid,"[MS == ", rowid, "]")

    })) %>%
    select(Output, declar, modif) -> liness


  ex$errorf <-   expr(feps_ODE_compiled <- function(model_switch,xt,parameters,epsi,poped.db){
    ## -- Residual Error function
    MS<-model_switch
    y <- model2(model_switch,xt,parameters,poped.db)[[1]]


    !!!parse_exprs(liness$declar)
    !!!parse_exprs(liness$modif)


    return(list(y=y,poped.db=poped.db))
  })


  #ex$dumb4<- expr(print(feps_ODE_compiled))
  OD_input_error %>%
    mutate(addprop = map2(add, prop, ~ c(.x, .y))) %>%
    pull(addprop) %>%
    reduce(c) -> sigmas


  OD_input_error %>%
    pull(Output) %>%
    map(~ paste0(.x,"_", c("add", "prop"))) %>%
    reduce(c) -> sigmas_names

  sgimasbool <- as.double(!grepl("F", sigmas))

  sigmas <- as.double(gsub("F", "", sigmas))

  sgimasbool <- sgimasbool[sigmas != 0]
  sigmas_names <- sigmas_names[sigmas != 0]
  sigmas <- sigmas[sigmas != 0]
  names(sgimasbool) <- sigmas_names
  names(sigmas) <- sigmas_names
  sigmas_names <- sigmas_names[as.logical(sgimasbool)]


  fError_fun="feps_ODE_compiled"
  sigmasbool <- expr(c(!!!sgimasbool))
  sigmas <- expr(c(!!!sigmas))

}


## The we need to handle time with
## matrix of time (as string or expression)
## number of output if several
OD_input %>%
  group_by(Output) %>%
  nest() %>%
  rowid_to_column("noutput") %>%
  unnest() %>%
  group_by(Group) %>%
  nest() %>%
  # {.[[1, 2]]} -> x
  mutate(time = map_chr(data, function(x){

    if(length(unique(x$Output)) == 1){
      return(unique(x$TimeSample))
    }else{
      return( paste0("c(", paste0(x$TimeSample, collapse = ", "),")"))
    }
  })) %>%
  mutate(output = map_chr(data, function(x){

    temp <- x %>%
      mutate(TimeSample = map(TimeSample,~ eval(parse_expr(.x)))) %>%
      unnest() %>%  pull(noutput)

    return(paste0("c(", paste0(temp, collapse = ", "),")"))

  })) -> tempTime

# xt is the list of all time
xt <- parse_expr(paste0("list(", paste(tempTime$time, collapse = ","),")"))
# model_switch is the list of output
model_switch <- parse_expr(paste0("list(", paste(tempTime$output, collapse = ","),")"))
# print("xt okay")
# by default group, but we must add covariate if we have some

# we need to compute the groups of covariates
OD_input %>%
  select(Group, cov) %>%
  group_by(Group) %>%
  slice(1) %>%
  mutate(a = map2_chr(Group, cov, function(Group, cov){

    if(gsub(" ", "", cov) == "") return(paste0("group = ", Group))

    return(paste0("c(group = ", Group, ", ", cov,")"))
  })) %>%
  pull(a) %>%
  paste0(collapse = ", ") -> a

# if(gsub(" ", "", OD_input$cov[[1]]) != "")
  a <- paste0("list(", a, ")")

a <- parse_expr(a)
# a <- parse_expr(paste0("list(", paste("group = ", unique(OD_input$Group), collapse = ","),")"))
# print("a okay")

# We need the groupsize
groupsize <- parse_expr(paste0("c(", paste(OD_input %>%
                                             group_by(Group) %>%
                                             slice(1) %>%
                                             pull(nidgroup), collapse = ","),")"))
# print("groupsize okay")

## And finally gather all this information into poped.db
## object (the most complex of PopEd)
ex$poped <- expr(poped.db <- create.poped.database(ff_fun="model2",
                                                    fError_fun=!!fError_fun,
                                                    fg_fun="fg",
                                                    groupsize=!!groupsize,
                                                    m=!!length(unique(OD_input$Group)),      #number of groups
                                                    sigma= !!sigmas,
                                                    notfixed_sigma = !!sigmasbool,
                                                    bpop= !!vectorPar,
                                                    notfixed_bpop= !!(parameters %>% filter(! Param %in% covariates) %>%  mutate(E = if_else(E == "Esti", T, F)) %>% pull(E) %>% as.double()),
                                                    d= !!diagOmega, #c(CL=0.1,V1=0.1,KA=0.1,Q= 0.1, V2= 0.1, Favail=0.1)
                                                    notfixed_d = !!estimomega,
                                                    xt= !!xt,
                                                    model_switch = !!model_switch,
                                                    a = !!a))
#ex$dumb5<- expr(print(poped.db))
ex$output <- expr(popedResult <- list())

ex$plot <- expr(popedResult$plot_OD <- plot_model_prediction(poped.db,model_num_points = 500, facet_scales = "free"))
#ex$dumb6<- expr(print(popedResult$plot_OD ))
ex$evaluate <- expr(popedResult$result_OD <- evaluate_design(poped.db))



parameters2 %>%
  filter(omega > 0 & Param %in% names(estimomega)[estimomega == 1] ) %>%
  pull(Param) %>%
  paste0("_omega") -> names_omega_est

if(length(names_omega_est) == 1 & names_omega_est[[1]] == "_omega"){

  names_omega_est <- character()
}

names <- c(names_par_est,
           names_omega_est, sigmas_names)
# ex$rename <- expr(names(result_OD$rse)[1:!!length(names)] <- !! names)
ex$names <- expr(names(popedResult$result_OD$rse) <- !! names)
print("aeza")

ex$final <- expr(popedResult)


# import line below to avoid any environment conflict after
ex <- map(ex, ~ parse_expr(deparse(.x, width.cutoff = 500) %>% paste0(collapse = "\n")))
ex <- expr({!!!ex})


if(outputExpr == T) return(ex)
# print(ex)
 # output$code <-ex
return(eval(ex, envir = globalenv()))

}



# output <- pecc_PopEd(model ,OD_input  , states , events , parameters )
#
# output
# # print("parameters")
# # print(parameters)
# #fg
#
# for(a in output){
#   print(a)
#   eval(a)
# }
#
# eval(parse_expr(paste0(deparse(output$model2), collapse = "\n")))
#
# eval(expr({!!!output}))
#
# map(output, function(x){
#
#   try(eval(parse_expr(paste0(deparse(x), collapse = "\n"))))
#
# })



# # tlag handling -----------------------------------------------------------
#
# if(sum(events$tlag) > 0){
#
#   listtlag <- list()
#
#   for(a in unique(events$ADM[events$tlag == T])){
#
#     tlagnamepar <- paste0("tlag_", a) %>% gsub(pattern = "_$", replacement = "")
#     listtlag[[a]] <-  expr(events$time[events$ADM == !!a] <- events$time[events$ADM == !!a] + parameter[[!!tlagnamepar]])
#
#   }
# }else{
#   listtlag <- NA
# }
# # Bioavailibilty handling -------------------------------------------------
#
#
# if(sum(events$F) > 0){
#
#   listBioAv <- list()
#
#   for(a in unique(events$ADM[events$F == T])){
#
#     bioavnamepar <- paste0("BioAv_", a) %>% gsub(pattern = "_$", replacement = "")
#     listBioAv[[a]] <-  expr(events$value[events$ADM == !!a] <- events$value[events$ADM == !!a] * parameter[[!!bioavnamepar]])
#
#   }
# }else{
#   listBioAv <- NA
# }


