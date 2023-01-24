# # model <- "dCentral <- -ke * Central"
# model <- "dGut <- - ka * Gut
# dCentral <- - Cl/V * Central  + ka * Gut
# conc_plot <- Central / V"
# parameters_dff <- tibble(Param = c("Cl", "ka", "V"), Value = c(15,3,"c(30,35)"), Distrib = c("logN", "logN", "logN"))
# matrix_etas <- matrix(c(0.3,0.05,0.05,0.05,0.3,0.05,0.05,0.05,0.3), ncol = 3, dimnames = list(c("Cl", "ka", "V"), c("Cl", "ka", "V")))
# matrix_etas <- matrix(c("0.3",0,0,"","0.2",0,"","","0.3F"), ncol = 3, dimnames = list(c("Cl",  "V","ka"), c("Cl",  "V","ka")))
# matrix_etas <- as.data.frame(matrix_etas)
# matrix_eta <- matrix_etas
# n <- 100
# parameters_dff <- tibble(Param = "ke", Value = 15, Distrib = "logN")
# matrix_etas <- c(ke = "0.1") -> matrix_eta
# n <- 10
## Role of this function: create set of parameters to inject into the model
#'@param n number of simulation. If n = 0, then no need to make simulations
#'@param parameter_df A tibble containing colonum Par, Value and Distrib
#'@param matrix_eta Variance covariance Matrix
#'@param sd booléen, wethere the matrix represent variance or standard error
#'@param returnExpr booléen, wethere whe sould just return the expression or evaluate it
#'@param matrixShiny is the matrix directly the matrix or just half of it
#'@export

random_etas <- function(n = 100, parameter_df = parameters_dff, matrix_eta = matrix_etas, sd = T, returnExpr = T, matrixShiny = T){



  # because now and endinf F mean "fixed", s we have to remove that !

  rownmatix <- rownames(matrix_eta)
  matrix_eta %>%
    map_dfr(~ gsub("F$","",.x)) -> matrix_eta
  # print("what")
  rownames(matrix_eta) <- rownmatix


### ex will be the output, a  list of expression
  ex <- list()
  # print("what")
# we will need to evaluate it. If the input is in double it will create problem after that
parameter_df$Value <- as.character(parameter_df$Value )





# Take the input
parameter_df %>%
  # Evaluate the colonne, each Value being an atomic vector
  mutate(Value = map(Value, ~ eval(parse_expr(.x)))) %>%
  # Compute the number of different value per paramater
  mutate(n = map_dbl(Value,~ length(.x))) -> parameter

# Get which paramater have several tried values
parameter %>% filter(n >1) %>% pull(Param) -> severalValue

# Get set
parameter %>%
{invoke(.fn = crossing, .$Value)} -> parameter
names(parameter) <- parameter_df$Param

# Compute name of each set
parameter$nameparset <- ""

if(length(severalValue) > 0){


  for(a in severalValue){

    parameter$nameparset <- paste0(parameter$nameparset, a, " = ", parameter[[a]], " \n ")
  }

  parameter$nameparset <- gsub(" \\\n $", "",  parameter$nameparset)


}

# print("bordel")

###### When n = 0, it means we do not need to makek simulation, nothing to do more
if(n == 0){

  ex$parameter <- expr( ind_param <- tibble(!!! map(parameter, ~ expr(!!.x))) )
  if(returnExpr == T)  return(ex)
  if(returnExpr == F)   return(eval(ex))

}


ex$parameter <- expr( parameters <- tibble(!!! map(parameter, ~ expr(!!.x))) )
############################# if n > 0
  ### If we do not have a complete matrix, recrete one
  if(matrixShiny == T){


  temp2 <-   matrix_eta %>%
      mutate(Par = colnames(matrix_eta)) %>%
      gather(-Par, key = Par2, value = value) %>%
      mutate(test = map2_chr(Par, Par2, ~ paste(c(.x, .y)[order(c(.x,.y))], collapse = "-")))


    temp2$value[temp2$value == ""] <- NA


    ## recreate a whole matrix
    for(a in which(is.na(temp2$value))){

      test2 <- temp2[[a, "test"]]

      value <- temp2 %>%
        filter(test == test2 & (!is.na(value) )) %>%
        pull(value)


      temp2[[a, "value"]] <- value

    }


    temp2 %>%
      select(-test) %>%
      spread(key = Par2, value = value) %>%
      select(-Par) %>%
      map_df(~ as.double(.x)) %>%
      as.matrix -> matrix_eta




  }

# switch from standard error to matrix
   if(sd == T ) matrix_eta <- matrix_eta ^ 2

   # remove value when diagonal = 0, for parameter without variability
      # which(diag(matrix_eta) == 0)
    # which(diag(matrix_eta) != 0)
diagneg <- which(diag(matrix_eta) == 0)
novar <- colnames(matrix_eta)[diagneg]

if(length(diagneg)>0){

  parameter_df$Distrib[ parameter_df$Param %in% novar] <- "No"
  matrix_eta <- matrix_eta[-diagneg,-diagneg]

}




# print("haaa")
   ex$matrix <-   expr(matrix_eta <- tibble(!!! map(matrix_eta %>% as.data.frame , ~ expr(!!.x))) %>% as.matrix)

# issue of name when we have only 1 parmater


### compute eta (nrow = number of simulation * number of parameter * number of set)
  ex$eta <-   expr(eta <- matrix(rnorm(!!n * nrow(matrix_eta)  ), ncol = nrow(matrix_eta)) %*% chol(matrix_eta))

# repete this matrix for each subset in case of multiple value for a same parameter
if(nrow(parameter) >1) ex$eta2 <- expr(eta <- apply(eta, 2, function(x) rep(x, nrow(parameters))))

  ### three compute final with thetas and etas
ex2 <- list()


    ex2$theta <- expr(ind_param <- map_dfc( parameters, ~ rep(.x, !!n))  %>% arrange(nameparset))

    ## if there is only one eta, problem with names, so we have to handle that separately
    if(length(grep("logN" , parameter_df$Distrib)) > 0 & NCOL(matrix_eta) > 1){
      ex2$logN0 <- expr( parameterslogN <- !!parameter_df$Param[parameter_df$Distrib == "logN" & !is.na(parameter_df$Distrib) ])
      ex2$logN <-  expr(ind_param[parameterslogN] <- ind_param[parameterslogN] * exp(eta[ ,parameterslogN]))
    }else if(length(grep("logN" , parameter_df$Distrib)) > 0 & NCOL(matrix_eta) == 1){

      ex2$logN0 <- expr( parameterslogN <- !!parameter_df$Param[parameter_df$Distrib == "logN" & !is.na(parameter_df$Distrib)])
      ex2$logN <-  expr(ind_param[parameterslogN] <- ind_param[parameterslogN] * exp(eta))

    }

    if(length(grep("Norm" , parameter_df$Distrib)) > 0 & NCOL(matrix_eta) > 1){
      ex2$Norm0 <- expr( parametersNorm <- !!parameter_df$Param[parameter_df$Distrib == "Norm" & !is.na(parameter_df$Distrib)])
      ex2$Norm <-  expr(ind_param[parametersNorm] <- ind_param[parametersNorm] + eta[ ,parametersNorm])
    }else if(length(grep("Norm" , parameter_df$Distrib)) > 0 & NCOL(matrix_eta) == 1){

      ex2$Norm0 <- expr( parametersNorm <- !!parameter_df$Param[parameter_df$Distrib == "Norm" & !is.na(parameter_df$Distrib)])
      ex2$Norm <-  expr(ind_param[parametersNorm] <- ind_param[parametersNorm] * exp(eta))

    }


#
#     if(nrow(parameter) > 1){
#
#
#       ex2 <- expr( ind_param <-
#
#         parameters %>%
#         rownames_to_column("setPara") %>%
#         group_by(setPara) %>%
#         nest() %>%
#         #  %>%
#           mutate(outp = map(data, function(parameters){
#
#           !!!ex2
#
#           return(ind_param)
#
#         })) %>%
#
#           unnest(outp)
#       )
#
#
#     }else{
#
#
#
#       ex2$setPara <- expr(ind_param <- ind_param %>% mutate(setPara = "1"))
#
#     }



    if(returnExpr == T){

      # return(expr({!!!c(ex, ex2, expr(ind_param))}))
      return(c(ex, ex2))

      } else{

        for(a in c(ex, ex2)) eval(a)

        return(ind_param)

    }

 }

# parameter <- random_etas(parameter_df = parameters_dff , matrix_eta = matrix_etas, n = 4, returnExpr = T)
# # # # # #
# # # # # #
# # # # # #
# parameter <- parameter %>% rename(BioAv = ka)
# events <- tibble(Proto = c("1", "1:2"), var = c("Gut", "Gut"), time = c("1:24", "0"), value = 5:6, method = c("add", "add"), use = c(T,T), delete = c(F,F), tlag = F, ADM = "", F = F,
# Perf = c("rate", "time"), filterPlot = c("Dose == 50", "Dose ==100"), Perf_num = c(5,5)) %>%
# slice(1:2)# time or rate
# model <- "dX <- X * Cl
# conc_plot <- X/V"
# states = tibble(Cmt = c("Gut", "Central"), t0 = "0")

# model <- "
# ke <- Cl/V
# Conc_plot <- ka/V * exp(- t * ke)"
# parameters_dff
# #
# states = tibble(Cmt = "X", t0 = "0.1")
# states = tibble(Cmt = "Central", t0 = "0")
# error <- tibble(output = "conc", err_add = "0.3F", err_prop = "0.1")

#'@param parameter output from random_etas function
#'@param model simplified model
#'@param states states, initial values of compartment
#'@param events deSolve events table
#'@param times  time sequence
#'@param timesF do we need to substitute time???
#'@export
make_simulations <- function(parameter, model, states, events, times = seq(0,100,1), timesF = F, error = NA, Progress = T, returnExpr = T){


  ## Perfusion handling

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
      states <- add_row(states, Cmt = paste0("Perf_", a), t0 = "0")

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


  transformeddmodel <- deSolve_pecc(model)
  # ind_param <- mtcars  %>% mutate(setPara = 1)
  # print("okay")




  if("data.frame" %in% class(parameter)){

    # add a column and nest of line per id
    # Important ! Here paramters will be a tibble with
    # id: number of the row
    # nameparset: name of the model (for plotting after)
    # paramter: dataframe of parameter
    ex <- list()

    ex$par <- expr(parameters_df <- !!substitute(parameter) %>% rownames_to_column("id") %>% group_by(id, nameparset) %>% nest(.key = "parameter"))

  }else{
    ex <- parameter
    ex$par <- expr(parameters_df <- ind_param %>% rownames_to_column("id") %>% group_by(id, nameparset) %>% nest(.key = "parameter"))

  }

  if(timesF == F) ex$time <- expr(times <- !!enexpr(times))
# eval(paramater)
  ### handling parameters

     # eval(ex$par)
  ### handling events
  # events <-
# print("par")

# if we have equation differntial
if(length(transformeddmodel$state)>0){
## events handling

  if("use" %in% names(events)){
  events <- events %>%
    filter(use == T) %>%
    select(-starts_with("use"), - starts_with("delete"))
}
  # print("events")
  # print(events)
# For label in plot
  # events <- map_if(events %>% select(-F, -tlag, -Perf), is.factor,~ as.character(.x)) %>%
  #          as.data.frame()


  ## tlag handling
  if(sum(events$tlag) > 0){

    listtlag <- list()

    for(a in unique(events$ADM[events$tlag == T])){

      tlagnamepar <- paste0("tlag_", a) %>% gsub(pattern = "_$", replacement = "")
      listtlag[[a]] <-  expr(events$time[events$ADM == !!a] <- events$time[events$ADM == !!a] + parameter[[!!tlagnamepar]])

    }
  }else{
     listtlag <- NA
  }

  ## Bioavailability handling

  if(sum(events$F) > 0){

    listBioAv <- list()

    for(a in unique(events$ADM[events$F == T])){

      bioavnamepar <- paste0("BioAv_", a) %>% gsub(pattern = "_$", replacement = "")
      listBioAv[[a]] <-  expr(events$value[events$ADM == !!a] <- events$value[events$ADM == !!a] * parameter[[!!bioavnamepar]])

    }
  }else{
    listBioAv <- NA
  }



  ## r

  events <- events %>%
    select(Proto, var, time, value, method, ADM) %>%
    mutate(var = as.character(var), method = as.character(method))


  events_expr <- expr(tibble(!!! map(events, ~ expr(!!.x))))

  # several time or proto
  if(evaltime == T | max((events %>%
     mutate(nmax = map_dbl(time, ~ length(eval(parse_expr(.x))))))$nmax) > 1){


    events_expr <- expr(!!events_expr %>%
           mutate(time = map(time, ~ eval(parse_expr(.x)))) %>%
           unnest(time))

  }





  if( max((events %>%
           mutate(nmax = map_dbl(Proto, ~ length(eval(parse_expr(.x))))))$nmax) > 1){


    events_expr <- expr(!!events_expr %>%
                          mutate(Proto = map(Proto, ~ eval(parse_expr(.x)))) %>%
                          unnest(Proto))

  }


  # For plot labeling
  if(length(unique(events$Proto)) == 1){
    events_expr  <- expr(!!events_expr  %>% mutate(Proto = ""))
  }else{
    events_expr  <-expr(!!events_expr  %>% mutate(Proto = paste0( "Prot ", Proto )))
  }



  ex$events <-  expr(events_df <- !!events_expr   %>%
         mutate(value = as.double(value), time = as.double(time)) %>%
         group_by(Proto) %>%
         nest(.key = "events"))




  ## states handling
  # do we need to evaluate the initial value (ex = "kin/kout")
  test_chara_states <- sum(is.na(as.double(states$t0)))

  # if there is characeter and several initial values
  # need to know how many parameters system we have we have...
  # eval(parse_expr(deparse(ex$parameter) %>% paste0(collapse = " ") %>% gsub(pattern = ".+<-", replacement = "") %>%
  # paste0(" %>% nrow"))) -> nparamatersset

  # if we need to compute and several paramters set
  if(test_chara_states > 0 ){ #& nparamatersset > 1

    ex$states <- expr(states <- tibble(!!!map(states, ~ expr(!!.x))))

    states2 <-   expr(states <- with(data = parameter, states %>%
                                       mutate(t0 = map_dbl(t0, ~ eval(parse_expr(.x))))) %>%
                        transmute(temp = paste0(Cmt, " = ", t0)) %>%
                        pull(temp) %>%
                        {paste0("c(",   reduce(., paste, sep = ", "), ")")} %>%
                        parse_expr %>%
                        eval)

    ### if there are all value -> put actual values
  }else if(test_chara_states == 0){ # & nparamatersset > 1

    # ex$states <- expr(states <- tibble(!!!map(states %>% mutate(t0 = as.double(t0)), ~ expr(!!.x))))
    ex$states <- expr( states <-   !!parse_expr(states %>%
                                                  transmute(temp = paste0(Cmt, " = ", t0)) %>%
                                                  pull(temp) %>%
                                                  {paste0("c(",   reduce(., paste, sep = ", "), ")")}))

    states2 <- NA
  }
} ### end if differntial equation


  ## handling model

  ex$model <- expr(model <- !!transformeddmodel$model)
# print("here")

# print("here")
  ##

  # & test_chara_states$test > 1
  #### If several sets of parameters
# print("okay on commence les siulations")

if(Progress == T & length(transformeddmodel$state) >0){
  ex$simulations <-  expr(
    withProgress(message = 'Making simulations', value = 0, {
      ysim <- 1
      simmax <- nrow(parameters_df) * nrow(events_df)
      simulations <- crossing(parameters_df, events_df ) %>%
        mutate(simul = map2( parameter, events, function(parameter, events){


          incProgress(1/simmax, detail = paste( ysim, "/", simmax))
          ysim <<- ysim + 1
          !!states2
          !!!listtlag
          !!!listBioAv



          as.data.frame(ode( events = list(data = as.data.frame(events)), y = states, times = times, func = model, parms = as.data.frame(parameter), method = "lsode"))



        })) %>%
        select(-parameter, -events) %>%
        unnest }) )

}else if(Progress == F & length(transformeddmodel$state) >0){

  # the exact same (copy/past) without the bar progression for shiny
  ex$simulations <-  expr(

      simulations <- crossing(parameters_df, events_df ) %>%
        mutate(simul = map2( parameter, events, function(parameter, events){


          !!states2
          !!!listtlag
          !!!listBioAv


          as.data.frame(ode( events = list(data = as.data.frame(events)), y = states, times = times, func = model, parms = as.data.frame(parameter)))


        })) %>%
        select(-parameter, -events) %>%
        unnest )


}else{

ex$simulations  <- expr(simulations <-  parameters_df %>%
       mutate(simul = map(parameter, ~  model(t = times, parameter = .x))) %>%
  select(-parameter)%>%
        unnest %>%
    mutate(Proto = "")
)
}

  ## adding error

  if(!is.na(error[[1]])){

    for(a in unique(error$output[error$output != ""])){
      errad <- gsub("F", "", error$err_add[error$output == a]) %>% as.double()
      errprop <- gsub("F", "",error$err_prop[error$output == a]) %>% as.double()


      if(errad > 0 & errprop == 0){

        ex[[paste0("error", a)]] <- expr(simulations[[!!a]] <-
                                           simulations[[!!a]] + rnorm(length(simulations[[!!a]]), 0 , !!errad) )
      }else if(errad == 0 & errprop > 0){

        ex[[paste0("error", a)]] <- expr(simulations[[!!a]] <-
                                           simulations[[!!a]] * (1 + rnorm(length(simulations[[!!a]]), 0 , !!errprop) ))


      }else if(errad  > 0 & errprop > 0){

        ex[[paste0("error", a)]] <- expr(simulations[[!!a]] <-
                                           simulations[[!!a]] * (1 + rnorm(length(simulations[[!!a]]), 0 , !!errprop) )+ rnorm(length(simulations[[!!a]]), 0 , !!errad))



      }


    }


  }

  ## Return

  if(returnExpr == T){

    # return(expr({!!!c(ex, ex2, expr(ind_param))}))
    return(ex)

  } else{

    for(a in ex) eval(a)

    return(simulations)

  }
  # expr({!!!c(ex)})


}

# make_simulations_rxode(parameter = paramInd,model =  model_demo, states = states, events = events, times = c(seq(0,100,1)), Progress = F, returnExpr =F)


#'@export
make_simulations_rxode <- function(parameter, model, states, events, times = seq(0,100,1), error = NA, Progress = T, returnExpr = T){


  # Add states different from 0 to model

  whichnot0 <- which(as.character(states$t0) != "0")
  for(a in whichnot0){

    model <-  paste0(states$Cmt[[a]],"(0) <- ", states$t0[[a]]  ,"\n", model)


  }

  # Analyse the model
  transformeddmodel <- deSolve_pecc(model)
  # ind_param <- mtcars  %>% mutate(setPara = 1)
  # print("okay")


  # number of proto

  if("use" %in% names(events)){
    events <- events %>%
      filter(use == T) %>%
      select(-starts_with("use"), - starts_with("delete"))
  }

  protonames <- map(events$Proto, ~eval(parse_expr(.x))) %>% reduce(c) %>% unique()
  protonames <- paste0("Prot ", protonames)
  if(length(unique(protonames)) == 1) protonames <- ""

  if("data.frame" %in% class(parameter)){

    # add a column and nest of line per id
    # Important ! Here paramters will be a tibble with
    # id: number of the row
    # nameparset: name of the model (for plotting after)
    # paramter: dataframe of parameter
    ex <- list()
    ex$par <- expr(parameters_df <- !!substitute(parameter) %>% crossing(Proto = !!protonames) %>% rownames_to_column("id")%>%
                     mutate(id = as.double(id)))

  }else{
    ex <- parameter
    ex$par <- expr(parameters_df <- ind_param %>% crossing(Proto = !!protonames) %>% rownames_to_column("id") %>%
                     mutate(id = as.double(id)))

  }

print("rezr")
  # if(length(unique(protonames)) == 1){
  #
  #   ex$par <- expr(!!ex$par %>%
  #                    mutate(id = as.double(id)))
  #
  # }else{
  #
  #   ex$par <- expr(!!ex$par %>%
  #                    mutate(id = as.double(id)))
  #
  # }



  ex$time <- expr(times <- !!enexpr(times))
  # eval(paramater)
  ### handling parameters

  # eval(ex$par)
  ### handling events
  # events <-
  # print("par")


  ## events handling
  events0 <- events
  events <-
    events %>%
    select(Proto, var, time, value, method, ADM, Perf, starts_with("Perf_num")) %>%
    mutate(var = as.character(var), method = as.character(method)) %>%
    # mutate(var = paste0("(", var,")")) %>%
    rename(amt = value, cmt = var) %>%
    mutate(evid = 1)

  if(sum("rate" %in% events$Perf) >0){

    events <- events %>%
      mutate(rate = if_else(Perf == "rate", as.double(Perf_num), NA_real_))

  }

  if(sum("time" %in% events$Perf) >0){

    events <- events %>%
      mutate(dur = if_else(Perf == "time", as.double(Perf_num), NA_real_))

  }




  events_expr <- expr(tibble(!!! map(events, ~ expr(!!.x))))

  if(length(unique(protonames)) == 1){
    events_expr <- expr(!!events_expr %>% mutate(Proto = ""))
  }else{

    events_expr <- expr(!!events_expr %>% mutate(Proto = paste0("Prot ", Proto)))
  }

  toobserv <- c(transformeddmodel$state, transformeddmodel$output_manual, transformeddmodel$toplot)
  toobserv <- toobserv[toobserv != ""]



  events_expr <- expr(!!events_expr %>%
                        mutate(time = map(time, ~ eval(parse_expr(.x)))) %>%
                        unnest(time))



  # to unnest
  if( max((events %>%
           mutate(nmax = map_dbl(Proto, ~ length(eval(parse_expr(.x))))))$nmax) > 1){


    events_expr <- expr(!!events_expr %>%
                          mutate(Proto = map(Proto, ~ eval(parse_expr(.x)))) %>%
                          unnest(Proto))

  }



  events_expr <- expr(!!events_expr %>%
                        bind_rows(


                          crossing(time = !!enexpr(times), evid = 0, cmt = !!toobserv[[1]], Proto = !!as.character(protonames))

                        ))




  events_expr <- expr(!!events_expr %>%
                        group_by(Proto) %>%
                        nest() %>%
                        full_join(parameters_df %>% select(id, Proto) %>% mutate(Proto = as.character(Proto))) %>%
                        unnest() %>%
                        mutate(amt = as.double(amt),   time = as.double(time)) %>% arrange(id, time) %>%
                      mutate(amt = as.double(amt), time = as.double(time)) )



  ## tlag handling
  if(sum(events0$tlag) > 0){


    tlags <- which(events0$tlag == T)

    for(a in tlags){

      adm <-  events0$ADM[[a]]
      adm <- if_else(adm == "", "1", as.character(adm))

      tlagnamepar <- paste0("tlag_",adm)

      events_expr <-  expr( !!events_expr %>%
              left_join(parameters_df %>% select(id, !!parse_expr(tlagnamepar )) %>% mutate(ADM = !!adm)) %>%
              mutate(time = time + if_else(is.na(!!parse_expr(tlagnamepar ) ),0,!!parse_expr(tlagnamepar )  )))



    }
  }

  ## Bioavailability handling

  if(sum(events0$F) > 0){



    bioav <- which(events0$F == T)

    for(a in bioav){

      adm <-  events0$ADM[[a]]
      adm <- if_else(adm == "", "1", as.character(adm))

      bioavnamepar <- paste0("BioAv_",adm)

      events_expr <-  expr( !!events_expr %>%
                              left_join(parameters_df %>% select(id, !!parse_expr(bioavnamepar )) %>% mutate(ADM = !!adm)) %>%
                              mutate(amt = amt * if_else(is.na(!!parse_expr(bioavnamepar ) ),0,!!parse_expr(bioavnamepar )  )))



    }
  }





  ex$events <-  expr(events_df <- !!events_expr%>%
                       arrange(id, time) )





    ## r












    # several time or proto



    # # For plot labeling
    # if(length(unique(events$Proto)) == 1){
    #   protolab = ""
    #   events_expr  <- expr(!!events_expr  %>% mutate(Proto = ""))
    # }else{
    #
    #   protolab <- invoke(map(parse_exprs(events$Proto), ~ eval(.x)),.fn = c) %>% unique
    #   protolab = paste0("Prot ", protolab)
    #   events_expr  <-expr(!!events_expr  %>% mutate(Proto = paste0( "Prot ", Proto )))
    # }



  ## Transformation into RxODE
# textmodel <- textmodel[-(1:2)]
# textmodel[-c(length(textmodel),length(textmodel) - ]

  ex$model <- expr(model <- !!transformeddmodel$modelrxode)



    ex$simulations <-  expr(simulations <- as_tibble(model$solve(parameters_df %>% select(-nameparset, -Proto), events_df)) %>%
                              mutate(securejoin = 1) %>%
                              left_join(parameters_df %>% mutate(securejoin = 1)))


    # %>%
      # left_join(parameters_df)
  ## adding error

  if(!is.na(error[[1]])){

    for(a in unique(error$output[error$output != ""])){
      errad <- gsub("F", "", error$err_add[error$output == a]) %>% as.double()
      errprop <- gsub("F", "",error$err_prop[error$output == a]) %>% as.double()


      if(errad > 0 & errprop == 0){

        ex[[paste0("error", a)]] <- expr(simulations[[!!a]] <-
                                           simulations[[!!a]] + rnorm(length(simulations[[!!a]]), 0 , !!errad) )
      }else if(errad == 0 & errprop > 0){

        ex[[paste0("error", a)]] <- expr(simulations[[!!a]] <-
                                           simulations[[!!a]] * (1 + rnorm(length(simulations[[!!a]]), 0 , !!errprop) ))


      }else if(errad  > 0 & errprop > 0){

        ex[[paste0("error", a)]] <- expr(simulations[[!!a]] <-
                                           simulations[[!!a]] * (1 + rnorm(length(simulations[[!!a]]), 0 , !!errprop) )+ rnorm(length(simulations[[!!a]]), 0 , !!errad))



      }


    }


  }

  ## Return

  if(returnExpr == T){

    # return(expr({!!!c(ex, ex2, expr(ind_param))}))
    return(ex)

  } else{

    for(a in ex) eval(a)

    return(simulations)

  }
  # expr({!!!c(ex)})


}


# model
# simulations <-
  # make_simulations(timesF = F, returnExpr = F,parameter = parameter,times = 1:100, model = model, states = states, events = events, Progress = F, error = error )
# returnExpr = T
# parameter = parameter
# times = 1:100
# model = model
# states = states
# events = events
# Progress = T

# simulations %>%
#   group_by(nameparset, Proto) %>%
#    slice(1:2)

# wrap <- "output"
# xlog <- F
# ylog <- T
# plot_table <- tibble(Plot = c(1,1 ), Todisplay = c("X", "conc_plot"), Check = c(T, T), Point = c(F, F), Filter_of_dataset = c("", "TIME < 20"), YTYPE = c(100,5))
# plot_table_cov  <- tibble(Plot = c(1,2), ylog = c(T, F), xlog = c(F,F), wrap = c("Output", "Output"))
# #
# test <- plot_simulations(simulations = simulations, plot_table =plot_table, plot_table_cov = plot_table_cov, n = 0)
# #
# plot_table <- tibble(Plot = c(1,1, 1), Todisplay = c("X", "output_conc", "Y"), Check = c(T, T, T), Point = c(F, F, F), Filter_of_dataset = c("", "TIME < 20", ""), YTYPE = c(100,5, NA))
# plot_simulations(simulations = simulations, plot_table =plot_table, plot_table_cov = plot_table_cov, n = 0, name_df = "mtcars", xpoint = "mpg", ypoint = "qsec")

# plot_table <- tibble(Plot = c(1,1 ), Todisplay = c("X", "conc_plot"), Check = c(T, T), Point = c(F, F), Filter_of_dataset = c("", "TIME < 20"), YTYPE = c(100,5))
# plot_table_cov  <- tibble(Plot = c(1,2), ylog = c(T, F), xlog = c(F,F), wrap = c("Output", "Output"))
# simtest <- simulations
# plot_simulations(returnExpr = F,simulations = simtest, plot_table = plot_table, plot_table_cov = plot_table_cov)
# # xpoint = "TIME"
# ypoint="DV"
#'@param simulations output from make_simulation function
#'@param plot_table shiny output, dataframe with Plot, Todisplay, Check, Point, Filter_of_dataset_ YTYPE
#'@param plot_table_cov shiny output, with for each plot the log, wrap...
#'@param xpoint name of x if point should be added
#'@param xpoint name of y if point should be added
#'@param n number of simulation (same as in random etas, just to know if we need to compute stat)
#'@param ymindisplayed to avoid stupid value of min display
#'@param ytype_header name of the column YTYPE
#'@param explo name of the dataset to use
#'@export
plot_simulations <- function(simulations, plot_table, plot_table_cov, xpoint = "", ypoint = "", n = 10, ymindisplayed = 0.0001, ymaxdisplayed = NA, ytype_header = character(), events = NULL, name_df = "explo", rmt0 = F, scalewrap = "free", returnExpr = T){

  # plottable <- mtcars
  # desolve <- expr(desolve)
  #
  # if(ymindisplayed != 0)
  #
  #   expr(!!desolve %>%
  #          filter())
###

####
  plot_table_cov$wrap <- as.character(plot_table_cov$wrap)


  if("data.frame" %in% class(simulations)){

    ex <- list()
    simulationname <- substitute(simulations)

  }else{
    ex <- simulations
    simulationname <- expr(simulations)
  }




  if(ymindisplayed > 0){

 ex$modification <- expr(for(a in names(!!simulationname)[! names(!!simulationname) %in% c("id", "nameparset", "time", "Proto")]){

   (!!simulationname)[[a]] <- if_else( (!!simulationname)[[a]]< !!ymindisplayed , 0,  (!!simulationname)[[a]])
   }
)
  }


  if(rmt0 == T) simulationname <- expr(!!simulationname %>% filter(time !=0))

  # print("ymaxdisplayed")
  # print(ymaxdisplayed)
  # print(is.na(ymaxdisplayed ))
  if(!is.na(ymaxdisplayed )){

    ex$modification2 <- expr(for(a in names(!!simulationname)[! names(!!simulationname) %in% c("id", "nameparset", "time", "Proto")]){

      (!!simulationname)[[a]] <- if_else( (!!simulationname)[[a]]>= !!ymaxdisplayed , NA_real_,  (!!simulationname)[[a]])
    }
    )

    print(ex$modification2)
  }

# if need to express quantiles
  if(n > 0){
    ex$quantiles <-  expr(quantiles <-  !!simulationname %>%
                             mutate(nameparset = fct_inorder(nameparset)) %>%
                             gather(-id, -nameparset, -Proto, - time, key = "key", value = "value") %>%
                            filter(key %in% c(!!!unique(plot_table$Todisplay))) %>%
                             group_by(time, key, nameparset, Proto) %>%
                             nest() %>%
                             crossing(tibble(q1 = seq(0.1,0.8,0.1),
                                             q2 = q1 + 0.1)) %>%
                             mutate(qx = map2_dbl(q1, data, function(q1,data){

                               unname(quantile(data$value, q1, na.rm = T))

                             })) %>%
                             mutate(qx2 = map2_dbl(q2, data, function(q2,data){

                               unname(quantile(data$value, q2, na.rm = T))

                             })) %>%
                             select(-data))
  }
  # print("izsdqsdqs")

  ## Plot handling
  ### IF YTYPE is mentionned, then recreate the whole filter

  if(length(ytype_header) >0){

    ## In order to recreate the whole filter !
    plot_table <- plot_table %>%
      # udpate Filter_of_dataset
      mutate(Filter_of_dataset = map2_chr(Filter_of_dataset, YTYPE, function(Filter_of_dataset,YTYPE){

        ## cleaning
        if(is.na(YTYPE) | is.null(YTYPE) | gsub(" *","",YTYPE) == "") YTYPE <- ""
        if(is.na(Filter_of_dataset)|is.null(Filter_of_dataset)|gsub(" *","",Filter_of_dataset) == "") Filter_of_dataset <- ""


        ### all possibiities
        if(Filter_of_dataset != "" & YTYPE != "") output <- paste0("(",Filter_of_dataset, ") & ", ytype_header," == ", YTYPE)
        if(Filter_of_dataset != "" & YTYPE == "") output <- Filter_of_dataset
        if(Filter_of_dataset == "" & YTYPE != "") output <- paste0(ytype_header," == ", YTYPE)
        if(Filter_of_dataset == "" & YTYPE == "") output <- ""

        output

      } ))
  }

  # remove plot we do not want
  if("Check" %in% names(plot_table))  plot_table <- plot_table %>% filter(Check == T)



  # print("izezeze")
  plot_all_info <- expr(tibble(!!!map(plot_table %>% left_join(plot_table_cov),~ expr(!!.x))))
  # print("test addition plotpoint")
  plot_all_info_eval <- eval(plot_all_info)
  ### if we need to add poinds
   if(sum(plot_table$Point) > 0){
  # print("test addition plotpoint")



     plottablepoint <- plot_table %>%
       filter(Point == T) %>%
       select(Plot, Todisplay, Filter_of_dataset)

    # handling protocol variation
     exprfilterPlotAdm <- expr(mutate())

    try({
# print("insidethetry")

      if(events %>%
         filter(use == T & filterPlot != "") %>%
         nrow > 1 ){
     events %>%
       filter(use == T & filterPlot != "") %>%
       mutate(Proto = paste0( filterPlot, " ~ \"Prot ", Proto,"\"" )) %>%
       pull(Proto) %>%
       parse_exprs %>%
       {expr(case_when(!!!.))} -> exprfilterPlotAdm
# print(exprfilterPlotAdm)
        exprfilterPlotAdm <- expr(mutate(Proto = !!exprfilterPlotAdm))
      }else if( events %>%
                filter(use == T & filterPlot != "") %>%
                nrow == 1){

        events %>%
          filter(use == T & filterPlot != "") %>%
          mutate(Proto = paste0( filterPlot, " ~ \"\"" )) %>%
          pull(Proto) %>%
          parse_exprs %>%
          {expr(case_when(!!!.))} -> exprfilterPlotAdm


        exprfilterPlotAdm <- expr(mutate(Proto = !!exprfilterPlotAdm))
      }
    })

    ##


     df_filtering <- expr(df_filtered <- tibble(!!!map(plottablepoint, ~ .x)) %>%
            # filter(Point == T) %>%
            mutate(dataset = map2(Filter_of_dataset,Todisplay, function(x,y){

              if(is.na(x) | x == ""){

                explo_temp <- crossing(!!parse_expr(name_df), key = y) %>%
                  !!exprfilterPlotAdm

                try( explo_temp <-  explo_temp %>% filter(!is.na(Proto)), silent = T)

              }else{

                explo_temp <- crossing( !!parse_expr(name_df) %>% filter_(x), key = y)%>%
                  !!exprfilterPlotAdm

                try( explo_temp <-  explo_temp %>% filter(!is.na(Proto)), silent = T)

              }

              explo_temp

            })))


    # eval(plotpoint_output)
     ex$df_filtering <- df_filtering

    plot_all_info <- expr(plot_all_info <- !!(plot_all_info) %>% left_join(df_filtered))

   }else{

     df_filtering <- NULL
     df_filtered <- NULL
   }

  # ex$plot_all_info <- plot_all_info


############################ PLot creation
   # print("plottable")
   # print(plottable)
  plot_all_info_eval  %>%
    mutate(wrap = as.character(wrap)) %>%
     # for each plot
     group_split(Plot) %>%
     # {.[[1]] }-> x
     map(function(x){
       # print("xwrap")
       # print(x$wrap)
      # print(expr(!!parse_expr(todispla)))
       col <- case_when(x$wrap == "None" ~ "mutate(color = paste0(key, \"\n\", nameparset, \"\n\", Proto, \"\n\") %>%
                        gsub(pattern = \"\n\n\n?\" , replacement = \"\n\"))",
                        x$wrap == "Output" ~  "mutate(color = paste0( nameparset, \"\n\", Proto, \"\n\")%>%
                        gsub(pattern = \"\n\n\n?\" , replacement = \"\n\"))",
                        x$wrap == "Param" ~  "mutate(color = paste0( key, \"\n\", Proto, \"\n\")%>%
                        gsub(pattern = \"\n\n\n?\" , replacement = \"\n\"))",
                        x$wrap == "Event" ~  "mutate(color = paste0( nameparset, \"\n\", key, \"\n\")%>%
                        gsub(pattern = \"\n\n\n?\" , replacement = \"\n\"))",
                        x$wrap %in% c("OP", "O|P") ~  "mutate(color =  Proto)",
                        x$wrap %in% c("OE", "O|E") ~  "mutate(color =  nameparset)",
                        x$wrap %in% c("PE", "P|E") ~  "mutate(color =  key)",
                        x$wrap == "OPE" ~  "mutate(color = \"\")")




        wrap_grid <- case_when(x$wrap == "None" ~ NA_character_,
                                                x$wrap == "Output" ~  paste0("facet_wrap(~key, scales =", scalewrap,")"),
                                                x$wrap == "Param" ~  paste0("facet_wrap(~fct_reorder(nameparset, as.double(gsub(\".+= *\", \"\", nameparset)) ), scales =", scalewrap,")"),
                                                x$wrap == "Event" ~  paste0("facet_wrap(~Proto, scales =", scalewrap,")"),
                                                x$wrap == "OP" ~  paste0("facet_wrap(~paste0(nameparset, \"\n\", key), scales =", scalewrap,")"),
                                                x$wrap == "OE" ~ paste0("facet_wrap(~paste0(Proto, \"\n\", key),  scales =", scalewrap,")"),
                                                x$wrap == "PE" ~  paste0("facet_wrap(~paste0(Proto, \"\n\", fct_reorder(nameparset, as.double(gsub(\".+= *\", \"\", nameparset)))), scales =", scalewrap,")"),
                                                x$wrap == "O|P" ~  paste0("facet_grid(key~fct_reorder(nameparset, as.double(gsub(\".+= *\", \"\", nameparset))),  scales =", scalewrap,")"),
                                                x$wrap == "O|E" ~  paste0("facet_grid(key~Proto,  scales =", scalewrap,")"),
                                                x$wrap == "P|E" ~  paste0("facet_grid(fct_reorder(nameparset, as.double(gsub(\".+= *\", \"\", nameparset)))~Proto,  scales =", scalewrap,")"),
                                                x$wrap == "OPE" ~  paste0("facet_wrap(~paste0(Proto, \"\n\", fct_reorder(nameparset, as.double(gsub(\".+= *\", \"\", nameparset))), \"\n\", key), scales =", scalewrap,")"))

        # wrap_grid <- case_when(x$wrap == "None" ~ NA_character_)

        # print("izezeze")
        # print(wrap_grid)
if(n > 0){

       plot_temp <- expr(quantiles %>%
            (!!parse_expr(col[1])) %>%
              filter(key %in% c(!!!x$Todisplay)) %>%
              mutate(forfct = as.double(nameparset)) %>%
              ggplot+
              geom_ribbon( aes(x = time, ymin = qx, ymax = qx2,  alpha = paste0(q1,"-", q2), fill = fct_reorder(color,forfct)))+
              geom_line(data = quantiles %>%  (!!parse_expr(col[1])) %>% filter(key %in% c(!!!x$Todisplay), q1 == 0.5) , aes(time,  qx, col =  fct_reorder(color, as.double(nameparset))))+
              labs(x = "Time", y = "", alpha = "quantiles", fill = "", col = "")+
              scale_alpha_manual(values = c(0.2, 0.4, 0.6, 0.8, 0.8, 0.6, 0.4, 0.2))+
              theme_bw())

}else{



  quantiles_expr <-  NULL

  plot_temp <- expr( !!simulationname %>%
                      gather(-id, -nameparset, -Proto, - time, key = "key", value = "value") %>%
                      (!!parse_expr(col[1])) %>%
                      filter(key %in% c(!!!x$Todisplay)) %>%
                      ggplot)

  if( sum(x$Point) >=1  & n == 0){ #x$Plot %in% plotpoint$Plot &


    plot_temp <- expr(!!plot_temp +
                        geom_point(data = reduce(df_filtered$dataset, bind_rows),
                                   aes(!!parse_expr(xpoint), !!parse_expr(ypoint))))

 }

  plot_temp <- expr(!!plot_temp +
                      geom_line(aes(time,  value, col = fct_reorder(color, as.double(id))), size = 1.5)+
                      labs(x = "Time", y = "", col = "")+
                      theme_bw())



}


        if( sum(x$Point) >=1  & n > 0){ #x$Plot %in% plotpoint$Plot &


          plot_temp <- expr(!!plot_temp +
                              geom_point(data = reduce(df_filtered$dataset, bind_rows),
                                         aes(!!parse_expr(xpoint), !!parse_expr(ypoint))))

        }


        # print("wrap_grid")
        # print(wrap_grid[1])
       if(!is.na(wrap_grid[1])) plot_temp <- expr(!!plot_temp + !!parse_expr(wrap_grid[1]) )
        if(x$ylog[[1]] == T) plot_temp <- expr(!!plot_temp + scale_y_log10() ) #limits = c(!!ymindisplayed , NA)
       if(x$xlog[[1]] == T) plot_temp <- expr(!!plot_temp + scale_x_log10() )
        # print("zzzz")
        # print( plotpoint$Plot)
        # print(x$Plot)
        # print(x$Plot %in% plotpoint$Plot )
          # print(x$Point)
          # print(x)




return(plot_temp)

     }) -> plots



   if(length(plots) == 1){

     plotfinal <- plots[[1]]

   }else{

     plotfinal <- expr(plot_grid(!!!plots, ncol = 1))

   }


ex$breaks2 <- expr( breaks2 <- map(-40:40,~ 1:9*10^.x) %>% reduce(c))
ex$labels2 <-  expr(labels2 <- as.character(breaks2))

ex$labels22 <-   expr(labels2[-seq(1,length(labels2), 9)] <- "")

ex$plot <- plotfinal
# return(expr({!!!c(simulation, quantiles_expr ,plot_temp)}))
 # simulations$simulations <- quantiles_expr
if(returnExpr == T){

  # return(expr({!!!c(ex, ex2, expr(ind_param))}))
  return(ex)

} else{

  for(a in ex) eval(a)

  return(eval(plotfinal))

}
}

# shinyApp(pecc_ui_reunif, pecc_server_reunif)
#
#
# shinyApp(pecc_ui, pecc_server)
#
# test <- plot_simulations(simulations = simulations, output = c("X", "output_conc"), wrap = "None",plot_table = plot_table, plot_table_cov = plot_table_cov, xpoint = "Test", ypoint = "lol")
#
#
# test %>%
#   map(function(x){
#
#     deparse(x, width.cutoff = 500) %>%.
#       paste(collapse = "\n") %>%
#       gsub(pattern = "%>%", replacement = "%>%\n" ) %>%
#       gsub(pattern = "\\+", replacement = "+\n" ) %>%
#       gsub(pattern = "  *", replacement = " " )
#
#   }) %>%
#   paste(collapse = "\n\n") %>%
#   cat
#
#
# map(c("X", "Y"), function(x) expr(!!x))
