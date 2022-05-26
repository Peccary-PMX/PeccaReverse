# ytype = tibble(output = 4L, YTYPE = NA_integer_, err_add = "0.1",
#                err_prop = "0.3", export = TRUE, rm = FALSE) %>% filter(export == T)
# add_param = ""
# y <- tibble(Cmt = c("A1", "A2", "A3"),
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
# parms  <- tibble(Param = c("CL", "KA", "Q",
#  "V1", "V2"), Value = c("0.1", "0.1",
# "5", "3", "5"), Distrib = c("logN", "logN", "logN", "logN",
# "Norm"), E = c("Esti", "Fix", "Input", "Esti", "Esti"))
#
#
#
# depot <- tibble(Proto = "1", var = "A1", time = "0",
#                       value = 50, method = "add", use = TRUE, delete = FALSE, ADM = "3",
#                       F = FALSE, tlag = FALSE, Perf = "None")
# mb_output <- tibble(output = "test", YTYPE = NA_integer_, err_add = "0.1",
#                     err_prop = "0.3", export = TRUE, rm = FALSE)

# ode_monolix(func = func, parms = parms, mb_output = mb_output, y = y, depot = depot, outputcat = F) %>% cat

#' ODE monolix
#'
#' @param func same as the ode
#' @param parms same as the ode
#' @param y same as the ode
#' @param datafile Either a .tab or an other cfl
#' @param ...
#'
#' @return
#' @export
ode_monolix <- function(func = Lorenz,mb_output, parms  = parameters, y = state, depot = "C",  datafile = "",add_param ="",ytype = data.frame(), outputcat = T,...){

  if(is.character(func)) func <- deSolve_pecc(model_demo)$model

  blocs <- list()

  init <- y

### adding biodisponibility if needed
  if(nrow(depot) > 0 ){

    depot %>%
      filter(F == T) -> test

    if(nrow(test) >0){
      newparms <- paste0("Fadm",test %>% pull(ADM))
      # previousnames <- names(parms)
      # parms <- c(parms, newparms)
      # names(parms) <- c(previousnames, newparms)

    }

    }
################


  if(add_param !=""){


    toadd <- init[gsub("0$","",add_param)]
    names(toadd) <- paste0(names(toadd),"0")

    parms <- c(parms, toadd)

    init[gsub("0$","",add_param)]<- names(toadd)

  }





  blocs$longitudinal <- "[LONGITUDINAL]\n"

  #bloc input

  blocs$input <- paste0("input ={", paste0(parms$Param, collapse = ", "),"}\n")

  regressor <- parms %>% filter(E == "Input") %>% pull(Param)

  if(length(regressor)>0){

    blocs$regressor <- map_chr(regressor, ~ paste0(.x, " = {use=regressor}\n")) %>%
      reduce(paste0)

  }
# print(depot)
  # bloc pk depot
  if(nrow(depot) > 0){
  depot %>%
      group_by(var) %>%
      slice(1) %>%
      mutate(test = if_else(F == FALSE,paste0("depot(target = ",var,  ", adm = ", ADM, ")"),
                            paste0("depot(target = ",var, ", adm = ", ADM,  ", p = Fadm", ADM,")"  )))%>%
    pull(test) %>%
    paste(collapse = "\n") %>%
      gsub(pattern = " *adm = *, *", replacement = " ")%>%
      gsub(pattern = ", *adm = *\\)", replacement = ")")-> depotemp
  }else{

    depotemp <- ""
}


  blocs$PK <- paste0("PK:
                     ", depotemp)

  #bloc EQUATION iniital condition

  blocs$equation1 <- paste0("\nEQUATION:
         ; Initial conditions
         t0 = 0
         ",
         map2(init$Cmt, init$t0, ~ paste0(.x, "_0 = ", .y)) %>%
           paste(collapse = "\n"))




  #bloc EQUATION dynamimcal model
  model_lines <- deparse(func,  width.cutoff = 500)[-c(1:3)] %>%
    gsub(pattern = "/",replacement =  " / ")

  model_lines <- model_lines[- c(grep("list\\(c\\(", model_lines): length(model_lines))]

  model_lines_2 <- model_lines %>%
    gsub(pattern = "<-", replacement = "=") %>%
    gsub(pattern = "\\^", replacement = " ^ ")




  for(x in init$Cmt){


    model_lines_2 <- gsub(paste0(" *d", x), paste0("ddt_", x), model_lines_2)

  }



  model_lines_3 <- model_lines_2 %>%
    paste0(" ") # in need this space for the next paste0(" ",x, " ")


  ## gestion of if else structure

  grep("if *\\(", model_lines_3) %>%
    walk(function(x) {

      # find the ending brackets of an if declaration
      number <- grep("\\}", model_lines_3)[grep("\\}", model_lines_3) >= x][[1]]

      model_lines_3[number] <<- gsub("\\}", "",model_lines_3[number])
      model_lines_3[x] <<- gsub("\\(|\\)|\\{", "",model_lines_3[x])

      # if( gsub("[[:blank:]]", "", model_lines_3[number]) == "")
        # model_lines_3 <<- model_lines_3[-number]

    })


  grep("else if \\(", model_lines_3) %>%
    walk(function(x) {
      # print(x)
      # find the ending brackets of an if declaration
      number <- grep("\\}", model_lines_3)[grep("\\}", model_lines_3) >= x][[1]]
      # print(number)
      model_lines_3[number] <<- gsub("\\}", "",model_lines_3[number]) # add <<-
      model_lines_3[x] <<- gsub("else if", "elseif",model_lines_3[x])# add <<-
      model_lines_3[x] <<- gsub("\\(|\\)|\\{", "",model_lines_3[x])
      # if( gsub("[[:blank:]]", "", model_lines_3[number]) == "")
      # model_lines_3 <- model_lines_3[-number]# add <<-

    })



  grep("else *\\{", model_lines_3) %>%
    walk(function(x){

      # find the ending brackets of an else declaration
      number <- grep("\\}", model_lines_3)[grep("\\}", model_lines_3) >= x][[1]]

      model_lines_3[number] <<- gsub("\\}", "end\n",model_lines_3[number])
      model_lines_3[x] <<- gsub("else *\\{", "else",model_lines_3[x])

    })

  model_lines_3 <- model_lines_3 %>%
    gsub(pattern = "^ *",replacement =  "")

  blocs$equation2 <- paste0( "\n; Dynamical model
                             ",model_lines_3 %>%
                             paste0("\n") %>%
                             reduce(paste0))


if(is.data.frame(mb_output)) mb_output <- mb_output$output
  blocs$OUTPUT <- paste0("\n\nOUTPUT: \noutput = {", paste0(mb_output, collapse = ","), "}")

blocs <- map(blocs, ~ gsub("^ *", "", .x) %>%
               gsub(pattern = "\n *", replacement = "\n")) %>%
  gsub(pattern = "output_", replacement = "")

  if(outputcat == T){
    blocs %>%
      reduce(c) %>%
      cat(sep = "")
  }else{


    return( blocs %>%
              paste(collapse = "\n\n") )
  }

}

# shinyApp(pecc_ui_reunif, pecc_server_reunif)
# shinyApp(pecc_ui, pecc_server)


# #
# #
# Lorenz <- function(t, state, parameters) {
#   with(as.list(c(state, parameters)), {
#
#     if(t < tlag_input ){
#
#       growth <- 0
#       elimination <- 0
#
#     }else{
#
#       growth =   kgrowth #* (1 - TIME ** hillt50/ (T50 ** hillt50 +TIME ** hillt50))
#       elimination = UCART * kdeg * Lym ** hillL50/ (1 + Lym ** hillL50 +L50 ** hillL50)
#
#     }
#
#
#
#
#
#     KE=CL_input/V1_input
#     K12=Q_input/V1_input
#     K21=Q_input/V2_input
#
#
#     CPAL <- Alen1 / V1_input
#     CUCART <- UCART / VCART
#
#     dUCART <-  growth - elimination
#     dAlen1 =  K21*Alen2 - K12*Alen1 - KE*Alen1
#     dAlen2 =  - K21*Alen2 + K12*Alen1
#     dLym= KIN_input - KOUT_input*Lym*(1+ EFF_input*CPAL)
#
#     list(c(dUCART, dAlen1, dAlen2, dLym), c(Alemtuzumab = CPAL, UCART19 = CUCART, growth = growth, elimination = elimination))
#   })
# }
# #
# parameters <- c(tlag_input =7.81 + 7, kgrowth = 1790, T50 = 18.6, kdeg = 0.453, VCART = 77, hillt50 = 10, L50 = 0.01, hillL50 = 3, CL_input = 0.37, V1_input = 3.46, V2_input = 2.26, Q_input = 0.34, KIN_input = 0.0677, KOUT_input = 0.0410, EFF_input = 6.07 )
#
# state      <- c( UCART = 0, Alen1 = 0, Alen2 = 0, Lym = 0.5)
#
# times      <- seq(0, 100, by = 0.1)
#
# #
# #
# eventdat <- data.frame(var = c("Alen1", "Alen1", "Alen1", "Alen1",  "Alen1"),
#                        time = c(0, 1, 2, 3, 4) ,
#                        value = c(14, 14, 14, 14,14),
#                        method = c("add", "add", "add", "add","add"))
#
# out <- ode(y = state, times = times, func = Lorenz, parms = parameters,  events = list(data = eventdat)) %>%
#   as_tibble()

#
#
# out %>%
#   mutate(eliminationbrut = elimination / UCART19) %>%
#   gather(Alemtuzumab, UCART19, Lym, growth, elimination,  eliminationbrut, key = "key", value ="value") %>%
#   ggplot()+
#   geom_line(aes(time, value))+
#   facet_wrap(~key, scales = "free")+
#   theme_bw()
#
#
# ode_nonmem(y = state, times = times, func = Lorenz, parms = parameters,  events = list(data = eventdat))
# # ode_nonmem(y = state, times = times, func = Lorenz, parms = parameters,  events = list(data = eventdat), datafile = "file:///Z:/ETUDES/SPK/CLSPKPOOL/ANACIN/USERS/TDPE_CB_2/UCART19/Modeling/Cytometry_plasma/4_time_delay_linearG.cfl")
# ode_monolix(y = state, times = times, func = Lorenz, parms = parameters,  events = list(data = eventdat), outputcat = F)
