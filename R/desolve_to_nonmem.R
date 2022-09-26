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
# "5", "3", "5"), Distrib = c("logN", "Norm", "logN", "logN",
# "fix"), E = c("Esti", "Fix", "Input", "Esti", "Esti"))
#
# omega <- diag(parms$Value) %>%
#   as.data.frame()
# names(omega) <- parms$Param
# rownames(omega) <- parms$Param
# # omega[[1,1]] <- "0.3L"
#
# depot <- tibble(Proto = "1", var = "A1", time = "0",
#                       value = 50, method = "add", use = TRUE, delete = FALSE, ADM = "",
#                       F = FALSE, tlag = FALSE, Perf = "None")
#
#
# mb_output <- tibble(output = "test", YTYPE = NA_integer_, err_add = "0.1",
#                     err_prop = "0.3", export = TRUE, rm = FALSE)
#
# ode_nonmem(func = func, omega = omega, parms = parms, y = y, mb_output = mb_output )
#' DeSolve pharmacometrics model creator
#'
#' @description Create a NONMEM code with a deSolve model.
#' @param func Function (same as deSolve in functions)
#' @param parms Parameters (same as deSolve in functions)
#' @param y same as the ode
#' @param output vector of the output names. In case of several outputs, use the format c("name1 = X", "name2 = Y",...), with X/Y being the YTYPE number.
#' @param add_param vector of parameter names to add, mostly usefull for adding baseline parameters.
#' If it has the format "X1_0", with  X1 being the name of a cmt, add it as a baseline for the concerned cmt.
#' If the format is "Z = 21", add a new parameter with the name Z and initial value 21.
#' @param mu_referencig booléen, declaration of parameters with/without mu referencing.
#' @param BLQ
#' @param datafile Either a .tab or an other cfl
#' @param ...
#' @author Thibaud Derippe
#'
#'@return
#'@export


# mega = omega
# mu_referencig = F
# BLQ = F
# add_param = ""
# mb_output = mb_output
# y = states
# times = times
# func = eval(deSolve_peccc[[1]])
# outputcat = F
# parms = parameters
# datafile = ""

ode_nonmem <- function(func = Lorenz, omega, parms  = parameters, y = state, mb_output , add_param ="",    mu_referencig = T, BLQ = F, datafile = "", outputcat = T,  ...){

  if(is.character(func)) func <- deSolve_pecc(model_demo)$model

  blocs <- list()

  parms <- parms %>%
    filter(E != "Input") %>%
    rowid_to_column("theta") %>%
    # compute n of eta
    left_join(
      parms %>%
        filter(E != "Input") %>%
        filter(Distrib != "NoVar") %>%
        rowid_to_column("eta")
    )


  init <- y # y to have the same paramter name as ode,
  #but I hate to work with such names (e.g because of .y in imap)



  # # For handling both Adrien and my ways of coding
  # if(is.null(names(init))){
  #   namesinit <- paste0("X", 1:length(init))
  # } else{
  #   namesinit <-names(init)
  # }
  #
  # names(init) <- namesinit


  # If manually adding a parameter

  # if(add_param != ""){
  #
  #   walk(add_param, function(x){
  #
  #     ## If it is an initial value of a compartment
  #     if(gsub("_0$", "", x) %in% namesinit){
  #
  #       # Take the value used in deSolve
  #       valueparam <- init[names(init) == gsub("_0$", "", x)]
  #
  #       # Replaced it by the name of the parameter
  #       init[names(init) == gsub("_0$", "", x)] <<- x
  #
  #       # Add the paramater in the list
  #       parms[x] <- valueparam
  #       parms <<- parms
  #
  #     }else{
  #       ### Else it the user just want to add an other parameter (with or without value)
  #
  #
  #       xname <- gsub(" *=.+", "", x)
  #       xvalue <- if_else(grep("=",x) %>% length > 0 , gsub(".+= *", "", x), "")
  #
  #       parms[xname] <- xvalue
  #       parms <<- parms
  #
  #     }
  #   })
  #
  # }

  # Bloc PROB

  blocs$prob <- paste0("$PROB name of the model and/or description \n\n")


  # If we provide a cfl model, copy the $INPUT, $DATA and $SUBROUTINE

  if(gsub(".+\\.", "", datafile) == "cfl"){

    cfl_reference <- suppressWarnings(readLines(datafile))

    input_ref <- grep("\\$INPUT", cfl_reference)
    model_ref <- grep("\\$MODEL", cfl_reference)


    temp <- cfl_reference[input_ref:(model_ref - 1)] %>%
      paste0("\n\n")

    blocs$copy1 <-  temp[temp !="\n\n"]

  }else{

    # Else we do it manually
    # Bloc input

    if(datafile != ""){

      blocs$input <- paste0("$INPUT ", readLines(datafile, n = 1),"\n\n")

    }else{

      blocs$input <- "$INPUT\n\n"

    }

    # Bloc data

    if(datafile != ""){


      blocs$data <- paste0("$DATA ", gsub("(.+/)|(.+\\\\)/", "", datafile)," IGNORE = ",  gsub("(\\$INPUT)| ", "", blocs$input) %>%
                             substr(0,1),"\n\n")

    }else{

      blocs$data <- "$DATA your_file.data IGNORE = I\n\n"

    }


    # Bloc subroutine

    blocs$subroutine <- "$SUBROUTINE ADVAN13 TOL9\n\n"
  } ### end else


  # Bloc model


  blocs$model <- map_chr(init$Cmt, ~ paste0("COMP = (", .x, ")\n")) %>%
    reduce(paste0) %>%
    {paste0("$MODEL\n", ., "\n")}



  # Bloc PK

  #part1 Parameter definition with theta and eta attribution


  param_estimated <- parms$Value[parms$E != "Input"]
  names(param_estimated) <- parms$Param[parms$E != "Input"]


  if(mu_referencig == T){

    parms %>%
      mutate(init = map2_chr(theta, Distrib, ~
        case_when(.y == "logN" ~ paste0("MU_",.x," = LOG(THETA(", .x,"))" ),
                  .y %in% c("Norm","NoVar") ~ paste0("MU_",.x," = THETA(", .x,")" ),
                  T == "NoVar" ~ "")
      )) %>%
      mutate(eta = if_else(!is.na(eta), paste0("+ ETA(", eta,")"), "")) %>%
      mutate(param = pmap_chr(list(Param, theta, eta, Distrib), function(Param, theta, eta, Distrib){

        case_when(Distrib == "logN" ~ paste0(Param, " = EXP(MU_", theta, eta,")" ),
                  Distrib == "Norm" ~ paste0(Param, " = MU_", theta, eta ),
                  Distrib == "NoVar" ~ paste0(Param, " = MU_", theta ),
                  T ~ "")

      })) -> temp

    lines <-  paste( paste0(temp$init[!is.na(temp$init)], collapse = "\n"),
                     paste0(temp$param, collapse = "\n"), sep  = "\n\n")


  }else{

    parms %>%
      mutate(etal = map2_chr(Distrib, eta, ~case_when(.x == "logN" ~ paste0(" * exp(ETA(", .y,"))"),
                                                  .x == "Norm" ~ paste0(" + ETA(", .y,")"),
                                                  T ~ ""))) %>%
      mutate(line = paste0(Param, " = THETA(",theta,")", etal)) %>%
      pull(line) -> lines

  }
  blocs$pk1 <- paste0("$PK\n\n", paste0(lines, collapse = "\n"),"\n\n")
  # part2: initialisation of compartment

  blocs$pk2 <-
  init %>%
    rowid_to_column() %>%
    mutate(line = paste0("A_0(", rowid,") = ", t0, " ; ", Cmt, " initialisation")) %>%
    pull(line) %>%
    paste0(collapse = "\n")

  # BLoc DES


  model_lines <- deparse(func,  width.cutoff = 500)%>%
    gsub(pattern = "/",replacement =  " / ")


  ##### replacement in case of PFIM / Adrien style of writing

  if(   grep("with\\( *as\\.list *\\(", model_lines)  %>% length == 0){



    model_lines <-  model_lines[- grep("<- *p\\[", model_lines)]

    model_lines <- model_lines %>%
      gsub(pattern = "y\\[", replacement = "X") %>%
      gsub(pattern = "\\]", replacement = "") %>%
      gsub(pattern = "yd", replacement = "dX")

    temp_to_sup <- grep("\\{|\\}", model_lines)
    model_lines <-  model_lines[-c(1:temp_to_sup[[1]], temp_to_sup[[length(temp_to_sup)]]:length(model_lines))]

  }else{

    temp_to_sup <- grep("\\{|\\}", model_lines)
    model_lines <-  model_lines[-c(1:temp_to_sup[[2]], temp_to_sup[[length(temp_to_sup)-1]]:length(model_lines))]


  }

  ######

  model_lines <- model_lines[- c(grep("list\\(c\\(", model_lines): length(model_lines))]

  model_lines_2 <-
    model_lines %>%
    gsub(pattern = "<-", replacement = "=") %>%
    gsub(pattern = " *<= *", replacement = ".LE.") %>%
    gsub(pattern = " *>= *", replacement = ".GE.") %>%
    gsub(pattern = " *< *", replacement = ".LT.") %>%
    gsub(pattern = " *> *", replacement = ".GT.") %>%
    gsub(pattern = " *& *", replacement = ".AND.") %>%
    gsub(pattern = " *\\| *", replacement = ".OR.") %>%
    gsub(pattern = "^ *if *\\(", replacement = "\nIF(") %>%
    gsub(pattern = "\\^", replacement = " ** ") %>%
    gsub(pattern = "(?<![[:alnum:]])t(?![[:alnum:]])", replacement = "T", perl = T)    ## if reference to t as time

  model_lines_3 <- model_lines_2 %>%
    paste0(" ") # I need this space for the next paste0(" ",x, " ")




  for(x in init$Cmt){

    number <- which(init$Cmt == x)

    model_lines_3 <- gsub(paste0(" *d", x), paste0("DADT(", number, ") "), model_lines_3) %>%
      gsub(pattern = paste0(" \\(?",x, "\\)? "), replacement = paste0(" A(", number, ") "))

  }

  ## gestion of if else structure
  # x = 14 #  33  48  57  78  90  96 104 119

  grep("IF\\(", model_lines_3) %>%
    walk(function(x) {
# print(x)
      # find the ending brackets of an if declaration
      number <- grep("\\}", model_lines_3)[grep("\\}", model_lines_3) >= x][[1]]
# print(number)
      model_lines_3[number] <<- gsub("\\}", "",model_lines_3[number]) # add <<-
      model_lines_3[x] <<- gsub("\\{", " THEN",model_lines_3[x])# add <<-

      # if( gsub("[[:blank:]]", "", model_lines_3[number]) == "")
        # model_lines_3 <- model_lines_3[-number]# add <<-

    })

  grep("else if \\(", model_lines_3) %>%
    walk(function(x) {
      # print(x)
      # find the ending brackets of an if declaration
      number <- grep("\\}", model_lines_3)[grep("\\}", model_lines_3) >= x][[1]]
      # print(number)
      model_lines_3[number] <<- gsub("\\}", "",model_lines_3[number]) # add <<-
      model_lines_3[x] <<- gsub("else if", "ELSEIF",model_lines_3[x])# add <<-
      model_lines_3[x] <<- gsub("\\{", " THEN",model_lines_3[x])
      # if( gsub("[[:blank:]]", "", model_lines_3[number]) == "")
      # model_lines_3 <- model_lines_3[-number]# add <<-

    })

  grep("else *\\{", model_lines_3) %>%
    walk(function(x){

      # find the ending brackets of an else declaration
      number <- grep("\\}", model_lines_3)[grep("\\}", model_lines_3) >= x][[1]]

      model_lines_3[number] <<- gsub("\\}", "ENDIF\n",model_lines_3[number])
      model_lines_3[x] <<- gsub("else *\\{", "ELSE",model_lines_3[x])

    })


  model_lines_3 <- gsub(" *== *", ".EQ.", model_lines_3)

  blocs$DES <- c("\n\n$DES\n\n",
                 model_lines_3 %>%
                   gsub(pattern = "^ *", replacement = "") %>%
                   paste0("\n")
  )

  # Bloc error



  #### First part: attributing theta and if needed output to ytype

  if(nrow (mb_output) == 1){

    nparm <- nrow(parms)
    # add parameter
    if(mb_output$err_add >0){

      addline <- paste0( "ADD = THETA(", nparm + 1,")")
      parms <- bind_rows(parms, tibble(theta = nparm + 1, Param = "err_add", Value = mb_output$err_add))

      nparm <- nparm + 1
      if(mb_output$err_prop  == 0) wline <- "W = ADD"
    }else{
      addline <- ""
    }


    # prop parameter
    if(mb_output$err_prop >0){

      propline <- paste0( "PROP = THETA(", nparm + 1,")")
      parms <- bind_rows(parms, tibble(theta = nparm + 1, Param = "err_prop", Value = mb_output$err_prop))
      if(mb_output$err_add  == 0) wline <- "W = (OUTPUT * OUTPUT * PROP ** 2)**0.5"
    }else{

      propline <- ""
    }

    # right equation
    if(mb_output$err_add > 0 & mb_output$err_prop > 0) wline <- "W = (OUTPUT * OUTPUT * PROP**2 + ADD**2)**0.5 + 0.0001"

    blocs$error1 <-    paste0("\n$ERROR
                              ", addline,"
                              ", propline,"
                              OUTPUT =", mb_output$output ,"
                              IPRED = OUTPUT
                              ", wline,"
                              IRES = DV - IPRED
                              IWRES = (DV - IPRED)/(W)\n\n") %>%
      gsub(pattern = "  *", replacement = " ") %>%
      gsub(pattern = "\n ", replacement = "\n")


  }else{
    # In case of several YTYPE
    nparm <- nrow(parms)
    blocs$error1 <- map_chr(mb_output$output, function(x){

      xname <- gsub(" *=.+", "", x)
      xvalue <- mb_output$YTYPE[mb_output$output == x]
      addvalue <- mb_output$err_add[mb_output$output == x]
      propvalue <- mb_output$err_prop[mb_output$output == x]

      if(addvalue >0){

        addline <- paste0( "ADD = THETA(", nparm + 1,")")
        parms <<- bind_rows(parms, tibble(theta = nparm + 1, Param = paste0("err_add", xname), Value = addvalue))

        nparm <<- nparm + 1

        if(propvalue  == 0) wline <- "W = ADD"
      }else{

        addline <- ""
      }


      # prop parameter
      if(propvalue >0){

        propline <- paste0( "PROP = THETA(", nparm + 1,")")
        parms <<- bind_rows(parms, tibble(theta = nparm + 1, Param = paste0("err_prop", xname), Value = propvalue))

        nparm <<- nparm + 1
        if(addvalue  == 0) wline <- "W = (OUTPUT * OUTPUT * PROP ** 2)**0.5"
      }else{

        propline <- ""
      }

      # right equation
      if(addvalue > 0 & propvalue > 0) wline <- "W = (OUTPUT * OUTPUT * PROP**2 + ADD**2)**0.5 + 0.0001"

      blocs$error1 <-    paste0("\n$ERROR
                              ", addline,"
                              ", propline,"
                              OUTPUT =", mb_output$output ,"
                              IPRED = OUTPUT
                              ", wline,"
                              IRES = DV - IPRED
                              IWRES = (DV - IPRED)/(W)\n\n") %>%
        gsub(pattern = "  *", replacement = " ") %>%
        gsub(pattern = "\n ", replacement = "\n")




      paste0("IF(YTYPE.EQ.", xvalue ,") THEN\n", addline,"\n",
             propline,"
             OUTPUT = ", xname,"\n",  wline, "\nENDIF\n\n") -> temp


      return(temp)

    }) %>%
      reduce(paste0) %>%
      {paste0("\n$ERROR\n\n", ., "IPRED = OUTPUT
              IRES = DV - IPRED
              IWRES = (DV - IPRED)/(W)\n\n")  } %>%
      gsub(pattern = "  *", replacement = " ") %>%
      gsub(pattern = "\n ", replacement = "\n")

      }

  ### Controling or not the BLQ

  if(BLQ == F){

    blocs$error2 <- "Y = OUTPUT +  W * EPS(1)"

  }else{

    blocs$error2 <-  ";------------------BLQ handling with M3 METHOD --------------;
    DUM =(LLOQ - IPRED)/W
    CUMD = PHI(DUM)

    IF(BLQ.EQ.0.OR.NPDE_MODE.EQ.1) THEN
    F_FLAG = 0
    Y = OUTPUT +  W * EPS(1)
    ENDIF

    IF(BLQ.EQ.1.AND.NPDE_MODE.EQ.0) THEN
    F_FLAG = 1
    Y = CUMD
    MDVRES = 1
    IRES = 0
    IWRES = 0
    ENDIF

    IF(BLQ.EQ.1) DV_LOQ = LLOQ"


  }




  # Bloc Theta
blocs$THETA1 <- parms %>%
  mutate(E = as.character(E)) %>%
  mutate(E = case_when(is.na(Distrib) & length(grep("F", Value)) == 1 ~ "Fix",
                       is.na(Distrib) & length(grep("F", Value)) == 0 ~ "Esti",
                       T ~ E)) %>%
  mutate(E = gsub("F", "", E)) %>%
  mutate(test = if_else(E == "Esti", paste0("(0,",Value,") ; ",Param, " - theta ",theta),
                        paste0(Value," FIX ; ",Param, " - theta ",theta))) %>%
  {paste0("\n\n$THETA\n", paste0(.$test, collapse = "\n"))}



# bloc Omega sigma and cov
blocs$OMEGA  <- parms %>%
  filter(!is.na(eta)) %>%
  left_join(
tibble(Param = names(omega), value = diag(omega %>% as.matrix))
) %>%
  mutate(value = gsub("F$", " FIX", value)) %>%
  mutate(value = gsub(" *0 *$", "0 FIX", value)) %>%
  mutate(test = paste0(value," ; ", Param, " - eta ", eta))%>%
  {paste0("\n\n$OMEGA\n", paste0(.$test, collapse = "\n"))}

# If we provide a cfl model, copy the $INPUT, $DATA and $SUBROUTINE

  if(gsub(".+\\.", "", datafile) == "cfl"){


    sigma_ref <- grep("\\$SIGMA", cfl_reference)

    temp <- cfl_reference[sigma_ref:length(cfl_reference)] %>%
      paste0("\n\n") %>%
      gsub(pattern = "FILE *= *.+\\.(TAB)|(tab)", replacement = "FILE = youroutputfile.TAB")

    blocs$copy2 <-  temp[temp !="\n\n"]


  }else{

    blocs$sigma_cov <- "\n\n$SIGMA 1 FIX \n\n$COV  MATRIX=S"

    # blocs$EST

    blocs$EST <- "\n\n$EST METHOD=SAEM LAPLACIAN INTERACTION NBURN = 1000 NITER = 300 ISAMPLE = 2 SIGL = 10 CTYPE = 3"


    # bloc table

    blocs$table <- "\n\n$TABLE ID	TIME DV IPRED ONEHEADER NOAPPEND NOPRINT FILE = youroutputfile.TAB"
  }#end else


  if(outputcat == T){
    blocs %>%
      reduce(c) %>%
      cat(sep = "")
  }else{


    return( blocs %>%
              reduce(c) %>%
              reduce(paste0))
  }
  }


