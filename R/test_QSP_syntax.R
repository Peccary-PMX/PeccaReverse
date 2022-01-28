# deg_param <- c("Bcl2", "Bclxl", "Mcl1", "BIM", "tBID", "PUMA", "NNOXA", "Bcl2-BIM", "Bclxl-BIM",
#                "Mcl1-BIM", "Bcl2-tBID", "Bclxl-tBID", "Mcl1-tBID",  "Mcl1-tBID", "Bck2-PUMA", "Bclxl-PUMA",
# "Mcl1-PUMA")
#
# params <- read.table("file:///D:/these/Second_project/model_lidner.csv", header = T, sep = ";",stringsAsFactors = F) %>%
#   as_tibble %>%
#   mutate(kdeg = kdeg * 3600, kforward = kforward * 3600, kbackward = kbackward * 3600, k = k * 3600)
# #
#
#
#
# test_per_bloc(params)

#' QSP to Peccary
#' @description
#'
#' @param params
#' @export
#'

QSP_to_Peccary <- function(params, outputExpr = T){

  model_bloc <- character()
  model_deriv <- character()
  parameters_output <- tibble(Parameter = character(), Value = character())
  compt <- compteur()
  init_value <- character()

  params$kdeg <- as.character(params$kdeg)
  # a = 35
  # a <- a + 1
  for(a in 1:nrow(params)){ #nrow(params)
    # print(a)
    line_input <- params %>% slice(a)

    compounds <- gsub("(^ *)|( *$)", "", line_input$Entity)

    # Handling initial value

    if(!is.na(line_input$init_value) & line_input$init_value != ""){

      init_value <- c(init_value, paste0("d", compounds,"_0 <- ", line_input$init_value))

    }


    # Handling degradation

    if(!is.na(line_input$kdeg) & line_input$kdeg !=""){ #length(str_split(compounds, " |\\+")[[1]]) == 1


      nline <- grep(paste0("^ *d", compounds," +"), model_deriv)
      newparam <- paste0("kdeg_", compounds)
      newbloc <-  paste0("- ", newparam, " * ", compounds)
      namenewbloc <- paste0("R", compt(), "_deg_",compounds)

      model_bloc <- c(model_bloc, paste0(namenewbloc, " <- ",newbloc ))

      if(length(nline) == 0){

        model_deriv <- c(model_deriv, paste0("d", compounds, " <- ", namenewbloc))

      }else{

        model_deriv[nline] <- paste0(model_deriv[nline], " + ", namenewbloc)

      }

      parameters_output <- bind_rows(parameters_output,
                                     tibble(Parameter = newparam, Value = line_input$kdeg))

    }

    # Handling production

    if(!is.na(line_input$kprod) & line_input$kprod !=""){ #length(str_split(compounds, " |\\+")[[1]]) == 1


      nline <- grep(paste0("^ *d", compounds," +"), model_deriv)
      newparam <- paste0("kprod_", compounds)

      if(is.null(line_input$order_kprod)) line_input$order_kprod <- 0

      if(line_input$order_kprod == 0){
      newbloc <-  paste0( newparam)
      }else{
        newbloc <-  paste0( newparam)
      }
      namenewbloc <- paste0("R", compt(), "_prod_",compounds)

      model_bloc <- c(model_bloc, paste0(namenewbloc, " <- ",newbloc ))

      if(length(nline) == 0){

        model_deriv <- c(model_deriv, paste0("d", compounds, " <- ", namenewbloc))

      }else{

        model_deriv[nline] <- paste0(model_deriv[nline]," + ", namenewbloc)

      }

      parameters_output <- bind_rows(parameters_output,
                                     tibble(Parameter = newparam, Value = line_input$kprod))

    }


    # Handling kforward and kbackword
    if(!is.na(line_input$kforward) & line_input$kforward != ""){


      cpd_split <- str_split(compounds, " ")[[1]]

      if(length(cpd_split) & grepl("_", cpd_split)){
        # to handle input such as "Bclxl_BAXc

        complexx <- compounds
        substrate <- str_split(complexx, "_")[[1]]
      }else{
        # to handle input such as "BAK2 + BAK2 -> BAK4"
        splitt <- str_split(compounds, " *-> *")[[1]]

        complexx <- splitt[[2]]
        substrate <- str_split(splitt[[1]], " *\\+ *")[[1]]

      }

      # handling complex lines
      nline <- grep(paste0("^ *d", complexx," +"), model_deriv)
      newparamforward <- paste0("kforward_", complexx)
      newparambackward <- paste0("kbackward_", complexx)
      newbloc <-  paste0(newparamforward, " * ", paste0(substrate, collapse = " * "),
                         " - ",  newparambackward," * ", complexx)
      namenewbloc <- paste0("R", compt(), "_complex_",complexx)
      model_bloc <- c(model_bloc, paste0(namenewbloc, " <- ",newbloc ))

      if(length(nline) == 0){

        model_deriv <- c(model_deriv, paste0("d", complexx, " <- ", namenewbloc))

      }else{

        model_deriv[nline] <- paste0(model_deriv[nline], " + ", namenewbloc)

      }

      # handling substrate lines
      for(b in substrate){

        nline <- grep(paste0("^ *d", b," +"), model_deriv)


        if(length(nline) == 0){

          model_deriv <- c(model_deriv, paste0("d", b, " <- - ", namenewbloc))

        }else{

          model_deriv[nline] <- paste0(model_deriv[nline]," - ", namenewbloc)

        }

      }

      parameters_output <- bind_rows(parameters_output,
                                     tibble(Parameter = c(newparamforward, newparambackward),
                                            Value = as.character(c(line_input$kforward, line_input$kbackward))))



    } # end handling kforward

    # Handle k

    if(!is.na(line_input$k) & line_input$k != ""){

      compounds <- gsub("(^ *)|( *$)", "", line_input$Entity)
      cpd_split <- str_split(compounds, " ")[[1]]


      # to handle input such as "BAK2 + BAK2 -> BAK4"
      splitt <- str_split(compounds, " *-> *")[[1]]

      left <- str_split(splitt[[1]], " *\\+ *")[[1]]
      right <- str_split(splitt[[2]], " *\\+ *")[[1]]



      newparam <- paste0("k_", paste0(gsub(" *\\+ *","", left),collapse = "_"))

      parameters_output <- bind_rows(parameters_output,
                                     tibble(Parameter = c(newparam),
                                            Value =  as.character(c(line_input$k))))

      # note: it won't work if we have a second order rate -> need to be improved in the futur
      newbloc <- paste0(" - ", newparam, " * ", left)
      namenewbloc <- paste0("R", compt(), "_disso_",left)
      model_bloc <- c(model_bloc, paste0(namenewbloc, " <- ",newbloc ))

      for(b in c(left, right)){

        nline <- grep(paste0("^ *d", b," +"), model_deriv)


        if(b %in% left){

          newbloc <-  paste0(" + ", namenewbloc)

        }else{

          newbloc <-  paste0(" - ",namenewbloc)

        }


        if(length(nline) == 0){

          model_deriv <- c(model_deriv, paste0("d", b, " <- ", newbloc))

        }else{

          model_deriv[nline] <- paste0(model_deriv[nline], newbloc)

        }

      }


    }



  } # end "for each line" loop

  model_bloc  <- model_bloc %>%
    gsub(pattern = "<- *\\+ *", replacement = "<- ")

  model_deriv <- model_deriv %>%
    gsub(pattern = "<- *\\+ *", replacement = "<- ")

  parameters_output %>%
    mutate(Value = Value) %>%
    mutate(out = paste0(Parameter, " <- ", Value)) %>%
    pull(out) %>%
    paste0(collapse = "\n") -> allvalues


  paste0("####### Fixed values #######\n",allvalues,
         "\n####### End fixed values #######\n\n\n####### Initial conditions #######\n",  paste0(init_value, collapse = "\n"),
         "\n####### End initial conditions #######\n\n\n####### All rates #######\n",
         paste0(model_bloc, collapse = "\n\n"),
         "\n####### End All rates #######\n\n\n####### All ODES #######\n",
         paste0(model_deriv, collapse = "\n\n")) -> output


  if(outputExpr)return(output %>% cat)

  return(output)

}



test <- function(params){

  model_deriv <- character()
  parameters_output <- tibble(Parameter = character(), Value = double())
a = 17
  a <- a + 1
  for(a in 1:nrow(params)){ #nrow(params)
# print(a)
  line_input <- params %>% slice(a)

  compounds <- gsub("(^ *)|( *$)", "", line_input$Entity)

  # Handling degradation

  if(!is.na(line_input$kdeg)){ #length(str_split(compounds, " |\\+")[[1]]) == 1

    nline <- grep(paste0("^ *d", compounds," +"), model_deriv)
    newparam <- paste0("kdeg_", compounds)

    newbloc <-  paste0("- ", newparam, " * ", compounds)
    if(length(nline) == 0){

      model_deriv <- c(model_deriv, paste0("d", compounds, " <- ", newbloc))

    }else{

      model_deriv[nline] <- paste0(model_deriv[nline], newbloc)

    }

    parameters_output <- bind_rows(parameters_output,
                                   tibble(Parameter = newparam, Value = line_input$kdeg))

  }

  # Handling kforward and kbackword
  if(!is.na(line_input$kforward)){


    cpd_split <- str_split(compounds, " ")[[1]]

    if(length(cpd_split) & grepl("_", cpd_split)){
      # to handle input such as "Bclxl_BAXc

      complexx <- compounds
      substrate <- str_split(complexx, "_")[[1]]
    }else{
      # to handle input such as "BAK2 + BAK2 -> BAK4"
      splitt <- str_split(compounds, " *-> *")[[1]]

      complexx <- splitt[[2]]
      substrate <- str_split(splitt[[1]], " *\\+ *")[[1]]

      }

    # handling complex lines
    nline <- grep(paste0("^ *d", complexx," +"), model_deriv)
    newparamforward <- paste0("kforward_", complexx)
    newparambackward <- paste0("kbackward_", complexx)
    newbloc <-  paste0(newparamforward, " * ", paste0(substrate, collapse = " * "),
                       " - ",  newparambackward," * ", complexx)

    if(length(nline) == 0){

      model_deriv <- c(model_deriv, paste0("d", complexx, " <- ", newbloc))

    }else{

      model_deriv[nline] <- paste0(model_deriv[nline], " + ", newbloc)

    }

    # handling substrate lines
    for(b in substrate){

      nline <- grep(paste0("^ *d", b," +"), model_deriv)

      newbloc <-  paste0(" - ", newparamforward, " * ", paste0(substrate, collapse = " * "),
                         " + ",  newparambackward," * ", complexx)

      if(length(nline) == 0){

        model_deriv <- c(model_deriv, paste0("d", b, " <- ", newbloc))

      }else{

        model_deriv[nline] <- paste0(model_deriv[nline], newbloc)

      }

    }

    parameters_output <- bind_rows(parameters_output,
                                   tibble(Parameter = c(newparamforward, newparambackward),
                                          Value = c(line_input$kforward, line_input$kbackward)))



  } # end handling kforward

  # Handle k

  if(!is.na(line_input$k)){

    compounds <- gsub("(^ *)|( *$)", "", line_input$Entity)
    cpd_split <- str_split(compounds, " ")[[1]]


      # to handle input such as "BAK2 + BAK2 -> BAK4"
      splitt <- str_split(compounds, " *-> *")[[1]]

      left <- str_split(splitt[[1]], " *\\+ *")[[1]]
      right <- str_split(splitt[[2]], " *\\+ *")[[1]]



    newparam <- paste0("k_", paste0(gsub(" *\\+ *","", left),collapse = "_"))

    parameters_output <- bind_rows(parameters_output,
                                   tibble(Parameter = c(newparam),
                                          Value = c(line_input$k)))

    # note: it won't work if we have a second order rate -> need to be improved in the futur

    for(b in c(left, right)){

      nline <- grep(paste0("^ *d", b," +"), model_deriv)


      if(b %in% left){

        newbloc <-  paste0(" - ", newparam, " * ", left)

      }else{

        newbloc <-  paste0(" + ", newparam, " * ", left)

      }


      if(length(nline) == 0){

        model_deriv <- c(model_deriv, paste0("d", b, " <- ", newbloc))

      }else{

        model_deriv[nline] <- paste0(model_deriv[nline], newbloc)

      }

    }


  }



  } # end "for each line" loop

  model_deriv[grep("dBIM <- ", model_deriv)]

  model_deriv <- model_deriv %>%
    gsub(pattern = "<- *\\+ *", replacement = "<- ")

  parameters_output %>%
    mutate(Value = Value * 3600) %>%
    mutate(out = paste0(Parameter, " <- ", Value)) %>%
    pull(out) %>%
    paste0(collapse = "\n") -> allvalues


  paste0(allvalues, "\n\n\n", paste0(model_deriv, collapse = "\n\n")) -> output

if(outputExpr)return(output %>% cat)

  return(output)

    }
# test(params)

