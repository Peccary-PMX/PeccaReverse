

# table <- output$res[1:2] %>%
#   left_join(thetas %>%
#               rename(Param = theta_names)  %>%
#               select(Param, distrib))

# Table des paramètres ----------------------------------------------------

#' nonmem_to_desolve
#'@export
table_param <- function(table = tibble(Param = character(), Value = character())){

names(table)[1:2] <- c("Param", "Value")

# Add Distrib
if(ncol(table) >=3){

  names(table)[3] <- "Distrib"

  table[[3]] <- factor(table[[3]], levels = c("logN", "Norm", "NoVar"))
  table[[3]][is.na(table[[3]])] <- "logN"
}else{

  table$Distrib <- factor("logN", levels = c("logN", "Norm", "NoVar"))

}

# Add Estimation
if(ncol(table) >=4){


  names(table)[4] <- "E"

  table[[4]] <- factor(table[[4]], levels = c("Esti", "Fix", "Input"))
  table[[4]][is.na(table[[4]])] <- "Esti"
}else{

table$E <- factor("Esti", levels = c("Esti", "Fix", "Input"))

}

return(table[1:4])

}




# Table des inputs --------------------------------------------------------

# var = factor(output$initial_values$Cmt[[1]], levels = output$initial_values$Cmt  )
# table = inputbase %>%
#   left_join(
#
#     df %>%
#       filter(eventidentifier == 1, identifier == reference) %>%
#       select(time, administration, amount) %>%
#       rename(adm = administration) %>%
#       mutate(adm  = as.character(adm))
#   );table

#' nonmem_to_desolve
#'@export
table_input <- function(table = NULL, var = NULL){


  if(is.null(table) & ! is.null(var)){

    return(
      tibble(Proto = c("1", "")) %>%
      mutate(var = var) %>%
      mutate(time = c("0", ""), #
             value = c(0.0, NA), #
             method = factor(c("add", NA), levels = c("add", "mult", "rep")),
             use = c(T,F),
             delete = c(F,F),
             ADM = c("1", ""),
             filterPlot = "",
             F = c(F,F),
             tlag = c(FALSE,FALSE),
             Perf = factor("None", levels = c("None", "rate","time"))
      )
      )

  }else if(!is.null(table) & is.null(var)){

    if(sum(names(table) %in% "Proto") == 0) table$Proto <-"1"
    names(table)[names(table) %in% c("target")] <- "var"

    names(table)[tolower(names(table)) %in% c("time", "t")] <- "time"
    if(sum(names(table) %in% "time") == 0) table$time <-"0"
    table$time <- as.character(table$time)

    names(table)[tolower(names(table)) %in% c("amt", "amount")] <- "value"
    if(sum(names(table) %in% "value") == 0) table$value <- 0

    if(sum(names(table) %in% "method") == 0){
      table$method <- factor("add", levels = c("add", "mult", "rep"))
    }else{
      table$method <- factor(table$method, levels = c("add", "mult", "rep"))
      table$method[is.na(table$method)] <- "add"
      }

    if(sum(names(table) %in% "use") == 0) table$use <- TRUE
    table$use <- as.logical(table$use)
    if(sum(names(table) %in% "delete") == 0) table$delete <- FALSE
    table$delete <- as.logical(table$delete)

  names(table)[tolower(names(table)) %in% c("adm")] <- "ADM"
  if(sum(names(table) %in% "ADM") == 0) table$ADM <- "1"
  table$ADM <- as.character(table$ADM)

  if(sum(names(table) %in% "F") == 0) table$F <- FALSE
  table$F <- as.logical(table$F)
  if(sum(names(table) %in% "tlag") == 0) table$tlag <- FALSE
  table$tlag <- as.logical(table$tlag)
  if(sum(names(table) %in% "Perf") == 0){
    table$Perf <- factor("None", levels = c("None", "rate","time"))
  }else{

    table$Perf <- factor(table$Perf, levels = c("None", "rate","time"))
    table$Perf[is.na(table$Perf)] <- "None"
  }

  if(! "filterPlot" %in% names(table)){

    table$filterPlot <- ""

  }
  table$filterPlot <- as.character(table$filterPlot)

 table <-  table %>%
    select(Proto, var, time, value, method, use, delete, ADM, filterPlot, F, tlag, Perf, contains("Perf_num"))

return(table)
  }
  #Perf_num

}

# Table des display --------------------------------------------------------
#' nonmem_to_desolve
#'@export
table_display <- function(table = NULL, possiblevalues = NA_character_){

  if(!is.null(table)){

    if(! "Plot" %in% names(table)) table$Plot <- 1L
    if(class(table$Plot) != "integer") table$Plot <- as.integer(table$Plot)
    if(! "Todisplay" %in% names(table)) table$Todisplay <- possiblevalues
    if(! "Check" %in% names(table)) table$Todisplay <- TRUE
    if(! "Point" %in% names(table)) table$Todisplay <- FALSE
    if(! "Filter_of_dataset" %in% names(table)) table$Todisplay <- NA_character_
    if(! "YTYPE" %in% names(table)) table$YTYPE <- NA_character_



    table %>%
      select(Plot, Todisplay, Check, Point, Filter_of_dataset, YTYPE)

  }else{
  tibble(Plot = 1L,
         Todisplay = possiblevalues,
         Check = T, Point = F, Filter_of_dataset = "", YTYPE = "")
  }
}
