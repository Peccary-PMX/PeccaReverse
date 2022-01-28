# path_nonmem = "file:///D:/Peccary/Exemple_demo/run_nonmem/2_4_StandardErrorSansW12.cfl"
# path_nonmem = "file:///D:/Peccary/Exemple_demo/run_nonmem/2_4_StandardErrorSansW12.res"
path_nonmem = "file:///D:/Peccary/Exemple_demo/run_nonmem/1_comp.res"
# path_nonmem = "file:///D:/Peccary/Exemple_demo/run_nonmem/cpt_sn38_1cpt.res"
# path_nonmem <- "file:///D:/Peccary/Exemple_demo/ivivc.res"
# output <-
# nonmem_to_desolve(path_nonmem)
# output$model %>% cat
#' nonmem_to_desolve
#'@export
nonmem_to_desolve <- function(path_nonmem){


  run <- try(createRun(path_nonmem))

  output <- list()
  ### Find model file from mxtran

  path <-  path_nonmem %>%
    gsub(pattern = "file:///", replacement = "")

  root <- str_split(path, "/")[[1]]
  root <- root[- length(root)] %>%
    paste0(collapse = .Platform$file.sep)
  # read files
  lines <- readLines(path)

  # Get df ------------------------------------------------------------------

  # Obtain name of dataset
  datasetname <-
    lines[grep("\\$DATA", lines)] %>%
    gsub(pattern = "\\$DATA *", replacement = "") %>%
    gsub(pattern = " .*", replacement = "")

  # Read dataset (TODO)

  pathdataset <- file.path(root, "DATA", datasetname)

  # Get input
  if(file.exists(pathdataset)){

    # guess the separator
    firstline <- readLines(pathdataset, n = 1)
    nvirgule <- length(str_split(firstline, ",")[[1]])
    npointvirgule <- length(str_split(firstline, ";")[[1]])

    sepp <- case_when(nvirgule > 1 ~ ",",
                    npointvirgule > 1 ~ ";",
                    T ~ " ")

    df <- read.table(pathdataset, header = T, sep = sepp, na.strings = ".") %>%
      as_tibble

    lines[grep("\\$INPUT", lines):(grep("\\$DATA", lines)-1)] %>%
      gsub(pattern = "(\\$INPUT)|(;.+)", replacement = "") %>%
      paste0(collapse = " ") %>%
      gsub(pattern = " *= *", replacement = "=") %>%
      gsub(pattern = "  *", replacement = " ") -> headerr

    str_split(headerr, pattern = " ")[[1]]-> headerr
    headerr <- headerr[- which(headerr == "")]

   gsub(".+=", "", headerr) -> headerr

   # if(length(df) > length(headerr)) headerr <- c(headerr, rep("DROP", length(df) - length(header)))

   # names(df) <- headerr

   df <- df[names(df) %in% headerr]
   # df <- df[- which(names(df) == "DROP")]

   names(df)[toupper(names(df)) == "ID"] <- "identifier"
   names(df)[toupper(names(df)) == "TIME"] <- "time"
   names(df)[toupper(names(df)) == "DV"] <- "observation"
   names(df)[toupper(names(df)) == "AMT"] <- "amount"
   names(df)[toupper(names(df)) == "CMT"] <- "administration"
   names(df)[toupper(names(df)) == "EVID"] <- "eventidentifier"
   names(df)[toupper(names(df)) == "YTYPE"] <- "observationtype"
   names(df)[toupper(names(df)) == "MDV"] <- "missingdependentvariable"

   output$df <- df
  }

# names eta and theta -----------------------------------------------------


  # get names (copy from setReplaceNonmemrun -> make a separate fuction?)--------------------------------

  # eta
  eta_name <- function(res_file = lines){

    #### F I R S T   M E T H O D:   A T   T I M E    O F    D E C L A R A T I O N
    eta_name_1st_mtd <- character()
    a <- 1
    temp <- "initiation"
    res_file_commentless <- gsub(";.+", "", res_file)


    while(class(temp) != "try-error"){
      temp <- try(res_file_commentless[grep(paste0("[^H]ETA(.",a,")"), toupper(res_file_commentless))[[1]]], silent = T)
      value <-  try(gsub(" ", "",str_split(temp, "=")[[1]][1]), silent = T)
      if(class(temp) != "try-error"){
        eta_name_1st_mtd<- c(eta_name_1st_mtd, value)
      }
      a <- a + 1
    }

    for (a in eta_name_1st_mtd){
      if(length(eta_name_1st_mtd[eta_name_1st_mtd == a]) > 1){
        eta_name_1st_mtd[eta_name_1st_mtd == a] <- ""
      }
    }
    ### S E C O N D   M E T H O D:  C O M M E N T S   N A M E

    eta_begin <- grep("^[[:blank:]]?\\$OMEGA", toupper(res_file))[[1]]
    eta_end <- grep("^[[:blank:]]?\\$SIGMA", toupper(res_file))[[1]]
    res_file_shorten <- res_file[eta_begin:(eta_end-1)]
    res_file_shorten <- res_file_shorten[ - (toupper(res_file_shorten) == "$OMEGA")]
    res_file_shorten <- res_file_shorten[-grep("^[[:space:]]*;", res_file_shorten)]
    res_file_shorten <- res_file_shorten[res_file_shorten != ""]
    res_file_shorten <- gsub(".+;", "", res_file_shorten)
    res_file_shorten <- gsub("^[[:space:]]*", "", res_file_shorten)
    res_file_shorten <- gsub("[[:space:]].*","",res_file_shorten)

    for (a in res_file_shorten){
      if(length(res_file_shorten[res_file_shorten == a]) > 1){
        res_file_shorten[res_file_shorten == a] <- ""
      }
    }

    ### C O M P A R I S O N    A N D    F I N A L     I N P U T

    eta_name <- character()
    for (a in 1:length(eta_name_1st_mtd)){


      test1 <- toupper(eta_name_1st_mtd[a])
      test2 <- toupper(res_file_shorten[a])

      if(is.na(test1) == T) test1 <- ""
      if(is.na(test2) == T) test2  <- ""

      if(test1 == test2){
        if(test1 == "" & test2 ==""){
          eta_name <- c(eta_name, a)
        }else{
          eta_name <- c(eta_name, eta_name_1st_mtd[a])
        }
      }else if(test1 == "" | test2 ==""){

        eta_name <- c(eta_name, paste0(test1,test2))

      }else{

        eta_name <- c(eta_name, paste0(test1, "/",test2))
      }


    }
    return(paste0("eta_",eta_name))

  }
  eta_names <- eta_name()#; eta_names

  theta_name <- function(res_file = lines){

    #### F I R S T   M E T H O D:   A T   T I M E    O F    D E C L A R A T I O N
    theta_name_1st_mtd <- character()
    a <- 1
    temp <- "initiation"
    res_file_commentless <- gsub(";.+", "", res_file)


    while(class(temp) != "try-error"){
      temp <- try(res_file_commentless[grep(paste0("THETA\\( *",a," *\\)"), toupper(res_file_commentless))[[1]]], silent = T)
      # Problem if for instance we have  "$THETA 10 FIX " in first estimation
      if(length(grep("^\\$THETA", temp)) > 0){

        temp <- try(ejrzerjezlckjlefserlksejresres) # create an error
      }
      value <-  try(gsub(" ", "",str_split(temp, "=")[[1]][1]), silent = T)

      if(length(grep("MU_", value) >0)){

        temp <-  res_file_commentless[grep(paste0(".+ = .+?",value), toupper(res_file_commentless))]
        value <-  try(gsub(" ", "",str_split(temp, "=")[[1]][1]), silent = T)
      }

      if(class(temp) != "try-error"){
        theta_name_1st_mtd<- c(theta_name_1st_mtd, value)
      }
      a <- a + 1
    }

    for (a in theta_name_1st_mtd){
      if(length(theta_name_1st_mtd[theta_name_1st_mtd == a]) > 1){
        theta_name_1st_mtd[theta_name_1st_mtd == a] <- ""
      }
    }
    ### S E C O N D   M E T H O D:  C O M M E N T S   N A M E

    theta_begin <- grep("^[[:blank:]]?\\$THETA", toupper(res_file))[[1]]
    theta_end <- grep("^[[:blank:]]?\\$OMEGA", toupper(res_file))[[1]]
    res_file_shorten <- res_file[theta_begin:(theta_end-1)]

    res_file_shorten <- gsub(" *\\$THETA *", "", res_file_shorten)
    res_file_shorten <- res_file_shorten[-grep("^[[:space:]]*;", res_file_shorten)]
    res_file_shorten <- res_file_shorten[res_file_shorten != ""]
    res_file_shorten <- gsub(".+;", "", res_file_shorten)
    res_file_shorten <- gsub("^[[:space:]]*", "", res_file_shorten)

    # then it is either the first word either if 1 = Cl the afther =
    res_file_shorten_test <- gsub("[[:space:]].*","",res_file_shorten)

    if(sum(!is.na(as.double(res_file_shorten_test))) == length(res_file_shorten)){

      res_file_shorten <- gsub(" *[[:digit:]] *=? *", "", res_file_shorten)

    }else{

      res_file_shorten <- res_file_shorten_test
    }



    for (a in res_file_shorten){
      if(length(res_file_shorten[res_file_shorten == a]) > 1){
        res_file_shorten[res_file_shorten == a] <- ""
      }
    }


    ### C O M P A R I S O N    A N D    F I N A L     I N P U T

    theta_name <- character()
    for (a in 1:length(theta_name_1st_mtd)){

      test1 <- toupper(theta_name_1st_mtd[a])
      test2 <- toupper(res_file_shorten[a])

      if(is.na(test1) == T) test1 <- ""
      if(is.na(test2) == T) test2  <- ""

      # slightly different from model analysis
      if(test1 != ""){
        theta_name <- c(theta_name, test1)
      }else if(test2 != ""){

        theta_name <- c(theta_name, test2)
      }else{
        theta_name <- c(theta_name, paste0("THETA_",a))
      }


    }
    return(theta_name)

  }
  theta_names <- theta_name()



  # Get values --------------------------------------------------------------

  # Get initial values (work for cfl and res files)
  # carefull, I removed theta used as error at the end of model section

  thetabalise <- grep("\\$THETA", lines)
  omegabalise <- grep("\\$OMEGA", lines)


  lines[(thetabalise):(omegabalise - 1)] %>%
    gsub(pattern = "(;.+)| |\\(|\\)|(\\$THETA)", replacement = "") -> temp

  str_split( temp[temp != ""], ",", simplify = T) %>%
    as_tibble %>%
    rowid_to_column() %>%
    group_split(rowid) %>%
    map_chr(function(x){

      if("V2" %in% names(x)){

        if(x$V2 == "") return(gsub("(FIX)| ", "", x$V1))

        return(gsub("(FIX)| ", "", x$V2))
      }
      return(gsub("(FIX)| ", "", x$V1))
    }) -> inivalu


  restemp <- tibble( Param = tolower(theta_names),
                     Initial = inivalu )

  # Get Populational Values (only .res, hence the try if only .cfl is provided)

  try({
    thetabaliseres <- grep("\\THETA - VECTOR OF FIXED EFFECTS PARAMETERS", lines)[1]
    omegabaliseres <- grep("\\OMEGA - COV MATRIX FOR RANDOM EFFECTS ", lines)[1]

    lines[(thetabaliseres + 1):(omegabaliseres - 1)] %>%
      gsub(pattern = "TH.{,3}", replacement = "") -> temp

    str_split(temp[gsub(" ", "", temp) != ""] %>%
                paste(collapse = "")," ")[[1]] -> temp

    restemp$Pop <- temp[temp != ""]
  }, silent = T)

  # Get individal values
  try({

    restemp

    for(a in 1:nrow(run@individualValues)){

      line <- run@individualValues[a, ]
      names(line) <- tolower(names(line))
      namenewcol <- paste0("id", line$id)
      restemp[namenewcol] <- restemp$Pop

      for(b in 2:ncol(line)){

        restemp[[namenewcol]][[which(restemp$Param == names(line)[b])]] <- line[[b]]

      }

    }

    # baliseTable <- grep("\\$TABLE", lines)
    # baliseEndcfl <-   grep("NM-TRAN MESSAGES", lines)
    # if(length(baliseEndcfl) == 0) baliseEndcfl <- length(lines)
    #
    # lines[baliseTable:baliseEndcfl] %>%
    #   paste(collapse = " ") %>%
    #   gsub(pattern = ".*FILE *= *", replacement = "") %>%
    #   gsub(pattern = "  *.+", replacement = "") -> tabfile
    #
    # # ridiculous but lets go with that for now
    # file.path(str_split(path, "/")[[1]][-length(str_split(path, "/")[[1]])] %>%
    #             paste0(collapse = "/"), tabfile) -> tabfile
    #
    # # (TODO) handle all the ONEHEADER/NOHEADER etc situations
    # tabfile <- read.table(tabfile, skip = 1, header = T)
    #
    # tabfile[theta_names]
    # # the user need to put the values of parameter otherwise it would be useless...
    # # restemp$Param[restemp$Param %in% names(tabfile)]
  }, silent = T) # end try

  output$res <- restemp

  # initial values ----------------------------------------------------------

  # some write everything in the same row
  # for instance $MODEL      COMP=(DEPOT) COMP=(CENTRAL)
  # hence this str_split/collapse process to fill every cases
  cmtnames <-
    str_split(lines[grep("COMP *=? *\\(", lines)] %>%
                paste(collapse = "") %>%
                gsub(pattern = "COMP=?", replacement = "\nCOMP="), "\n")[[1]] %>%
    gsub(pattern = "(COMP *= *\\()|(\\).*| )", replacement = "")


  cmtnames <- cmtnames[-(cmtnames %in% c("", "$MODEL"))] %>% tolower()

  tibble(Cmt = cmtnames) %>%
    rownames_to_column() %>%
    mutate(t0 = map_chr(rowname, function(x){

      grep(paste0("A_0\\(",x,"\\) *="), lines) -> temp

      if(length(temp) == 0) return("0")

      lines[temp] %>%
        gsub(pattern = "(.+=)| ", replacement = "")

    })) %>%
    select(Cmt, t0) -> output$initial_values


  # model -------------------------------------------------------------------

  # get the advan

  lines[grep("\\$SUBROUTINE", lines):(grep("\\$MODEL", lines) - 1)] %>%
    gsub(pattern = ";.+", replacement = "") %>%
    paste0(collapse = " ") %>%
    gsub(pattern = "\\$SUBROUTINE *A", replacement = "A") %>%
    gsub(pattern = " ", replacement = "") -> subroutine

  balisePK <- grep("\\$PK",lines)
  baliseerror <- grep("\\$ERROR",lines)

  linesModel <- lines[(balisePK+1):(baliseerror-1)]

  # if advan5, we have to recreate the explicit lines
  if(subroutine == "ADVAN5"){


    linesModel[grep("^ *k[[:digit:]]", tolower(linesModel))] %>%
      tolower %>%
      as_tibble() %>%
      mutate(formul = gsub(".+= *", "", value)) %>%
      mutate(value = gsub(" *=.+", "", value)) %>%
      mutate(from = gsub("k|([[:digit:]]$)", "", value)) %>%
      mutate(value = gsub("..", "", value)) %>%
      rename(to = value) %>%
      arrange(from)-> equ

    tibble(cmt = 1:length(unique(equ$from))) %>%
      mutate(line = map_chr(cmt, function(a){



        paste0(equ$formul[equ$from == a]) -> out

        out <- case_when(length(out) == 0 ~ "",
                         length(out) == 1 ~ paste0("-A(",a,") * ", out ),
                         T ~ paste0("-A(",a,") * (", paste0(out, collapse = " + " ),")")) %>%
          unique


        inn <- paste0(equ$formul[equ$to == a])

        equ %>%
          filter(to == a) %>%
          mutate(inn = paste0("A(",from,") * ", formul )) %>%
          pull(inn) -> inn

        inn <- case_when(length(inn) == 0 ~ "",
                         T ~ paste0(inn,collapse = " + " ))

        if(out != "") inn <- paste0(" + ", inn)

        paste0("DADT(",a,") = ", out, inn)

      }))-> equ2


    linesModel <- c(linesModel, equ2$line)

    linesModel <- linesModel[-grep("^ *k[[:digit:]]", tolower(linesModel))]

  }

  linesModel <- tolower(linesModel)

  linesModel <- linesModel %>%
    gsub(pattern = "\\*", replacement = " * ") %>%
    gsub(pattern = "\\+", replacement = " + ") %>%
    gsub(pattern = "\\/", replacement = " / ") %>%
    gsub(pattern = "\\(", replacement = "( ") %>%
    gsub(pattern = "\\)", replacement = " )") %>%
    gsub(pattern = "\\/", replacement = " / ") %>%
    gsub(pattern = "\\-", replacement = " - ") %>%
    gsub(pattern = "=", replacement = " = ") %>%
    gsub(pattern = " \\* *\\* ", replacement = " ** ") %>%
    gsub(pattern = " \\-  \\-", replacement = " -- ")%>%
    gsub(pattern = " = ", replacement = " <- ") %>%
    gsub(pattern = "\\.eq\\.", replacement = " == ") %>%
    gsub(pattern = "\\.gt\\.", replacement = " > ") %>%
    gsub(pattern = "\\.lt\\.", replacement = " < ") %>%
    gsub(pattern = "\\.ge\\.", replacement = " >= ") %>%
    gsub(pattern = "\\.le\\.", replacement = " < ") %>%
    gsub(pattern = " then *$", replacement = "{") %>%
    gsub(pattern = " *endif *$", replacement = "}")


 # handle mtime( 1 ) <- lagt

for(a in grep("mtime *\\(", linesModel )){

analyse <- str_split(linesModel[a], pattern  = "<-")[[1]] %>%
  gsub(pattern = "( $)|(^ )",replacement = "")


# linesModel[a]
linesModel[a] <- paste0("if(t <", analyse[[2]],"){\n", analyse[[2]], "_bool = 0\n}else{\n",analyse[[2]], "_bool = 1\n}\n" )

linesModel <- gsub(paste0("mpast\\( ", gsub("(mtime *\\( *)|( *\\))", "", analyse[[1]]) ," \\)"), paste0(analyse[[2]], "_bool"), linesModel)

}

  #analyse THETA
  tibble(theta_names = tolower(theta_names)) %>%
    mutate(line = map_chr(theta_names,function(x){

      grep(paste0(tolower(x), " *="), tolower(lines)) -> temp
      # what hapen if several lines? or if/else structure
      # could be very complex to fit everything....
      if(length(temp) > 0) return(gsub("(;.+)", "", lines[temp])) #|(.*= *)

      return("")

    })) %>%
    mutate(distrib = map_chr(line, function(line){

      if(length(grep("exp\\(eta\\(", tolower(line))) == 1) return("logN")

      return("NoVar")
    })) %>%
    mutate(neta = map2_chr(line, distrib, function(line,distrib){

      if(distrib == "NoVar") return(NA)
      gsub(".*eta\\(", "", tolower(line)) %>%
        gsub(  pattern = "\\).+", replacement = "")

    })) -> thetas

  # if there is extra random effect -> create a parameter with fixed value 1?

  # remove eta
  linesModel %>%
    gsub(pattern = " *\\* *exp\\( *eta\\( *[[:digit:]]* *\\) *\\)",replacement =  "") %>%
    gsub(pattern = " *\\+ *eta\\( *[[:digit:]]* *\\)",replacement =  "") -> linesModel


  # replace (THETA(1)) etc by the name of the parameter
  # remove stupid line  like Vd = Vd? might be usefull in if else structure...
  # a = "lambd0411"
  for(a in tolower(theta_names)){

    nt <- which(tolower(theta_names) == a)

    linesModel <- gsub(paste0("theta\\( *", nt," *\\)"), a,linesModel )

    # remove line like VD = VD but NOT if inside a if else structure
    ngrep <- grep(paste0(a, " *<-"),linesModel)
    # print(a)
    # print(ngrep)
    linesModel[ngrep] -> temp
    gsub(paste0(tolower(a), "|(<-)| |(;.+)"),"",tolower(temp)) -> temp

    # print(temp)
    if(length(ngrep[temp == ""]) > 0){
      linesModel <- linesModel[- ngrep[temp == ""]]
    }

  }

  # linesModel[1:10]
  # modify cmtnames
  # a = "A"
  cmtnames <- gsub(",.+", "", cmtnames)
  for(a in tolower(cmtnames)){

    nt <- which(tolower(cmtnames) == a)

    linesModel <- gsub(paste0("adt\\( *", nt," *\\)"), a,linesModel )

    linesModel <- gsub(paste0("a\\( *", nt," *\\)"), a,linesModel )

  }

  # Some other element
  linesModel <- gsub(paste0(";.+"), "",linesModel )
  linesModel <- linesModel[linesModel !=""]
  linesModel <- linesModel[linesModel !=" "]
  linesModel <- linesModel[linesModel !="$des"]

  if(length(grep("obs *<- *dv" , linesModel)) > 0) linesModel <- linesModel[-grep("obs *<- *dv" , linesModel)]

 try({ # we need to add the output
  errorLine <- lines[baliseerror:(thetabalise-1)] %>% tolower
  # get value of output throug ipred
  errorLine[grep("ipred *=", errorLine)] %>%
    gsub(pattern = "(ipred *= *)|(;.+)",replacement = "") -> firstoutput


  if(firstoutput == "f"){
    # is there scaling factor ?
    nsoutput <- character()
    svalues <- character()
    try({
    linesModel[grep("s[[:digit:]] *<-", linesModel)] %>%
      gsub(pattern = "( *<-.+)|s", replacement = "") %>%
      as.double -> nsoutput

    linesModel[grep("s[[:digit:]] *<-", linesModel)] %>%
      gsub(pattern = "(.+<-)|(;.+)| ", replacement = "") -> svalues

    }, silent = T)

    if(length(nsoutput) > 0 ){
      outputtemp <- paste0(paste0("out", 1:length(nsoutput)), "_plot <- ", cmtnames[nsoutput], " / ",  svalues)
      linesModel <- linesModel[-grep("s[[:digit:]] *<-", linesModel)]
      }else{
      outputtemp <- paste0("out_plot <- ",cmtnames[length(cmtnames)])
    }

    linesModel <- c(linesModel, outputtemp)

  }
 }, silent = T)
  output$model <- linesModel %>% paste(collapse = "\n")
  # remove theta used as error (to keep as error)
 testthetaerror <- try(str_split(errorLine[grep("theta\\(", errorLine )], "theta\\("))

 if(class(testthetaerror) != "try-error"){

  thetatoremove <- gsub(pattern = "\\).*", replacement = "", testthetaerror[[1]][-1])
 output$res %>%
   slice(- as.double(thetatoremove) ) -> output$res


 }

# values ------------------------------------------------------------------

# just a simple rearrangement of values in order to directly inject it into Shiny app
 try({output$res[1:2] %>%
   left_join(thetas %>%
               rename(Param = theta_names)  %>%
               select(Param, distrib)) %>%
   table_param() -> output$values #to improve
 })

# matrix  ------------------------------------------------------------------
 rm(temp)
 #if .res
try({
 baliseomega <- grep("OMEGA - COV MATRIX FOR RANDOM EFFECTS", lines)[1]
 balisesigma <- grep("SIGMA - COV MATRIX FOR RANDOM EFFEC",lines)[1]

 str_split(lines[(baliseomega+1):(balisesigma-1)] %>%
   paste(collapse = ""), " ")[[1]] -> temp

 temp <- temp[temp != ""]
 temp <- temp[-grep("^ETA", temp)] %>% as.double
}, silent = T)
 #if .cfl
 if(!exists("temp")){



   temp <-  map(1:length(eta_names), function(x){

       temp <- 1:x * 0
       temp[length(temp)] <- 0.3
       return(temp)
     }) %>% reduce(c)
 }

matrixx <- matrix(NA, nrow = length(eta_names),length(eta_names),dimnames = list(eta_names,eta_names))
print(temp)

# start <- 1
for(a in 1:length(eta_names)){


  matrixx[a, 1:a] <- temp[1:a]
  temp <- temp[- (1:a)]

}
matrixx <- matrixx %>% as.data.frame()

namesma <- rownames(matrixx) %>% gsub(pattern = "eta_", replacement = "") %>% tolower()

rownames(matrixx) <- namesma
colnames(matrixx) <- namesma
# matrixx
# print(matrixx)
output$matrix <- matrixx %>%
  map_dfr(~ as.character(.x))



# Todisplay ---------------------------------------------------------------

output$Todisplay <- table_display(output$initial_values$Cmt)

# Input -------------------------------------------------------------------

output$input <- table_input(var = factor(output$initial_values$Cmt[[1]], levels = output$initial_values$Cmt ))

  return(output)
}
