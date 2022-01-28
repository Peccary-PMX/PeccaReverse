
# text <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_21_01_11_5ytype/cov_analsysis/Ref_without_cov_no_growth_estimElim_IL7onExpanHillIL7_10HigherEff4_IL750free2.mlxtran"
# text <- "file:///D:/Peccary/Exemple_demo/ivivc.res"
# text <- "file:///D:/these/Second_project/model_lidner.csv"

# read.table("file:///D:/these/Second_project/model_lidner.csv", header = T, sep = ";",stringsAsFactors = F)
# text <- "file:///D:/Peccary/Exemple_demo/run_nonmem/1_comp.res"
# pecc_import_model(text)

#' pecc_import_model
#' @export
pecc_import_model <- function(text){


  config <- file.path(find.package("peccary"),  "Peccary_configuration.R")
  config <- readLines(config)  # keyword_NONMEM from configuration file
  eval(parse_expr(config[grep("ext_nm_res_file", config)]))
  pre_keyword_NONMEM <- paste0("\\", ext_nm_res_file, "$")
  # pre_keyword_NONMEM <- c("\\.res$","\\.lst$")
  keyword_NONMEM <- paste0("(", pre_keyword_NONMEM, ")", collapse = "|")


 firstline <-  gsub("(\n|#).+", "", text) %>%
    gsub(pattern = "(^ *)|( *$)", replacement =  "")

# if csv -> QSP

 if(grepl("csv$", text)){


first_line <- readLines(text, n = 1)

sep = ";"
test <- str_split(first_line, pattern = ";")[[1]]
if(length(test) < 3){
  test <- str_split(first_line, pattern = ";")[[1]]
    sep = ","
}


tabl <- read.table(text, sep = sep, header = T, stringsAsFactors = F)

a <- list()

a$model <- QSP_to_Peccary(tabl, outputExpr = F)

return(a)
 }

## if mlxtran
 if(grepl("\\.mlxtran", firstline)){

  return(monolix_to_desolve(text))

## if an other file -> just read it first
 }else if(grepl(keyword_NONMEM, firstline)){

   return(nonmem_to_desolve(firstline))

 }else if(file.exists(firstline)){

   text <- readLines(firstline) %>% paste0(collapse = "\n")


   # if it is a Monolix code


 }

 if(grepl("ddt_",text) %>% max ){

   return(monolix_to_desolve_str_model(text))

 }else if(grepl("\n* *XP\\([[:digit:]]\\) *=", text) %>% max){

   return(adapt_to_desolve(text))

 }

 return(NULL)

}
