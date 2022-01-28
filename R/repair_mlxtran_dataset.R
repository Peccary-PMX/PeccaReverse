# path_mxtran <- "file:///D:/these/Pecc_test/3_Models/1_Models/lymphodepletion2/endogeneousProlif.mlxtran"
# path <- "D:/these/Pecc_test/3_Models/1_Models/lymphodepletion2"
# path <- "file:///D:/these/Pecc_test/3_Models/1_Models/000_20_04_24"
#' @export
pecc_repair_mlxtran_header <- function(path){

  path <-  path %>%
    gsub(pattern = "file:///", replacement = "")

if(length(grep(".mlxtran", path)) == 1){

  path_mxtran <- path

}else{

  path_mxtran <-  list.files(path)
  path_mxtran <- path_mxtran[grep(".mlxtran", path_mxtran)]
  path_mxtran <- paste(path,path_mxtran,sep = "/"  )
}
  individual_run <- function(path_mxtran){
  path <-  path_mxtran %>%
    gsub(pattern = "file:///", replacement = "")


  root <- gsub( gsub(".+/", "", path), "", path)
  lines <- readLines(path)

  ### Find dataset file from mxtran


  dataset <- lines[grep("file = ", lines)]; dataset <- dataset[1] %>%
    gsub(pattern = "(file =)|(')|( )",replacement =  "")

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

  #

  df <- read.table(path_dataset, header = T, sep = delimiter, na = ".", nrows = 1) %>%
    as_tibble

  balise_header <- grep("header = ", lines)



  fileConn<-file(gsub(".mlxtran", "_peccbckp.mlxtran",path_mxtran))
  writeLines(lines, fileConn)
  close(fileConn)


  lines[balise_header] <-  paste0("header = {", paste0(names(df), collapse = ", "),"}")



  fileConn<-file(path_mxtran)
  writeLines(lines, fileConn)
  close(fileConn)
return(paste0(path_mxtran, " - done"))
}

  map(path_mxtran, ~ individual_run(path_mxtran = .x) )

}
