# model <- "C----------------------------------------------------------------------C
# C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
# C----c-----------------------------------------------------------------C
# kg=P(1)
# Kmax_b=P(2)
# KC50_b = P(3)
# Psi = P(4)         !Interaction parameter
# KKill_v=P(5)
# H = P(6)
# Ktr_b=P(7)
# Ktr_v=P(8)
#
#
# Cb=3.0  !Conc of BTZ
# Cv= 0.0005  ! Conc of Vor 0.5nM
#
# Killv= KKill_v*Cv
# killb= (Kmax_b*Cb**H)/((KC50_b*Psi)**H+Cb**H)
#
# !Vor 1 compt Transit
# Vor1=X(1)     ! Transit compartment 1
# XP(1)=(Killv-Vor1)/Ktr_v
#
# BTZ1=X(2)     !BTZ 2nM Transit compartment 1
# BTZ2=X(3)     !BTZ 2nM transit compartment 2
# BTZ3=X(4)     !BTZ 2nM transit compartment 3
# BTZ4=X(5)     !BTZ 2nM tansit compartment 4
#
# !BTZ 4 compt Transit
# XP(2)=1.0/Ktr_b*(Killb-BTZ1)       !BTZ 2nM trasnit 1
# XP(3)=1.0/Ktr_b*(BTZ1-BTZ2)      !BTZ 2nM trasnit 2
# XP(4)=1.0/Ktr_b*(BTZ2-BTZ3)      !BTZ 2nM trasnit 3
# XP(5)=1.0/Ktr_b*(BTZ3-BTZ4)      !BTZ 2nM trasnit 4
#
# !Double combination Psi in KC50 of BTZ
# N=X(6)    !Cell Number
# XP(6)= Kg*N-(BTZ4+Vor1)*N
# "


#' adapt_to_desolve
#' @export
adapt_to_desolve <- function(model){

  output <- list()
  model <- str_split(model, "\n")[[1]]
  model <- gsub("^ *", "", model)
  if(grepl("^C--", model[[3]])) model <- model[4:length(model)]

### remove model declaration
  model <- model[- grep("= *P\\([[:digit:]]*\\)", model)]
  model <- gsub(" *\\!.+", "", model)

### Find names of parameter
  tibble(cmt = model[grep("= *X\\([[:digit:]]*\\)", model)]) %>%
    mutate(number = gsub("(.+= *X\\()|(\\) *)", "", cmt)) %>%
    mutate(cmt = gsub(" *=.+", "", cmt)) -> cmt

  for(a in 1:nrow(cmt)){

    model <- gsub(paste0("XP\\(", cmt$number[[a]], "\\)"), paste0("d", cmt$cmt[[a]]), model)
    model <- model[-grep(paste0("= *X\\(", cmt$number[[a]]), model)]
  }


  model <- gsub(" *= *", " <- ", model)
  model <- gsub(" <- <- ", " == ", model)

  output$model <- paste0(model, collapse = "\n")



  return(output)
}
