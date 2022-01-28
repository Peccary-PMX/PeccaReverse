#' @description azeazezaezae
#'
#' @param model
#' @export

# library(diagram)
#
# library(peccary)
# model <- "##Pk.1comp.Cl.Vd.difEq
# ke <- Cl/Vd
# dCentral <- -ke * Central
# Conc  <- Central/ Vd
#
# ##PD.Simeoni
# tumVol_output <- X1 + X2 + X3 + X4
# growth <- lambda0 * tumVol /( ( (1 + lambda0 * X1 /lambda1)**psi)**(1/psi))
# dX1 <- X1 * (growth - Conc * k2)
# dX2 <- X1 * Conc * k2 - k1 * X2
# dX3 <- k1 * (X2 - X3)
# dX4 <- k1 * (X3 - X4)
# "
#
#
# model <- "##Pk.2comp.ke.k12.k21.difEq
# #k12 <- Q/V1
# #k21 <- Q/V2
# #ke <- Cl/V1
# dCentral <- k21 * Periph - k12 * Central - ke * Central
# dPeriph <- - k21 * Periph + k12 * Central
# Conc <- Central / V1
#
# ##PD.Simeoni
# tumVol_output <- X1 + X2 + X3 + X4
# growth <- lambda0 * tumVol /( ( (1 + lambda0 * X1 /lambda1)**psi)**(1/psi))
# dX1 <- X1 * (growth - Central * k2)
# dX2 <- X1 * Central * k2 - k1 * X2
# dX3 <- k1 * (X2 - X3)
# dX4 <- k1 * (X3 - X4)
# "
#
# model <- "##Pk.2comp.ke.k12.k21.difEq
# #k12 <- Q/V1
# #k21 <- Q/V2
# #ke <- Cl/V1
# dCentral <- k21 * Periph - k12 * Central - ke * Central
# dPeriph <- - k21 * Periph + k12 * Central
# Conc <- Central / V1
#
# dIL7 <- kinIL7 -koutIL7 * IL7
# "
#
# # developement
# x <- "X1 * (growth - Conc * k2) + X2 * (aha - lol)"
# x <- " - X1 * ke"
# y <- "X1"
#
# x <- "k1 * (X2 - X3)"
# y <- "X2"
# develo(x)
# par(mar=c(5.1, 4.1, 4.1, 2.1))
# pecc_scheme(model, relsize = 1)
pecc_scheme <- function(model, relsize = 1){

develo <- function(x , y = NULL){


  # recursif system
  if(is.null(y)){

    totest <- str_split(gsub("[[:punct:]]"," ", x) %>%
      gsub(pattern = "  ", replacement = " "), pattern = " ")[[1]]

    totest <- totest[totest != ""]

    finaloutput <- x

    for(a in totest){

    finaloutput <- develo(finaloutput, a)

    }
    return(finaloutput)
  }


  # ou se trouve X1
  whereiscomp <- str_locate(x, y)
  # start[[1,1]]

  # Remove all before (for now you need to do x * (a + b) and not (a + b) * x...)
  xtemp <-
    substr(x,whereiscomp[[1,1]] - 1,  nchar(x) )

  # remove end of factorisation
  str_locate_all(xtemp, "\\(|\\)")[[1]] %>%
    as_tibble() %>%
    mutate(which = map_chr(start, ~ substr(xtemp, .x, .x))) %>%
    mutate(test = 0) -> parenAna
# print(parenAna)
# print(nrow(parenAna))
  if(nrow(parenAna %>% filter(!is.na(end))) == 0) return(x)

  parenAna[["test"]][[1]] <- 1

  for(a in 2:nrow(parenAna)){

    if(parenAna[["which"]][[a]] == "(")  parenAna[["test"]][[a]] <- parenAna[["test"]][[a - 1]] + 1
    if(parenAna[["which"]][[a]] == ")")  parenAna[["test"]][[a]] <- parenAna[["test"]][[a - 1]] - 1

  }

  parenAna %>%
    mutate(whereiscomp =  str_locate(xtemp, y)[[2]]) %>%
    filter(whereiscomp < end) %>%
    filter(test == 0) %>%
    slice(1) %>%
    pull(end) -> endbloc

  parenAna %>%
    filter(start < endbloc) %>%
    pull(1) -> starparen

  xtemp2 <- substr(xtemp,0, endbloc + 1 )

  tocheckfactorisation <- substr(xtemp,0 , starparen -1) %>% gsub(pattern = " ",replacement =  "")

  if(substr(tocheckfactorisation, nchar(tocheckfactorisation), nchar(tocheckfactorisation)) == "*"){

    xtemp2inspar <- substr(xtemp,starparen +1, endbloc-1) # inspar = inside parenthesis
    allblocs <- str_locate_all(xtemp2inspar, "\\+|-")[[1]]

    finalline <- ""

    for(a in 1:(nrow(allblocs) + 1)){

      if(a == 1){

        finalline <- paste0(y, " * ", substr(xtemp2inspar,0, allblocs[[a,1]])  )

      }else if(a != nrow(allblocs) + 1){

        finalline <- paste0(finalline, y, " * ", substr(xtemp2inspar,allblocs[[a -1,1]] + 1, allblocs[[a,1]])  )
      }else{

        finalline <- paste0(finalline," ", y, " * ", substr(xtemp2inspar,allblocs[[a -1,1]] + 1,nchar(xtemp2inspar))  )

      }

    } # end for loop

    toreplace <- xtemp2 %>%
      gsub(pattern = "\\*", replacement = "\\\\*") %>%
      gsub(pattern = "\\(", replacement = "\\\\(") %>%
      gsub(pattern = "\\)", replacement = "\\\\)")


    return(    gsub(toreplace, finalline, x))

  }

return(x)
}







deSModel <- deSolve_pecc(model)

box <- deSModel$state
# model <- deSModel$model

model <- str_split(model, "\n")[[1]]
# target = "IL7"
# from = "from"



crossing(target = box, from = box) %>%
    mutate(equation = map2(target, from, function(target, from){

      # print("new")
      # print(target)
      # print(from)
    test <- try({
      line <-gsub(".+<-", "",  model[grepl(paste0("d",target, " *<-"),model)])
      line <- develo(line)

      line <- gsub("\\+ *", "+ ", line)
      line <- gsub("- *", "- ", line)

      splitline <- str_split(line, pattern = "\\+|\\-")[[1]]

      # finalines <-
      tibble(equation = splitline[grep(from, splitline)]) %>%
        mutate(equation = map_chr(equation, function(x){

        start <-   str_locate(string = line, pattern = fixed(x))[1]

        if(start > 1){
        # readd "-" if needed
          substr(line,0, start) %>%
            gsub(pattern = " ", replacement = "") %>%
            grepl(pattern = "-$") -> testminus

        if(testminus) x <- paste0("- ", x)

        }
        x
        }))



    })

    if(class(test) == "try-error") return(tibble(equation= ""))
    if(length(test) == 0) return(tibble(equation= ""))
    return(test)

    })) %>%
    unnest() %>%
    filter(equation != "") %>%
    mutate(equation = gsub("(^ *)|( *$)", "", equation)) %>%
    mutate(equation = gsub("  *", " ", equation)) %>%
  mutate(order = map_dbl(equation, function(x){

    sum(map_lgl(box, ~ grepl(paste0(" ", .x, " "), paste0(" ",x," "))))

  })) -> allblocs ; allblocs




# Handling second order elimination

allblocs %>%
  filter(order == 2, target == from) %>%
  mutate(secorder = map2_chr(equation, target, function(x,y){

    temp <- box[map_lgl(box, ~ grepl(.x, x) )]
    temp[temp != y]
  })) -> secorder

if(nrow(secorder) >0){
for(a in 1:nrow(secorder)){


  secorder %>%
    slice(a) -> secordertemp

  # rm the two first lines
  allblocs <- allblocs %>%
    filter(!(order == 2 & from == secordertemp$secorder))

  # add
}
}


  #remove elimination that are input

for(a in box){

  allblocs %>%
    filter(target == a & from == a) %>%
    pull(equation) -> bloc

  for(b in bloc[grepl("^-", bloc)]){

gsub("^- *", "", b) -> totest


  allblocs %>%
    mutate(test1 = if_else(equation == totest, 1, 0)) %>%
    mutate(test2 = if_else(from == a, 1, 0)) %>%
    mutate(testfinal = test1 + test2) %>%
    filter(testfinal == 2) %>%
    nrow -> finaltest

  if(finaltest > 0) allblocs <- allblocs %>%
    filter(!(target == a & from == a &  equation == b))


  }

}


## now creation !!
library(diagram)
par(ask=F)
elpos <- coordinates (N = length(box), mx = 0, my = 0, relsize = relsize)

# elpos2 <-coordinates (N = length(box), mx = 0, my = 0, relsize = 0.8 ); elpos2

# elpos <- coordinates (pos = c(2, 4), mx = 0, my = 0.1)
# a <- 1
# a <- 2
# a <- 3
# a <- 4
# a <- 5

text <- T
a <- 5
a <- a + 1


allblocs$done <- 0
openplotmat(main = "test")
for(a in 1:nrow(allblocs)){

  allblocs %>%
    slice(a) -> line;line

  nbox <- which(box == line$target)

  texttemp <- gsub("- *","",line$equation) %>%
    gsub(pattern = paste0(" *\\* *",line$from ), replacement = "") %>%
    gsub(pattern = paste0(line$from , " *\\* *"), replacement = "")

  # direct elimimnation
  if(line$target == line$from & grepl("^ *-", line$equation)){




    straightarrow (from = elpos[nbox, ] * c(0.9,1), to = elpos[nbox, ] * c(0.9,0.7), path = "right", arr.pos= 1)


  if(text == T)  text(elpos[nbox, 1] * 0.8, elpos[nbox, 2] * 0.8, texttemp)


  }else if(line$target == line$from & !grepl("^ *-", line$equation)){

    # self growth

    selfarrow(pos =  elpos[nbox, ], arr.pos = 0.25, path = "U")
    if(text == T) text(elpos[nbox, 1], elpos[nbox, 2] +0.15, texttemp)

  }else{

    # transfert

    # does a reverse transfert also exist?
    allblocs %>%
      filter(target == line$from, from == line$target) -> testreverse


    nboxtarget <- which(box == line$target)
    nboxfrom <- which(box == line$from)
    fromtemp <- elpos[nboxfrom, ]
    totemp <- elpos[nboxtarget, ]
    xtemp <-   mean(c(elpos[nboxfrom, 1], elpos[nboxtarget, 1]))
    ytemp <- mean(c(elpos[nboxfrom, 2], elpos[nboxtarget, 2]))


    if(round(totemp[[1]],3) != round(fromtemp[[1]],3)){
      ytemp <- ytemp * 1.2
    }else{

      xtemp <- xtemp * 0.9
    }
    ##
    if(nrow(testreverse) > 0){

      # if not done yet -> slide to the left or bottom

      if(testreverse$done == 0){

        fromtemp[1] <-  fromtemp[1] * 0.95
        totemp[1] <-  totemp[1] * 0.95
        xtemp <- xtemp

      }else{
      # otherwise -> slide to the right or upper
        fromtemp[1] <-  fromtemp[1] * 1.02
        totemp[1] <-  totemp[1] * 1.02
        xtemp <- xtemp /0.95 * 1.14
      }



    }



      straightarrow (from = fromtemp, to = totemp)

      if(text == T)  text(xtemp, ytemp ,texttemp)



  }

  # handle 2nd order

  if(line$order == 2){



    fromtemp <- box[map_lgl(box, ~ grepl(.x, line$equation))]
    fromtemp <- fromtemp[fromtemp != line$from]
    nboxfrom <- which(box == fromtemp)
    fromtemp <- elpos[nboxfrom, ]

    nboxto1 <- which(box == line$from)
    nboxto2 <- which(box == line$target)
    xtemp <-   mean(c(elpos[nboxto1, 1], elpos[nboxto2, 1]))
    ytemp <- mean(c(elpos[nboxto1, 2], elpos[nboxto2, 2]))

    totemp <- c(xtemp, ytemp)

    straightarrow (from = fromtemp, to = totemp, lty =1, lcol = "red", arr.pos = 0.8,
                   segment = c(0,0.8), arr.type = "ellipse")

  }


 # update done

  allblocs[["done"]][[a]] <- 1
}

for(a in box){


  textrect(elpos[which(box == a), ], 0.1, 0.1, lab = a, shadow.size = 0)
}


}



