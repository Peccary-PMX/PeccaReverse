# func <- "
#
# Cb <- 3.0
# Cv <- 0.0005
#
# Killv <- KKill_v*Cv
# killb <- (Kmax_b*Cb**H)/((KC50_b*Psi)**H+Cb**H)
#
#
# dVor1 <- (Killv-Vor1)/Ktr_v
#
#
# dBTZ1 <- 1.0/Ktr_b*(Killb-BTZ1)
# dBTZ2 <- 1.0/Ktr_b*(BTZ1-BTZ2)
# dBTZ3 <- 1.0/Ktr_b*(BTZ2-BTZ3)
# dBTZ4 <- 1.0/Ktr_b*(BTZ3-BTZ4)
#
#
# dN <- Kg*N-(BTZ4+Vor1)*N
#
# "
#' desolve_to_adapt
#' @export
desolve_to_adapt <- function(func){

  bloclist <- list()

  analysis  <-  deSolve_pecc(func)

# Bloc Memo ---------------------------------------------------------------


  bloclist$memo <-
"**********************************************************************
C                           ADAPT                                     *
C                         Version 5                                   *
C**********************************************************************
C                                                                     *
C                           MODEL                                     *
C                                                                     *
C    This file contains Fortran subroutines into which the user       *
C    must enter the relevant model equations and constants.           *
C    Consult the User's Guide for details concerning the format for   *
C    entered equations and definition of symbols.                     *
C                                                                     *
C       1. Symbol-  Parameter symbols and model constants             *
C       2. DiffEq-  System differential equations                     *
C       3. Output-  System output equations                           *
C       4. Varmod-  Error variance model equations                    *
C       5. Covmod-  Covariate model equations (ITS,MLEM)              *
C       6. Popinit- Population parameter initial values (ITS,MLEM)    *
C       7. Prior -  Parameter mean and covariance values (ID,NPD,STS) *
C       8. Sparam-  Secondary parameters                              *
C       9. Amat  -  System state matrix                               *
C                                                                     *
C**********************************************************************

C######################################################################C
"


# bloc Symbol -------------------------------------------------------------

  bloclist$symbol <-  paste0("C######################################################################C

        Subroutine SYMBOL
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'
",
   "
C----------------------------------------------------------------------C
C   Enter as Indicated                                                 C
C----c-----------------------------------------------------------------C


      NDEqs   =  ", length(analysis$state),"   ! Enter # of Diff. Eqs.
      NSParam =  ", length(analysis$parameter) ,"   ! Enter # of System Parameters.
      NVparam =  2   ! Enter # of Variance Parameters.
      NSecPar =  0   ! Enter # of Secondary Parameters.
      NSecOut =  0  ! Enter # of Secondary Outputs (not used).
      Ieqsol  =  1  ! Model type: 1 - DIFFEQ, 2 - AMAT, 3 - OUTPUT only.
      Descr   = ' Name '

C----------------------------------------------------------------------C
C   Enter Symbol for Each System Parameter (eg. Psym(1)='Kel')         C
C----c-----------------------------------------------------------------C
\n",
paste0(imap_chr(analysis$parameter, ~ paste0("Psym(", .y, ") = ", .x)), collapse =  "\n"),
"\n
C----------------------------------------------------------------------C
C   Enter Symbol for Each Variance Parameter {eg: PVsym(1)='Sigma'}    C
C----c-----------------------------------------------------------------C
            PVsym(1)='Slope'
            PVsym(2)='intercept'

C----------------------------------------------------------------------C
C   Enter Symbol for Each Secondary Parameter {eg: PSsym(1)='CLt'}     C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C

Return
End
"
    )


# Bloc diffeq -------------------------------------------------------------

bloclist$diffeq <- paste0("C######################################################################C

                          Subroutine DIFFEQ(T,X,XP)
                          Implicit None

                          Include 'globals.inc'
                          Include 'model.inc'

                          Real*8 T,X(MaxNDE),XP(MaxNDE)
                          Real*8 ",paste0(c(analysis$parameter, analysis$state), collapse = ", ") , "
C----------------------------------------------------------------------C
C   Enter Differential Equations Below  {e.g.  XP(1) = -P(1)*X(1) }    C
C----c-----------------------------------------------------------------C")


  model <- deparse(analysis$model)
  model <- model[3:(length(model)-3)]

  for(a in length(analysis$parameter):1){

model <- c(paste0(analysis$parameter[[a]], " = P(", a,")"), model)

  }

  for(a in  length(analysis$state):1){

    model <- c(paste0(analysis$state[[a]], " = X(", a,")"), model)
    model <- gsub(paste0("d", analysis$state[[a]], " *<-"), paste0("XP(", a,") <-"), model)
  }

  model <- gsub(" *<- *", " = ", model)
  model <- gsub("^ *", "", model)

  bloclist$diffeq2 <- paste0(paste0(model, collapse = "\n"), "
                          Return
        End")


# bloc output -------------------------------------------------------------

  bloclist$output <- paste0(


"
        Subroutine OUTPUT(Y,T,X)
        Implicit None

        Include 'globals.inc'
        Include 'model.inc'

        Real*8 Y(MaxNOE),T,X(MaxNDE)
        Real*8 ???


CC
C----------------------------------------------------------------------C
C   Enter Output Equations Below   {e.g.  Y(1) = X(1)/P(2) }           C
C----c-----------------------------------------------------------------C

????

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
        Return
        End",


"C######################################################################C

Subroutine VARMOD(V,T,X,Y)
Implicit None

Include 'globals.inc'
Include 'model.inc'

Real*8 V(MaxNOE),T,X(MaxNDE),Y(MaxNOE)

CC
C----------------------------------------------------------------------C
C   Enter Variance Model Equations Below                               C
C         {e.g. V(1) = (PV(1) + PV(2)*Y(1))**2 }                       C
C----c-----------------------------------------------------------------C
V(1)=(PV(2)+PV(1)*Y(1))**2.0

C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
Return
End

C######################################################################C

Subroutine COVMOD(Pmean, ICmean, PC)
C  Defines any covariate model equations (MLEM, ITS)
Implicit None

Include 'globals.inc'
Include 'model.inc'

Real*8 PC(MaxNCP)
Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)

CC
C----------------------------------------------------------------------C
C     Enter # of Covariate Parameters                                  C
C----c-----------------------------------------------------------------C

NCparam = 0    ! Enter # of Covariate Parameters.

CC
C----------------------------------------------------------------------C
C   Enter Symbol for Covariate Params {eg: PCsym(1)='CLRenal'}         C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   For the Model Params. that Depend on Covariates Enter the Equation C
C         {e.g. Pmean(1) =  PC(1)*R(2) }                               C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
Return
End

C######################################################################C

Subroutine POPINIT(PmeanI,ICmeanI,PcovI,ICcovI, PCI)
C  Initial parameter values for population program parameters (ITS, MLEM)

Implicit None

Include 'globals.inc'
Include 'model.inc'

Integer I,J
Real*8 PmeanI(MaxNSP+MaxNDE), ICmeanI(MaxNDE)
Real*8 PcovI(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcovI(MaxNDE,MaxNDE)
Real*8 PCI(MaxNCP)

CC
C----------------------------------------------------------------------C
C  Enter Initial Values for Population Means                           C
C          {  e.g. PmeanI(1) = 10.0    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Initial Values for Pop. Covariance Matrix (Lower Triang.)    C
C         {  e.g. PcovI(2,1) = 0.25    }                               C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Values for Covariate Model Parameters                        C
C         {  e.g. PCI(1) = 2.0    }                                    C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
Return
End

C######################################################################C

Subroutine PRIOR(Pmean,Pcov,ICmean,ICcov)
C  Parameter mean and covariance values for MAP estimation (ID,NPD,STS)
Implicit None

Include 'globals.inc'
Include 'model.inc'

Integer I,J
Real*8 Pmean(MaxNSP+MaxNDE), ICmean(MaxNDE)
Real*8 Pcov(MaxNSP+MaxNDE,MaxNSP+MaxNDE), ICcov(MaxNDE,MaxNDE)

CC
C----------------------------------------------------------------------C
C  Enter Nonzero Elements of Prior Mean Vector                         C
C          {  e.g. Pmean(1) = 10.0    }                                C
C----c-----------------------------------------------------------------C


CC
C----------------------------------------------------------------------C
C   Enter Nonzero Elements of Covariance Matrix (Lower Triang.)       C
C         {  e.g. Pcov(2,1) = 0.25    }                                C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
Return
End

C######################################################################C

Subroutine SPARAM(PS,P,IC)
Implicit None

Include 'globals.inc'

Real*8 PS(MaxNSECP), P(MaxNSP+MaxNDE), IC(MaxNDE)

CC
C----------------------------------------------------------------------C
C   Enter Equations Defining Secondary Paramters                       C
C           {  e.g.  PS(1) = P(1)*P(2)   }                             C
C----c-----------------------------------------------------------------C




C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
Return
End

C######################################################################C

Subroutine AMAT(A)
Implicit None

Include 'globals.inc'
Include 'model.inc'

Integer I,J
Real*8 A(MaxNDE,MaxNDE)

DO I=1,Ndeqs
Do J=1,Ndeqs
A(I,J)=0.0D0
End Do
End Do

CC
C----------------------------------------------------------------------C
C   Enter non zero elements of state matrix  {e.g.  A(1,1) = -P(1) }   C
C----c-----------------------------------------------------------------C


C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C
Return
End

C######################################################################C"



  )


  return( paste0(bloclist, collapse = "\n") %>%
            gsub(pattern = "\n *",replacement =  "\n"))
}
