        !COMPILER-GENERATED INTERFACE MODULE: Tue Nov 28 23:35:29 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITDRV__genmod
          INTERFACE 
            SUBROUTINE NITDRV(N,XCUR,FCUR,XPLS,FPLS,STEP,F,JACV,RPAR,   &
     &IPAR,FTOL,STPTOL,NNIMAX,IJACV,IKRYSL,KDMAX,IRPRE,IKSMAX,IRESUP,   &
     &IFDORD,IBTMAX,IETA,ITERM,NFE,NJVE,NRPRE,NLI,NNI,NBT,RWORK,DINPR,  &
     &DNORM)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XCUR(N)
              REAL(KIND=8) :: FCUR(N)
              REAL(KIND=8) :: XPLS(N)
              REAL(KIND=8) :: FPLS(N)
              REAL(KIND=8) :: STEP(N)
              EXTERNAL F
              EXTERNAL JACV
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              REAL(KIND=8) :: FTOL
              REAL(KIND=8) :: STPTOL
              INTEGER(KIND=4) :: NNIMAX
              INTEGER(KIND=4) :: IJACV
              INTEGER(KIND=4) :: IKRYSL
              INTEGER(KIND=4) :: KDMAX
              INTEGER(KIND=4) :: IRPRE
              INTEGER(KIND=4) :: IKSMAX
              INTEGER(KIND=4) :: IRESUP
              INTEGER(KIND=4) :: IFDORD
              INTEGER(KIND=4) :: IBTMAX
              INTEGER(KIND=4) :: IETA
              INTEGER(KIND=4) :: ITERM
              INTEGER(KIND=4) :: NFE
              INTEGER(KIND=4) :: NJVE
              INTEGER(KIND=4) :: NRPRE
              INTEGER(KIND=4) :: NLI
              INTEGER(KIND=4) :: NNI
              INTEGER(KIND=4) :: NBT
              REAL(KIND=8) :: RWORK(*)
              REAL(KIND=8) :: DINPR
              EXTERNAL DINPR
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
            END SUBROUTINE NITDRV
          END INTERFACE 
        END MODULE NITDRV__genmod
