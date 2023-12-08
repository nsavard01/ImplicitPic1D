<<<<<<< HEAD
        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec  7 16:42:37 2023
=======
        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec  7 17:55:47 2023
>>>>>>> temp
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITGM__genmod
          INTERFACE 
            SUBROUTINE NITGM(N,XCUR,FCUR,FCNRM,STEP,ETA,F,JACV,RPAR,IPAR&
     &,IJACV,IRPRE,IKSMAX,IRESUP,IFDORD,NFE,NJVE,NRPRE,NLI,KDMAX,KDMAXP1&
     &,VV,RR,SVBIG,SVSML,W,RWORK,RSNRM,DINPR,DNORM,ITRMKS)
              INTEGER(KIND=4) :: KDMAXP1
              INTEGER(KIND=4) :: KDMAX
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XCUR(N)
              REAL(KIND=8) :: FCUR(N)
              REAL(KIND=8) :: FCNRM
              REAL(KIND=8) :: STEP(N)
              REAL(KIND=8) :: ETA
              EXTERNAL F
              EXTERNAL JACV
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              INTEGER(KIND=4) :: IJACV
              INTEGER(KIND=4) :: IRPRE
              INTEGER(KIND=4) :: IKSMAX
              INTEGER(KIND=4) :: IRESUP
              INTEGER(KIND=4) :: IFDORD
              INTEGER(KIND=4) :: NFE
              INTEGER(KIND=4) :: NJVE
              INTEGER(KIND=4) :: NRPRE
              INTEGER(KIND=4) :: NLI
              REAL(KIND=8) :: VV(N,KDMAXP1)
              REAL(KIND=8) :: RR(KDMAX,KDMAX)
              REAL(KIND=8) :: SVBIG(KDMAX)
              REAL(KIND=8) :: SVSML(KDMAX)
              REAL(KIND=8) :: W(KDMAX)
              REAL(KIND=8) :: RWORK(N)
              REAL(KIND=8) :: RSNRM
              REAL(KIND=8) :: DINPR
              EXTERNAL DINPR
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
              INTEGER(KIND=4) :: ITRMKS
            END SUBROUTINE NITGM
          END INTERFACE 
        END MODULE NITGM__genmod
