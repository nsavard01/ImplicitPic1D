        !COMPILER-GENERATED INTERFACE MODULE: Tue Dec  5 23:58:19 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITSTB__genmod
          INTERFACE 
            SUBROUTINE NITSTB(N,XCUR,FCUR,FCNRM,STEP,ETA,F,JACV,RPAR,   &
     &IPAR,IJACV,IRPRE,IKSMAX,IFDORD,NFE,NJVE,NRPRE,NLI,R,RTIL,P,PHAT,V,&
     &T,RWORK1,RWORK2,RSNRM,DINPR,DNORM,ITRMKS)
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
              INTEGER(KIND=4) :: IFDORD
              INTEGER(KIND=4) :: NFE
              INTEGER(KIND=4) :: NJVE
              INTEGER(KIND=4) :: NRPRE
              INTEGER(KIND=4) :: NLI
              REAL(KIND=8) :: R(N)
              REAL(KIND=8) :: RTIL(N)
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: PHAT(N)
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: T(N)
              REAL(KIND=8) :: RWORK1(N)
              REAL(KIND=8) :: RWORK2(N)
              REAL(KIND=8) :: RSNRM
              REAL(KIND=8) :: DINPR
              EXTERNAL DINPR
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
              INTEGER(KIND=4) :: ITRMKS
            END SUBROUTINE NITSTB
          END INTERFACE 
        END MODULE NITSTB__genmod
