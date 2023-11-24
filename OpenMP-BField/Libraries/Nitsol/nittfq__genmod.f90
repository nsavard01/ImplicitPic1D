        !COMPILER-GENERATED INTERFACE MODULE: Fri Nov 24 20:04:09 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITTFQ__genmod
          INTERFACE 
            SUBROUTINE NITTFQ(N,XCUR,FCUR,FCNRM,STEP,ETA,F,JACV,RPAR,   &
     &IPAR,IJACV,IRPRE,IKSMAX,IFDORD,NFE,NJVE,NRPRE,NLI,R,RCGS,RTIL,D,P,&
     &Q,U,V,Y,RWORK1,RWORK2,RSNRM,DINPR,DNORM,ITRMKS)
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
              REAL(KIND=8) :: RCGS(N)
              REAL(KIND=8) :: RTIL(N)
              REAL(KIND=8) :: D(N)
              REAL(KIND=8) :: P(N)
              REAL(KIND=8) :: Q(N)
              REAL(KIND=8) :: U(N)
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: RWORK1(N)
              REAL(KIND=8) :: RWORK2(N)
              REAL(KIND=8) :: RSNRM
              REAL(KIND=8) :: DINPR
              EXTERNAL DINPR
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
              INTEGER(KIND=4) :: ITRMKS
            END SUBROUTINE NITTFQ
          END INTERFACE 
        END MODULE NITTFQ__genmod
