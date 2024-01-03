        !COMPILER-GENERATED INTERFACE MODULE: Wed Jan  3 22:40:05 2024
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITJV__genmod
          INTERFACE 
            SUBROUTINE NITJV(N,XCUR,FCUR,F,JACV,RPAR,IPAR,IJACV,IFDORD, &
     &ITASK,NFE,NJVE,NRPRE,V,Z,RWORK1,RWORK2,DNORM,ITRMJV)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XCUR(N)
              REAL(KIND=8) :: FCUR(N)
              EXTERNAL F
              EXTERNAL JACV
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              INTEGER(KIND=4) :: IJACV
              INTEGER(KIND=4) :: IFDORD
              INTEGER(KIND=4) :: ITASK
              INTEGER(KIND=4) :: NFE
              INTEGER(KIND=4) :: NJVE
              INTEGER(KIND=4) :: NRPRE
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: Z(N)
              REAL(KIND=8) :: RWORK1(N)
              REAL(KIND=8) :: RWORK2(N)
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
              INTEGER(KIND=4) :: ITRMJV
            END SUBROUTINE NITJV
          END INTERFACE 
        END MODULE NITJV__genmod
