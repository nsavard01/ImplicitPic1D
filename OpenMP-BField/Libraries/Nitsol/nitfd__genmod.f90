<<<<<<< HEAD
<<<<<<< HEAD
        !COMPILER-GENERATED INTERFACE MODULE: Wed Jan 17 21:00:38 2024
=======
        !COMPILER-GENERATED INTERFACE MODULE: Sat Jan 20 00:25:38 2024
>>>>>>> RF_test
=======
        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 16 19:14:11 2024
>>>>>>> MCC
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITFD__genmod
          INTERFACE 
            SUBROUTINE NITFD(N,XCUR,FCUR,F,RPAR,IPAR,IJACV,IFDORD,NFE,V,&
     &Z,RWORK,DNORM,ITRMJV)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XCUR(N)
              REAL(KIND=8) :: FCUR(N)
              EXTERNAL F
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              INTEGER(KIND=4) :: IJACV
              INTEGER(KIND=4) :: IFDORD
              INTEGER(KIND=4) :: NFE
              REAL(KIND=8) :: V(N)
              REAL(KIND=8) :: Z(N)
              REAL(KIND=8) :: RWORK(N)
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
              INTEGER(KIND=4) :: ITRMJV
            END SUBROUTINE NITFD
          END INTERFACE 
        END MODULE NITFD__genmod
