<<<<<<< HEAD
<<<<<<< HEAD
        !COMPILER-GENERATED INTERFACE MODULE: Wed Jan 17 21:00:37 2024
=======
        !COMPILER-GENERATED INTERFACE MODULE: Sat Jan 20 00:25:37 2024
>>>>>>> RF_test
=======
        !COMPILER-GENERATED INTERFACE MODULE: Fri Feb 16 19:14:11 2024
>>>>>>> MCC
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NITSOL__genmod
          INTERFACE 
            SUBROUTINE NITSOL(N,X,F,JACV,FTOL,STPTOL,INPUT,INFO,RWORK,  &
     &RPAR,IPAR,ITERM,DINPR,DNORM)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              EXTERNAL F
              EXTERNAL JACV
              REAL(KIND=8) :: FTOL
              REAL(KIND=8) :: STPTOL
              INTEGER(KIND=4) :: INPUT(10)
              INTEGER(KIND=4) :: INFO(6)
              REAL(KIND=8) :: RWORK(*)
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              INTEGER(KIND=4) :: ITERM
              REAL(KIND=8) :: DINPR
              EXTERNAL DINPR
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
            END SUBROUTINE NITSOL
          END INTERFACE 
        END MODULE NITSOL__genmod
