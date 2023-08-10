        !COMPILER-GENERATED INTERFACE MODULE: Thu Aug 10 00:01:37 2023
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
