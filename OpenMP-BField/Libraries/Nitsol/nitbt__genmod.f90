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
        MODULE NITBT__genmod
          INTERFACE 
            SUBROUTINE NITBT(N,XCUR,FCNRM,STEP,ETA,XPLS,FPLS,FPNRM,OFTJS&
     &,REDFAC,NFE,IBT,IBTMAX,F,RPAR,IPAR,DNORM,ITRMBT)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: XCUR(N)
              REAL(KIND=8) :: FCNRM
              REAL(KIND=8) :: STEP(N)
              REAL(KIND=8) :: ETA
              REAL(KIND=8) :: XPLS(N)
              REAL(KIND=8) :: FPLS(N)
              REAL(KIND=8) :: FPNRM
              REAL(KIND=8) :: OFTJS
              REAL(KIND=8) :: REDFAC
              INTEGER(KIND=4) :: NFE
              INTEGER(KIND=4) :: IBT
              INTEGER(KIND=4) :: IBTMAX
              EXTERNAL F
              REAL(KIND=8) :: RPAR(*)
              INTEGER(KIND=4) :: IPAR(*)
              REAL(KIND=8) :: DNORM
              EXTERNAL DNORM
              INTEGER(KIND=4) :: ITRMBT
            END SUBROUTINE NITBT
          END INTERFACE 
        END MODULE NITBT__genmod
