        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 13:27:03 2016
        MODULE CALCFC__genmod
          INTERFACE 
            SUBROUTINE CALCFC(N,M,X,F,CON)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: M
              REAL(KIND=8), INTENT(IN) :: X(:)
              REAL(KIND=8), INTENT(OUT) :: F
              REAL(KIND=8), INTENT(OUT) :: CON(:)
            END SUBROUTINE CALCFC
          END INTERFACE 
        END MODULE CALCFC__genmod
