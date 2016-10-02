        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 13:27:01 2016
        MODULE BISEC__genmod
          INTERFACE 
            SUBROUTINE BISEC(NV,VMIN,VMAX,NPIXEL,CONVERGENCE,USPEC)
              INTEGER(KIND=4), INTENT(IN) :: NV
              REAL(KIND=8), INTENT(IN) :: VMIN(NV)
              REAL(KIND=8), INTENT(IN) :: VMAX(NV)
              INTEGER(KIND=4), INTENT(IN) :: NPIXEL
              REAL(KIND=8), INTENT(IN) :: CONVERGENCE
              INTERFACE 
                FUNCTION USPEC(NV,V)
                  INTEGER(KIND=4) :: NV
                  REAL(KIND=8) :: V(NV)
                  REAL(KIND=8) :: USPEC
                END FUNCTION USPEC
              END INTERFACE 
            END SUBROUTINE BISEC
          END INTERFACE 
        END MODULE BISEC__genmod
