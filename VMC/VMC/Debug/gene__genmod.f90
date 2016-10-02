        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 13:27:01 2016
        MODULE GENE__genmod
          INTERFACE 
            SUBROUTINE GENE(NSEG,NHEAD,NGEN,ISEED,PMUTE,CONVERGENCE,    &
     &USPEC)
              INTEGER(KIND=4) :: NHEAD
              INTEGER(KIND=4) :: NSEG
              INTEGER(KIND=4) :: NGEN
              INTEGER(KIND=4) :: ISEED
              REAL(KIND=8) :: PMUTE
              REAL(KIND=8) :: CONVERGENCE
              INTERFACE 
                SUBROUTINE USPEC(NSEG,NHEAD,CODE,F)
                  INTEGER(KIND=4) :: NHEAD
                  INTEGER(KIND=4) :: NSEG
                  INTEGER(KIND=4) :: CODE(NSEG,NHEAD)
                  REAL(KIND=8) :: F(NHEAD)
                END SUBROUTINE USPEC
              END INTERFACE 
            END SUBROUTINE GENE
          END INTERFACE 
        END MODULE GENE__genmod
