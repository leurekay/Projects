        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 20:18:05 2016
        MODULE BANDSTRUCTURE__genmod
          INTERFACE 
            SUBROUTINE BANDSTRUCTURE(ND,NA,GM,GX,UNITCELL)
              USE CRYSTAL
              INTEGER(KIND=4) :: NA
              INTEGER(KIND=4) :: ND
              REAL(KIND=8) :: GM(ND)
              REAL(KIND=8) :: GX(ND)
              TYPE (UNITCELLCONFIG) ,TARGET :: UNITCELL
            END SUBROUTINE BANDSTRUCTURE
          END INTERFACE 
        END MODULE BANDSTRUCTURE__genmod