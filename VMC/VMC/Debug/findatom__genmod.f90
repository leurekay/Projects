        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 20:18:06 2016
        MODULE FINDATOM__genmod
          INTERFACE 
            FUNCTION FINDATOM(RV,UNITCELL,RT)
              USE CRYSTAL
              REAL(KIND=8) :: RV(3)
              TYPE (UNITCELLCONFIG) ,TARGET :: UNITCELL
              REAL(KIND=8) :: RT(3)
              INTEGER(KIND=4) :: FINDATOM
            END FUNCTION FINDATOM
          END INTERFACE 
        END MODULE FINDATOM__genmod
