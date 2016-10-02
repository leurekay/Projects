        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 20:18:05 2016
        MODULE PAIRING__genmod
          INTERFACE 
            SUBROUTINE PAIRING(N,GK,ND,KV,UNITCELL)
              USE CRYSTAL
              INTEGER(KIND=4) :: ND
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: GK(N,N)
              REAL(KIND=8) :: KV(ND)
              TYPE (UNITCELLCONFIG) ,TARGET :: UNITCELL
            END SUBROUTINE PAIRING
          END INTERFACE 
        END MODULE PAIRING__genmod
