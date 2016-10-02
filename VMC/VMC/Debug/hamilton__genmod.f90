        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 20:18:05 2016
        MODULE HAMILTON__genmod
          INTERFACE 
            SUBROUTINE HAMILTON(N,HK,ND,KV,UNITCELL,SPIN)
              USE CRYSTAL
              INTEGER(KIND=4) :: ND
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: HK(N,N)
              REAL(KIND=8) :: KV(ND)
              TYPE (UNITCELLCONFIG) ,TARGET :: UNITCELL
              INTEGER(KIND=4) :: SPIN
            END SUBROUTINE HAMILTON
          END INTERFACE 
        END MODULE HAMILTON__genmod
