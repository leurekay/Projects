        !COMPILER-GENERATED INTERFACE MODULE: Thu Sep 29 20:18:06 2016
        MODULE CONNECTION__genmod
          INTERFACE 
            SUBROUTINE CONNECTION(RI,IA,BOND,RJ,JA,INFO)
              USE VMCPLACE
              INTEGER(KIND=4) :: RI(ND)
              INTEGER(KIND=4) :: IA
              REAL(KIND=8) :: BOND(3)
              INTEGER(KIND=4) :: RJ(ND)
              INTEGER(KIND=4) :: JA
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE CONNECTION
          END INTERFACE 
        END MODULE CONNECTION__genmod
