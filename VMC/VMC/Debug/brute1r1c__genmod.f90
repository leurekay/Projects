        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:03 2016
        MODULE BRUTE1R1C__genmod
          INTERFACE 
            SUBROUTINE BRUTE1R1C(N,ITHROW,ROW,JTHCOL,COL,GINV,DET,RATIO,&
     &INFO)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ITHROW
              COMPLEX(KIND=8) :: ROW(N)
              INTEGER(KIND=4) :: JTHCOL
              COMPLEX(KIND=8) :: COL(N)
              COMPLEX(KIND=8) :: GINV(N,N)
              COMPLEX(KIND=8) :: DET
              COMPLEX(KIND=8) :: RATIO
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE BRUTE1R1C
          END INTERFACE 
        END MODULE BRUTE1R1C__genmod
