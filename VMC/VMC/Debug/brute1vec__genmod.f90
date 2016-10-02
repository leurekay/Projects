        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:03 2016
        MODULE BRUTE1VEC__genmod
          INTERFACE 
            SUBROUTINE BRUTE1VEC(N,GINV,ITH,VEC,ROWCOL,DET,RATIO,INFO)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: GINV(N,N)
              INTEGER(KIND=4) :: ITH
              COMPLEX(KIND=8) :: VEC(N)
              INTEGER(KIND=4) :: ROWCOL
              COMPLEX(KIND=8) :: DET
              COMPLEX(KIND=8) :: RATIO
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE BRUTE1VEC
          END INTERFACE 
        END MODULE BRUTE1VEC__genmod
