        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:03 2016
        MODULE DYSON__genmod
          INTERFACE 
            SUBROUTINE DYSON(N,GINV,G,ITH,VEC,ROWCOL,RATIO,INFO)
              INTEGER(KIND=4) :: N
              COMPLEX(KIND=8) :: GINV(N,N)
              COMPLEX(KIND=8) :: G(N,N)
              INTEGER(KIND=4) :: ITH
              COMPLEX(KIND=8) :: VEC(N)
              INTEGER(KIND=4) :: ROWCOL
              COMPLEX(KIND=8) :: RATIO
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DYSON
          END INTERFACE 
        END MODULE DYSON__genmod
