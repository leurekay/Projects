        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:03 2016
        MODULE ROWCOL__genmod
          INTERFACE 
            SUBROUTINE ROWCOL(N,ITHROW,ROW,JTHCOL,COL,GINV,G,RATIO,INFO)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ITHROW
              COMPLEX(KIND=8) :: ROW(N)
              INTEGER(KIND=4) :: JTHCOL
              COMPLEX(KIND=8) :: COL(N)
              COMPLEX(KIND=8) :: GINV(N,N)
              COMPLEX(KIND=8) :: G(N,N)
              COMPLEX(KIND=8) :: RATIO
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE ROWCOL
          END INTERFACE 
        END MODULE ROWCOL__genmod
