        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:04 2016
        MODULE ZSVD__genmod
          INTERFACE 
            SUBROUTINE ZSVD(M,N,A,U,S,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              COMPLEX(KIND=8) :: A(M,N)
              COMPLEX(KIND=8) :: U(M,M)
              REAL(KIND=8) :: S(*)
              COMPLEX(KIND=8) :: V(N,N)
            END SUBROUTINE ZSVD
          END INTERFACE 
        END MODULE ZSVD__genmod
