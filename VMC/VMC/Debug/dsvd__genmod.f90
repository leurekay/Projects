        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:04 2016
        MODULE DSVD__genmod
          INTERFACE 
            SUBROUTINE DSVD(M,N,A,U,S,V)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(M,N)
              REAL(KIND=8) :: U(M,M)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: V(N,N)
            END SUBROUTINE DSVD
          END INTERFACE 
        END MODULE DSVD__genmod
