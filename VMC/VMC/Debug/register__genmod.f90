        !COMPILER-GENERATED INTERFACE MODULE: Wed Sep 28 19:05:05 2016
        MODULE REGISTER__genmod
          INTERFACE 
            SUBROUTINE REGISTER(N,TABLE,V,NVAR,VNAME,VMIN,VMAX,VAR,MODEL&
     &,OPTIMIZE)
              INTEGER(KIND=4) :: N
              CHARACTER(LEN=6) :: TABLE(*)
              REAL(KIND=8) :: V(*)
              INTEGER(KIND=4) :: NVAR
              CHARACTER(LEN=6) :: VNAME(*)
              REAL(KIND=8) :: VMIN(*)
              REAL(KIND=8) :: VMAX(*)
              REAL(KIND=8) :: VAR(*)
              CHARACTER(LEN=12) :: MODEL
              LOGICAL(KIND=4) :: OPTIMIZE
            END SUBROUTINE REGISTER
          END INTERFACE 
        END MODULE REGISTER__genmod
