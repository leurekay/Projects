!**************************************************************************************************
subroutine broadcast_cmplx(seed)
  use MPI_place
  implicit none
  complex*16 seed
  call MPI_BCAST(seed, 1, mpi_double_complex, 0, mpi_comm_world, ierr )          !MPI statement
  return
end subroutine broadcast_cmplx

!**************************************************************************************************
subroutine broadcast_dble(seed)
  use MPI_place
  implicit none
  real*8 seed
  call MPI_BCAST(seed, 1, mpi_double_precision, 0, mpi_comm_world, ierr )          !MPI statement
  return
end subroutine broadcast_dble

!**************************************************************************************************
subroutine broadcast_integer(seed)
  use MPI_place
  implicit none
  integer seed
  call MPI_BCAST(seed, 1, mpi_integer, 0, mpi_comm_world, ierr )          !MPI statement
  return
end subroutine broadcast_integer
!**************************************************************************************************
subroutine broadcast_integer_vector(N, A)
  use MPI_place
  implicit none
  integer N, Newtype
  integer A(N)

  !define a structure, Newtype, with N elements of required precision
  call MPI_TYPE_CONTIGUOUS(N, mpi_integer, Newtype, ierr)          !MPI statement
  call MPI_TYPE_COMMIT(Newtype, ierr)                              !MPI statement

  !broadcast A as 1 data of Newtype
  call MPI_BCAST(A, 1, Newtype, 0, mpi_comm_world, ierr )          !MPI statement

  !free the structure declaration
  call MPI_TYPE_FREE(Newtype, ierr)                                !MPI statement
  return
end subroutine broadcast_integer_vector


!**************************************************************************************************
subroutine broadcast_cmplx_vector(N, A)
  use MPI_place
  implicit none
  integer N, Newtype
  complex*16 A(N)

  !define a structure, Newtype, with N elements of required precision
  call MPI_TYPE_CONTIGUOUS(N, mpi_double_complex, Newtype, ierr)     !MPI statement
  call MPI_TYPE_COMMIT(Newtype, ierr)                                !MPI statement

  !broadcast A as 1 data of Newtype
  call MPI_BCAST(A, 1, Newtype, 0, mpi_comm_world, ierr )            !MPI statement

  !free the structure declaration
  call MPI_TYPE_FREE(Newtype, ierr)                                  !MPI statement
  return
end subroutine broadcast_cmplx_vector


!**************************************************************************************************
subroutine broadcast_dble_vector(N, A)
  use MPI_place
  implicit none
  integer N, Newtype
  real*8 A(N)

  !define a structure, Newtype, with N elements of required precision
  call MPI_TYPE_CONTIGUOUS(N, mpi_double_precision, Newtype, ierr)   !MPI statement
  call MPI_TYPE_COMMIT(Newtype, ierr)                                !MPI statement

  !broadcast A as 1 data of Newtype
  call MPI_BCAST(A, 1, Newtype, 0, mpi_comm_world, ierr )            !MPI statement

  !free structure declaration
  call MPI_TYPE_FREE(Newtype, ierr)                                  !MPI statement
  return
end subroutine broadcast_dble_vector

!**************************************************************************************************

subroutine broadcast_cmplx_matrix(N, M, A)
  use MPI_place
  implicit none
  integer N, M, Newtype
  complex*16 A(N, M)

  !define a structure, Newtype, with N elements of required precision
  call MPI_TYPE_CONTIGUOUS(N*M, mpi_double_complex, Newtype, ierr)   !MPI statement
  call MPI_TYPE_COMMIT(Newtype, ierr)                                !MPI statement

  !broadcast A as 1 data of Newtype
  call MPI_BCAST(A, 1, Newtype, 0, mpi_comm_world, ierr )            !MPI statement

  !free structure declaration
  call MPI_TYPE_FREE(Newtype, ierr)                                  !MPI statement
  return
end subroutine broadcast_cmplx_matrix


!**************************************************************************************************
subroutine broadcast_dble_matrix(N, M, A)
  use MPI_place
  implicit none
  integer N, M, Newtype
  real*8 A(N, M)
  
  !define a structure, Newtype, with N elements of required precision  
  call MPI_TYPE_CONTIGUOUS(N*M, mpi_double_precision, Newtype, ierr)   !MPI statement
  call MPI_TYPE_COMMIT(Newtype, ierr)                                  !MPI statement

  !broadcast A as 1 data of Newtype
  call MPI_BCAST(A, 1, Newtype, 0, mpi_comm_world, ierr )              !MPI statement

  !free structure declaration
  call MPI_TYPE_FREE(Newtype, ierr)                                    !MPI statement
  return
end subroutine broadcast_dble_matrix


!**************************************************************************************************
subroutine broadcast_integer_matrix(N, M, A)
  use MPI_place
  implicit none
  integer N, M, Newtype, A(N, M)

  !define a structure, Newtype, with N elements of required precision
  call MPI_TYPE_CONTIGUOUS(N*M, mpi_integer, Newtype, ierr)          !MPI statement
  call MPI_TYPE_COMMIT(Newtype, ierr)                                !MPI statement

  !broadcast A as 1 data of Newtype
  call MPI_BCAST(A, 1, Newtype, 0, mpi_comm_world, ierr )            !MPI statement

  !free structure declaration
  call MPI_TYPE_FREE(Newtype, ierr)                                  !MPI statement
  return
end subroutine broadcast_integer_matrix

!**************************************************************************************************
subroutine reduce_cmplx_matrix(N, M, A)
  use MPI_place
  implicit none
  integer N, M
  complex*16 A(N, M), B(N, M)
  B = 0; call MPI_REDUCE(A, B, N*M, mpi_double_complex, MPI_SUM, 0, mpi_comm_world, ierr )   !MPI statement
  A = B; call MPI_BCAST(A, N*M, mpi_double_complex, 0, mpi_comm_world, ierr )                !MPI statement
  return
end subroutine reduce_cmplx_matrix    

!**************************************************************************************************
subroutine reduce_dble_matrix(N, M, A)
  use MPI_place
  implicit none
  integer N, M
  real*8 A(N, M), B(N, M)
  B = 0; call MPI_REDUCE(A, B, N*M, mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr )   !MPI statement
  A = B; call MPI_BCAST(A, N*M, mpi_double_precision, 0, mpi_comm_world, ierr )                !MPI statement
  return
end subroutine reduce_dble_matrix    
!**************************************************************************************************

subroutine reduce_integer_matrix(N, M, A)
  use MPI_place
  implicit none
  integer N, M
  integer A(N, M), B(N, M)
  B = 0; call MPI_REDUCE(A, B, N*M, mpi_integer, MPI_SUM, 0, mpi_comm_world, ierr )            !MPI statement
  A = B; call MPI_BCAST(A, N*M, mpi_integer, 0, mpi_comm_world, ierr )                         !MPI statement
  return
end subroutine reduce_integer_matrix    
!**************************************************************************************************

subroutine reduce_cmplx_vector(N, A)
  use MPI_place
  implicit none
  integer N
  complex*16 A(N), B(N)
  B = 0; call MPI_REDUCE(A, B, N, mpi_double_complex, MPI_SUM, 0, mpi_comm_world, ierr )     !MPI statement
  A = B; call MPI_BCAST(A, N, mpi_double_complex, 0, mpi_comm_world, ierr )                  !MPI statement
  return
end subroutine reduce_cmplx_vector    

!**************************************************************************************************
subroutine reduce_dble_vector(N, A)
  use MPI_place
  implicit none
  integer N
  real*8 A(N), B(N)
  B = 0; call MPI_REDUCE(A, B, N, mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr )    !MPI statement
  A = B; call MPI_BCAST(A, N, mpi_double_precision, 0, mpi_comm_world, ierr )                 !MPI statement
  return
end subroutine reduce_dble_vector    

!**************************************************************************************************
subroutine reduce_integer_vector(N, A)
  use MPI_place
  implicit none
  integer N
  integer A(N), B(N)
  B = 0; call MPI_REDUCE(A, B, N, mpi_integer, MPI_SUM, 0, mpi_comm_world, ierr )      !MPI statement
  A = B; call MPI_BCAST(A, N, mpi_integer, 0, mpi_comm_world, ierr )                   !MPI statement
  return
end subroutine reduce_integer_vector


subroutine reduce_cmplx(A)
  use MPI_place
  implicit none
  complex*16 A, B
  B = 0; call MPI_REDUCE(A, B, 1, mpi_double_complex, MPI_SUM, 0, mpi_comm_world, ierr )     !MPI statement
  A = B; call MPI_BCAST(A, 1, mpi_double_complex, 0, mpi_comm_world, ierr )                  !MPI statement
  return
end subroutine reduce_cmplx   

!**************************************************************************************************
subroutine reduce_dble(A)
  use MPI_place
  implicit none
  real*8 A, B
  B = 0; call MPI_REDUCE(A, B, 1, mpi_double_precision, MPI_SUM, 0, mpi_comm_world, ierr )    !MPI statement
  A = B; call MPI_BCAST(A, 1, mpi_double_precision, 0, mpi_comm_world, ierr )                 !MPI statement
  return
end subroutine reduce_dble   

!**************************************************************************************************
subroutine reduce_integer(A)
  use MPI_place
  implicit none
  integer A, B
  B = 0; call MPI_REDUCE(A, B, 1, mpi_integer, MPI_SUM, 0, mpi_comm_world, ierr )      !MPI statement
  A = B; call MPI_BCAST(A, 1, mpi_integer, 0, mpi_comm_world, ierr )                   !MPI statement
  return
end subroutine reduce_integer
