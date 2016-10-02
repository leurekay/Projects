subroutine cobyla_opt()
  use cobyla2;  use vmcplace
  implicit none

  integer, parameter :: iprint  = 3
  integer maxfun/100/
  real*8 rhobeg, rhoend
  real*8, allocatable :: x(:)

  integer n, m, i

  !caution: for the time being, cobyla_opt can only be called in the main program because of 
  !cobyla itself is not sufficiently user-friendly

  LMoptimize = 0;  optimize = 1;  optinfo = 0;  call definejob( )

  n = nvar;  if(n > 5)stop 'too many vars to be optimized by cobyla'
  allocate( x(n) )

  !cobyla parameters
  m = 2 * n;  rhobeg = 0.5d0;  rhoend = 1.d-3
  
  !initial value of x, corresponding to mid points of vars in bounds user specifed
  x = 0.5

  if(ithread == 0)open(17, file = 'histo.'//model)
  
  !cobyla  optimization
  call cobyla (n, m, x, rhobeg, rhoend, iprint, maxfun)

  if(ithread /= 0)return

  close(17)
  
  !tranport x to var used in MCthread
  var(1:n) =  vmin(1:n) + (  vmax(1:n) -  vmin(1:n) ) * x

  open(10, file = 'optimized.'//model)

  write(10, *)'optimized variables:'
  do i = 1, n
	 write(10, *) vname(i), ' = ',  var(i)
  end do
  close(10)

  return
end subroutine cobyla_opt


subroutine calcfc(n, m, x, f, con)
  use vmcplace
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: n, m
  REAL*8, INTENT(IN)  :: x(:)
  REAL*8, INTENT(OUT) :: f
  REAL*8, INTENT(OUT) :: con(:)

  integer i
  
  call MPI_barrier(mpi_comm_world, ierr)

  !tranport x to var used in MCthread
  var(1:n) =  vmin(1:n) + (  vmax(1:n) -  vmin(1:n) ) * x

  !count times of optimization performed
  optinfo = optinfo + 1
  !if optinfo > 1, dynamic memories are no longer allocated but only updated, e.g., in unitcell

  !define the complete model and its infrastructure  
  call definejob();  call infrastructure( )

  call MCthread();  call MPI_average( );  f = Eav

  !constrains for vars normalized by their spans
  con(1:n) = x;   con(n+1:m) = 1 - x 

  if(ithread == 0)write(17, 100) var(1:n), f, delE, rhoav, delrho

100 format(1x, 9f12.5)
  return
end subroutine calcfc



