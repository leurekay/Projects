subroutine r2xr2Kondo( )
  use vmcplace
  implicit none

  integer nb, ib, a
  real*8  Vnn, Jnn, JK, pi, bond(2)
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom


  integer i, nx, ny
  real*8, pointer :: mz, mzK, hopK, gapK, swave, dwave
  character*6 name, vtable(20)
  character*20 explain
  real*8, target :: v(20)


  !============================================================================
  !                   constants and mc controls
  model = 'r2xr2Kondo'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'control.'//model)
     read(10, *)explain;  read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain;  read(10, *)bruteforce;   close(10)
  end do

  !============================================================================
  !               lattice, electron filling, and phyical parameters

  nd = 2;         ndim = nd + 1;         na1uc = 4;     nambu = na1uc * 2
  U = 0;          Vnn = 0;               Jnn = 0
  JK = 0.3;       mott = 0

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'lattice.'//model)
     read(10, *)explain;  read(10, *)nx, ny
     read(10, *)explain;  read(10, *)antiperiodic;  close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2)
  naa = nuc * na1uc;         npair = naa / 2
  
  !npair = naa / 2 sets the model at p-h symmetric point

  !===============================================================================
  !                variational vars 

  vtable(1:6) = (/'mz', 'mzK', 'hopK', 'gapK', 'swave', 'dwave'/)
  mz => v(1);   mzK => v(2);  hopK => v(3);  gapK => v(4);  swave => v(5);  dwave => v(6)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     call register(6, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do  
  
  if(optimize .and. optinfo == 0) return


  !=============================================================================
  !              automated inputs

  unitcell%natom = na1uc

  abc => unitcell%abc;  dualabc => unitcell%dualabc
  abc(:, 1) = (/1, 0, 0/)
  abc(:, 2) = (/0, 1, 0/)
  abc(:, 3) = (/0, 0, 1/)
  call dualvector(abc(:, 2), abc(:, 3), abc(:, 1), dualabc(:, 1))  
  call dualvector(abc(:, 3), abc(:, 1), abc(:, 2), dualabc(:, 2))  
  call dualvector(abc(:, 1), abc(:, 2), abc(:, 3), dualabc(:, 3))  

  if(.not. optimize .or. optinfo == 1) allocate( unitcell%atom(na1uc) )
  unitcell%atom(1)%r = 0
  unitcell%atom(2)%r = (/0.5, 0.5, 0./)
  unitcell%atom(3)%r = (/0., 0., .25/)
  unitcell%atom(4)%r = (/0.5, 0.5, .25/)

  unitcell%atom%nbond = (/3, 3, 0, 0/)
  unitcell%atom%kondo = (/0, 0, 1, 1/)
  unitcell%atom%E = 0
  unitcell%atom%B = 0
  unitcell%atom%muloc = 0
  unitcell%atom%pairloc = 0
  unitcell%atom%mzloc = (/mz, -mz, -mzK, mzK/)

  do a = 1, na1uc;  atom => unitcell%atom(a)

     nb = atom%nbond;  if(nb == 0)cycle

     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )

     bond = (/0.5, 0.5/)
	 do ib = 1, 2;  if(ib > 1)call rotate90(bond)
	    atom%bond(1:2, ib) = bond
     end do     
	 atom%bond(:, 3) = (/0., 0., 0.25/)
	    
     atom%tz = 0;  atom%t = (/1, 1, 0/);	 atom%J = (/Jnn, Jnn, JK/);	 atom%V = 0;  atom%V(3) = JK/4
	 atom%hopz = 0;  atom%hop(1:2) = 1;  atom%hop(3) = hopK
	 atom%pair = 0;  atom%pair(1:2) = (/swave+dwave, swave-dwave/); atom%pair(3) = gapK
	 atom%pairz = 0
  end do

  if(.not. optimize .or. optinfo == 1) call neighbors(unitcell)

  return
end subroutine r2xr2Kondo