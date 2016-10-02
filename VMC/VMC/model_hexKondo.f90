subroutine hexKondo( )
  use vmcplace
  implicit none

  integer nb, ib, a, i
  real*8 JK, Vnn, bond(2)
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom

  integer nx, ny
  character*6 vtable(20)
  character*20 explain
  real*8, target :: v(20)
  real*8, pointer :: mz, mzK, hopK, gapK


  !============================================================================
  !                    constants  and mc controls
  model = 'hexKondo'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'control.'//model)
     read(10, *)explain;  read(10, *)nbath, nsample, nmc1sample
     read(10, *)explain;  read(10, *)bruteforce;  close(10)
  end do

  !==============================================================================
  !                lattice, electron filling and physical parameters

  nd = 2;  ndim = nd + 1;  na1uc = 4;   nambu = na1uc * 2
  U = 0;   Vnn = 0;        mott = 0

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'lattice.'//model)
     read(10, *)explain;  read(10, *)nx, ny
     read(10, *)explain;  read(10, *)antiperiodic
     read(10, *)explain;  read(10, *)JK;  close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2)
  naa = nuc * na1uc;           npair = naa / 2

  if(.not. allocated(occup) )then
    allocate( occup(naa, 2) ); occup  = 0
    do i = 1, nuc
       occup(i, 1) = i; occup(i + nuc, 2) = i                     !spin up/dn on A/B sublattice
       occup(i + 2*nuc, 2) = i + nuc; occup(i+3*nuc, 1) = i + nuc !spin dn/up on C/D sublattice 
    end do
  end if  

  !===============================================================================
  !                variational

  vtable(1:4) = (/'mz', 'mzK', 'hopK', 'gapK'/)
  mz => v(1);   mzK => v(2);  hopK => v(3);  gapK => v(4)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(4, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do

  if(.not. allocated(sdwav) ) allocate( sdwav(naa), delsdw(naa) )
  if(.not. allocated(cdwav) ) allocate( cdwav(naa), delcdw(naa) )

  if(optimize .and. optinfo == 0) return

  !================================================================================
  !              automated inputs

  unitcell%natom = na1uc

  abc => unitcell%abc;  dualabc => unitcell%dualabc
  abc(:, 1) = (/1, 0, 0/)
  abc(:, 2) = (/1., sqrt(3.), 0./)/2
  abc(:, 3) = (/0, 0, 1/)
  call dualvector(abc(:, 2), abc(:, 3), abc(:, 1), dualabc(:, 1))  
  call dualvector(abc(:, 3), abc(:, 1), abc(:, 2), dualabc(:, 2))  
  call dualvector(abc(:, 1), abc(:, 2), abc(:, 3), dualabc(:, 3))  

  if(.not. optimize .or. optinfo == 1) allocate( unitcell%atom(na1uc) )
  unitcell%atom(1)%r = 0
  unitcell%atom(2)%r = (/0.5, 0.5/sqrt(3.), 0./)
  unitcell%atom(3)%r = (/0., 0., .5/)
  unitcell%atom(4)%r = (/0.5, 0.5/sqrt(3.), .5/)

  unitcell%atom%nbond = (/4, 1, 0, 0/)
  unitcell%atom%kondo = (/0, 0, 1, 1/)
  unitcell%atom%E = 0
  unitcell%atom%B = 0
  unitcell%atom%muloc = 0
  unitcell%atom%pairloc = 0
  unitcell%atom%mzloc = (/mz, -mz, -mzK, mzK/)

  do a = 1, na1uc;  atom => unitcell%atom(a)

     nb = atom%nbond;  if(nb == 0)cycle

     if(.not. optimize .or. optinfo == 1)allocate( atom%bond(3, nb) );  atom%bond = 0
     if(.not. optimize .or. optinfo == 1)allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1)allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )

     if(a == 2)then
	   atom%bond(:, 1) = (/0., 0., 0.5/)
	   atom%t = 0;    atom%tz = 0;  atom%J = JK;  atom%V = 0
	   atom%hop = hopK;  atom%hopz = 0;  atom%pair = gapK;  atom%pairz  = 0
	 else if(a == 1)then

       bond = (/0.5, 0.5/sqrt(3.)/)
	   do ib = 1, 3;  if(ib > 1)call rotate120(bond)
	      atom%bond(1:2, ib) = bond
       end do     
	   atom%bond(:, 4) = (/0., 0., 0.5/)
	    
       atom%tz = 0;         atom%t = (/1, 1, 1, 0/)
	   atom%J(1:3) = 0;     atom%J(4) = JK;	     atom%V(1:3) = Vnn;   atom%V(4) = 0
	   atom%hopz = 0;       atom%hop(1:3) = 1;   atom%hop(4) = hopK
	   atom%pairz = 0;      atom%pair = 0;       atom%pair(4) = gapK
     end if
  end do

  if(.not. optimize .or. optinfo == 1)call neighbors(unitcell)

  return
end subroutine hexKondo