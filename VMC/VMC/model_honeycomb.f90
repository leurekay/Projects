
subroutine honeycomb( )
  use vmcplace
  implicit none

  integer nb, ib, a
  real*8  pi, Vnn, Jnn, bond(2)
  complex*16 one
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom

  integer nx, ny, i
  character*6 vtable(20)
  character*20 explain
  real*8, target :: v(20)
  real*8, pointer :: mu, mz, did, gU_, Vdh_, Vzz_, Vcc_
  logical iEf/.true./

  !=======================================================================
  !                   constants
  one = dcmplx(0.d0, 1.d0);   pi = asin(1.d0)
  model = 'honeycomb'

  !============================================================================
  !                   constants and mc controls
  
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     open(10, file = 'control.'//model)
     read(10, *)explain;  read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain;  read(10, *)bruteforce;   close(10)
  end do

  !============================================================================
  !               lattice, electron filling, and phyical parameters

  nd = 2;      ndim = nd + 1;      na1uc = 2;    nambu = na1uc * 2

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'lattice.'//model)
     read(10, *)explain;  read(10, *)nx, ny, npair
     read(10, *)explain;  read(10, *)antiperiodic
     read(10, *)explain;  read(10, *)U, Vnn, Jnn
     read(10, *)explain;  read(10, *) mott; close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2);   naa = nuc * na1uc
  
  !van hove point: npair = naa * 3 / 8         


  !===============================================================================
  !                variational vars 

  vtable(1:7) = (/'mu', 'mz', 'did', 'gU', 'Vdh', 'Vzz', 'Vcc'/)
  mu => v(1);  mz => v(2);  did => v(3);  gU_ => v(4);  Vdh_ => v(5);  Vzz_ => V(6); Vcc_ => v(7)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(7, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do  
  
  if(optimize .and. optinfo == 0) return

  !notice:  mz represents afm order

  !==========================================================================
  !                     automated inputs
  gU = gU_;  Vdh = Vdh_;  Vzz = Vzz_;   Vcc =  Vcc_

  unitcell%natom = na1uc

  abc => unitcell%abc;  dualabc => unitcell%dualabc
  abc(:, 1) = (/1, 0, 0/)
  abc(:, 2) = (/1., sqrt(3.), 0./)/2
  abc(:, 3) = (/0, 0, 1/)
  call dualvector(abc(:, 2), abc(:, 3), abc(:, 1), dualabc(:, 1))  
  call dualvector(abc(:, 3), abc(:, 1), abc(:, 2), dualabc(:, 2))  
  call dualvector(abc(:, 1), abc(:, 2), abc(:, 3), dualabc(:, 3))  

  if(.not. optimize .or. optinfo ==1 ) allocate( unitcell%atom(na1uc) )
  unitcell%atom(1)%r = 0
  unitcell%atom(2)%r = (/0.5, 0.5/sqrt(3.), 0./)
  unitcell%atom%nbond = (/3, 0/)
  unitcell%atom%kondo = 0
  unitcell%atom%E = 0
  unitcell%atom%B = 0
  unitcell%atom%muloc = mu
  unitcell%atom%mzloc = (/1, -1/) * mz

  do a = 1, na1uc;  atom => unitcell%atom(a)

     nb = atom%nbond;  if(nb == 0)cycle

     if(.not. optimize .or. optinfo ==1 ) allocate( atom%bond(3, nb) );  atom%bond = 0
     if(.not. optimize .or. optinfo ==1 ) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo ==1 ) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )

     bond = (/0.5, 0.5/sqrt(3.)/)
	 do ib = 1, nb;  if(ib > 1)call rotate120(bond)
	    atom%bond(1:2, ib) = bond
     end do     
	    
     atom%t=1;  atom%tz = 0; atom%J=0;  atom%V = Vnn  
	 atom%hop = 1;  atom%hopz = 0;      atom%pairz = 0
	 atom%pair = exp( one * (/0, 2, 4/) * pi / 3 ) * did 
  
  end do

  if(.not. optimize  .or. optinfo == 1) call neighbors(unitcell)

  !=====================================================================================
  !             adjust fermilevel if applicable

  do i = 1, nvar; if( vname(i) == 'mu' ) iEf = .false.;  end do
  !if mu is to be optimized, do not adjust it automatically

  if(.not. iEf )return
  !even if mu is not present in vars to be optimized, or vmc is not in the optimization mode
  !one can still decide whether fermilevel is to be adjusted. 

  call adjustFermilevel( )

  if(ithread /= 0)return

  open(10, file = 'warning.'//model)
  write(10, *)'Warning: muloc adjusted automatically as'
  write(10, *)unitcell%atom%muloc; close(10)
  !if fermilevel adjusted automatically, a warning message is saved for user

  return
end subroutine honeycomb



