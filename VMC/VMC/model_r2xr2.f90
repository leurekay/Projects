
subroutine r2xr2_mannual( )
  use vmcplace
  implicit none

  integer nb, ib, a
  complex*16 one, unity
  real*8 hopnn, Vnn, Jnn, tnn, t2nd
  character*20 explain

  integer i, nx, ny, nhole
  character*6 vtable(20)
  real*8, target, dimension(20) :: v
  real*8, pointer :: mu, mz, hop2nd, S0, TQ, SQ, T0, ddw
  logical iEf /.true./    !switch to adjusted Ef automatically
  
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom
  logical, external :: active, inactive


  !============================================================================
  !                     constants
  one = dcmplx(0.d0, 1.d0);   unity = 1.d0
  model = 'r2xr2'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'control.'//model)
     read(10, *)explain
     read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain
     read(10, *)bruteforce;  close(10)
  end do

  !==============================================================================
  !                 parameters for lattice and physical hamiltonian 
  
  nd = 2;     ndim = nd + 1
  na1uc = 2;  nambu = na1uc*2
  
  !import nx, ny, nhole and boundary conditions from file
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     open(10, file = 'lattice.'//model)
     read(10, *)explain; read(10, *)nx, ny, nhole
     read(10, *)explain; read(10, *)antiperiodic; close(10)
  end do

  limits = (/nx, ny, na1uc/)
  nuc = limits(1) * limits(2)
  naa = nuc * na1uc
  Npair = (naa - nhole) / 2

  U = 0; mott = .true.

  tnn = 0.4;      t2nd = - 0.3d0 * tnn;  Jnn = 0.13;   Vnn = -Jnn / 4
  tnn = tnn/Jnn;  t2nd = t2nd / Jnn;    Vnn = Vnn / Jnn;  Jnn = 1

  if(.not. allocated(sdwav) ) allocate( sdwav(naa), delsdw(naa) )

  !============================================================================
  !                  variational vars
  
  !vars fixed in the present model
  gU = 1;   Vdh = 0;   Vzz = 0;   Vcc = 0;   hopnn = 1

  !vars, or part of them, to be optimized (names and pointers must match)
  vtable(1:8) = (/ 'mu', 'mz',  'hop2nd', 'S0',  'TQ',  'SQ', 'T0',  'ddw' /)
  mu => v(1);  mz => v(2);  hop2nd => v(3);  S0 => v(4)
  TQ => v(5);  SQ => v(6);  T0 => v(7);      ddw => v(8)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(8, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do

  if(optimize .and. optinfo == 0) return


  !===========================================================================
  !                 automated inputs

  unitcell%natom = na1uc

  abc => unitcell%abc;  dualabc => unitcell%dualabc
  abc(:, 1) = (/1, 0, 0/)
  abc(:, 2) = (/0, 1, 0/)
  abc(:, 3) = (/0, 0, 1/)
  call dualvector(abc(:, 2), abc(:, 3), abc(:, 1), dualabc(:, 1))  
  call dualvector(abc(:, 3), abc(:, 1), abc(:, 2), dualabc(:, 2))  
  call dualvector(abc(:, 1), abc(:, 2), abc(:, 3), dualabc(:, 3))  

  if(.not. optimize .or. optinfo == 1)allocate( unitcell%atom(na1uc) )
  unitcell%atom(1)%r = 0
  unitcell%atom(2)%r = (/0.5, 0.5, 0./)

  unitcell%atom%nbond = 4
  unitcell%atom%kondo = 0
  unitcell%atom%E = 0
  unitcell%atom%B = 0
  unitcell%atom%muloc = mu 
  unitcell%atom%mzloc = (/1, -1/) * mz      !fix the zeeman potential at atom 1 as positive for definiteness
  unitcell%atom%pairloc = one * 1.d-5

  do a = 1, na1uc;  atom => unitcell%atom(a)

     nb = atom%nbond

     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )

     atom%bond(:, 1) = (/0.5, 0.5, 0./);  atom%bond(:, 2) = (/-0.5, 0.5, 0./)
	 atom%bond(:, 3) = (/1, 0, 0/);       atom%bond(:, 4) = (/0, 1, 0/)
	    
     atom%t = (/tnn, tnn, t2nd, t2nd/);  atom%tz = 0
	 atom%J = (/1, 1, 0, 0/) * Jnn;      atom%V = (/1, 1, 0, 0/) * Vnn

	 atom%hop(1:2) = hopnn;  atom%hop(3:4) = hop2nd;   atom%hopz = 0
	 if(a == 1)then
	   atom%hop(1:2) = atom%hop(1:2) + one * ddw * (/1, -1/)
	 else
	   atom%hop(1:2) = atom%hop(1:2) - one * ddw * (/1, -1/)
     end if

	 !Q-pair first
	 atom%pair = 0;  atom%pair(1:2) = (/unity, one/) * SQ;  if(a == 2) atom%pair = - atom%pair
	 atom%pairz = 0; atom%pairz(1:2) = (/TQ, -TQ/);         if(a == 2) atom%pairz = - atom%pairz
     
	 !0-pair second
	 atom%pair(1:2) = atom%pair(1:2) + (/1, -1/) * S0 
     atom%pairz(1:2) = atom%pairz(1:2) + (/unity, one /) * T0

     atom%pairz = - atom%pairz
  end do


  if(.not. optimize .or. optinfo == 1) call neighbors(unitcell)


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
end subroutine r2xr2_mannual



subroutine r2xr2_ionic( )
  use vmcplace
  implicit none

  integer nb, ib, a
  complex*16 one, unity
  real*8 Vnn, Vion
  character*20 explain

  integer i, nx, ny
  character*6 vtable(20)
  real*8, target, dimension(20) :: v
  real*8, pointer :: mu, mz, cdw, tv2, tv3, dsc, gU_, Vdh_, Vzz_, Vcc_
  logical iEf /.false./    !switch to adjusted Ef automatically
  
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom
  logical, external :: active, inactive


  !============================================================================
  !                     constants
  one = dcmplx(0.d0, 1.d0);   unity = 1.d0
  model = 'r2xr2_ionic'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'control.'//model)
     read(10, *)explain
     read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain
     read(10, *)bruteforce;  close(10)
  end do

  !==============================================================================
  !                 parameters for lattice and physical hamiltonian 
  
  nd = 2;     ndim = nd + 1;   na1uc = 2;  nambu = na1uc*2
  
  !import nx, ny, nhole and boundary conditions from file
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     open(10, file = 'lattice.'//model)
     read(10, *)explain; read(10, *)nx, ny, npair
     read(10, *)explain; read(10, *)antiperiodic
     read(10, *)explain; read(10, *)Vion, U, Vnn;  close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2)
  naa = nuc * na1uc;           mott = .false.

  !============================================================================
  !                  variational vars
  
  !vars, or part of them, to be optimized (names and pointers must match)
  vtable(1:10) = (/ 'mu', 'mz', 'cdw', 'tv2', 'tv3', 'dsc', 'gU', 'Vdh', 'Vzz', 'Vcc'/)
  mu => v(1);   mz => v(2);   cdw => v(3);    tv2 => v(4);   tv3 => v(5)
  dsc => v(6);  gU_ => v(7);  Vdh_ => v(8);   Vzz_ => v(9);  Vcc_ => v(10)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(10, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do

  if(optimize .and. optinfo == 0) return

  !===========================================================================
  !                 automated inputs

  gU = gU_;  Vdh = Vdh_;  Vcc = Vcc_;  Vzz = Vzz_

  unitcell%natom = na1uc

  abc => unitcell%abc;  dualabc => unitcell%dualabc
  abc(:, 1) = (/1, 0, 0/)
  abc(:, 2) = (/0, 1, 0/)
  abc(:, 3) = (/0, 0, 1/)
  call dualvector(abc(:, 2), abc(:, 3), abc(:, 1), dualabc(:, 1))  
  call dualvector(abc(:, 3), abc(:, 1), abc(:, 2), dualabc(:, 2))  
  call dualvector(abc(:, 1), abc(:, 2), abc(:, 3), dualabc(:, 3))  

  if(.not. optimize .or. optinfo == 1)allocate( unitcell%atom(na1uc) )
  unitcell%atom(1)%r = 0
  unitcell%atom(2)%r = (/0.5, 0.5, 0./)

  unitcell%atom%nbond = 6
  unitcell%atom%kondo = 0
  unitcell%atom%E = (/1, -1/) * Vion
  unitcell%atom%B = 0
  unitcell%atom%muloc = (/mu-cdw, mu+cdw/)
  unitcell%atom%mzloc = (/1, -1/) * mz      
  unitcell%atom%pairloc = dcmplx(0.d0, 1.d-5)

  do a = 1, na1uc;  atom => unitcell%atom(a)

     nb = atom%nbond

     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )

     atom%bond(:, 1) = (/0.5, 0.5, 0./);  atom%bond(:, 2) = (/-0.5, 0.5, 0./)
     atom%bond(:, 3) = (/1, 0/);          atom%bond(:, 4) = (/0, 1/)
     atom%bond(:, 5) = (/1, 1/);          atom%bond(:, 6) = (/-1, 1/)
	    
     atom%t = 0;  atom%t(1:2) = 1;  atom%tz = 0;  atom%J = 0;  atom%V = Vnn

	 atom%hop = 0;  atom%hop(1:2) = 1;  atom%hopz = 0
     atom%hop(3:6) = (/tv2, tv2, tv3, tv3/)
     if(a == 2)atom%hop(3:6) = -(/tv2, tv2, tv3, tv3/)
    
     atom%pair = 0;	 atom%pair(1:2) = (/1, -1/) * dsc;  atom%pairz = 0
  end do

  if(.not. optimize .or. optinfo == 1) call neighbors(unitcell)


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
end subroutine r2xr2_ionic


