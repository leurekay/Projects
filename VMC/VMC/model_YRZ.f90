subroutine YRZ
  use vmcplace
  implicit none

  integer nb, ib, a
  real*8 Vnn, Jnn, tnn, t2nd
  character*20 explain

  integer i, j, k, nx, ny, nhole, iup, idn, iuc
  character*6 vtable(20)
  real*8, target, dimension(20) :: v
  real*8, pointer :: mu, h2nd, hperp, hresv
  logical iEf /.true./    !switch to adjust Ef automatically
  
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom
  logical, external :: active, inactive


  !============================================================================
  !                     constants
  if( model /= 'YRZ' ) stop 'model mismatch'

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

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2)
  naa = nuc * na1uc;           Npair = (naa - nhole) / 2

  U = 0; mott = .true.

  tnn = 0.4;      t2nd = - 0.3d0 * tnn;  Jnn = 0.13;   Vnn = -Jnn / 4
  tnn = tnn/Jnn;  t2nd = t2nd / Jnn;    Vnn = Vnn / Jnn;  Jnn = 1

  if(.not. allocated(sdwav) ) allocate( sdwav(naa), delsdw(naa) )
  if(.not. allocated(cdwav) ) allocate( cdwav(naa), delcdw(naa) )  

  if(.not. allocated(occup) )then
    allocate( occup(naa, 2) ); occup = 0; iup = 0;  idn = 0

    !staggered spin on reservior layer
    do i = 1, nx; do j = 1, ny;  iuc = i + (j-1) * nx
       k = mod(i + j, 2)
       if(k == 1)then
         idn = idn + 1;  occup(iuc+nuc, 2) = idn
       else
         iup = iup + 1;  occup(iuc+nuc, 1) = iup
       end if       
    end do;  end do

    !staggered spin on active layer, with nholes
    do i = 1, nx; do j = 1, ny;  iuc = i + (j-1) * nx
       k = mod(i + j, 2)
       if(k == 1)then
         iup = iup + 1;  if(iup <= npair) occup(iuc, 1) = iup
       else
         idn = idn + 1;  if(idn <= npair) occup(iuc, 2) = idn
       end if       
    end do;  end do
  end if  

  !============================================================================
  !                  variational vars
  
  !vars fixed in the present model
  gU = 1;   Vdh = 0;   Vzz = 0;   Vcc = 0

  !vars, or part of them, to be optimized (names and pointers must match)
  vtable(1:4) = (/ 'mu', 'h2nd', 'hperp', 'hresv' /)
  mu => v(1);  h2nd => v(2);  hperp => v(3);  hresv => v(4)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(4, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
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
  unitcell%atom(2)%r = (/0, 0, 1/) * 0.5d0

  unitcell%atom%nbond = (/8, 2/)
  unitcell%atom%kondo = (/0, 1/)
  unitcell%atom%E = 0
  unitcell%atom%B = 0
  unitcell%atom%muloc = (/1, 0/) * mu
  unitcell%atom%mzloc = 0
  unitcell%atom%pairloc = dcmplx(0.d0, 1.d0) * 1.d-5

  do a = 1, na1uc;  atom => unitcell%atom(a)

     nb = atom%nbond

     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )
     atom%pair = 0;  atom%pairz = 0

     select case (a)
       case (1)
       atom%bond(:, 1) = (/1, 0, 0/);  atom%bond(:, 2) = (/0, 1, 0/)
       atom%bond(:, 3) = (/1, 1, 0/);  atom%bond(:, 4) = (/-1, 1, 0/)
       atom%bond(:, 5) = (/1.d0, 0.d0, 0.5d0/);  atom%bond(:, 6) = (/-1.d0, 0.d0, 0.5d0/)
       atom%bond(:, 7) = (/0.d0, 1.d0, 0.5d0/);  atom%bond(:, 8) = (/0.d0, -1.d0, 0.5d0/)	    
       atom%t = (/tnn, tnn, t2nd, t2nd, 0.d0, 0.d0, 0.d0, 0.d0/);  atom%tz = 0
	   atom%J = (/1, 1, 0, 0, 0, 0, 0, 0/) * Jnn;      atom%V = (/1, 1, 0, 0, 0, 0, 0, 0/) * Vnn

	   atom%hop(1:2) = 1;  atom%hop(3:4) = h2nd
       atom%hop(5:8) = (/1, 1, -1, -1/) * hperp; atom%hopz = 0
       case(2)
       atom%bond(:, 1) = (/1, 0, 0/);  atom%bond(:, 2) = (/0, 1, 0/)
       atom%t = 0;  atom%tz = 0
	   atom%J = 0;  atom%V = 0
	   atom%hop(1:2) = -hresv; atom%hopz = 0
     end select
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
end subroutine YRZ