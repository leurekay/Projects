subroutine kagomeJ( )
  use vmcplace
  implicit none

  integer i, j, k, ib, nb, ia, a
  
  integer nx, ny
  real*8  Vnn, Jnn
  real*8  bond(2), shift(3)
  character*6 vtable(97)
  real*8, target, dimension(97) :: v
  real*8, pointer :: dr(:), di(:), hr(:), hi(:), mu
  character*20 explain
  logical iEf/.true./

  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom, btom


  !============================================================================
  !                    constants and mc controllers
  model = 'kagomeJ'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'control.'//model)
     read(10, *)explain;  read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain;  read(10, *)bruteforce;  close(10)
  end do

  !==============================================================================
  !                lattice, npair and physical parameters

  nd = 2;  ndim = nd + 1;  na1uc = 12;  nambu = na1uc * 2
  !Npair = (naa * 5) / 12
  !U = 0;  mott = 0;   Vnn = 2;  Jnn = 0

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'lattice.'//model)
     read(10, *)explain;  read(10, *)nx, ny, npair
     read(10, *)explain;  read(10, *)antiperiodic
     read(10, *)explain;  read(10, *)U, Vnn, Jnn
     read(10, *)explain;  read(10, *)mott;   close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2);  naa = nuc * na1uc

  !===============================================================================
  !                variational vars 

  vtable(1:24) = (/'dr11', 'dr12', 'dr21', 'dr22', 'dr31', 'dr32', 'dr41', 'dr42', 'dr51', 'dr52', 'dr61', 'dr62', 'dr71', 'dr72', 'dr81', 'dr82', 'dr91', 'dr92', 'dr101', 'dr102', 'dr111', 'dr112', 'dr121', 'dr122'/)
  vtable(25:48) = (/'di11', 'di12', 'di21', 'di22', 'di31', 'di32', 'di41', 'di42', 'di51', 'di52', 'di61', 'di62', 'di71', 'di72', 'di81', 'di82', 'di91', 'di92', 'di101', 'di102', 'di111', 'di112', 'di121', 'di122'/)
  vtable(49:72) = (/'hr11', 'hr12', 'hr21', 'hr22', 'hr31', 'hr32', 'hr41', 'hr42', 'hr51', 'hr52', 'hr61', 'hr62', 'hr71', 'hr72', 'hr81', 'hr82', 'hr91', 'hr92', 'hr101', 'hr102', 'hr111', 'hr112', 'hr121', 'hr122'/)
  vtable(73:96) = (/'hi11', 'hi12', 'hi21', 'hi22', 'hi31', 'hi32', 'hi41', 'hi42', 'hi51', 'hi52', 'hi61', 'hi62', 'hi71', 'hi72', 'hi81', 'hi82', 'hi91', 'hi92', 'hi101', 'hi102', 'hi111', 'hi112', 'hi121', 'hi122'/)
  vtable(97) = 'mu'
  dr => v(1:24);  di => v(25:48);  hr => v(49:72);  hi => v(73:96); mu => v(97)
  !import var bounds and determine vars from file (and external optimizer)

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     call register(96, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do

  if(optimize .and. optinfo == 0) return


  !==========================================================================
  !       automated inputs

  unitcell%natom = na1uc

  abc => unitcell%abc;  dualabc => unitcell%dualabc
  abc(:, 1) = (/1, 0, 0/)
  abc(:, 2) = (/1.d0, sqrt(3.d0), 0.d0/)/2
  abc(:, 3) = (/0, 0, 1/)
  call dualvector(abc(:, 2), abc(:, 3), abc(:, 1), dualabc(:, 1))  
  call dualvector(abc(:, 3), abc(:, 1), abc(:, 2), dualabc(:, 2))  
  call dualvector(abc(:, 1), abc(:, 2), abc(:, 3), dualabc(:, 3))  

  if(.not. optimize .or. optinfo == 1) allocate( unitcell%atom(na1uc) )
  unitcell%atom(1)%r = 0
  unitcell%atom(2)%r = (/1.d0, sqrt(3.d0), 0.d0/)/8
  unitcell%atom(3)%r = (/-1.d0, sqrt(3.d0), 0.d0/)/8
  unitcell%atom%kondo = 0

  nb = 2
  do ia = 1, 3;  atom => unitcell%atom(ia)

     atom%nbond = nb;  atom%E = 0;       atom%B = 0
     atom%muloc = mu;  atom%mzloc = 0;    atom%pairloc = 0
	    
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0

     atom%t = 1;      atom%tz = 0;    atom%J = Jnn;   atom%V = 0 
	 atom%hop = 0;    atom%hopz = 0;  atom%pair = 0;  atom%pairz = 0

     select case (ia) 
	   case (1);  bond = abc(1:2, 2) / 4
	   case (2);  bond = abc(1:2, 1) / 4
	   case (3);  bond = abc(1:2, 1) / 4
	 end select
	 atom%bond(1:2, 1) = bond
	 call rotate60(bond);  if(ia == 3)call rotate60(bond)
	 atom%bond(1:2, 2) = bond     
  
  end do

  !copy info in the first primitive cell to the other three primitive cells in the same unitcell 
  !direct assignment btom = atom does not seem to work correctly, so this is done de forte
  do j = 1, 2;  do i = 1, 2
     k = i + (j-1) * 2;  if(k == 1)cycle
	 shift = (i-1) * abc(:, 1) / 2 + (j-1) * abc(:, 2) / 2
	 do a = 1, 3;  atom => unitcell%atom(a)
	    ia = a + (k-1) * 3;  btom => unitcell%atom(ia)

        if(.not. optimize .or. optinfo == 1) allocate( btom%t(nb), btom%tz(nb), btom%J(nb), btom%V(nb) )
        if(.not. optimize .or. optinfo == 1) allocate( btom%hop(nb), btom%hopz(nb), btom%pair(nb), btom%pairz(nb) )
        if(.not. optimize .or. optinfo == 1) allocate( btom%bond(3, nb) )

	    btom%r = atom%r + shift;  btom%nbond = nb;  btom%E = 0;  btom%B = 0
		btom%muloc = atom%muloc;  btom%mzloc = atom%mzloc;  btom%pairloc = atom%pairloc
		
		btom%bond = atom%bond;  btom%t = atom%t;  btom%tz = atom%tz;  btom%J = atom%J;  btom%V = atom%V
		btom%hop = atom%hop; btom%hopz = atom%hopz;  btom%pair = atom%pair;  btom%pairz = atom%pairz
	 end do
  end do;  end do

  do ia = 1, 12;  atom => unitcell%atom(ia);  do ib = 1, 2;  j = ib + (ia-1) * 2
     atom%pair(ib) = cmplx( dr(j), di(j) );  atom%pairz(ib) = cmplx(hr(j), hi(j))  !atom%hop(ib) = cmplx(hr(j), hi(j))
  end do;  end do


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

end subroutine kagomeJ


subroutine kagome( )
  use vmcplace
  implicit none

  integer i, nb, ia
  complex*16 one
  
  integer nx, ny
  real*8  Vnn, Jnn
  real*8  pi, bond(2)
  real*8, pointer :: mu, cdw, did, fbond, sloc, sbond, gU_, Vdh_, Vcc_, Vzz_
  character*6 vtable(20)
  real*8, target, dimension(20) :: v
  character*20 explain
  logical iEf/.true./

  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom


  !============================================================================
  !                    constants and mc controllers
  one = dcmplx(0.d0, 1.d0);  pi = asin(1.d0) * 2;  model = 'kagome'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle
     open(10, file = 'control.'//model)
     read(10, *)explain;  read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain;  read(10, *)bruteforce;  close(10)
  end do

  !==============================================================================
  !                lattice, npair and physical parameters

  nd = 2;  ndim = nd + 1;  na1uc = 3;  nambu = na1uc * 2
  !Npair = (naa * 5) / 12
  !U = 0;  mott = 0;   Vnn = 2;  Jnn = 0

  do i = 0, nthread - 1
    call MPI_barrier( mpi_comm_world, ierr )  
    if(ithread /= i)cycle
    open(10, file = 'lattice.'//model)
    read(10, *)explain;  read(10, *)nx, ny, npair
    read(10, *)explain;  read(10, *)antiperiodic
    read(10, *)explain;  read(10, *)U, Vnn, Jnn
    read(10, *)explain;  read(10, *)mott;   close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2);  naa = nuc * na1uc

  !===============================================================================
  !                variational vars 

  vtable(1:10) = (/'mu', 'cdw', 'sloc', 'sbond', 'did', 'fbond', 'gU', 'Vdh', 'Vzz', 'Vcc'/)
  mu => v(1);   cdw => v(2);  sloc => v(3);   sbond => v(4);  did => v(5);  fbond => v(6)
  gU_ => v(7);  Vdh_ => v(8);  Vzz_ => v(9);  Vcc_ => v(10)
  
  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(10, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do  
  if(optimize .and. optinfo == 0) return

  
  !==========================================================================
  !       automated inputs

  gU = gU_;  Vdh = Vdh_;  Vzz = Vzz_;  Vcc = Vcc_

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
  unitcell%atom(2)%r = (/1., sqrt(3.), 0./)/4
  unitcell%atom(3)%r = (/-1., sqrt(3.), 0./)/4
  unitcell%atom%kondo = 0

  nb = 2
  do ia = 1, 3;  atom => unitcell%atom(ia)

     atom%nbond = nb;  atom%E = 0;       atom%B = 0
     atom%muloc = mu;  atom%mzloc = 0;  atom%pairloc = sloc
	    
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0

     atom%t=1;      atom%tz = 0;    atom%J = Jnn;   atom%V = Vnn 
	 atom%hop = 1;  atom%hopz = 0;  atom%pair = 0;  atom%pairz = 0

     select case (ia) 
	   case (1);  bond = abc(1:2, 2) / 2;  atom%muloc = atom%muloc + cdw
	              atom%pair = sbond + did * exp( one * (/1, 3/) * 2 * pi / 3) + fbond * (/1, 1/)
	   case (2);  bond = abc(1:2, 1) / 2;  atom%muloc = atom%muloc - cdw * 0.5
	              atom%pair = sbond + did * exp( one * (/2, 4/) * 2 * pi / 3) + fbond * (/-1, -1/)
	   case (3);  bond = abc(1:2, 1) / 2;  atom%muloc = atom%muloc - cdw * 0.5
	              atom%pair = sbond + did * exp( one * (/5, 6/) * 2 * pi / 3) + fbond * (/1, -1/)
	 end select
	 atom%bond(1:2, 1) = bond
	 call rotate60(bond);  if(ia == 3)call rotate60(bond)
	 atom%bond(1:2, 2) = bond     
  
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

end subroutine kagome


subroutine kagome2x2( )
  use vmcplace
  implicit none

  integer i, j, k, nb, ia, a
  complex*16 one
  
  integer nx, ny
  real*8  Vnn, Jnn
  real*8  pi, bond(2), shift(3), hop(2, 12)
  real*8, pointer :: mu, cdw, did, david, sloc, sbond, gU_, Vdh_, Vcc_, Vzz_
  character*6 vtable(20)
  real*8, target, dimension(20) :: v
  character*20 explain
  logical iEf/.true./

  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom, btom


  !============================================================================
  !                    constants and mc controllers
  one = dcmplx(0.d0, 1.d0);  pi = asin(1.d0) * 2;   model = 'kagome2x2'

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     open(10, file = 'control.'//model)
     read(10, *)explain;  read(10, *) nbath, nsample, nmc1sample
     read(10, *)explain;  read(10, *)bruteforce;  close(10)
  end do

  !==============================================================================
  !                lattice, npair and physical parameters

  nd = 2;  ndim = nd + 1;  na1uc = 12;  nambu = na1uc * 2
  !Npair = (naa * 5) / 12
  !U = 0;  mott = 0;   Vnn = 2;  Jnn = 0

  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle 
     open(10, file = 'lattice.'//model)
     read(10, *)explain;  read(10, *)nx, ny, npair
     read(10, *)explain;  read(10, *)antiperiodic
     read(10, *)explain;  read(10, *)U, Vnn, Jnn
     read(10, *)explain;  read(10, *)mott;   close(10)
  end do

  limits = (/nx, ny, na1uc/);  nuc = limits(1) * limits(2);  naa = nuc * na1uc

  !===============================================================================
  !                variational vars 

  vtable(1:10) = (/'mu', 'cdw', 'david', 'sloc', 'sbond', 'did', 'gU', 'Vdh', 'Vzz', 'Vcc'/)
  mu => v(1);  cdw => v(2);  david => v(3);  sloc => v(4);    sbond => v(5);  did => v(6)
  gU_ => v(7);    Vdh_ => v(8);    Vzz_ => v(9);   Vcc_ => v(10)

  !import var bounds and determine vars from file (and external optimizer)
  do i = 0, nthread - 1
     call MPI_barrier( mpi_comm_world, ierr )  
     if(ithread /= i)cycle  
     call register(10, vtable, v, nvar, vname, vmin, vmax, var, model, optimize)
  end do    
  if(optimize .and. optinfo == 0) return

  
  !==========================================================================
  !       automated inputs

  gU = gU_;  Vdh = Vdh_;  Vzz = Vzz_;  Vcc = Vcc_

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
  unitcell%atom(2)%r = (/1., sqrt(3.), 0./)/8
  unitcell%atom(3)%r = (/-1., sqrt(3.), 0./)/8
  unitcell%atom%kondo = 0

  nb = 2
  do ia = 1, 3;  atom => unitcell%atom(ia)

     atom%nbond = nb;  atom%E = 0;       atom%B = 0
     atom%muloc = mu;  atom%mzloc = 0;   atom%pairloc = sloc * one
	    
     if(.not. optimize .or. optinfo == 1) allocate( atom%t(nb), atom%tz(nb), atom%J(nb), atom%V(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%hop(nb), atom%hopz(nb), atom%pair(nb), atom%pairz(nb) )
     if(.not. optimize .or. optinfo == 1) allocate( atom%bond(3, nb) );  atom%bond = 0

     atom%t = 1;      atom%tz = 0;    atom%J = Jnn;   atom%V = Vnn 
	 atom%hop = 1;    atom%hopz = 0;  atom%pair = 0;  atom%pairz = 0

     select case (ia) 
	   case (1);  bond = abc(1:2, 2) / 4;  atom%muloc = atom%muloc + cdw
	              atom%pair = did * exp( one * (/1, 3/) * 2 * pi / 3)
	   case (2);  bond = abc(1:2, 1) / 4;  atom%muloc = atom%muloc - cdw * 0.5
	              atom%pair = did * exp( one * (/2, 4/) * 2 * pi / 3)
	   case (3);  bond = abc(1:2, 1) / 4;  atom%muloc = atom%muloc - cdw * 0.5
	              atom%pair = did * exp( one * (/5, 6/) * 2 * pi / 3)
	 end select
	 atom%bond(1:2, 1) = bond
	 call rotate60(bond);  if(ia == 3)call rotate60(bond)
	 atom%bond(1:2, 2) = bond     
  
  end do

  !copy info in the first primitive cell to the other three primitive cells in the same unitcell 
  !direct assignment btom = atom does not seem to work correctly, so this is done de forte
  do j = 1, 2;  do i = 1, 2
     k = i + (j-1) * 2;  if(k == 1)cycle
	 shift = (i-1) * abc(:, 1) / 2 + (j-1) * abc(:, 2) / 2
	 do a = 1, 3;  atom => unitcell%atom(a)
	    ia = a + (k-1) * 3;  btom => unitcell%atom(ia)

        if(.not. optimize .or. optinfo == 1) allocate( btom%t(nb), btom%tz(nb), btom%J(nb), btom%V(nb) )
        if(.not. optimize .or. optinfo == 1) allocate( btom%hop(nb), btom%hopz(nb), btom%pair(nb), btom%pairz(nb) )
        if(.not. optimize .or. optinfo == 1) allocate( btom%bond(3, nb) )

	    btom%r = atom%r + shift;  btom%nbond = nb;  btom%E = 0;  btom%B = 0
		btom%muloc = atom%muloc;  btom%mzloc = atom%mzloc;  btom%pairloc = atom%pairloc
		
		btom%bond = atom%bond;  btom%t = atom%t;  btom%tz = atom%tz;  btom%J = atom%J;  btom%V = atom%V
		btom%hop = atom%hop; btom%hopz = atom%hopz;  btom%pair = atom%pair;  btom%pairz = atom%pairz
	 end do
  end do;  end do

  hop(1, :) = 0.5 - (/0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1/) 
  hop(2, :) = 0.5 - (/1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0/) 

  do ia = 1, 12;  atom => unitcell%atom(ia)
     atom%hop = 1 + hop(:, ia) * david
     atom%pair = atom%pair * hop(:, ia) * 2                  !multiply staggered sign for d+id assigned previously         
	 atom%pair = atom%pair + hop(:, ia) * sbond              !this would assign a david-star pairing using sbond as the amplitude
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

end subroutine kagome2x2






subroutine optimized_var( )
  use vmcplace
  implicit none

  integer i
  character*6 vtable(48)
  real*8 mu
  real*8, target, dimension(48) :: v
  real*8, pointer :: dr(:), di(:), hr(:), hi(:)
  character*20 explain

  vtable(1:24) = (/'dr11', 'dr12', 'dr21', 'dr22', 'dr31', 'dr32', 'dr41', 'dr42', 'dr51', 'dr52', 'dr61', 'dr62', 'dr71', 'dr72', 'dr81', 'dr82', 'dr91', 'dr92', 'dr101', 'dr102', 'dr111', 'dr112', 'dr121', 'dr122'/)
  vtable(25:48) = (/'di11', 'di12', 'di21', 'di22', 'di31', 'di32', 'di41', 'di42', 'di51', 'di52', 'di61', 'di62', 'di71', 'di72', 'di81', 'di82', 'di91', 'di92', 'di101', 'di102', 'di111', 'di112', 'di121', 'di122'/)
  vtable(49:72) = (/'hr11', 'hr12', 'hr21', 'hr22', 'hr31', 'hr32', 'hr41', 'hr42', 'hr51', 'hr52', 'hr61', 'hr62', 'hr71', 'hr72', 'hr81', 'hr82', 'hr91', 'hr92', 'hr101', 'hr102', 'hr111', 'hr112', 'hr121', 'hr122'/)
  vtable(73:96) = (/'hi11', 'hi12', 'hi21', 'hi22', 'hi31', 'hi32', 'hi41', 'hi42', 'hi51', 'hi52', 'hi61', 'hi62', 'hi71', 'hi72', 'hi81', 'hi82', 'hi91', 'hi92', 'hi101', 'hi102', 'hi111', 'hi112', 'hi121', 'hi122'/)

  dr => v(1:24);  di => v(25:48);  hr => v(49:72);  hi => v(73:96)
  !import var bounds and determine vars from file (and external optimizer)

  open(10, file = 'optimized.kagomeJ')
  read(10, *)dr, di, hr, hi, mu;  close(10)

  open(10, file = 'var.kagomeJ')
  do i = 1, 24
     write(10, 100)vtable(i), dr(i), dr(i)
  end do

  do i = 1, 24
     write(10, 100)vtable(i+24), di(i), di(i)
  end do

  do i = 1, 24
     write(10, 100)vtable(i+48), hr(i), hr(i)
  end do

  do i = 1, 24
     write(10, 100)vtable(i+72), hi(i), hi(i)
  end do

  write(10, 100)'mu', mu, mu
  close(10)


100 format(1x, a6, 2f20.6)
  return
end subroutine optimized_var