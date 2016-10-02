subroutine adjustFermilevel( )
  use vmcplace
  implicit none

  integer spin, ik, i, ia, ja
  real*8, dimension ( nd,  nuc) :: Ru, kv
  real*8 ek( na1uc ), eval( naa ), Ef(2)
  real*8, dimension ( nd ) :: vk
  complex*16, dimension (na1uc,  na1uc) :: hk
  real*8, external :: fermilevel

  !purpose: adjust muloc and mzloc so that in the free system the
  !         filling matches numer of electrons per spin, n1spin = Npair.
  !         this is useful to get a good guess of muloc and mzloc.

  call getMesh(Ru, kv)

  do spin = 1, 2    
     do ik = 1, nuc;  vk = kv(:, ik)
        call hamilton(na1uc, hk, nd, vk, unitcell, spin)	 
	    call ZHEIGEN(na1uc, hk, ek)
		eval( (ik-1)*na1uc+1 : ik*na1uc ) = ek
	 end do
	 Ef(spin) = fermilevel( naa, eval, npair )	 
  end do

  unitcell%atom%muloc = unitcell%atom%muloc + sum(Ef) / 2
  unitcell%atom%mzloc = unitcell%atom%mzloc + ( Ef(1) - Ef(2) ) / 2

  !print*, 'muloc and mzloc adjusted as follows:'
  !print*, unitcell%atom%muloc, unitcell%atom%mzloc

  return
end subroutine adjustFermilevel


function element( Ri, Rj )
  use vmcplace
  implicit none
  integer, dimension ( ndim) :: Ri, Rj
  complex*16 element

  integer a, b, ir2r, sign, id
  integer r2r( nd ), rt(  nd )
  integer, external :: typewriter

  a = Ri(ndim);  b = Rj(ndim)
  r2r = Ri(1:nd) - Rj(1:nd)
  rt = mod( r2r + limits(1:nd), limits(1:nd) )
  
  sign = 1
  do id = 1, nd; if(  antiperiodic(id) .and. r2r(id) /= rt(id) )sign = -sign; end do

  rt = rt + 1;  ir2r = typewriter(nd, limits(1:nd), rt)

  if( whichtable == 0)then
    element =  table(a, b, ir2r) * sign
  else
    element =  tables(a, b, ir2r, whichtable) * sign
  end if

  return
end function element

subroutine newvector(spin, Rf, vec)
  use vmcplace
  implicit none

  integer spin
  integer Rf( ndim)
  complex*16 vec( Npair)

  integer ia, i, ja, j
  complex*16 element

  if(spin == 1)then
    do ja = 1, naa; j =  occup(ja, 2); if(j == 0)cycle
	   vec(j) = element( Rf,  A2R(:, ja) )
	end do
  else
    do ia = 1, naa; i =  occup(ia, 1); if(i == 0)cycle
	   vec(i) = element(  A2R(:, ia), Rf )
	end do
  end if

  return
end subroutine newvector


subroutine tabulate( )
  use vmcplace
  implicit none
  
  integer idim, ik, i, ia, ja
  real*8, dimension (nd,  nuc) :: Ru, kv
  real*8 tmax, ek(nambu)
  complex*16 one
  real*8, dimension (nd) :: vk
  complex*16 ak(nambu, nambu)
  complex*16, dimension ( na1uc,  na1uc) :: hk, pk, work, VL, VR, zek(na1uc)
  logical icheck/.true./

  one = dcmplx(0.d0, 1.d0)

  call getMesh(Ru, kv)

  if(.not. allocated(table) ) allocate(  table(na1uc, na1uc, nuc) )
  table = 0

  do ik = 1, nuc;  vk = kv(:, ik)
     call hamilton(na1uc, hk, nd, vk, unitcell, 1)   !spin up h(k)
	 ak(1:na1uc, 1:na1uc) = hk
     
     if(icheck)then
       call ZHEIGEN( na1uc, hk, ek(1:na1uc) )
       if( count( abs(ek(1:na1uc)) < 1.d-10 ) > 0 )stop 'open shell found in hk_up @ tabulate'
    end if

	 call pairing(na1uc, pk, nd, vk, unitcell)
	 ak(1:na1uc, na1uc+1:nambu) = pk
	 ak(na1uc+1:nambu, 1:na1uc) = transpose( conjg(pk) )
	 	 
	 call hamilton(na1uc, hk, nd, -vk, unitcell, 2)  !spin down h(-k)
	 ak(na1uc+1:nambu, na1uc+1:nambu) = - conjg(hk)

     if(icheck)then
       call ZHEIGEN( na1uc, hk, ek(1:na1uc) )
       if( count( abs(ek(1:na1uc)) < 1.d-10 ) > 0 ) stop 'open shell found in hk_down @ tabulate'
	 end if 

	 call ZHEIGEN(nambu, ak, ek)
     if( count( ek < 0 ) /= na1uc )stop 'p-h symmetry fails in bdg hamiltonian @ tabulate'
	 
	 hk = ak(1:na1uc, 1:na1uc)         !acts as Uk = {uk} for negative-energy BdG eigenstates
	 pk = ak(na1uc+1:nambu, 1:na1uc)   !acts as Vk = {vk} for negative-energy BdG eigenstates

     if(icheck)then	 
       work = pk; call ZGEIGEN(na1uc, work, VL, VR, zek)
       if( count( abs(zek) < 1.d-10 ) > 0 )stop 'singular kernel matrix found @ tabulate'
     end if

	 call ZINVERT(na1uc, pk)          !inverse of Vk
     pk = matmul(hk, pk)              !pk = Uk * inv(Vk)
	 
	 call gaugephase()

	 do i = 1, nuc
	    table(:, :, i) = table(:, :, i) + pk * exp( one*sum(vk*Ru(:, i)) )
	 end do
  end do

  table = table / nuc 

  !tmax = 0
  !do i = 1, nuc; do ia = 1, na1uc; do ja = 1, na1uc
  !   tmax = max( tmax, abs(table(ia, ja, i)) )
  !end do; end do; end do

  !table = table / tmax
  
  return
  contains
    subroutine gaugephase()
	  integer ia, ja
	  real*8 rv(nd)
	  do ja = 1, na1uc; do ia = 1, na1uc
	     rv =  unitcell%atom(ia)%r(1:nd) -  unitcell%atom(ja)%r(1:nd)
		 pk(ia, ja) = pk(ia, ja) * exp( one * sum(vk*rv) )
	  end do; end do
	  return
	end subroutine gaugephase	
end subroutine tabulate

subroutine getMesh(Ru, kv)
  use vmcplace
  implicit none

  real*8 Ru( nd,  nuc), kv( nd,  nuc)
  
  integer i, id
  real*8 twopi, offset
  integer Rint( nd)

  twopi = asin(1.d0) * 4

  do i = 1, nuc
     call typewriter2coord(i, nd, limits(1:nd), Rint)    
	 Ru(:, i) = 0
	 do id = 1, nd
	    Ru(:, i) = Ru(:, i) + (Rint(id) - 1) * unitcell%abc(1:nd, id)
	 end do
  end do

  do i = 1, nuc
     call typewriter2coord(i, nd, limits(1:nd), Rint)    
	 kv(:, i) = 0
	 do id = 1, nd; offset = 0; if(  antiperiodic(id) ) offset = 0.5d0
	    kv(:, i) = kv(:, i) + (Rint(id) - 1 + offset) * unitcell%dualabc(1:nd, id) / limits(id)
	 end do
  end do

  kv = kv * twopi

  return
end subroutine getMesh





