subroutine normal_band(unitcell)
  use crystal

  implicit none
  type (unitcellconfig), target :: unitcell

  integer spin, nd, ndim 
  real*8, dimension(2) :: G, X, M, dk
  logical, external :: active
  real*8, parameter :: pi = 3.14159265358979d0

  nd = 2;  ndim = unitcell%natom
  
  G = 0;  M = (/1, 1/);  X = (/1, 0/)
  if( active( sum( unitcell%abc(:, 1) * unitcell%abc(:, 2) ) ) )then
    M=(/2.d0/3, 2.d0/dsqrt(3.d0)/); X = (/1.d0, 1.d0/dsqrt(3.d0)/)
  end if
  M = M * pi;  X = X * pi

  open(10, file = 'normalband.dat')

  do spin = 1, 2
     dk = (M - G) / 100;  call linecut(G, dk, 99)
     dk = (X-M) / 100;  call linecut(M, dk, 99)
     dk = (G-X) / 100;  call linecut(X, dk, 100)
  end do;  close(10)

  return
  contains
    subroutine linecut(k0, dk, n)
      implicit none
      integer i
      complex*16 hk(ndim, ndim)
      real*8 ek(ndim)
      real*8 k0(2), dk(2), kv(2)
      integer n
      do i = 0, n;  kv = k0 + dk * i
         call hamilton(ndim, hk, nd, kv, unitcell, spin)
         call ZHEIGEN(ndim, hk, ek);  write(10, 100)ek
      end do
100   format(1x, f15.6)
      return
    end subroutine linecut
end subroutine normal_band



subroutine bdg_band(unitcell)
  use crystal
  implicit none
  type (unitcellconfig), target :: unitcell

  integer nd, ndim, nambu
  real*8, dimension(2) :: G, X, M, dk
  logical, external :: active
  real*8, parameter :: pi = 3.14159265358979d0

  nd = 2;  ndim = unitcell%natom;  nambu = ndim * 2
  
  G = 0;  M = (/1, 1/);  X = (/1, 0/)
  if( active( sum( unitcell%abc(:, 1) * unitcell%abc(:, 2) ) ) )then
    M=(/2.d0/3, 2.d0/dsqrt(3.d0)/); X = (/1.d0, 1.d0/dsqrt(3.d0)/)
  end if
  M = M * pi;  X = X * pi

  open(10, file = 'bdgband.dat')

  dk = (M-G) / 100;  call linecut(G, dk, 99)
  dk = (X-M) / 100;  call linecut(M, dk, 99)
  dk = (G-X) / 100;  call linecut(X, dk, 100)
  
  close(10)

  return
  contains
    subroutine linecut(k0, dk, n)
      implicit none
      integer i
      complex*16 hk(ndim, ndim), gk(ndim, ndim), ak(nambu, nambu)
      real*8 ek(nambu)
      real*8 k0(2), dk(2), kv(2)
      integer n
      do i = 0, n;  kv = k0 + dk * i
         call hamilton(ndim, hk, nd, kv, unitcell, 1);  ak(1:ndim, 1:ndim) = hk
         call hamilton(ndim, hk, nd, -kv, unitcell, 2); ak(ndim+1:nambu, ndim+1:nambu) = - conjg(hk)
         call pairing(ndim, gk, nd, kv, unitcell);  ak(1:ndim, ndim+1:nambu) = gk;  ak(ndim+1:nambu, 1:ndim) = transpose( conjg(gk) )
         call ZHEIGEN(nambu, ak, ek);  write(10, 100)ek
      end do
100   format(1x, f15.6)
      return
    end subroutine linecut
end subroutine bdg_band

  

subroutine hamilton(n, hk, nd, kv, unitcell, spin)
  use crystal
  implicit none
  integer n, nd, spin
  real*8 kv(nd)
  complex*16 hk(n, n)
  type (unitcellconfig), target :: unitcell

  complex*16 one, phase
  integer ia, ja, ib
  type (atomconfig), pointer :: atom, btom

  if(n /= unitcell%natom) stop 'error in dim of hk'

  one = dcmplx(0.d0, 1.d0);  hk = 0

  do ia = 1, n;  atom => unitcell%atom(ia)   
	 do ib = 1, atom%nbond; ja = atom%nbatom(ib)
		phase = exp( one * sum( kv * atom%bond(1:nd, ib) ) )
		hk(ia, ja) = hk(ia, ja) - ( atom%hop(ib) + atom%hopz(ib) * (3-2*spin) ) * phase
     end do
  end do
  hk = hk + transpose( conjg(hk) )

  do ia = 1, n;  atom => unitcell%atom(ia)
	 hk(ia, ia) = hk(ia, ia) - atom%muloc - atom%mzloc * (3-2*spin)
  end do

  return
end subroutine hamilton


subroutine pairing(n, gk, nd, kv, unitcell)
  use crystal
  implicit none
  integer n, nd
  real*8 kv(nd)
  complex*16 gk(n, n)
  type (unitcellconfig), target :: unitcell

  complex*16 one, phase
  integer ia, ib, ja
  type (atomconfig), pointer :: atom, btom

  if(n /= unitcell%natom) stop 'error in dim of gk'

  one = dcmplx(0.d0, 1.d0);  gk = 0

  do ia = 1, n;  atom => unitcell%atom(ia)     
	 do ib = 1, atom%nbond; ja = atom%nbatom(ib)
		phase = exp( one * sum( kv * atom%bond(1:nd, ib) ) )
		gk(ia, ja) = gk(ia, ja) - ( atom%pair(ib) + atom%pairz(ib) ) * phase
		gk(ja, ia) = gk(ja, ia) - ( atom%pair(ib) - atom%pairz(ib) ) * conjg(phase)     
	 end do
  end do

  do ia = 1, n;  atom => unitcell%atom(ia)
	 gk(ia, ia) = gk(ia, ia) - atom%pairloc
  end do

  return
end subroutine pairing









subroutine bandstructure(nd, na, GM, GX, unitcell)
  use crystal
  implicit none
  integer nd, na
  type (unitcellconfig), target :: unitcell
  real*8, dimension(nd) :: GX, GM

  integer nk, nambu, nx, ny
  integer idim, ik, i
  real*8 pi, kx, ky, eta, w, ek(2*na)
  complex*16 one
  real*8, dimension (nd) :: vk, dk, kinit, kfinal
  complex*16, dimension(2*na, 2*na) :: ak, gk
  complex*16, dimension (na, na) :: hk, pk
  complex*16 bra(2*na, 1), ket(1, 2*na)

  one = dcmplx(0.d0, 1.d0); pi = asin(1.d0) * 2; eta = 0.05
  nambu = na*2

  open(10, file = 'band.dat');  open(11, file='aup.dat');  open(12, file='adn.dat')
  kinit = 0;   kfinal = GM;  call linecut(0)
  kinit = GM;  kfinal = GX;  call linecut(0)
  kinit = GX;  kfinal = 0;   call linecut(1)  
  close(10); close(11);  close(12)

  return

  contains
  subroutine linecut(info)
  integer info, m
  dk = (kfinal - kinit) / 100
  do i = 1, 100+info;  vk = ( kinit + dk * (i-1) ) * pi
     call hamilton(na, hk, nd, vk, unitcell, 1)   !spin up h(k)
	 ak(1:na, 1:na) = hk
	 
	 call pairing(na, pk, nd, vk, unitcell)
	 ak(1:na, na+1:nambu) = pk
	 ak(na+1:nambu, 1:na) = transpose( conjg(pk) )
	 	 
	 call hamilton(na, hk, nd, -vk, unitcell, 2)  !spin down h(-k)
	 ak(na+1:nambu, na+1:nambu) = - conjg(hk)
	 
	 call ZHEIGEN(nambu, ak, ek);	 write(10, 100)ek
	 do m = 1, nambu
	    write(11, *)abs( sum( ak(1:na, m) ) )**2
		write(12, *)abs( sum(ak(na+1:nambu, m) ) )**2
	 end do
  end do
100 format(1x, e15.6)
  return
  end subroutine linecut
end subroutine bandstructure