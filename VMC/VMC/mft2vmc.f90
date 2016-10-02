subroutine mft2vmc( )
  use vmcplace
  implicit none

  integer a, nb, ib
  complex*16 one
  real*8 pairloc, mu, mz, Vnn, bond(2)
  real*8, pointer, dimension(:, :) :: abc, dualabc
  type (atomconfig), pointer :: atom
  character*20 sample, explain

  open(10, file='vmcinput.txt')

  read(10, *)explain
  read(10, *)sample
  
  read(10, *)explain
  read(10, *)nd;  ndim = nd + 1

  read(10, *)explain
  read(10, *)na1uc, nambu;  unitcell%natom = na1uc

  read(10, *)explain
  read(10, *)unitcell%abc

  allocate( unitcell%atom(na1uc) )

  do a = 1, na1uc;  atom => unitcell%atom(a)
     read(10, *)explain
     read(10, *)atom%r

     read(10, *)explain
     read(10, *)nb;  atom%nbond = nb

     allocate(atom%bond(3, nb));  atom%bond = 0

     read(10, *)explain
	 read(10, *)atom%bond(1:nd, :)

     allocate(atom%t(nb), atom%J(nb), atom%V(nb) )

     read(10, *)explain
     read(10, *)atom%t, atom%J;  atom%V = 0

     read(10, *)explain
     read(10, *)atom%E

     read(10, *)explain
     read(10, *)atom%B
 
     read(10, *)explain
     read(10, *)atom%muloc

     read(10, *)explain
	 read(10, *)atom%mzloc

     allocate( atom%hop(nb), atom%pair(nb), atom%pairz(nb) )
     read(10, *)explain
  
     do ib = 1, nb
	    read(10, *)atom%hop(ib), atom%pair(ib), atom%pairz(ib)
	 end do
  end do
  close(10)


  return
end subroutine mft2vmc