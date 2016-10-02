subroutine InfraStructure(  )
  use vmcplace
  implicit none
  
  integer i
  real*8, parameter :: gU_max = 0.9999d0

  if(optimize .and. optinfo > 1)then
    call tabulate( );  return       
	!the rest infras already established
  end if

  !overwrite user's spec if any
  ndim = nd + 1;   na1uc = unitcell%natom;   nambu = na1uc * 2;    limits(ndim) = na1uc
  nuc = 1;  do i = 1, nd;  nuc = nuc * limits(i);  end do;         naa = nuc * na1uc

  if(npair <= 0)stop 'npair <= 0 is invalid'
  if(npair * 2 > naa)stop 'npair * 2 > naa is invalid'
  
  !mapping from atoms on lattice to coordinates
  if(.not. allocated( a2R ) ) allocate( a2R(ndim, naa) )
  if(.not. allocated( kondo ) ) allocate( kondo(naa) )

  kondolattice = .false.;   kondosites = 0
  do i = 1, naa
     call typewriter2coord(i, ndim, limits, a2R(:, i))
	 kondo(i) =  unitcell%atom( a2R(ndim, i) )%kondo
	 if( kondo(i) )then
	    kondolattice = .true.
	    kondosites =  kondosites + 1
	 end if
  end do

  !if( kondolattice .and. mott)mott = .false.

  if(mott)then;  gU = 1;  else;  gU = min(gU, gU_max); end if

  call tabulate(  )

  !nearest neighbors of all atoms on lattice (used in penalty)
  call includenn(  )

  !all independent neighbors (used in updateConfig and measure)
  call coordination(  )
  return
end subroutine InfraStructure





