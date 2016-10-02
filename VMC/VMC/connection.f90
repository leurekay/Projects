subroutine connection(Ri, ia, bond, Rj, ja, info)
  use vmcplace
  implicit none

  integer ia, ja, info
  integer, dimension (nd) :: Ri, Rj    !integer coord in unit of unitcells
  real*8 bond(3)

  real*8 Rv(3), Rt(3), vecs(3, 3)
  integer id, Rjj(nd)


  vecs = unitcell%abc; Rv = 0
  do id = 1, nd; Rv = Rv + Ri(id) * vecs(:, id); end do              
  !absolute position of unitcell Ri

  Rv = Rv + unitcell%atom(ia)%r + bond - unitcell%atom(ja)%r        
  !absolute position of a candidate unitcell holding ja


  vecs = unitcell%dualabc
  do id = 1, nd; Rj(id) = nint( sum(Rv * vecs(:, id)) ); end do   
  !translate Rv into integer coordinate Rj in terms of a, b, c

  vecs = unitcell%abc;  Rt = 0
  do id = 1, nd; Rt = Rt + Rj(id) * vecs(:, id); end do           
  !vector corresponds to unitcell Rj

  info = 0
  if( sum( abs(Rt-Rv) ) > 1.d-5 )return                         
  !ja is not connected to (ia, Ri) by the bond
  
  Rjj = mod( Rj-1 + 2 * limits(1:nd), limits(1:nd) ) + 1        
  !truncate Rj into positive labels
  
  info = 1
  do id = 1, nd
     if( Rjj(id) /= Rj(id) .and. antiperiodic(id) ) info = - info
  end do

  Rj = Rjj

  return
end subroutine connection

subroutine coordination()
  use vmcplace
  implicit none
  integer ia, a, ib, b, info
  integer, dimension (ndim) :: Ri, Rj
  type (atomconfig), pointer :: atom
  integer, external :: typewriter

  mbond = 0
  do a = 1, na1uc;  atom => unitcell%atom(a)
     mbond = max(mbond, atom%nbond)
  end do

  if(.not. allocated(nghbr) ) allocate( nghbr(mbond, naa) )

  do ia = 1, naa;   Ri = A2R(:, ia);  a = Ri(ndim);  atom => unitcell%atom(a)
     do ib = 1, atom%nbond;  b = atom%nbatom(ib);  Rj(ndim) = b
		call connection( Ri(1:nd), a, atom%bond(:, ib), Rj(1:nd), b, info)
		if(info == 0)stop 'connection error'
		nghbr(ib, ia) = typewriter(ndim, limits, Rj) * info         !info encodes boundary condition
     end do
  end do
  
  return
end subroutine coordination

subroutine includenn( )
  use vmcplace
  implicit none

  integer ia, ib, a, b, info
  integer, dimension (ndim) :: Ri, Rj
  real*8 bond(3), shortbond(3)
  real*8 lsq, shortlsq
  logical hexagonal, active
  integer typewriter
  type (atomconfig), pointer :: atom

  hexagonal = active( sum( unitcell%abc(:, 1) * unitcell%abc(:, 2) ) )
  mnn = 4; if(hexagonal) mnn = 6
  if(nd == 1) mnn = 2
    
  shortbond = 100
  shortlsq = sum(shortbond * shortbond)
  
  do ia = 1, unitcell%natom;  atom => unitcell%atom(ia)
     do ib = 1, atom%nbond
	    bond = atom%bond(:, ib);  lsq = sum(bond*bond)
		if( lsq < shortlsq )then; shortlsq = lsq; shortbond = bond; end if
	 end do
  end do

  if( sum( abs(shortbond(1:2)) ) < 1.d-5 )mnn = 2    !vertical shortest bond

  if(.not. allocated(nn) ) allocate( nn(mnn, naa) );  nn = 0

  do ia = 1, naa
     bond = shortbond; Ri = a2R(:, ia); a = Ri(ndim)
	 do ib = 1, mnn;  if(ib>1)call rotate(bond)
		do b = 1, unitcell%natom
		   call connection(Ri(1:nd), a, bond, Rj(1:nd), b, info);  if(info == 0)cycle
		   Rj(ndim) = b;  nn(ib, ia) = typewriter(ndim, limits, Rj)
		end do
	  end do
   end do

   return
   contains
     subroutine rotate(bond)
	   real*8 bond(3)
	   if(mnn == 2)then; bond = - bond; return; end if
	   if(hexagonal)then; call rotate60(bond(1:2)); else; call rotate90(bond(1:2)); end if
       return
     end subroutine rotate
end subroutine includenn

subroutine neighbors(unitcell)
  use crystal
  implicit none
  
  type (unitcellconfig), target :: unitcell

  type (atomconfig), pointer :: atom, btom
  integer na, nb, ia, ja, ib
  real*8, dimension (3) :: rv, rt
  integer, external :: findatom

  na = unitcell%natom
  do ia = 1, na;  atom => unitcell%atom(ia)
     nb = atom%nbond;  if(nb == 0)cycle
	 allocate( atom%nbatom(nb) )
	 do ib = 1, nb;  rv = atom%r + atom%bond(:, ib)
	    ja = findatom(rv, unitcell, rt)
		if(ja <= 0)stop 'nbatom not found @ neighbors'
        atom%nbatom(ib) = ja
	 end do
  end do     

  return
end subroutine neighbors

function findatom(rv, unitcell, rt)
  use crystal
  implicit none
  integer findatom
  real*8, dimension (3) :: rv, rt
  type (unitcellconfig), target :: unitcell
  type (atomconfig), pointer :: atom

  real*8 dr(3)
  real*8, pointer, dimension (:, :) :: abc, dualabc
  integer natom, i, ia, ib, ic

  natom=unitcell%natom
  if(natom == 1)then; findatom = 1; rt = rv; return; end if

  abc => unitcell%abc;  dualabc => unitcell%dualabc

  findatom = -1
  do i = 1, natom; atom => unitcell%atom(i)
     dr = rv - atom%r
     ia = nint( sum( dr * dualabc(:, 1) ) )
	 ib = nint( sum( dr * dualabc(:, 2) ) )
     ic = nint( sum( dr * dualabc(:, 3) ) )
     dr = dr - ia * abc(:, 1) - ib * abc(:, 2) - ic * abc(:, 3)
	 if( sum( abs(dr) ) < 1.d-5 )then
		findatom = i;  rt = ia * abc(:, 1) + ib * abc(:, 2) + ic * abc(:, 3); return
	 end if
  end do
  return
end function findatom