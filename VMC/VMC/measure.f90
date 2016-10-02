  
subroutine measure()
  use vmcplace
  implicit none

  real*8 Eonsite, Ebond, E_s
  integer ia, nkondo
  integer, external :: moment
  real*8, external :: NeelMoment
  
  if(LMoptimize)call getdlogs( )
   
  E_s = ( Eonsite() + Ebond() ) / naa
  Eav = Eav + E_s;  delE = delE + E_s * E_s

  rho_s = rho_s / naa
  rhoav = rhoav + rho_s;  delrho = delrho + rho_s * rho_s

  call density_distribution( )

  isample =  isample + 1

  return
end subroutine measure

subroutine getdlogs( )
  use vmcplace
  implicit none
  
  complex*16, dimension (npair, npair) :: Dwork, bra(LM, 1), ket(1, LM)
  integer i, j, ia, ja, iv
  complex*16, external :: element, trace

  if(.not. allocated(dlogs) ) allocate( dlogs(nvar) )

  if(conjugate .and. .not. allocated(ratios) ) allocate( ratios(LM) )

  if(bruteforce)then
    if(.not. allocated(G) ) allocate( G(npair, npair) )
	G = D;  call ZINVERT(npair, G)
  end if

  do iv = 1, nfield
     whichtable = iv;   Dwork = 0
     do ja = 1, naa; j = occup(ja, 2);  if(j == 0)cycle
     do ia = 1, naa; i = occup(ia, 1);  if(i == 0)cycle
        Dwork(i, j) = element(  a2R(:, ia),  a2R(:, ja) )
     end do; end do  
	 if(conjugate) Dv(:, :, iv) = Dwork
     dlogs(iv) = trace(npair, G, Dwork)
  end do

  if(nfield < nvar) call dlogPenalty( )

  if(newton)then
    Avdlog = Avdlog + dlogs
  else
    bra(1, 1) = 1;  bra(2:LM, 1) = conjg(dlogs)
    ket(1, 1) = 1;  ket(1, 2:LM) = dlogs
    Svv = Svv + matmul(bra, ket)
  end if 

  return
end subroutine getdlogs

subroutine density_distribution( )
  use vmcplace
  implicit none

  integer ia, m, c
  integer, external :: moment, charge
  
  if( allocated(sdwav) )then
    do ia = 1, naa;  m = moment( occup(ia, :) )
       sdwav(ia) = sdwav(ia) + m;  delsdw(ia) = delsdw(ia) + m * m
    end do
  end if

  if( allocated(cdwav) )then
    do ia = 1, naa;    c = charge( occup(ia, :) )
       cdwav(ia) = cdwav(ia) + c;  delcdw(ia) = delcdw(ia) + c * c
    end do
  end if

  return
end subroutine density_distribution


function Eonsite()
  use vmcplace
  implicit none
  real*8 Eonsite

  integer ia
  integer, dimension (2) :: si
  integer, external :: doublon, charge, moment
  complex*16 bra(LM, 1), ket(1, LM)
  type (atomconfig), pointer :: atom

  Eonsite = 0

  do ia = 1, naa; si =  occup(ia, :)

     atom => unitcell%atom( a2R(ndim, ia) )
     Eonsite = Eonsite + atom%E * charge(si) - atom%B * moment(si)

	 if( mott .and. doublon(si) )stop 'doublon found for mott system'

	 if(.not. mott .and. .not. kondo(ia) ) Eonsite = Eonsite + doublon(si) *  U

  end do

  if(.not. LMoptimize) return
  
  if(newton)then
    Edlog = Edlog + Eonsite * dlogs / naa
  else
    bra(1, 1) = 1;  bra(2:LM, 1) = conjg(dlogs)
	ket(1, 1) = 1;  ket(1, 2:LM) = dlogs
	Evv = Evv + matmul(bra, ket) * Eonsite / naa
  end if

  return
end function Eonsite

function Ebond()
  use vmcplace
  implicit none
  real*8 Ebond

  integer ia, ja, a, b, ib, spin, info
  character*10 mode
  integer, dimension (2) :: si, sj
  logical superExchangeable, Exchangeable, inactive, active
  integer moment, charge, doublon, typewriter
  integer, dimension( ndim) :: Ri, Rj
  logical kondobond
  real*8 add, bond(2)
  complex*16 teff, Evbond(LM), bra(LM, 1), ket(1, LM)
  type (atomconfig), pointer :: atom

  Ebond = 0;  mode = 'measure'
  rho_s = 0    !prepare to measure rho in the present config

  if(LMoptimize .and. conjugate)then
    Evbond = 0
    bra(1, 1) = 1;  bra(2:LM, 1) = conjg(dlogs)
	ket(1, 1) = 1;  ket(1, 2:LM) = dlogs
  end if

  do ia = 1, naa; si = occup(ia, :)
     Ri = A2R(:, ia);  a = Ri(ndim);  atom => unitcell%atom(a)
     do ib = 1, atom%nbond; b = atom%nbatom(ib);  bond = atom%bond(1:2, ib)

		!Rj(ndim) = b
        !call connection(Ri(1:nd), a, atom%bond(:, ib), Rj(1:nd), b, info, vmcin)
		!if(info==0)cycle     !Ri + ra is not connected via the bond to an orbital b in some other unitcell
		!ja = typewriter(ndim, limits, Rj);  sj = occup(ja, :)
		
		ja = nghbr(ib, ia);  info = 1
		if(ja < 0)then;  ja = -ja;  info = -1;  end if
		sj = occup(ja, :)

		kondobond = ( kondo(ia) .or. kondo(ja) )
        
		if( active( abs(atom%t(ib)) + abs(atom%tz(ib)) ) )then
		   do spin = 1, 2
		      
			  if( kondobond )cycle
			  !hop to/from a kondo site is forbidden

		      if(.not. exchangeable(spin, si, sj, mott) )cycle
		      call hopping(spin, ia, ja, mode)    

		      ratio = conjg( ratio ) * info                  
			  !the complex conjugation is needed by definition, and info accounts for boundary condition
			  
			  teff = atom%t(ib) + atom%tz(ib) * (3-2*spin)            !for C_i^dag C_j
		      if( si(spin) /= 0)teff = conjg(teff)                    !for C_j^dag C_i

		      Ebond = Ebond - real( teff * ratio ) 
              rho_s = rho_s + real( teff * ratio ) * bond * bond
			  if(LMoptimize .and. conjugate) Evbond = Evbond - teff * conjg(ratios) * info
			  			         
		   end do
		end if

		if( active( atom%J(ib) ) .or. active( atom%V(ib) ) )then
          add = atom%V(ib) * charge(si) * charge(sj) + 0.25 * atom%J(ib) * moment(si) * moment(sj)
          Ebond = Ebond + add
          if(LMoptimize .and. conjugate) Evbond = Evbond + add * bra(:, 1)
		end if

        if( inactive( atom%J(ib) ) )cycle
		if( .not. superExchangeable(si, sj) )cycle

        !spin flip translated as superexchange, no boundary condition is needed here
		!since each kind of sign change occurs twice
		call superExchange(ia, ja, mode)
	    Ebond = Ebond - 0.5 * atom%J(ib) * real(ratio)
		if(LMoptimize .and. conjugate) Evbond = Evbond - 0.5 * atom%J(ib) * conjg(ratios) 

     end do
  end do

  if(.not. LMoptimize)return

  if(newton)then
    Edlog = Edlog + Ebond * dlogs / naa
  else
    bra(:, 1) = Evbond / naa;  Evv = Evv + matmul(bra, ket)
  end if

  return
end function Ebond



