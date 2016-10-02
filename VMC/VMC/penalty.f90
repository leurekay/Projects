
function penalty(ia, ja, info)
  use vmcplace
  implicit none

  integer ia, ja, info
  real*8 penalty

  !info = 0:  spin exchange if moment(ia) * moment(ja) /= 0,  or pair hopping if moment(ia) = moment(ja) = 0
  !info /= 0: electron with spin = info hops between isps and jsps, both spin and charge configs are changed

  integer nu, inn, k, momenti, momentj, chrgi, chrgj, deld
  integer, dimension (2) :: si, sj
  real*8 delm, delc, deldh, factor
  integer, external :: moment, charge, doublon, doublonholon
  logical, external :: inactive

  penalty = 1
  
  if(  mott .and. inactive( abs( Vzz) + abs( Vcc) ) )return    
  !in mott system, Vdh is inactive since no doublon occurs

  factor = abs( Vzz) + abs( Vcc) + abs( Vdh) + abs( gU) 
  if(.not. mott .and. inactive(factor) )return

  si = occup(ia, :);  sj = occup(ja, :)
  
  delm = 0;  delc = 0;  deldh = 0;  deld = 0

  do nu = -1, 1, 2
     if( nu == 1)then
       !exchange states on isps and jsps according to info
       if( info == 0 )then
         !spin exchange on isps and jsps
         call exchange( occup(ia, 1), occup(ja, 1) ); call exchange( occup(ia, 2), occup(ja, 2) )
       else
         !electron with spin = info hops between isps and jsps
         call exchange( occup(ia, info), occup(ja, info) )
       end if
     end if

     deld = deld + nu * ( doublon( occup(ia, :) ) + doublon( occup(ja, :) ) )

	 if( inactive( abs( Vzz) + abs( Vdh) + abs( Vcc) ) ) cycle

     momenti = moment( occup(ia, :) );  momentj = moment( occup(ja, :) )
     chrgi = charge( occup(ia, :) );    chrgj = charge( occup(ja, :) )

     do inn = 1,  mnn

	    !nn of isps
	    k = nn(inn, ia)
		if(k > 0)then       !make sure that k is an nn of isps
          factor = 1;	  if(k == ja)factor = 0.5        
		  !the factor of 1/2 reduces double counting
	      
		  delm = delm + nu * moment(occup(k, :)) * momenti * factor
		  delc = delc + nu * charge(occup(k, :)) * chrgi * factor
		  deldh = deldh + nu * doublonholon( charge(occup(k, :)), chrgi ) * factor
        end if

        !nn of jsps
        k = nn(inn, ja)
        if(k > 0)then       !make sure that k is an nn of jsps
		  factor = 1;	  if(k == ia)factor = 0.5        
		  !the factor of 1/2 reduces double counting
	      
		  delm = delm + nu * moment(occup(k, :)) * momentj * factor
		  delc = delc + nu * charge(occup(k, :)) * chrgj * factor
		  deldh = deldh + nu * doublonholon( charge(occup(k, :)), chrgj ) * factor
        end if
	 end do
  end do
    
  penalty = exp( - Vzz * delm +  Vdh * deldh -  Vcc * delc)
  if(.not.  mott) penalty = penalty * (1 -  gU) ** deld 

  occup(ia, :) = si;  occup(ja, :) = sj

  return
end function penalty

subroutine dlogPenalty( )
  use vmcplace
  implicit none

  integer ia, ja, inn, k
  integer momenti, momentj, chrgi, chrgj
  integer Dtot, DHtot, zztot, cctot 
  integer, dimension (2) :: si, sj
  integer, external :: moment, charge, doublon, doublonholon

  Dtot = 0;  DHtot = 0;  zztot = 0;  cctot = 0

  do ia = 1, naa;  si = occup(ia, :)

     Dtot = Dtot + doublon( si )

     momenti = moment( si );  chrgi = charge( si )
	 
     do inn = 1,  mnn

	    !nn of isps
	    k = nn(inn, ia);  if( k <= 0 ) cycle

		sj = occup(k, :);  chrgj = charge( sj );  momentj = moment( sj )

		DHtot = DHtot + doublonholon( chrgi, chrgj )
		zztot = zztot + momenti * momentj
		cctot = cctot + chrgi * chrgj
     end do
  end do
    
  where( vname(1:nvar) == 'gU') dlogs = - Dtot / (1-gU)
  where( vname(1:nvar) == 'Vdh') dlogs = DHtot * 0.5
  where( vname(1:nvar) == 'Vzz') dlogs = - zztot * 0.5
  where( vname(1:nvar) == 'Vcc') dlogs = - cctot * 0.5

  return
end subroutine dlogPenalty



subroutine left_dlogpenalty(ia, ja, info)
  use vmcplace
  implicit none

  integer ia, ja, info

  !info = 0:  spin exchange if moment(ia) * moment(ja) /= 0,  or pair hopping if moment(ia) = moment(ja) = 0
  !info /= 0: electron with spin = info hops between isps and jsps, both spin and charge configs are changed

  integer nu, inn, k, momenti, momentj, chrgi, chrgj
  integer, dimension (2) :: si, sj
  real*8 deld, delm, delc, deldh, factor
  integer, external :: moment, charge, doublon, doublonholon

  si = occup(ia, :);  sj = occup(ja, :)
  
  delm = 0;  delc = 0;  deldh = 0;  deld = 0

  do nu = -1, 1, 2
     if( nu == 1)then
       !exchange states on isps and jsps according to info
       if( info == 0 )then
         !spin exchange on isps and jsps
         call exchange( occup(ia, 1), occup(ja, 1) ); call exchange( occup(ia, 2), occup(ja, 2) )
       else
         !electron with spin = info hops between isps and jsps
         call exchange( occup(ia, info), occup(ja, info) )
       end if
     end if

     deld = deld + nu * ( doublon( occup(ia, :) ) + doublon( occup(ja, :) ) )

     momenti = moment( occup(ia, :) );  momentj = moment( occup(ja, :) )
     chrgi = charge( occup(ia, :) );    chrgj = charge( occup(ja, :) )

     do inn = 1,  mnn

	    !nn of isps
	    k = nn(inn, ia)
		if(k > 0)then       !make sure that k is an nn of isps
          factor = 1;	  if(k == ja)factor = 0.5        
		  !the factor of 1/2 reduces double counting
	      
		  delm = delm + nu * moment(occup(k, :)) * momenti * factor
		  delc = delc + nu * charge(occup(k, :)) * chrgi * factor
		  deldh = deldh + nu * doublonholon( charge(occup(k, :)), chrgi ) * factor
        end if

        !nn of jsps
        k = nn(inn, ja)
        if(k > 0)then       !make sure that k is an nn of jsps
		  factor = 1;	  if(k == ia)factor = 0.5        
		  !the factor of 1/2 reduces double counting
	      
		  delm = delm + nu * moment(occup(k, :)) * momentj * factor
		  delc = delc + nu * charge(occup(k, :)) * chrgj * factor
		  deldh = deldh + nu * doublonholon( charge(occup(k, :)), chrgj ) * factor
        end if
	 end do
  end do
    
  where( vname(1:nvar) == 'gU') ratios(2:LM) = dlogs - deld / (1-gU)
  where( vname(1:nvar) == 'Vdh') ratios(2:LM) = dlogs + deldh
  where( vname(1:nvar) == 'Vzz') ratios(2:LM) = dlogs - delm
  where( vname(1:nvar) == 'Vcc') ratios(2:LM) = dlogs - delc

  occup(ia, :) = si;  occup(ja, :) = sj

  return
end subroutine left_dlogpenalty