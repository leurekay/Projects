subroutine hopping(spin, ia, ja, mode)
  use vmcplace
  implicit none

  integer spin, ia, ja          !|ia> and |ja> are two single-particle basis states
  character*10 mode             !mode = 'bath', 'sampling', 'measure'

  real*8 punish, p                                 !metropolis probability
  complex*16 detratio                              !ratio between new and old determinants
  complex*16 vec( Npair)                           !row or column updated by relocating an electron
  complex*16, dimension (npair, npair) :: Dwork, Gwork

  integer iv, who                 !who (with 'spin') is going to hop
  integer Rf( ndim)               !|Rf> is the final sps after fermion hopping
  logical exchangeable            !user specified logical func, judge whether exchange of sps is valid
  real*8, external :: penalty     !user specified penalty func for hop (which could change both spin and charge configs)
  real*8, external :: rand
  complex*16, external :: trace

  if(mode /= 'measure') nattempt =   nattempt + 1  

  punish = penalty(ia, ja, spin)

  if(occup(ia, spin) /= 0 )then
     Rf = A2R(:, ja); who = occup(ia, spin)
  else
     Rf = A2R(:, ia); who = occup(ja, spin)
  end if
  
  whichtable = 0
        
  !calculate new row or column due to electron relocation
  call newvector(spin, Rf, vec)  
    
  if(mode == 'measure' .and. LMoptimize .and. conjugate)then   

    ratios(1) = 1;  Dwork = D;  Gwork = G

	call dyson(Npair, Dwork, Gwork, who, vec, spin, detratio, 1)
	ratio = detratio * punish
	
	do iv = 1, nfield; whichtable = iv;  Dwork = Dv(:, :, iv)
	   call newvector(spin, Rf, vec)
	   if(spin == 1)then; Dwork(who, :) = vec;  else; Dwork(:, who) = vec; end if
	   ratios(iv+1) = trace(npair, Gwork, Dwork)
	end do
	
	if(nfield < nvar) call left_dlogpenalty(ia, ja, spin)
	
	ratios = ratios * ratio;  return
  end if

  !calculate ratio between new and old determinants
  if( bruteforce)then
     call brute1vec(nPair, D, who, vec, spin, det, detratio, 0)
  else
     call dyson(Npair, D, G, who, vec, spin, detratio, 0)
  end if

  !combine penalty ratio 		
  ratio = punish * detratio 
    
  if(mode=='measure')return

  !metropolis probability
  p = abs(ratio)**2; if(mode == 'bath')p = p/(1+p)
  if(rand(iseed) > p)return

  !update p2p and G
  if( bruteforce)then
    call brute1vec(nPair, D, who, vec, spin, det, detratio, 1)
  else
    call dyson(Npair, D, G, who, vec, spin, detratio, 1)
  end if

  !update coordinates and state
  call exchange(occup(ia, spin), occup(ja, spin))

   naccepted =    naccepted + 1

  return
end subroutine hopping


subroutine superexchange(ia, ja, mode)
  use vmcplace
  implicit none

  integer iv, ia, ja               !|ia> and |ja> are two single-particle basis states
  character*10 mode           !mode = 'bath', 'sampling', 'measure'

  complex*16 vec( Npair, 2)

  real*8 punish, p                                        !metropolis probability
  complex*16 detratio                             !ratio between new and old determinants
  real*8, external :: penalty
  logical, external :: active

  integer spin
  integer Rf(ndim, 2), who(2)
  complex*16, dimension (npair, npair) :: Dwork, Gwork
  real*8, external :: rand
  complex*16, external :: trace

  if(mode/='measure')  nattempt =   nattempt + 1

  punish = penalty(ia, ja, 0)
  !penalty func must be calculated first since occup is changed in the interdiate stage of 
  !the following spin-exchange attempt 

  whichtable = 0

  do spin = 1, 2
     if( occup(ia, spin) /= 0 )then
       Rf(:, spin) = A2R(:, ja); who(spin) = occup(ia, spin)
     else
       Rf(:, spin) = A2R(:, ia); who(spin) = occup(ja, spin)
     end if
     call newvector( spin, Rf(:, spin), vec(:, spin) )
	 
	 !assuming spin-up already hopped, change occup(:, spinup)
	 if(spin==1)call exchange( occup(ia, 1), occup(ja, 1) )
  end do  

  if(mode == 'measure' .and. LMoptimize .and. conjugate)then   
    
	ratios(1) = 1;  Dwork = D;  Gwork = G
	call rowcol(Npair, who(1), vec(:, 1), who(2), vec(:, 2), Dwork, Gwork, detratio, 1)
	ratio = detratio * punish
	
	do iv = 1, nfield; whichtable = iv;  Dwork = Dv(:, :, iv)
	   call exchange( occup(ia, 1), occup(ja, 1) )                   !restore occup
       do spin = 1, 2
	      call newvector( spin, Rf(:, spin), vec(:, spin) )
		  if(spin == 1) call exchange( occup(ia, 1), occup(ja, 1) )  !change occup for spin up
	   end do
	   Dwork(who(1), :) = vec(:, 1);  Dwork(:, who(2)) = vec(:, 2)
	   ratios(iv+1) = trace(npair, Gwork, Dwork)
	end do

	call exchange( occup(ia, 1), occup(ja, 1) )    !restore occup
	if(nfield < nvar) call left_dlogpenalty(ia, ja, 0)

	ratios = ratios * ratio;  return

  end if

  if( bruteforce)then
    call brute1r1c(nPair, who(1), vec(:, 1), who(2), vec(:, 2), D, det, detratio, 0)
  else
    call rowcol(Npair, who(1), vec(:, 1), who(2), vec(:, 2), D, G, detratio, 0)
  end if
  ratio = punish * detratio 


  !restore spin-up occup if the superexchange is only for measurement purpose 
  if(mode == 'measure')then;  call exchange( occup(ia, 1), occup(ja, 1) );  return; end if
    
  !metropolis probability
  p=abs(ratio)**2; if(mode == 'bath') p = p/(1+p)

  !restore spin-up occup if attempted superexchange is not accepted
  if(rand(iseed) > p)then; call exchange(occup(ia, 1), occup(ja, 1)); return; end if

  !update p2p, G and spin-down occup (spin-up occup already changed in the attempt)
  if( bruteforce)then
    call brute1r1c(nPair, who(1), vec(:, 1), who(2), vec(:, 2), D, det, detratio, 1)
  else
    call rowcol(Npair, who(1), vec(:, 1), who(2), vec(:, 2), D, G, detratio, 1)
  end if
  call exchange( occup(ia, 2), occup(ja, 2) )

  naccepted =   naccepted + 1

  return
end subroutine superexchange


subroutine update(mode)
  use vmcplace
  implicit none
  character*10 mode

  call updateOnsite(mode)
  call updateBond(mode)
  
  return
end subroutine update  

subroutine updateOnsite(mode)
  use vmcplace
  implicit none
  character*10 mode

  integer ia
  integer, dimension(2) :: si
  integer, external :: doublon, moment

  !in principle, local interactions are well defined only for Hubbard-type models
  !there is nothing to be done if an atom contains one orbital only, except for a check of doublon

  !CAUTION: this subroutine needs to be significantly modified if an atom contains many orbitals, whence 
  !spin exchange and pair hopping are allowed between orbitals, or orbital hybridization when soc is present

  if(.not.  mott)return

  do ia = 1,  naa;  si =  occup(ia, :)
	 if(  mott .and. doublon(si) )stop 'doublon found for mott system'
	 if( moment(si) == 0 .and.  kondo(ia) )stop 'no spin found on kondo site'
  end do

  return
end subroutine updateOnsite

subroutine updateBond(mode)
  use vmcplace
  implicit none
  character*10 mode

  integer ia, ja, a, b, spin, m, mb, ib, info
  integer, dimension(2) :: si, sj
  logical superExchangeable, Exchangeable, pairExchangeable
  integer typewriter
  integer, dimension( ndim) :: Ri, Rj
  integer, external :: irand
  type (atomconfig), pointer :: atom, btom
  logical kondobond

  do m = 1, naa; ia = irand(iseed, naa)
     Ri = A2R(:, ia);  a = Ri(ndim);  atom => unitcell%atom(a)
     do mb = 1, atom%nbond; ib = irand(iseed, atom%nbond);  b = atom%nbatom(ib)

		!Rj(ndim) = b
		!call connection( Ri(1:nd), a, atom%bond(:, ib), Rj(1:nd), b, info, vmcin)
		!if(info == 0)stop 'connection error found'      !Ri + ra is not connected via the bond to an atom b in some other unitcell

		!ja = typewriter(ndim, limits, Rj)

		ja = abs( nghbr(ib, ia) )
		kondobond = ( kondo(ia) .or. kondo(ja) )
		
		!hopping on bond
		do spin = 1, 2

		   if( kondobond ) cycle
		   !hopping to/from a kondo site is invalid
           
		   si = occup(ia, :); sj = occup(ja, :)                
		   !si and sj should be assigned each time in use, since they may have been modified implicitly 
		   
		   if( exchangeable(spin, si, sj, mott) )call hopping(spin, ia, ja, mode)
        end do

        !spin exhcange on bond
		si = occup(ia, :); sj = occup(ja, :)
	    !si and sj should be assigned each time in use, since they may have been modified implicitly 
		if( superExchangeable(si, sj))call superExchange(ia, ja, mode)

        !pair hopping on bond
		
		if( kondobond )cycle
		!pairing hopping to/from a kondo site is forbidden

        si = occup(ia, :); sj = occup(ja, :)
 	    !si and sj should be assigned each time in use, since they may have been modified implicitly 
		if( pairExchangeable(si, sj) )call superExchange(ia, ja, mode)

     end do
  end do

  return
end subroutine updateBond