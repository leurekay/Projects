subroutine MCthread(  )
  use vmcplace
  implicit none
  character*10 mode
  integer mc

  call initconfig(  )
  
  !initial zeros
  rhoav = 0;      delrho = 0
  Eav = 0;      Ekav = 0;       Epav = 0;      delE = 0  
  isample = 0;  nattempt = 0;   naccepted = 0

  if( allocated(sdwav) )then;  sdwav = 0;  delsdw = 0;  end if
  if( allocated(cdwav) )then;  cdwav = 0;  delcdw = 0;  end if

  if(LMoptimize)then
    if(newton)then; Edlog = 0;  Avdlog = 0; end if
	if(conjugate)then; Svv = 0;  Evv = 0;  end if
  end if

  mc = 0
  do while ( isample < nsample )

     mc = mc + 1

     mode = 'bath'; if(mc > Nbath) mode = 'sampling'
	 !adjust in Metropolis update mode
     
	 call update(mode)

     !if(nattempt == 0)stop 'vmc failed to find new config'
	 !something is wrong if no attempt was made in MC

     if(mc > Nbath .and. mod(mc, Nmc1sample) == 0)call measure()
	 !do measurement every Nmc1sample scans

     if(bruteforce) cycle

	 call checkG(npair, D, G)          
	 !check accuracy of G obtained by Dyson in MC

  end do

  acceptance =  naccepted * 1./ nattempt

  rhoav= rhoav / nsample;       delrho = delrho / nsample - rhoav * rhoav
  where(delrho > 0) delrho = sqrt( delrho / (nsample - 1) )

  Eav =  Eav / nsample;      delE =  delE / nsample -  Eav * Eav
  if(delE > 0) delE = sqrt(  delE / (nsample-1) )

  if( allocated(sdwav) )then
    sdwav = sdwav / nsample;  delsdw = delsdw / nsample - sdwav * sdwav
    where(delsdw > 0) delsdw = sqrt( delsdw / (nsample-1) )
  end if
  
  if( allocated(cdwav) )then
    cdwav = cdwav / nsample;  delcdw = delcdw / nsample - cdwav * cdwav
    where(delcdw > 0) delcdw = sqrt( delcdw / (nsample-1) )
  end if

  if(.not. LMoptimize) return

  if(newton)then; Edlog = Edlog / nsample;  Avdlog = Avdlog / nsample; end if
  if(conjugate)then; Evv = Evv / nsample;  Svv = Svv / nsample;  end if
    
  return
end subroutine MCthread



