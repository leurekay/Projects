subroutine initconfig( )
  use vmcplace
  implicit none

  integer i, j, ia, ja, spin, isps, jsps, icount, jcount
  integer milisecond
  logical exchangeable, superExchangeable
  real*8, external :: rand
  complex*16, external :: element, ZDET
  logical mott0

  !get seed number
  iseed = milisecond() + ithread * 10000

  mott0 = mott
  do i = 1, naa; if(  kondo(i) ) mott = .true.; end do
  !make sure that initially no double occupancy is allowed if there are kondo impurities

  !prepare occup (append history if applicable)
  if(.not. allocated(occup) )then
    allocate(  occup(naa, 2) );    occup = 0
	do i = 1, Npair; occup(i, 1) = i; occup(naa-i+1, 2) = i; end do

    !randomize electron positions via random hopping
    do spin = 1, 2; do ia = 1, naa;  do ja = 1, naa
       if( .not. exchangeable( spin, occup(ia, :), occup(ja, :), mott ) )cycle
       if(rand(iseed) > 0.5)call exchange(occup(ia, spin), occup(ja, spin))
    end do; end do; end do
    
    !randomize electron positions via random spin exchange
    do ia = 1, naa; do ja = 1, naa
       if( .not. superExchangeable( occup(ia, :), occup(ja, :) ) ) cycle
       if(rand(iseed) > 0.5) cycle
	   call exchange(occup(ia, 1), occup(ja, 1))
	   call exchange(occup(ia, 2), occup(ja, 2))
    end do; end do
  end if

  if(.not. allocated(D) )allocate( D(Npair, Npair) )
  if(.not. allocated(G) )allocate( G(Npair, Npair) )

  !initial p2p matrix and its inverse  

  D = 0;  G =0;  whichtable = 0
  do ja = 1, naa; j = occup(ja, 2);  if(j == 0)cycle
  do ia = 1, naa; i = occup(ia, 1);  if(i == 0)cycle
     D(i, j) = element(  a2R(:, ia),  a2R(:, ja) )
  end do; end do  
  
  if(bruteforce)then; det = ZDET(Npair, D); else; G = D; call ZINVERT(Npair, G); end if
  
  mott = mott0
  return
end subroutine initconfig