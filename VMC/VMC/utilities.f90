function trace(n, A, B)
  implicit none
  integer n
  complex*16 A(n, n), B(n, n), trace

  integer i
  
  trace = 0
  do i = 1, n
     trace = trace + sum( A(i, :) * B(:, i) )
  end do

  return
end function trace

function Jastraw( vname )
  implicit none
  character*6 vname
  logical Jastraw

  Jastraw = 0
  if( vname == 'gU' )Jastraw = 1
  if( vname == 'Vdh' )Jastraw = 1
  if( vname == 'Vzz' )Jastraw = 1
  if( vname == 'Vcc' )Jastraw = 1

  return
end function Jastraw

subroutine register(n, table, v, nvar, vname, vmin, vmax, var, model, optimize)
  implicit none
  integer n, nvar, info
  logical optimize
  character*6 table(*), vname(*)
  character*12 model
  real*8 v(*), var(*), vmin(*), vmax(*)

  integer i
  character*6 name
  real*8 lower, upper
  logical, external :: inactive

  open(10, file = 'var.'//model)
  
  nvar = 0     !counter of vars to be optimized

  if(optimize)then
    do i = 1, n	 	    
 	   read(10, *)name, lower, upper
	   if(lower > upper)stop 'invalid var bound found in file'
	   if( name /= table(i) )stop 'vname mismatch in file'

	   v(i) = lower;  if( inactive(upper - lower) ) cycle            ! v = lower = upper
		
	   nvar = nvar + 1;     v(i) = var(nvar)                         ! v determined outside by user
	   vname(nvar) = name;  vmin(nvar) = lower;  vmax(nvar) = upper  ! register name and bounds
    end do
    if(nvar == 0)stop 'no vars in file are to be optimized'
  else
    do i = 1, n
	   read(10, *)name, v(i)
	   if( name /= table(i) )stop 'vname mismatch in file'
	end do  
  end if
    
  close(10)

  return
end subroutine register 

function fermi(x)
  implicit none
  real*8 x, fermi

  fermi = exp( -abs(x) )
  fermi = 1.d0 / ( 1.d0 + fermi )
  if(x > 0) fermi = 1 - fermi

  return
end function fermi


function fermilevel(nval, eval, n)
  implicit none
  integer nval, n
  real*8 fermilevel, eval(nval)

  !purpose: determine the fermi level if n particles occupy nval states of energy eval

  integer imap(nval)
  real*8 bandwidth
  logical iwrite/.false./

  if(n < 0 .or. n > nval ) stop 'invalid filling @ fermilevel'

  call ascendingorder(nval, eval, imap)

  bandwidth = eval( imap(nval) ) - eval( imap(1) )

  if( n == 0 )then
    fermilevel = eval( imap(1) ) - bandwidth * 1.d-3
  else if ( n == nval )then
    fermilevel = eval( imap(nval) ) + bandwidth * 1.d-3
  else
    fermilevel = 0.5d0 * ( eval( imap(n) ) + eval( imap(n+1) ) )
    if( abs( fermilevel - eval( imap(n+1) ) ) < bandwidth * 1.d-10 )then
       !print*, eval(imap);  print*, 'fermilevel =', fermilevel
	   if(iwrite)print*, 'open shell found @ fermilevel'
    end if
  end if

  return
end function fermilevel


function typewriter(ndim, limits, coord)
  implicit none
  integer ndim, typewriter
  integer limits(ndim), coord(ndim)

  !mapping between coordinates and typewriter index
  integer idim, layer(ndim)
  
  layer(1)=1; do idim=2, ndim; layer(idim)=layer(idim-1)*limits(idim-1); end do       
  !layer(:) = (/1, 1*N1, 1*N1*N2, .../)

  typewriter=sum((coord-1)*layer)+1  
  
  return
end function typewriter

subroutine typewriter2coord( i, ndim, limits, coord )
  implicit none
  integer i, ndim
  integer limits(ndim), coord(ndim)

  !convert the typewriter index i to coordinates

  integer k, idim
  integer layer(ndim)

  layer(1)=1; do idim=2, ndim; layer(idim)=layer(idim-1)*limits(idim-1); end do       
  !layer(:) = (/1, 1*N1, 1*N1*N2, .../)

  k=i
  do idim=ndim, 1, -1     
     coord(idim)=(k-1)/layer(idim)+1; k=k-(coord(idim)-1)*layer(idim)
  end do

  return
end subroutine typewriter2coord


function exchangeable(spin, si, sj, mott)
  implicit none
  integer spin, si(2), sj(2)
  logical exchangeable, mott

  integer fi(2), fj(2)
  integer, external :: doublon

  fi = si;  where(fi > 0) fi = 1
  fj = sj;  where(fj > 0) fj = 1

  exchangeable = ( fi(spin) /= fj(spin) )
  if(.not. exchangeable)return
  if( .not. mott )return
  	 
  if(doublon(fi) .and. mott)stop 'doublon found before exchange'
  if(doublon(fj) .and. mott)stop 'doublon found before exchange'

  call exchange( fi(spin), fj(spin) )

  if(doublon(fi) .and. mott)exchangeable = .false.
  if(doublon(fj) .and. mott)exchangeable = .false.
       
  return
end function exchangeable

function superexchangeable(si, sj)
  implicit none
  integer si(2), sj(2)
  logical superexchangeable

  integer, external :: moment

  superexchangeable = ( moment(si) * moment(sj) < 0 )
       
  return
end function superexchangeable

function pairExchangeable(si, sj)
  logical pairExchangeable
  integer, dimension (2) :: si, sj
  integer charge

  pairExchangeable = ( (charge(si) * charge(sj) == 0) .and. (charge(si)+charge(sj) == 2) )

  return
end function PairExchangeable

function inactive(x)
  logical inactive
  real*8 x
  inactive = ( abs(x) < 1.d-5 )
  return
end function inactive

function active(x)
  logical active
  real*8 x
  active = ( abs(x) >= 1.d-5 )
  return
end function active


function onsite( nd, ra, rb )
  integer nd
  real*8, dimension(nd) :: ra, rb
  logical onsite
  onsite = ( sum( abs(ra-rb) ) < 1.d-5 )
  return
end function onsite

subroutine exchange(i, j)
  implicit none
  integer i, j, k
  k=i; i=j; j=k
  return
end subroutine exchange

function doublon(s)
  implicit none
  integer doublon, s(2)
  doublon = 0
  if(s(1)>0 .and. s(2)>0)doublon = 1
  return
end function doublon

function holon(s)
  implicit none
  integer holon, s(2)
  holon = 0
  if(s(1) == 0 .and. s(2) == 0)holon = 1
  return
end function holon

function charge(s)
  implicit none
  integer charge, s(2)
  integer i
  charge = 0
  do i=1, 2; if(s(i)>0)charge = charge + 1; end do
  return
end function charge

function moment(s)
  implicit none
  integer moment, s(2)
  moment = 0
  if(s(1)>0)moment = 1
  if(s(2)>0)moment = moment - 1
  return
end function moment

function doublonholon(chi, chj)
  implicit none
  integer chi, chj, doublonholon
  doublonholon = 0
  if( (chi == 0 .and. chj == 2) .or. (chi == 2 .and. chj == 0) ) doublonholon = 1
  return
end function doublonholon

subroutine rotate90(r)
  implicit none
  real*8 r(2), x
  x=r(1)
  r(1)=-r(2)
  r(2)=x
  return
end subroutine rotate90

subroutine rotate45(r)
  implicit none
  real*8 r(2), work(2)
  work=r
  r(1)=work(1)-work(2)
  r(2)=work(1)+work(2)
  r=r/sqrt(2.d0)
  return
end subroutine rotate45

subroutine rotate60(r)
  implicit none
  real*8 r(2), work(2), rt3
  work=r; rt3=sqrt(3.d0)
  r(1)=0.5d0*( work(1)- work(2)*rt3 )
  r(2)=0.5d0*( rt3*work(1) + work(2) )
  return
end subroutine rotate60

subroutine rotate120(r)
  implicit none
  real*8 r(2)
  integer i
  do i=1, 2; call rotate60(r); end do
  return
end subroutine rotate120


subroutine dualvector(a, b, c, ab)
  implicit none
  real*8, dimension (3) :: a, b, c, ab

  ab(1)=a(2)*b(3)-a(3)*b(2)
  ab(2)=a(3)*b(1)-a(1)*b(3)
  ab(3)=a(1)*b(2)-a(2)*b(1)

  ab=ab/sum(ab*c)
  return
end subroutine dualvector


function milisecond()
   implicit none
   integer milisecond

   integer date_time(8)
   character*10 b(3)

   call date_and_time(b(1), b(2), b(3), date_time)
   milisecond = date_time(8) + date_time(7)*1000

   !print *, 'date_time    array values:'
   !print *, 'year=', date_time(1)
   !print *, 'month_of_year=', date_time(2)
   !print *, 'day_of_month=', date_time(3)
   !print *, 'time difference in minutes=', date_time(4)
   !print *, 'hour of day=', date_time(5)
   !print *, 'minutes of hour=', date_time(6)
   !print *, 'seconds of minute=', date_time(7)
   !print *, 'milliseconds of second=', date_time(8)
   !print *, 'DATE=', b(1)
   !print *, 'TIME=', b(2)
   !print *, 'ZONE=', b(3)

   return
 end function milisecond


FUNCTION RAND(ISEED)
  IMPLICIT NONE
  INTEGER ISEED
  REAL*8 RAND
  !RAND = RAN(ISEED)
  DO WHILE (1)
     ISEED = MOD(16807*ISEED, 2147483647)
     RAND = 0.5D0 * ( 1 + REAL(ISEED)*4.6566128752458E-10  )
     IF(RAND > 0.D0 .and. RAND < 1.D0 )EXIT
  END DO
  RETURN
END FUNCTION RAND

FUNCTION IRAND(ISEED, L)
  IMPLICIT NONE
  INTEGER ISEED, L, IRAND
  REAL*8, EXTERNAL :: RAND
  IF(L == 0)STOP 'L = 0 INVALIDE @ IRAND'
  IF(L == 1)THEN; IRAND = 1; RETURN; END IF
  DO WHILE (1)
     IRAND = INT( L*RAND(ISEED) ) + 1
	 IF(IRAND >=1 .and. IRAND <= L)EXIT
  END DO
  RETURN
END FUNCTION IRAND


SUBROUTINE ASCENDINGLIMITS(M, S, SX, LL, UL)
  IMPLICIT NONE

  INTEGER M, LL, UL, ML
  REAL*8 SX, S(*)
  
  !PURPOSE: FIND THE UPPER AND LOWER LIMITS OF SX IN A REAL SEQUENCE OF ASCENDING ORDER

  IF( SX >= S(M) )THEN; UL=M+1; LL=M+1; RETURN; END IF
  IF( SX <= S(1) )THEN; UL=1; LL=1; RETURN; END IF

  UL = M; LL = 1
  DO WHILE (UL-LL > 1)
     ML = (UL+LL)/2
     IF( SX > S(ML) )THEN; LL = ML; CYCLE; END IF
 	 IF( SX < S(ML) )THEN; UL = ML; CYCLE; END IF
	 UL = ML; LL = ML
   END DO
   RETURN
END SUBROUTINE ASCENDINGLIMITS



SUBROUTINE ASCENDINGORDER(N, S, IMAP)
  IMPLICIT NONE
  INTEGER N, IMAP(*)
  REAL*8 S(*)
  
  REAL*8 SS(N)

  !PURPOSE: BRING A REAL SEQUENCE S INTO ASCENDING ORDER, S(IMAP)
  !NOTICE:  THE SEQUENCE S IS UNCHANGED ON EXIT, ONLY THE MAPPING IS OUTPUTED

  INTEGER I, LL, UL
  REAL*8 SI

  SS = S(1:N);  IMAP(1) = 1

  DO I = 2, N
     IF( SS(I) > SS(I-1) )THEN
	   IMAP(I) = I
	 ELSE
	   CALL ASCENDINGLIMITS( I-1, SS, SS(I), LL, UL)
	   SI = SS(I); SS(UL+1:I) = SS(UL:I-1); SS(UL) = SI
	   IMAP(UL+1:I) = IMAP(UL:I-1);  IMAP(UL) = I
     END IF
  END DO

  DO I = 2, N;  IF( SS(I) < SS(I-1) )STOP 'IASCENDING ORDER FAILED'; END DO
  DO I = 2, N;  IF( S( IMAP(I) ) < S( IMAP(I-1) ) ) STOP 'IMAP ERROR'; END DO

  RETURN
END SUBROUTINE ASCENDINGORDER
