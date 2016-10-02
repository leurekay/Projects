
subroutine gene_opt()
  use vmcplace
  implicit none

  integer, parameter :: ngen = 1000, nhead = 80
  real*8, parameter :: pmute = 0.15, convergence = 1.d-5

  integer nseg, jseed, i
  integer, external :: milisecond
  
  LMoptimize = 0;  optimize = 1;  optinfo = 0;  call definejob( )
  
  nseg = nvar;  if(nvar > 100)stop 'too many vars to be optimized by gene'

  if(ithread == 0)open(17, file = 'histo.'//model)
  
  if(ithread == 0) jseed = milisecond( );  call MPI_barrier( mpi_comm_world, ierr);  call broadcast_integer(jseed)

  call gene(nseg, nhead, ngen, jseed, pmute, convergence, fitness)

  if(ithread /= 0)return

  close(17)
  
  open(16, file = 'optimized.'//model);  write(16, *)'optimized variables:'
  do i = 1, nvar; write(16, *) vname(i), ' = ',  var(i);  end do;  close(16)

  open(16, file = 'rho.'//model)
  write(16, *)'rho = ', rhoav;  write(16, *)'delrho = ', delrho; close(16)

  return
  contains

  subroutine fitness(nseg, nhead, code, f)
    implicit none
    integer nseg, nhead
    integer code(nseg, nhead)
    real*8 f(nhead)

    integer i, iv, ibest, imap(nhead)
    real*8 ceil, delf(nhead)
    integer, external :: maxint
    logical iwrite/.false./

    call MPI_barrier(mpi_comm_world, ierr)

    ceil = dble( maxint( ) )
  
    do i = 1, nhead

       var(1:nseg) =  vmin(1:nseg) + (  vmax(1:nseg) -  vmin(1:nseg) ) * code(1:nseg, i) / ceil

       optinfo = optinfo + 1
       call definejob( );  call infrastructure( )
       call MCthread( );   call MPI_average( )
	 
	   f(i) = Eav;      delf(i) = delE

       if(.not. iwrite)cycle;   if(ithread /= 0)cycle

       print*, i
	   do iv = 1, nseg;  print*, vname(iv), '  ', var(iv);  end do
	   print*, 'Eav, delE =', Eav, delE

    end do

    call ascendingorder(nhead, f, imap);    ibest = imap(1)

    if(ithread == 0)then
	  do iv = 1, nseg; print*, vname(iv), '= ', var(iv);  end do
	  print*, 'best energy =', f(ibest);  print*, ' '
	  !output best var and the associated results
      var(1:nseg) =  vmin(1:nseg) + (  vmax(1:nseg) -  vmin(1:nseg) ) * code(:, ibest) / ceil    
	  if(nseg <= 4) write(17, 100) var(1:nseg), f(ibest), delf(ibest)
	  if(nseg > 4) write(17, 200) var(1:nseg), f(ibest), delf(ibest)
	end if

    !transform energy into fitness so that lowest energy corresponds to best fitness
	f = f(imap);    f = f(nhead) - f;     code = code(:, imap)

100 format(1x, 4f9.4, 4f10.5)
200 format(1x, e20.10)
    return
  end subroutine fitness

end subroutine gene_opt



!===========================================================================
!          portable gene algorithm

subroutine gene(nseg, nhead, ngen, iseed, pmute, convergence, uspec) 
  implicit none

  integer nseg, nhead, ngen, iseed
  real*8 pmute, convergence

  !arguments:
  !nseg, user specified number of segments in a gene code of an individual. each segment is represented by an integer
  !nhead, user specified population of the herd, which should be sufficiently large to garantee diversity and correct convergence
  !ngen,  user specified number of evolution iterations. the evolution may terminate earlier when convergence is reached
  !iseed, user specified seed for random number generator
  !pmute, user specified probability of gene mutation
  !cnvergence, user specified criterion for convergence of diversity of the fitness of the herd
  !uspec, dummy subroutine calculating the fitness of each individual in the herd, to be specified by user
  interface
    subroutine uspec(nseg, nhead, code, f)
	  integer nseg, nhead
	  integer code(nseg, nhead)
	  real*8 f(nhead)
	end subroutine uspec
  end interface

  integer nmax, bestcode(nseg), code(nseg, nhead), goodcode(nseg, nhead)
  real*8 bestf, worstf, f(nhead), psel(nhead), diversity
  logical selectbest

  integer, external :: irand, milisecond, maxint
  real*8, external :: rand

  integer i, k, ig, isg, ibt, lr

  !maximal positive integer of 4 byts ( 4*8 bits )
  nmax = maxint( )
  
  !initial herd with random gene codes
  do i = 1, nhead;  do k = 1, nseg
     code(k, i) = irand(iseed, nmax)
  end do;  end do

  do ig = 1, ngen  	 
	 
	 !fitness of the herd (arranged in descending order)
	 call uspec(nseg, nhead, code, f)	 	 
     bestf = f(1);  bestcode = code(:, 1)

     !check diversity of the herd in terms of fitness, and quit if converged
     worstf = f(nhead);  diversity = bestf - worstf

	 if( diversity <= convergence ) exit

	 !prob for an individual to survive to the next generation
	 psel = f / sum(f)  

     !select good individuals to survive
     k = 0;  selectbest = 0
	 do i = 1, nhead;  if( rand(iseed) > psel(i) )cycle
	    k = k + 1;  goodcode(:, k) = code(:, i)
		if(i == 1)selectbest = 1
	 end do

	 !priority for the best code
	 if(.not. selectbest)then;  k = k + 1; goodcode(:, k) = bestcode; end if

     !gene crossing between two neighboring (fitness ordered) individuals
	 do i = 1, nhead-1, 2

	    !select a crossing point
		isg = irand(iseed, nseg);  ibt = irand(iseed, 31) - 1
		
		!select a direction along which to exchange bits in codes
		lr = 1;  if(rand(iseed) > 0.5) lr = -1

		!exchange bits in codes i and i+1, along lr direction and starting from the crossing point (isg, ibt)
		call genecross( nseg, code(:, i), code(:, i+1), isg, ibt, lr )
	 end do

	 !priority to the survived good individuals
	 code(:, nhead-k+1:nhead) = goodcode(:, 1:k)

	 !gene mutating
	 do i = 1, nhead;  if( rand(iseed) > pmute)cycle
		isg = irand(iseed, nseg);  ibt = irand(iseed, 31) -1

		!flip ibt-th bit of isg-th segment in a code
	 	call iflip( code(isg, i), ibt )
	 end do

  end do

  return
end subroutine gene


subroutine genecross(nseg, code1, code2, isg, ibt, lr)
  implicit none
  integer nseg, isg, ibt, lr
  integer, dimension (nseg) :: code1, code2
  
  integer i, work

  if(lr == 1)then
    do i = isg + 1, nseg
       work = code1(i);  code1(i) = code2(i);  code2(i) = work
    end do

    do i = ibt + 1, 30
       call bitcross( code1(isg), i, code2(isg), i )
    end do
  else
    do i = 1, isg - 1
	   work = code1(i);  code1(i) = code2(i);  code2(i) = work
	end do

	do i = 0, ibt - 1
	   call bitcross( code1(isg), i, code2(isg), i )
	end do
  end if

  return
end subroutine genecross 


function maxint( )
  integer maxint

  integer i

  maxint = 1
  do i = 0, 30
     maxint = IBSET(maxint, i)
  end do
	
  return
end function maxint


subroutine iFlip(k, i)
  integer k, i

  !purpose: flip the i-th bit of a field k
  !notice: the 1st bit site is defined as 0-th bit
    
  if( BTEST(k, i) )then; k = IBCLR(k, i); else; k = IBSET(k, i); end if

  return
end subroutine iFlip


subroutine bitcross(k1, i1, k2, i2)
  integer k1, k2, i1, i2
	
  !purpose: exchange i1-th bit in k1 and i2-th bit in k2

  if( BTEST(k1, i1) )then
    if(.not. BTEST(k2, i2) )then
      k1 = IBCLR(k1, i1);  k2 = IBSET(k2, i2)
    end if
  else
    if( BTEST(k2, i2) )then
      k1 = IBSET(k1, i1);  k2 = IBCLR(k2, i2)
    end if
  end if
  return
end subroutine bitcross

