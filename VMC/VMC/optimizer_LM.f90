subroutine LMoptimizer( )
  use vmcplace
  implicit none

  real*8, allocatable, dimension (:) :: dvar, step

  integer nter, i, iv, iter, jseed
  logical iwrite/.true./, irandom/.false./
  logical converged/.false./, auto_stop/.false./, switch_method/.false./
  real*8 scale, Eavold

  integer, external :: milisecond
  real*8, external :: rand

  !conjugate must be specified by user in job.txt
  newton = .not. conjugate

  LMoptimize = 1;   optimize = 1;   optinfo = 0;  call definejob( )
  if(nvar > 100)stop 'too many vars to be optimized by LMoptimizer'

  LM = nvar + 1
  allocate( Edlog(nvar), Avdlog(nvar), dvar(nvar), step(nvar), Svv(LM, LM), Evv(LM, LM) )
  !all relevant matrices are allocated to enable method switching

  jseed = milisecond()
  call MPI_barrier(mpi_comm_world, ierr)
  call broadcast_integer(jseed)

  call initVar( )

  if(ithread == 0) open(17, file = 'histo.'//model)

  nter = 100;  Eavold = 1.d10
  do iter = 1, nter
    call MPI_barrier(mpi_comm_world, ierr)

    !switch to conjugate after some steps of newton 
    if(iter > 10) switch_method = newton

    optinfo = optinfo + 1;  call definejob;     call infrastructure( )
    call maketables( ); 	call MCThread( );   call MPI_average( )

    if(Eav > Eavold .and. iter > 1) switch_method = newton
    converged = ( abs(Eav - Eavold) < delE ) 

    if(ithread == 0) then
	   if(nvar <= 6) then; write(17, 100) var(1:nvar), Eav, delE; else; write(17, 100)Eav, delE; end if
       if(Eav < Eavold)then;  Eavold = Eav; call saveVar( );  call normal_band(unitcell);  end if
       if(iwrite)then
         do i = 1, nvar; print*, vname(i), ' = ', var(i);  end do
         print*, 'Eav = ', Eav;  print*, 'delE = ', delE
         print*, 'rho = ', rhoav;  print*, 'delrho = ', delrho;  print*, ' '
	   end if
    else
       if(Eav < Eavold)Eavold = Eav
       !make sure this is done on 0th and all of the other nodes since Eold is used later to judge 
       !convergence on each node. (otherwise convergence must be broadcasted instead.)
    end if

    if(converged .and. auto_stop)exit

	if(newton)then
	  dvar = - 2 * real( Edlog - Avdlog * Eav ) * 100
	  do i = 1, nvar;  if( abs( dvar(i) ) > step(i) ) dvar = dvar * step(i) / abs( dvar(i) );  end do
      scale = rand(jseed);  dvar = dvar * scale
	else
	  call conjugate_search(dvar, jseed)
	end if

	var(1:nvar) = var(1:nvar) + dvar 
    where(var(1:nvar) > vmax(1:nvar)) var(1:nvar) = vmax(1:nvar)
    where(var(1:nvar) < vmin(1:nvar)) var(1:nvar) = vmin(1:nvar)

    if(.not. switch_method) cycle

    conjugate = .true.;  newton = .false.
    irandom = .false.;   call initVar( )        !read optimized vars to start with
  end do
  if(ithread == 0) close(17)

100 format(1x, 8f12.5)
200 format(1x, e20.10)

  return
  contains
      
	 subroutine initVar( )
       implicit none   
       character*6 name; character*1 equal;  integer inode
       do inode = 0, nthread - 1
          call MPI_barrier(mpi_comm_world, ierr);  if(inode /= ithread)cycle
          open(16, file = 'optimized.'//model)   
          do i = 1, nvar
             step(i) = 0.1 * ( vmax(i) - vmin(i) )
             var(i) = vmin(i) + ( vmax(i) - vmin(i) ) * 0.5 * rand(jseed)
	         if(.not. eof(16) )then
	           read(16, *) name, equal, var(i)
	           if( name /= vname(i) )stop 'vname mismatch in optimized input file'
			   if(irandom)var(i) = var(i) + 0.02 * (rand(jseed)-0.5)
            end if
       end do;  close(16); end do
	   return
	end subroutine initVar

    subroutine saveVar( )
	  implicit none
      open(16, file = 'optimized.'//model)
      do i = 1, nvar;  write(16, *)vname(i), ' = ', var(i);  end do;  close(16)

      if( allocated(sdwav) )then
        open(16, file = 'sdw.'//model)
        do i = 1, naa;  write(16, 100) sdwav(i), delsdw(i); end do;  close(16)
      end if

      if( allocated(cdwav) )then
        open(16, file = 'cdw.'//model)
        do i = 1, naa;  write(16, 100) cdwav(i), delcdw(i); end do;  close(16)
      end if
100   format(1x, 2f12.6)

      open(16, file = 'rho.'//model)
      write(16, *)'rhoav = ', rhoav;  write(16, *)'delrho = ', delrho; close(16)
      
	  return
	end subroutine saveVar
end subroutine LMoptimizer


subroutine conjugate_search(dVar, jseed)
  use vmcplace
  implicit none  
  integer jseed
  real*8 dVar(nvar)

  integer i, iv, best, imap(LM)
  real*8 scale, vari, stabilizer, val(LM)
  complex*16 lambda
  logical valid
  complex*16, dimension (LM, LM) :: VL, VR, work(LM)
  real*8, external :: rand

  stabilizer = 1.d-3
  do i = 1, LM;  Svv(i, i) = Svv(i, i) + stabilizer;  end do

  VL = Svv;  call ZINVERT(LM, VL);  Evv = matmul(VL, Evv)
  call ZGEIGEN(LM, Evv, VL, VR, work);  val = real(work)
  call ASCENDINGORDER(LM, val, imap)

  best = imap(1);  work = VR(:, best) / VR(1, best)

  lambda = sum( work(2:LM) * Svv(1, 2:LM) ) / Svv(1, 1)

  !lambda determined and used as follows:
  !write the eigen state in the variational hilbert space as
  !          psi = (1 + lambda) psi_0 + (a_i psi_i - lambda psi_0).
  !we want to make the second part orthogonal to the first part via a 
  !parameter lambda,
  !          lambda = sum( a_i Svv(0, i) ) / Svv(0, 0). 
  !if the change is too big, one may want to scale down the second part
  !    psi => (1 + lambda) psi_0 + s (a_i psi_i - lambda psi_0),
  !which is rewritten as
  !    psi => (1 + lambda - lambda s) psi_0 + s (a_i psi_i).
  !this leads to
  !    psi => psi_0 + s (a_i psi_i) / (1 + lambda - lambda s),
  !so that the variational parameters are updated as
  !    x_i => x_i + s a_i / (1 + lambda - lambda s)

  scale = 1.d-1;  valid = .false.
  
  do while (.not. valid)

    dVar = real( work(2:LM) * scale / (1 + lambda - lambda * scale) )

    !truncation in case of overflow from user specified boundaries    
    valid = .true.
    do iv = 1, nvar;  vari = var(iv) + dvar(iv)
       if( vari <= Vmax(iv) .and. vari >= Vmin(iv) )cycle
       valid = .false.;  scale = scale / 2;  exit
	end do

    if(scale < 1.d-8) stop 'boundary of var approached'

  end do

  return
end subroutine conjugate_search


subroutine maketables( )
  use vmcplace
  implicit none

  integer i, nfields
  complex*16, dimension (na1uc, na1uc, nuc) :: twork
  real*8 vwork(nvar), delv
  logical, external :: Jastraw

  nfield = 0;  do i = 1, nvar;  if( Jastraw( vname(i) ) )cycle;  nfield = nfield + 1; end do
  do i = 1, nfield
     if( Jastraw( vname(i) ) )stop 'Jastraw vars should be positioned after mean field vars'  
  end do

  if( nfield == 0 )return

  if(.not. allocated(tables) ) allocate( tables(na1uc, na1uc, nuc, nfield) )
  if(conjugate .and. .not. allocated(Dv) ) allocate( Dv(npair, npair, nfield) )

  twork = table;  vwork = var(1:nfield)
  do i = 1, nfield
     delv = ( vmax(i) - vmin(i) ) * 1.d-5
     var(1:nfield) = vwork;  var(i) = var(i) + delv
	 optinfo = optinfo + 1;	 call definejob( );  call infrastructure( )
	 tables(:, :, :, i) = ( table - twork ) / delv
  end do
  table = twork;  var(1:nfield) = vwork
  
  return
end subroutine maketables


 