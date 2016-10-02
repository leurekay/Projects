subroutine vpoint()
  use vmcplace
  implicit none
  integer i

  LMoptimize = 0;  optimize = 0
  
  call definejob(  );    call infrastructure(  )

  if(ithread == 0)then
    call normal_band(unitcell);  call bdg_band(unitcell)
  end if

  call MCthread( );      call MPI_average( )

  if(ithread /= 0)return
  !tell MPI to write in the main thread only

  open(10, file = 'vpoint.'//model)
  write(10, *)'acceptance = ', acceptance
  write(10, *)'Eav, delE =', Eav, delE
  write(10, *)'rhoav = ', rhoav
  write(10, *)'delrho = ', delrho; close(10)

  if( allocated(sdwav) )then
    open(10, file = 'sdw.'//model)
    do i = 1, naa; write(10, 100)sdwav(i), delsdw(i); end do; close(10)
  end if

  if( allocated(cdwav) )then
    open(10, file = 'cdw.'//model)
    do i = 1, naa; write(10, 100)cdwav(i), delcdw(i); end do; close(10)
  end if

100 format(1x, 6f12.6)    
  return
end subroutine vpoint

subroutine vline()
  use vmcplace
  implicit none

  integer i
  real*8 x
  logical iwrite/.true./

  LMoptimize = 0;  optimize = 1;  optinfo = 0;  call definejob( )
  
  if( nvar /= 1 )stop 'line scan not well defined'

  if(ithread == 0)then
    open(11, file = 'vline.'//model)
    open(12, file = 'sdw.'//model)
    open(13, file = 'cdw.'//model)
  end if

  do x = 0.d0, 1.d0, 5.d-2
     call MPI_barrier(mpi_comm_world, ierr)

     var(1) = vmin(1) + ( vmax(1) - vmin(1) ) * x

	 optinfo = optinfo + 1
	 
	 call definejob(  );     call infrastructure(  )
	 call MCthread( );  	 call MPI_average(  )

     if(ithread /= 0) cycle
	 write(11, 100) var(1), Eav, delE, acceptance, rhoav, delrho
	 if(iwrite) write(*, 100) var(1), Eav, delE, acceptance, rhoav, delrho

     if( allocated(sdwav) )then
       do i = 1, naa; write(12, 100) sdwav(i), delsdw(i); end do
     end if

     if( allocated(cdwav) )then
       do i = 1, naa;  write(13, 100) cdwav(i), delcdw(i); end do
     end if

  end do
100 format(1x, 8f12.6)

  if(ithread == 0)then
    close(11);  close(12);  close(13)
  end if

  return
end subroutine vline


subroutine MPI_average( )
  use vmcplace
  implicit none

  if(nthread == 1)return

  delE = Eav * Eav;  delrho = rhoav * rhoav
  if( allocated(sdwav) ) delsdw = sdwav * sdwav
  if( allocated(cdwav) ) delcdw = cdwav * cdwav

  call MPI_barrier( mpi_comm_world, ierr )  
  
  call reduce_dble( Eav );      Eav = Eav / nthread
  call reduce_dble( delE );     delE = delE / nthread
  delE = delE - Eav * Eav;      if(delE > 0)delE = sqrt( delE /(nthread-1) )

  call reduce_dble_vector(2, rhoav);    rhoav = rhoav / nthread
  call reduce_dble_vector(2, delrho);   delrho = delrho / nthread
  delrho = delrho - rhoav * rhoav
  where(delrho > 0) delrho = sqrt( delrho / (nthread - 1) )

  if( allocated(sdwav) )then
    call reduce_dble_vector(naa, sdwav);  sdwav = sdwav / nthread
    call reduce_dble_vector(naa, delsdw); delsdw = delsdw / nthread - sdwav * sdwav
    where(delsdw > 0) delsdw = sqrt( delsdw / (nthread -1) )
  end if

  if( allocated(cdwav) )then
    call reduce_dble_vector(naa, cdwav);  cdwav = cdwav / nthread
    call reduce_dble_vector(naa, delcdw);  delcdw = delcdw / nthread - cdwav * cdwav
    where(delcdw > 0) delcdw = sqrt( delcdw / (nthread -1) )
  end if

  if(LMoptimize .and. newton)then
    call reduce_cmplx_vector(nvar, Edlog);  Edlog = Edlog / nthread
    call reduce_cmplx_vector(nvar, Avdlog); Avdlog = Avdlog / nthread
  end if
  
  if(LMoptimize .and. conjugate)then
    call reduce_cmplx_matrix(LM, LM, Evv);  Evv = Evv / nthread
    call reduce_cmplx_matrix(LM, LM, Svv);  Svv = Svv / nthread
  end if

  return
end subroutine MPI_average






























subroutine MPI_average_obsolete( )
  use vmcplace
  implicit none

  real*8 Eavs(nthread)
  complex*16, dimension(nvar, nthread) :: Edlogs, Avdlogs
  complex*16, dimension(LM, LM, nthread) :: Evvs, Svvs
  integer i, nelement

  integer status( MPI_status_size )

  !==============================================================
  !                       collect data 

  call MPI_barrier( mpi_comm_world, ierr )            
  !tell MPI to wait for data communication

  if(LMoptimize .and. conjugate) nelement = LM * LM

  if(ithread /= 0)then 

    !if ithread is not the main thread, ask MPI to send data in other threads to the main thread
    call MPI_send( Eav, 1, mpi_double_precision, 0, 88, mpi_comm_world, ierr )
    
	if(LMoptimize)then
	  if(newton)then
	    call MPI_send( Edlog, nvar, mpi_double_complex, 0, 68, mpi_comm_world, ierr )
        call MPI_send( Avdlog, nvar, mpi_double_complex, 0, 78, mpi_comm_world, ierr )
      else
	    call MPI_send( Evv, nelement, mpi_double_complex, 0, 68, mpi_comm_world, ierr )
        call MPI_send( Svv, nelement, mpi_double_complex, 0, 78, mpi_comm_world, ierr )
	  end if
	end if
  else

    !if ithread is the main thread, collect data in all threads and analyze data

	!data in the main thread
    Eavs(1) = Eav              
	if(LMoptimize .and. newton)then;  Edlogs(:, 1) = Edlog;  Avdlogs(:, 1) = Avdlog; end if
	if(LMoptimize .and. conjugate)then;  Evvs(:, :, 1) = Evv;  Svvs(:, :, 1) = Svv; end if

    !receive data in other i-thread's and save as dat(i+1)
    do i = 1, nthread - 1
	   call MPI_recv( Eavs(i+1), 1, mpi_double_precision, i, 88, mpi_comm_world, status, ierr )
	   if(.not.LMoptimize)cycle
	   if(newton)then
	     call MPI_recv( Edlogs(:, i+1), nvar, mpi_double_complex, i, 68, mpi_comm_world, ierr )
         call MPI_recv( Avdlogs(:, i+1), nvar, mpi_double_complex, i, 78, mpi_comm_world, ierr )
       else
	     call MPI_recv( Evvs(:, :, i+1), nelement, mpi_double_complex, i, 68, mpi_comm_world, ierr )
         call MPI_recv( Svvs(:, :, i+1), nelement, mpi_double_complex, i, 78, mpi_comm_world, ierr )
       end if
    end do
	
	!average and standard error over the threads
	if(nthread > 1)then
	  Eav = sum( Eavs ) / nthread
      delE = sum( Eavs**2 ) / nthread - Eav**2
	  if(delE > 0) delE = sqrt(delE) 
	  if(nthread > 2)delE = delE/sqrt(nthread-1.d0)

	  if(LMoptimize)then
	     if(newton)then
	        Edlog = 0;  do i = 1, nthread;  Edlog = Edlog + Edlogs(:, i);  end do
		    Avdlog = 0;  do i = 1, nthread;  Avdlog = Avdlog + Avdlogs(:, i);  end do
		    Edlog = Edlog / nthread;  Avdlog = Avdlog / nthread
	     else
		    Evv = 0;  do i = 1, nthread;  Evv = Evv + Evvs(:, :, i); end do
			Svv = 0;  do i = 1, nthread;  Svv = Svv + Svvs(:, :, i); end do
			Evv = Evv / nthread;   Svv = Svv / nthread
		 end if
	  end if
	end if

  end if
  !=====================================================================


  !=====================================================================
  !              publish Eav in the main thread to other threads

  call MPI_barrier( mpi_comm_world, ierr )
  !tell MPI to wait for data communication again

  if(ithread /= 0)then 

    !if ithread is not the main thread, let it receive data from the main thread
    call MPI_recv( Eav, 1, mpi_double_precision, 0, 89, mpi_comm_world, status, ierr )
    
	if(LMoptimize)then
	  if(newton)then
	    call MPI_recv( Edlog, nvar, mpi_double_complex, 0, 69, mpi_comm_world, status, ierr )
        call MPI_recv( Avdlog, nvar, mpi_double_complex, 0, 79, mpi_comm_world, status, ierr )
      else
	    call MPI_recv( Evv, nelement, mpi_double_complex, 0, 69, mpi_comm_world, status, ierr )
        call MPI_recv( Svv, nelement, mpi_double_complex, 0, 79, mpi_comm_world, status, ierr )
	  end if
	end if
  else

    !if ithread is the main thread, send its data to other threads
    do i = 1, nthread - 1
	  call mpi_send( Eav, 1, mpi_double_precision, i, 89, mpi_comm_world, ierr )
      if(.not.LMoptimize)cycle
	  if(newton)then
	    call mpi_send( Edlog, nvar, mpi_double_complex, i, 69, mpi_comm_world, ierr )
	    call mpi_send( Avdlog, nvar, mpi_double_complex, i, 79, mpi_comm_world, ierr )
	  else
	    call mpi_send( Evv, nelement, mpi_double_complex, i, 69, mpi_comm_world, ierr )
	    call mpi_send( Svv, nelement, mpi_double_complex, i, 79, mpi_comm_world, ierr )
      end if
	end do
  end if
  !==========================================================================


  call MPI_barrier( mpi_comm_world, ierr )  
  !tell MPI to wait 

  return
end subroutine MPI_average_obsolete