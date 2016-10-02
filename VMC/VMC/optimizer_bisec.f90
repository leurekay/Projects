subroutine bisec_opt()
  use vmcplace
  implicit none

  integer, parameter :: npixel = 5
  real*8, parameter :: convergence = 1.d-3

  integer i

  LMoptimize = 0;  optimize = 1;  optinfo = 0;  call definejob( )
  !get info of var table

  if( nvar > 5 )stop 'too many vars to be optimized by bisec'

  if(ithread == 0)open(17, file = 'histo.'//model)
  !open output file in main thread only

  call bisec(nvar, vmin, vmax, npixel, convergence, energy)
  
  if(ithread /= 0)return
  
  close(17)

  open(11, file = 'optimized.'//model)
  do i = 1, nvar;  write(11, *)vname(i), ' = ', var(i);  end do
  close(11)

  return
  contains 

  function Energy(nv, v)
    implicit none

    integer nv
    real*8 Energy, v(nv)

    call MPI_barrier(mpi_comm_world, ierr)

    !transport v to var to be visible in vmcplace
    var(1:nv) = v

    optinfo = optinfo + 1
    call definejob( );  call infrastructure( )
    call MCthread( );   call MPI_average( )

    Energy = Eav

    if(ithread == 0) write(17, 100)v, Eav, delE, rhoav, delrho
100 format(1x, 7f11.5)

    return
  end function Energy

end subroutine bisec_opt









!===========================================================================
!    portable bisec subroutine 

subroutine bisec(nv, vmin, vmax, npixel, convergence, uspec)
  implicit none

  integer, intent(in) :: nv, npixel
  real*8, intent(in) :: convergence, vmin(nv), vmax(nv)

  !interface to user specified function to be minimized
  interface
    function uspec(nv, v)
	  integer nv
	  real*8 uspec, v(nv)
	end function uspec
  end interface

  real*8, dimension(nv) :: lower, upper, var
  real*8 x, dx, vspan, E, vopt, Eopt
  integer iv, ix
  logical converged

  !scanning step in the affine path (for all vars)
  dx = 1.d0 /(npixel - 1.d0)

  !copy of the initial bounds
  lower = vmin;  upper = vmax

  !set inital vars at the midpoints within their bounds
  var = (lower + upper) / 2

  do while (1)    
     
	 converged = .true.
	 !assume convergence, and see whether it stays
     
	 do iv = 1, nv

        if( abs( upper(iv) - lower(iv) ) < convergence )cycle
		!if the bounds converged along the direction, no scan is performed further along this direction

		converged = .false.      
		!if any direction needs to be scanned and optimized, convergence is not yet reached
          
        Eopt = 1.d10
        !infinite upper bound for E to be optimized

        vspan = upper(iv) - lower(iv)
		!span of var(iv)

	    do ix = 1, npixel;  x = (ix-1) * dx
     	       
		   var(iv) = lower(iv) + vspan * x;   E = uspec(nv, var)
		   !E for var obtained by user specified energy function

		   if(E < Eopt)then; Eopt = E; vopt = var(iv); end if
		   !update Eopt and vopt
        
	    end do
		
		var(iv) = vopt
		!record vopt in the scanning direction
		 	 
        lower(iv) = ( lower(iv) + vopt ) / 2;   upper(iv) = (upper(iv) + vopt) / 2 
	    !renew the bounds of var by bisection

      end do  !end of iv
	  if(converged) exit
  end do      !end of optimization

  return
end subroutine bisec