subroutine definejob(  )
  use vmcplace
  implicit none

  !default parameters (can be reassigned by user)  
  bruteforce = .false.;  whichtable = 0  
  gU = 0;    Vdh = 0;       Vzz = 0;      Vcc = 0
     
  select case (model)
	case ( 'hexKondo' );   call hexKondo( )
	case ( 'honeycomb' );  call honeycomb( )
    case ( 'kagome' );     call kagome(  )
    case ( 'kagomeJ' );    call kagomeJ(  )
    case ( 'kagome2x2' );  call kagome2x2( )
	case ( 'r2xr2' );      call r2xr2_mannual(  )
	case ( 'r2xr2_ionic' ); call r2xr2_ionic(  )
	case ( 'r2xr2Kondo' ); call r2xr2Kondo( )
    case ('YRZ');          call YRZ( )
    case default;          stop 'unknown model'
  end select

  return
end subroutine definejob
