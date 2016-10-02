program a
    implicit none
    real(kind=8) :: n,m
	real(kind=8)::o
    n=4
    call f(n,m)
    print *,"m=",m
    call b(o)
    print *,"dot_product=",o
end program a
subroutine f(x,y)
    implicit none
    real(kind=8),intent(in) ::x
    real(kind=8), intent(out)::y
    y=x*x
end subroutine f

    
    
