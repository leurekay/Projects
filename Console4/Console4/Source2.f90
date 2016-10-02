subroutine b(o)
    implicit none
    real(kind=8),intent(out)::o
    real(kind=8),dimension(3)::x,y
    x=(/1.,2.,3./)
    y=(/4.,5.,6./)
    o=dot_product(x,y)
    
end subroutine b