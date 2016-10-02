program array
    implicit none 
    real*8,dimension(3,2)::a
    real*8,dimension(2,3)::b
    b=reshape((/1,2,3,4,5,6/),(/2,3/))
    print *,b
end program array