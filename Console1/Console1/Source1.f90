program loop
implicit none
integer,parameter :: n=100
real (kind=8), dimension(n) :: x,y
integer :: i

do i=1,n
    if (i==1) then
        x(i)=1
    else
        x(i)=i*x(i-1)
    endif
enddo

print *,x(n)
end program loop

