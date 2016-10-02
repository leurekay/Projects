program Fcode_cn
  integer::i=0, j=1, k=2, m(2)
  !interface  ! ½Ó¿Ú¿é
   ! subroutine sub1(i,j,k)
    !  integer,optional::k
     ! integer i,j
    !end subroutine
    !function func1(j,k)
     ! integer j, k
      !integer func1(2)
    !end function
  !end interface
  m = func1(j,k)
  print*, m  ! m=3,-1
  call sub1(i,j,k)
  print*, i  ! i=3
  call sub1(i,j)
  print*, i  ! i=1
  pause
end

function func1(j,k)
  integer j, k
  integer func1(2)
  func1(1) = j + k
  func1(2) = j - k
end function
subroutine sub1(i,j,k)
  integer,optional::k
  integer i,j
  if( present(k) ) then
    i = j + k
  else
    i = j
  end if
end subroutine
