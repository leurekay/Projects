subroutine checkdyson(N)
  implicit none
  integer N, iseed, ith, ispin, i, j, k, ithrow, jthcol
  complex*16 ratio
  complex*16, dimension (N) :: vec, col, row
  complex*16, dimension (N, N) :: G, Ginv, work
  complex*16  det, det1
  complex*16 ZDET

  iseed=97235

  do j=1, N; do i=1, N
     Ginv(i, j)=cmplx( ran(iseed)-0.5, ran(iseed)-0.5 )
  end do; end do
  G=Ginv; call ZINVERT(N, G); det=ZDET(N, Ginv)

  do k=1, 100
  do ith=1, N
     row = cmplx( ran(iseed)-0.5, ran(iseed)-0.5 )
	 col = cmplx( ran(iseed)-0.5, ran(iseed)-0.5 )
     ithrow=int(N*ran(iseed))+1; jthcol=int(N*ran(iseed))+1
	 call rowcol(N, ithrow, row, jthcol, col, Ginv, G, ratio, 1)
	 Det1=ZDET(N, Ginv)
	 if(abs(det1/det-ratio)>1.d-6)then
	   print*, ratio, det1/det
	   stop 'rowcol failed'
	 end if
	 Det = Det1
	 work=matmul(G, Ginv)
	 do i=1, N; do j=1, N
	    if(i==j)then; if(abs(work(i, i)-1)>1.d-7)print*, work(i, i); end if
		if(i/=j)then; if(abs(work(i, j))>1.d-7)print*, work(i, j); end if
     end do; end do
  end do
  end do

100 continue

  do ith=1, N; do ispin=1, 2
	 do i=1, N; vec(i) = cmplx( ran(iseed)-0.5, ran(iseed)-0.5 ); end do
     !call dyson(N, Ginv, G, ith, vec, ispin, ratio, 0)
	 call dyson(N, Ginv, G, ith, vec, ispin, ratio, 1)
	 Det1=ZDET(N, Ginv)
	 if(abs(det1/det-ratio)>1.d-6)then
	   print*, ratio, det1/det
	   stop 'dyson failed'
	 end if
	 det=det1
  end do; end do
  print*, 'dyson ok'
  return
end subroutine checkdyson

subroutine dyson(N, Ginv, G, ith, vec, rowcol, ratio, info)
  implicit none
  integer N, ith, rowcol, info
  complex*16 ratio, vec(N), Ginv(N, N), G(N, N)
  complex*16 work(N, N), Gii

  logical iwrite/.false./
  logical icheck/.true./

  !if rowcol=1, the ith-row of Ginv => vec
  !if rowcol=2, the ith-column of Ginv => vec
  !ratio = det(new Ginv) / det(old Ginv)

  !options:
  !  if info=0, only ratio is calculated
  !  if info=1, Ginv and G are also updated

  if(rowcol==1)then
    ratio = sum( vec*G(:, ith) )
  else
    ratio = sum( G(ith, :)*vec )
  end if

  if(info==0)return

  !update Ginv
  if(rowcol==1)then; Ginv(ith, :)=vec; else; Ginv(:, ith)=vec; end if
  
  !prepare auxiliary green's func: work
  Gii=G(ith, ith)

  !in case Gii=0 or ratio is too small, G is renewed de fort
  if(abs(Gii)<1.d-10 .or. abs(ratio)<1.d-10)then
    G = Ginv; call ZINVERT(N, G); if(iwrite)print*, 'band G or Ginv encountered @ dyson'; return
  end if

  work = G - matmul( G(:, ith:ith), G(ith:ith, :) ) / Gii
  !notice that work(ith, :)=work(:, ith)=0

  !save corrected G(ith, ith)
  Gii = Gii/ratio

  !calculate corrected G(:, ith) and G(ith, :) using work and Ginv
  G(:, ith:ith) = - matmul(work, Ginv(:, ith:ith))*Gii
  G(ith:ith, :) = - matmul(Ginv(ith:ith, :), work)*Gii

  !recover G(ith, ith), which was erased by work in the above 
  G(ith, ith) = Gii

  !get all elements of corrected G
  G = matmul(G(:, ith:ith), G(ith:ith, :))/Gii + work     
  
  if(icheck)call checkG(N, Ginv, G)
  return
end subroutine dyson

subroutine rowcol(N, ithrow, row, jthcol, col, Ginv, G, ratio, info)
  implicit none
  integer N, ithrow, jthcol, info
  complex*16 ratio
  complex*16, dimension (N) :: row, col
  complex*16, dimension (N, N) :: Ginv, G

  integer, dimension (N) :: maprow, mapcol
  integer i, j, i0, j0, irow, jcol, icol, jrow
  complex*16, dimension (N) :: G0i, Gi0
  complex*16 Mij, Mji, G00

  logical iwrite/.false./
  logical icheck/.true./

  !Ginv(ithrow, :) is to be replaced by row(:)
  !Ginv(:, jthcol) is to be replaced by col(:)
  !it is assumed that col overwrites row at the crossing point (ithrow, jthcol)
  !ratio = det(new Ginv) / det(old Ginv)

  !options:
  !  if info=0, only ratio is calculated
  !  if info=1, Ginv and G are also updated

  if(ithrow>N .or. ithrow<1)stop 'invalid row label'
  if(jthcol>N .or. jthcol<1)stop 'invalid col label'

  !row labels are mapped so that ithrow --> first row
  maprow(1) = ithrow 
  do i=1, ithrow-1; maprow(i+1) = i; end do
  do i=ithrow+1, N; maprow(i) = i; end do

  !col labels are mapped so that jthcol --> first column
  mapcol(1) = jthcol
  do i=1, jthcol-1; mapcol(i+1) = i; end do
  do i=jthcol+1, N; mapcol(i) = i; end do

  !(row, col) mapped to (0, 0)
  i0 = maprow(1); j0 = mapcol(1)

  !notice: row and col obey different mapping in G and Ginv:
  !for G(i, j), i <--> mapcol(i),  j <--> maprow(j), 
  !for Ginv(i, j), i <--> maprow(i), j <--> mapcol(j).

  ratio = col(ithrow) * G(j0, i0)
  !assuming col overwrites row at the crossing position
   
  do j = 2, N;  jrow = maprow(j)
  do i = 2, N;  icol = mapcol(i)
     Mij =  G(icol, jrow)*G(j0, i0) - G(icol, i0) * G(j0, jrow)
     ratio = ratio - row(icol) * Mij * col(jrow)
  end do; end do

  if(info == 0)return

  !update Ginv
  Ginv(ithrow, :) = row;  Ginv(:, jthcol) = col

  if(abs(ratio) < 1.d-10 .or. abs(G(j0, i0)) < 1.d-10 )then
    G = Ginv; call ZINVERT(N, G); if(iwrite)print*, 'bad G or Ginv encountered @ rowcol'; return
  end if
  
  G00 = G(j0, i0) / ratio
  Gi0 = G(:, i0); G0i = G(j0, :)        

  !update G(:, i0) and G(i0, :)
  do j = 2, N; jrow = maprow(j); jcol = mapcol(j)  
     G(jcol, i0) = 0; G(j0, jrow) = 0
     do i = 2, N;  irow = maprow(i); icol = mapcol(i)
	    Mij = G(icol, jrow)*G(j0, i0) - Gi0(icol) * G0i(jrow) 
	    Mji = G(jcol, irow)*G(j0, i0) - Gi0(jcol) * G0i(irow) 
	    G(jcol, i0) = G(jcol, i0) - Mji*col(irow)/ratio
	    G(j0, jrow) = G(j0, jrow) - row(icol)*Mij/ratio
	 end do
  end do

  !update G(i, j) for i/=j0 and j/=i0
  do j = 2, N;   jrow = maprow(j)
  do i = 2, N;   icol = mapcol(i)
     Mij = G(icol, jrow)*G(j0, i0) - Gi0(icol) * G0i(jrow)      
	 G(icol, jrow) = Mij/G(j0, i0) + G(icol, i0) * G(j0, jrow)/G00
  end do; end do

  G(j0, i0) = G00

  if(icheck)call checkG(N, Ginv, G)
  return
end subroutine rowcol


subroutine checkG(N, Ginv, G)
  integer N
  complex*16, dimension (N, N) :: Ginv, G

  real*8, parameter :: tolerance = 1.d-4
  real*8 error
  integer i
  logical quiet/.true./

  error = 0
  do i = 1, N; error = max( error, abs( sum( Ginv(i, :) * G(:, i) ) -1 ) ); end do

  if(error > tolerance)then
    if(.not. quiet)print*, 'G * Ginv error =', error
    G = Ginv; call ZINVERT(N, G)
  else
    if(.not. quiet)print*, 'G is fine'
  end if

  return
end subroutine checkG



subroutine brute1vec(N, Ginv, ith, vec, rowcol, det, ratio, info)
  implicit none
  integer N, ith, rowcol, info
  complex*16 det, ratio, vec(N), Ginv(N, N)
  complex*16 work(N, N)
  complex*16, external :: ZDET

  if(info /= 0)then
    det = det * ratio
	if(rowcol == 1)then; Ginv(ith, :) = vec; else; Ginv(:, ith) = vec; end if
  else
    work = Ginv
	if(rowcol == 1)then; work(ith, :) = vec; else; work(:, ith) = vec; end if
	ratio = ZDET(N, work) / det
  end if

  return
end subroutine brute1vec


subroutine brute1r1c(N, ithrow, row, jthcol, col, Ginv, det, ratio, info)
  implicit none
  integer N, ithrow, jthcol, info
  complex*16 det, ratio
  complex*16, dimension (N) :: row, col
  complex*16, dimension (N, N) :: Ginv, work
  complex*16, external :: ZDET

  if(info /= 0)then
    det = det * ratio
	Ginv(ithrow, :) = row;  Ginv(:, jthcol) = col
  else
    work = Ginv;  work(ithrow, :) = row;  work(:, jthcol) = col
	ratio = ZDET(n, work) / det
  end if

  return
end subroutine brute1r1c


SUBROUTINE BRUTEnVECS(N, GINV, NV, V, LV, COL, DET, RATIO, INFO)
  implicit none
  INTEGER, INTENT(IN) :: N, NV, INFO
  INTEGER, INTENT(IN) :: LV(NV)          
  LOGICAL, INTENT(IN) :: COL(NV)         
  COMPLEX*16, INTENT(IN) :: V(N, NV)
  COMPLEX*16 DET, RATIO
  COMPLEX*16, DIMENSION(N, N) :: GINV, WORK

  INTEGER I
  COMPLEX*16, EXTERNAL :: ZDET
  
  IF( INFO /= 0 )THEN
    DO I = 1, NV
	   IF( COL(I) )THEN; GINV(:, LV(I) ) = V(:, I); ELSE; GINV( LV(I), :) = V(:, I); END IF
	END DO
	DET = DET * RATIO
  ELSE
    WORK = GINV
    DO I = 1, NV
	   IF( COL(I) )THEN; WORK(:, LV(I) ) = V(:, I); ELSE; WORK( LV(I), :) = V(:, I); END IF
	END DO
	RATIO = ZDET(N, WORK) / DET
  END IF
  RETURN
END SUBROUTINE BRUTEnVECS






 