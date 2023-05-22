module util

implicit none

real(kind=8),parameter  :: pi=3.141592653589793238d0
complex(kind=8),parameter :: i_imag=cmplx(0,1,kind=8)
real(kind=8),parameter :: eleVolt = 1.609d-19
real(kind=8),parameter :: mass_proton = 1.6726219d-27
real(kind=8),parameter :: kb = 1.380649d-23
real(kind=8),parameter :: hbar = 1.0545718d-34
 

contains
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  complex(kind=8), dimension(:,:), intent(in) :: A
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Ainv

  complex(kind=8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: i, n, info, tag

  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  tag = ifnan(A)
  if (tag .eq. 1) then
      stop
  end if

  Ainv(:,:) = A(:,:)
  n = size(A,1)
  !do i = 1,n
  !   Ainv(i,i) = Ainv(i,i) + 1.0d-16
  !end do

  ! ZGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     ! = 0:  successful exit
     !     < 0:  if INFO = -i, the i-th argument had an illegal value
     !     > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
     !           has been completed, but the factor U is exactly
     !           singular, and division by zero will occur if it is used
     !           to solve a system of equations.
      write(*,*) 'info'
      write(*,*) info
      write(*,*) Ainv(info,info)
      write(*,*) 'A'
      write(*,*) A(1,:)
     write(*,*) 'Dimension of matrix:', size(A,1)
     stop 'Matrix is numerically singular!'
  end if

  ! ZGETRI computes the inverse of a matrix using the LU factorization
  ! computed by ZGETRF.
  call ZGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  end function inv

  function psinv(A) result(Ainv)

  complex(kind=8), dimension(:,:), intent(in) :: A
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Ainv
  integer(kind=4) :: LDA, LDU, LDVT
  integer(kind=4) :: LWMAX
  integer(kind=4) :: INFO,LWORK
  complex(kind=8) :: S(size(A,1)), Rwork(5*size(A,1))
  complex(kind=8) ::U(size(A,1),size(A,1)), VT(size(A,2),size(A,2)),&
                    WORK(1000)
  complex(kind=8) :: SS(size(A,1))
  integer(kind=4) :: m,n,i,incx
  complex(kind=8) :: alpha=1.0d0
  complex(kind=8) :: beta=0.0d0
 
  external ZGESVD
  external zgemm
  external zscal

 
  incx = 1 ! INCX is INTEGER storage spacing between elements of ZX
  lwmax = 1000

  m = size(A,1)
  n = size(A,2)
  LDA = m 
  LDU = m 
  LDVT = n

  LWORK = -1
  CALL ZGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,&
             WORK, LWORK, RWORK, INFO )
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

!     Compute SVD.
  CALL ZGESVD( 'A', 'A', M, N, A, LDA, S, U, LDU, VT, LDVT,&
             WORK, LWORK, RWORK, INFO )
  
  SS = 0.0d0
  do i = 1,m
      if (abs(s(i)) > 1.0d-17) then
          SS(i) = 1.0d0/s(i) 
          call zscal(m,ss(i),u(i,:),incx)
      end if
  end do

  call zgemm('T','T',n,m,m,alpha,vt,ldvt,&
        u,ldu,beta,Ainv,n)

!     Check for convergence.

  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm computing SVD failed to converge.'
!     STOP
  END IF
  end function psinv

  ! matrix inversion for real matrix
  function inv_real(A) result(Ainv)
  real(kind=8), dimension(:,:), intent(in) :: A
  real(kind=8), dimension(size(A,1),size(A,2)) :: Ainv

  real(kind=8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  end function inv_real

  function diagmat(eigs) result(mat)
  complex(kind=8),intent(in) :: eigs(:)
  complex(kind=8) :: mat(size(eigs,1),size(eigs,1))
  integer(kind=4) :: i

  mat = 0.0d0
  do i = 1,size(eigs,1)
      mat(i,i) = eigs(i)
  end do

  end function diagmat

  function eyemat(n) result(mat)
  integer(kind=4),intent(in) :: n
  real(kind=8) :: mat(n,n)
  integer(kind=4) :: i
  mat = 0.0d0
  do i = 1, n
    mat(i,i) = 1.0d0
  end do
  end function eyemat

  function eyemat_cmplx(n) result(mat)
  integer(kind=4),intent(in) :: n
  complex(kind=8) :: mat(n,n)
  integer(kind=4) :: i
  mat = 0.0d0
  do i = 1, n
    mat(i,i) = 1.0d0
  end do
  end function eyemat_cmplx


  subroutine eigenH(A,w,vecr)

  integer(kind=4) :: nwork=1
  integer(kind=4) :: info
  real(kind=8),allocatable :: rwork(:)
  complex(kind=8),allocatable :: work(:)
  complex(kind=8),intent(in) :: A(:,:)
  complex(kind=8),dimension(size(A,1),size(A,2)) :: At
  real(kind=8),dimension(size(A,1)),intent(out) :: w
  complex(kind=8),dimension(size(A,1),size(A,2)),&
               intent(out) :: vecr
  integer(kind=4) :: n

  n = size(A,1)
  At = A
  allocate(work(nwork))
  allocate(rwork(3*n))
    
  ! compute the eigenvalue and eigenvector
  ! of Hermitian matrix
  ! info = 0, the eigenvalues in ascending order
  call zheev("V","U",n,At(:,:),n,w,work,-1,rwork,info)

  if(real(work(1)).gt.nwork) then
        nwork=idnint(2*real(work(1)))
        deallocate(work)
        allocate(work(nwork))
  end if

  call zheev("V","U",n,At(:,:),n,w,work,nwork,rwork,info)
  vecr = At
  end subroutine

  subroutine eigen(A,w,vecr)
  ! compute the eigenvalue and eigenvectors for a
  ! n by n nonsymmetric complex matrix 
  complex(kind=8), dimension(:,:), intent(in) :: A
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Atemp
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Asq
  complex(kind=8), intent(out) :: w(size(A,1))
  complex(kind=8), allocatable :: work(:)
  real(kind=8), allocatable :: rwork(:)
  complex(kind=8), dimension(size(A,1),size(A,2)) :: vecl
  complex(kind=8),dimension(size(A,1),size(A,2)),intent(out) :: vecr
  integer(kind=4) :: lda,  info, nwork = 1

  Atemp = A
  
  if (size(A,1) /= size(A,2)) then
    stop 'Not a square matrix!'
  end if
  lda = max(size(A,1),1)
  allocate(rwork(max(1,2*lda)))
  allocate(work(nwork))
  
  call ZGEEV("N", "V", lda, Atemp, lda, w,&
   vecl, lda, vecr, lda,work,-1,&
   rwork, info)

  if(real(work(1)).gt.nwork) then
        nwork=idnint(2*real(work(1)))
        deallocate(work)
        allocate(work(nwork))
  end if

  call ZGEEV("N", "V", lda, Atemp, lda, w,&
   vecl, lda, vecr, lda,work,nwork,&
   rwork, info) 

  end subroutine eigen

  subroutine eigen_sorting(A,w,vecr)
  ! compute the eigenvalue and eigenvectors for a
  ! n by n nonsymmetric complex matrix 
  complex(kind=8), dimension(:,:), intent(in) :: A
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Atemp
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Asq
  complex(kind=8), intent(out) :: w(size(A,1))
  integer(kind=4) :: sort_idx(size(A,1))
  complex(kind=8), allocatable :: work(:)
  real(kind=8), allocatable :: rwork(:)
  complex(kind=8), dimension(size(A,1),size(A,2)) :: vecl
  complex(kind=8), dimension(size(A,1),size(A,2)) :: vecp
  complex(kind=8),dimension(size(A,1),size(A,2)),intent(out) :: vecr
  integer(kind=4) :: lda,  info, nwork = 1
  integer(kind=4) :: i,nprop

  Atemp = A
  
  if (size(A,1) /= size(A,2)) then
    stop 'Not a square matrix!'
  end if
  lda = max(size(A,1),1)
  allocate(rwork(max(1,2*lda)))
  allocate(work(nwork))
  
  call ZGEEV("N", "V", lda, Atemp, lda, w,&
   vecl, lda, vecr, lda,work,-1,&
   rwork, info)

  if(real(work(1)).gt.nwork) then
        nwork=idnint(2*real(work(1)))
        deallocate(work)
        allocate(work(nwork))
  end if

  call ZGEEV("N", "V", lda, Atemp, lda, w,&
   vecl, lda, vecr, lda,work,nwork,&
   rwork, info) 

   sort_idx = sorting_index_descent(abs(w))
   w = w(sort_idx)         ! order the eigenvalue 
   vecr = vecr(:,sort_idx) ! according to their amplitude

   nprop = 0
   vecp = 0.0d0

   do i = 1,size(A,1)
        if ((abs(abs(w(i))-1.0d0) .gt. 5.0d-3).and.(i.gt.1)) then
           go to 102 
        end if
   end do

   102 nprop = i-1 
   if (nprop .gt. 0) then
       vecp(:,1:nprop) = vecr(:,1:nprop)
       vecr(:,1:nprop) = gram_schmidt(vecp(:,1:nprop))
   end if

  end subroutine eigen_sorting

  LOGICAL FUNCTION DUMMY_SELECT(ARG) 
  IMPLICIT NONE


  COMPLEX(kind=8),intent(in) :: ARG

  DUMMY_SELECT = (ARG .EQ. ARG)
 
  end function DUMMY_SELECT

  function ifnan(A) result(idx)

  complex(kind=8),dimension(:,:),intent(in) :: A
  integer(kind=4) :: idx,m,n,i,j
  real(kind=8) :: infinity

  infinity = huge(8)
  m = size(A,1)
  n = size(A,2)
  idx = 0
  do i = 1,m
      do j = 1,n
          if (isnan(abs(A(i,j)))) then
              write(*,*)'nan', i,j
              idx = 1
              go to 3233
          end if
      end do
  end do

  3233 continue

  end function ifnan


  function sqrtm(A) result(QUQH)
  ! Firstly, computes for an N-by-N complex nonsymmetric
  ! matrix A, the eigenvalues, the Schur form T.
  
  complex(kind=8), dimension(:,:), intent(in) :: A
  complex(kind=8), dimension(size(A,1),size(A,2)) :: S
  complex(kind=8) :: w(size(A,1))
  complex(kind=8), allocatable :: work(:)
  real(kind=8), allocatable :: rwork(:)
  logical,allocatable :: bwork(:)
  complex(kind=8), dimension(size(A,1),size(A,2)) :: q
  complex(kind=8), dimension(size(A,1),size(A,2)) :: U
  complex(kind=8), dimension(size(A,1),size(A,2)) :: QUQH
  complex(kind=8), dimension(size(A,1),size(A,2)) :: QU
  integer(kind=4) :: N,info,nwork = 1
  integer(kind=4) :: i,j,k,sdim, sdiag
  integer(kind=4) :: tag

  
  sdim = 0
  tag = ifnan(A)
  if (tag .eq. 1) then
      stop
  end if

  S(:,:) = A(:,:)
  q = 0.0d0
  U = 0.0d0
  QUQH = 0.0d0
  QU = 0.0d0

  if (size(A,1) /= size(A,2)) then
    stop 'Not a square matrix!'
  end if
  N = max(size(A,1),1)
  allocate(rwork(max(1,N)))
  allocate(bwork(N))
  allocate(work(nwork))

  call ZGEES('V', 'N',DUMMY_SELECT, N, &
   S, N, sdim, w,&
   q, N,work,-1,&
   rwork, bwork, info)

  if(real(work(1)).gt.nwork) then
        nwork=idnint(2*real(work(1)))
        deallocate(work)
        allocate(work(nwork))
  end if

  call ZGEES('V', 'N',DUMMY_SELECT, N, &
  S, N, sdim, w,&
  q, N,work,nwork,&
  rwork, bwork, info)

  do j = 1,N-1              ! set the lower triangle to zero
        do i = j+1,N
           U(i,j) = 0.0d0
        enddo
  enddo

  do i = 1,N                ! set the diagonal elements
        U(i,i) = sqrt(S(i,i))
  enddo
  
  do sdiag = 1,N-1          ! loop over the N-1 super-diagonals
     do i = 1,N-sdiag
         j = i+sdiag
         U(i,j) = S(i,j)
         do k = i+1,j-1
             U(i,j) = U(i,j) - U(i,k)*U(k,j)
         enddo
         if ((U(i,i) + U(j,j)) .ne. cmplx(0.0,0.0,kind=8))then
             U(i,j) = U(i,j) / (U(i,i) + U(j,j))
         else
             U(i,j) = 0.0d0
         end if
     enddo
  enddo
  call ZGEMM('N', 'N', N, N, N, cmplx(1.0,0.0,kind=8), Q, N, U, N, cmplx(0.0,0.0,kind=8), QU, N)
  call ZGEMM('N', 'C', N, N, N, cmplx(1.0,0.0,kind=8), QU, N,Q, N, cmplx(0.0,0.0,kind=8), QUQH, N)
  if(info /= 0) then
       write(*,*) 'S'
       write(*,*) S
       stop
       write(*,*) 'U'
       write(*,*) U
       !write(*,*) U
       !write(*,*) QU
       !write(*,*) QUQH
       write(*,*) "zgemm returned info = ", info
       write(*,*) "the ", info, "-th argument had an illegal value"
       stop
  end if

  end function

  function lstsq(A, b) result(x)
  ! compute least square solution to A x = B for cmplx A, b
  complex(kind=8), intent(in) :: A(:,:), B(:,:)
  complex(kind=8), allocatable :: x(:,:)
  ! LAPACK variables:
  integer :: info, lda, ldb, lwork, m, n, nrhs, rank
  real(kind=8) :: rcond
  real(kind=8),allocatable :: rwork(:) 
  complex(kind=8), allocatable :: work(:), At(:,:), Bt(:,:)
  integer, allocatable :: jpvt(:)

  m = size(A,1)
  n = size(A,2)
  nrhs = size(B,2)
  lda = max(1,m)
  ldb = max(1,m,n)
  allocate(x(n,nrhs), At(m,n), Bt(m,nrhs), &
            jpvt(n), work(1),rwork(2*n))
  rcond = 0.0
  jpvt(:) = 0
  Bt(:,:) = B(:,:)  ! only one right-hand side
  At(:,:) = A(:,:)
  call zgelsy(m, n, nrhs, At, lda, Bt, ldb,&
        jpvt, rcond, rank, work, &
         -1, rwork, info)  ! query optimal workspace size
  lwork = int(real(work(1)))
  deallocate(work)
  allocate(work(lwork))  ! allocate with ideal size

  call zgelsy(m, n, nrhs, At, lda, Bt, ldb,&
         jpvt, rcond, rank, work, &
         lwork, rwork, info)
  if(info /= 0) then
       write(*,*) "zgelsy returned info = ", info
       write(*,*) "the ", info, "-th argument had an illegal value"
  end if
  x(:,:) = Bt(:,:)
  end function lstsq

  function trace(A) result(X)
  complex(kind=8),intent(in) :: A(:,:)
  complex(kind=8) :: X
  integer(kind=4) :: i, n
  
  if (size(A,1) /= size(A,2)) then
    stop 'Not a square matrix!'
  end if
 
  n = size(A,1)
  
  X = 0.0d0

  do i = 1,n
    X = A(i,i) + X
  end do
  end function trace

  subroutine find_degenerate(eigval,start,idx_list,next_start,case_no,&
             deg_all)

  implicit none

  complex(kind=8),dimension(:),intent(in) :: eigval
  complex(kind=8) :: edeg
  integer(kind=4),intent(in) :: start
  integer(kind=4),dimension(:),intent(in) :: deg_all
  integer(kind=4),intent(out) :: next_start
  integer(kind=4),intent(out) :: case_no
  integer(kind=4),allocatable,intent(out) :: idx_list(:)
  integer(kind=4) :: ne, i, j, edeg_idx,&
                     n, itemp
  complex(kind=8),allocatable :: temp(:)
  
  itemp = 0
  next_start = 0
  case_no = 0 ! case number is zero meaning there are no degenerate states
  ne = size(eigval,1)
  outer: do i = start, ne-1
      do j = i+1,ne
          if ((abs(eigval(i)-eigval(j)) .lt. 1.0d-4).and.&
           (abs(abs(eigval(j))-1.0d0) .lt. 1.0d-3)) then 
              case_no = 1
              edeg_idx = j
              itemp = i
              edeg = eigval(edeg_idx)
          end if
          if( case_no .eq. 1) exit outer
      end do
  end do outer

  next_start = itemp+1

  n = 0
  if (ALLOCATED(idx_list)) then
      deallocate(idx_list)
  end if
  allocate(idx_list(1))
  idx_list = 0

  allocate(temp(1))
 
  if (case_no .eq. 1) then
     do i = 1, ne
         if (abs(eigval(i)-edeg) .lt. 1.0d-4) then
             n = n + 1
             if (n.eq.1) then
                 idx_list(1) = i
                 temp(1) = i
             else
                 deallocate(idx_list)
                 allocate(idx_list(n))
                 idx_list(1:n-1) = temp(:)
                 idx_list(n) = i ! append new element

                 deallocate(temp)
                 allocate(temp(n))
                 temp = idx_list
             end if
         end if
     end do
     do i =1,size(idx_list,1)
         do j = 1,size(deg_all,1)
             if (idx_list(i).eq.deg_all(j)) then
                 deallocate(idx_list)
                 allocate(idx_list(1))
                 idx_list(1) = 0
             end if
         end do
     end do 
  end if

  end subroutine find_degenerate

  function in_ws(a_vectors,k_pos) result(output)

    implicit none
   
    real(kind=8),intent(in) :: a_vectors(3,3), k_pos(3)
    real(kind=8) :: dist(125),dist_min,dist_min2
    real(kind=8) :: dist_remain(124)
    integer(kind=4) :: min_idx
    integer(kind=4) :: icnt,i1,i2,i3
    real(kind=8) :: ndiff(3)
    integer(kind=4) :: output

    dist = 0
    icnt = 0
    do i1 =-2,2
        do i2 =-2,2
            do i3 =-2,2
                ndiff(1) = dble(i1)
                ndiff(2) = dble(i2)
                ndiff(3) = dble(i3)

                icnt =icnt+ 1

                dist(icnt) = dot_product(k_pos-matmul(ndiff,a_vectors),&
                                         k_pos-matmul(ndiff,a_vectors)) 
            end do
        end do
    end do
    dist_min = minval(dist)

    ! find the second smallest element
    min_idx = minloc(dist,1)
    if (min_idx .eq. 1) then
        dist_remain = dist(2:125)
    else if (min_idx .eq. 125) then
        dist_remain = dist(1:124)
    else
        dist_remain(1:min_idx-1) = dist(1:min_idx-1)
        dist_remain(min_idx:124) = dist(min_idx+1:125)
    end if
    dist_min2 = minval(dist_remain)


    if (abs(dist(63)-dist_min) .lt. 1d-3) then
        if (abs(dist_min2-dist_min).lt.1.0d-3) then
            output = 2 ! zone boundary of FBZ
        else
            output = 1 ! inside FBZ
        end if
    else
        output = 0
    end if
    end function in_ws 

    function move_in_ws3(k_pos,a_vectors) result(output)

    implicit none
   
    real(kind=8),intent(in) :: a_vectors(3,3), k_pos(3)
    real(kind=8) :: dist(125)
    integer(kind=4) :: i_min
    integer(kind=4) :: icnt,i1,i2,i3
    real(kind=8) :: ndiff(3)
    real(kind=8) :: output(3)
    real(kind=8) :: k_all(125,3)

    dist = 0
    icnt = 0
    k_all(:,:) = 0.0d0

    do i1 =-2,2
        do i2 =-2,2
            do i3 =-2,2
                ndiff(1) = dble(i1)
                ndiff(2) = dble(i2)
                ndiff(3) = dble(i3)

                icnt =icnt+ 1
                k_all(icnt,:) = ndiff
                dist(icnt) = dot_product(k_pos-matmul(ndiff,a_vectors),&
                                         k_pos-matmul(ndiff,a_vectors)) 
            end do
        end do
    end do
    i_min = minloc(dist,1)
    output = k_pos - matmul(k_all(i_min,:),a_vectors)

    end function move_in_ws3 

    function move_in_ws2(k_pos,a_vectors) result(output)

    implicit none
   
    real(kind=8),intent(in) :: a_vectors(3,3), k_pos(3)
    real(kind=8) :: dist(25)
    integer(kind=4) :: i_min
    integer(kind=4) :: icnt,i1,i2,i3
    real(kind=8) :: ndiff(3)
    real(kind=8) :: output(3)
    real(kind=8) :: k_all(25,3)

    dist = 0
    icnt = 0
    k_all(:,:) = 0.0d0

    do i1 =-2,2
        do i2 =-2,2
                ndiff(1) = dble(i1)
                ndiff(2) = dble(i2)
                ndiff(3) = 0.0d0

                icnt =icnt+ 1
                k_all(icnt,:) = ndiff
                dist(icnt) = dot_product(k_pos-matmul(ndiff,a_vectors),&
                                         k_pos-matmul(ndiff,a_vectors)) 
        end do
    end do
    i_min = minloc(dist,1)
    output = k_pos - matmul(k_all(i_min,:),a_vectors)

    end function move_in_ws2

    function get_pos_index_in_pc(pos_sc,pos_pc,recivec) result(idx)

    implicit none

    real(kind=8),intent(in) :: pos_sc(3)
    real(kind=8),intent(in) :: pos_pc(:,:)
    real(kind=8),intent(in) :: recivec(3,3)
    integer(kind=4) :: natm_pc, i, idx
    real(kind=8) :: coord_crys_temp(3)

    natm_pc = size(pos_pc,1)
    do i = 1,natm_pc
        coord_crys_temp(:) = matmul(&
        pos_sc-pos_pc(i,1:3),recivec)
        if ((abs(coord_crys_temp(1)-idnint(coord_crys_temp(1)))<1d-4).and.&
            (abs(coord_crys_temp(2)-idnint(coord_crys_temp(2)))<1d-4).and.&
            (abs(coord_crys_temp(3)-idnint(coord_crys_temp(3)))<1d-4)) then
            idx = i    
        end if  
    end do


    end function get_pos_index_in_pc

    function fourier_transform(Ham,n,nb,G,R) result(Ham_f)

    ! Fourier transform of a block circualr matrix
    
    implicit none

    integer(kind=4),intent(in) :: n
    integer(kind=4),intent(in) :: nb
    real(kind=8),intent(in) :: G(3)
    real(kind=8),intent(in) :: R(3)
    complex(kind=8),intent(in) :: Ham(:,:)
    complex(kind=8) :: Ham_f(size(Ham,1),size(Ham,2))
    integer(kind=4) :: i,j,k,l
    integer(kind=4) :: na, mm

    na = n/nb
    Ham_f = 0.0d0
      
    do i = 1, nb
        do j = 1,nb
            Ham_f((i-1)*na+1:i*na,(i-1)*na+1:i*na) = &
            Ham_f((i-1)*na+1:i*na,(i-1)*na+1:i*na)+&
            Ham(1:na,(j-1)*na+1:j*na)*&
            exp(i_imag*dot_product(&
            G*dble(i-1)/dble(nb),&
            R*dble(j-1)))
        end do
    end do     
    end function fourier_transform

    function decimation(Es_in0,e_in1,alpha_in,&
             beta_in,convergence) result(g)
        
    implicit none

    complex(kind=8),intent(in) :: Es_in0(:,:)
    complex(kind=8),intent(in) :: e_in1(:,:)
    complex(kind=8),intent(in) :: alpha_in(:,:)
    complex(kind=8),intent(in) :: beta_in(:,:)
    real(kind=8),intent(in) :: convergence
    complex(kind=8) :: Es0(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: e1(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: es(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: alpha(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: beta(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: ive(size(Es_in0,1),size(Es_in0,2))
    real(kind=8) :: maxdiff
    complex(kind=8) :: g(size(Es0,1),size(Es0,2))

    Es0 = Es_in0
    e1 = e_in1
    es = e_in1
    alpha = alpha_in
    beta = dconjg(transpose(alpha_in)) ! beta_in
    maxdiff = maxval(abs(alpha))
    do while( maxdiff .gt. convergence)
        ive = inv(Es0-e1)
        e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                matmul(matmul(beta,ive),alpha)
        es = es + matmul(matmul(alpha,ive),beta)
        alpha = matmul(matmul(alpha,ive),alpha)
        beta = matmul(matmul(beta,ive),beta)
        maxdiff = maxval(abs(alpha))
    end do
    g = inv(Es0-es)
    
    end function decimation

    function gram_schmidt(U1) result(U2)

    implicit none
    complex(kind=8),intent(in) :: U1(:,:)
    complex(kind=8),dimension(size(U1,1),size(U1,2)) :: U2
    integer(kind=4) :: i,j,n


    U2 = U1
    n = size(U1,2)
    do i = 2,n
        do j=1,i-1
            U2(:,i) = U2(:,i) - &
             dot_product(U1(:,i),U2(:,j))/&
             dot_product(U2(:,j),U2(:,j))*&
             U2(:,j)
        end do
    end do
    
    do i = 1,n
        U2(:,i) = U2(:,i)/sqrt(dot_product(U2(:,i),U2(:,i)))
    end do

    end function
    
    function log_fbz(lambda,az) result(kz)

    implicit none
    complex(kind=8),intent(in) :: lambda
    real(kind=8),intent(in) :: az
    real(kind=8) :: kz, Gz

    Gz = 2*pi/az

    kz = modulo(real(log(lambda)/(i_imag*az)),Gz)
    if (kz.gt.Gz*0.5)then
        kz = kz - Gz ! make sure kz is inside fbz 
    end if
                               
    end function log_fbz

    function sorting(A) result(B)

    implicit none

    real(kind=8) :: A(:)
    real(kind=8) :: B(size(A,1))
    integer(kind=4) :: irow, nsize, krow
    real(kind=8) :: buf

    nsize = size(A,1)
    B = A
    do irow = 1, nsize
        krow = minloc(B(irow:nsize), dim=1) + irow - 1
        buf  = B(irow)
        B(irow) = B(krow)
        B(krow) = buf
    enddo
    end function sorting

    function sorting_index(A) result(I)

    implicit none

    real(kind=8) :: A(:)
    integer(kind=4) :: I(size(A,1))
    real(kind=8) :: B(size(A,1))
    real(kind=8) :: buf
    integer(kind=8) :: idx
    integer(kind=4) :: irow, nsize, krow

    nsize = size(A,1)
    B = A
    do irow = 1, nsize
        I(irow) = irow
    end do
    do irow = 1, nsize
        krow = minloc(B(irow:nsize), dim=1) + irow - 1
        buf  = B(irow)
        B(irow) = B(krow)
        B(krow) = buf
        idx  = I(irow)
        I(irow) = I(krow)
        I(krow) = idx
    enddo

    end function sorting_index

    function sorting_index_descent(A) result(I)

    implicit none

    real(kind=8) :: A(:)
    integer(kind=4) :: I(size(A,1))

    I = sorting_index(A)
    I = I(size(A,1):1:-1)

    end function sorting_index_descent


    function rotate_z(pos,ang) result(pos1)

    implicit none

    real(kind=8),intent(in)  :: pos(:,:)
    real(kind=8),intent(in)  :: ang
    integer(kind=4) :: nv
    real(kind=8)  :: pos1(size(pos,1),size(pos,2))
    integer(kind=4) :: i
 
    nv = size(pos,1)

    do i = 1,nv
        pos1(i,1) = pos(i,1)*cos(ang)-pos(i,2)*sin(ang)
        pos1(i,2) = pos(i,1)*sin(ang)+pos(i,2)*cos(ang)
        pos1(i,3) = pos(i,3)
    end do

    end function

    function rotate_z1(pos,ang) result(pos1)

    implicit none

    real(kind=8),intent(in)  :: pos(:)
    real(kind=8),intent(in)  :: ang
    integer(kind=4) :: nv
    real(kind=8)  :: pos1(size(pos,1))
    integer(kind=4) :: i
 
    nv = size(pos,1)

    pos1(1) = pos(1)*cos(ang)-pos(2)*sin(ang)
    pos1(2) = pos(1)*sin(ang)+pos(2)*cos(ang)
    pos1(3) = pos(3)

    end function


    function hermitian(A) result(B)

    implicit none

    complex(kind=8) :: A(:,:)
    complex(kind=8) :: B(size(A,1),size(A,2))

    B = 0.5*(A+transpose(dconjg(A)))

    end function hermitian

    function my_sqrt(M) result(N)

    implicit none

    real(kind=8) :: M(:)
    real(kind=8) :: N(size(M,1))
    integer(kind=4) :: i

    N = 0.0d0
    do i = 1,size(M,1)
       if (M(i) .lt.0.0d0) then
           N(i) = -sqrt(-M(i))
       else
           N(i) = sqrt(M(i))
       end if    
    end do
    end function my_sqrt

    function find_min_dist(kpoint,n) result(kmin) 

    implicit none

    real(kind=8) :: kpoint(:,:)
    integer(kind=4) :: n,i
    real(kind=8) :: kabs(n)
    real(kind=8) :: kmin(3)
    integer(kind=4) :: minid

    do i = 1,n
        kabs(i) = sqrt(dot_product(kpoint(i,:),kpoint(i,:)))
    end do

    minid = minloc(kabs,1)
    kmin = kpoint(minid,:)
    end function

    function bose_einstein(freq,temperature) result(fbe)

    implicit none

    real(kind=8) :: fbe
    real(kind=8) :: freq,temperature
    real(kind=8) :: expf

    expf = exp(hbar*freq*1.0d12*2*pi/kb/temperature)

    fbe = 1.0d0/(expf-1.0d0)

    end function

    function dbosedT(freq,temperature) result(fbe)

    implicit none

    real(kind=8) :: fbe
    real(kind=8) :: freq,temperature
    real(kind=8) :: expf

    expf = exp(hbar*freq*1.0d12*2*pi/kb/temperature)

    fbe = expf/(expf-1.0d0)**2*hbar*freq*1.0d12*2*pi/kb/temperature**2

    end function

    function simpson(x,y) result(f)

    implicit none

    real(kind=8) :: x(:),dx
    complex(kind=8) :: y(:)
    complex(kind=8) :: f
    integer(kind=4) :: n,i,mm

    n = size(x,1)
    dx = x(n)-x(n-1)
    f = 0.0d0
    if (mod(n-1,2) .eq.0) then
        mm = (n-1)/2
        f = f + y(1)
        do i = 1,mm
            f = f + 4.0d0*y(2*i)
        end do
        do i = 1,mm-1
            f = f + 2.0d0*y(2*i+1)
        end do
        f = f + y(n)
    else
        write(*,*) "ne should be an odd number" 
        stop
        
    end if
    f = f*dx/3.0d0
    end function

end module util
