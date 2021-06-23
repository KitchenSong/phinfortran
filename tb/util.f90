module util

implicit none

real(kind=8),parameter  :: pi=3.141592653589793238d0
complex(kind=8),parameter :: i_imag=cmplx(0,1,kind=8)

contains
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  complex(kind=8), dimension(:,:), intent(in) :: A
  complex(kind=8), dimension(size(A,1),size(A,2)) :: Ainv

  complex(kind=8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external ZGETRF
  external ZGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! ZGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call ZGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     write(*,*) size(A,1)

     stop 'Matrix is numerically singular!'
  end if

  ! ZGETRI computes the inverse of a matrix using the LU factorization
  ! computed by ZGETRF.
  call ZGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
  end function inv

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


  function eyemat(n) result(mat)
  integer(kind=4),intent(in) :: n
  real(kind=8) :: mat(n,n)
  integer(kind=4) :: i
  mat = 0.0d0
  do i = 1, n
    mat(i,i) = 1.0d0
  end do
  end function eyemat

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
        nwork=nint(2*real(work(1)))
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
        nwork=nint(2*real(work(1)))
        deallocate(work)
        allocate(work(nwork))
  end if

  call ZGEEV("N", "V", lda, Atemp, lda, w,&
   vecl, lda, vecr, lda,work,nwork,&
   rwork, info) 

  end subroutine eigen

  LOGICAL FUNCTION DUMMY_SELECT(ARG) 
  IMPLICIT NONE


  COMPLEX(kind=8),intent(in) :: ARG

  DUMMY_SELECT = (ARG .EQ. ARG)
 
  end function DUMMY_SELECT


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
  sdim = 0

  S = A

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
        nwork=nint(2*real(work(1)))
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
         if (U(i,i) + U(j,j) .ne. cmplx(0.0,0.0,kind=8))then
             U(i,j) = U(i,j) / (U(i,i) + U(j,j))
         else
             U(i,j) = 0.0d0
         end if
     enddo
  enddo
  call ZGEMM('N', 'N', N, N, N, cmplx(1.0,0.0,kind=8), Q, N, U, N, cmplx(0.0,0.0,kind=8), QU, N)
  call ZGEMM('N', 'C', N, N, N, cmplx(1.0,0.0,kind=8), QU, N,Q, N, cmplx(0.0,0.0,kind=8), QUQH, N)
  if(info /= 0) then
       write(*,*) "zgemm returned info = ", info
       write(*,*) "the ", info, "-th argument had an illegal value"
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

    if (abs(dist_min2-dist_min).lt.1.0d-5) then
        output = 2
    else if (abs(dist(63)-dist_min) .lt. 1d-5) then
        output = 1
    else
        output = 0
    end if

    end function in_ws 

    function in_ws_misori(a_vectors,k_pos,n1,n2,n3) result(output)

    implicit none
   
    real(kind=8),intent(in) :: a_vectors(3,3), k_pos(3)
    real(kind=8) :: dist((2*n1+1)*(2*n2+1)*(2*n3+1)),dist_min,dist_min2
    real(kind=8) :: dist_remain((2*n1+1)*(2*n2+1)*(2*n3+1)-1)
    integer(kind=4) :: min_idx
    integer(kind=4) :: icnt,i1,i2,i3
    integer(kind=4),intent(in) :: n1,n2,n3
    real(kind=8) :: ndiff(3)
    integer(kind=4) :: output
    integer(kind=4) :: ntot

    ntot = (2*n1+1)*(2*n2+1)*(2*n3+1)

    dist = 0
    icnt = 0
    do i1 =-n1,n1
        do i2 =-n2,n2
            do i3 =-n3,n3
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
        dist_remain = dist(2:ntot)
    else if (min_idx .eq. ntot) then
        dist_remain = dist(1:ntot-1)
    else
        dist_remain(1:min_idx-1) = dist(1:min_idx-1)
        dist_remain(min_idx:ntot-1) = dist(min_idx+1:ntot)
    end if
    dist_min2 = minval(dist_remain)

    if (abs(dist_min2-dist_min).lt.1.0d-5) then
        output = 2
    else if (abs(dist((ntot+1)/2)-dist_min) .lt. 1d-5) then
        output = 1
    else
        output = 0
    end if

    end function in_ws_misori


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
        if ((abs(coord_crys_temp(1)-nint(coord_crys_temp(1)))<1d-4).and.&
            (abs(coord_crys_temp(2)-nint(coord_crys_temp(2)))<1d-4).and.&
            (abs(coord_crys_temp(3)-nint(coord_crys_temp(3)))<1d-4)) then
            idx = int(pos_pc(i,4))
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
    beta = beta_in
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

    function find_min_dist(kpoint,n) result(kmin) 

    implicit none

    real(kind=8) :: kpoint(:,:)
    integer(kind=4) :: n,i
    real(kind=8) :: kabs(n)
    real(kind=8) :: kmin(3)
    integer(kind=4) :: minid

    do i = 1,n
        kabs(i) = dsqrt(dot_product(kpoint(i,:),kpoint(i,:)))
    end do

    minid = minloc(kabs,1)
    kmin = kpoint(minid,:)
    end function

    function hermitian(A) result(B)

    implicit none

    complex(kind=8) :: A(:,:)
    complex(kind=8) :: B(size(A,1),size(A,2))

    B = 0.5*(A+transpose(dconjg(A)))

    end function hermitian

    function antihermitian(A) result(B)

    implicit none

    complex(kind=8) :: A(:,:)
    complex(kind=8) :: B(size(A,1),size(A,2))

    B = 0.5*(A-transpose(dconjg(A)))

    end function antihermitian



end module util
