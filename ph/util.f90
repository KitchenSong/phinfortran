module util

implicit none

real(kind=8),parameter  :: pi=3.141592653589793238d0
real(kind=8),parameter  :: angbohr = 1.889725989d0
real(kind=8),parameter  :: bohr2ang = 0.529177d0 
real(kind=8),parameter  :: mass_proton = 1.6726219d-27
real(kind=8),parameter  :: eleVolt = 1.602176634d-19
complex(kind=8),parameter :: i_imag=cmplx(0,1,kind=8)
contains

function linspace(k1,k2,nk,endpoint) result(k)

real(kind=8) :: k1,k2
integer(kind=4) :: nk
integer(kind=4) :: i
real(kind=8),allocatable :: k(:)
integer(kind=4) :: endpoint

allocate(k(nk))
if (endpoint .ne. 0) then
    do i = 1,nk
        k(i) = dble(i-1)/dble(nk-1)*(k2-k1) + k1 
    end do
else
    do i = 1,nk
        k(i) = dble(i-1)/dble(nk)*(k2-k1) + k1 
    end do
end if

end function

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
FUNCTION cross(a, b)
  real(kind=8) :: cross(3)
  real(kind=8),dimension(3),INTENT(IN) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross
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

function hermitian(A) result(B)

    implicit none

    complex(kind=8) :: A(:,:)
    complex(kind=8) :: B(size(A,1),size(A,2))

    B = 0.5*(A+transpose(dconjg(A)))

end function hermitian

function participation(A) result(B)

    implicit none

    complex(kind=8) :: A(:)
    integer(kind=4) :: n,i,j
    real(kind=8)    :: B
    complex(kind=8) :: C,D

    n = size(A,1)
    C = 0.0d0
    D = 0.0d0

    do i = 1,n/3
        C = C + dot_product(A(3*(i-1)+1:3*i),A(3*(i-1)+1:3*i))
        D = D + dot_product(A(3*(i-1)+1:3*i),A(3*(i-1)+1:3*i))**2
    end do
    C = C**2
    D = D* n/3
    B = C/D

   
end function participation

function taumn(Temp,taum,taun,wm,wn) result(tmn)
    
    implicit none

    real(kind=8) :: temp,tmn,taum,taun,wm,wn

    tmn = (wm+wn)**2/4/wm/wn*(taum+taun)/((taum+taun)**2+(wm-wn)**2)+&
        (wm-wn)**2/4/wm/wn*(taum+taun)/((taum+taun)**2+(wm+wn)**2)



end function






end module

