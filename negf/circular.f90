module circular

    implicit none

    complex(kind=8),allocatable :: P1(:,:), P2(:,:), &
                                   Pinv1(:,:), Pinv2(:,:)
contains

    subroutine set_pmatrix(n,nb,G,R,P,Pinv)
    
    use config
    use util

    implicit none

    integer(kind=4),intent(in) :: n
    integer(kind=4),intent(in) :: nb
    real(kind=8),intent(in) :: G(3)
    real(kind=8),intent(in) :: R(3)
    integer(kind=4) :: na
    integer(kind=4) :: i,j,k
    complex(kind=8),intent(out) :: P(n,n),Pinv(n,n)
    na = n/nb
    P = 0.0d0
    Pinv = 0.0d0
    do i = 1,nb
        do j = 1,nb
             do k = 1,na
                 P((i-1)*na+k,&
                   (j-1)*na+k) = &
                 exp(i_imag*dot_product(&
                 G*dble(j-1)/dble(nb),&
                 dble(i-1)*R))/sqrt(dble(nb))
            end do
        end do
    end do
    
    ! P is unitary matrix
    Pinv = dconjg(transpose(P))     
    end subroutine set_pmatrix
    

    function circular_g(Es_in0,e_in1,alpha_in,&
             beta_in,convergence,&
             n,nb,nc,G1,R1,G2,R2) result(green1)

    ! We have a difference of 1d-10 compared with
    ! direct matrix inversion

    use util    
 
    implicit none

    complex(kind=8),intent(in) :: Es_in0(:,:)
    complex(kind=8),intent(in) :: e_in1(:,:)
    complex(kind=8),intent(in) :: alpha_in(:,:)
    complex(kind=8),intent(in) :: beta_in(:,:)
    real(kind=8),intent(in) :: convergence
    integer(kind=4),intent(in) :: n, nb, nc
    real(kind=8),intent(in) :: G1(3),G2(3)
    real(kind=8),intent(in) :: R1(3),R2(3)
    integer(kind=4) :: i,j,na,nac, k
    complex(kind=8) :: Es0(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: e1(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: alpha(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: beta(size(Es_in0,1),size(Es_in0,2))
    complex(kind=8) :: green1(size(Es0,1),size(Es0,2))
    complex(kind=8),allocatable :: green2(:,:)
    complex(kind=8),allocatable :: Esy0(:,:)
    complex(kind=8),allocatable :: ey1(:,:)
    complex(kind=8),allocatable :: alphay(:,:)
    complex(kind=8),allocatable :: betay(:,:)

    na = n/nb
    allocate(green2(na,na),Esy0(na,na),ey1(na,na),&
             alphay(na,na),betay(na,na))

    nac = na/nc
    green1 = 0.0d0

    Es0 = fourier_transform(Es_in0,n,nb,G1,R1) 
    e1 = fourier_transform(e_in1,n,nb,G1,R1) 
    alpha = fourier_transform(alpha_in,n,nb,G1,R1) 
    beta = fourier_transform(beta_in,n,nb,G1,R1) 
 
    do i = 1,nb
        Esy0 = fourier_transform(Es0((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
               na,nc,G2,R2) 
        ey1 = fourier_transform(e1((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
              na,nc,G2,R2) 
        alphay = fourier_transform(alpha((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
                 na,nc,G2,R2) 
        betay = fourier_transform(beta((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
                na,nc,G2,R2) 
        green2 = 0.0d0
        do j = 1,nc
            green2((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac) = &
                 decimation(Esy0((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac),&
                 ey1((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac),&
                 alphay((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac),&
                 betay((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac),convergence)
!            do k = 1, 12
!                write(*,*) k
!                write(*,*) ey1((j-1)*nac+k,(j-1)*nac+1:j*nac)
!             write(*,*) maxval(abs(ey1((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac)-dconjg(transpose(ey1((j-1)*nac+1:j*nac,(j-1)*nac+1:j*nac)))))
!            end do
!             stop
        end do
        green1((i-1)*na+1:i*na,(i-1)*na+1:i*na) = &
            matmul(matmul(P2,green2),Pinv2)
    end do
    green1 = matmul(matmul(P1,green1),Pinv1)

    end function circular_g

    

    function circular_inv(H,n,nb,nc,&
             G1,R1,G2,R2) result(Hinv1)

    use util

    implicit none

    integer(kind=4),intent(in) :: n, nb, nc
    real(kind=8),intent(in) :: G1(3),G2(3)
    real(kind=8),intent(in) :: R1(3),R2(3)
    complex(kind=8),intent(in) :: H(:,:)
    complex(kind=8) :: Hinv1(size(H,1),size(H,2))
    complex(kind=8) :: Hinv2(size(H,1)/nb,size(H,2)/nb)
    complex(kind=8) :: H1(size(H,1),size(H,2))
    complex(kind=8) :: H2(size(H,1)/nb,size(H,2)/nb)
    integer(kind=4) :: i,j, na, nac
    
    na = n/nb
    nac = na/nc
    Hinv1 = 0.0d0

    H1 = fourier_transform(H,n,nb,G1,R1)
    do i = 1,nb
        H2 = fourier_transform(H1((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
            na,nc,G2,R2)
        Hinv2 = 0.0d0
        do j = 1,nc
            Hinv2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac) = &
            inv(H2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac))
        end do
        Hinv1((i-1)*na+1:i*na,(i-1)*na+1:i*na)&
        = matmul(matmul(P2,Hinv2),Pinv2)
    end do
    Hinv1 = matmul(matmul(P1,Hinv1),Pinv1)

    end function circular_inv

    subroutine circular_eigen_axb(Ha,Hb,n,nb,nc,&
                 G1,R1,G2,R2,eigs1_diag,evs1)

    use util

    implicit none

    integer(kind=4),intent(in) :: n, nb, nc
    real(kind=8),intent(in) :: G1(3),G2(3)
    real(kind=8),intent(in) :: R1(3),R2(3)
    complex(kind=8),intent(in) :: Ha(:,:),Hb(:,:) ! Ha, Hb shares same periodicity
    complex(kind=8)            :: Ha1(size(Ha,1),size(Ha,2)),Hb1(size(Hb,1),size(Hb,2))
    complex(kind=8)            :: Ha2(size(Ha,1)/nb,size(Ha,2)/nb),Hb2(size(Hb,1)/nb,size(Hb,2)/nb)
    complex(kind=8) :: Hd2(size(Ha,1)/nb,size(Ha,2)/nb)
    complex(kind=8),intent(out) :: evs1(size(Ha,1),size(Ha,2))
    complex(kind=8) :: evs2(size(Ha,1)/nb,size(Ha,2)/nb)
    complex(kind=8),intent(out) :: eigs1_diag(size(Ha,1))
    complex(kind=8) :: eigs1(size(Ha,1),size(Ha,2))
    complex(kind=8) :: eigs2(size(Ha,1)/nb,size(Ha,2)/nb)
    complex(kind=8) :: eigs3(size(Ha,1)/nb/nc)
    complex(kind=8) :: evs3(size(Ha,1)/nb/nc,size(Ha,1)/nb/nc)
    integer(kind=4) :: i,j, na, nac
    
    na = n/nb
    nac = na/nc

    Ha1 = fourier_transform(Ha,n,nb,G1,R1)
    Hb1 = fourier_transform(Hb,n,nb,G1,R1)

    eigs1 = 0.0d0
    evs1 = 0.0d0

    do i = 1,nb
        Ha2 = fourier_transform(Ha1((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
            na,nc,G2,R2)
        Hb2 = fourier_transform(Hb1((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
            na,nc,G2,R2)
        Hd2 = matmul(Ha2,Hb2)

        eigs2 = 0.0d0
        evs2 = 0.0d0
        do j = 1,nc
            call eigen_sorting(Hd2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac),eigs3,evs3)
            eigs2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac) = diagmat(eigs3) 
            evs2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac) = evs3
        end do

        eigs1((i-1)*na+1:i*na,(i-1)*na+1:i*na) = eigs2
        evs1((i-1)*na+1:i*na,(i-1)*na+1:i*na) = matmul(P2,evs2)      
    end do

    evs1 = matmul(P1,evs1)

    do i = 1,n
        eigs1_diag(i) = eigs1(i,i) 
    end do

    end subroutine circular_eigen_axb

    subroutine circular_eigen(H,n,nb,nc,&
                 G1,R1,G2,R2,eigs1_diag,evs1)

    use util

    implicit none

    integer(kind=4),intent(in) :: n, nb, nc
    real(kind=8),intent(in) :: G1(3),G2(3)
    real(kind=8),intent(in) :: R1(3),R2(3)
    complex(kind=8),intent(in) :: H(:,:)
    complex(kind=8)            :: H1(size(H,1),size(H,2))
    complex(kind=8) :: H2(size(H,1)/nb,size(H,2)/nb)
    complex(kind=8),intent(out) :: evs1(size(H,1),size(H,2))
    complex(kind=8) :: evs2(size(H,1)/nb,size(H,2)/nb)
    complex(kind=8),intent(out) :: eigs1_diag(size(H,1))
    complex(kind=8) :: eigs1(size(H,1),size(H,2))
    complex(kind=8) :: eigs2(size(H,1)/nb,size(H,2)/nb)
    complex(kind=8) :: eigs3(size(H,1)/nb/nc)
    complex(kind=8) :: evs3(size(H,1)/nb/nc,size(H,1)/nb/nc)
    integer(kind=4) :: i,j, na, nac
    
    na = n/nb
    nac = na/nc

    H1 = fourier_transform(H,n,nb,G1,R1)

    eigs1 = 0.0d0
    evs1 = 0.0d0

    do i = 1,nb
        H2 = fourier_transform(H1((i-1)*na+1:i*na,(i-1)*na+1:i*na),&
            na,nc,G2,R2)

        eigs2 = 0.0d0
        evs2 = 0.0d0

        do j = 1,nc
            call eigen_sorting(H2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac),eigs3,evs3)
            eigs2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac) = diagmat(eigs3) 
            evs2((j-1)*nac+1:j*nac,&
            (j-1)*nac+1:j*nac) = evs3
        end do

        eigs1((i-1)*na+1:i*na,(i-1)*na+1:i*na) = eigs2
        evs1((i-1)*na+1:i*na,(i-1)*na+1:i*na) = matmul(P2,evs2)      
    end do

    evs1 = matmul(P1,evs1)

    do i = 1,n
        eigs1_diag(i) = eigs1(i,i) 
    end do

    end subroutine circular_eigen

end module circular
