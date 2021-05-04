module device
    
    implicit none

    integer(kind=4),allocatable :: layerstart(:),layerend(:)
    integer(kind=4),allocatable :: layer_list(:)
    integer(kind=4)   :: nlayer
    integer(kind=4) :: nd, maxlen

contains
    subroutine blockwise_setup(Ham,norb)
    ! this subroutine make the hamiltonian for the device
    ! region blockwise

    use surface, only : nl,nr,n_buffer_l,n_buffer_r
    use config
    use input    
    use write_data

    implicit none

    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(in) :: Ham(norb,norb)
    integer(kind=4) :: start, layer1, layer2
    integer(kind=4) :: start0, cont
    integer(kind=4) :: i,j, maxij1, maxij3, maxij2
    integer(kind=4),allocatable :: tmp(:)
    integer(kind=4),allocatable :: tmp1(:)
    complex(kind=8),allocatable :: VLD(:,:), VRD(:,:),HD(:,:)
    complex(kind=8),allocatable :: Vtemp(:,:), Vtemp_l(:),Vtemp_b(:)
    integer(kind=4) :: lay1,lay2, lastblock
    real(kind=8) :: positions_all(norb,5)


    nlayer = 1
    allocate(layer_list(nlayer))
    allocate(tmp(nlayer+1))

    ! dimensions of Hamiltonian for device region
    nd = norb - n_buffer_l - n_buffer_r - 2*nl - 2*nr

    ! Ham_LD
    allocate(VLD(nl,nd))
    VLD(:,:) = Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
                   n_buffer_l+2*nl+1:n_buffer_l+2*nl+nd)

    maxij1 = 0
    do i = 1,nl
        do j = 1,nd-1
            if ((abs(VLD(i,j)) .gt. 1.0d-7) &
            .and.(maxval(abs(VLD(i,j+1:nd))).lt.1.0d-7)&
            .and.(j.ge.maxij1)) then
                maxij1 = j
            end if
        end do
    end do

    ! Ham_RD
    maxij3 = 0
    allocate(VRD(nr,nd))
    VRD(:,:) = Ham(norb-n_buffer_r-2*nr+1:norb-n_buffer_r-nr,&
    norb-n_buffer_r-2*nr-nd+1:norb-n_buffer_r-2*nr)
    do i = 1, nr
       do j = 2,nd
           if ((abs(VRD(i,j)) .gt. 1.0d-7) &
           .and.(maxval(abs(VRD(i,1:j-1))).lt.1.0d-7)&
           .and.(nd-j+1.gt.maxij3)) then
               maxij3 = nd-j+1
            end if
       end do
    end do
    ! Hamiltonian for device region
    allocate(HD(nd,nd))
    HD(:,:) = Ham(n_buffer_l+2*nl+1:n_buffer_l+2*nl+nd,&
    n_buffer_l+2*nl+1:n_buffer_l+2*nl+nd)

    start = 1

    layer1 = 1
    layer2 = 1
    cont = 0
    start0 = 1

    allocate(Vtemp(1,1))
    allocate(Vtemp_l(1))
    allocate(Vtemp_b(1))
    do while ((start .lt. nd-1).and. (cont .lt. 2).and.(start .ge. start0))
        maxij2 = 0
        start0 = start
        if (cont .eq. 0) then
            do i = maxij1,nd-2
                do j = i+1,nd-1
                    layer1 = i+1-start ! dimension of first block
                    layer2 = j-i
                    ! coupling between i and i+2 layer
                    deallocate(Vtemp)
                    allocate(Vtemp(layer1,nd-start-layer1-layer2+1))
                    Vtemp(:,:) = HD(start:start+layer1-1,&
                    start+layer1+layer2:nd)
                    ! Elements to the left of V should be non-zero
                    deallocate(Vtemp_l)
                    allocate(Vtemp_l(layer1))
                    Vtemp_l(:) = HD(start:start+layer1-1,start+layer1+layer2-1)
                    ! Elements to the bottom of V should also be non-zero
                    deallocate(Vtemp_b)
                    allocate(Vtemp_b(nd-start-layer1-layer2+1))
                    Vtemp_b(:) = HD(start+layer1,start+layer1+layer2:nd)
                    if ((maxval(abs(Vtemp_l)).gt.1.0d-7).and.(maxval(abs(Vtemp_b)).gt.1.0d-7) &
                    .and. (maxval(abs(Vtemp)).lt.1.0d-7)&
                    .and.( size(Vtemp,2).gt. maxij2)) then
                        maxij2 = size(Vtemp,2)
                        lay1 = layer1
                        lay2 = layer2
                        if (nlayer .eq. 1)then
                            layer_list(1) = lay1
                            tmp(1) = lay1
                        else
                            tmp(nlayer) = lay1
                            deallocate(layer_list)
                            allocate(layer_list(nlayer))
                            layer_list(:) = tmp(:)
                            deallocate(tmp)
                            allocate(tmp(nlayer+1))
                            tmp(1:nlayer) = layer_list(:)
                        end if
                        nlayer = nlayer + 1

                    end if
                end do
            end do
            cont = cont+ 1
            start = start + lay1
        else if (cont .ne. 0) then
            layer1 = lay2
            cont = 2
            do j = start+layer1,nd-1
                layer2 = j+1-start-layer1
                ! coupling between i and i+2 layer
                deallocate(Vtemp)
                allocate(Vtemp(layer1,nd-start-layer1-layer2+1))
                Vtemp(:,:) = HD(start:start+layer1-1,&
                start+layer1+layer2:nd)
                ! Elements to the left of V should be non-zero
                deallocate(Vtemp_l)
                allocate(Vtemp_l(layer1))
                Vtemp_l(:) = HD(start:start+layer1-1,start+layer1+layer2-1)
                ! Elements to the bottom of V should also be non-zero
                deallocate(Vtemp_b)
                allocate(Vtemp_b(nd-start-layer1-layer2+1))
                Vtemp_b(:) = HD(start+layer1,start+layer1+layer2:nd)
                if ((maxval(abs(Vtemp_l)).gt.1.0d-7).and.(maxval(abs(Vtemp)).lt.1.0d-7)&
                .and.(size(Vtemp,2).gt.maxij2)) then
                    maxij2 = size(Vtemp,2)
                    lay1 = layer1
                    lay2 = layer2
                    if (nlayer .eq. 1)then
                        layer_list(1) = lay1
                        tmp(1) = lay1
                    else
                        tmp(nlayer) = lay1
                        deallocate(layer_list)
                        allocate(layer_list(nlayer))
                        layer_list(:) = tmp(:)
                        deallocate(tmp)
                        allocate(tmp(nlayer+1))
                        tmp(1:nlayer) = layer_list(:)
                    end if
                    nlayer = nlayer + 1
                    cont = 1
                    start = start + lay1
                    if (start+layer1+layer2.eq.nd+1)then
                        if (nlayer .eq. 1)then
                            layer_list(1) = lay1
                            tmp(1) = lay1
                        else
                            tmp(nlayer) = lay1
                            deallocate(layer_list)
                            allocate(layer_list(nlayer))
                            layer_list(:) = tmp(:)
                            deallocate(tmp)
                            allocate(tmp(nlayer+1))
                            tmp(1:nlayer) = layer_list(:)
                        end if
                        nlayer = nlayer + 1
                    end if
                end if
            end do
        end if
    end do

    tmp(nlayer) = nd - sum(layer_list)
    deallocate(layer_list)
    allocate(layer_list(nlayer))
    layer_list(:) = tmp(:)
    deallocate(tmp)
    allocate(tmp(nlayer+1))
    tmp(1:nlayer) = layer_list(:)

    ! Merge last two blocks if the last block is too small
    lastblock = size(layer_list,1)

    allocate(tmp1(1))
    do while (layer_list(lastblock) .lt. maxij3)
        deallocate(tmp1)
        allocate(tmp1(size(layer_list,1)))
        tmp1(:) = layer_list(:)
        deallocate(layer_list)
        lastblock = lastblock-1
        allocate(layer_list(lastblock))
        layer_list(1:lastblock) = tmp1(1:lastblock)
        layer_list(lastblock) = layer_list(lastblock) + tmp1(lastblock+1)
    end do
    nlayer = size(layer_list,1)

    ! the index for the start and the end of each block
    allocate(layerstart(nlayer))
    allocate(layerend(nlayer))

    layerstart(1) = 1
    layerend(1) = layer_list(1)
    do i = 2,nlayer
        layerstart(i) = layerend(i-1) + 1
        layerend(i) = layerstart(i) + layer_list(i) -1
    end do

    deallocate(tmp)
    deallocate(tmp1)

    maxlen = maxval(layer_list) ! largest block

    start = 0
    positions_all = 0.0d0
    do i =1,natoms
        do j = 1,3
            start = start + 1
            positions_all(start,1:4) = positions(i,1:4)
        end do
    end do
    do i = 1,nlayer
        ! fifth dimension is the index of layer
        do j = n_buffer_l+2*nl+layerstart(i),&
               n_buffer_l+2*nl+layerend(i)
            positions_all(j,5) = dble(i)
        end do
    end do
    
    call write_positions(positions_all)

    end subroutine blockwise_setup

    subroutine G_block_mpi(ie,ik_core,Ham,norb,gl,gr,gml,gmr,&
               sl,sr,E_k,transl,transr,G_D_RL,G_D_LR,&
               G_L,G_R)

    use config
    use surface, only : nl, nr

    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: ie, ik_core
    integer(kind=4),intent(in) :: norb
    integer(kind=4) :: mm, nn
    complex(kind=8),intent(in) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(in) :: gl(numprocs,nl,nl),gml(numprocs,nl,nl)
    complex(kind=8),intent(in) :: gr(numprocs,nr,nr),gmr(numprocs,nr,nr)
    complex(kind=8),intent(in) :: sl(numprocs,nl,nl),sr(numprocs,nr,nr)
    real(kind=8),intent(in) :: E_k(ne,nk)
    real(kind=8) :: transr(ne,nk),transl(ne,nk)
    complex(kind=8),intent(out) :: G_D_RL(numprocs,nr,nl)
    complex(kind=8),intent(out) :: G_D_LR(numprocs,nl,nr)
    complex(kind=8),intent(out) :: G_L(numprocs,nl,nl)
    complex(kind=8),intent(out) :: G_R(numprocs,nr,nr)


    mm=myid+1
    nn = k_start_core(mm)+ik_core-1
    call G_block(mm,Ham(mm,:,:),norb,&
         transpose(dconjg(gl(mm,:,:))),gr(mm,:,:),&
         -gml(mm,:,:),gmr(mm,:,:),&
         transpose(dconjg(sl(mm,:,:))),sr(mm,:,:),E_k(ie,nn),&
         transl(ie,nn),transr(ie,nn),&
         G_D_RL(mm,:,:),G_D_LR(mm,:,:),&
         G_L(mm,:,:),G_R(mm,:,:))

    end subroutine G_block_mpi

    subroutine G_block(mm,Ham,norb,gl,gr,gml,gmr,sl,sr,&
                       E,transl,transr,G_D_RL,G_D_LR,&
                       G_L,G_R)

    use util
    use config
    use surface , only : nl, nr, n_buffer_l, n_buffer_r

    implicit none

    integer(kind=4),intent(in) :: mm
    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(in) :: Ham(norb,norb)
    complex(kind=8),intent(in) :: gl(nl,nl)
    complex(kind=8),intent(in) :: gr(nr,nr)
    complex(kind=8),intent(in) :: gml(nl,nl)
    complex(kind=8),intent(in) :: gmr(nr,nr)
    complex(kind=8),intent(in) :: sl(nl,nl)
    complex(kind=8),intent(in) :: sr(nr,nr)
    real(kind=8),intent(in) :: E
    real(kind=8) :: transl
    real(kind=8) :: transr
    complex(kind=8),intent(out) :: G_D_RL(nr,nl), G_D_LR(nl,nr)
    complex(kind=8),intent(out) :: G_L(nl,nl), G_R(nr,nr)

    complex(kind=8) :: GGGGl(nl,nl), GGGGr(nr,nr)
    complex(kind=8),allocatable :: Es(:,:),Vcoup(:,:),&
                                   H(:,:)
    complex(kind=8),allocatable :: Gblock_NN(:,:),&
                                 Gblock_oN(:,:),Gblock_No(:,:)
    complex(kind=8),allocatable :: Gblock_NN0(:,:),&
                            Gblock_oN0(:,:),Gblock_No0(:,:)
    integer(kind=4) :: i, device_start, j
    complex(kind=8),allocatable :: GRinv(:,:), GLinv(:,:)


    ! index that marks the start of device region
    device_start = n_buffer_l+2*nl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Left to right sweep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(Es(layer_list(1),layer_list(1)),&
            Vcoup(nl,layer_list(1)),&
            H(layer_list(1),layer_list(1)),&
            Gblock_NN(layer_list(1),layer_list(1)),&
            Gblock_oN(nl,layer_list(1)),&
            Gblock_NN0(layer_list(1),layer_list(1)),&
            Gblock_oN0(nl,layer_list(1)))

    ! Ham_LD
    Vcoup = Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
        device_start+layerstart(1):&
        device_start+layerend(1))
    H = Ham(device_start + layerstart(1):&
            device_start + layerend(1),&
            device_start + layerstart(1):&
            device_start + layerend(1))
    Es = (E+10.0d0*(E*eta+eta0)*i_imag)*eyemat(layer_list(1))
    ! G_1, 1
    Gblock_NN = inv(Es-H-matmul(matmul(transpose(dconjg(Vcoup)),gl(:,:)),Vcoup))
    ! G_1, 1
    Gblock_oN = matmul(matmul(gl(:,:),Vcoup),Gblock_NN)
    Gblock_NN0 = Gblock_NN
    Gblock_oN0 = Gblock_oN

    do i = 2, nlayer
        deallocate(Vcoup)
        allocate(Vcoup(layer_list(i-1),layer_list(i)))
        Vcoup = Ham(device_start+layerstart(i-1):&
                device_start+layerend(i-1),&
                device_start+layerstart(i):&
                device_start+layerend(i)) ! V_n_n+1

        deallocate(Es)
        allocate(Es(layer_list(i),layer_list(i)))
        Es = (E+10.0d0*(E*eta+eta0)*i_imag)*eyemat(layer_list(i))

        deallocate(H)
        allocate(H(layer_list(i),layer_list(i)))
        H = Ham(device_start+layerstart(i):&
                device_start+layerend(i),&
                device_start+layerstart(i):&
                device_start+layerend(i)) ! H_n+1_n+1

        deallocate(Gblock_NN)
        allocate(Gblock_NN(layer_list(i),layer_list(i)))
        Gblock_NN = inv(Es - H -&
        matmul(matmul(transpose(dconjg(Vcoup)),Gblock_NN0),&
        Vcoup)) ! G_n+1_n+1

        deallocate(Gblock_NN0)
        allocate(Gblock_NN0(layer_list(i),layer_list(i)))
        Gblock_NN0 = Gblock_NN

        deallocate(Gblock_oN0)
        allocate(Gblock_oN0(nl,layer_list(i-1)))
        Gblock_oN0 = Gblock_oN ! G_0_n

        deallocate(Gblock_oN)
        allocate(Gblock_oN(nl,layer_list(i)))
        Gblock_oN = matmul(matmul(Gblock_oN0,Vcoup),Gblock_NN) ! G_0_n+1
    end do

    deallocate(Gblock_NN0)
    allocate(Gblock_NN0(layer_list(nlayer),layer_list(nlayer)))
    Gblock_NN0 = Gblock_NN

    deallocate(Es)
    allocate(Es(nr,nr))
    Es = (E+10.0d0*(E*eta+eta0)*i_imag)*eyemat(nr)

    allocate(GRinv(nr,nr))
    GRinv = Es - Ham(norb-n_buffer_r-2*nr+1:&
    norb-n_buffer_r-nr,norb-n_buffer_r-2*nr+1:&
    norb-n_buffer_r-nr)- sr(:,:)

    deallocate(Vcoup)
    allocate(Vcoup(layer_list(nlayer),nr))

    Vcoup(:,:) = Ham(device_start+layerstart(nlayer):&
                     device_start+layerend(nlayer),&
                     norb-n_buffer_r-2*nr+1:&
                     norb-n_buffer_r-nr)
    deallocate(Gblock_NN)
    allocate(Gblock_NN(nr,nr))
    Gblock_NN = inv(GRinv -&
                matmul(matmul(transpose(dconjg(Vcoup)),&
                Gblock_NN0),Vcoup))

    deallocate(Gblock_oN0)
    allocate(Gblock_oN0(nl,layer_list(nlayer)))
    Gblock_oN0 = Gblock_oN ! G_0_n

    deallocate(Gblock_oN)
    allocate(Gblock_oN(nl,nr))
    Gblock_oN = matmul(matmul(Gblock_oN0,Vcoup),Gblock_NN) ! G_0_n+1

    G_R =  Gblock_NN ! right green function in full green function
    G_D_LR = Gblock_oN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Right to left sweep
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    deallocate(Vcoup,H,Es,Gblock_NN,Gblock_NN0)
    allocate(Vcoup(nr,layer_list(nlayer)),&
            H(layer_list(nlayer),layer_list(nlayer)),&
            Es(layer_list(nlayer),layer_list(nlayer)),&
            Gblock_NN(layer_list(nlayer),layer_list(nlayer)),&
            Gblock_NN0(layer_list(nlayer),layer_list(nlayer)),&
            Gblock_No(nr,layer_list(nlayer)),&
            Gblock_No0(nr,layer_list(nlayer)))

    Vcoup = Ham(norb-n_buffer_r-2*nr+1:&
                norb-n_buffer_r-nr,&
                device_start+layerstart(nlayer):&
                device_start+layerend(nlayer)) ! V_N+1_N

    H = Ham(device_start+layerstart(nlayer):&
            device_start+layerend(nlayer),&
            device_start+layerstart(nlayer):&
            device_start+layerend(nlayer)) ! H_N_N
    Es = (E+10.0d0*(E*eta+eta0)*i_imag)*eyemat(layer_list(nlayer))

    ! G_N, N
    Gblock_NN = inv(Es-H-&
                matmul(matmul(transpose(dconjg(Vcoup)),gr(:,:)),Vcoup))
    Gblock_No = matmul(matmul(gr(:,:),Vcoup),Gblock_NN)
    Gblock_NN0 = Gblock_NN
    Gblock_No0 = Gblock_No
    
    do i = nlayer,2,-1
        deallocate(Vcoup)
        allocate(Vcoup(layer_list(i),layer_list(i-1)))
        Vcoup = Ham(device_start+layerstart(i):&
                    device_start+layerend(i),&
                    device_start+layerstart(i-1):&
                    device_start+layerend(i-1)) ! V_n+1_n

        deallocate(Es)
        allocate(Es(layer_list(i-1),layer_list(i-1)))
        Es = (E + 10.0d0*(E*eta+eta0)*i_imag) * eyemat(layer_list(i-1))

        deallocate(H)
        allocate(H(layer_list(i-1),layer_list(i-1)))
        H = Ham(device_start+layerstart(i-1):&
                device_start+layerend(i-1),&
                device_start+layerstart(i-1):&
                device_start+layerend(i-1)) ! H_n_n

        deallocate(Gblock_NN)
        allocate(Gblock_NN(layer_list(i-1),layer_list(i-1)))
        Gblock_NN = inv(Es-H-&
               matmul(matmul(transpose(dconjg(Vcoup)),&
               Gblock_NN0),Vcoup)) ! G_n_n
        if (abs(E - 31.2953890185127).lt. 0.001) then
            Gblock_NN = matmul(Es-H-&
                        matmul(matmul(transpose(dconjg(Vcoup)),&
                        Gblock_NN0),Vcoup),Gblock_NN)
                        do j = 1,size(Gblock_NN,1)
                            write(*,*) Gblock_NN(j,j)
                        end do
!            stop
        end if

        deallocate(Gblock_NN0)
        allocate(Gblock_NN0(layer_list(i-1),layer_list(i-1)))
        Gblock_NN0 = Gblock_NN

        deallocate(Gblock_No0)
        allocate(Gblock_No0(nr,layer_list(i)))
        Gblock_No0 = Gblock_No ! G_N+1_n+1

        deallocate(Gblock_No)
        allocate(Gblock_No(nr,layer_list(i-1)))
        Gblock_No = matmul(matmul(Gblock_No0,Vcoup)&
        ,Gblock_NN) ! G_N+1_n
    end do

    deallocate(Gblock_NN0)
    allocate(Gblock_NN0(layer_list(1),layer_list(1)))
    Gblock_NN0 = Gblock_NN

    deallocate(Es)
    allocate(Es(nl,nl))
    Es = (E+10.0d0*(E*eta+eta0)*i_imag)*eyemat(nl)

    allocate(GLinv(nl,nl))
    GLinv = Es - Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
       n_buffer_l+nl+1:n_buffer_l+2*nl) - sl(:,:) ! G_1_1

    deallocate(Vcoup)
    allocate(Vcoup(nl,layer_list(1)))
    Vcoup = Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
               device_start+layerstart(1):&
               device_start+layerend(1)) ! V_0_1

    deallocate(Gblock_NN)
    allocate(Gblock_NN(nl,nl))
    Gblock_NN = inv(GLinv - matmul(matmul(Vcoup,Gblock_NN0),&
              transpose(dconjg(Vcoup)))) ! G_0_0

    deallocate(Gblock_No0)
    allocate(Gblock_No0(nr,layer_list(1)))
    Gblock_No0 = Gblock_No ! G_N_1

    deallocate(Gblock_No)
    allocate(Gblock_No(nr,nl))
    Gblock_No = matmul(matmul(Gblock_No0,&
           transpose(dconjg(Vcoup))),Gblock_NN) ! G_R_N+1_0

    G_L =  Gblock_NN ! left green function in full green function
    G_D_RL = Gblock_No

    ! G_0_N+1 Gamma_R G_N+1_0 Gamma_L
    GGGGl = matmul(matmul(Gblock_oN,gmr(:,:)),&
                   matmul(transpose(dconjg(Gblock_oN)),gml(:,:)))
    ! G_N+1_0 Gamma_L G_0_N+1 Gamma_R
    GGGGr = matmul(matmul(Gblock_No,gml(:,:)),&
                   matmul(transpose(dconjg(Gblock_No)),gmr(:,:)))
    transl = real(trace(GGGGl))
    transr = real(trace(GGGGr))

    if (transl .lt. 1.0d-60) then
        transl = 0.0d0
    end if
    if (transr .lt. 1.0d-60) then
        transr = 0.0d0
    end if

    end subroutine G_block

end module device
