module lesser

    implicit none
      
    real(kind=8),allocatable    :: temperature_i(:)  
    complex(kind=8),allocatable    :: sigma_w(:,:)
    complex(kind=8),allocatable    :: dsigma_wdt(:,:)
    complex(kind=8),allocatable  :: jacobian(:,:,:,:)

contains

subroutine set_sigma_bp_w(temp_i)
    
    use config
    use util

    implicit none

    integer(kind=4) :: i,j,n
    real(kind=8),allocatable  :: tau(:)
    real(kind=8) :: temp_i(:)
    real(kind=8)    :: B,T,C

    n = size(temp_i,1)
    allocate(tau(ne))
    T = 300.0d0
    B = B_probe !5.0d-20*(2.0d0*3.1415926d0*1.0d12)
    C = T_probe !430.0d0

    do j = 1,n
        T = temp_i(j)
        tau = 1.0d0/(B*Egrid(:,1)*T*exp(-C/T))
        do i = 1,ne
            sigma_w(j,i) = -2.0d0*i_imag*sqrt(Egrid(i,1))/tau(i)
            dsigma_wdt(j,i) =  -2.0d0*i_imag*sqrt(Egrid(i,1))*&
             (B*Egrid(i,1)*exp(-C/T)*(1.0d0+C/T))
            ! dsigma_wdt = sigma_w/T*(1+C/T)
        end do
    end do
    if (lifetime .eq. 0) then
        sigma_w = 0.0d0
    end if

    
end subroutine

subroutine get_G_full_mpi(ie,ik_core,Ham,norb,gl,gr,gml,gmr,&
               sl,sr,E_k,GLeft,G_full,&
               G_qp_q,G_q_qp,G_N_Np,G_Np_N,&
               G_zp_z,G_z_zp,&
               transl,transr,G_D_RL,G_D_LR,&
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
    complex(kind=8),intent(out) :: GLeft(:,:,:,:)
    complex(kind=8),intent(out) :: G_full(:,:,:,:)
    complex(kind=8),intent(out) :: G_qp_q(:,:,:,:),G_q_qp(:,:,:,:)
    complex(kind=8),intent(out) :: G_Np_N(:,:,:),G_N_Np(:,:,:)
    complex(kind=8),intent(out) :: G_zp_z(:,:,:),G_z_zp(:,:,:)
    real(kind=8) :: transr(ne,nk),transl(ne,nk)
    complex(kind=8),intent(out) :: G_D_RL(numprocs,nr,nl)
    complex(kind=8),intent(out) :: G_D_LR(numprocs,nl,nr)
    complex(kind=8),intent(out) :: G_L(numprocs,nl,nl)
    complex(kind=8),intent(out) :: G_R(numprocs,nr,nr)



    mm=myid+1
    nn = k_start_core(mm)+ik_core-1

    call get_G_full(ie,mm,Ham(mm,:,:),norb,&
         transpose(dconjg(gl(mm,:,:))),gr(mm,:,:),&
         -gml(mm,:,:),gmr(mm,:,:),&
         transpose(dconjg(sl(mm,:,:))),sr(mm,:,:),E_k(ie,nn),&
         GLeft(mm,:,:,:),G_full(mm,:,:,:),&
         G_qp_q(mm,:,:,:),G_q_qp(mm,:,:,:),&
         G_Np_N(mm,:,:),G_N_Np(mm,:,:),&
         G_zp_z(mm,:,:),G_z_zp(mm,:,:),&
         transl(ie,nn),transr(ie,nn),&
         G_D_RL(mm,:,:),G_D_LR(mm,:,:),&
         G_L(mm,:,:),G_R(mm,:,:))
 
 
end subroutine

subroutine get_G_full_tdb_mpi(ie,ik_core,gl,gr,gml,gmr,&
               sl,sr,E_k,GLeft,G_full,&
               G_qp_q,G_q_qp,G_N_Np,G_Np_N,&
               G_zp_z,G_z_zp,&
               transl,transr,G_D_RL,G_D_LR,&
               G_L,G_R,&
               hamd,hamdv,hamld,hamrd,&
               haml,hamr)

    use config
    use surface, only : nl, nr

    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: ie, ik_core
    integer(kind=4) :: mm, nn
    complex(kind=8),intent(in) :: gl(numprocs,nl,nl),gml(numprocs,nl,nl)
    complex(kind=8),intent(in) :: gr(numprocs,nr,nr),gmr(numprocs,nr,nr)
    complex(kind=8),intent(in) :: sl(numprocs,nl,nl),sr(numprocs,nr,nr)
    real(kind=8),intent(in) :: E_k(ne,nk)
    complex(kind=8),intent(out) :: GLeft(:,:,:,:)
    complex(kind=8),intent(out) :: G_full(:,:,:,:)
    complex(kind=8),intent(out) :: G_qp_q(:,:,:,:),G_q_qp(:,:,:,:)
    complex(kind=8),intent(out) :: G_Np_N(:,:,:),G_N_Np(:,:,:)
    complex(kind=8),intent(out) :: G_zp_z(:,:,:),G_z_zp(:,:,:)
    real(kind=8) :: transr(ne,nk),transl(ne,nk)
    complex(kind=8),intent(out) :: G_D_RL(numprocs,nr,nl)
    complex(kind=8),intent(out) :: G_D_LR(numprocs,nl,nr)
    complex(kind=8),intent(out) :: G_L(numprocs,nl,nl)
    complex(kind=8),intent(out) :: G_R(numprocs,nr,nr)
    complex(kind=8),intent(in) :: hamd(:,:,:,:),hamdv(:,:,:,:)
    complex(kind=8),intent(in) :: hamld(:,:,:),hamrd(:,:,:)
    complex(kind=8),intent(in) :: haml(:,:,:),hamr(:,:,:)



    mm=myid+1
    nn = k_start_core(mm)+ik_core-1

    call get_G_full_tdb(ie,mm,&
         transpose(dconjg(gl(mm,:,:))),gr(mm,:,:),&
         -gml(mm,:,:),gmr(mm,:,:),&
         transpose(dconjg(sl(mm,:,:))),sr(mm,:,:),E_k(ie,nn),&
         GLeft(mm,:,:,:),G_full(mm,:,:,:),&
         G_qp_q(mm,:,:,:),G_q_qp(mm,:,:,:),&
         G_Np_N(mm,:,:),G_N_Np(mm,:,:),&
         G_zp_z(mm,:,:),G_z_zp(mm,:,:),&
         transl(ie,nn),transr(ie,nn),&
         G_D_RL(mm,:,:),G_D_LR(mm,:,:),&
         G_L(mm,:,:),G_R(mm,:,:),&
         hamd(mm,:,:,:),hamdv(mm,:,:,:),&
         hamld(mm,:,:),hamrd(mm,:,:),&
         haml(mm,:,:),hamr(mm,:,:))
 
 
end subroutine


subroutine get_G_full(ie,mm,Ham,norb,gl,gr,gml,gmr,sl,sr,&
                     E,GLeft,G_full,G_qp_q,G_q_qp,&
                     G_Np_N,G_N_Np,G_zp_z,G_z_zp,&
                     transl,transr,G_D_RL,G_D_LR,&
                     G_L,G_R)

    use util
    use config
    use surface , only : nl, nr, n_buffer_l, n_buffer_r
    use device  , only : nlayer,layer_list,layerstart,layerend

    implicit none

    integer(kind=4),intent(in) :: ie
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
    complex(kind=8),allocatable :: Es(:,:),Vcoup(:,:),&
                                   H(:,:)
    complex(kind=8),allocatable :: Gblock_NN(:,:),&
                                 Gblock_oN(:,:),Gblock_No(:,:)
    complex(kind=8),allocatable :: Gblock_NN0(:,:),&
                            Gblock_oN0(:,:),Gblock_No0(:,:) 
    integer(kind=4) :: i, device_start, j
    complex(kind=8),allocatable :: GRinv(:,:), GLinv(:,:)
    complex(kind=8) :: GGGGl(nl,nl), GGGGr(nr,nr)
    complex(kind=8) :: GLeft(:,:,:),G_full(:,:,:)
    complex(kind=8) :: G_qp_q(:,:,:),G_q_qp(:,:,:)
    complex(kind=8) :: G_Np_N(:,:),G_N_Np(:,:)
    complex(kind=8) :: G_zp_z(:,:),G_z_zp(:,:)
    real(kind=8) :: transl
    real(kind=8) :: transr
    complex(kind=8),intent(out) :: G_D_RL(nr,nl), G_D_LR(nl,nr)
    complex(kind=8),intent(out) :: G_L(nl,nl), G_R(nr,nr)


    complex(kind=8),allocatable :: eyed(:,:,:)
    complex(kind=8),allocatable :: Gblock_NN1(:,:),vcoup1(:,:)

    ! Algorithm ref: Modeling of nanoscale devices. Proceedings of the IEEE,
    ! 96(9), 1511-1550

    ! prepare self energy matrix due to butikker probe
    allocate(eyed(nlayer,layer_list(1),layer_list(1)))
    do i = 1,nlayer
        eyed(i,:,:) = sigma_w(i,ie)*eyemat(layer_list(1))
    end do
   
    ! index that marks the start of device region
    device_start = n_buffer_l+2*nl
    ! First, left sweep to get left-connected Green function
    allocate(Es(layer_list(1),layer_list(1)),&
            Vcoup(nl,layer_list(1)),&
            H(layer_list(1),layer_list(1)),&
            Gblock_NN(layer_list(1),layer_list(1)),&
            Gblock_oN(nl,layer_list(1)),&
            Gblock_NN0(layer_list(1),layer_list(1)),&
            Gblock_oN0(nl,layer_list(1)))
    allocate(Gblock_NN1(layer_list(nlayer),layer_list(nlayer)))
    ! Ham_LD
    Vcoup = Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
            device_start+layerstart(1):&
            device_start+layerend(1))
    H = Ham(device_start + layerstart(1):&
            device_start + layerend(1),&
            device_start + layerstart(1):&
            device_start + layerend(1))
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(layer_list(1))
    ! G_1, 1 with Buttiker self-energy
    Gblock_NN = inv(Es-H-matmul(matmul(transpose(dconjg(Vcoup)),gl(:,:)),Vcoup)-eyed(1,:,:))
    ! left-connected Green function
    GLeft(1,:,:) = Gblock_NN
    ! G_1, 1
    Gblock_oN = matmul(matmul(gl(:,:),Vcoup),Gblock_NN)
    Gblock_NN0 = Gblock_NN
    Gblock_oN0 = Gblock_oN

    deallocate(Vcoup)
    deallocate(Es)
    deallocate(H)
    deallocate(Gblock_NN)
    allocate(Vcoup(layer_list(1),layer_list(1)))
    allocate(Es(layer_list(1),layer_list(1)))
    allocate(H(layer_list(1),layer_list(1)))
    allocate(Gblock_NN(layer_list(1),layer_list(1)))

    do i = 2, nlayer 
        Vcoup = Ham(device_start+layerstart(i-1):&
                device_start+layerend(i-1),&
                device_start+layerstart(i):&
                device_start+layerend(i)) ! V_n_n+1
        Es = (E+(E*eta+eta0)*i_imag)*eyemat(layer_list(i))
        H = Ham(device_start+layerstart(i):&
                device_start+layerend(i),&
                device_start+layerstart(i):&
                device_start+layerend(i)) ! H_n+1_n+1
        Gblock_NN = inv(Es - H -&
        matmul(matmul(transpose(dconjg(Vcoup)),Gblock_NN0),&
        Vcoup)-eyed(i,:,:)) ! G_n+1_n+1
        if (i.eq. nlayer) then
            allocate(Vcoup1(layer_list(nlayer),nr))

            Vcoup1(:,:) = Ham(device_start+layerstart(nlayer):&
                         device_start+layerend(nlayer),&
                         norb-n_buffer_r-2*nr+1:&
                         norb-n_buffer_r-nr)
     
            Gblock_NN1 = inv(Es - H -&
                        matmul(matmul(transpose(dconjg(Vcoup)),Gblock_NN0),&
                            Vcoup)-eyed(i,:,:)-matmul(matmul(Vcoup1,gr),transpose(dconjg(Vcoup1))))
        end if

        ! left-connected Green function
        GLeft(i,:,:) = Gblock_NN
        Gblock_NN0 = Gblock_NN
        Gblock_oN0 = Gblock_oN ! G_0_n
        Gblock_oN = matmul(matmul(Gblock_oN0,Vcoup),Gblock_NN) ! G_0_n+1        
    end do
    ! Then, get the retarded full Green function G_n+1_n+1
    Gblock_NN0 = Gblock_NN

    deallocate(Es)
    allocate(Es(nr,nr))
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(nr)

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
    ! G_N+1_N+1
    allocate(Gblock_NN(nr,nr))
    Gblock_NN = inv(GRinv -&
                matmul(matmul(transpose(dconjg(Vcoup)),&
                Gblock_NN0),Vcoup))
    ! G_N_N
    G_full(nlayer,:,:) = GLeft(nlayer,:,:) + &
    matmul(matmul(matmul(matmul(GLeft(nlayer,:,:),Vcoup),Gblock_NN),&
    transpose(dconjg(Vcoup))),GLeft(nlayer,:,:)) ! Eq. B3
    GLeft(nlayer,:,:) = Gblock_NN1
    G_full(nlayer,:,:) = GLeft(nlayer,:,:)
    ! G_N+1_N
    G_Np_N = matmul(matmul(Gblock_NN,transpose(dconjg(Vcoup))),GLeft(nlayer,:,:))
    ! G_N_N+1
    G_N_Np = matmul(matmul(GLeft(nlayer,:,:),Vcoup),Gblock_NN)

    deallocate(Gblock_oN0)
    allocate(Gblock_oN0(nl,layer_list(nlayer)))
    Gblock_oN0 = Gblock_oN ! G_0_n

    deallocate(Gblock_oN)
    allocate(Gblock_oN(nl,nr))
    Gblock_oN = matmul(matmul(Gblock_oN0,Vcoup),Gblock_NN) ! G_0_n+1

    G_R =  Gblock_NN ! right green function in full green function
    G_D_LR = Gblock_oN
    
    ! Next, right sweep to obtain the full retarded Green function
    deallocate(Vcoup)
    allocate(Vcoup(layer_list(nlayer),layer_list(nlayer)))


    do i = nlayer-1,1,-1
        ! -A_q_q+1
        Vcoup = Ham(device_start+layerstart(i):&
                device_start+layerend(i),&
                device_start+layerstart(i+1):&
                device_start+layerend(i+1)) 
        ! G_q+1_q
        G_qp_q(i,:,:) = matmul(matmul(G_full(i+1,:,:),&
                        transpose(dconjg(Vcoup))),&
                        GLeft(i,:,:))
        ! G_q_q+1
        G_q_qp(i,:,:) = matmul(matmul(GLeft(i,:,:),&
                        Vcoup),&
                        G_full(i+1,:,:))
        ! G_q_q
        G_full(i,:,:) = GLeft(i,:,:)+&
                  matmul(matmul(GLeft(i,:,:),Vcoup)&
                 ,G_qp_q(i,:,:))
    end do

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
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(layer_list(nlayer))

    ! G_N, N
    Gblock_NN = inv(Es-H-&
                matmul(matmul(transpose(dconjg(Vcoup)),gr(:,:)),Vcoup)-eyed(nlayer,:,:))
    Gblock_No = matmul(matmul(gr(:,:),Vcoup),Gblock_NN)
    Gblock_NN0 = Gblock_NN
    Gblock_No0 = Gblock_No

    deallocate(Vcoup)
    allocate(Vcoup(layer_list(1),layer_list(1)))
    deallocate(Es)
    allocate(Es(layer_list(1),layer_list(1)))
    deallocate(H)
    allocate(H(layer_list(1),layer_list(1)))
    deallocate(Gblock_NN)
    allocate(Gblock_NN(layer_list(1),layer_list(1)))

    do i = nlayer,2,-1
        Vcoup = Ham(device_start+layerstart(i):&
                    device_start+layerend(i),&
                    device_start+layerstart(i-1):&
                    device_start+layerend(i-1)) ! V_n+1_n

        Es = (E + (E*eta+eta0)*i_imag) * eyemat(layer_list(i-1))

        H = Ham(device_start+layerstart(i-1):&
                device_start+layerend(i-1),&
                device_start+layerstart(i-1):&
                device_start+layerend(i-1)) ! H_n_n

        Gblock_NN = inv(Es-H-&
               matmul(matmul(transpose(dconjg(Vcoup)),&
               Gblock_NN0),Vcoup)-eyed(i-1,:,:)) ! G_n_n
   
        Gblock_NN0 = Gblock_NN

        Gblock_No0 = Gblock_No ! G_N+1_n+1

        Gblock_No = matmul(matmul(Gblock_No0,Vcoup)&
        ,Gblock_NN) ! G_N+1_n
    end do

    deallocate(Gblock_NN0)
    allocate(Gblock_NN0(layer_list(1),layer_list(1)))
    Gblock_NN0 = Gblock_NN

    deallocate(Es)
    allocate(Es(nl,nl))
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(nl)

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


end subroutine

subroutine get_G_full_tdb(ie,mm,gl,gr,gml,gmr,sl,sr,&
                     E,GLeft,G_full,G_qp_q,G_q_qp,&
                     G_Np_N,G_N_Np,G_zp_z,G_z_zp,&
                     transl,transr,G_D_RL,G_D_LR,&
                     G_L,G_R,&
                     hamd,hamdv,&
                     hamld,hamrd,&
                     haml,hamr)

    use util
    use config
    use surface , only : nl, nr, n_buffer_l, n_buffer_r
    use device  , only : nlayer,layer_list,layerstart,layerend

    implicit none

    integer(kind=4),intent(in) :: ie
    integer(kind=4),intent(in) :: mm
    complex(kind=8),intent(in) :: gl(nl,nl)
    complex(kind=8),intent(in) :: gr(nr,nr)
    complex(kind=8),intent(in) :: gml(nl,nl)
    complex(kind=8),intent(in) :: gmr(nr,nr)
    complex(kind=8),intent(in) :: sl(nl,nl)
    complex(kind=8),intent(in) :: sr(nr,nr)
    real(kind=8),intent(in) :: E
    complex(kind=8),allocatable :: Es(:,:),Vcoup(:,:),&
                                   H(:,:)
    complex(kind=8),allocatable :: Gblock_NN(:,:),&
                                 Gblock_oN(:,:),Gblock_No(:,:)
    complex(kind=8),allocatable :: Gblock_NN0(:,:),&
                            Gblock_oN0(:,:),Gblock_No0(:,:) 
    integer(kind=4) :: i, device_start, j
    complex(kind=8),allocatable :: GRinv(:,:), GLinv(:,:)
    complex(kind=8) :: GGGGl(nl,nl), GGGGr(nr,nr)
    complex(kind=8) :: GLeft(:,:,:),G_full(:,:,:)
    complex(kind=8) :: G_qp_q(:,:,:),G_q_qp(:,:,:)
    complex(kind=8) :: G_Np_N(:,:),G_N_Np(:,:)
    complex(kind=8) :: G_zp_z(:,:),G_z_zp(:,:)
    real(kind=8) :: transl
    real(kind=8) :: transr
    complex(kind=8),intent(out) :: G_D_RL(nr,nl), G_D_LR(nl,nr)
    complex(kind=8),intent(out) :: G_L(nl,nl), G_R(nr,nr)
    complex(kind=8),intent(in) :: hamd(:,:,:),hamdv(:,:,:)
    complex(kind=8),intent(in) :: hamld(:,:),hamrd(:,:)
    complex(kind=8),intent(in) :: haml(:,:),hamr(:,:)


    complex(kind=8),allocatable :: eyed(:,:,:)
    complex(kind=8),allocatable :: Gblock_NN1(:,:),vcoup1(:,:)

    ! Algorithm ref: Modeling of nanoscale devices. Proceedings of the IEEE,
    ! 96(9), 1511-1550

    ! prepare self energy matrix due to butikker probe
    allocate(eyed(nlayer,layer_list(1),layer_list(1)))
    do i = 1,nlayer
        eyed(i,:,:) = sigma_w(i,ie)*eyemat(layer_list(1))
    end do
    
    ! index that marks the start of device region
    device_start = n_buffer_l+2*nl
    ! First, left sweep to get left-connected Green function
    allocate(Es(layer_list(1),layer_list(1)),&
            Vcoup(nl,layer_list(1)),&
            H(layer_list(1),layer_list(1)),&
            Gblock_NN(layer_list(1),layer_list(1)),&
            Gblock_oN(nl,layer_list(1)),&
            Gblock_NN0(layer_list(1),layer_list(1)),&
            Gblock_oN0(nl,layer_list(1)))
    allocate(Gblock_NN1(layer_list(nlayer),layer_list(nlayer)))
    ! Ham_LD
    Vcoup = hamld 
    !Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
    !        device_start+layerstart(1):&
    !        device_start+layerend(1))
    H = hamd(1,:,:) 
    !Ham(device_start + layerstart(1):&
    !        device_start + layerend(1),&
    !        device_start + layerstart(1):&
    !        device_start + layerend(1))
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(layer_list(1))
    ! G_1, 1 with Buttiker self-energy
    ! gl is the advanced left surface GF
    Gblock_NN = inv(Es-H-matmul(matmul(transpose(dconjg(Vcoup)),gl(:,:)),Vcoup)-eyed(1,:,:))
    ! left-connected Green function
    GLeft(1,:,:) = Gblock_NN
    ! G_1, 1
    Gblock_oN = matmul(matmul(gl(:,:),Vcoup),Gblock_NN)
    Gblock_NN0 = Gblock_NN
    Gblock_oN0 = Gblock_oN

    deallocate(Vcoup)
    deallocate(Es)
    deallocate(H)
    deallocate(Gblock_NN)
    allocate(Vcoup(layer_list(1),layer_list(1)))
    allocate(Es(layer_list(1),layer_list(1)))
    allocate(H(layer_list(1),layer_list(1)))
    allocate(Gblock_NN(layer_list(1),layer_list(1)))

    do i = 2, nlayer 
        Vcoup = hamdv(i-1,:,:) 
        !Ham(device_start+layerstart(i-1):&
        !        device_start+layerend(i-1),&
        !        device_start+layerstart(i):&
        !        device_start+layerend(i)) ! V_n_n+1
        Es = (E+(E*eta+eta0)*i_imag)*eyemat(layer_list(i))
        H = hamd(i,:,:)
        !Ham(device_start+layerstart(i):&
        !        device_start+layerend(i),&
        !        device_start+layerstart(i):&
        !        device_start+layerend(i)) ! H_n+1_n+1
        Gblock_NN = inv(Es - H -&
        matmul(matmul(transpose(dconjg(Vcoup)),Gblock_NN0),&
        Vcoup)-eyed(i,:,:)) ! G_n+1_n+1
        if (i.eq. nlayer) then
            allocate(Vcoup1(layer_list(nlayer),nr))

            Vcoup1(:,:) = transpose(dconjg(hamrd))
            !Ham(device_start+layerstart(nlayer):&
            !             device_start+layerend(nlayer),&
            !             norb-n_buffer_r-2*nr+1:&
            !             norb-n_buffer_r-nr)
     
            Gblock_NN1 = inv(Es - H -&
                        matmul(matmul(transpose(dconjg(Vcoup)),Gblock_NN0),&
                            Vcoup)-eyed(i,:,:)-matmul(matmul(Vcoup1,gr),transpose(dconjg(Vcoup1))))
        end if

        ! left-connected Green function
        GLeft(i,:,:) = Gblock_NN
        Gblock_NN0 = Gblock_NN
        Gblock_oN0 = Gblock_oN ! G_0_n
        Gblock_oN = matmul(matmul(Gblock_oN0,Vcoup),Gblock_NN) ! G_0_n+1        
    end do
    ! Then, get the retarded full Green function G_n+1_n+1
    Gblock_NN0 = Gblock_NN

    deallocate(Es)
    allocate(Es(nr,nr))
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(nr)

    allocate(GRinv(nr,nr))
    GRinv = Es - hamr - sr(:,:) 

    deallocate(Vcoup)
    allocate(Vcoup(layer_list(nlayer),nr))

    Vcoup(:,:) = transpose(dconjg(hamrd))
    !Ham(device_start+layerstart(nlayer):&
    !                 device_start+layerend(nlayer),&
    !                 norb-n_buffer_r-2*nr+1:&
    !                 norb-n_buffer_r-nr)
    deallocate(Gblock_NN)
    ! G_N+1_N+1
    allocate(Gblock_NN(nr,nr))
    Gblock_NN = inv(GRinv -&
                matmul(matmul(transpose(dconjg(Vcoup)),&
                Gblock_NN0),Vcoup))
    ! G_N_N
    G_full(nlayer,:,:) = GLeft(nlayer,:,:) + &
    matmul(matmul(matmul(matmul(GLeft(nlayer,:,:),Vcoup),Gblock_NN),&
    transpose(dconjg(Vcoup))),GLeft(nlayer,:,:)) ! Eq. B3
    GLeft(nlayer,:,:) = Gblock_NN1
    G_full(nlayer,:,:) = GLeft(nlayer,:,:)
    ! G_N+1_N
    G_Np_N = matmul(matmul(Gblock_NN,transpose(dconjg(Vcoup))),GLeft(nlayer,:,:))
    ! G_N_N+1
    G_N_Np = matmul(matmul(GLeft(nlayer,:,:),Vcoup),Gblock_NN)

    deallocate(Gblock_oN0)
    allocate(Gblock_oN0(nl,layer_list(nlayer)))
    Gblock_oN0 = Gblock_oN ! G_0_n

    deallocate(Gblock_oN)
    allocate(Gblock_oN(nl,nr))
    Gblock_oN = matmul(matmul(Gblock_oN0,Vcoup),Gblock_NN) ! G_0_n+1

    G_R =  Gblock_NN ! right green function in full green function
    G_D_LR = Gblock_oN
    
    ! Next, right sweep to obtain the full retarded Green function
    deallocate(Vcoup)
    allocate(Vcoup(layer_list(nlayer),layer_list(nlayer)))


    do i = nlayer-1,1,-1
        ! -A_q_q+1
        Vcoup = hamdv(i,:,:)
        !Ham(device_start+layerstart(i):&
        !        device_start+layerend(i),&
        !        device_start+layerstart(i+1):&
        !        device_start+layerend(i+1)) 
        ! G_q+1_q
        G_qp_q(i,:,:) = matmul(matmul(G_full(i+1,:,:),&
                        transpose(dconjg(Vcoup))),&
                        GLeft(i,:,:))
        ! G_q_q+1
        G_q_qp(i,:,:) = matmul(matmul(GLeft(i,:,:),&
                        Vcoup),&
                        G_full(i+1,:,:))
        ! G_q_q
        G_full(i,:,:) = GLeft(i,:,:)+&
                        matmul(matmul(GLeft(i,:,:),Vcoup)&
                        ,G_qp_q(i,:,:))
    end do

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

    Vcoup = hamrd
    !Ham(norb-n_buffer_r-2*nr+1:&
    !            norb-n_buffer_r-nr,&
    !            device_start+layerstart(nlayer):&
    !            device_start+layerend(nlayer)) ! V_N+1_N

    H = hamd(nlayer,:,:)
    !Ham(device_start+layerstart(nlayer):&
    !        device_start+layerend(nlayer),&
    !        device_start+layerstart(nlayer):&
    !        device_start+layerend(nlayer)) ! H_N_N
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(layer_list(nlayer))

    ! G_N, N
    Gblock_NN = inv(Es-H-&
                matmul(matmul(transpose(dconjg(Vcoup)),gr(:,:)),Vcoup)-eyed(nlayer,:,:))
    Gblock_No = matmul(matmul(gr(:,:),Vcoup),Gblock_NN)

    Gblock_NN0 = Gblock_NN
    Gblock_No0 = Gblock_No

    deallocate(Vcoup)
    allocate(Vcoup(layer_list(1),layer_list(1)))
    deallocate(Es)
    allocate(Es(layer_list(1),layer_list(1)))
    deallocate(H)
    allocate(H(layer_list(1),layer_list(1)))
    deallocate(Gblock_NN)
    allocate(Gblock_NN(layer_list(1),layer_list(1)))

    do i = nlayer,2,-1
        Vcoup = transpose(dconjg(hamdv(i-1,:,:))) 
        !Ham(device_start+layerstart(i):&
        !            device_start+layerend(i),&
        !            device_start+layerstart(i-1):&
        !            device_start+layerend(i-1)) ! V_n+1_n

        Es = (E + (E*eta+eta0)*i_imag) * eyemat(layer_list(i-1))

        H = hamd(i-1,:,:) 
        !Ham(device_start+layerstart(i-1):&
        !        device_start+layerend(i-1),&
        !        device_start+layerstart(i-1):&
        !        device_start+layerend(i-1)) ! H_n_n

        Gblock_NN = inv(Es-H-&
               matmul(matmul(transpose(dconjg(Vcoup)),&
               Gblock_NN0),Vcoup)-eyed(i-1,:,:)) ! G_n_n
   
        Gblock_NN0 = Gblock_NN

        Gblock_No0 = Gblock_No ! G_N+1_n+1

        Gblock_No = matmul(matmul(Gblock_No0,Vcoup)&
        ,Gblock_NN) ! G_N+1_n
    end do

    deallocate(Gblock_NN0)
    allocate(Gblock_NN0(layer_list(1),layer_list(1)))
    Gblock_NN0 = Gblock_NN

    deallocate(Es)
    allocate(Es(nl,nl))
    Es = (E+(E*eta+eta0)*i_imag)*eyemat(nl)

    allocate(GLinv(nl,nl))
    GLinv = Es - haml - sl(:,:) ! G_1_1

    deallocate(Vcoup)
    allocate(Vcoup(nl,layer_list(1)))
    Vcoup = Hamld
    !Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
            !   device_start+layerstart(1):&
            !   device_start+layerend(1)) ! V_0_1

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

    if (ie.eq.1) then
        transl = 0.0d0
    end if
    if (ie.eq.1) then
        transr = 0.0d0
    end if

end subroutine


subroutine get_G_lesser_mpi(ie,ik_core,Ham,norb,gl,gr,gml,gmr,&
               sl,sr,E_k,GLeft,G_full,&
               G_qp_q,G_q_qp,G_N_Np,G_Np_N,&
               G_zp_z,G_z_zp,GnL,Gn,dGnLdT,dGndT,&
               flux_i,Gn_i,A_i)

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
    complex(kind=8)            :: GLeft(:,:,:,:)
    complex(kind=8)            :: G_full(:,:,:,:)
    complex(kind=8)            :: G_qp_q(:,:,:,:),G_q_qp(:,:,:,:)
    complex(kind=8)            :: G_Np_N(:,:,:),G_N_Np(:,:,:)
    complex(kind=8)            :: G_zp_z(:,:,:),G_z_zp(:,:,:)
    complex(kind=8)            :: GnL(:,:,:,:),Gn(:,:,:,:)
    complex(kind=8)            :: dGnLdT(:,:,:,:,:),dGndT(:,:,:,:,:)
    complex(kind=8)         :: flux_i(:,:,:),Gn_i(:,:,:),A_i(:,:,:)

    mm=myid+1
    nn = k_start_core(mm)+ik_core-1

    call get_G_lesser(ie,mm,nn,Ham(mm,:,:),norb,&
         transpose(dconjg(gl(mm,:,:))),gr(mm,:,:),&
         -gml(mm,:,:),gmr(mm,:,:),&
         transpose(dconjg(sl(mm,:,:))),sr(mm,:,:),E_k(ie,nn),&
         GLeft(mm,:,:,:),G_full(mm,:,:,:),&
         G_qp_q(mm,:,:,:),G_q_qp(mm,:,:,:),&
         G_Np_N(mm,:,:),G_N_Np(mm,:,:),&
         G_zp_z(mm,:,:),G_z_zp(mm,:,:),&
         GnL(mm,:,:,:),Gn(mm,:,:,:),&
         dGnLdT(mm,:,:,:,:),dGndT(mm,:,:,:,:),&
         flux_i(nn,ie,:),Gn_i(nn,ie,:),A_i(nn,ie,:))
 
 
end subroutine

subroutine get_G_lesser_tdb_mpi(ie,ik_core,gl,gr,gml,gmr,&
               sl,sr,E_k,GLeft,G_full,&
               G_qp_q,G_q_qp,G_N_Np,G_Np_N,&
               G_zp_z,G_z_zp,GnL,Gn,dGnLdT,dGndT,&
               flux_i,Gn_i,A_i,&
               hamd,hamdv,hamld,hamrd,&
               flux_left,flux_right)

    use config
    use surface, only : nl, nr

    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: ie, ik_core
    integer(kind=4) :: mm, nn
    complex(kind=8),intent(in) :: gl(numprocs,nl,nl),gml(numprocs,nl,nl)
    complex(kind=8),intent(in) :: gr(numprocs,nr,nr),gmr(numprocs,nr,nr)
    complex(kind=8),intent(in) :: sl(numprocs,nl,nl),sr(numprocs,nr,nr)
    real(kind=8),intent(in) :: E_k(ne,nk)
    complex(kind=8)            :: GLeft(:,:,:,:)
    complex(kind=8)            :: G_full(:,:,:,:)
    complex(kind=8)            :: G_qp_q(:,:,:,:),G_q_qp(:,:,:,:)
    complex(kind=8)            :: G_Np_N(:,:,:),G_N_Np(:,:,:)
    complex(kind=8)            :: G_zp_z(:,:,:),G_z_zp(:,:,:)
    complex(kind=8)            :: GnL(:,:,:,:),Gn(:,:,:,:)
    complex(kind=8)            :: dGnLdT(:,:,:,:,:),dGndT(:,:,:,:,:)
    complex(kind=8)         :: flux_i(:,:,:),Gn_i(:,:,:),A_i(:,:,:)
    complex(kind=8)         :: flux_left(:,:),flux_right(:,:)
    complex(kind=8)          :: hamd(:,:,:,:),hamdv(:,:,:,:)
    complex(kind=8)          :: hamld(:,:,:),hamrd(:,:,:)

    mm=myid+1
    nn = k_start_core(mm)+ik_core-1

    call get_G_lesser_tdb(ie,mm,nn,&
         transpose(dconjg(gl(mm,:,:))),gr(mm,:,:),&
         -gml(mm,:,:),gmr(mm,:,:),&
         transpose(dconjg(sl(mm,:,:))),sr(mm,:,:),E_k(ie,nn),&
         GLeft(mm,:,:,:),G_full(mm,:,:,:),&
         G_qp_q(mm,:,:,:),G_q_qp(mm,:,:,:),&
         G_Np_N(mm,:,:),G_N_Np(mm,:,:),&
         G_zp_z(mm,:,:),G_z_zp(mm,:,:),&
         GnL(mm,:,:,:),Gn(mm,:,:,:),&
         dGnLdT(mm,:,:,:,:),dGndT(mm,:,:,:,:),&
         flux_i(nn,ie,:),Gn_i(nn,ie,:),A_i(nn,ie,:),&
         hamd(mm,:,:,:),hamdv(mm,:,:,:),&
         hamld(mm,:,:),hamrd(mm,:,:),&
         flux_left(nn,ie),flux_right(nn,ie))
 
 
end subroutine


subroutine get_G_lesser(ie,mm,nn,Ham,norb,gl,gr,gml,gmr,sl,sr,&
                        E,GLeft,G_full,G_qp_q,G_q_qp,&
                        G_Np_N,G_N_Np,G_zp_z,G_z_zp,&
                        GnL,Gn,dGnLdT,dGndT,&
                        flux_i,Gn_i,A_i_sum)

    use util
    use config
    use surface , only : nl, nr, n_buffer_l, n_buffer_r
    use device  , only : nlayer,layer_list,layerstart,layerend

    implicit none

    integer(kind=4),intent(in) :: ie
    integer(kind=4),intent(in) :: mm
    integer(kind=4),intent(in) :: nn
    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(in) :: Ham(norb,norb)
    complex(kind=8),intent(in) :: gl(nl,nl)
    complex(kind=8),intent(in) :: gr(nr,nr)
    complex(kind=8),intent(in) :: gml(nl,nl)
    complex(kind=8),intent(in) :: gmr(nr,nr)
    complex(kind=8),intent(in) :: sl(nl,nl)
    complex(kind=8),intent(in) :: sr(nr,nr)
    real(kind=8),intent(in) :: E
    complex(kind=8),allocatable :: Es(:,:),Vcoup(:,:),&
                                   H(:,:),Vcoup1(:,:)
    complex(kind=8),allocatable :: Gblock_NN(:,:),&
                                 Gblock_oN(:,:),Gblock_No(:,:)
    complex(kind=8),allocatable :: Gblock_NN0(:,:),&
                            Gblock_oN0(:,:),Gblock_No0(:,:) 
    integer(kind=4) :: i, device_start, j
    complex(kind=8),allocatable :: GRinv(:,:), GLinv(:,:)
    complex(kind=8) :: GLeft(:,:,:),G_full(:,:,:)
    complex(kind=8) :: GnL(:,:,:),Gn(:,:,:)
    complex(kind=8) :: G_qp_q(:,:,:),G_q_qp(:,:,:)
    complex(kind=8) :: G_Np_N(:,:),G_N_Np(:,:)
    complex(kind=8) :: G_zp_z(:,:),G_z_zp(:,:)
    complex(kind=8) :: dGnLdT(:,:,:,:),dGndT(:,:,:,:)
    complex(kind=8) :: flux_i(:),Gn_i(:),A_i_sum(:)
    complex(kind=8),allocatable :: eyed(:,:,:)
    complex(kind=8),allocatable :: deyedt(:,:,:)
    complex(kind=8),allocatable :: sigma_in(:,:,:)
    complex(kind=8),allocatable :: gamma_i(:,:,:)
    complex(kind=8),allocatable :: dgamma_idt(:,:,:)
    complex(kind=8),allocatable :: A_i(:,:,:)
    complex(kind=8),allocatable :: sigma_first(:,:),sigma_prev(:,:),sigma_last(:,:)
    complex(kind=8),allocatable :: gamma_first(:,:),gamma_last(:,:)
    complex(kind=8),allocatable :: dsigma_dt(:,:)
    real(kind=8)    :: freq

    ! frequency
    freq = sqrt(E)

    ! prepare self energy matrix due to butikker probe
    allocate(eyed(nlayer,layer_list(1),layer_list(1)))
    allocate(deyedt(nlayer,layer_list(1),layer_list(1)))
    do i = 1,nlayer
        eyed(i,:,:) = sigma_w(i,ie)*eyemat(layer_list(1))
        deyedt(i,:,:) = dsigma_wdt(i,ie)*eyemat(layer_list(1))
    end do
    ! prepare self energy
    allocate(sigma_in(nlayer,layer_list(1),layer_list(1)))
    sigma_in = 0.0d0
    ! prepare escape rate
    allocate(gamma_i(nlayer,layer_list(1),layer_list(1)))
    gamma_i = 0.0d0
    allocate(dgamma_idt(nlayer,layer_list(1),layer_list(1)))
    dgamma_idt = 0.0d0
    ! prepare spectral function
    allocate(A_i(nlayer,layer_list(1),layer_list(1)))
    A_i = 0.0d0
    A_i_sum = 0.0d0
    Gn_i = 0.0d0
    ! index that marks the start of device region
    device_start = n_buffer_l+2*nl
    ! First, left sweep to get left-connected Green function
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
    ! self-energy due to left lead in the first block
    allocate(sigma_first(layer_list(1),layer_list(1)))
    allocate(gamma_first(layer_list(1),layer_list(1)))
    sigma_first = matmul(matmul(transpose(dconjg(Vcoup)),gl),&
                             Vcoup) 
    gamma_first = i_imag*(sigma_first-transpose(dconjg(sigma_first)))
    ! spectral function
    A_i(1,:,:) = i_imag*(G_full(1,:,:)-transpose(dconjg(G_full(1,:,:))))
    ! lesser self-energy (in-scattering) for first block
    gamma_i(1,:,:) =  i_imag*(eyed(1,:,:)-transpose(dconjg(eyed(1,:,:))))
    dgamma_idt(1,:,:) =  i_imag*(deyedt(1,:,:)-transpose(dconjg(deyedt(1,:,:))))
    sigma_in(1,:,:) = gamma_first*bose_einstein(freq,T_left)&
                      +gamma_i(1,:,:)&
                      *bose_einstein(freq,temperature_i(1))
    ! left-connected lesser Green function in first block
    GnL(1,:,:) = matmul(matmul(GLeft(1,:,:),sigma_in(1,:,:)),&
                 transpose(dconjg(GLeft(1,:,:))))    
    ! dsigma_dt for the first block
    allocate(dsigma_dt(layer_list(1),layer_list(1)))
    dsigma_dt(:,:) =  gamma_i(1,:,:)&
                      *dbosedT(freq,temperature_i(1))&
                    + dgamma_idt(1,:,:)&
                      *bose_einstein(freq,temperature_i(1))
    ! temperature derivative of left-connected lesser GF
    dGnLdT(1,1,:,:) = matmul(matmul(GLeft(1,:,:),dsigma_dt),&
                      transpose(dconjg(GLeft(1,:,:))))
    !do j = 2,nlayer
    dGnLdT(1,2:nlayer,:,:) = 0.0d0 ! dG_nL_1_1/dT_j
    !end do

    deallocate(Vcoup)
    allocate(Vcoup(layer_list(1),layer_list(1)))
    allocate(sigma_prev(layer_list(1),layer_list(1)))

    do i = 2,nlayer-1
        Vcoup = Ham(device_start+layerstart(i-1):&
                device_start+layerend(i-1),&
                device_start+layerstart(i):&
                device_start+layerend(i)) ! V_n_n+1
        ! sigma_in_i+1_i+1
        sigma_prev = matmul(matmul(transpose(dconjg(Vcoup)),GnL(i-1,:,:)),&
                             Vcoup)

        ! spectral function A_i
        A_i(i,:,:) = i_imag*(G_full(i,:,:)-transpose(dconjg(G_full(i,:,:))))
 
        ! Gamma_i
        gamma_i(i,:,:) = i_imag*(eyed(i,:,:)-transpose(dconjg(eyed(i,:,:))))
        dgamma_idt(i,:,:) = i_imag*(deyedt(i,:,:)-transpose(dconjg(deyedt(i,:,:))))
        ! Sigma_in_i+1_i+1
        sigma_in(i,:,:) = gamma_i(i,:,:)&
                         *bose_einstein(freq,temperature_i(i))
        ! GnL_i+1_i+1
        GnL(i,:,:) = sigma_in(i,:,:) + sigma_prev
        GnL(i,:,:) = matmul(matmul(GLeft(i,:,:),GnL(i,:,:)),&
                     transpose(dconjg(GLeft(i,:,:)))) ! Eq. B8             
        do j = 1,i-1
            ! j<i+1
            ! d sigma_in_i+1_i+1 dT_j
            dsigma_dt(:,:) = matmul(matmul(transpose(dconjg(Vcoup)),dGnLdT(i-1,j,:,:)),&
                        Vcoup)
            ! dGnL_i+1_i+1 dT_j
            dGnLdT(i,j,:,:) =  matmul(matmul(GLeft(i,:,:),dsigma_dt),&
                               transpose(dconjg(GLeft(i,:,:))))
        end do
        ! j = i+1
        dsigma_dt(:,:) = gamma_i(i,:,:)&
                         *dbosedT(freq,temperature_i(i))&
                        +dgamma_idt(i,:,:)&
                         *bose_einstein(freq,temperature_i(i))
        ! dGnL_i+1_i+1 dT_i+1
        dGnLdT(i,i,:,:) = matmul(matmul(GLeft(i,:,:),dsigma_dt),&
                          transpose(dconjg(GLeft(i,:,:))))
        ! j>i+1
        dGnLdT(i,i+1:nlayer,:,:) = 0.0d0
    end do
    ! GnL_N_N
    Vcoup = Ham(device_start+layerstart(nlayer-1):&
                device_start+layerend(nlayer-1),&
                device_start+layerstart(nlayer):&
                device_start+layerend(nlayer)) ! V_n_n+1
 
    sigma_prev = matmul(matmul(transpose(dconjg(Vcoup)),GnL(nlayer-1,:,:)),&
                                 Vcoup)

    allocate(sigma_last(layer_list(nlayer),layer_list(nlayer)),&
             Vcoup1(layer_list(nlayer),nr))
    allocate(gamma_last(layer_list(nlayer),layer_list(nlayer)))
    Vcoup1(:,:) = Ham(device_start+layerstart(nlayer):&
                      device_start+layerend(nlayer),&
                      norb-n_buffer_r-2*nr+1:&
                      norb-n_buffer_r-nr)
    ! self energy due to right lead
    sigma_last = matmul(matmul(Vcoup1,gr),&
                 transpose(dconjg(Vcoup1)))
    gamma_last = i_imag*(sigma_last - transpose(dconjg(sigma_last)))

    ! spectral function A_N
    A_i(nlayer,:,:) = i_imag*(G_full(nlayer,:,:)-transpose(dconjg(G_full(nlayer,:,:))))
 
    ! Gamma_i
    gamma_i(nlayer,:,:) = i_imag*(eyed(nlayer,:,:)-transpose(dconjg(eyed(nlayer,:,:))))
    dgamma_idt(nlayer,:,:) = i_imag*(deyedt(nlayer,:,:)-transpose(dconjg(deyedt(nlayer,:,:))))


    sigma_in(nlayer,:,:) =  gamma_i(nlayer,:,:)&
                           *bose_einstein(freq,temperature_i(nlayer))&
                           +gamma_last&
                           *bose_einstein(freq,T_right)
    
    GnL(nlayer,:,:) = sigma_in(nlayer,:,:) + sigma_prev
    ! G_full(nlayer,:,:) is not GLeft(nlayer) in current implementation
    GnL(nlayer,:,:) = matmul(matmul(GLeft(nlayer,:,:),GnL(nlayer,:,:)),&
                      transpose(dconjg(GLeft(nlayer,:,:))))
    do j = 1,nlayer-1
        ! d sigma_in_N_N dT_j
        dsigma_dt(:,:) = matmul(matmul(transpose(dconjg(Vcoup)),dGnLdT(nlayer-1,j,:,:)),&
                         Vcoup)
        ! dGnL_N_N dT_j
        dGnLdT(nlayer,j,:,:) =  matmul(matmul(GLeft(nlayer,:,:),dsigma_dt),&
                                transpose(dconjg(GLeft(nlayer,:,:))))
    end do
    ! j = N
    dsigma_dt(:,:) = gamma_i(nlayer,:,:)&
                     *dbosedT(freq,temperature_i(nlayer))&
                    +dgamma_idt(nlayer,:,:)&
                     *bose_einstein(freq,temperature_i(nlayer))
    ! dGnL_N_N dT_N
    dGnLdT(nlayer,nlayer,:,:) = matmul(matmul(GLeft(nlayer,:,:),dsigma_dt),&
                          transpose(dconjg(GLeft(nlayer,:,:))))
 
 
    ! Now, we use left connected lesser GF to compute full lesser GF
    Gn(nlayer,:,:) = GnL(nlayer,:,:)
    dGndT(nlayer,:,:,:) = dGnLdT(nlayer,:,:,:)


    deallocate(Vcoup)
    allocate(Vcoup(layer_list(nlayer),layer_list(nlayer)))
    do i = nlayer-1,1,-1
        Vcoup = Ham(device_start+layerstart(i):&
                    device_start+layerend(i),&
                    device_start+layerstart(i+1):&
                    device_start+layerend(i+1)) ! V_n_n+1
        
        Gn(i,:,:) = GnL(i,:,:) + matmul(matmul(matmul(matmul(GLeft(i,:,:),&
                    Vcoup),Gn(i+1,:,:)),transpose(dconjg(Vcoup))),&
                    transpose(dconjg(GLeft(i,:,:)))) &
                    + matmul(matmul(GnL(i,:,:),Vcoup),&
                    transpose(dconjg(G_q_qp(i,:,:))))&
                    + matmul(matmul(G_q_qp(i,:,:),transpose(dconjg(Vcoup))),&
                    GnL(i,:,:)) 
        do j = 1,nlayer
            dGndT(i,j,:,:) = dGnLdT(i,j,:,:) + matmul(matmul(matmul(matmul(GLeft(i,:,:),&
                   Vcoup),dGndT(i+1,j,:,:)),transpose(dconjg(Vcoup))),&
                   transpose(dconjg(GLeft(i,:,:)))) &
                 + matmul(matmul(dGnLdT(i,j,:,:),Vcoup),&
                   transpose(dconjg(G_q_qp(i,:,:))))&
                 + matmul(matmul(G_q_qp(i,:,:),transpose(dconjg(Vcoup))),&
                   dGnLdT(i,j,:,:))
        end do
    end do
   
    ! zero frequency limit
    if (ie.eq.1) then
        Gn(:,:,:) = 0.0d0
        dGndT(:,:,:,:) = 0.0d0
        A_i(:,:,:) = 0.0d0
    end if

    do i = 1,nlayer        
        flux_i(i) = -1.0d0*trace(matmul(gamma_i(i,:,:),Gn(i,:,:)))&
                    + trace(matmul(gamma_i(i,:,:),A_i(i,:,:)))&
                    * bose_einstein(freq,temperature_i(i))
        do j = 1,nlayer
            jacobian(nn,ie,i,j) = -1.0d0*trace(matmul(gamma_i(i,:,:),dGndT(i,j,:,:)))
            if (i.eq.j) then
                jacobian(nn,ie,i,j) = jacobian(nn,ie,i,j) + &
                                  trace(matmul(gamma_i(i,:,:),A_i(i,:,:))&
                                  *dbosedT(freq,temperature_i(j))) + &
                                  trace(matmul(dgamma_idt(i,:,:),A_i(i,:,:))&
                                  *bose_einstein(freq,temperature_i(j)))
            end if
        end do
    end do

    do i = 1,nlayer
        A_i_sum(i) =  trace(A_i(i,:,:))
        Gn_i(i) =  trace(Gn(i,:,:))
    end do


    deallocate(sigma_first,sigma_last,sigma_prev,Vcoup,Vcoup1,sigma_in)

end subroutine

subroutine get_G_lesser_tdb(ie,mm,nn,gl,gr,gml,gmr,sl,sr,&
                        E,GLeft,G_full,G_qp_q,G_q_qp,&
                        G_Np_N,G_N_Np,G_zp_z,G_z_zp,&
                        GnL,Gn,dGnLdT,dGndT,&
                        flux_i,Gn_i,A_i_sum,&
                        hamd,hamdv,hamld,hamrd,&
                        flux_left,flux_right)

    use util
    use config
    use surface , only : nl, nr, n_buffer_l, n_buffer_r
    use device  , only : nlayer,layer_list,layerstart,layerend

    implicit none

    integer(kind=4),intent(in) :: ie
    integer(kind=4),intent(in) :: mm
    integer(kind=4),intent(in) :: nn
    complex(kind=8),intent(in) :: gl(nl,nl)
    complex(kind=8),intent(in) :: gr(nr,nr)
    complex(kind=8),intent(in) :: gml(nl,nl)
    complex(kind=8),intent(in) :: gmr(nr,nr)
    complex(kind=8),intent(in) :: sl(nl,nl)
    complex(kind=8),intent(in) :: sr(nr,nr)
    real(kind=8),intent(in) :: E
    complex(kind=8),allocatable :: Es(:,:),Vcoup(:,:),&
                                   H(:,:),Vcoup1(:,:)
    complex(kind=8),allocatable :: Gblock_NN(:,:),&
                                 Gblock_oN(:,:),Gblock_No(:,:)
    complex(kind=8),allocatable :: Gblock_NN0(:,:),&
                            Gblock_oN0(:,:),Gblock_No0(:,:) 
    integer(kind=4) :: i, device_start, j
    complex(kind=8),allocatable :: GRinv(:,:), GLinv(:,:)
    complex(kind=8) :: GLeft(:,:,:),G_full(:,:,:)
    complex(kind=8) :: GnL(:,:,:),Gn(:,:,:)
    complex(kind=8) :: G_qp_q(:,:,:),G_q_qp(:,:,:)
    complex(kind=8) :: G_Np_N(:,:),G_N_Np(:,:)
    complex(kind=8) :: G_zp_z(:,:),G_z_zp(:,:)
    complex(kind=8) :: dGnLdT(:,:,:,:),dGndT(:,:,:,:)
    complex(kind=8) :: flux_i(:),Gn_i(:),A_i_sum(:)
    complex(kind=8) :: flux_left,flux_right
    complex(kind=8),allocatable :: eyed(:,:,:)
    complex(kind=8),allocatable :: deyedt(:,:,:)
    complex(kind=8),allocatable :: sigma_in(:,:,:)
    complex(kind=8),allocatable :: gamma_i(:,:,:)
    complex(kind=8),allocatable :: dgamma_idt(:,:,:)
    complex(kind=8),allocatable :: A_i(:,:,:)
    complex(kind=8),allocatable :: sigma_first(:,:),sigma_prev(:,:),sigma_last(:,:)
    complex(kind=8),allocatable :: gamma_first(:,:),gamma_last(:,:)
    complex(kind=8),allocatable :: dsigma_dt(:,:)
    real(kind=8)    :: freq
    complex(kind=8) :: hamd(:,:,:),hamdv(:,:,:)
    complex(kind=8) :: hamld(:,:),hamrd(:,:)

    ! frequency
    freq = sqrt(E)

    ! prepare self energy matrix due to butikker probe
    allocate(eyed(nlayer,layer_list(1),layer_list(1)))
    allocate(deyedt(nlayer,layer_list(1),layer_list(1)))
    do i = 1,nlayer
        eyed(i,:,:) = sigma_w(i,ie)*eyemat(layer_list(1))
        deyedt(i,:,:) = dsigma_wdt(i,ie)*eyemat(layer_list(1))
    end do
    ! prepare self energy
    allocate(sigma_in(nlayer,layer_list(1),layer_list(1)))
    sigma_in = 0.0d0
    ! prepare escape rate
    allocate(gamma_i(nlayer,layer_list(1),layer_list(1)))
    gamma_i = 0.0d0
    allocate(dgamma_idt(nlayer,layer_list(1),layer_list(1)))
    dgamma_idt = 0.0d0 
    ! prepare spectral function
    allocate(A_i(nlayer,layer_list(1),layer_list(1)))
    A_i = 0.0d0
    A_i_sum = 0.0d0
    Gn_i = 0.0d0
    ! index that marks the start of device region
    device_start = n_buffer_l+2*nl
    ! First, left sweep to get left-connected Green function
    allocate(Es(layer_list(1),layer_list(1)),&
             Vcoup(nl,layer_list(1)),&
             H(layer_list(1),layer_list(1)),&
             Gblock_NN(layer_list(1),layer_list(1)),&
             Gblock_oN(nl,layer_list(1)),&
             Gblock_NN0(layer_list(1),layer_list(1)),&
             Gblock_oN0(nl,layer_list(1)))
    ! Ham_LD
    Vcoup = hamld
    !Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
    !            device_start+layerstart(1):&
    !            device_start+layerend(1))
    ! self-energy due to left lead in the first block
    allocate(sigma_first(layer_list(1),layer_list(1)))
    allocate(gamma_first(layer_list(1),layer_list(1)))
    sigma_first = matmul(matmul(transpose(dconjg(Vcoup)),gl),&
                             Vcoup) 
    gamma_first = i_imag*(sigma_first-transpose(dconjg(sigma_first)))
    ! spectral function
    A_i(1,:,:) = i_imag*(G_full(1,:,:)-transpose(dconjg(G_full(1,:,:))))
    ! lesser self-energy (in-scattering) for first block
    gamma_i(1,:,:) =  i_imag*(eyed(1,:,:)-transpose(dconjg(eyed(1,:,:))))
    dgamma_idt(1,:,:) =  i_imag*(deyedt(1,:,:)-transpose(dconjg(deyedt(1,:,:))))
    ! left flux
    flux_left = trace(matmul(gamma_first,A_i(1,:,:)))*bose_einstein(freq,T_left)

    sigma_in(1,:,:) = gamma_first*bose_einstein(freq,T_left)&
                      +gamma_i(1,:,:)&
                      *bose_einstein(freq,temperature_i(1))
    ! left-connected lesser Green function in first block
    GnL(1,:,:) = matmul(matmul(GLeft(1,:,:),sigma_in(1,:,:)),&
                 transpose(dconjg(GLeft(1,:,:))))    
    ! dsigma_dt for the first block
    allocate(dsigma_dt(layer_list(1),layer_list(1)))
    dsigma_dt(:,:) =  gamma_i(1,:,:)&
                      *dbosedT(freq,temperature_i(1))&
                    + dgamma_idt(1,:,:)&
                      *bose_einstein(freq,temperature_i(1))
    ! temperature derivative of left-connected lesser GF
    dGnLdT(1,1,:,:) = matmul(matmul(GLeft(1,:,:),dsigma_dt),&
                      transpose(dconjg(GLeft(1,:,:))))
    !do j = 2,nlayer
    dGnLdT(1,2:nlayer,:,:) = 0.0d0 ! dG_nL_1_1/dT_j
    !end do

    deallocate(Vcoup)
    allocate(Vcoup(layer_list(1),layer_list(1)))
    allocate(sigma_prev(layer_list(1),layer_list(1)))

    do i = 2,nlayer-1
        Vcoup = hamdv(i-1,:,:)
        !Ham(device_start+layerstart(i-1):&
        !        device_start+layerend(i-1),&
        !        device_start+layerstart(i):&
        !        device_start+layerend(i)) ! V_n_n+1
        ! sigma_in_i+1_i+1
        sigma_prev = matmul(matmul(transpose(dconjg(Vcoup)),GnL(i-1,:,:)),&
                             Vcoup)

        ! spectral function A_i
        A_i(i,:,:) = i_imag*(G_full(i,:,:)-transpose(dconjg(G_full(i,:,:))))
 
        ! Gamma_i
        gamma_i(i,:,:) = i_imag*(eyed(i,:,:)-transpose(dconjg(eyed(i,:,:))))
        dgamma_idt(i,:,:) = i_imag*(deyedt(i,:,:)-transpose(dconjg(deyedt(i,:,:))))
        ! Sigma_in_i+1_i+1
        sigma_in(i,:,:) = gamma_i(i,:,:)&
                         *bose_einstein(freq,temperature_i(i))
        ! GnL_i+1_i+1
        GnL(i,:,:) = sigma_in(i,:,:) + sigma_prev
        GnL(i,:,:) = matmul(matmul(GLeft(i,:,:),GnL(i,:,:)),&
                     transpose(dconjg(GLeft(i,:,:)))) ! Eq. B8             
        do j = 1,i-1
            ! j<i+1
            ! d sigma_in_i+1_i+1 dT_j
            dsigma_dt(:,:) = matmul(matmul(transpose(dconjg(Vcoup)),dGnLdT(i-1,j,:,:)),&
                        Vcoup)
            ! dGnL_i+1_i+1 dT_j
            dGnLdT(i,j,:,:) =  matmul(matmul(GLeft(i,:,:),dsigma_dt),&
                               transpose(dconjg(GLeft(i,:,:))))
        end do
        ! j = i+1
        dsigma_dt(:,:) = gamma_i(i,:,:)&
                         *dbosedT(freq,temperature_i(i))&
                        +dgamma_idt(i,:,:)&
                         *bose_einstein(freq,temperature_i(i))
        ! dGnL_i+1_i+1 dT_i+1
        dGnLdT(i,i,:,:) = matmul(matmul(GLeft(i,:,:),dsigma_dt),&
                          transpose(dconjg(GLeft(i,:,:))))
        ! j>i+1
        dGnLdT(i,i+1:nlayer,:,:) = 0.0d0
    end do
    ! GnL_N_N
    Vcoup = hamdv(nlayer-1,:,:)
    !Ham(device_start+layerstart(nlayer-1):&
            !    device_start+layerend(nlayer-1),&
            !    device_start+layerstart(nlayer):&
            !    device_start+layerend(nlayer)) ! V_n_n+1
 
    sigma_prev = matmul(matmul(transpose(dconjg(Vcoup)),GnL(nlayer-1,:,:)),&
                                 Vcoup)

    allocate(sigma_last(layer_list(nlayer),layer_list(nlayer)),&
             Vcoup1(layer_list(nlayer),nr))
    allocate(gamma_last(layer_list(nlayer),layer_list(nlayer)))
    Vcoup1(:,:) = transpose(dconjg(Hamrd))
    ! Ham(device_start+layerstart(nlayer):&
                   !   device_start+layerend(nlayer),&
                   !   norb-n_buffer_r-2*nr+1:&
                   !   norb-n_buffer_r-nr)
    ! self energy due to right lead
    sigma_last = matmul(matmul(Vcoup1,gr),&
                 transpose(dconjg(Vcoup1)))
    gamma_last = i_imag*(sigma_last - transpose(dconjg(sigma_last)))

    ! spectral function A_N
    A_i(nlayer,:,:) = i_imag*(G_full(nlayer,:,:)-transpose(dconjg(G_full(nlayer,:,:))))
    ! right flux
    flux_right = trace(matmul(gamma_last,A_i(nlayer,:,:)))*bose_einstein(freq,T_right)
 
    ! Gamma_i
    gamma_i(nlayer,:,:) = i_imag*(eyed(nlayer,:,:)-transpose(dconjg(eyed(nlayer,:,:))))
    dgamma_idt(nlayer,:,:) = i_imag*(deyedt(nlayer,:,:)-transpose(dconjg(deyedt(nlayer,:,:))))
    sigma_in(nlayer,:,:) =  gamma_i(nlayer,:,:)&
                           *bose_einstein(freq,temperature_i(nlayer))&
                           +gamma_last&
                           *bose_einstein(freq,T_right)
    
    GnL(nlayer,:,:) = sigma_in(nlayer,:,:) + sigma_prev
    ! G_full(nlayer,:,:) is not GLeft(nlayer) in current implementation
    GnL(nlayer,:,:) = matmul(matmul(GLeft(nlayer,:,:),GnL(nlayer,:,:)),&
                      transpose(dconjg(GLeft(nlayer,:,:))))
    do j = 1,nlayer-1
        ! d sigma_in_N_N dT_j
        dsigma_dt(:,:) = matmul(matmul(transpose(dconjg(Vcoup)),dGnLdT(nlayer-1,j,:,:)),&
                         Vcoup)
        ! dGnL_N_N dT_j
        dGnLdT(nlayer,j,:,:) =  matmul(matmul(GLeft(nlayer,:,:),dsigma_dt),&
                                transpose(dconjg(GLeft(nlayer,:,:))))
    end do
    ! j = N
    dsigma_dt(:,:) = gamma_i(nlayer,:,:)&
                     *dbosedT(freq,temperature_i(nlayer))&
                    +dgamma_idt(nlayer,:,:)&
                     *bose_einstein(freq,temperature_i(nlayer))
    ! dGnL_N_N dT_N
    dGnLdT(nlayer,nlayer,:,:) = matmul(matmul(GLeft(nlayer,:,:),dsigma_dt),&
                          transpose(dconjg(GLeft(nlayer,:,:))))
 
 
    ! Now, we use left connected lesser GF to compute full lesser GF
    Gn(nlayer,:,:) = GnL(nlayer,:,:)
    ! right flux
    flux_right = (flux_right - trace(matmul(gamma_last,Gn(nlayer,:,:))))*freq

    dGndT(nlayer,:,:,:) = dGnLdT(nlayer,:,:,:)


    deallocate(Vcoup)
    allocate(Vcoup(layer_list(nlayer),layer_list(nlayer)))
    do i = nlayer-1,1,-1
        Vcoup = Hamdv(i,:,:)
        ! Ham(device_start+layerstart(i):&
                 !   device_start+layerend(i),&
                 !   device_start+layerstart(i+1):&
                 !   device_start+layerend(i+1)) ! V_n_n+1
        
        Gn(i,:,:) = GnL(i,:,:) + matmul(matmul(matmul(matmul(GLeft(i,:,:),&
                    Vcoup),Gn(i+1,:,:)),transpose(dconjg(Vcoup))),&
                    transpose(dconjg(GLeft(i,:,:)))) &
                    + matmul(matmul(GnL(i,:,:),Vcoup),&
                    transpose(dconjg(G_q_qp(i,:,:))))&
                    + matmul(matmul(G_q_qp(i,:,:),transpose(dconjg(Vcoup))),&
                    GnL(i,:,:)) 
        do j = 1,nlayer
            dGndT(i,j,:,:) = dGnLdT(i,j,:,:) + matmul(matmul(matmul(matmul(GLeft(i,:,:),&
                   Vcoup),dGndT(i+1,j,:,:)),transpose(dconjg(Vcoup))),&
                   transpose(dconjg(GLeft(i,:,:)))) &
                 + matmul(matmul(dGnLdT(i,j,:,:),Vcoup),&
                   transpose(dconjg(G_q_qp(i,:,:))))&
                 + matmul(matmul(G_q_qp(i,:,:),transpose(dconjg(Vcoup))),&
                   dGnLdT(i,j,:,:))
        end do
    end do

    ! left flux
    flux_left = (flux_left - trace(matmul(gamma_first,Gn(1,:,:))))*freq

   
    ! zero frequency limit
    if (ie.eq.1) then
        Gn(:,:,:) = 0.0d0
        dGndT(:,:,:,:) = 0.0d0
        A_i(:,:,:) = 0.0d0
    end if

    ! flux from left and right contact

    
    ! flux through each probe 
    do i = 1,nlayer        
        flux_i(i) = -1.0d0*trace(matmul(gamma_i(i,:,:),Gn(i,:,:)))&
                    + trace(matmul(gamma_i(i,:,:),A_i(i,:,:)))&
                    * bose_einstein(freq,temperature_i(i))
        do j = 1,nlayer
            jacobian(nn,ie,i,j) = -1.0d0*trace(matmul(gamma_i(i,:,:),dGndT(i,j,:,:)))
            if (i.eq.j) then
                jacobian(nn,ie,i,j) = jacobian(nn,ie,i,j) + &
                                  trace(matmul(gamma_i(i,:,:),A_i(i,:,:))&
                                  *dbosedT(freq,temperature_i(j)))+ &
                                  trace(matmul(dgamma_idt(i,:,:),A_i(i,:,:))&
                                  *bose_einstein(freq,temperature_i(j)))
            end if
        end do
    end do

    do i = 1,nlayer
        A_i_sum(i) =  trace(A_i(i,:,:))
        Gn_i(i) =  trace(Gn(i,:,:))
    end do


    deallocate(sigma_first,sigma_last,sigma_prev,Vcoup,Vcoup1,sigma_in)

end subroutine


subroutine integrate_jacobian(jacobian_ij,j_all,flux_i,flux)
    
    use config 
    use device  , only : nlayer
    use util
    
    implicit none

    complex(kind=8) :: jacobian_ij(:,:),j_all(:,:,:,:)
    complex(kind=8) :: flux_i(:,:,:),flux(:)
    integer(kind=4) :: i,j 
    complex(kind=8),allocatable :: integrand(:)

    allocate(integrand(ne))
    do i = 1,nlayer
        integrand(:) = sum(flux_i(:,:,i),1)*sqrt(Egrid(:,1))/dble(nk)
        flux(i) = simpson(sqrt(Egrid(:,1)),integrand)
        do j = 1,nlayer
            integrand(:) = sum(j_all(:,:,i,j),1)*sqrt(Egrid(:,1))/dble(nk)
            jacobian_ij(i,j) = simpson(sqrt(Egrid(:,1)),integrand)
        end do
    end do

end subroutine

subroutine update_temperature(jacobian_ij,flux,tdiff)

    use util
    use device  , only : nlayer

    implicit none

    complex(kind=8) :: jacobian_ij(:,:)
    complex(kind=8) :: flux(:)
    real(kind=8)    :: tdiff
    real(kind=8),allocatable    :: temp_old(:)


    allocate(temp_old(nlayer))
    temp_old = temperature_i

    temperature_i(:) = temperature_i(:) - real(matmul(inv(jacobian_ij),flux))
    if (minval(temperature_i) .lt. 0.0d0 ) then
        temperature_i = 1.0d0
    end if

   
    tdiff = sqrt(dot_product(temperature_i-temp_old,temperature_i-temp_old)/dble(nlayer))
    write(*,*) '|| T_old - T_new || = ', tdiff


end subroutine

subroutine write_flux(flux_all,gn_all,a_all,flux_left,flux_right)

    implicit none

    complex(kind=8) :: flux_all(:,:,:) ! flux through probes
    complex(kind=8) :: gn_all(:,:,:)
    complex(kind=8) :: a_all(:,:,:)
    complex(kind=8) :: flux_left(:,:)
    complex(kind=8) :: flux_right(:,:)


    open(unit=4,file="flux.dat",status="UNKNOWN",action="write",form="unformatted")
    open(unit=5,file="gn.dat",status="UNKNOWN",action="write",form="unformatted")
    open(unit=6,file="A.dat",status="UNKNOWN",action="write",form="unformatted")
    open(unit=7,file="flux_left.dat",status="UNKNOWN",action="write",form="unformatted")
    open(unit=8,file="flux_left.dat",status="UNKNOWN",action="write",form="unformatted")
    write(4) real(flux_all)
    write(5) real(gn_all)
    write(6) real(a_all)
    write(7) real(flux_left)
    write(8) real(flux_right)
    close(4)
    close(5)
    close(6)
    close(7)
    close(8)

end subroutine

subroutine get_ph_temperature(A_i_sum,Gn_i,temperature_ph)

    use device  , only : nlayer
    use config
    use util

    implicit none

    complex(kind=8),intent(in) :: gn_i(:,:)
    complex(kind=8),intent(in) :: a_i_sum(:,:)
    real(kind=8)    :: temperature_ph(:)
    real(kind=8),allocatable  :: temperature_old(:)
    real(kind=8),allocatable  :: dflayer(:),flayer(:)
    real(kind=8)    :: diff,jdiff
    integer(kind=4)  :: i,nn

    allocate(dflayer(ne),flayer(ne),temperature_old(nlayer))

    temperature_old(:) = -100.0d0
    
    do nn = 1,nlayer
        do while (abs(temperature_old(nn)-temperature_ph(nn)).gt.1.0e-4)
            temperature_old(nn) = temperature_ph(nn)
            do i = 1,ne
                dflayer(i) = dbosedT(sqrt(Egrid(i,1)),temperature_ph(nn))
                flayer(i) = bose_einstein(sqrt(Egrid(i,1)),temperature_ph(nn))
            end do
            jdiff = real(simpson(sqrt(Egrid(:,1)),&
                 -dflayer(:)*Egrid(:,1)*a_i_sum(:,nn)))
            diff = real(simpson(sqrt(Egrid(:,1)),Gn_i(:,nn)*Egrid(:,1)&
                   -flayer(:)*Egrid(:,1)*a_i_sum(:,nn)))
            temperature_ph(nn) = temperature_old(nn) - diff/jdiff
        end do   
    end do

end subroutine

subroutine write_temperature(temperature_i,temperature_ph)

    implicit none

    real(kind=8) :: temperature_i(:),temperature_ph(:)

    open(unit=4,file="temperature_i.dat",status="UNKNOWN",action="write",form="unformatted")
    write(4) temperature_i
    close(4)

    open(unit=4,file="temperature_ph.dat",status="UNKNOWN",action="write",form="unformatted")
    write(4) temperature_ph
    close(4)
 
    
end subroutine



end module
