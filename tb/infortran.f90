program infortran

    use config
    use input
    use param
    use hamiltonian
    use surface
    use device
    use velocity_mat 
    use write_data
    use circular
    use util

    implicit none

    include "mpif.h"

    integer(kind=4) :: ll, ie, i, ik, j, k, l, nb, m
    real(kind=8) :: kcry(3),b_sc(3,3),ktp(3)
    complex(kind=8),allocatable :: Ham(:,:,:), eigvec(:,:,:)
    complex(kind=8),allocatable :: Px_L(:,:,:), Py_L(:,:,:)
    complex(kind=8),allocatable :: Px_R(:,:,:), Py_R(:,:,:)
    complex(kind=8),allocatable :: gl_adv(:,:,:), gr_ret(:,:,:)
    complex(kind=8),allocatable :: gl_plus(:,:,:), gr_minus(:,:,:)
    complex(kind=8),allocatable :: gml_adv(:,:,:), gmr_ret(:,:,:)
    complex(kind=8),allocatable :: sl_adv(:,:,:), sr_ret(:,:,:)
    complex(kind=8),allocatable :: G_D_RL(:,:,:),G_D_LR(:,:,:)
    complex(kind=8),allocatable :: G_L(:,:,:),G_R(:,:,:)
    complex(kind=8),allocatable :: tlr(:,:,:,:),trl(:,:,:,:)
    real(kind=8),allocatable :: transl(:,:),transr(:,:)
    real(kind=8),allocatable :: transl_s(:,:),transr_s(:,:)
    real(kind=8),allocatable :: transl_ns(:,:),transr_ns(:,:)
    real(kind=8),allocatable :: transl_reduce(:,:),transr_reduce(:,:)
    real(kind=8),allocatable :: transl_s_reduce(:,:),transr_s_reduce(:,:)
    real(kind=8),allocatable :: transl_ns_reduce(:,:),transr_ns_reduce(:,:)
    real(kind=8),allocatable :: refll_s(:,:),reflr_s(:,:)
    real(kind=8),allocatable :: refll_ns(:,:),reflr_ns(:,:)
    real(kind=8),allocatable :: refll_ns_reduce(:,:),reflr_ns_reduce(:,:)
    real(kind=8),allocatable :: refll_s_reduce(:,:),reflr_s_reduce(:,:)
    real(kind=8),allocatable :: surface_dos_l(:,:),surface_dos_r(:,:)
    real(kind=8),allocatable :: surface_dos_l_reduce(:,:),surface_dos_r_reduce(:,:)
    real(kind=8),allocatable :: device_dos_reduce(:,:),device_dos(:,:)
    real(kind=8),allocatable :: ldos(:,:,:,:)
    real(kind=8),allocatable :: eig(:,:)
    complex(kind=8),allocatable :: eigv_pc(:,:)
    complex(kind=8),allocatable :: Ham_f(:,:,:)
    complex(kind=8),allocatable :: Ham_uc_l(:,:),Ham_uc_r(:,:)
    real(kind=8),allocatable    :: eig_uc_l(:),eig_uc_r(:)
    complex(kind=8),allocatable :: evs_uc_l(:,:),evs_uc_r(:,:)
    complex(kind=8),allocatable :: Ham_pc_l(:,:),Ham_pc_r(:,:)
    real(kind=8),allocatable    :: eig_pc_l(:),eig_pc_r(:)
    complex(kind=8),allocatable :: evs_pc_l(:,:),evs_pc_r(:,:)
    ! slab calculation for surface band structure
    complex(kind=8),allocatable :: Ham_slab_l(:,:),Ham_slab_r(:,:)
    real(kind=8),allocatable    :: eig_slab_l(:),eig_slab_r(:)
    complex(kind=8),allocatable :: evs_slab_l(:,:),evs_slab_r(:,:)
    real(kind=8),allocatable ::  pos_slab_l(:,:),pos_slab_r(:,:)
    real(kind=8) ::  latvec_slab(3,3),recvec_slab(3,3)
    integer(kind=4) :: natoms_slab_l,natoms_slab_r
    integer(kind=4) :: norb_slab_l,norb_slab_r
    integer(kind=4),allocatable :: atom_orb_num_slab_l(:)
    integer(kind=4),allocatable :: atom_orb_num_slab_r(:)
    integer(kind=4)             :: norb ! total number of orbitals
    integer(kind=4) :: norb_pc_l, norb_uc_l, norb_pc_r, norb_uc_r          
    integer(kind=4) :: temp_count,ntemp
    real(kind=8) :: p_crys_temp(2)
    ! search k point for misorientation unit cell
    real(kind=8),allocatable :: k_misori(:,:)
    integer(kind=4) :: n_k_misori
    real(kind=8) :: k_mis_temp(3),dk_mis_temp(3)
    real(kind=8),allocatable :: kuc_misori(:,:)

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

    call load_configure()
    call load_input()
    call k_grid_in_ws()
    call load_param()
   
    ! obtain the size of Hamiltonian 
    norb = sum(atom_orb_num) ! norb can be known only after loading params
    norb_pc_l = sum(atom_orb_num_pc_l)
    norb_pc_r = sum(atom_orb_num_pc_r)
    norb_uc_l = sum(atom_orb_num_uc_l)
    norb_uc_r = sum(atom_orb_num_uc_r)

    allocate(pos_orb_pc_l(norb_pc_l,4))
    allocate(pos_orb_uc_l(norb_uc_l,4))
    allocate(pos_orb_pc_r(norb_pc_r,4))
    allocate(pos_orb_uc_r(norb_uc_r,4))

    ! get the dimensions for Hamiltonian
    ! of primitive cell
    call get_pos_orb(norb_pc_l,norb_uc_l,&
                     norb_pc_r,norb_uc_r)
    
    ! calculate the eigenvales on k points
    ! along k path
    if (path_mode.eq.1)then
        if (side(1:1) .eq. 'l') then
            if (allocated(Egrid)) deallocate(Egrid)
            allocate(Egrid(norb_pc_l,n_kpath_points),&
                 eigv_pc(norb_pc_l,norb_pc_l))
            ne = norb_pc_l ! here, we overwrite ne
            nk = n_kpath_points ! and nk
            nkcore = nk/numprocs
            Egrid = 0.0d0
            eigv_pc = 0.0d0
            nb = max(n_bloch_x,n_bloch_y)
            deallocate(k_grid)
            allocate(k_grid(nk,3))
            k_grid = 100.0
            b_sc = recivec
            b_sc(3,3) = 2*pi/az1
            do i = 1,n_kpath_points
                !call gen_ham_pc_l(norb_pc_l,&
                !matmul(kpath_points(i,:),recivec_pc_l),&
                !Egrid(:,i),eigv_pc)
                do j = -nb,nb
                    do k = -nb,nb
                        do l = -nb,nb
                            ktp = matmul(kpath_points(i,:),&
                            recivec_pc_l)+&
                            matmul((/dble(j),dble(k),dble(l)/),&
                            b_sc)
                            kcry = matmul(ktp,inv_real(b_sc))
                            !write(*,*) kpath_points(i,:)
                            if (abs(kcry(1).lt.0.5).and.&
                                abs(kcry(2).lt.0.5).and.&
                                abs(kcry(3).lt.0.5)) then
                                k_grid(i,:) = ktp    
                            end if
                        end do
                    end do
                end do
               ! write(*,*)'uc', Egrid(:,1)
            end do
        end if
    end if


    allocate(Ham(numprocs,norb,norb))
    allocate(eigvec(numprocs,norb,norb))
    allocate(eig(numprocs,norb))

    Ham = 0.0d0
    eigvec = 0.0d0
    eig = 0.0d0
       
    ! i th k point on a core  
    call gen_ham_mpi(1,1,Ham,eigvec,eig,norb)   

    call get_periodic_lead_orb_num()

    if (path_mode.eq.1)then
        ! calculate the uc band
!        allocate(Ham_uc_l(norb_uc_l,norb_uc_l),&
!                 eig_uc_l(norb_uc_l),&
!                 evs_uc_l(norb_uc_l,norb_uc_l))
!        allocate(Ham_uc_r(norb_uc_r,norb_uc_r),&
!                 eig_uc_r(norb_uc_r),&
!                 evs_uc_r(norb_uc_r,norb_uc_r))

!        do i = 1,n_kpath_points
!            call gen_ham_uc_l(norb_uc_l,&
!                           matmul(kpath_points(i,:),&
!                           reci_uc_l),Ham_uc_l,&
!                           eig_uc_l,evs_uc_l)
!            call write_eig('eigenl',eig_uc_l)
!            call gen_ham_uc_r(norb_uc_r,&
!                           matmul(kpath_points(i,:),&
!                           reci_uc_r),Ham_uc_r,&
!                           eig_uc_r,evs_uc_r)
!            call write_eig('eigenr',eig_uc_r)
!        end do
        ! the unit cell for misorientation
        uc_misori = 2*pi*transpose(inv_real(uc_misori))
        allocate(k_misori((2*n_uc_misori1+1)*(2*n_uc_misori2+1)*(2*n_uc_misori2+1),3))
        allocate(kuc_misori(size(kpath_points,1),3))
        kuc_misori = 0.0d0
        ! Calculate the pc band
        allocate(Ham_pc_l(norb_pc_l,norb_pc_l),&
                 eig_pc_l(norb_pc_l),&
                 evs_pc_l(norb_pc_l,norb_pc_l))
        allocate(Ham_pc_r(norb_pc_r,norb_pc_r),&
                 eig_pc_r(norb_pc_r),&
                 evs_pc_r(norb_pc_r,norb_pc_r))

        do i = 1,n_kpath_points
            k_misori = 0.0d0 
            n_k_misori = 0
            k_mis_temp(:) = matmul(kpath_points(i,:),&
                                recivec_pc_l)
            do j = -n_uc_misori1,n_uc_misori1
                do k = -n_uc_misori2,n_uc_misori2
                    do l = -n_uc_misori3,n_uc_misori3
                        dk_mis_temp(:) = matmul((/dble(j),dble(k),dble(l)/),&
                                         uc_misori)
                        if (in_ws(uc_misori,k_mis_temp+dk_mis_temp)>0) then
                            n_k_misori = n_k_misori + 1
                            k_misori(n_k_misori,:) = k_mis_temp+dk_mis_temp
                        end if
                    end do
                end do
            end do
            kuc_misori(i,:) = find_min_dist(k_misori(1:n_k_misori,:),n_k_misori)
            call gen_ham_pc_l(norb_pc_l,&
                           matmul(kpath_points(i,:),&
                           recivec_pc_l),Ham_pc_l)
            call eigenH(Ham_pc_l,eig_pc_l,evs_pc_l)
            call write_eig('eigenl',eig_pc_l)
            call gen_ham_pc_r(norb_pc_r,&
                           matmul(kpath_points(i,:),&
                           recivec_pc_r),Ham_pc_r)
            call eigenH(Ham_pc_r,eig_pc_r,evs_pc_r)
            call write_eig('eigenr',eig_pc_r)
        end do
        open(unit=15, file=trim(adjustl(output_dir))//"/kuc_misori.dat",status="REPLACE",action="write")
        m = size(kuc_misori,1)
        write(15,*) m
        do i = 1,m
             WRITE(15,1015) kuc_misori(i,:)
             1015 format(3f20.10)
        end do
        close(unit=15)
        stop
    end if

    if (path_mode .eq. 4)  then
        latvec_slab = 0.0d0
        recvec_slab = 0.0d0
        latvec_slab(1:2,1:2) = latvec(1:2,1:2) 
        latvec_slab(3,3) = latvec(3,3)*5 ! large vaccuum
        recvec_slab = 2.0d0* pi * transpose(inv_real(latvec_slab))
        natoms_slab_l = 2*period_left+buffer_left
        natoms_slab_r = 2*period_right+buffer_right
        norb_slab_l = norb_uc_l/period_left*natoms_slab_l
        norb_slab_r = norb_uc_r/period_right*natoms_slab_r
        allocate(Ham_slab_l(norb_slab_l,norb_slab_l),&
                 eig_slab_l(norb_slab_l),&
                 evs_slab_l(norb_slab_l,norb_slab_l),&
                 pos_slab_l(natoms_slab_l,4),&
                 atom_orb_num_slab_l(natoms_slab_l))
        allocate(Ham_slab_r(norb_slab_r,norb_slab_r),&
                 eig_slab_r(norb_slab_r),&
                 evs_slab_r(norb_slab_r,norb_slab_r),&
                 pos_slab_r(natoms_slab_r,4),&
                 atom_orb_num_slab_r(natoms_slab_r))
        pos_slab_l = 0.0d0
        pos_slab_l(1:natoms_slab_l,:) = positions(1:natoms_slab_l,:) 
        atom_orb_num_slab_l = atom_orb_num(1:natoms_slab_l)
        pos_slab_r = 0.0d0
        pos_slab_r(1:natoms_slab_r,:) = positions(natoms-natoms_slab_r+1:natoms,:) 
        atom_orb_num_slab_r = atom_orb_num(natoms-natoms_slab_r+1:natoms)
        do i = 1,n_kpath_points
            call gen_ham_slab(norb_slab_l,&
                          matmul(kpath_points(i,:),recvec_slab),&
                          Ham_slab_l,&
                          eig_slab_l,evs_slab_l,pos_slab_l,&
                          latvec_slab,atom_orb_num_slab_l)
            call gen_ham_slab(norb_slab_r,&
                          matmul(kpath_points(i,:),recvec_slab),&
                          Ham_slab_r,&
                          eig_slab_r,evs_slab_r,pos_slab_r,&
                          latvec_slab,atom_orb_num_slab_r)
            call write_eig('slabl',eig_slab_l)
            call write_eig('slabr',eig_slab_r)
        end do
        stop
    end if

    if (path_mode.eq.5) then
        ! calculate the uc band
        allocate(Ham_uc_l(norb_uc_l,norb_uc_l),&
                 eig_uc_l(norb_uc_l),&
                 evs_uc_l(norb_uc_l,norb_uc_l))
        allocate(Ham_uc_r(norb_uc_r,norb_uc_r),&
                 eig_uc_r(norb_uc_r),&
                 evs_uc_r(norb_uc_r,norb_uc_r))

        ntemp = 100

        do i = 1,n_kpath_points
            do j = -ntemp, ntemp
                call gen_ham_uc_l(norb_uc_l,&
                           matmul(kpath_points(i,:)&
                           +(/0.0d0,0.0d0,dble(j)/dble(ntemp)/),&
                           reci_uc_l),Ham_uc_l,&
                           eig_uc_l,evs_uc_l)
                call write_eig('eigen_uc_l',eig_uc_l)
                call gen_ham_uc_r(norb_uc_r,&
                           matmul(kpath_points(i,:)&
                           +(/0.0d0,0.0d0,dble(j)/dble(ntemp)/),&
                           reci_uc_r),Ham_uc_r,&
                           eig_uc_r,evs_uc_r)
                call write_eig('eigen_uc_r',eig_uc_r)
            end do
        end do
        stop
    end if


    allocate(gl_adv(numprocs,nl,nl))
    allocate(gr_ret(numprocs,nr,nr))
    allocate(gl_plus(numprocs,nl,nl))
    allocate(gr_minus(numprocs,nl,nl))
    allocate(gml_adv(numprocs,nl,nl))
    allocate(gmr_ret(numprocs,nr,nr))
    allocate(sl_adv(numprocs,nl,nl))
    allocate(sr_ret(numprocs,nr,nr))
    allocate(transl(ne,nk))
    allocate(transr(ne,nk))
    allocate(transl_s(ne,nk))
    allocate(transr_s(ne,nk))
    allocate(transl_ns(ne,nk))
    allocate(transr_ns(ne,nk))
    allocate(transl_reduce(ne,nk))
    allocate(transr_reduce(ne,nk))
    allocate(transl_s_reduce(ne,nk))
    allocate(transr_s_reduce(ne,nk))
    allocate(transl_ns_reduce(ne,nk))
    allocate(transr_ns_reduce(ne,nk))
    allocate(refll_s(ne,nk))
    allocate(reflr_s(ne,nk)) 
    allocate(refll_ns(ne,nk))
    allocate(reflr_ns(ne,nk))
    allocate(refll_s_reduce(ne,nk))
    allocate(reflr_s_reduce(ne,nk)) 
    allocate(refll_ns_reduce(ne,nk))
    allocate(reflr_ns_reduce(ne,nk))
    allocate(surface_dos_l(ne,nk))
    allocate(surface_dos_r(ne,nk))
    allocate(surface_dos_l_reduce(ne,nk))
    allocate(surface_dos_r_reduce(ne,nk))
    allocate(device_dos(ne,nk))
    allocate(device_dos_reduce(ne,nk))
    allocate(G_D_RL(numprocs,nr,nl))
    allocate(G_D_LR(numprocs,nl,nr))
    allocate(G_L(numprocs,nl,nl))
    allocate(G_R(numprocs,nr,nr))
    allocate(P1(nl,nl),Pinv1(nl,nl))
    allocate(P2(nl/n_bloch_y,nl/n_bloch_y),&
             Pinv2(nl/n_bloch_y,nl/n_bloch_y))

    gl_adv = 0.0d0
    gr_ret = 0.0d0
    gl_plus = 0.0d0
    gr_minus = 0.0d0
    gml_adv = 0.0d0
    gmr_ret = 0.0d0
    sl_adv = 0.0d0
    sr_ret = 0.0d0
    G_D_RL = 0.0d0
    G_D_LR = 0.0d0
    G_L = 0.0d0
    G_R = 0.0d0
    P1 = 0.0d0
    Pinv1 = 0.0d0
    P2 = 0.0d0
    Pinv2 = 0.0d0
    ! set up the permutation matrix for circular matrix inversion
    call set_pmatrix(nl,n_bloch_y,&
         reci_uc_l(2,:),latvec_uc_l(2,:),&
         P1,Pinv1)
    call set_pmatrix(nl/n_bloch_y,n_bloch_x,&
         reci_uc_l(1,:),latvec_uc_l(1,:),&
         P2,Pinv2)

    transl = 0.0d0
    transr = 0.0d0
    transl_s = 0.0d0
    transr_s = 0.0d0
    transl_ns = 0.0d0
    transr_ns = 0.0d0
    transl_reduce = 0.0d0
    transr_reduce = 0.0d0
    transl_s_reduce = 0.0d0
    transr_s_reduce = 0.0d0
    transl_ns_reduce = 0.0d0
    transr_ns_reduce = 0.0d0
    refll_s = 0.0d0
    reflr_s = 0.0d0
    refll_ns = 0.0d0
    reflr_ns = 0.0d0
    refll_s_reduce = 0.0d0
    reflr_s_reduce = 0.0d0
    refll_ns_reduce = 0.0d0
    reflr_ns_reduce = 0.0d0
    surface_dos_l = 0.0d0
    surface_dos_r = 0.0d0
    surface_dos_l_reduce = 0.0d0
    surface_dos_r_reduce = 0.0d0
    device_dos = 0.0d0
    device_dos_reduce = 0.0d0

    if (myid .eq. 0) then
        call blockwise_setup(Ham(1,:,:),norb)
        call transverse_uc_setup()
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
 
    ! broadcast non-array variable   
    call MPI_BCAST(nlayer,1,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(nd,1,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(maxlen,1,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(natm_l_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(natm_r_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nl_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(nr_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr) 

    ! still we need to run allocation on other cores.
    if (myid .gt. 0) then
        allocate(layer_list(nlayer),layerstart(nlayer),&
        layerend(nlayer),pidxL(nl),pidxR(nr),&
        L_orb_num_list(natm_l_uc),R_orb_num_list(natm_r_uc),&
        L_list(natm_l_uc),R_list(natm_r_uc),&
        pos_l(nl,3),pos_r(nr,3))
    end if

    ! finally we broadcast array variables
    call MPI_BCAST(layer_list,nlayer,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(layerstart,nlayer,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(layerend,nlayer,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(pidxL, nl, MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr) 
    call MPI_BCAST(pidxR, nr, MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(L_orb_num_list,natm_l_uc,MPI_INTEGER,0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(R_orb_num_list,natm_r_uc,MPI_INTEGER,0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(L_list,natm_l_uc,MPI_INTEGER,0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(R_list,natm_r_uc,MPI_INTEGER,0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(pos_l,nl*3,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(pos_r,nr*3,MPI_DOUBLE_PRECISION,0, &
         MPI_COMM_WORLD, ierr)

    allocate(orb_uc_l_idx(norb_uc_l),&
             orb_uc_r_idx(norb_uc_r))
   
    temp_count = 0
    do i = 1,nl
        p_crys_temp = matmul(pos_l(i,1:2),&
                      inv_real(latvec_uc(1:2,1:2)))
        if ((p_crys_temp(1).gt.-0.001).and.&
            (p_crys_temp(1).lt. 0.999).and.&
            (p_crys_temp(2).gt.-0.001).and.&
            (p_crys_temp(2).lt. 0.999)) then
            temp_count = temp_count + 1
            orb_uc_l_idx(temp_count) = i
        end if
    end do
    temp_count = 0
    do i = 1,nr
        p_crys_temp = matmul(pos_r(i,1:2),&
                      inv_real(latvec_uc(1:2,1:2)))
        if ((p_crys_temp(1).gt.-0.001).and.&
            (p_crys_temp(1).lt. 0.999).and.&
            (p_crys_temp(2).gt.-0.001).and.&
            (p_crys_temp(2).lt. 0.999)) then
            temp_count = temp_count + 1
            orb_uc_r_idx(temp_count) = i
        end if
    end do

    ! below is the start of the (heavy) calculation

    allocate(ldos(ne,numprocs,nlayer,maxlen))
    ldos = 0.0d0

    if (path_mode.ne.6) then
        do ik = 1,nkcore
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call gen_ham_mpi(1,ik,Ham,eigvec,eig,norb)        
            do ie = 1,ne
        
                call surface_green_mpi(ie,ik,Ham,norb,&
                                   gl_adv,gr_ret,gl_plus,gr_minus,&
                                   gml_adv,gmr_ret,&
                                   sl_adv,sr_ret,Egrid)
                call surface_ldos_mpi(ie,ik,&
                                   gl_adv,gr_ret,&
                                   surface_dos_l_reduce,surface_dos_r_reduce)

                   
                call G_block_mpi(ie,ik,Ham,norb,&
                           gl_adv,gr_ret,&
                           gml_adv,gmr_ret,sl_adv,sr_ret,Egrid,&
                           transl_reduce(:,:),transr_reduce(:,:),&
                           G_D_RL(:,:,:),G_D_LR(:,:,:),&
                           G_L(:,:,:),G_R(:,:,:))
               
                call get_ldos_mpi(ie,ik,Ham,norb,gl_adv,gr_ret,&
                     gml_adv,gmr_ret,&
                     sl_adv,sr_ret,Egrid,ldos(:,:,:,:),device_dos_reduce(:,:))       

                call vel_mat_mpi(ie,ik,Egrid,Ham,norb,&
                     gl_adv,gr_ret,gl_plus,gr_minus,&
                     gml_adv,gmr_ret,&
                     G_D_RL(:,:,:),G_D_LR(:,:,:),&
                     G_L(:,:,:),G_R(:,:,:),&
                     transl_s_reduce(:,:),transr_s_reduce(:,:),&
                     transl_ns_reduce(:,:),transr_ns_reduce(:,:),&
                     refll_s_reduce(:,:),reflr_s_reduce(:,:),&
                     refll_ns_reduce(:,:),reflr_ns_reduce(:,:),&
                     norb_pc_l,norb_pc_r,&
                     norb_uc_l,norb_uc_r)

                if (myid.eq.0) then
                    write(*,23,ADVANCE='NO') dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0
                    do i = 1,50
                        if ((i-1)*100/50 .le.&
                            int(dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0)) then
                            write(*,"(A)",ADVANCE='NO')'='
                        else 
                            write(*,"(A)",ADVANCE='NO')' '
                        end if
                    end do
                    write(*,"(A)")']'
                    23  format ('Process: ',f6.2,' % [')  
                end if
            end do
            call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
        call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(device_dos_reduce,device_dos,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)


         if (myid.eq.0) then
             call write_device_dos(device_dos) 
             call write_surface_ldos(surface_dos_l,surface_dos_r) 
         end if
    else
        do ik = 1,nkcore
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call gen_ham_mpi(1,ik,Ham,eigvec,eig,norb)        
            do ie = 1,ne
        
                call surface_green_mpi(ie,ik,Ham,norb,&
                                   gl_adv,gr_ret,gl_plus,gr_minus,&
                                   gml_adv,gmr_ret,&
                                   sl_adv,sr_ret,Egrid)                   
                call surface_ldos_mpi(ie,ik,&
                                   gl_adv,gr_ret,&
                                   surface_dos_l_reduce,surface_dos_r_reduce)
                if (myid.eq.0) then
                    write(*,24,ADVANCE='NO') dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0
                    do i = 1,50
                        if ((i-1)*100/50 .le.&
                            int(dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0)) then
                            write(*,"(A)",ADVANCE='NO')'='
                        else 
                            write(*,"(A)",ADVANCE='NO')' '
                        end if
                    end do
                    write(*,"(A)")']'
                    24  format ('Process: ',f6.2,' % [')  
                end if
            end do
            call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

        call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    
         if (myid.eq.0) then
             call write_surface_ldos(surface_dos_l,surface_dos_r) 
         end if
         call MPI_FINALIZE(ierr)
         stop
    end if

    call MPI_ALLREDUCE(transl_reduce,transl,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(transr_reduce,transr,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(transl_s_reduce,transl_s,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(transr_s_reduce,transr_s,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(transl_ns_reduce,transl_ns,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(transr_ns_reduce,transr_ns,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(refll_s_reduce,refll_s,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(reflr_s_reduce,reflr_s,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(refll_ns_reduce,refll_ns,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(reflr_ns_reduce,reflr_ns,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)



    call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

    if (myid.eq.0) then
        if (path_mode.ne.1)then
            call write_trans(Egrid,&
            transl_s,transl_ns,&
            transr_s,transr_ns,&
            transl,transr,&
            refll_s,refll_ns,&
            reflr_s,reflr_s)
        end if
    end if

    call MPI_FINALIZE(ierr)

end program infortran
