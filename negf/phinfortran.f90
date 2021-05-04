program phinfortran
     
    use config 
    use input
    use fc
    use dyn
    use circular
    use surface
    use write_data
    use device
    use lesser
    use velocity_mat
    use util

    implicit none

    include "mpif.h"

    integer(kind=4) :: ll, ie, i, ik, j, k, l, nb, iit, niit
    real(kind=8),allocatable :: force_constant(:,:)
    real(kind=8),allocatable :: force_constant_l(:,:)
    real(kind=8),allocatable :: force_constant_r(:,:)
    integer(kind=4),allocatable :: idx_mat(:,:,:),idx_mat_l(:,:,:),idx_mat_r(:,:,:)
    integer(kind=4),allocatable :: idx_mat_l_v(:,:,:),idx_mat_r_v(:,:,:)
    integer(kind=4),allocatable :: idx_mat_red(:),idx_mat_l_red(:),idx_mat_r_red(:)
    integer(kind=4),allocatable :: idx_mat_l_v_red(:),idx_mat_r_v_red(:)
    integer(kind=4),allocatable :: num_mat(:,:),num_mat_l(:,:),num_mat_r(:,:)
    integer(kind=4),allocatable :: num_mat_l_v(:,:),num_mat_r_v(:,:)
    real(kind=8),allocatable :: dist_mat(:,:,:,:),dist_mat_l(:,:,:,:),dist_mat_r(:,:,:,:)
    real(kind=8),allocatable :: dist_mat_l_v(:,:,:,:),dist_mat_r_v(:,:,:,:)
    real(kind=8),allocatable :: dist_mat_red(:,:),dist_mat_l_red(:,:),dist_mat_r_red(:,:)
    real(kind=8),allocatable :: dist_mat_l_v_red(:,:),dist_mat_r_v_red(:,:)

    complex(kind=8),allocatable :: Ham(:,:,:),eigvec(:,:,:)
    complex(kind=8),allocatable :: Ham_l(:,:,:),Ham_r(:,:,:)
    complex(kind=8),allocatable :: Hamd(:,:,:,:)
    complex(kind=8),allocatable :: Hamdv(:,:,:,:)
    complex(kind=8),allocatable :: Haml(:,:,:),Hamr(:,:,:)
    complex(kind=8),allocatable :: Hamlv(:,:,:),Hamrv(:,:,:)
    complex(kind=8),allocatable :: Hamld(:,:,:),Hamrd(:,:,:)
    complex(kind=8),allocatable :: Ham_l_v(:,:,:),Ham_r_v(:,:,:)
    complex(kind=8),allocatable :: eigvec_l(:,:,:),eigvec_r(:,:,:)
    real(kind=8),allocatable :: eig(:,:)
    real(kind=8),allocatable :: eig_l(:,:),eig_r(:,:)
    real(kind=8),allocatable    :: pos_left_lead(:,:),pos_right_lead(:,:)
    real(kind=8),allocatable    :: pos_left_lead_uc(:,:),pos_right_lead_uc(:,:)
    real(kind=8),allocatable    :: pos_left_lead_sc(:,:),pos_right_lead_sc(:,:)
    complex(kind=8),allocatable :: gl_adv(:,:,:), gr_ret(:,:,:)
    complex(kind=8),allocatable :: gl_plus(:,:,:), gr_minus(:,:,:)
    complex(kind=8),allocatable :: gml_adv(:,:,:), gmr_ret(:,:,:)
    complex(kind=8),allocatable :: sl_adv(:,:,:), sr_ret(:,:,:)
    complex(kind=8),allocatable :: G_D_RL(:,:,:),G_D_LR(:,:,:)
    complex(kind=8),allocatable :: GnL(:,:,:,:),Gn(:,:,:,:)
    complex(kind=8),allocatable :: GL(:,:,:,:),G_full(:,:,:,:)
    complex(kind=8),allocatable :: dGnLdT(:,:,:,:,:),dGndT(:,:,:,:,:)
    complex(kind=8),allocatable :: G_q_qp(:,:,:,:),G_qp_q(:,:,:,:)
    complex(kind=8),allocatable :: G_N_Np(:,:,:),G_Np_N(:,:,:)
    complex(kind=8),allocatable :: G_z_zp(:,:,:),G_zp_z(:,:,:)
    complex(kind=8),allocatable :: jacobian_all(:,:,:,:),jacobian_ij(:,:)
    complex(kind=8),allocatable :: flux_all(:,:,:),flux_i(:,:,:),flux(:)

    complex(kind=8),allocatable :: flux_left(:,:),flux_left_reduce(:,:)
    complex(kind=8),allocatable :: flux_right(:,:),flux_right_reduce(:,:)
    complex(kind=8),allocatable :: A_all(:,:,:),A_i(:,:,:)
    complex(kind=8),allocatable :: Gn_all(:,:,:),Gn_i(:,:,:)
    real(kind=8)                :: tdiff
    real(kind=8),allocatable    :: temperature_ph(:)
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
    integer(kind=4) :: norb, norb_sc, norb_uc, norb_layer
    integer(kind=4) ::  norb_uc_l,  norb_uc_r
    integer(kind=4) :: temp_sum, temp_count
    complex(kind=8),allocatable :: dyn_uc_l(:,:),dyn_uc_r(:,:) 
    real(kind=8),allocatable    :: eig_uc_l(:),eig_uc_r(:)
    complex(kind=8),allocatable :: evs_uc_l(:,:),evs_uc_r(:,:) 
    complex(kind=8),allocatable :: dyn_pc_l(:,:),dyn_pc_r(:,:) 
    real(kind=8),allocatable    :: eig_pc_l(:),eig_pc_r(:)
    complex(kind=8),allocatable :: evs_pc_l(:,:),evs_pc_r(:,:) 
    real(kind=8),allocatable    :: vels_pc_l(:,:),vels_pc_r(:,:) 
    real(kind=8),allocatable    :: vels_uc_l(:,:),vels_uc_r(:,:) 
    complex(kind=8),allocatable :: dyn_slab_l(:,:),dyn_slab_r(:,:)
    real(kind=8),allocatable    :: eig_slab_l(:),eig_slab_r(:)
    complex(kind=8),allocatable :: evs_slab_l(:,:),evs_slab_r(:,:)
    real(kind=8),allocatable ::  pos_slab_l(:,:),pos_slab_r(:,:)
    real(kind=8) ::  latvec_slab(3,3),recvec_slab(3,3)
    integer(kind=4) :: natoms_slab_l,natoms_slab_r
    integer(kind=4) :: norb_slab_l,norb_slab_r
    real(kind=8) :: kpath3(3)
    real(kind=8) :: cell_qe(3,3),cell_sc_qe(3,3)
    complex(kind=8),allocatable :: dyn_qe(:,:),dyn_all(:,:),dyn_pc_qe(:,:)
    real(kind=8),allocatable    :: fc_qe(:,:,:,:,:,:,:),pos_qe(:,:),mass_qe(:)
    real(kind=8),allocatable :: pos_sc_qe(:,:),mass_sc_qe(:)
    integer(kind=4) :: nx_qe,ny_qe,nz_qe
    integer(kind=4),allocatable :: idx_sc2uc(:),idx_all(:)
    real(kind=8)    :: rws(124,3),rd(124)
    real(kind=8)    :: cell_all(3,3)
    real(kind=8),allocatable    :: pos_all(:,:),mass_all(:)
    integer(kind=4),allocatable :: mass_id(:)
    real(kind=8) :: intprogress
    integer(kind=4) :: intstep
    integer(kind=4),allocatable :: neighborlist(:,:),nnbrs(:)


    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

    call load_configure()

    if (qe_fc .ne. 0) then
        call read_qe_fc(cell_qe,cell_sc_qe,pos_qe,pos_sc_qe,&
                       idx_sc2uc,idx_all,mass_qe,mass_sc_qe,fc_qe,&
                       nx_qe,ny_qe,nz_qe,rws,rd,cell_all,pos_all,mass_id) 
        natoms = size(pos_all,1)
        natoms_uc = size(pos_qe,1)

        allocate(positions(size(pos_all,1),size(pos_all,2)))
        positions = pos_all
         
        norb = natoms*3 
        norb_uc = natoms_uc*3 
        norb_layer = natoms_uc*3*n_bloch_x*n_bloch_y
        norb_uc_l = natoms_uc*3
        norb_uc_r = natoms_uc*3

        latvec = cell_all
        latvec_left = cell_all
        latvec_right = cell_all
        latvec_uc = cell_qe
        latvec_uc_l = cell_qe
        latvec_uc_r = cell_qe


        recivec = 2*pi*transpose(inv_real(latvec))
        reci_uc_l = 2*pi*transpose(inv_real(latvec_uc_l))
        reci_uc_r = 2*pi*transpose(inv_real(latvec_uc_r))
        left_buffer_orb =  n_bloch_x*n_bloch_y*natoms_uc*3
        right_buffer_orb = n_bloch_x*n_bloch_y*natoms_uc*3
        nl = n_bloch_x*n_bloch_y*natoms_uc*3
        nr = n_bloch_x*n_bloch_y*natoms_uc*3
        n_buffer_l = n_bloch_x*n_bloch_y*natoms_uc*3
        n_buffer_r = n_bloch_x*n_bloch_y*natoms_uc*3
        a_z_L = latvec_uc_l(3,3)
        a_z_R = latvec_uc_r(3,3)

        nd = norb-2*nl-2*nr-left_buffer_orb-right_buffer_orb
        nlayer = nz - 6

        call read_pc_qe()
        allocate(neighborlist(natoms,200))
        allocate(nnbrs(natoms))
        neighborlist(:,:)=0
        nnbrs(:)=0
        ! find neighboring list
        if (update_nbr.eq.1) then
            call find_neighbors(neighborlist,nnbrs,pos_all,cell_all,natoms) 
        end if
        
        ! read mass
        allocate(mass_all(size(pos_all,1)))
        do i = 1,size(pos_all,1)
            mass_all(i) = mass_no(mass_id(i))
        end do
        if (myid.eq.0) then
            allocate(layer_list(nlayer),layerstart(nlayer),layerend(nlayer))
            maxlen = n_bloch_x*n_bloch_y*natoms_uc*3
            layer_list(:) = maxlen 
            layerstart(1) = 1
            layerend(1) = layer_list(1)
            do i = 2,nlayer
                layerstart(i) = layerend(i-1) + 1
                layerend(i) = layerstart(i) + layer_list(i) -1
            end do

            ! initialize probe temperature at different blocks
            allocate(temperature_i(nlayer))
            if (nlayer.gt.1) then
                do i = 1,nlayer
                    temperature_i(i) = t_left + 0.25d0*(t_right-t_left) &
                                    + 0.5d0* dble(i-1)*(t_right-t_left)/dble(nlayer-1)
                end do
            else
                temperature_i(1) = t_left
            end if
            ! initialize the temperature of the lattice
            allocate(temperature_ph(nlayer))
            temperature_ph = temperature_i
                   
            allocate(dyn_qe(norb_uc,norb_uc))
            allocate(dyn_all(norb,norb))
            call gen_dyn_uc_qe(dyn_qe,cell_qe,cell_sc_qe,pos_qe,pos_sc_qe,&
                           idx_sc2uc,mass_qe,mass_sc_qe,fc_qe,&
                           nx_qe,ny_qe,nz_qe,rws,rd)
    !        call gen_dyn_qe(dyn_all,cell_qe,cell_sc_qe,cell_all,pos_qe,pos_all,&
    !                 idx_sc2uc,idx_all,mass_qe,mass_sc_qe,mass_all,fc_qe,nx_qe,ny_qe,nz_qe,&
    !                         1,1,nz,rws,rd)
            allocate(dyn_pc_qe(norb_pc_l,norb_pc_l))
            call gen_dyn_pc_qe(dyn_pc_qe,cell_qe,pos_qe,mass_qe,fc_qe,nx_qe,ny_qe,nz_qe,&
                                        rws,rd)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! broadcast non-array variable
        call MPI_BCAST(nlayer,1,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nd,1,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)
        call MPI_BCAST(maxlen,1,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)


        ! initialize Jacobian matrix
        allocate(jacobian(nk,ne,nlayer,nlayer))
        jacobian = 0.0d0
        allocate(jacobian_all(nk,ne,nlayer,nlayer))
        jacobian_all = 0.0d0
        allocate(jacobian_ij(nlayer,nlayer))
        jacobian_ij = 0.0d0
        ! initial flux array
        allocate(flux_i(nk,ne,nlayer))
        flux_i = 0.0d0 
        allocate(flux_all(nk,ne,nlayer))
        flux_all = 0.0d0 
        allocate(flux(nlayer))
        flux = 0.0d0
        allocate(flux_left(nk,ne))
        flux_left = 0.0d0 
        allocate(flux_right(nk,ne))
        flux_right = 0.0d0 
        allocate(flux_left_reduce(nk,ne))
        flux_left_reduce = 0.0d0 
        allocate(flux_right_reduce(nk,ne))
        flux_right_reduce = 0.0d0 
 
 
 
        ! spectral density 
        allocate(A_i(nk,ne,nlayer),A_all(nk,ne,nlayer))
        A_i = 0.0d0 
        A_all = 0.0d0
        allocate(Gn_i(nk,ne,nlayer),Gn_all(nk,ne,nlayer))
        Gn_i = 0.0d0
        Gn_all = 0.0d0 


        ! these are the matrices that consumes most memory 
        ! initialize left-connected Green function
        allocate(GL(numprocs,nlayer,norb_layer,norb_layer))
        ! initialize full Green function
        allocate(G_full(numprocs,nlayer,norb_layer,norb_layer))

        if (buttiker .eq. 1) then
            ! initialize left-connected lesser Green function GnL
            allocate(GnL(numprocs,nlayer,norb_layer,norb_layer))
            ! initialize lesser Green function Gn
            allocate(Gn(numprocs,nlayer,norb_layer,norb_layer))
            ! initialize temperature derivative 
            ! of left-connected lesser Green function GnL
            allocate(dGnLdT(numprocs,nlayer,nlayer,norb_layer,norb_layer))
            ! initialize temperature derivative of
            ! lesser Green function Gn
            allocate(dGndT(numprocs,nlayer,nlayer,norb_layer,norb_layer))
        else
            allocate(GnL(numprocs,1,1,1))
            allocate(Gn(numprocs,1,1,1))
            allocate(dGnLdT(numprocs,1,1,1,1))
            allocate(dGndT(numprocs,1,1,1,1))
        end if
        GnL = 0.0d0
        Gn = 0.0d0
        GL = 0.0d0
        G_full = 0.0d0
        dGnLdT = 0.0d0
        dGndT = 0.0d0

        ! initialize off-diagonal full Green function 
        ! G_q+1_q
        allocate(G_qp_q(numprocs,nlayer,norb_layer,norb_layer))
        ! G_q_q+1
        allocate(G_q_qp(numprocs,nlayer,norb_layer,norb_layer))
        ! G_N+1_N
        allocate(G_Np_N(numprocs,nr,norb_layer))
        ! G_N_N+1
        allocate(G_N_Np(numprocs,norb_layer,nr))
        ! G_1_0
        allocate(G_zp_z(numprocs,norb_layer,nl))
        ! G_0_1
        allocate(G_z_zp(numprocs,nl,norb_layer))
        G_qp_q = 0.0d0
        G_q_qp = 0.0d0
        G_Np_N = 0.0d0
        G_N_Np = 0.0d0
        G_zp_z = 0.0d0
        G_z_zp = 0.0d0
        ! tri-diagonal block device dyn
        allocate(Hamd(numprocs,nlayer,norb_layer,norb_layer))
        allocate(Hamdv(numprocs,nlayer,norb_layer,norb_layer))
        ! tri-diagonal block lead dyn
        allocate(Haml(numprocs,nl,nl),Hamr(numprocs,nr,nr))
        allocate(Hamlv(numprocs,nl,nl),Hamrv(numprocs,nr,nr))
        allocate(Hamld(numprocs,nl,norb_layer))
        allocate(Hamrd(numprocs,nr,norb_layer))

        Hamd = 0.0d0
        Hamdv = 0.0d0
        Haml = 0.0d0
        Hamr = 0.0d0
        Hamlv = 0.0d0
        Hamrv = 0.0d0
        Hamrd = 0.0d0
        Hamld = 0.0d0
    else
        call load_input()
        norb = natoms * 3
        norb_sc = natoms_sc * 3
        norb_uc = natoms_uc * 3
        allocate(force_constant(norb_uc,norb_sc))
        call read_fc(force_constant,norb_uc,norb_sc,filename_sandwich)
        call asr_cutoff(force_constant,positions_sc,r_cutoff,&
             n_fc_uc_x,n_fc_uc_y,n_fc_uc_z) 

    end if

    call k_grid_in_ws()

    if (qe_fc.eq. 1) then
        allocate(Ham(numprocs,1,1))
        allocate(eigvec(numprocs,1,1))
        allocate(eig(numprocs,1))
    else
        allocate(Ham(numprocs,norb,norb))
        allocate(eigvec(numprocs,norb,norb))
        allocate(eig(numprocs,norb))
    end if

    Ham = 0.0d0
    eigvec = 0.0d0
    eig = 0.0d0
    
    if (qe_fc .eq. 0) then

        call get_periodic_lead_orb_num() 

        norb_uc_l = nl/(n_bloch_x*n_bloch_y)
        norb_uc_r = nr/(n_bloch_x*n_bloch_y)

        allocate(force_constant_l(norb_uc_l,&
                 norb_uc_l*norb_sc/norb_uc))
        allocate(force_constant_r(norb_uc_r,&
                 norb_uc_r*norb_sc/norb_uc))

        call read_fc(force_constant_l,norb_uc_l,norb_uc_l*norb_sc/norb_uc,filename_left_lead)
        call read_fc(force_constant_r,norb_uc_r,norb_uc_r*norb_sc/norb_uc,filename_right_lead)
    end if

    if (qe_fc.ne.0) then
        allocate(Ham_l(1,1,1))
        allocate(Ham_r(1,1,1))
        allocate(Ham_l_v(1,1,1)) ! coupling between layers
        allocate(Ham_r_v(1,1,1)) ! coupling between layers
        allocate(eigvec_l(1,1,1))
        allocate(eigvec_r(1,1,1))
        allocate(eig_l(1,1))
        allocate(eig_r(1,1))
    else
        allocate(Ham_l(numprocs,nl,nl))
        allocate(Ham_r(numprocs,nr,nr))
        allocate(Ham_l_v(numprocs,nl,nl)) ! coupling between layers
        allocate(Ham_r_v(numprocs,nr,nr)) ! coupling between layers
        allocate(eigvec_l(numprocs,nl,nl))
        allocate(eigvec_r(numprocs,nr,nr))
        allocate(eig_l(numprocs,nl))
        allocate(eig_r(numprocs,nr))
    end if

    Ham_l = 0.0d0
    Ham_r = 0.0d0
    Ham_l_v = 0.0d0
    Ham_r_v = 0.0d0
    eigvec_l = 0.0d0
    eigvec_r = 0.0d0
    eig_l = 0.0d0
    eig_r = 0.0d0

    period_left_uc = period_left/n_bloch_x/n_bloch_y
    period_right_uc = period_right/n_bloch_x/n_bloch_y
    buffer_left_uc = buffer_left/n_bloch_x/n_bloch_y
    buffer_right_uc = buffer_right/n_bloch_x/n_bloch_y

    allocate(pos_left_lead(period_left,4),&
             pos_right_lead(period_right,4))
    allocate(pos_left_lead_uc(period_left_uc,4),&
             pos_left_lead_sc(period_left_uc*norb_sc/norb_uc,4))
    allocate(pos_right_lead_uc(period_right_uc,4),&
             pos_right_lead_sc(period_right_uc*norb_sc/norb_uc,4))
    if (qe_fc.eq.0) then

        pos_left_lead(:,:) = positions(buffer_left+1:buffer_left+period_left,:)
        pos_right_lead(:,:) = positions(natoms-buffer_right-period_right+1:natoms-buffer_right,:)
        pos_left_lead(:,3) = pos_left_lead(:,3) - pos_left_lead(1,3)
        pos_right_lead(:,3) = pos_right_lead(:,3) - pos_right_lead(1,3)

        pos_left_lead_uc(:,:) = positions_uc(buffer_left_uc+1:buffer_left_uc+period_left_uc,:)
        pos_right_lead_uc(:,:) = positions_uc(natoms_uc-buffer_right_uc-period_right_uc+1:natoms_uc-buffer_right_uc,:)
        pos_left_lead_uc(:,3) = pos_left_lead_uc(:,3) - pos_left_lead_uc(1,3)
        pos_right_lead_uc(:,3) = pos_right_lead_uc(:,3) - pos_right_lead_uc(1,3)

        call fc_sc_positions(pos_left_lead_uc,n_fc_uc_x,n_fc_uc_y,n_fc_uc_z,&
             pos_left_lead_sc,latvec_uc_l)
        call fc_sc_positions(pos_right_lead_uc,n_fc_uc_x,n_fc_uc_y,n_fc_uc_z,&
             pos_right_lead_sc,latvec_uc_r)

        call asr_cutoff(force_constant_l,pos_left_lead_sc,r_cutoff,&
             n_fc_uc_x,n_fc_uc_y,n_fc_uc_z) 
        call asr_cutoff(force_constant_r,pos_right_lead_sc,r_cutoff,&
             n_fc_uc_x,n_fc_uc_y,n_fc_uc_z) 


        ! primitive cell
        allocate(force_constant_pc_l(norb_pc_l,norb_pc_l*27))
        allocate(force_constant_pc_r(norb_pc_r,norb_pc_r*27))
       
        call extract_fc_pc(force_constant,natoms_pc_l,positions_pc_l,pos_pc_l_fc_sc,&
             latvec_pc_l,positions_sc,force_constant_pc_l,1)
        call asr_cutoff(force_constant_pc_l,pos_pc_l_fc_sc,r_cutoff,&
                1,1,1)
        call extract_fc_pc(force_constant,natoms_pc_r,positions_pc_r,pos_pc_r_fc_sc,&
             latvec_pc_r,positions_sc,force_constant_pc_r,2)
        call asr_cutoff(force_constant_pc_r,pos_pc_r_fc_sc,r_cutoff,&
                1,1,1)

        if ((path_mode.eq.1).or.(path_mode.eq.5))then
            allocate(dyn_uc_l(norb_uc_l,norb_uc_l),&
                     eig_uc_l(norb_uc_l),&
                     evs_uc_l(norb_uc_l,norb_uc_l))
            allocate(dyn_uc_r(norb_uc_r,norb_uc_r),&
                     eig_uc_r(norb_uc_r),&
                     evs_uc_r(norb_uc_r,norb_uc_r))
            allocate(dyn_pc_l(norb_pc_l,norb_pc_l),&
                     eig_pc_l(norb_pc_l),&
                     evs_pc_l(norb_pc_l,norb_pc_l))
            allocate(dyn_pc_r(norb_pc_r,norb_pc_r),&
                     eig_pc_r(norb_pc_r),&
                     evs_pc_r(norb_pc_r,norb_pc_r))
            allocate(vels_pc_l(norb_pc_l,3),&
                     vels_pc_r(norb_pc_r,3),&
                     vels_uc_l(norb_uc_l,3),&
                     vels_uc_r(norb_uc_r,3))


            do i = 1,n_kpath_points
                if (path_mode.eq.5) then

                ! This session writes eigenvalues of unitcell
                    if (crystal .eq. 1) then
                        kpath3(:) =  matmul(kpath_points(i,:),&
                                     reci_uc_l)
                    else
                        kpath3(:) = rotate_z1(kpath_points(i,:),pi/4.0d0)
                    end if

                    call gen_dyn_pc(norb_uc_l,&
                              matmul(kpath_points(i,:),&
                              reci_uc_l),&
                              dyn_uc_l,&
                              force_constant_l,positions_uc_l,&
                              pos_left_lead_sc)
                    call eigenH(dyn_uc_l,eig_uc_l,evs_uc_l)
                    call write_eig('eigenl',eig_uc_l)
                    call group_velocity_pc(norb_uc_l,& 
                                    kpath3,dyn_uc_l,&               
                                    force_constant_l,positions_uc_l,&
                                    pos_left_lead_sc,&      
                                    eig_uc_l,evs_uc_l,vels_uc_l) 
                    call write_vel('velsl',vels_uc_l)

                    if (crystal .eq. 1) then
                        kpath3(:) =  matmul(kpath_points(i,:),&
                                     reci_uc_r)
                    else
                        kpath3(:) = rotate_z1(kpath_points(i,:),pi/4.0d0)
                    end if

                    call gen_dyn_pc(norb_uc_r,&
                                  matmul(kpath_points(i,:),&
                                  reci_uc_r),&
                                  dyn_uc_r,&
                                  force_constant_r,positions_uc_r,&
                                  pos_right_lead_sc)
                    call eigenH(dyn_uc_r,eig_uc_r,evs_uc_r)
                    call write_eig('eigenr',eig_uc_r)
                    call group_velocity_pc(norb_uc_r,& 
                                    kpath3,dyn_uc_r,&               
                                    force_constant_r,positions_uc_r,&
                                    pos_right_lead_sc,&      
                                    eig_uc_r,evs_uc_r,vels_uc_r) 

                    call write_vel('velsr',vels_uc_r)

                end if

                ! This session writes eigenvalues of primitive cell
                if (path_mode.eq.1) then
     
                    if (crystal .eq. 1) then
                        kpath3(:) =  matmul(kpath_points(i,:),&
                                     recivec_pc_l)
                    else
                        kpath3(:) = rotate_z1(kpath_points(i,:),pi/4.0d0)
                    end if
                    call gen_dyn_pc(norb_pc_l,&
                                  kpath3,&
                                  dyn_pc_l,&
                                  force_constant_pc_l,positions_pc_l,&
                                  pos_pc_l_fc_sc)
                    call eigenH(dyn_pc_l,eig_pc_l,evs_pc_l)
                    call write_eig('eigenl',eig_pc_l)
                    call group_velocity_pc(norb_pc_l,& 
                                    kpath3,dyn_pc_l,&               
                                    force_constant_pc_l,positions_pc_l,&
                                    pos_pc_l_fc_sc,&      
                                    eig_pc_l,evs_pc_l,vels_pc_l) 
                    call write_vel('velsl',vels_pc_l)

                    
                    if (crystal .eq. 1) then
                        kpath3(:) =  matmul(kpath_points(i,:),&
                                     recivec_pc_r)
                    else
                        kpath3(:) = rotate_z1(kpath_points(i,:),pi/4.0d0)
                    end if
                    call gen_dyn_pc(norb_pc_r,&
                                  kpath3,&
                                  dyn_pc_r,&
                                  force_constant_pc_r,positions_pc_r,&
                                  pos_pc_r_fc_sc)
                    call eigenH(dyn_pc_r,eig_pc_r,evs_pc_r)
                    call write_eig('eigenr',eig_pc_r)
                    call group_velocity_pc(norb_pc_r,& 
                                    kpath3,dyn_pc_r,&               
                                    force_constant_pc_r,positions_pc_r,&
                                    pos_pc_r_fc_sc,&      
                                    eig_pc_r,evs_pc_r,vels_pc_r) 
                    call write_vel('velsr',vels_pc_r)
                end if
            end do
            stop
        end if
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

    if (qe_fc .eq. 0) then
        allocate(num_mat(norb/3,norb/3))
        allocate(idx_mat(norb/3,norb/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)),&
        dist_mat(norb/3,norb/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1),3))

        call gen_fc_to_dyn_idx(idx_mat,dist_mat,num_mat,force_constant,positions_sc,positions,&
                              latvec,n_bloch_x,n_bloch_y,norb)

        temp_sum = sum(sum(num_mat,1),1)
     
        allocate(idx_mat_red(temp_sum))
        allocate(dist_mat_red(temp_sum,3))
        idx_mat_red = 0
        dist_mat_red = 0.0d0

        temp_count = 1
        do i = 1,norb/3
            do j = 1,norb/3
                if (num_mat(i,j).gt.0)then
                    do k = 1, num_mat(i,j)
                        dist_mat_red(temp_count,1:3) = dist_mat(i,j,k,1:3)
                        idx_mat_red(temp_count) = idx_mat(i,j,k)
                        temp_count = temp_count + 1
                    end do
                end if
            end do
        end do
        deallocate(dist_mat,idx_mat)


        call gen_dyn_mpi(1,1,force_constant,Ham,eigvec,eig,norb,&
             positions_sc,positions,latvec,idx_mat_red,dist_mat_red,&
             num_mat)

        allocate(num_mat_l(nl/3,nl/3),num_mat_l_v(nl/3,nl/3))
        allocate(idx_mat_l(nl/3,nl/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z)),&
        idx_mat_l_v(nl/3,nl/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z)),&
        dist_mat_l(nl/3,nl/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1),3),&
        dist_mat_l_v(nl/3,nl/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1),3))


        call gen_lead_fc_to_dyn_idx(idx_mat_l,idx_mat_l_v,dist_mat_l,dist_mat_l_v,num_mat_l,&
                              num_mat_l_v,&
                              force_constant_l,pos_left_lead_sc,pos_left_lead,&
                              latvec_left,n_bloch_x,n_bloch_y,nl)

        temp_sum = sum(sum(num_mat_l,1),1)
        allocate(idx_mat_l_red(temp_sum))
        allocate(dist_mat_l_red(temp_sum,3))
        temp_sum = sum(sum(num_mat_l_v,1),1)
        allocate(idx_mat_l_v_red(temp_sum))
        allocate(dist_mat_l_v_red(temp_sum,3))

        idx_mat_l_red = 0
        idx_mat_l_v_red = 0
        dist_mat_l_red = 0.0d0
        dist_mat_l_v_red = 0.0d0

        temp_count = 1
        do i = 1,nl/3
            do j = 1,nl/3
                if (num_mat_l(i,j).gt.0)then
                    do k = 1, num_mat_l(i,j)
                        dist_mat_l_red(temp_count,1:3) = dist_mat_l(i,j,k,1:3)
                        idx_mat_l_red(temp_count) = idx_mat_l(i,j,k)
                        temp_count = temp_count + 1
                    end do
                end if
            end do
        end do
     
        temp_count = 1
        do i = 1,nl/3
            do j = 1,nl/3
                if (num_mat_l_v(i,j).gt.0)then
                    do k = 1, num_mat_l_v(i,j)
                        dist_mat_l_v_red(temp_count,1:3) = dist_mat_l_v(i,j,k,1:3)
                        idx_mat_l_v_red(temp_count) = idx_mat_l_v(i,j,k)
                        temp_count = temp_count + 1
                    end do
                end if
            end do
        end do
      
        deallocate(dist_mat_l,dist_mat_l_v,idx_mat_l,idx_mat_l_v)
     
        call gen_dyn_lead_mpi(1,1,force_constant_l,Ham_l,Ham_l_v,eigvec_l,eig_l,nl,&
             pos_left_lead_sc,pos_left_lead,latvec_left,idx_mat_l_red,idx_mat_l_v_red,&
             dist_mat_l_red,dist_mat_l_v_red,num_mat_l,num_mat_l_v)
        

        allocate(num_mat_r(nr/3,nr/3),num_mat_r_v(nr/3,nr/3))
        allocate(idx_mat_r(nr/3,nr/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z)),&
        idx_mat_r_v(nr/3,nr/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z)),&
        dist_mat_r(nr/3,nr/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1),3),&
        dist_mat_r_v(nr/3,nr/3,(2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1),3))


        call gen_lead_fc_to_dyn_idx(idx_mat_r,idx_mat_r_v,dist_mat_r,dist_mat_r_v,num_mat_r,&
                              num_mat_r_v,&
                              force_constant_r,pos_right_lead_sc,pos_right_lead,&
                              latvec_right,n_bloch_x,n_bloch_y,nr)

        temp_sum = sum(sum(num_mat_r,1),1)
        allocate(idx_mat_r_red(temp_sum))
        allocate(dist_mat_r_red(temp_sum,3))
        temp_sum = sum(sum(num_mat_r_v,1),1)
        allocate(idx_mat_r_v_red(temp_sum))
        allocate(dist_mat_r_v_red(temp_sum,3))

        idx_mat_r_red = 0
        idx_mat_r_v_red = 0
        dist_mat_r_red = 0.0d0
        dist_mat_r_v_red = 0.0d0

        temp_count = 1
        do i = 1,nr/3
            do j = 1,nr/3
                if (num_mat_r(i,j).gt.0)then
                    do k = 1, num_mat_r(i,j)
                        dist_mat_r_red(temp_count,1:3) = dist_mat_r(i,j,k,1:3)
                        idx_mat_r_red(temp_count) = idx_mat_r(i,j,k)
                        temp_count = temp_count + 1
                    end do
                end if
            end do
        end do
        
        temp_count = 1
        do i = 1,nr/3
            do j = 1,nr/3
                if (num_mat_r_v(i,j).gt.0)then
                    do k = 1, num_mat_r_v(i,j)
                        dist_mat_r_v_red(temp_count,1:3) = dist_mat_r_v(i,j,k,1:3)
                        idx_mat_r_v_red(temp_count) = idx_mat_r_v(i,j,k)
                        temp_count = temp_count + 1
                    end do
                end if
            end do
        end do
     
        deallocate(dist_mat_r,dist_mat_r_v,idx_mat_r,idx_mat_r_v)

        call gen_dyn_lead_mpi(1,1,force_constant_r,Ham_r,Ham_r_v,eigvec_r,eig_r,nr,&
             pos_right_lead_sc,pos_right_lead,latvec_right,idx_mat_r_red,idx_mat_r_v_red,&
             dist_mat_r_red,dist_mat_r_v_red,num_mat_r,num_mat_r_v)

        call swap_lead_region_dyn_mpi(1,1,Ham,Ham_l,Ham_r,Ham_l_v,Ham_r_v)

        if (myid .eq. 0) then
            call blockwise_setup(Ham(1,:,:),norb)
            allocate(temperature_i(nlayer))
            temperature_i(:) = t_left
        end if
    end if
    
    
    if (myid .eq. 0) then
        call transverse_uc_setup()
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! broadcast non-array variable
    call MPI_BCAST(natm_l_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(natm_r_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nl_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(nr_uc,1,MPI_INTEGER, 0, &
         MPI_COMM_WORLD, ierr)

    if (qe_fc.eq.0) then
        call MPI_BCAST(nlayer,1,MPI_INTEGER,0,&
        MPI_COMM_WORLD, ierr) 
    end if

    ! still we need to run allocation on other cores.
    if (myid .gt. 0) then
        allocate(layer_list(nlayer),layerstart(nlayer),&
        layerend(nlayer),pidxL(nl),pidxR(nr),&
        L_orb_num_list(natm_l_uc),R_orb_num_list(natm_r_uc),&
        L_list(natm_l_uc),R_list(natm_r_uc),&
        pos_l(nl,3),pos_r(nr,3),temperature_i(nlayer))
    end if

    ! finally we broadcast array variables
    call MPI_BCAST(layer_list,nlayer,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(layerstart,nlayer,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(layerend,nlayer,MPI_INTEGER,0,&
         MPI_COMM_WORLD, ierr)
    call MPI_BCAST(temperature_i,nlayer,MPI_DOUBLE_PRECISION,0,&
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

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
   
    ! initialize buttiker self-energy
    allocate(sigma_w(nlayer,ne))
    sigma_w = 0.0d0
    allocate(dsigma_wdt(nlayer,ne))
    dsigma_wdt = 0.0d0
 
    call set_sigma_bp_w(temperature_i)

    if (path_mode .eq. 4)  then
        latvec_slab = 0.0d0
        recvec_slab = 0.0d0
        latvec_slab(1:2,1:2) = latvec(1:2,1:2) 
        latvec_slab(3,3) = latvec(3,3)*5 ! large vaccuum
        recvec_slab = 2.0d0* pi * transpose(inv_real(latvec_slab))
        natoms_slab_l = 2*period_left+buffer_left
        natoms_slab_r = 2*period_right+buffer_right
        norb_slab_l = norb_uc_l*nthick*n_bloch_x*n_bloch_y
        norb_slab_r = norb_uc_r*nthick*n_bloch_x*n_bloch_y
        allocate(dyn_slab_l(norb_slab_l,norb_slab_l),&
                 eig_slab_l(norb_slab_l),&
                 evs_slab_l(norb_slab_l,norb_slab_l),&
                 pos_slab_l(natoms_slab_l,4))
        allocate(dyn_slab_r(norb_slab_r,norb_slab_r),&
                 eig_slab_r(norb_slab_r),&
                 evs_slab_r(norb_slab_r,norb_slab_r),&
                 pos_slab_r(natoms_slab_r,4))
        pos_slab_l = 0.0d0
        pos_slab_l(1:natoms_slab_l,:) = positions(1:natoms_slab_l,:) 
        pos_slab_r = 0.0d0
        pos_slab_r(1:natoms_slab_r,:) = positions(natoms-natoms_slab_r+1:natoms,:) 
        do i = 1,n_kpath_points
            call gen_ham_slab(norb_slab_l,&
                           matmul(kpath_points(i,:),recvec_slab),&
                           dyn_slab_l,&
                           eig_slab_l,evs_slab_l,&
                           latvec_slab,force_constant_l,&
                           pos_left_lead_sc,pos_left_lead,&
                           latvec_left,idx_mat_l_red,idx_mat_l_v_red,dist_mat_l_red,&
                           dist_mat_l_v_red,num_mat_l,num_mat_l_v,&
                           nl) 
            call gen_ham_slab(norb_slab_r,&
                           matmul(kpath_points(i,:),recvec_slab),&
                           dyn_slab_r,&
                           eig_slab_r,evs_slab_r,&
                           latvec_slab,force_constant_r,&
                           pos_right_lead_sc,pos_right_lead,&
                           latvec_right,idx_mat_r_red,idx_mat_r_v_red,dist_mat_r_red,&
                           dist_mat_r_v_red,num_mat_r,num_mat_r_v,&
                           nr) 
            call write_eig('slabl',eig_slab_l)
            call write_eig('slabr',eig_slab_r)
        end do
        stop
    end if


    if (path_mode .ne. 6) then
        if (qe_fc.eq.0) then
            do ik = 1,nkcore
                call gen_dyn_mpi(1,ik,force_constant,Ham,eigvec,eig,norb,&
                     positions_sc,positions,latvec,idx_mat_red,dist_mat_red,&
                     num_mat)

                call gen_dyn_lead_mpi(1,ik,force_constant_l,Ham_l,Ham_l_v,eigvec_l,eig_l,nl,&
                     pos_left_lead_sc,pos_left_lead,latvec_left,idx_mat_l_red,idx_mat_l_v_red,&
                     dist_mat_l_red,dist_mat_l_v_red,num_mat_l,num_mat_l_v)
         
                call gen_dyn_lead_mpi(1,ik,force_constant_r,Ham_r,Ham_r_v,eigvec_r,eig_r,nr,&
                     pos_right_lead_sc,pos_right_lead,latvec_right,idx_mat_r_red,idx_mat_r_v_red,&
                     dist_mat_r_red,dist_mat_r_v_red,num_mat_r,num_mat_r_v)

                !call swap_lead_region_dyn_mpi(1,ik,Ham,Ham_l,Ham_r,Ham_l_v,Ham_r_v)

                do ie = 1,ne
                    call surface_green_mpi(ie,ik,Ham,norb,&
                                          gl_adv,gr_ret,gl_plus,gr_minus,&
                                          gml_adv,gmr_ret,&
                                          sl_adv,sr_ret,Egrid)
                    call surface_ldos_mpi(ie,ik,&
                                       gl_adv,gr_ret,&
                                       surface_dos_l_reduce,surface_dos_r_reduce,Egrid)

                    call G_block_mpi(ie,ik,Ham,norb,&
                               gl_adv,gr_ret,&
                               gml_adv,gmr_ret,sl_adv,sr_ret,Egrid,&
                               transl_reduce(:,:),transr_reduce(:,:),&
                               G_D_RL(:,:,:),G_D_LR(:,:,:),&
                               G_L(:,:,:),G_R(:,:,:))

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
                    end if
                end do
            end do
            call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        
            if (myid.eq.0) then
                call write_surface_ldos(surface_dos_l,surface_dos_r) 
            end if
        else if ((qe_fc .eq. 2).and.(buttiker.eq.1)) then
            ! qe force constant
            ! temperature difference compared with last loop
            tdiff = 100.0d0
            intstep = ne/10

            outer: do iit = 1, 100
                jacobian = 0.0d0
                jacobian_all = 0.0d0
                jacobian_ij = 0.0d0
                flux_i = 0.0d0
                flux_all = 0.0d0
                flux = 0.0d0
                Gn_i = 0.0d0
                Gn_all = 0.0d0
                A_i = 0.0d0
                A_all = 0.0d0
                Ham = 0.0d0
                GnL = 0.0d0
                Gn = 0.0d0
                GL = 0.0d0
                G_full = 0.0d0
                transl_s_reduce = 0.0d0
                transl_ns_reduce = 0.d0
                transr_s_reduce = 0.0d0
                transr_ns_reduce = 0.0d0
                transl_reduce = 0.0d0
                transr_reduce = 0.0d0
                refll_s_reduce = 0.0d0
                refll_ns_reduce = 0.0d0
                reflr_s_reduce = 0.0d0
                reflr_ns_reduce = 0.0d0
                do ik = 1,nkcore
                    call gen_dyn_qe_mpi(1,ik,norb,Ham,cell_qe,cell_sc_qe,cell_all,pos_qe,pos_all,&
                     idx_sc2uc,idx_all,mass_qe,mass_sc_qe,mass_all,fc_qe,nx_qe,ny_qe,nz_qe,&
                             n_bloch_x,n_bloch_y,nz,rws,rd,neighborlist,nnbrs)
                    do ie = 1,ne
                        call surface_green_mpi(ie,ik,Ham,norb,&
                                              gl_adv,gr_ret,gl_plus,gr_minus,&
                                              gml_adv,gmr_ret,&
                                              sl_adv,sr_ret,Egrid)
                        call surface_ldos_mpi(ie,ik,&
                                           gl_adv,gr_ret,&
                                           surface_dos_l_reduce,surface_dos_r_reduce,Egrid)

                        !call G_block_mpi(ie,ik,Ham,norb,&
                        !           gl_adv,gr_ret,&
                        !           gml_adv,gmr_ret,sl_adv,sr_ret,Egrid,&
                        !           transl_reduce(:,:),transr_reduce(:,:),&
                        !           G_D_RL(:,:,:),G_D_LR(:,:,:),&
                        !           G_L(:,:,:),G_R(:,:,:))
                        call get_G_full_mpi(ie,ik,Ham,norb,&
                                   gl_adv,gr_ret,&
                                   gml_adv,gmr_ret,sl_adv,sr_ret,&
                                   Egrid,GL,G_full,&
                                   G_qp_q,G_q_qp,G_Np_N,G_N_Np,&
                                   G_zp_z,G_z_zp,&
                                   transl_reduce(:,:),transr_reduce(:,:),&
                                   G_D_RL(:,:,:),G_D_LR(:,:,:),&
                                   G_L(:,:,:),G_R(:,:,:))
                        call get_G_lesser_mpi(ie,ik,Ham,norb,&
                                   gl_adv,gr_ret,&
                                   gml_adv,gmr_ret,sl_adv,sr_ret,&
                                   Egrid,GL,G_full,&
                                   G_qp_q,G_q_qp,G_Np_N,G_N_Np,&
                                   G_zp_z,G_z_zp,&
                                   GnL,Gn,dGnLdT,dGndT,flux_i,&
                                   Gn_i,A_i)

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
                 !          intprogress = dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0
                 !          if (abs(int(intprogress)-intprogress).lt. 1.0d-2) then
                 !               write(*,23,ADVANCE='NO') intprogress
                 !               do i = 1,50
                 !                   if ((i-1)*100/50 .le.&
                 !                       int(intprogress)) then
                 !                    if (mod(int(intprogress),5).eq.0) then
                                    if (mod(ie,intstep).eq. 0) then
                                        write(*,"(A)",ADVANCE='NO')'.'
                                    end if
                 !                    end if
                 !                   else
                 !                       write(*,"(A)",ADVANCE='NO')' '
                  !                  end if
                  !              end do
                  !              write(*,"(A)")']'
                                23  format ('Process: ',f6.2,' % [')
                 !           end if
                            if ((ie.eq.ne).and.(ik.eq.nkcore)) then
                                write(*,"(A)") ' '
                            end if
                        end if
                    end do
                end do

                call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
                
                call MPI_ALLREDUCE(jacobian,jacobian_all,ne*nk*nlayer*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(flux_i,flux_all,ne*nk*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(Gn_i,Gn_all,ne*nk*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(A_i,A_all,ne*nk*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
     
                call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            
                if (myid.eq.0) then
                     ! write sdos
                     !call write_surface_ldos(surface_dos_l,surface_dos_r) 
                     call integrate_jacobian(jacobian_ij,jacobian_all,flux_all,flux)
                     call write_flux(flux_all,gn_all,a_all,flux_left,flux_right)
                     call get_ph_temperature(sum(A_all,1),sum(Gn_all,1),&
                          temperature_ph)
                     call update_temperature(jacobian_ij,flux,tdiff)
                     call write_temperature(temperature_i,temperature_ph)
                 end if

                 call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
                 
                 ! broadcast to other cores
                 call MPI_BCAST(temperature_i,nlayer,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD, ierr)
                 call MPI_BCAST(tdiff,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD, ierr)
                 ! update scattering rate
                 call set_sigma_bp_w(temperature_i)


                if (tdiff .lt. 1.0d-4) exit outer

             end do outer
         else 
            ! qe force constant
            ! temperature difference compared with last loop
            tdiff = 100.0d0
            intstep = (nkcore*ne)/5
            if (buttiker.eq.0) then
                niit = 1
            else
                niit = 100
            end if
            outerr: do iit = 1, niit
                jacobian = 0.0d0
                jacobian_all = 0.0d0
                jacobian_ij = 0.0d0
                flux_i = 0.0d0
                flux_all = 0.0d0
                flux = 0.0d0
                flux_left = 0.0d0
                flux_left_reduce = 0.0d0
                flux_right = 0.0d0
                flux_right_reduce = 0.0d0  
                Gn_i = 0.0d0
                Gn_all = 0.0d0
                A_i = 0.0d0
                A_all = 0.0d0
                Ham = 0.0d0
                GnL = 0.0d0
                Gn = 0.0d0
                GL = 0.0d0
                G_full = 0.0d0
                transl_s_reduce = 0.0d0
                transl_ns_reduce = 0.d0
                transr_s_reduce = 0.0d0
                transr_ns_reduce = 0.0d0
                transl_reduce = 0.0d0
                transr_reduce = 0.0d0
                refll_s_reduce = 0.0d0
                refll_ns_reduce = 0.0d0
                reflr_s_reduce = 0.0d0
                reflr_ns_reduce = 0.0d0
                do ik = 1,nkcore
                      call gen_dyn_qe_tdb_mpi(1,ik,cell_qe,cell_sc_qe,cell_all,pos_qe,pos_all,&
                      idx_sc2uc,idx_all,mass_qe,mass_sc_qe,mass_all,fc_qe,nx_qe,ny_qe,nz_qe,&
                             n_bloch_x,n_bloch_y,nz,rws,rd,neighborlist,nnbrs,&
                             Hamd,Hamdv,Haml,Hamr,&
                             Hamlv,Hamrv,&
                             Hamld,Hamrd,&
                             nl,nr,nd,nlayer,norb_layer)        
                    do ie = 1,ne
                        call surface_green_tdb_mpi(ie,ik,&
                                              gl_adv,gr_ret,gl_plus,gr_minus,&
                                              gml_adv,gmr_ret,&
                                              sl_adv,sr_ret,Egrid,&
                                              Haml,Hamr,&
                                              Hamlv,Hamrv)
                        !call surface_ldos_mpi(ie,ik,&
                        !                   gl_adv,gr_ret,&
                        !                   surface_dos_l_reduce,surface_dos_r_reduce,Egrid)

                        !call G_block_mpi(ie,ik,Ham,norb,&
                        !           gl_adv,gr_ret,&
                        !           gml_adv,gmr_ret,sl_adv,sr_ret,Egrid,&
                        !           transl_reduce(:,:),transr_reduce(:,:),&
                        !           G_D_RL(:,:,:),G_D_LR(:,:,:),&
                        !           G_L(:,:,:),G_R(:,:,:))
                        call get_G_full_tdb_mpi(ie,ik,&
                                   gl_adv,gr_ret,&
                                   gml_adv,gmr_ret,sl_adv,sr_ret,&
                                   Egrid,GL,G_full,&
                                   G_qp_q,G_q_qp,G_Np_N,G_N_Np,&
                                   G_zp_z,G_z_zp,&
                                   transl_reduce(:,:),transr_reduce(:,:),&
                                   G_D_RL(:,:,:),G_D_LR(:,:,:),&
                                   G_L(:,:,:),G_R(:,:,:),&
                                   hamd,hamdv,hamld,hamrd,&
                                   haml,hamr)
                        if (buttiker.eq.1) then
                            call get_G_lesser_tdb_mpi(ie,ik,&
                                   gl_adv,gr_ret,&
                                   gml_adv,gmr_ret,sl_adv,sr_ret,&
                                   Egrid,GL,G_full,&
                                   G_qp_q,G_q_qp,G_Np_N,G_N_Np,&
                                   G_zp_z,G_z_zp,&
                                   GnL,Gn,dGnLdT,dGndT,flux_i,&
                                   Gn_i,A_i,&
                                   hamd,hamdv,hamld,hamrd,&
                                   flux_left_reduce,&
                                   flux_right_reduce)
                        end if

                        call vel_mat_tdb_mpi(ie,ik,Egrid,&
                             gl_adv,gr_ret,gl_plus,gr_minus,&
                             gml_adv,gmr_ret,&
                             G_D_RL(:,:,:),G_D_LR(:,:,:),&
                             G_L(:,:,:),G_R(:,:,:),&
                             transl_s_reduce(:,:),transr_s_reduce(:,:),&
                             transl_ns_reduce(:,:),transr_ns_reduce(:,:),&
                             refll_s_reduce(:,:),reflr_s_reduce(:,:),&
                             refll_ns_reduce(:,:),reflr_ns_reduce(:,:),&
                             norb_pc_l,norb_pc_r,&
                             norb_uc_l,norb_uc_r,&
                             hamd,hamdv,haml,hamr,&
                             hamld,hamrd,hamlv,hamrv)

                        if (myid.eq.0) then
                 !          intprogress = dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0
                 !          if (abs(int(intprogress)-intprogress).lt. 1.0d-2) then
                 !               write(*,23,ADVANCE='NO') intprogress
                 !               do i = 1,50
                 !                   if ((i-1)*100/50 .le.&
                 !                       int(intprogress)) then
                 !                    if (mod(int(intprogress),5).eq.0) then
                                    if (mod((ik-1)*ne+ie,intstep).eq. 0) then
                                        write(*,"(A)",ADVANCE='YES')'.'
                                    end if
                 !                    end if
                 !                   else
                 !                       write(*,"(A)",ADVANCE='NO')' '
                  !                  end if
                  !              end do
                  !              write(*,"(A)")']'
                 !           end if
                            if ((ie.eq.ne).and.(ik.eq.nkcore)) then
                                write(*,"(A)") ' '
                            end if
                        end if
                    end do
                end do

                call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
                
                call MPI_ALLREDUCE(jacobian,jacobian_all,ne*nk*nlayer*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(flux_i,flux_all,ne*nk*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(flux_left_reduce,flux_left,ne*nk,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
                    
                call MPI_ALLREDUCE(flux_right_reduce,flux_right,ne*nk,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(Gn_i,Gn_all,ne*nk*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(A_i,A_all,ne*nk*nlayer,&
                    MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
     
                call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

                call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
            
                if ((myid.eq.0) .and. (buttiker.eq.1)) then
                     ! write sdos
                     !call write_surface_ldos(surface_dos_l,surface_dos_r) 
                     call integrate_jacobian(jacobian_ij,jacobian_all,flux_all,flux) 


                     call write_flux(flux_all,gn_all,a_all,flux_left,flux_right)
                    
                    
                     call get_ph_temperature(sum(A_all,1)/dble(nk),sum(Gn_all,1)/dble(nk),&
                          temperature_ph)
 
                     call update_temperature(jacobian_ij,flux,tdiff)
                     
                     call write_temperature(temperature_i,temperature_ph)

                 end if

                 call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
                 ! broadcast to other cores
                 call MPI_BCAST(temperature_i,nlayer,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD, ierr)
                 call MPI_BCAST(tdiff,1,MPI_DOUBLE_PRECISION,0,&
                      MPI_COMM_WORLD, ierr)
                 ! update scattering rate
                 call set_sigma_bp_w(temperature_i)


                if (tdiff .lt. 1.0d-4) exit outerr

             end do outerr
!        else
!            intstep = ne/10
!                do ik = 1,nkcore
!                    call gen_dyn_qe_mpi(1,ik,norb,Ham,cell_qe,cell_sc_qe,cell_all,pos_qe,pos_all,&
!                     idx_sc2uc,idx_all,mass_qe,mass_sc_qe,mass_all,fc_qe,nx_qe,ny_qe,nz_qe,&
!                             n_bloch_x,n_bloch_y,nz,rws,rd,neighborlist,nnbrs)
!
!                    do ie = 1,ne
!                        call surface_green_mpi(ie,ik,Ham,norb,&
!                                              gl_adv,gr_ret,gl_plus,gr_minus,&
!                                              gml_adv,gmr_ret,&
!                                              sl_adv,sr_ret,Egrid)
!                        call surface_ldos_mpi(ie,ik,&
!                                           gl_adv,gr_ret,&
!                                           surface_dos_l_reduce,surface_dos_r_reduce,Egrid)
!
!                        call G_block_mpi(ie,ik,Ham,norb,&
!                                   gl_adv,gr_ret,&
!                                   gml_adv,gmr_ret,sl_adv,sr_ret,Egrid,&
!                                   transl_reduce(:,:),transr_reduce(:,:),&
!                                   G_D_RL(:,:,:),G_D_LR(:,:,:),&
!                                   G_L(:,:,:),G_R(:,:,:))
!                        !call get_G_full_mpi(ie,ik,Ham,norb,&
!                        !           gl_adv,gr_ret,&
!                        !           gml_adv,gmr_ret,sl_adv,sr_ret,&
!                        !           Egrid,GL,G_full,&
!                        !           G_qp_q,G_q_qp,G_Np_N,G_N_Np,&
!                        !           G_zp_z,G_z_zp,&
!                        !           transl_reduce(:,:),transr_reduce(:,:),&
!                        !           G_D_RL(:,:,:),G_D_LR(:,:,:),&
!                        !           G_L(:,:,:),G_R(:,:,:))
!                        !call get_G_lesser_mpi(ie,ik,Ham,norb,&
!                        !           gl_adv,gr_ret,&
!                        !           gml_adv,gmr_ret,sl_adv,sr_ret,&
!                        !           Egrid,GL,G_full,&
!                        !           G_qp_q,G_q_qp,G_Np_N,G_N_Np,&
!                        !           G_zp_z,G_z_zp,&
!                        !           GnL,Gn,dGnLdT,dGndT,flux_i,&
!                        !           Gn_i,A_i)
!
!                        call vel_mat_mpi(ie,ik,Egrid,Ham,norb,&
!                             gl_adv,gr_ret,gl_plus,gr_minus,&
!                             gml_adv,gmr_ret,&
!                             G_D_RL(:,:,:),G_D_LR(:,:,:),&
!                             G_L(:,:,:),G_R(:,:,:),&
!                             transl_s_reduce(:,:),transr_s_reduce(:,:),&
!                             transl_ns_reduce(:,:),transr_ns_reduce(:,:),&
!                             refll_s_reduce(:,:),reflr_s_reduce(:,:),&
!                             refll_ns_reduce(:,:),reflr_ns_reduce(:,:),&
!                             norb_pc_l,norb_pc_r,&
!                             norb_uc_l,norb_uc_r)
!
!
!                        if (myid.eq.0) then
!                 !          intprogress = dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0
!                 !          if (abs(int(intprogress)-intprogress).lt. 1.0d-2) then
!                 !               write(*,23,ADVANCE='NO') intprogress
!                 !               do i = 1,50
!                 !                   if ((i-1)*100/50 .le.&
!                 !                       int(intprogress)) then
!                 !                    if (mod(int(intprogress),5).eq.0) then
!                                    if (mod(ie,intstep).eq. 0) then
!                                        write(*,"(A)",ADVANCE='NO')'.'
!                                    end if
!                 !                    end if
!                 !                   else
!                 !                       write(*,"(A)",ADVANCE='NO')' '
!                  !                  end if
!                  !              end do
!                  !              write(*,"(A)")']'
!                 !           end if
!                            if ((ie.eq.ne).and.(ik.eq.nkcore)) then
!                                write(*,"(A)") ' '
!                            end if
!                        end if
!                    end do
!                end do
!
!                call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!                     
!                call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
!                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!
!                call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
!                 MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
!            
!                if (myid.eq.0) then
!                     ! write sdos
!                     call write_surface_ldos(surface_dos_l,surface_dos_r) 
!                 end if
!                 ! broadcast to other cores
        end if
    else
        do ik = 1,nkcore
            call gen_dyn_mpi(1,ik,force_constant,Ham,eigvec,eig,norb,&
                 positions_sc,positions,latvec,idx_mat_red,dist_mat_red,&
                 num_mat)

            call gen_dyn_lead_mpi(1,ik,force_constant_l,Ham_l,Ham_l_v,eigvec_l,eig_l,nl,&
                 pos_left_lead_sc,pos_left_lead,latvec_left,idx_mat_l_red,idx_mat_l_v_red,&
                 dist_mat_l_red,dist_mat_l_v_red,num_mat_l,num_mat_l_v)
     
            call gen_dyn_lead_mpi(1,ik,force_constant_r,Ham_r,Ham_r_v,eigvec_r,eig_r,nr,&
                 pos_right_lead_sc,pos_right_lead,latvec_right,idx_mat_r_red,idx_mat_r_v_red,&
                 dist_mat_r_red,dist_mat_r_v_red,num_mat_r,num_mat_r_v)

            call swap_lead_region_dyn_mpi(1,ik,Ham,Ham_l,Ham_r,Ham_l_v,Ham_r_v)

            do ie = 1,ne
                call surface_green_mpi(ie,ik,Ham,norb,&
                                      gl_adv,gr_ret,gl_plus,gr_minus,&
                                      gml_adv,gmr_ret,&
                                      sl_adv,sr_ret,Egrid)
                call surface_ldos_mpi(ie,ik,&
                                   gl_adv,gr_ret,&
                                   surface_dos_l_reduce,surface_dos_r_reduce,Egrid)

                call G_block_mpi(ie,ik,Ham,norb,&
                           gl_adv,gr_ret,&
                           gml_adv,gmr_ret,sl_adv,sr_ret,Egrid,&
                           transl_reduce(:,:),transr_reduce(:,:),&
                           G_D_RL(:,:,:),G_D_LR(:,:,:),&
                           G_L(:,:,:),G_R(:,:,:))

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
                    write(*,24,ADVANCE='NO') dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0
                    do i = 1,50
                        if ((i-1)*100/50 .le.&
                            int(dble((ik-1)*ne+ie)/dble(nkcore*ne)*100.0d0)) then
                            write(*,"(A)",ADVANCE='NO')'.'
                        else
                            write(*,"(A)",ADVANCE='NO')' '
                        end if
                    end do
                    write(*,"(A)")']'
                    24  format ('Process: ',f6.2,' % [')
                end if
            end do
        end do
        call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

        call MPI_ALLREDUCE(surface_dos_l_reduce,surface_dos_l,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(surface_dos_r_reduce,surface_dos_r,ne*nk,&
         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    
         if (myid.eq.0) then
             call write_surface_ldos(surface_dos_l,surface_dos_r) 
         end if
    end if    

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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
            reflr_s,reflr_ns)
        end if
    end if

    call MPI_FINALIZE(ierr)

end program phinfortran
