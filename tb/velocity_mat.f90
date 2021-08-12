module velocity_mat
      
    implicit none
    
    integer(kind=4),allocatable :: L_list(:),R_list(:)
    integer(kind=4),allocatable :: L_orb_num_list(:), R_orb_num_list(:)
    real(kind=8),allocatable :: pos_l(:,:), pos_r(:,:)
    integer(kind=4),allocatable :: pidxL(:),pidxR(:)
    integer(kind=4) :: natm_l_uc, natm_r_uc
    integer(kind=4) :: nl_uc,nr_uc

contains

    subroutine transverse_uc_setup()

    use surface, only : nl,nr,n_buffer_l,n_buffer_r  
    use param
    use input
    use util
    use config

    implicit none
        
    integer(kind=4) :: i,j,k,l,countL,countR,temp
    real(kind=8) :: n(3)
    real(kind=8) :: pos(3),p(3),dp(3)
    real(kind=8) :: poszl, poszr
    integer(kind=4),allocatable :: tempL(:),tempR(:),&
                                   tempL1(:),tempR1(:)
    allocate(L_list(1),R_list(1),&
             L_orb_num_list(1),R_orb_num_list(1),&
             tempL(1),tempR(1),tempL1(1),tempR1(1))
    
    L_list = 0 ! the indices for transverse uc for left lead
    L_orb_num_list = 0 ! corresponding number of orbs
    R_list = 0 ! the indices for transverse uc for right lead
    R_orb_num_list = 0 ! corresponding number of orbs
    tempL = 0
    tempR = 0
    tempL1 = 0
    tempR1 = 0

    countL = 1
    nl_uc = 0
    do i = buffer_left+period_left+1,buffer_left+2*period_left
        pos = positions(i,1:3)
        n = matmul(pos,inv_real(latvec_uc))
        if ((n(1) .lt. 0.999d0).and.(n(2) .lt. 0.999d0) .and.&
            (n(1) .gt. -0.001d0).and.(n(2) .gt. -0.001d0)) then
            deallocate(L_list)
            deallocate(L_orb_num_list)
            allocate(L_list(countL))
            allocate(L_orb_num_list(countL))
            L_list(countL) = i
            L_orb_num_list(countL) = atom_orb_num(i)
            if (countL .gt. 1) then
                L_list(1:countL-1) = tempL(:)
                L_orb_num_list(1:countL-1) = tempL1(:)
            end if
            deallocate(tempL)
            allocate(tempL(countL))
            tempL = L_list
            deallocate(tempL1)
            allocate(tempL1(countL))
            tempL1 = L_orb_num_list

            countL = countL + 1

            nl_uc = nl_uc + atom_orb_num(i)
        end if
    end do
    natm_l_uc = countL-1 ! number of atoms in transverse uc

    countR = 1
    nr_uc = 0
    do i = natoms-buffer_right-2*period_right+1,&
           natoms-buffer_right-period_right
        pos = positions(i,1:3)
        n = matmul(pos,inv_real(latvec_uc))
        if ((n(1) .lt. 0.999d0).and.(n(2) .lt. 0.999d0) .and.&
            (n(1) .gt. -0.001d0).and.(n(2) .gt. -0.001d0)) then
            deallocate(R_list)
            deallocate(R_orb_num_list)
            allocate(R_list(countR))
            allocate(R_orb_num_list(countR))
            R_list(countR) = i
            R_orb_num_list(countR) = atom_orb_num(i)
            if (countR .gt. 1) then
                R_list(1:countR-1) = tempR(:)
                R_orb_num_list(1:countR-1) = tempR1(:)
            end if
            deallocate(tempR)
            allocate(tempR(countR))
            tempR = R_list
            deallocate(tempR1)
            allocate(tempR1(countR))
            tempR1 = R_orb_num_list
            
            countR = countR + 1

            nr_uc = nr_uc + atom_orb_num(i)
        end if
    end do
    natm_r_uc = countR-1 ! number of atoms in transverse uc

    allocate(pos_l(nl,4))
    allocate(pos_r(nr,4))
    allocate(pidxL(nl))
    allocate(pidxR(nr))
    pos_l = 0.0d0
    pos_r = 0.0d0
    poszl = minval(positions(buffer_left+period_left+1:&
                  buffer_left+2*period_left,3)) 
    poszr = minval(positions(&
                  natoms-buffer_right-2*period_right+1:&
                  natoms-buffer_right-period_right,3)) 

    countL = 1
    do i = buffer_left+period_left+1,buffer_left+2*period_left
        pos = positions(i,1:3)
        do j = 1,natm_l_uc
            p = positions(L_list(j),1:3)
            do k = 1,n_bloch_x
                do l = 1,n_bloch_y
                    dp(:) = pos(:)-p(:)&
                    -matmul((/dble(k-1),dble(l-1),0.0d0/),latvec_uc(:,:))
                    if (sqrt(dp(1)**2+dp(2)**2+dp(3)**2).lt.1.0d-3) then
                        temp = j 
                    end if
                end do
            end do
        end do
        ! index with same orbitals but different positions
        do j = 1,atom_orb_num(i)
            pos_l(countL+j-1,1:3) = pos(:)
            pos_l(countL+j-1,4) = &
            get_pos_index_in_pc(&
            (/pos(1),pos(2),pos(3)-poszl/),&
            positions_pc_l(:,:),&
            transpose(recivec_pc_l)/2/pi)
        end do
        pidxL(countL:countL+atom_orb_num(i)-1) = temp
        countL = countL + atom_orb_num(i)
    end do

    countR = 1
    do i = natoms-buffer_right-2*period_right+1,&
           natoms-buffer_right-period_right
        pos = positions(i,1:3)
        do j = 1,natm_r_uc
            p = positions(R_list(j),1:3)
            do k = 1,n_bloch_x
                do l = 1,n_bloch_y
                    dp(:) = pos(:)-p(:)&
                    -matmul((/dble(k-1),dble(l-1),0.0d0/),latvec_uc(:,:))
                    if (sqrt(dp(1)**2+dp(2)**2+dp(3)**2).lt.1.0d-5) then
                        temp = j 
                    end if
                end do
            end do
        end do
        do j = 1,atom_orb_num(i)
            pos_r(countR+j-1,1:3) = pos(:)
            pos_r(countR+j-1,4) = &
            get_pos_index_in_pc(&
            (/pos(1),pos(2),pos(3)-poszr/),&
            positions_pc_r(:,:),&
            transpose(recivec_pc_r)/2/pi)
        end do
        pidxR(countR:countR+atom_orb_num(i)-1) = temp
        countR = countR + atom_orb_num(i)
    end do  

    end subroutine transverse_uc_setup

    subroutine vel_mat_mpi(ie,ik_core,E_k,Ham,norb,gl,gr,glp,grm,gml,gmr,&
                  G_D_RL,G_D_LR,G_L,G_R,&
                  transl_s,transr_s,transl_ns,transr_ns,&
                  refll_s,reflr_s,refll_ns,reflr_ns,&
                  norb_pc_l,norb_pc_r,&
                  norb_uc_l,norb_uc_r)

    use config
    use surface, only : nl, nr, n_buffer_l, n_buffer_r

    implicit none

    include "mpif.h"

    real(kind=8),intent(in) :: E_k(ne,nk)
    integer(kind=4),intent(in) :: ik_core, ie
    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(in) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(in) :: gl(numprocs,nl,nl),gr(numprocs,nr,nr)
    complex(kind=8),intent(in) :: glp(numprocs,nl,nl),grm(numprocs,nr,nr)
    complex(kind=8),intent(in) :: gml(numprocs,nl,nl),gmr(numprocs,nr,nr)
    complex(kind=8),intent(in) :: G_D_RL(numprocs,nr,nl),G_D_LR(numprocs,nl,nr)
    complex(kind=8),intent(in) :: G_L(numprocs,nl,nl),G_R(numprocs,nr,nr)
    integer(kind=4),intent(in) :: norb_pc_l,norb_pc_r
    integer(kind=4),intent(in) :: norb_uc_l,norb_uc_r
    real(kind=8),intent(out) :: transl_s(ne,nk),transr_s(ne,nk)
    real(kind=8),intent(out) :: transl_ns(ne,nk),transr_ns(ne,nk)
    real(kind=8),intent(out) :: refll_s(ne,nk),reflr_s(ne,nk)
    real(kind=8),intent(out) :: refll_ns(ne,nk),reflr_ns(ne,nk)
    integer(kind=4) :: mm,nn

    mm=myid+1
    nn = k_start_core(mm)+ik_core-1
    call vel_mat(mm,nn,E_k(ie,nn),Ham(mm,:,:),norb,&
        gl(mm,:,:),gr(mm,:,:),glp(mm,:,:),grm(mm,:,:),&
        gml(mm,:,:),gmr(mm,:,:),G_D_RL(mm,:,:),G_D_LR(mm,:,:),&
        G_L(mm,:,:), G_R(mm,:,:), &
        k_grid(nn,:),&
        transl_s(ie,nn),transr_s(ie,nn),&
        transl_ns(ie,nn),transr_ns(ie,nn),&
        refll_s(ie,nn),reflr_s(ie,nn),&
        refll_ns(ie,nn),reflr_ns(ie,nn),&
        norb_pc_l,norb_pc_r,norb_uc_l,norb_uc_r) 

    end subroutine vel_mat_mpi

    subroutine vel_mat(mm,nn,E,Ham,norb,gl,gr,glp,grm,&
                       gml,gmr,G_D_RL,G_D_LR,&
                       G_L,G_R,&
                       ksc,&
                       transl_s,transr_s,&
                       transl_ns,transr_ns,&
                       reflectl_s,reflectr_s,&
                       reflectl_ns,reflectr_ns,&
                       norb_pc_l,norb_pc_r,&
                       norb_uc_l,norb_uc_r)

    use config
    use input, only : latvec, positions, recivec,&
                      orb_uc_l_idx,orb_uc_r_idx,&
                      latvec_pc_l,latvec_pc_r,&
                      latvec_uc_l,latvec_uc_r,&
                      reci_uc_l,reci_uc_r,&
                      recivec_pc_l,recivec_pc_r,&
                      positions_pc_l,positions_pc_r
    use util
    use surface, only : nl, nr, n_buffer_l, n_buffer_r, &
                        a_z_L, a_z_R
    use write_data
    use hamiltonian
    use circular


    implicit none

    integer(kind=4),intent(in) :: mm,nn
    real(kind=8),intent(in) :: E
    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(in) :: Ham(norb,norb)
    complex(kind=8),intent(in) :: gl(nl,nl),gr(nr,nr) ! g_adv_L_minus, g_ret_R_plus
    complex(kind=8),intent(in) :: glp(nl,nl),grm(nr,nr) ! g_ret_L_plus, g_ret_R_minus
    complex(kind=8),intent(in) :: gml(nl,nl),gmr(nr,nr)
    complex(kind=8),intent(in) :: G_D_RL(nr,nl),G_D_LR(nl,nr)
    complex(kind=8),intent(in) :: G_L(nl,nl),G_R(nr,nr)
    integer(kind=4),intent(in) :: norb_pc_l,norb_pc_r
    integer(kind=4),intent(in) :: norb_uc_l,norb_uc_r
    real(kind=8),intent(in)    :: ksc(3)
    real(kind=8)   :: vel(3)
    real(kind=8)   :: ksc_cart(3)
    real(kind=8)   :: transl_s,transr_s
    real(kind=8)   :: transl_ns,transr_ns
    real(kind=8) :: reflectr_s, reflectr_ns
    real(kind=8) :: reflectl_s, reflectl_ns

    complex(kind=8) :: F_adv_L_m(nl,nl),lambda_adv_L_m(nl),&
            U_adv_L_m(nl,nl),V_adv_L_m(nl,nl),inv_U_adv_L_m(nl,nl),&
            V_adv_L_m_half(nl,nl)
    complex(kind=8) :: F_ret_R_p(nr,nr),lambda_ret_R_p(nr),&
            U_ret_R_p(nr,nr),V_ret_R_p(nr,nr),inv_U_ret_R_p(nr,nr),&
            V_ret_R_p_half(nr,nr)
    complex(kind=8) :: F_ret_L_m(nl,nl),lambda_ret_L_m(nl),&
            U_ret_L_m(nl,nl),V_ret_L_m(nl,nl),inv_U_ret_L_m(nl,nl),&
            V_ret_L_m_half(nl,nl)
    complex(kind=8) :: F_adv_R_p(nr,nr),lambda_adv_R_p(nr),&
            U_adv_R_p(nr,nr),V_adv_R_p(nr,nr),inv_U_adv_R_p(nr,nr),&
            V_adv_R_p_half(nr,nr)
    complex(kind=8) :: UL(nl,nl),UR(nr,nr)

    complex(kind=8) :: Q_L(nl,nl),Q_R(nr,nr)

    complex(kind=8) :: trl(nr,nl)
    complex(kind=8) :: tlr(nl,nr)
    complex(kind=8) :: rrr(nr,nr)
    complex(kind=8) :: rll(nl,nl)
    integer(kind=4) :: num_prop_l_adv, num_prop_r_ret
    integer(kind=4) :: num_prop_l_ret, num_prop_r_adv
    complex(kind=8),allocatable :: S(:,:)

    real(kind=8) :: t_mode_trace_L, t_mode_trace_R
    real(kind=8) :: t_mode_L(nl), t_mode_R(nr)
    real(kind=8) :: k_z_L_adv(nl), k_z_R_ret(nr)
    real(kind=8) :: k_z_L_ret(nl), k_z_R_adv(nr)

    real(kind=8) :: kpoint(3), pos_uc(3), pos(3)

    integer(kind=4) :: i, j, n, temp_idx, idx, temp_num_orb,&
                       it, ir, countL, countR, count1, count2
    complex(kind=8) :: temp_sum

    real(kind=8),allocatable :: weightL_adv(:,:),weightR_ret(:,:)
    integer(kind=4) :: bz_idx_L_adv(nl),bz_idx_R_ret(nr)
    real(kind=8),allocatable :: weightL_ret(:,:),weightR_adv(:,:)
    integer(kind=4) :: bz_idx_L_ret(nl),bz_idx_R_adv(nr)
    real(kind=8) :: k_xyz_L_ret(nl,3),k_xyz_L_adv(nl,3)
    real(kind=8) :: k_xyz_R_ret(nr,3),k_xyz_R_adv(nr,3)
    real(kind=8) :: k_xyz_L_ret_pc(nl,3),k_xyz_L_adv_pc(nl,3)
    real(kind=8) :: k_xyz_R_ret_pc(nr,3),k_xyz_R_adv_pc(nr,3)
    real(kind=8) :: vel_xyz_L_ret_pc(nl,3),vel_xyz_L_adv_pc(nl,3)
    real(kind=8) :: vel_xyz_R_ret_pc(nr,3),vel_xyz_R_adv_pc(nr,3)
    integer(kind=4),allocatable :: deg_list(:),deg_all(:),temp(:)
    integer(kind=4) :: start, next_start, case_no, ndeg
    complex(kind=8) :: temp_HL(nl,nl),temp_HR(nr,nr),&
                       temp_evL(nl,nl), temp_evR(nr,nr),&
                       temp_evL1(nl,nl), temp_evR1(nr,nr),&
                       tf_l(nl,nl),tf_r(nr,nr),tf1(nr,nr)
    real(kind=8) :: temp_eL(nl),temp_eR(nr)
    real(kind=8) :: temp_eL1(nl),temp_eR1(nr)
    complex(kind=8),allocatable :: temp_e(:), temp_ev(:,:)
    integer(kind=4) :: l_to_propl(nl), r_to_propr(nr) 
    integer(kind=4) :: l_adv_to_propl(nl), r_ret_to_propr(nr) 
    integer(kind=4) :: l_ret_to_propl(nl), r_adv_to_propr(nr) 
    integer(kind=4),allocatable :: propl_to_l(:), propr_to_r(:) 
    integer(kind=4),allocatable :: ptl_temp(:), ptr_temp(:) 

    complex(kind=8),allocatable :: A(:,:),B(:,:), X(:,:)
    
    real(kind=8) :: ktemp(3),vtemp(3)
    real(kind=8) :: deltakz,deltavz
    real(kind=8) :: kpc_unfold_l(125,3),kpc_unfold_r(125,3)
    real(kind=8) :: vel_unfold(125,3)
    integer(kind=4) :: kpc_deg_num_l,kpc_deg_num_r
    complex(kind=8) :: test(nl,nl),errors 
    integer(kind=4) :: tags = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Left advanced propagating from left to right
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    F_adv_L_m = matmul(gl,&
                 Ham(n_buffer_l+1:n_buffer_l+nl,&
                 n_buffer_l+nl+1:n_buffer_l+2*nl))
    
    call eigen(F_adv_L_m,lambda_adv_L_m,U_adv_L_m)

    ! the momentum in transport direction for left lead
    k_z_L_adv = 100.0
    num_prop_l_adv = 0
    UL = 0.0d0
    do i = 1, nl
        if (abs(abs(lambda_adv_L_m(i))-1.0) < 1.0d-3) then 
            num_prop_l_adv = num_prop_l_adv + 1
            l_adv_to_propl(num_prop_l_adv) = i
            UL(:,num_prop_l_adv) = U_adv_L_m(:,i)
            k_z_L_adv(i) = -log_fbz(lambda_adv_L_m(i),a_z_L)
            !write(*,*) k_z_L_adv(i)
            !write(*,*) i,k_z_L_adv(i)/(2*pi)*7.8177197260163
        end if
    end do 
    if (path_mode .ne. tags) then
        UL(:,1:num_prop_l_adv) = gram_schmidt(UL(:,1:num_prop_l_adv))
        
        do i = 1,num_prop_l_adv
            U_adv_L_m(:,l_adv_to_propl(i)) = UL(:,i)
        end do

    else
        UL(:,1:num_prop_l_adv) = gram_schmidt(UL(:,1:num_prop_l_adv))
        
        do i = 1,num_prop_l_adv
            U_adv_L_m(:,l_adv_to_propl(i)) = UL(:,i)
        end do

    ! There exist degenerate states and we need to apply
    ! rotation for each degenerate eigenvector to make them
    ! orthonormal to each other.
    
    tf_l = 0.0d0
    do i = 1,nl
        tf_l(i,i) = 1.0d0
    end do
    start = 1
    case_no = 1
    n = 0
    allocate(deg_all(1),temp(1))
    deg_all = 0
    do while (case_no .eq. 1)
        call find_degenerate(lambda_adv_L_m,start,deg_list,&
                        next_start,case_no,deg_all)
        ndeg = size(deg_list,1)
        if (ndeg .gt. 1) then
           ! first record the index of degenerate states
            n = n+ndeg
            if (n.eq.ndeg) then
                deallocate(deg_all,temp)
                allocate(deg_all(ndeg),temp(ndeg))
                deg_all(1:ndeg) = deg_list
                temp = deg_all 
            else
                deallocate(deg_all)
                allocate(deg_all(n))
                deg_all(1:n-ndeg) = temp
                deg_all(n-ndeg+1:n) = deg_list
                deallocate(temp)
                allocate(temp(n))
                temp = deg_all
            end if
            temp_HL = Ham(n_buffer_l+1:n_buffer_l+nl,&
                  n_buffer_l+1:n_buffer_l+nl)+ & 
                  Ham(n_buffer_l+1:n_buffer_l+nl,&
                  n_buffer_l+nl+1:n_buffer_l+2*nl)/&
                  lambda_adv_L_m(deg_list(1))+&
                  Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
                  n_buffer_l+1:n_buffer_l+nl)*&
                  lambda_adv_L_m(deg_list(1))
            
            call eigenH(temp_HL,temp_eL,temp_evL)

            allocate(A(nl,ndeg),B(nl,ndeg),X(ndeg,ndeg))        
            X = 0.0d0
            ! compute A
            do i = 1,ndeg
                A(:,i) = U_adv_L_m(:,deg_list(i))
            end do

    !        ! compute B
            countL = 1
            do i = 1,nl
                if (abs(E-temp_eL(i)).lt.1.0d-4) then
                    B(:,countL) = temp_evL(:,i)/&
                     lambda_adv_L_m(deg_list(1))
                    countL = countL +1
                end if
            end do

            ! compute X
            X = lstsq(A,B)
            ! collect all coefficients for degenerate states
            do i = 1,ndeg
                do j = 1,ndeg
                    tf_l(deg_list(i),deg_list(j))=X(i,j)
                end do
            end do
            if ( maxval(abs(X)) > 2) then
                write(*,*) 'X',abs(X)
                write(*,*) 'X'
!                stop
            end if
            errors = maxval(abs(matmul(A,X)-B))
!            write(*,*) 'errors',errors
            deallocate(A,B,X)
        end if
        start = next_start
    end do

    U_adv_L_m = matmul(U_adv_L_m,tf_l)
    do i = 1,num_prop_l_adv
        U_adv_L_m(:,l_adv_to_propl(i)) = U_adv_L_m(:,l_adv_to_propl(i))&
         /sqrt(dot_product(U_adv_L_m(:,l_adv_to_propl(i)),U_adv_L_m(:,l_adv_to_propl(i))))
    end do
    end if
     
    V_adv_L_m = -matmul(matmul(transpose(dconjg(U_adv_L_m)),&
                gml),U_adv_L_m)

    ! we need to correct V_adv_L_m as well
    tf_l = 0.0d0
    do i = 1,nl
        tf_l(i,i) = 1.0d0
    end do
    ndeg = size(deg_all,1)
    !if (deg_all(1).gt.0) then
    !    allocate(B(ndeg,ndeg),temp_e(ndeg),temp_ev(ndeg,ndeg))
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            B(i,j) = V_adv_L_m(deg_all(i),deg_all(j))
    !        end do
    !    end do
    !    call eigen(B,temp_e,temp_ev)
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            tf_l(deg_all(i),deg_all(j)) = temp_ev(i,j)
    !        end do
    !    end do
    !    deallocate(B,temp_e,temp_ev)
    !end if

    V_adv_L_m = matmul(matmul(transpose(dconjg(tf_l)),&
                V_adv_L_m),tf_l)

    inv_U_adv_L_m = inv(transpose(dconjg(U_adv_L_m)))
    !U_adv_L_m =circular_inv(transpose(dconjg(U_adv_L_m)),nl,n_bloch_y,n_bloch_x,reci_uc_l(1,:),-latvec_uc_l(2,:),&
    !reci_uc_l(2,:),-latvec_uc_l(1,:))
    !write(*,*) maxval(abs(inv_U_adv_L_m-U_adv_L_m))
    !stop
           
    V_adv_L_m_half = sqrtm(V_adv_L_m)

    ! exeute the unfolding process
    allocate(weightL_adv(nl,n_bloch_x*n_bloch_y))
    allocate(weightR_ret(nr,n_bloch_x*n_bloch_y))
    bz_idx_L_adv = 0
    bz_idx_R_ret = 0
    k_xyz_L_adv = 100.0
    k_xyz_R_ret = 100.0
    k_xyz_L_adv_pc = 100.0
    k_xyz_R_ret_pc = 100.0
    vel_xyz_L_adv_pc = 0.0d0
    vel_xyz_R_ret_pc = 0.0d0
    weightL_adv = 0.0d0
    weightR_ret = 0.0d0

    do i = 1,nl ! i th column (i th eigenvector)
        if (k_z_L_adv(i).lt.100.0) then
            countL = 1
            do j = 1, n_bloch_x*n_bloch_y
                kpoint(:) = move_in_ws2(matmul(k_shift(j,:) +&
                ksc(:),recivec),reci_uc_l)
                n = 1
                temp_idx = 1
                !write(*,*) 'kz'
                !write(*,*) abs(lambda_adv_L_m(i)),k_z_L_adv(i),recivec(3,3)

                do while (n.le.natm_l_uc)
                    idx = pidxL(temp_idx)
                    temp_num_orb = L_orb_num_list(idx)
                    pos_uc = positions(L_list(idx),:)
                    do it = 1,temp_num_orb
                        temp_sum = 0.0d0
                        do ir = 1,n_bloch_x*n_bloch_y
                            pos = pos_l(temp_idx+it-1+& 
                                       (ir-1)*nl/n_bloch_x/n_bloch_y,1:3)
                            temp_sum = temp_sum + &
                          exp(-i_imag*(-dot_product(matmul(ksc,recivec),pos)+&
                          dot_product(kpoint,pos-pos_uc)))*&
                          U_adv_L_m(temp_idx+it-1+&
                                       (ir-1)*nl/n_bloch_x/n_bloch_y,i)*&
                          lambda_adv_L_m(i)
                        end do
                        weightL_adv(i,countL) = weightL_adv(i,countL) +&
                              (abs(temp_sum))**2/dble(n_bloch_x*n_bloch_y)
                    end do
                    n = n+1 
                    temp_idx = temp_idx+ L_orb_num_list(idx)
                end do
                countL = countL+1
            end do

            !write(*,*) weightL_adv(i,:)
            bz_idx_L_adv(i) = maxloc(weightL_adv(i,:),1)
            k_xyz_L_adv(i,:) = move_in_ws2(matmul(ksc(:) + &
                               k_shift(maxloc(weightL_adv(i,:),1),:),&
                               recivec),reci_uc_l)
            k_xyz_L_adv(i,3) = k_z_L_adv(i)
            ksc_cart = matmul(ksc,recivec)
            ksc_cart(3) = k_z_L_adv(i)

            !write(*,*) weightL_adv(i,:)
            !write(*,*) k_shift(maxloc(weightL_adv(i,:),1),:)
            !stop



            !write(*,*) lambda_adv_L_m(i)
            !write(*,*) dot_product(u_adv_L_m(:,i),u_adv_L_m(:,i))
            !write(*,*) maxval(weightL_adv(i,:))
            !temp_HL = Ham(n_buffer_l+1:n_buffer_l+nl,&
            !      n_buffer_l+1:n_buffer_l+nl)+ &
            !      Ham(n_buffer_l+1:n_buffer_l+nl,&
            !      n_buffer_l+nl+1:n_buffer_l+2*nl)/&
            !      lambda_adv_L_m(i)+&
            !      Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
            !      n_buffer_l+1:n_buffer_l+nl)*&
            !      lambda_adv_L_m(i) 

            !call eigenH(temp_HL,temp_eL1,temp_evL1)
            !write(*,*) matmul(temp_HL,U_adv_L_m(:,i))-E*U_adv_L_m(:,i)
            !stop
            !do j = 1, size(temp_eL1,1)
            !write(*,*) temp_eL1(j)
            !if (abs(temp_eL1(j)-E).lt.0.001) then
            !    write(*,*)j
               !write(*,*)  maxval(abs(temp_evL1(:,j)-U_adv_L_m(:,i)*lambda_adv_L_m(i)))
            !end if
            !end do 
            !stop

            vel_unfold = 0.0d0
            call find_k_pc_l(E,norb_pc_l,&
                  ksc_cart,&
                  k_xyz_L_adv(i,:),&
                  norb_uc_l,nl,&
                  U_adv_L_m(:,i),&
                  pos_l,matmul(ksc,recivec),&
                  kpc_unfold_l,kpc_deg_num_l,vel_unfold)
!            write(*,*) 'degnum_adv_l',kpc_deg_num_l,matmul(k_xyz_L_adv(i,:),inv_real(reci_uc_l)),ksc_cart,matmul(k_xyz_L_adv(i,:),inv_real(recivec_pc_l)) 
            k_xyz_L_adv_pc(i,:) = find_min_dist(kpc_unfold_l(1:kpc_deg_num_l,:),kpc_deg_num_l)
            vel_xyz_L_adv_pc(i,:) = vel_unfold(1,:)
        end if
    end do   


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Left retarded propagating from right to left
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    F_ret_L_m = matmul(transpose(conjg(gl)),&
                   Ham(n_buffer_l+1:n_buffer_l+nl,&
                   n_buffer_l+nl+1:n_buffer_l+2*nl))


    call eigen(F_ret_L_m,lambda_ret_L_m,U_ret_L_m)

    ! the momentum in transport direction for left lead
    k_z_L_ret = 100.0
    num_prop_l_ret = 0 ! count the number of propagating modes
    l_ret_to_propl = 0
    UL = 0.0d0 ! keep the propagating modes
    allocate(ptl_temp(1),propl_to_l(1))
    j = 0
    do i = 1, nl
        if (abs(abs(lambda_ret_L_m(i))-1.0) < 1.0d-3) then 
            k_z_L_ret(i) = -log_fbz(lambda_ret_L_m(i),a_z_L)
            num_prop_l_ret = num_prop_l_ret + 1
            l_ret_to_propl(num_prop_l_ret) = i
            UL(:,num_prop_l_ret) = U_ret_L_m(:,i)

            j = j + 1
            if (j.eq.1) then
                propl_to_l(j) = i
                ptl_temp(j) = i
            else
                deallocate(propl_to_l)
                allocate(propl_to_l(j))
                propl_to_l(1:j-1) = ptl_temp
                propl_to_l(j) = i
                deallocate(ptl_temp)
                allocate(ptl_temp(j))
                ptl_temp = propl_to_l
            end if
        end if
    end do
    

    if (path_mode .ne. tags) then
        UL(:,1:num_prop_l_ret) = gram_schmidt(UL(:,1:num_prop_l_ret))
        do i = 1,num_prop_l_ret
            U_ret_L_m(:,l_ret_to_propl(i)) = UL(:,i)
        end do
    else

    ! There exist degenerate states and we need to apply
    ! rotation for each degenerate eigenvector to make them
    ! orthonormal to each other.
        UL(:,1:num_prop_l_ret) = gram_schmidt(UL(:,1:num_prop_l_ret))
        do i = 1,num_prop_l_ret
            U_ret_L_m(:,l_ret_to_propl(i)) = UL(:,i)
        end do

    tf_l = 0.0d0
    do i = 1,nl
        tf_l(i,i) = 1.0d0
    end do
    start = 1
    case_no = 1
    n = 0
    if (allocated(deg_all)) deallocate(deg_all)
    if (allocated(temp)) deallocate(temp)
    allocate(deg_all(1),temp(1))

    deg_all = 0
    do while (case_no .eq. 1)
        call find_degenerate(lambda_ret_L_m,start,deg_list,&
                        next_start,case_no,deg_all)
        ndeg = size(deg_list,1)

        if (ndeg .gt. 1) then
            ! first record the index of degenerate states
            n = n+ndeg
            if (n.eq.ndeg) then
                deallocate(deg_all,temp)
                allocate(deg_all(ndeg),temp(ndeg))
                deg_all(1:ndeg) = deg_list
                temp = deg_all 
            else
                deallocate(deg_all)
                allocate(deg_all(n))
                deg_all(1:n-ndeg) = temp
                deg_all(n-ndeg+1:n) = deg_list
                deallocate(temp)
                allocate(temp(n))
                temp = deg_all
            end if
            temp_HL = Ham(n_buffer_l+1:n_buffer_l+nl,&
                  n_buffer_l+1:n_buffer_l+nl)+ & 
                  Ham(n_buffer_l+1:n_buffer_l+nl,&
                  n_buffer_l+nl+1:n_buffer_l+2*nl)/&
                  lambda_ret_L_m(deg_list(1))+&
                  Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
                  n_buffer_l+1:n_buffer_l+nl)*&
                  lambda_ret_L_m(deg_list(1))
            
            call eigenH(temp_HL,temp_eL1,temp_evL1)

            allocate(A(nl,ndeg),B(nl,ndeg),X(ndeg,ndeg))        
            X = 0.0d0
            ! compute A
            do i = 1,ndeg
                A(:,i) = U_ret_L_m(:,deg_list(i))
            end do
            ! compute B
            countL = 1
            do i = 1,nl
                if (abs(E-temp_eL1(i)).lt.1.0d-4) then
                    B(:,countL) = temp_evL1(:,i)/&
                     lambda_ret_L_m(deg_list(1))
                    countL = countL +1
                end if
            end do
            ! compute X
            X = lstsq(A,B)
            ! collect all coefficients for degenerate states
            do i = 1,ndeg
                do j = 1,ndeg
                    tf_l(deg_list(i),deg_list(j))=X(i,j)
                end do
            end do
            deallocate(A,B,X)
        end if
        start = next_start
    end do

    U_ret_L_m = matmul(U_ret_L_m,tf_l)

    do i = 1,num_prop_l_ret
        U_ret_L_m(:,l_ret_to_propl(i)) = U_ret_L_m(:,l_ret_to_propl(i))&
         /sqrt(dot_product(U_ret_L_m(:,l_ret_to_propl(i)),U_ret_L_m(:,l_ret_to_propl(i))))
    end do

    end if
     
    V_ret_L_m = -matmul(matmul(transpose(dconjg(U_ret_L_m)),&
                transpose(dconjg(gml))),U_ret_L_m)
  
    ! we need to correct V_adv_L_m as well
    !tf_l = 0.0d0
    !do i = 1,nl
    !    tf_l(i,i) = 1.0d0
    !end do
    ndeg = size(deg_all,1)
    !if (deg_all(1).gt.0) then
    !    allocate(B(ndeg,ndeg),temp_e(ndeg),temp_ev(ndeg,ndeg))
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            B(i,j) = V_adv_L_m(deg_all(i),deg_all(j))
    !        end do
    !    end do
    !    call eigen(B,temp_e,temp_ev)
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            tf_l(deg_all(i),deg_all(j)) = temp_ev(i,j)
    !        end do
    !    end do
    !    deallocate(B,temp_e,temp_ev)
    !end if

    !V_ret_L_m = matmul(matmul(transpose(dconjg(tf_l)),&
    !            V_ret_L_m),tf_l)

    inv_U_ret_L_m = inv(U_ret_L_m)
           
    V_ret_L_m_half = sqrtm(V_ret_L_m) 

    ! exeute the unfolding process for reflected processes
    allocate(weightL_ret(nl,n_bloch_x*n_bloch_y))
    allocate(weightR_adv(nr,n_bloch_x*n_bloch_y))
    bz_idx_L_ret = 0
    bz_idx_R_adv = 0
    k_xyz_L_ret = 100.0
    k_xyz_R_adv = 100.0
    k_xyz_L_ret_pc = 100.0
    k_xyz_R_adv_pc = 100.0
    vel_xyz_L_ret_pc = 0.0d0
    vel_xyz_R_adv_pc = 0.0d0
    weightL_ret = 0.0d0
    weightR_adv = 0.0d0

    do i = 1,nl ! i th column (i th eigenvector)
        if (k_z_L_ret(i).lt.100.0) then
            countL = 1
            do j = 1, n_bloch_x*n_bloch_y
                kpoint(:) = move_in_ws2(matmul(k_shift(j,:) +&
                ksc(:),recivec),reci_uc_l)
                n = 1
                temp_idx = 1
                do while (n.le.natm_l_uc)
                    idx = pidxL(temp_idx)
                    temp_num_orb = L_orb_num_list(idx)
                    pos_uc = positions(L_list(idx),:)
                    do it = 1,temp_num_orb
                        temp_sum = 0.0d0
                        do ir = 1,n_bloch_x*n_bloch_y
                            pos = pos_l(temp_idx+it-1+(ir-1)&
                                       *nl/n_bloch_x/n_bloch_y,1:3)
                            temp_sum = temp_sum + &
                          exp(-i_imag*(-dot_product(matmul(ksc,recivec),pos)+&
                          dot_product(kpoint,pos-pos_uc)))*&
                          U_ret_L_m(temp_idx+it-1+(ir-1)&
                          *nl/n_bloch_x/n_bloch_y,i)*&
                          lambda_ret_L_m(i)
                        end do
                        weightL_ret(i,countL) = weightL_ret(i,countL) +&
                              (abs(temp_sum))**2/dble(n_bloch_x*n_bloch_y)
                    end do
                    n = n+1 
                    temp_idx = temp_idx + L_orb_num_list(idx)
                end do
                countL = countL+1
            end do
            !write(*,*) weightL_ret(i,:)

            bz_idx_L_ret(i) = maxloc(weightL_ret(i,:),1)
            k_xyz_L_ret(i,:) = move_in_ws2(matmul(ksc(:) + &
                               k_shift(maxloc(weightL_ret(i,:),1),:),&
                               recivec),reci_uc_l) 
            k_xyz_L_ret(i,3) = k_z_L_ret(i)

            ksc_cart = matmul(ksc,recivec)
            ksc_cart(3) = k_z_L_ret(i)

            vel_unfold = 0.0d0
            call find_k_pc_l(E,norb_pc_l,&
                  ksc_cart,&
                  k_xyz_L_ret(i,:),&
                  norb_uc_l,nl,&
                  U_ret_L_m(:,i),&
                  pos_l,matmul(ksc,recivec),&
                  kpc_unfold_l,kpc_deg_num_l,vel_unfold)

            k_xyz_L_ret_pc(i,:)=find_min_dist(kpc_unfold_l(1:kpc_deg_num_l,:),kpc_deg_num_l) 
            vel_xyz_L_ret_pc(i,:) = vel_unfold(1,:)
        end if 
    end do

    ! align positive kz and negative kz modes' indices
    !do i = 1,num_prop_l_ret
    !    deltakz = k_xyz_L_ret(l_ret_to_propl(i),3)&
    !            + k_xyz_L_adv(l_adv_to_propl(i),3)
    !    deltavz = vel_xyz_L_ret_pc(l_ret_to_propl(i),3)&
    !                + vel_xyz_L_adv_pc(l_adv_to_propl(i),3)
    !    if ((deltakz .gt. 1.0d-5).or.(deltavz .gt. 1.0d-5)) then
    !        write(*,*) i,deltakz,deltavz
    !        write(*,*) k_xyz_L_ret(l_ret_to_propl(i),3)/(2*pi/a_z_R),k_xyz_L_adv(l_adv_to_propl(i),3)/(2*pi/a_z_R)
    !    end if
    !end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Right retarded propagating from left to right
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    F_ret_R_p = matmul(gr,Ham(norb-n_buffer_r-nr+1:&
              norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
                norb-n_buffer_r-nr))    
    
    call eigen(F_ret_R_p,lambda_ret_R_p,U_ret_R_p)

    ! the momentum in transport direction for left lead
    k_z_R_ret = 100.0
    num_prop_r_ret = 0 ! count the number of propagating modes
    r_ret_to_propr = 0
    UR = 0.0d0 ! keep the propagating modes
    allocate(ptr_temp(1),propr_to_r(1))
    j = 0

    do i = 1, nr
        if (abs(abs(lambda_ret_R_p(i))-1.0) < 1.0d-3) then 
            k_z_R_ret(i) = log_fbz(lambda_ret_R_p(i),a_z_R)
            !write(*,*) k_z_R_ret(i)
            num_prop_r_ret    = num_prop_r_ret + 1
            r_ret_to_propr(num_prop_r_ret) = i
            UR(:,num_prop_r_ret) = U_ret_R_p(:,i)

            j = j + 1
            if (j.eq.1) then
                propr_to_r(j) = i
                ptr_temp(j) = i
            else
                deallocate(propr_to_r)
                allocate(propr_to_r(j))
                propr_to_r(1:j-1) = ptr_temp
                propr_to_r(j) = i
                deallocate(ptr_temp)
                allocate(ptr_temp(j))
                ptr_temp = propr_to_r
            end if
        end if
    end do
    
    if (path_mode .ne. tags) then
        UR(:,1:num_prop_r_ret) = gram_schmidt(UR(:,1:num_prop_r_ret))
        do i = 1,num_prop_r_ret
            U_ret_R_p(:,r_ret_to_propr(i)) = UR(:,i)
        end do
    else
    ! There exist degenerate states and we need to apply
    ! rotation for each degenerate eigenvector to make them
    ! orthonormal to each other.
         UR(:,1:num_prop_r_ret) = gram_schmidt(UR(:,1:num_prop_r_ret))
        do i = 1,num_prop_r_ret
            U_ret_R_p(:,r_ret_to_propr(i)) = UR(:,i)
        end do
   
    tf_r = 0.0d0
    do i = 1,nr
        tf_r(i,i) = 1.0d0
    end do
    start = 1
    case_no = 1
    n = 0
    if (allocated(deg_all)) deallocate(deg_all)
    if (allocated(temp)) deallocate(temp)
    allocate(deg_all(1),temp(1))
    deg_all(1) = 0
    do while (case_no .eq. 1)
        call find_degenerate(lambda_ret_R_p,start,deg_list,&
                        next_start,case_no,deg_all)
        ndeg = size(deg_list,1)
        if (ndeg .gt. 1) then
            ! first record the index of degenerate states
            n = n+ndeg
            if (n.eq.ndeg) then
                deallocate(deg_all,temp)
                allocate(deg_all(ndeg),temp(ndeg))
                deg_all(1:ndeg) = deg_list
                temp = deg_all 
            else
                deallocate(deg_all)
                allocate(deg_all(n))
                deg_all(1:n-ndeg) = temp
                deg_all(n-ndeg+1:n) = deg_list
                deallocate(temp)
                allocate(temp(n))
                temp = deg_all
            end if
            if (myid.eq.0)then
            end if
            temp_HR = Ham(norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r,norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r)+&
                  Ham(norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
                  norb-n_buffer_r-nr)/&
                  lambda_ret_R_p(deg_list(1))+&
                  Ham(norb-n_buffer_r-2*nr+1:&
                  norb-n_buffer_r-nr,norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r)*&
                  lambda_ret_R_p(deg_list(1))            
            call eigenH(temp_HR,temp_eR,temp_evR)
            allocate(A(nr,ndeg),B(nr,ndeg),X(ndeg,ndeg))        
            X = 0.0d0
            ! compute A
            do i = 1,ndeg
                A(:,i) = U_ret_R_p(:,deg_list(i))
            end do
            ! compute B
            countR = 1
            do i = 1,nr
                if (abs(E-temp_eR(i)).lt.1.0d-4) then
                    B(:,countR) = temp_evR(:,i)/&
                     lambda_ret_R_p(deg_list(1))
                    countR = countR +1
                end if
            end do
            ! compute X
            X = lstsq(A,B)
            ! collect all coefficients for degenerate states
            do i = 1,ndeg
                do j = 1,ndeg
                    tf_r(deg_list(i),deg_list(j))=X(i,j)
                end do
            end do
            deallocate(A,B,X)
        end if
        start = next_start
    end do
!      if (myid.eq.0) then
!             do i=1,size(deg_all,1)
!             do j=1,size(deg_all,1)
!            if (abs(tf_r(deg_all(i),deg_all(j))).gt.0) then
!            write(*,*) 'retR',deg_all(i),deg_all(j)
!            write(*,*) tf_r(deg_all(i),deg_all(j))
!            end if
!            end do
!            end do
            !stop
!     end if

    tf1 = tf_r
    U_ret_R_p = matmul(U_ret_R_p,tf_r) 
    do i = 1,num_prop_r_ret
        U_ret_R_p(:,r_ret_to_propr(i)) = U_ret_R_p(:,r_ret_to_propr(i))&
         /sqrt(dot_product(U_ret_R_p(:,r_ret_to_propr(i)),U_ret_R_p(:,r_ret_to_propr(i))))
    end do

    end if

    V_ret_R_p = matmul(matmul(transpose(dconjg(U_ret_R_p)),&
                gmr),U_ret_R_p)

    ! we need to correct V_ret_R_p as well
    tf_r = 0.0d0
    do i = 1,nr
        tf_r(i,i) = 1.0d0
    end do
    ndeg = size(deg_all,1)
    !if (deg_all(1).gt.0) then
    !    allocate(B(ndeg,ndeg),temp_e(ndeg),temp_ev(ndeg,ndeg))
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            B(i,j) = V_ret_R_p(deg_all(i),deg_all(j))
    !        end do
    !    end do
    !    call eigen(B,temp_e,temp_ev)
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            tf_r(deg_all(i),deg_all(j)) = temp_ev(i,j)
    !        end do
    !    end do
    !    deallocate(B,temp_e,temp_ev)
    !end if
    V_ret_R_p = matmul(matmul(transpose(dconjg(tf_r)),&
                V_ret_R_p),tf_r)
    
    inv_U_ret_R_p = inv(U_ret_R_p)

    V_ret_R_p_half = sqrtm(V_ret_R_p)

    do i = 1,nr ! i th column (i th eigenvector)
        if (k_z_R_ret(i)<100.0) then
            countR = 1
            do j = 1, n_bloch_x*n_bloch_y
                kpoint(:) = move_in_ws2(matmul(k_shift(j,:) +&
                ksc(:),recivec),reci_uc_r)
                n = 1
                temp_idx = 1
                do while (n.le.natm_r_uc)
                    idx = pidxR(temp_idx)
                    temp_num_orb = R_orb_num_list(idx)
                    pos_uc = positions(R_list(idx),:)
                    do it = 1,temp_num_orb
                        temp_sum = 0.0d0
                        do ir = 1,n_bloch_x*n_bloch_y
                            pos = pos_r(temp_idx+it-1+(ir-1)&
                                       *nr/n_bloch_x/n_bloch_y,1:3)
                            temp_sum = temp_sum + &
                         exp(-i_imag*(-dot_product(matmul(ksc,recivec),pos)+&
                          dot_product(kpoint,pos-pos_uc)))*& 
                          U_ret_R_p(temp_idx+it-1+(ir-1)&
                          *nr/n_bloch_x/n_bloch_y,i)*&
                          lambda_ret_R_p(i)
                        end do
                        weightR_ret(i,countR) = weightR_ret(i,countR) +&
                              abs(temp_sum)**2/dble(n_bloch_x*n_bloch_y)
                    end do
                    n = n+1 
                    temp_idx = temp_idx+ R_orb_num_list(idx)
                end do
                countR = countR+1
            end do
            bz_idx_R_ret(i) = maxloc(weightR_ret(i,:),1)
            k_xyz_R_ret(i,:) = move_in_ws2(matmul(ksc(:) + &
                               k_shift(maxloc(weightR_ret(i,:),1),:),&
                               recivec),reci_uc_r)
            k_xyz_R_ret(i,3) = k_z_R_ret(i)

            ksc_cart = matmul(ksc,recivec)
            ksc_cart(3) = k_z_R_ret(i)
!            write(*,*)weightR_ret(i,:)

            vel_unfold = 0.0d0
            call find_k_pc_r(E,norb_pc_r,&
                  ksc_cart,&
                  k_xyz_R_ret(i,:),&
                  norb_uc_r,nr,&
                  u_ret_R_p(:,i),&
                  pos_r,matmul(ksc,recivec),&
                  kpc_unfold_r,kpc_deg_num_r,vel_unfold)

            k_xyz_R_ret_pc(i,:) = find_min_dist(kpc_unfold_r(1:kpc_deg_num_r,:),kpc_deg_num_r)
            vel_xyz_R_ret_pc(i,:) = vel_unfold(1,:)
        end if
    end do  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Right advanced propagating from right to left
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    F_adv_R_p = matmul(transpose(dconjg(gr)),&
              Ham(norb-n_buffer_r-nr+1:&
              norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
                norb-n_buffer_r-nr))    
    
    call eigen(F_adv_R_p,lambda_adv_R_p,U_adv_R_p)

    ! the momentum in transport direction for left lead
    k_z_R_adv = 100.0
    num_prop_r_adv = 0
    r_adv_to_propr = 0
    UR = 0.0
    do i = 1, nr
        if (abs(abs(lambda_adv_R_p(i))-1.0) < 1.0d-3) then 
            num_prop_r_adv   = num_prop_r_adv + 1
            r_adv_to_propr(num_prop_r_adv) = i
            UR(:,num_prop_r_adv) = U_adv_R_p(:,i)
            k_z_R_adv(i) = log_fbz(lambda_adv_R_p(i),a_z_R)
        end if
    end do

    if (path_mode .ne. tags) then
        UR(:,1:num_prop_r_adv) = gram_schmidt(UR(:,1:num_prop_r_adv))
        do i = 1,num_prop_r_adv
            U_adv_R_p(:,r_adv_to_propr(i)) = UR(:,i)
        end do
    else
        UR(:,1:num_prop_r_adv) = gram_schmidt(UR(:,1:num_prop_r_adv))
        do i = 1,num_prop_r_adv
            U_adv_R_p(:,r_adv_to_propr(i)) = UR(:,i)
        end do

    ! There exist degenerate states and we need to apply
    ! rotation for each degenerate eigenvector to make them
    ! orthonormal to each other.
    tf_r = 0.0d0
    do i = 1,nr
        tf_r(i,i) = 1.0d0
    end do
    start = 1
    case_no = 1
    n = 0
    if (allocated(deg_all)) deallocate(deg_all)
    if (allocated(temp)) deallocate(temp)
    allocate(deg_all(1),temp(1))
    deg_all(1) = 0
    do while (case_no .eq. 1)
        call find_degenerate(lambda_adv_R_p,start,deg_list,&
                        next_start,case_no,deg_all)
        ndeg = size(deg_list,1)
        if (ndeg .gt. 1) then
            ! first record the index of degenerate states
            n = n+ndeg
            if (n.eq.ndeg) then
                deallocate(deg_all,temp)
                allocate(deg_all(ndeg),temp(ndeg))
                deg_all(1:ndeg) = deg_list
                temp = deg_all 
            else
                deallocate(deg_all)
                allocate(deg_all(n))
                deg_all(1:n-ndeg) = temp
                deg_all(n-ndeg+1:n) = deg_list
                deallocate(temp)
                allocate(temp(n))
                temp = deg_all
            end if
            temp_HR = Ham(norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r,norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r)+&
                  Ham(norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
                  norb-n_buffer_r-nr)/&
                  lambda_adv_R_p(deg_list(1))+&
                  Ham(norb-n_buffer_r-2*nr+1:&
                  norb-n_buffer_r-nr,norb-n_buffer_r-nr+1:&
                  norb-n_buffer_r)*&
                  lambda_adv_R_p(deg_list(1))            
            call eigenH(temp_HR,temp_eR1,temp_evR1)
            allocate(A(nr,ndeg),B(nr,ndeg),X(ndeg,ndeg))        
            X = 0.0d0
            ! compute A
            do i = 1,ndeg
                A(:,i) = U_adv_R_p(:,deg_list(i))
            end do
            ! compute B
            countR = 1
            do i = 1,nr
                if (abs(E-temp_eR1(i)).lt.1.0d-4) then
                    B(:,countR) = temp_evR1(:,i)/&
                     lambda_adv_R_p(deg_list(1))
                    countR = countR +1
                end if
            end do
            ! compute X
            X = lstsq(A,B)
            ! collect all coefficients for degenerate states
            do i = 1,ndeg
                do j = 1,ndeg
                    tf_r(deg_list(i),deg_list(j))=X(i,j)
                end do
            end do
            deallocate(A,B,X)
        end if
        start = next_start
    end do

    U_adv_R_p = matmul(U_adv_R_p,tf_r)

    do i = 1,num_prop_r_adv
        U_adv_R_p(:,r_adv_to_propr(i)) = U_adv_R_p(:,r_adv_to_propr(i))&
         /sqrt(dot_product(U_adv_R_p(:,r_adv_to_propr(i)),U_adv_R_p(:,r_adv_to_propr(i))))
    end do

    end if
!    U_adv_R_p = matmul(U_adv_R_p,tf1)

!    if (myid.eq.0)then
!             write(*,*) dot_product(U_adv_R_p(:,125),U_adv_R_p(:,126))
            !stop
! end if

    V_adv_R_p = matmul(matmul(transpose(dconjg(U_adv_R_p)),&
                transpose(dconjg(gmr))),U_adv_R_p)


    ! we need to correct V_adv_R_p as well
    tf_r = 0.0d0
    do i = 1,nr
        tf_r(i,i) = 1.0d0
    end do
    ndeg = size(deg_all,1)
    !if (deg_all(1).gt.0) then
    !    allocate(B(ndeg,ndeg),temp_e(ndeg),temp_ev(ndeg,ndeg))
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            B(i,j) = V_adv_R_p(deg_all(i),deg_all(j))
    !        end do
    !    end do
    !    call eigen(B,temp_e,temp_ev)
    !    do i = 1,ndeg
    !        do j = 1,ndeg
    !            tf_r(deg_all(i),deg_all(j)) = temp_ev(i,j)
    !        end do
    !    end do
    !    deallocate(B,temp_e,temp_ev)
    !end if
    V_adv_R_p = matmul(matmul(transpose(dconjg(tf_r)),&
                V_adv_R_p),tf_r)
    
    inv_U_adv_R_p = inv(transpose(dconjg(U_adv_R_p)))

    V_adv_R_p_half = sqrtm(V_adv_R_p)
   
    do i = 1,nr ! i th column (i th eigenvector)
        if (k_z_R_adv(i)<100.0) then
            countR = 1
            do j = 1, n_bloch_x*n_bloch_y
                kpoint(:) = move_in_ws2(matmul(k_shift(j,:) +&
                ksc(:),recivec),reci_uc_r)
                n = 1
                temp_idx = 1
                do while (n.le.natm_r_uc)
                    idx = pidxR(temp_idx)
                    temp_num_orb = R_orb_num_list(idx)
                    pos_uc = positions(R_list(idx),:)
                    do it = 1,temp_num_orb
                        temp_sum = 0.0d0
                        do ir = 1,n_bloch_x*n_bloch_y
                            pos = pos_r(temp_idx+it-1+(ir-1)&
                                       *nr/n_bloch_x/n_bloch_y,:)
                            temp_sum = temp_sum + &
                          exp(-i_imag*(-dot_product(matmul(ksc,recivec),pos)+&
                          dot_product(kpoint,pos-pos_uc)))*&  
                          U_adv_R_p(temp_idx+it-1+(ir-1)&
                           *nr/n_bloch_x/n_bloch_y,i)*&
                          lambda_adv_R_p(i)
                        end do
                        weightR_adv(i,countR) = weightR_adv(i,countR) +&
                              abs(temp_sum)**2/dble(n_bloch_x*n_bloch_y)
                    end do
                    n = n+1 
                    temp_idx = temp_idx+ R_orb_num_list(idx)
                end do
                countR = countR+1
            end do

            bz_idx_R_adv(i) = maxloc(weightR_adv(i,:),1)
            k_xyz_R_adv(i,:) = move_in_ws2(matmul(ksc(:) + &
                               k_shift(maxloc(weightR_adv(i,:),1),:),&
                               recivec),reci_uc_r)
            k_xyz_R_adv(i,3) = k_z_R_adv(i)

            ksc_cart = matmul(ksc,recivec)
            ksc_cart(3) = k_z_R_adv(i)
            !write(*,*) weightR_adv(i,:)

            vel_unfold = 0.0d0
            call find_k_pc_r(E,norb_pc_r,&
                  ksc_cart,&
                  k_xyz_R_adv(i,:),&
                  norb_uc_r,nr,&
                  u_adv_R_p(:,i),&
                  pos_r,matmul(ksc,recivec),&
                  kpc_unfold_r,kpc_deg_num_r,vel_unfold)

            k_xyz_R_adv_pc(i,:) = find_min_dist(kpc_unfold_r(1:kpc_deg_num_r,:),kpc_deg_num_r)
            vel_xyz_R_adv_pc(i,:) = vel_unfold(1,:)

        end if
    end do  

    ! retarded bulk Green's function
    Q_L = 0.0d0
    Q_L = (E+eta*i_imag)*eyemat(nl) - &
      Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
      n_buffer_l+nl+1:n_buffer_l+2*nl) - &
      matmul(matmul(Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
      n_buffer_l+1:n_buffer_l+nl) ,&
      transpose(dconjg(gl))),&
      Ham(n_buffer_l+1:n_buffer_l+nl,&
      n_buffer_l+nl+1:n_buffer_l+2*nl))-&
      matmul(matmul(Ham(n_buffer_l+1:n_buffer_l+nl,&
      n_buffer_l+nl+1:n_buffer_l+2*nl),&
      glp),&
      Ham(n_buffer_l+nl+1:n_buffer_l+2*nl,&
      n_buffer_l+1:n_buffer_l+nl))
    
    Q_R = 0.0d0
    Q_R = (E+eta*i_imag)*eyemat(nr) - &
      Ham(norb-n_buffer_r-2*nr+1:&
      norb-n_buffer_r-nr,norb-n_buffer_r-2*nr+1:&
      norb-n_buffer_r-nr)-&
      matmul(matmul(Ham(norb-n_buffer_r-2*nr+1:&
      norb-n_buffer_r-nr,norb-n_buffer_r-nr+1:&
      norb-n_buffer_r),gr),&
      Ham(norb-n_buffer_r-nr+1:&
      norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
      norb-n_buffer_r-nr))-&
      matmul(matmul(Ham(norb-n_buffer_r-nr+1:&
      norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
      norb-n_buffer_r-nr),dconjg(transpose(grm))),&
      Ham(norb-n_buffer_r-2*nr+1:&
      norb-n_buffer_r-nr,norb-n_buffer_r-nr+1:&
      norb-n_buffer_r))

    trl = 0.0d0
    tlr = 0.0d0
    rrr = 0.0d0
    rll = 0.0d0

    trl = matmul(inv_U_adv_L_m,V_adv_L_m_half)
    trl = matmul(G_D_RL,trl)
    trl = matmul(inv_U_ret_R_p,trl)
    trl = i_imag*(matmul(V_ret_R_p_half,trl))

    tlr = matmul(inv_U_adv_R_p,V_adv_R_p_half)
    tlr = matmul(G_D_LR,tlr)
    tlr = matmul(inv_U_ret_L_m,tlr)
    tlr = i_imag*matmul(V_ret_L_m_half,tlr)

    rrr = matmul(inv_U_adv_R_p,V_adv_R_p_half)
    rrr = matmul(G_R-inv(Q_R),rrr)
    rrr = matmul(inv_U_ret_R_p,rrr)
    rrr = i_imag*(matmul(V_ret_R_p_half,rrr)) ! ret_R to adv_R
 
    rll = matmul(inv_U_adv_L_m,V_adv_L_m_half)
    rll = matmul(G_L-inv(Q_L),rll)
    rll = matmul(inv_U_ret_L_m,rll)
    rll = i_imag*matmul(V_ret_L_m_half,rll) ! ret_L to adv_L

    allocate(S(num_prop_l_adv+num_prop_r_ret,num_prop_l_adv+num_prop_r_ret))
    S = 0.0d0
    count1 = 1
    if (myid.eq.0) then
    do i = 1, nl
        if (k_z_L_ret(i).lt.100.0) then
            count2 = 1
            do j = 1, nl
                if (k_z_L_adv(j).lt.100.0) then
                   S(count1,count2) = rll(i,j)
                   count2 = count2 + 1
                end if
            end do
            do j = 1, nr
                if (k_z_R_adv(j).lt.100.0) then
                   S(count1,count2) = tlr(i,j)
                   count2 = count2 + 1
                end if
            end do
            count1 = count1 +1
        end if
    end do
    do i = 1, nr
        if (k_z_R_ret(i).lt.100.0) then
            count2 = 1
            do j = 1, nl
                if (k_z_L_adv(j).lt.100.0) then
                   S(count1,count2) = trl(i,j)
                   count2 = count2 + 1
                end if
            end do
            do j = 1, nr
                if (k_z_R_adv(j).lt.100.0) then
                   S(count1,count2) = rrr(i,j)
                   count2 = count2 + 1
                end if
            end do
            count1 = count1 +1
        end if
    end do
    end if
    reflectr_s = 0.0
!    S = matmul(transpose(dconjg(S)),S)
    !write(*,*) 'n'
    !stop
 !    end if
 !      write(*,*) S(j,j)
 !   end do
!     write(*,*) num_prop_l
!     write(*,*)reflectr_s 
!    stop


  ! end if

    !S(1:nl,1:nl) = rll
    !S(1:nl,nl+1:nl+nr) = tlr
    !S(nl+1:nl+nr,1:nl) = trl
    !S(nl+1:nl+nr,nl+1:nl+nr) = rrr

  !  if (myid.eq.0)then
!            write(*,*) abs(S(5,:))
  !          S = matmul(transpose(dconjg(S)),S)
            !do i = 1,num_prop_l+num_prop_r
  !          write(*,*) S(5,:)
            !end do
  !          stop
  !  end if 
    

    t_mode_trace_L = real(trace(matmul(transpose(dconjg(trl)),trl)))

    t_mode_trace_R = real(trace(matmul(trl,transpose(dconjg(trl)))))
    
    t_mode_L = 0.0d0 
    t_mode_R = 0.0d0

    do i = 1, nl
        do j = 1, nr
            t_mode_L(i) = t_mode_L(i) + abs(trl(j,i))**2
            t_mode_R(j) = t_mode_R(j) + abs(trl(j,i))**2
            if ((bz_idx_L_adv(i).eq.bz_idx_R_ret(j)).and.&
                (bz_idx_L_adv(i).gt.0)) then
                transr_s = transr_s + abs(trl(j,i))**2
                !if (abs(trl(j,i))>0.01)then
                !write(*,*) j,i,abs(trl(j,i))**2
                !end if
            else if ((bz_idx_L_adv(i).ne.bz_idx_R_ret(j)).and.&
                (bz_idx_R_ret(j).gt.0).and.&
                (bz_idx_L_adv(i).gt.0)) then
                transr_ns = transr_ns + abs(trl(j,i))**2
            end if
        end do
    end do
    do i = 1, nl
        do j = 1, nr
            if ((bz_idx_L_ret(i).eq.bz_idx_R_adv(j)).and.&
                (bz_idx_L_ret(i).gt.0)) then
                transl_s = transl_s + abs(tlr(i,j))**2
            else if ((bz_idx_L_ret(i).ne.bz_idx_R_adv(j)).and.&
                (bz_idx_R_adv(j).gt.0).and.&
                (bz_idx_L_ret(i).gt.0)) then
                transl_ns = transl_ns + abs(tlr(i,j))**2
            end if
        end do
    end do
 
    reflectr_s = 0.0d0
    reflectr_ns = 0.0d0 

    do i = 1, nr
        do j = 1, nr
            if ((bz_idx_R_ret(i).eq.bz_idx_R_adv(j)).and.&
                (bz_idx_R_ret(i).gt.0)) then
                reflectr_s = reflectr_s + abs(rrr(i,j))**2
            else if ((bz_idx_R_ret(i).ne.bz_idx_R_adv(j)).and.&
                     (bz_idx_R_ret(i).gt.0).and.&
                     (bz_idx_R_adv(j).gt.0)) then
                reflectr_ns = reflectr_ns + abs(rrr(i,j))**2
            end if
        end do
    end do

    reflectl_s = 0.0d0
    reflectl_ns = 0.0d0 

    do i = 1, nl
        do j = 1, nl
            if ((bz_idx_L_ret(i).eq.bz_idx_L_adv(j)).and.&
                (bz_idx_L_ret(i).gt.0)) then
                reflectl_s = reflectl_s + abs(rll(i,j))**2
            else if ((bz_idx_L_ret(i).ne.bz_idx_L_adv(j)).and.&
                     (bz_idx_L_ret(i).gt.0).and.&
                     (bz_idx_L_adv(j).gt.0)) then
                reflectl_ns = reflectl_ns + abs(rll(i,j))**2
            end if
        end do
    end do

        call write_specular(E,nn,rll,k_xyz_L_ret,k_xyz_L_adv,&
                            k_xyz_L_ret_pc,k_xyz_L_adv_pc,&
                            vel_xyz_L_ret_pc,vel_xyz_L_adv_pc,'left')
        call write_specular(E,nn,rrr,k_xyz_R_ret,k_xyz_R_adv,&
                            k_xyz_R_ret_pc,k_xyz_R_adv_pc,&
                            vel_xyz_R_ret_pc,vel_xyz_R_adv_pc,'right')
        call write_specular(E,nn,trl,k_xyz_R_ret,k_xyz_L_adv,&
                            k_xyz_R_ret_pc,k_xyz_L_adv_pc,&
                            vel_xyz_R_ret_pc,vel_xyz_L_adv_pc,'trl')
        call write_specular(E,nn,tlr,k_xyz_L_ret,k_xyz_R_adv,&
                            k_xyz_L_ret_pc,k_xyz_R_adv_pc,&
                            vel_xyz_L_ret_pc,vel_xyz_R_adv_pc,'tlr')
        call write_md(E,nn,trl,rll,k_xyz_R_ret,k_xyz_L_ret,k_xyz_L_adv,&
                            k_xyz_R_ret_pc,k_xyz_L_ret_pc,k_xyz_L_adv_pc,&
                            vel_xyz_R_ret_pc,vel_xyz_L_ret_pc,vel_xyz_L_adv_pc,'mdl')
        call write_md(E,nn,tlr,rrr,k_xyz_L_ret,k_xyz_R_ret,k_xyz_R_adv,&
                            k_xyz_L_ret_pc,k_xyz_R_ret_pc,k_xyz_R_adv_pc,&
                            vel_xyz_L_ret_pc,vel_xyz_R_ret_pc,vel_xyz_R_adv_pc,'mdr') 

    if (verbosity .eq. 1) then
        call write_md2md(E,nn,tlr,k_xyz_L_ret,k_xyz_R_adv,&
                           k_xyz_L_ret_pc,k_xyz_R_adv_pc,&
                           vel_xyz_L_ret_pc,vel_xyz_R_adv_pc,'tlr')
        call write_md2md(E,nn,trl,k_xyz_R_ret,k_xyz_L_adv,&
                           k_xyz_R_ret_pc,k_xyz_L_adv_pc,&
                           vel_xyz_R_ret_pc,vel_xyz_L_adv_pc,'trl')
        call write_md2md(E,nn,rrr,k_xyz_R_ret,k_xyz_R_adv,&
                           k_xyz_R_ret_pc,k_xyz_R_adv_pc,&
                           vel_xyz_R_ret_pc,vel_xyz_R_adv_pc,'rrr')
        call write_md2md(E,nn,rll,k_xyz_L_ret,k_xyz_L_adv,&
                           k_xyz_L_ret_pc,k_xyz_L_adv_pc,&
                           vel_xyz_L_ret_pc,vel_xyz_L_adv_pc,'rll')
    end if

!   if (myid.eq.0) then
!           write(*,*) k_xyz_L_ret(:,3)
!           stop
!       write(*,*) transr_ns,transr_s
!       write(*,*) reflectr_ns,reflectr_s
!       write(*,*) reflectl_ns,reflectl_s
!       rrr = matmul(transpose(dconjg(rrr)),rrr)
!       write(*,*) maxval(abs(tlr))
!       write(*,*) num_prop_l,num_prop_r
!       write(*,*) '**', trace(matmul(transpose(dconjg(trl)),trl)),trace(matmul(transpose(dconjg(tlr)),tlr))

!       stop
!       stop
!  end if    

    end subroutine vel_mat
end module velocity_mat
