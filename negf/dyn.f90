module dyn

implicit none

contains

    subroutine gen_dyn_mpi(ie,ik_core,force_constant,Ham,eigvec,&
                           eig,norb,pos,pos_s,cell,idx_mat,dist_mat,&
                           num_mat)

    use config

    implicit none

    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ie,ik_core
    integer(kind=4),intent(in) :: norb
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    complex(kind=8),intent(out) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(out) :: eigvec(numprocs,norb,norb)
    real(kind=8),intent(out) :: eig(numprocs,norb)
    integer(kind=4) :: idx_mat(:)
    integer(kind=4) :: num_mat(:,:)
    real(kind=8)    :: dist_mat(:,:)

    mm = myid+1
    nn = k_start_core(mm)+ik_core-1
    call gen_dyn(ie,k_grid(nn,:),force_constant,&
         Ham(mm,:,:),eigvec(mm,:,:),eig(mm,:),norb,pos,pos_s,cell,&
         idx_mat(:),dist_mat(:,:),num_mat(:,:))

    end subroutine gen_dyn_mpi

    subroutine gen_dyn(ie,kpoint,force_constant,Ham,vec,eig,norb,pos,pos_s,cell,&
                       idx_mat,dist_mat,num_mat)

    use config
    use input
    use util

    implicit none

    integer(kind=4),intent(in) :: ie
    real(kind=8),intent(in) :: kpoint(3)
    real(kind=8) :: kpoint_cart(3)
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    integer(kind=4) :: norb, norb_uc
    complex(kind=8) :: Ham(norb,norb)
    complex(kind=8),intent(out) :: vec(norb,norb)
    real(kind=8),intent(out) :: eig(norb)
    complex(kind=8) :: eigc(norb)
    integer(kind=4) :: n_uc, n_sc
    integer(kind=4) :: istart
    integer(kind=4) :: i,j,atm_idx,atm_jdx
    integer(kind=4) :: jdx_mod
    integer(kind=4) :: ii,jj,kk,i_uc,idx,nn
    real(kind=8)  :: m1,m2
    integer(kind=4) :: ndiag,ll
    real(kind=8) :: dr(3)
    real(kind=8) :: r_shift(3)
    integer(kind=4) :: ix,iy,iz
    integer(kind=4) :: idx_s, is, jdx_s
    integer(kind=4) :: atm_idx_s,atm_jdx_s
    real(kind=8) :: pos_temp(3),r_temp(3)
    integer(kind=4) :: dx,dy,dz
    integer(kind=4) :: idx_mat(:)
    real(kind=8) ::    dist_mat(:,:)
    integer(kind=4) :: num_mat(:,:)
    integer(kind=4) :: icart,jcart
    integer(kind=4) :: temp_count,temp_idx,temp_i
    real(kind=8) :: tempsum,tempsum1,tempsum2

    n_sc = size(force_constant,2)/3
    norb_uc = norb/n_bloch_x/n_bloch_y
    n_uc = norb_uc / 3
    
    istart = n_uc*(n_fc_uc_x*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)&
            +n_fc_uc_y*(2*n_fc_uc_z+1)+n_fc_uc_z)+1

    Ham = 0.0d0 
    eig = 0.0d0
    eigc = 0.0d0
    vec = 0.0d0

    kpoint_cart = matmul(kpoint,recivec) ! new to convert

    temp_count = 1
    do atm_idx_s = 1,norb/3
        iz = (atm_idx_s-1)/(n_bloch_y*n_bloch_x*n_uc/nz)
        iy = (atm_idx_s - iz*(n_bloch_y*n_bloch_x*n_uc/nz) - 1)/(n_bloch_x*n_uc/nz)
        ix = (atm_idx_s - iz*(n_bloch_y*n_bloch_x*n_uc/nz) - iy*(n_bloch_x*n_uc/nz) - 1)&
             /(n_uc/nz)
        is = atm_idx_s - iz*(n_bloch_y*n_bloch_x*n_uc/nz) - iy*(n_bloch_x*n_uc/nz) &
            -ix * (n_uc/nz)
        i = istart-1 + iz*n_uc/nz + is
                    m1 = mass_no(int(pos_s(atm_idx_s,4)))

                    do nn = 1, norb/3
                        atm_jdx_s = nn
                        m2 = mass_no(int(pos_s(atm_jdx_s,4)))
!                    write(*,*) atm_idx_s,atm_jdx_s
                        if (num_mat(atm_idx_s,atm_jdx_s).gt.0) then
                            do temp_i = 1,num_mat(atm_idx_s,atm_jdx_s) 
                                dr = dist_mat(temp_count,1:3)
                                temp_idx = idx_mat(temp_count)
                                temp_count = temp_count + 1 
                                if (temp_idx.gt.0)then
                                    do icart = 1,3
                                        idx_s = (atm_idx_s-1)*3 + icart
                                        do jcart = 1,3
                                            jdx_s = (atm_jdx_s-1)*3 + jcart
                                            Ham(idx_s,jdx_s) = Ham(idx_s,jdx_s) + &
                                               force_constant((i-istart)*3+icart,(temp_idx-1)*3+jcart)*&
                                               exp(i_imag*dot_product(dr,kpoint_cart))&
                                               /sqrt(m1*m2)
!                         if((atm_idx_s.eq.7).and.(atm_jdx_s.eq.9).and.(icart.eq.1).and.(jcart.eq.1))then
!                          write(*,*)dr,force_constant((i-istart)*3+icart,(temp_idx-1)*3+jcart) 
!                end if
                                            end do
                                    end do
                                end if
                            end do
                        end if
               !          if ((atm_idx_s.eq.7).and.(atm_jdx_s.eq.9))then
               !     stop
               ! end if
                    end do
    end do

    Ham(:,:) = Ham(:,:)*eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2 ! THz^2
    Ham(:,:) = hermitian(Ham)
    
    !call eigenH(Ham,eig,vec)
    !write(*,*)  my_sqrt(eig(sorting_index(eig)))
    !write(*,*) kpoint_cart
    !stop

    !call eigen(Ham,eigc,vec)
    !write(*,*)  my_sqrt(abs(eigc(sorting_index(abs(eigc)))))


    end subroutine

    subroutine gen_dyn_lead_mpi(ie,ik_core,force_constant,Ham,Ham_v,eigvec,&
                           eig,norb,pos,pos_s,cell,&
                           idx_mat,idx_mat_v,dist_mat,&
                           dist_mat_v,num_mat,num_mat_v)


    use config

    implicit none

    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ie,ik_core
    integer(kind=4),intent(in) :: norb
    integer(kind=4) :: norb_uc
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    complex(kind=8),intent(out) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(out) :: Ham_v(numprocs,norb,norb)
    complex(kind=8),intent(out) :: eigvec(numprocs,norb,norb)
    real(kind=8),intent(out) :: eig(numprocs,norb)
    integer(kind=4),intent(in) :: idx_mat(:),idx_mat_v(:)
    real(kind=8),intent(in) :: dist_mat(:,:),dist_mat_v(:,:)
    integer(kind=4),intent(in) :: num_mat(:,:),num_mat_v(:,:)

    mm = myid+1
    nn = k_start_core(mm)+ik_core-1

    call gen_dyn_lead(ie,k_grid(nn,:),force_constant,&
         Ham(mm,:,:),Ham_v(mm,:,:),eigvec(mm,:,:),eig(mm,:),norb,pos,pos_s,&
         cell,idx_mat,idx_mat_v,dist_mat,dist_mat_v,num_mat,num_mat_v)

    end subroutine gen_dyn_lead_mpi


    subroutine gen_dyn_lead(ie,kpoint,force_constant,Ham,Ham_v,vec,eig,norb,pos,pos_s,&
                            cell,idx_mat,idx_mat_v,dist_mat,dist_mat_v,num_mat,num_mat_v)

    use config
    use input
    use util

    implicit none

    integer(kind=4),intent(in) :: ie
    real(kind=8),intent(in) :: kpoint(3)
    real(kind=8) :: kpoint_cart(3)
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    integer(kind=4),intent(in) :: idx_mat(:),idx_mat_v(:)
    real(kind=8),intent(in) :: dist_mat(:,:),dist_mat_v(:,:)
    integer(kind=4),intent(in) :: num_mat(:,:),num_mat_v(:,:)
    integer(kind=4) :: norb,norb_uc
    complex(kind=8) :: Ham(norb,norb)
    complex(kind=8) :: Ham_v(norb,norb)
    complex(kind=8),intent(out) :: vec(norb,norb)
    real(kind=8),intent(out) :: eig(norb)
    integer(kind=4) :: n_uc, n_sc
    integer(kind=4) :: istart
    integer(kind=4) :: i,j,atm_idx,atm_jdx
    integer(kind=4) :: ii,jj,kk,i_uc,idx,nn
    real(kind=8)  :: m1,m2
    real(kind=8) :: dr(3)
    integer(kind=4) :: ix,iy,iz
    integer(kind=4) :: idx_s, is, jdx_s
    integer(kind=4) :: atm_idx_s,atm_jdx_s
    real(kind=8) :: pos_temp(3),r_temp(3)
    integer(kind=4) :: dx,dy,dz
    integer(kind=4) :: icart,jcart
    integer(kind=4) :: temp_count,temp_count_v,temp_i
    integer(kind=4) :: temp_idx,temp_idx_v
   

    n_sc = size(force_constant,2)/3
    norb_uc = norb/n_bloch_x/n_bloch_y
    n_uc = norb_uc / 3
    
    istart = n_uc*(n_fc_uc_x*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)&
         +n_fc_uc_y*(2*n_fc_uc_z+1)+n_fc_uc_z)+1

    Ham = 0.0d0 
    Ham_v = 0.0d0
    eig = 0.0d0
    vec = 0.0d0

    kpoint_cart = matmul(kpoint,recivec) ! new to convert
    temp_count = 1
    temp_count_v = 1

    do atm_idx_s = 1,norb/3
        iy = (atm_idx_s - 1)/(n_bloch_x*n_uc)
        ix = (atm_idx_s - iy*(n_bloch_x*n_uc) - 1)&
             /(n_uc)
        is = atm_idx_s  - iy*(n_bloch_x*n_uc)&
            -ix * (n_uc)
        i = istart-1 + is
                m1 = mass_no(int(pos_s(atm_idx_s,4)))
                do nn = 1, norb/3
                    atm_jdx_s = nn
                    m2 = mass_no(int(pos_s(atm_jdx_s,4)))
                    if (num_mat(atm_idx_s,atm_jdx_s).gt.0) then
                          do temp_i = 1,num_mat(atm_idx_s,atm_jdx_s) 
                            dr = dist_mat(temp_count,1:3)
                            temp_idx = idx_mat(temp_count)
                            temp_count = temp_count + 1
                            if (temp_idx .gt. 0) then
                                do icart = 1,3
                                    idx_s = (atm_idx_s-1)*3 + icart
                                    do jcart = 1,3
                                        jdx_s = (atm_jdx_s-1)*3 + jcart               
                                        Ham(idx_s,jdx_s) = Ham(idx_s,jdx_s) + &
                                           force_constant((i-istart)*3+icart,(temp_idx-1)*3+jcart)*&
                                           exp(i_imag*dot_product(dr,kpoint_cart))&
                                           /sqrt(m1*m2)
!                                           write(*,*)icart,jcart,dr,force_constant((i-istart)*3+icart,(temp_idx-1)*3+jcart)
                                    end do
                                end do
                            end if
                        end do
                    end if
                    if (num_mat_v(atm_idx_s,atm_jdx_s).gt.0) then
                        do temp_i = 1,num_mat_v(atm_idx_s,atm_jdx_s)
                            dr = dist_mat_v(temp_count_v,1:3)
                            temp_idx_v = idx_mat_v(temp_count_v)
                            temp_count_v = temp_count_v + 1
                            if (temp_idx_v .gt. 0) then
                                do icart = 1,3
                                    idx_s = (atm_idx_s-1)*3 + icart
                                    do jcart = 1,3
                                        jdx_s = (atm_jdx_s-1)*3 + jcart               
                                        Ham_v(idx_s,jdx_s) = Ham_v(idx_s,jdx_s) + &
                                           force_constant((i-istart)*3+icart,(temp_idx_v-1)*3+jcart)*&
                                           exp(i_imag*dot_product(dr,kpoint_cart))&
                                           /sqrt(m1*m2)
!                if((atm_idx_s.eq.3).and.(atm_jdx_s.eq.1).and.(icart.eq.1).and.(jcart.eq.1))then
!                          write(*,*)'lead',dr,force_constant((i-istart)*3+icart,(temp_idx-1)*3+jcart)
!                end if
                                    end do
                                end do
                            end if
                        end do
                    end if
!                if ((atm_idx_s.eq.3).and.(atm_jdx_s.eq.1))then
!                    stop
!                end if

                end do
    end do
!    stop
   

    Ham(:,:) = Ham(:,:)*eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2 ! THz^2
    Ham(:,:) = hermitian(Ham)
    Ham_v(:,:) = Ham_v(:,:)*eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2 ! THz^2

    end subroutine gen_dyn_lead

    subroutine gen_ham_slab(norb_slab,&
                           k_uc_cart,&
                           dyn_slab,&
                           eigs,evs,&
                           latvec_slab,force_constant,&
                           pos,pos_s,&
                           cell,idx_mat,idx_mat_v,dist_mat,&
                           dist_mat_v,num_mat,num_mat_v,&
                           norb)
    
    use config
    use input
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3
    real(kind=8) :: r(3),dr(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    real(kind=8)  :: m1,m2
    integer(kind=4) :: l=1
    real(kind=8),intent(in)  :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_slab
    real(kind=8),intent(in)  :: force_constant(:,:)
    complex(kind=8),intent(out) :: dyn_slab(norb_slab,norb_slab)
    complex(kind=8) :: Ham(norb,norb)
    complex(kind=8) :: Ham_v(norb,norb)
    real(kind=8),intent(out) :: eigs(norb_slab)
    complex(kind=8),intent(out) :: evs(norb_slab,norb_slab)
    real(kind=8),intent(in) :: latvec_slab(3,3)
    integer(kind=4) :: natoms_slab
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    integer(kind=4),intent(in) :: idx_mat(:),idx_mat_v(:)
    real(kind=8),intent(in) :: dist_mat(:,:),dist_mat_v(:,:)
    integer(kind=4),intent(in) :: num_mat(:,:),num_mat_v(:,:)
    integer(kind=4) :: ix,iy,iz,is
    integer(kind=4) :: norb,norb_uc
    integer(kind=4) :: n_sc,n_uc,nn
    integer(kind=4) :: temp_count,temp_count_v
    integer(kind=4) :: atm_idx,atm_jdx,idx_s,jdx_s
    integer(kind=4) :: atm_idx_s,atm_jdx_s
    integer(kind=4) :: icart,jcart,temp_idx,temp_idx_v
    integer(kind=4) :: istart,temp_i

    dyn_slab = 0.0d0
    Ham = 0.0d0
    Ham_v = 0.0d0
    eigs = 0.0d0
    evs = 0.0d0

    n_sc = size(force_constant,2)/3
    norb_uc = norb/n_bloch_x/n_bloch_y
    n_uc = norb_uc / 3

    istart = n_uc*(n_fc_uc_x*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)&
         +n_fc_uc_y*(2*n_fc_uc_z+1)+n_fc_uc_z)+1

    temp_count = 1
    temp_count_v = 1

    do atm_idx_s = 1,norb/3
        iy = (atm_idx_s - 1)/(n_bloch_x*n_uc)
        ix = (atm_idx_s - iy*(n_bloch_x*n_uc) - 1)&
             /(n_uc)
        is = atm_idx_s  - iy*(n_bloch_x*n_uc)&
            -ix * (n_uc)
        i = istart-1 + is
        m1 = mass_no(int(pos_s(atm_idx_s,4)))
                do nn = 1, norb/3
                    atm_jdx_s = nn
                    m2 = mass_no(int(pos_s(atm_jdx_s,4)))
                    if (num_mat(atm_idx_s,atm_jdx_s).gt.0) then
                          do temp_i = 1,num_mat(atm_idx_s,atm_jdx_s) 
                            dr = dist_mat(temp_count,1:3)
                            temp_idx = idx_mat(temp_count)
                            temp_count = temp_count + 1
                            if (temp_idx .gt. 0) then
                                do icart = 1,3
                                    idx_s = (atm_idx_s-1)*3 + icart
                                    do jcart = 1,3
                                        jdx_s = (atm_jdx_s-1)*3 + jcart               
                                        Ham(idx_s,jdx_s) = Ham(idx_s,jdx_s) + &
                                           force_constant((i-istart)*3+icart,(temp_idx-1)*3+jcart)*&
                                           exp(i_imag*dot_product(dr,k_uc_cart))&
                                           /sqrt(m1*m2)
                                    end do
                                end do
                            end if
                        end do
                    end if
                    if (num_mat_v(atm_idx_s,atm_jdx_s).gt.0) then
                        do temp_i = 1,num_mat_v(atm_idx_s,atm_jdx_s)
                            dr = dist_mat_v(temp_count_v,1:3)
                            temp_idx_v = idx_mat_v(temp_count_v)
                            temp_count_v = temp_count_v + 1
                            if (temp_idx_v .gt. 0) then
                                do icart = 1,3
                                    idx_s = (atm_idx_s-1)*3 + icart
                                    do jcart = 1,3
                                        jdx_s = (atm_jdx_s-1)*3 + jcart               
                                        Ham_v(idx_s,jdx_s) = Ham_v(idx_s,jdx_s) + &
                                           force_constant((i-istart)*3+icart,(temp_idx_v-1)*3+jcart)*&
                                           exp(i_imag*dot_product(dr,k_uc_cart))&
                                           /sqrt(m1*m2)
                                    end do
                                end do
                            end if
                        end do
                    end if

                end do
    end do
    Ham(:,:) = Ham(:,:)*eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2 ! THz^2
    Ham(:,:) = hermitian(Ham)
    Ham_v(:,:) = Ham_v(:,:)*eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2 ! THz^2
    
    dyn_slab(1:norb,1:norb) = Ham(:,:)
    dyn_slab(1:norb,norb+1:2*norb) = Ham_v(:,:)
    dyn_slab((nthick-1)*norb+1:nthick*norb,&
    (nthick-1)*norb+1:nthick*norb) = Ham(:,:)
    dyn_slab((nthick-1)*norb+1:nthick*norb,&
    (nthick-2)*norb+1:(nthick-1)*norb) = transpose(dconjg(Ham_v))
 
    do i = 1,nthick-2
        dyn_slab(i*norb+1:(i+1)*norb,i*norb+1:(i+1)*norb) = Ham
        dyn_slab(i*norb+1:(i+1)*norb,(i+1)*norb+1:(i+2)*norb) = Ham_v
        dyn_slab(i*norb+1:(i+1)*norb,(i-1)*norb+1:i*norb) = &
        transpose(dconjg(Ham_v))
    end do

    dyn_slab = hermitian(dyn_slab)


    call eigenH(dyn_slab,eigs,evs)

    end subroutine gen_ham_slab

    subroutine gen_dyn_pc(norb_pc,&
                          kpoint_cart,dyn_pc,&
                          force_constant,pos_pc,&
                          pos_pc_fc_sc)

    use config, only : r_cutoff
    use input
    use util

    implicit none

    include "mpif.h" 

    integer(kind=4) :: norb_pc
    real(kind=8) :: kpoint_cart(3)
    real(kind=8),intent(in) :: force_constant(:,:)
    complex(kind=8),intent(out) :: dyn_pc(norb_pc,norb_pc)
    real(kind=8),intent(in) :: pos_pc(:,:)
    real(kind=8),intent(in) :: pos_pc_fc_sc(:,:)
    integer(kind=4) :: i,j,icart,jcart
    integer(kind=4) :: nx,ny,nz
    integer(kind=4) :: ix,iy,iz,idx
    integer(kind=4) :: natm_pc
    real(kind=8) :: pos_diff(3)
    real(kind=8) :: m1,m2

    dyn_pc(:,:) = 0.0d0
    nx = 1
    ny = 1
    nz = 1
    natm_pc = norb_pc/3
!    istart = norb_pc*(nx*(2*ny+1)*(2*nz+1)&
!             +ny*(2*nz+1)+nz)+1
    do i = 1,natm_pc
        m1 = mass_no(int(pos_pc(i,4))) 
        do j = 1,natm_pc
            m2 = mass_no(int(pos_pc(j,4)))
            do ix = -nx,nx
                do iy = -ny,ny
                    do iz = -nz,nz
                        idx = ((ix+nx)*(2*ny+1)*(2*nz+1)&
                          +(iy+ny)*(2*nz+1)+iz+nz)*natm_pc+j
                        pos_diff(:) = pos_pc_fc_sc(idx,1:3) - pos_pc(i,1:3) 
                        if (sqrt(dot_product(pos_diff,pos_diff)).lt.r_cutoff) then
                            do icart = 1,3
                                do jcart = 1,3
                                    dyn_pc(3*(i-1)+icart,3*(j-1)+jcart) = dyn_pc(3*(i-1)+icart,3*(j-1)+jcart) &
                                    + force_constant(3*(i-1)+icart,3*(idx-1)+jcart) &
                                    * exp(i_imag*dot_product(pos_diff,kpoint_cart))/sqrt(m1*m2) 
!                                    write(*,*)force_constant(3*(i-1)+icart,3*(idx-1)+jcart)
                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do
    end do
    dyn_pc(:,:) = dyn_pc(:,:)*eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2
    dyn_pc(:,:) = hermitian(dyn_pc)

    end subroutine gen_dyn_pc

    subroutine group_velocity_pc(norb_pc,&
                          kpoint_cart,dyn_pc,&
                          force_constant,pos_pc,&
                          pos_pc_fc_sc,&
                          eigs,evs,vels)

    use config, only : r_cutoff
    use input
    use util

    implicit none

    include "mpif.h" 

    integer(kind=4) :: norb_pc
    real(kind=8) :: kpoint_cart(3)
    real(kind=8),intent(in) :: force_constant(:,:)
    complex(kind=8),intent(out) :: dyn_pc(norb_pc,norb_pc)
    complex(kind=8) :: dyn_pc_p(norb_pc,norb_pc)
    complex(kind=8) :: dyn_pc_m(norb_pc,norb_pc)      
    complex(kind=8) :: dHdk(norb_pc,norb_pc)
    complex(kind=8) :: Vmat(norb_pc,norb_pc)
    real(kind=8) :: devi = 1.0d-5
    real(kind=8),intent(in) :: pos_pc(:,:)
    real(kind=8),intent(in) :: pos_pc_fc_sc(:,:)
    real(kind=8),intent(out) :: vels(norb_pc,3)
    integer(kind=4) :: i,j,icart,jcart
    integer(kind=4) :: nx,ny,nz
    integer(kind=4) :: ix,iy,iz,idx
    integer(kind=4) :: natm_pc
    real(kind=8) :: pos_diff(3)
    real(kind=8) :: eigs(norb_pc)
    complex(kind=8) :: evs(norb_pc,norb_pc)
    real(kind=8) :: ktemp1(3),ktemp2(3),delta_k
   
    call gen_dyn_pc(norb_pc,&
                    kpoint_cart,dyn_pc,&
                    force_constant,pos_pc,&
                    pos_pc_fc_sc)
!    write(*,*) maxval(abs(dyn_pc - transpose(dconjg(dyn_pc))))
!    stop
    call eigenH(dyn_pc,eigs,evs)
!    write(*,*)'eigs', eigs
    dyn_pc_p = 0.0d0
    dyn_pc_m = 0.0d0
    dHdk = 0.0d0

    do i = 1,3
        delta_k = devi

        ktemp1 = kpoint_cart
        ktemp1(i) = ktemp1(i)-delta_k
        ktemp2 = kpoint_cart
        ktemp2(i) = ktemp2(i)+delta_k

        call gen_dyn_pc(norb_pc,&
                      ktemp1,dyn_pc_m,&
                      force_constant,pos_pc,&
                      pos_pc_fc_sc)
        call gen_dyn_pc(norb_pc,&
                      ktemp2,dyn_pc_p,&
                      force_constant,pos_pc,&
                      pos_pc_fc_sc)

        dHdk = (dyn_pc_p-dyn_pc_m)/(2.0*delta_k)
        Vmat = matmul(matmul(&
               transpose(dconjg(evs)),dHdk),evs)
        do j = 1,norb_pc
            vels(j,i) = real(Vmat(j,j))/(2*sqrt(abs(eigs(j)))) ! v_x/y/z
        end do
    end do

    end subroutine group_velocity_pc

    subroutine group_velocity_pc_qe(norb_pc,&
                          kpoint_cart,dyn_pc,&
                          force_constant,pos_pc,&
                          pos_pc_fc_sc,&
                          eigs,evs,vels,mass_uc)

    use config, only : r_cutoff
    use input
    use util

    implicit none

    include "mpif.h" 

    integer(kind=4) :: norb_pc
    real(kind=8) :: kpoint_cart(3)
    real(kind=8),intent(in) :: force_constant(:,:)
    complex(kind=8),intent(out) :: dyn_pc(norb_pc,norb_pc)
    complex(kind=8) :: dyn_pc_p(norb_pc,norb_pc)
    complex(kind=8) :: dyn_pc_m(norb_pc,norb_pc)      
    complex(kind=8) :: dHdk(norb_pc,norb_pc)
    complex(kind=8) :: Vmat(norb_pc,norb_pc)
    real(kind=8) :: devi = 1.0d-5
    real(kind=8),intent(in) :: pos_pc(:,:)
    real(kind=8),intent(in) :: pos_pc_fc_sc(:,:)
    real(kind=8),intent(in) :: mass_uc(:)
    real(kind=8),intent(out) :: vels(norb_pc,3)
    integer(kind=4) :: i,j,icart,jcart
    integer(kind=4) :: nx,ny,nz
    integer(kind=4) :: ix,iy,iz,idx
    integer(kind=4) :: natm_pc
    real(kind=8) :: pos_diff(3)
    real(kind=8) :: eigs(norb_pc)
    complex(kind=8) :: evs(norb_pc,norb_pc)
    real(kind=8) :: ktemp1(3),ktemp2(3),delta_k
   
    call gen_dyn_pc_qe_k(kpoint_cart,dyn_pc,mass_uc)
!    write(*,*) maxval(abs(dyn_pc - transpose(dconjg(dyn_pc))))
!    stop
    call eigenH(dyn_pc,eigs,evs)
!    write(*,*)'eigs', eigs
    dyn_pc_p = 0.0d0
    dyn_pc_m = 0.0d0
    dHdk = 0.0d0

    do i = 1,3
        delta_k = devi

        ktemp1 = kpoint_cart
        ktemp1(i) = ktemp1(i)-delta_k
        ktemp2 = kpoint_cart
        ktemp2(i) = ktemp2(i)+delta_k

        call gen_dyn_pc_qe_k(ktemp1,dyn_pc_m,mass_uc)
        call gen_dyn_pc_qe_k(ktemp2,dyn_pc_p,mass_uc)

        dHdk = (dyn_pc_p-dyn_pc_m)/(2.0*delta_k)
        Vmat = matmul(matmul(&
               transpose(dconjg(evs)),dHdk),evs)
        do j = 1,norb_pc
            vels(j,i) = real(Vmat(j,j))/(2*sqrt(abs(eigs(j)))) ! v_x/y/z
        end do
    end do

    end subroutine group_velocity_pc_qe


    subroutine compare_evs(E,norb_pc,n,&
                          kpoint_cart,norb_uc,u,pos,&
                          k_sc_cart,k_uc_cart,ksc,i_deg_pc,&
                          weight,vxyz,force_constant,&
                          pos_pc_fc_sc,pos_pc,&
                          natoms_uc,natoms_pc,mass_uc,&
                          branch_idx_pc)

    use config
    use util

    implicit none

    integer(kind=4),intent(in) :: norb_pc
    real(kind=8),intent(in) :: E
    real(kind=8),intent(in) :: kpoint_cart(3)
    real(kind=8),intent(in) :: pos(:,:)
    integer(kind=4),intent(in) :: n
    complex(kind=8),intent(in):: u(n)
    real(kind=8) :: pos_shiftz(n,4)
    real(kind=8),intent(in) :: k_sc_cart(3)
    real(kind=8),intent(in) :: k_uc_cart(3)
    real(kind=8),intent(in) :: ksc(3)
    complex(kind=8)         :: dyn_pc(norb_pc,norb_pc)
    real(kind=8),intent(in) :: pos_pc(:,:)
    integer(kind=4),intent(in) :: norb_uc
    integer(kind=4),intent(out) :: i_deg_pc
    real(kind=8),intent(out) :: weight
    real(kind=8),intent(out) :: vxyz(3)
    integer(kind=4),intent(out) :: branch_idx_pc
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos_pc_fc_sc(:,:)
    integer(kind=4),intent(in) :: natoms_uc,natoms_pc 
    real(kind=8),intent(in)    :: mass_uc(:)
    real(kind=8) :: eigs(norb_pc)
    complex(kind=8) :: evs(norb_pc,norb_pc)
    complex(kind=8) :: u_pc(norb_pc)
    complex(kind=8) :: u_sc(n)
    real(kind=8) :: vels(norb_pc,3)
    complex(kind=8),allocatable :: evs_uc(:)
    complex(kind=8),allocatable :: evs_lead(:)
    complex(kind=8) :: temp_sum
    integer(kind=4) :: npc
    real(kind=8) :: Eq
    integer(kind=4) :: i,j,j1,l,k,id,idx

    Eq = sqrt(E)

    if (qe_fc .eq. 0) then
        call group_velocity_pc(norb_pc,&
                       kpoint_cart,dyn_pc,&
                        force_constant,pos_pc,&
                        pos_pc_fc_sc,&
                        eigs,evs,vels)
    else
        call group_velocity_pc_qe(norb_pc,&
                       kpoint_cart,dyn_pc,&
                        force_constant,pos_pc,&
                        pos_pc_fc_sc,&
                        eigs,evs,vels,mass_uc)
    end if
    npc = natoms_uc/natoms_pc
    allocate(evs_uc(norb_pc*npc))
    allocate(evs_lead(norb_uc*n_bloch_x*n_bloch_y))
    evs_lead = 0.0d0
    evs_uc = 0.0d0
    u_pc = 0.0d0
    u_sc = 0.0d0

    pos_shiftz  = pos
    pos_shiftz(:,3) = pos_shiftz(:,3) - minval(pos(:,3))
    i_deg_pc = 0
    weight = 0.0d0
    vxyz = 0.0d0
!    write(*,*)'eig',norb_pc, kpoint_cart
!    write(*,*) minval(abs(Eq-sqrt(abs(eigs(:)))))

    do i =1,norb_pc
        if (abs(Eq-sqrt(abs(eigs(i)))).lt.1.0d-2)then
            vxyz = vels(i,:)
            i_deg_pc = 1
            branch_idx_pc = i

            do j = 1,n_bloch_y
                do j1 = 1,n_bloch_x
                    do l = 1,npc
                        do k = 1,norb_pc
                            id = (((j-1)*n_bloch_x+j1-1)*npc+l-1)*norb_pc+k
                            idx = int(pos_shiftz(id,4))
                            evs_lead(id) =&
                            evs(k,i)*exp(i_imag*(dot_product(pos_shiftz(id,1:3), &
                            kpoint_cart)))
                            u_sc(id) = u(id) *&
                            exp(i_imag*(dot_product(pos_shiftz(id,1:3),ksc)))
                        end do
                    end do
                 end do
            end do
            evs_lead = evs_lead / sqrt(dble(n_bloch_x*n_bloch_y*npc))

            temp_sum = 0.0
            temp_sum = dot_product(u_sc(:),evs_lead(:))
            weight = abs(temp_sum)**2
            !if (abs(temp_sum) .lt. 0.8d0) then
            !    i_deg_pc = i_deg_pc - 1
            !end if
        end if
    end do            

    end subroutine compare_evs

    subroutine find_k_pc(E,norb_pc,&
                          k_sc_cart,&
                          k_uc_cart,&
                          norb_uc,n,u,pos,&
                          ksc,kpc_unfold,&
                          i_deg_pc_tot,vpc_unfold,&
                          reci_uc,recivec_pc,&
                          force_constant,&
                          pos_pc_fc_sc,pos_pc,&
                          natoms_uc,natoms_pc,mass_uc,&
                          branch_idx_pc)

    use config
    use util

    implicit none

    integer(kind=4),intent(in) :: norb_pc 
    real(kind=8),intent(in) :: E
    integer(kind=4) :: l,i
    real(kind=8),intent(in) :: k_sc_cart(3)
    real(kind=8),intent(in) :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_uc
    integer(kind=4),intent(in) :: n
    real(kind=8),intent(in) :: pos(n,3)
    complex(kind=8),intent(in) :: u(n)
    complex(kind=8) :: u1(norb_uc,1)
    complex(kind=8) :: u2(norb_uc,1)
    complex(kind=8) :: u3(norb_uc,1)
    real(kind=8),intent(in) :: ksc(3)
    real(kind=8),intent(in) :: mass_uc(:)
    real(kind=8)    :: force_constant(:,:) 
    integer(kind=4) :: i1,i2,i3,ibz
    real(kind=8) :: kpc(3),weight
    integer(kind=4) :: i_deg_pc = 0
    integer(kind=4),intent(out) :: i_deg_pc_tot
    real(kind=8),intent(out) :: kpc_unfold(125,3)
    real(kind=8),intent(out) :: vpc_unfold(125,3)
    integer(kind=4),intent(out) :: branch_idx_pc
    real(kind=8) :: vxyz(3)
    real(kind=8) :: v_pc(norb_pc,3)
    real(kind=8),intent(in)  :: pos_pc_fc_sc(:,:)
    real(kind=8),intent(in)  :: pos_pc(:,:)
    complex(kind=8) :: dyn_uc(norb_uc,norb_uc)
    real(kind=8) :: eigs(norb_uc)
    complex(kind=8) :: evs(norb_uc,norb_uc)
    complex(kind=8) :: dyn_pc(norb_pc,norb_pc)
    real(kind=8) :: eigs_pc(norb_pc)
    complex(kind=8) :: evs_pc(norb_pc,norb_pc)

    real(kind=8) :: reci_uc(3,3),recivec_pc(3,3)
    integer(kind=4) :: natoms_uc,natoms_pc

    natoms_pc = norb_pc/3

    l = max(n_bloch_x,n_bloch_y)
    i_deg_pc_tot = 0
    kpc_unfold(:,:) = 0.0d0
    vpc_unfold(:,:) = 0.0d0
    weight = 0.0
    vxyz = 0.0d0
    
    ibz = in_ws(reci_uc,k_uc_cart)
!    write(*,*)matmul(k_uc_cart,inv_real(reci_uc)),ibz 
    if (ibz.eq.2) then
        i_deg_pc = 1
        i_deg_pc_tot = i_deg_pc + i_deg_pc_tot
        kpc_unfold(i_deg_pc_tot,:) = k_uc_cart
        if (qe_fc.eq. 0) then
            call group_velocity_pc(norb_pc,&
                        k_uc_cart,dyn_pc,&
                        force_constant,pos_pc,&
                        pos_pc_fc_sc,&
                        eigs_pc,evs_pc,v_pc)
        else
            call group_velocity_pc_qe(norb_pc,&
                        k_uc_cart,dyn_pc,&
                        force_constant,pos_pc,&
                        pos_pc_fc_sc,&
                        eigs_pc,evs_pc,v_pc,mass_uc)
        end if
        do i1 = 1,norb_pc
            if(abs(sqrt(E)-sqrt(abs(eigs_pc(i1)))).lt.1.0d-2)then
                vxyz = v_pc(i1,:)
                vpc_unfold(1,:) = vxyz(:)
                branch_idx_pc = i1
            end if
        end do
    else
        do i1 = -l,l
            do i2 = -l,l
                do i3 = -l,l
                    kpc(:) = k_uc_cart(:) &
                           + matmul((/dble(i1),dble(i2),dble(i3)/),reci_uc)
                    ibz = in_ws(recivec_pc,kpc)
                    if (ibz.ne.0) then
                        call compare_evs(E,norb_pc,n,&
                              kpc,norb_uc,u,pos,&
                              k_sc_cart,k_uc_cart,ksc,i_deg_pc,&
                              weight,vxyz,force_constant,&
                              pos_pc_fc_sc,pos_pc,&
                              natoms_uc,natoms_pc,mass_uc,&
                              branch_idx_pc)
                        if (i_deg_pc .gt. 0) then
                            i_deg_pc_tot = i_deg_pc + i_deg_pc_tot
                            kpc_unfold(i_deg_pc_tot,:) = kpc(:)
                            vpc_unfold(i_deg_pc_tot,:) = vxyz(:)
                        end if 
                    end if
                end do
            end do
        end do
    end if
!    if (i_deg_pc_tot .eq.0 )then
!        write(*,*) 'one K point does not find the corresponding k in FBZ of left PC'
!        write(*,*) 'sc',ksc
!        write(*,*) 'uc',k_uc_cart
!        write(*,*) 'uc_cry',sqrt(E),matmul(k_uc_cart,inv_real(reci_uc))
!        write(*,*)inv_real(reci_uc)
!        write(*,*) sqrt(E),'vel',i_deg_pc_tot,vxyz,kpc
!        stop
!    end if

    end subroutine find_k_pc
    
    subroutine swap_lead_region_dyn_mpi(ie,ik_core,Ham,Ham_l,Ham_r,Ham_l_v,Ham_r_v)

    use config

    implicit none

    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ie,ik_core
    complex(kind=8) :: Ham(:,:,:)
    complex(kind=8) :: Ham_l(:,:,:)
    complex(kind=8) :: Ham_r(:,:,:)
    complex(kind=8) :: Ham_l_v(:,:,:)
    complex(kind=8) :: Ham_r_v(:,:,:)


    mm = myid+1
    nn = k_start_core(mm)+ik_core-1

    call swap_lead_region_dyn(Ham(mm,:,:),&
                             Ham_l(mm,:,:),Ham_r(mm,:,:),&
                             Ham_l_v(mm,:,:),Ham_r_v(mm,:,:))

    end subroutine swap_lead_region_dyn_mpi
    
    subroutine swap_lead_region_dyn(Ham,Ham_l,Ham_r,Ham_l_v,Ham_r_v)

    use config
    use input
    use util

    implicit none

    complex(kind=8) :: Ham(:,:)
    complex(kind=8) :: Ham_l(:,:)
    complex(kind=8) :: Ham_r(:,:)
    complex(kind=8) :: Ham_l_v(:,:)
    complex(kind=8) :: Ham_r_v(:,:)
    integer(kind=4) :: i, lefta,leftb, leftc,&
                       righta, rightb, rightc


    lefta = 3*buffer_left+1
    leftc = 3*buffer_left+6*period_left
    righta = 3*natoms-6*period_right-3*buffer_right+1
    rightc = 3*natoms-3*buffer_right

!    Ham(lefta:leftc,lefta:leftc) = 0.0d0
!    Ham(righta:rightc,righta:rightc) = 0.0d0

!    write(*,*) size(Ham(lefta+3*(i-1)*period_left:lefta-1+3*i*period_left,&
!                lefta+3*(i-1)*period_left:lefta-1+3*i*period_left),1)
!    write(*,*) minval(abs(Ham(lefta+3*(i-1)*period_left:lefta-1+3*i*period_left,&
!                lefta+3*(i-1)*period_left:lefta-1+3*i*period_left)))
!    write(*,*) minval(abs(Ham_l))
!    write(*,*) maxval(abs(Ham(lefta+3*(i-1)*period_left:lefta-1+3*i*period_left,&
!                lefta+3*(i-1)*period_left:lefta-1+3*i*period_left)-Ham_l))
    !write(*,*) maxval(abs( Ham(lefta:leftb,lefta:leftb) - Ham_l))
    !write(*,*)  maxloc(abs( Ham(lefta:leftb,lefta:leftb) - Ham_l))


    do i = 1,2
        Ham(lefta+3*(i-1)*period_left:lefta-1+3*i*period_left,&
            lefta+3*(i-1)*period_left:lefta-1+3*i*period_left)&
        = Ham_l
        Ham(rightc+1-3*i*period_right:rightc-3*(i-1)*period_right,&
            rightc+1-3*i*period_right:rightc-3*(i-1)*period_right) &
        = Ham_r
    end do

    leftb = 3*buffer_left+3*period_left

!    write(*,*) lefta
!    write(*,*) maxloc(abs(Ham(lefta:leftb,leftb+1:leftc) - Ham_l_v))
!    stop
    Ham(lefta:leftb,leftb+1:leftc) = Ham_l_v
    Ham(leftb+1:leftc,lefta:leftb) = transpose(dconjg(Ham_l_v)) 
    ! ideally, the device should be thick enough so that ASR are
    ! satisfied, for electron the onsite energy is physical so device can be
    ! shorter

    rightb = 3*natoms-3*period_right-3*buffer_right 

    Ham(righta:rightb,rightb+1:rightc) = Ham_r_v
    Ham(rightb+1:rightc,righta:rightb) = transpose(dconjg(Ham_r_v))
    if (maxval(abs(Ham( 3*buffer_left+3*period_left+1: 3*buffer_left+2*3*period_left,&
            3*buffer_left+1:3*buffer_left+3*period_left)-transpose(dconjg(Ham(3*buffer_left+1:3*buffer_left+3*period_left,3*buffer_left+3*period_left+1:3*buffer_left+2*3*period_left))))) .gt.0.00001) then
            stop
        end if

 

    Ham(:,:) = hermitian(Ham)

    end subroutine swap_lead_region_dyn

    subroutine  gen_fc_to_dyn_idx(idx_mat,dist_mat,num_mat,force_constant,&
                                  pos,pos_s,cell,nbx,nby,norb)

    use config
    use input
    use util


    implicit none

    include "mpif.h"

    integer(kind=4),intent(out) :: idx_mat(:,:,:)
    real(kind=8),intent(out) :: dist_mat(:,:,:,:)
    integer(kind=4),intent(out) :: num_mat(:,:)
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    integer(kind=4),intent(in) :: nbx,nby
    integer(kind=4) :: norb, norb_uc
    integer(kind=4) :: n_uc, n_sc
    integer(kind=4) :: istart
    integer(kind=4) :: i,j,atm_idx,atm_jdx
    integer(kind=4) :: jdx_mod
    integer(kind=4) :: ii,jj,kk,i_uc,idx,nn
    real(kind=8)  :: m1,m2
    integer(kind=4) :: ndiag,ll
    real(kind=8) :: dr(3),dr1(3),dr2(3)
    real(kind=8) :: r_shift(3)
    integer(kind=4) :: ix,iy,iz
    integer(kind=4) :: idx_s, is, jdx_s
    integer(kind=4) :: atm_idx_s,atm_jdx_s
    real(kind=8) :: pos_temp(3),r_temp(3)
    integer(kind=4) :: dx,dy,dz
    real(kind=8) :: i_cry(3)
    integer(kind=4) :: nb_no

    nb_no = 0
    num_mat = 0
    idx_mat = 0
    dist_mat = 0.0d0
  

    n_sc = size(force_constant,2)/3
    norb_uc = norb/n_bloch_x/n_bloch_y
    n_uc = norb_uc / 3

    
    istart = n_uc*(n_fc_uc_x*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)&
         +n_fc_uc_y*(2*n_fc_uc_z+1)+n_fc_uc_z)+1

    do iy = 1,n_bloch_y
        do ix = 1,n_bloch_x
            do iz = 1,nz 
                do is = 1,n_uc/nz ! current code is not general enough
                    i = istart-1 + (iz-1)*n_uc/nz + is
                    atm_idx_s = (iz-1)*n_bloch_y*n_bloch_x*n_uc/nz+&
                                (iy-1)          *n_bloch_x*n_uc/nz+&
                                (ix-1)                    *n_uc/nz+&
                                is

                    do nn = 1, norb/3
                        atm_jdx_s = nn
                        nb_no = 0
                        do dx = -1,1
                            do dy = -1,1
                                do dz = -1,1
                                    r_shift(:) =  matmul((/dble(dx),&
                                                  dble(dy),dble(dz)/),&
                                                  cell)
                                    dr = pos_s(atm_jdx_s,1:3) + r_shift &
                                       - pos_s(atm_idx_s,1:3)
                                    if (sqrt(dot_product(dr,dr)) .lt. r_cutoff) then
                                        do ii = 1,2*n_fc_uc_x+1
                                            do jj = 1,2*n_fc_uc_y+1
                                                do kk = 1,2*n_fc_uc_z+1
                                                    i_uc = (ii-1)*(2*n_fc_uc_y+1)&
                                                                 *(2*n_fc_uc_z+1)+&
                                                           (jj-1)*(2*n_fc_uc_z+1)+kk
                                                    do j = 1, n_uc ! in fc coordinates
                                                        idx = (i_uc-1)*n_uc+j
                                                        atm_jdx = (i_uc-1)*n_uc+j
                                                        pos_temp = pos(atm_jdx,1:3)
                                                        dr2 = pos_temp - pos(i,1:3)                 
                                                        if (sqrt(dot_product(dr2-dr,dr2-dr)).lt.1.0d-3)then
                                                            nb_no = nb_no + 1 ! equivalent basis
                                                            dist_mat(atm_idx_s,atm_jdx_s,nb_no,1:3) = dr
                                                            idx_mat(atm_idx_s,atm_jdx_s,nb_no) = idx
                                                        end if
                                                    end do
                                                end do
                                            end do
                                        end do
                                    end if
                                end do
                            end do
                        end do
                        num_mat(atm_idx_s,atm_jdx_s) = nb_no
                    end do
                end do
            end do
        end do
    end do

    end subroutine gen_fc_to_dyn_idx
   
    subroutine  gen_lead_fc_to_dyn_idx(idx_mat,idx_mat_v,dist_mat,dist_mat_v,&
                num_mat,num_mat_v,&
                force_constant,pos,pos_s,cell,nbx,nby,norb)

    use config
    use input
    use util


    implicit none

    include "mpif.h"

    integer(kind=4),intent(out) :: idx_mat(:,:,:)
    integer(kind=4),intent(out) :: idx_mat_v(:,:,:)
    real(kind=8),intent(out) :: dist_mat(:,:,:,:)
    real(kind=8),intent(out) :: dist_mat_v(:,:,:,:)
    integer(kind=4),intent(out) :: num_mat(:,:)
    integer(kind=4),intent(out) :: num_mat_v(:,:)
    real(kind=8),intent(in) :: force_constant(:,:)
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8),intent(in) :: pos_s(:,:)
    real(kind=8),intent(in) :: cell(3,3)
    integer(kind=4),intent(in) :: nbx,nby
    integer(kind=4) :: norb, norb_uc
    integer(kind=4) :: n_uc, n_sc
    integer(kind=4) :: istart
    integer(kind=4) :: i,j,atm_idx,atm_jdx
    integer(kind=4) :: jdx_mod
    integer(kind=4) :: ii,jj,kk,i_uc,idx,nn
    real(kind=8)  :: m1,m2
    real(kind=8) :: dr(3),dr1(3),dr2(3)
    real(kind=8) :: r_shift(3)
    integer(kind=4) :: ix,iy
    integer(kind=4) :: idx_s, is, jdx_s
    integer(kind=4) :: atm_idx_s,atm_jdx_s
    real(kind=8) :: pos_temp(3),r_temp(3)
    integer(kind=4) :: dx,dy,dz
    integer(kind=4) :: nb_no,nb_no_v

    idx_mat = 0
    idx_mat_v = 0
    dist_mat = 0.0d0
    dist_mat_v = 0.0d0
    num_mat = 0
    num_mat_v = 0

    n_sc = size(force_constant,2)/3
    norb_uc = norb/n_bloch_x/n_bloch_y
    n_uc = norb_uc / 3
    
    istart = n_uc*(n_fc_uc_x*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)&
         +n_fc_uc_y*(2*n_fc_uc_z+1)+n_fc_uc_z)+1

    do iy = 1,n_bloch_y
        do ix = 1,n_bloch_x
            do is = 1,n_uc ! current code is not g
                i = istart-1 + is
                atm_idx_s = (iy-1)*n_bloch_x*n_uc+&
                            (ix-1)*n_uc+&
                            is
                do nn = 1, norb/3
                    atm_jdx_s = nn
                    nb_no = 0
                    nb_no_v = 0
                    do dx = -1,1
                        do dy = -1,1
                            do dz = -1,1
                                r_shift(:) =  matmul((/dble(dx),&
                                              dble(dy),dble(dz)/),&
                                              cell)
                                dr = pos_s(atm_jdx_s,1:3) + r_shift&
                                   - pos_s(atm_idx_s,1:3)
                                if (sqrt(dot_product(dr,dr)) .lt. r_cutoff) then
                                    do ii = 1,2*n_fc_uc_x+1
                                        do jj = 1,2*n_fc_uc_y+1
                                            kk = n_fc_uc_z + 1 

                                            i_uc = (ii-1)*(2*n_fc_uc_y+1)&
                                                         *(2*n_fc_uc_z+1)+&
                                                   (jj-1)*(2*n_fc_uc_z+1)+kk

                                            do j = 1, n_uc ! in fc coordinates
                                                idx = (i_uc-1)*n_uc+j
                                                atm_jdx = (i_uc-1)*n_uc+j
                                                pos_temp = pos(atm_jdx,1:3)
                                                dr2 = pos_temp - pos(i,1:3)
                                                if (sqrt(dot_product(dr2-dr,dr2-dr)).lt.1.0d-3)then
                                                    nb_no = nb_no + 1 ! equivalent basis
                                                    dist_mat(atm_idx_s,atm_jdx_s,nb_no,1:3) = dr
                                                    idx_mat(atm_idx_s,atm_jdx_s,nb_no) = idx
                                                end if
                                            end do

                                            kk = 2*n_fc_uc_z + 1 
                                            i_uc = (ii-1)*(2*n_fc_uc_y+1)&
                                                         *(2*n_fc_uc_z+1)+&
                                                   (jj-1)*(2*n_fc_uc_z+1)+kk
                                            do j = 1, n_uc ! in fc coordinates
                                                idx = (i_uc-1)*n_uc+j  
                                                atm_jdx = (i_uc-1)*n_uc+j
                                                pos_temp = pos(atm_jdx,1:3)
                                                dr2 = pos_temp - pos(i,1:3) 
                                                if (sqrt(dot_product(dr2-dr,dr2-dr)).lt.1.0d-3)then
                                                    nb_no_v = nb_no_v + 1 ! equivalent basis
                                                    dist_mat_v(atm_idx_s,atm_jdx_s,nb_no_v,1:3) = dr
                                                    idx_mat_v(atm_idx_s,atm_jdx_s,nb_no_v) = idx
                                                end if
                                            end do
                                        end do
                                    end do
                                end if
                            end do
                        end do ! jj
                    end do ! ii
                    num_mat(atm_idx_s,atm_jdx_s) = nb_no
                    num_mat_v(atm_idx_s,atm_jdx_s) = nb_no_v
                end do
            end do
        end do
    end do
    end subroutine gen_lead_fc_to_dyn_idx

    subroutine asr_cutoff(force_constant,pos,rcut,nx,ny,nz)

    implicit none

    include "mpif.h"
    
    real(kind=8) :: force_constant(:,:)
    real(kind=8) :: force_constant_new(size(force_constant,1),size(force_constant,2))
    real(kind=8) :: rcut,temp_sum
    real(kind=8) :: pos(:,:)
    real(kind=8) :: dr(3)
    integer(kind=4) :: m,n,i,j,icart,jcart
    integer(kind=4) :: n_uc,istart
    integer(kind=4) :: nx,ny,nz

    m = size(force_constant,1)/3
    n_uc = m
    istart = n_uc*(nx*(2*ny+1)*(2*nz+1)&
            +ny*(2*nz+1)+nz)+1
    n = size(force_constant,2)/3

    force_constant_new = force_constant

    do i = 1,m
        do j = 1,n
            dr = pos(j,:) - pos(i-1+istart,:) 
            if (sqrt(dot_product(dr,dr)).gt.rcut)then
                do icart = 1,3
                    do jcart = 1,3
!                        write(*,*) force_constant_new((i-1)*3+icart,(j-1)*3+jcart)
                        force_constant_new((i-1)*3+icart,(j-1)*3+jcart) = 0.0d0
                    end do
                end do
!            else
!                        if (i.eq.1)then
!                            write(*,*)i,j,force_constant_new((i-1)*3+2,(j-1)*3+2),(j-1)*3+2
!                        end if

            end if
        end do
    end do
    do i = 1,m
        do icart = 1,3
            do jcart = 1,3
                temp_sum = sum(force_constant_new((i-1)*3+icart,jcart:n*3:3))-&
                force_constant_new((i-1)*3+icart,(i-1+istart-1)*3+jcart)
                force_constant_new((i-1)*3+icart,(i-1+istart-1)*3+jcart) = -temp_sum
            end do
        end do
    end do

    force_constant = force_constant_new

    end subroutine asr_cutoff


    ! quantum espresso ifc
    subroutine gen_dyn_qe_mpi(ie,ik_core,norb,dyn,cell,cell_sc,cell_all,pos,pos_all,&
                 idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nxall,nyall,nzall,&
                                      nx_sc,ny_sc,nz_sc,rws,rd,neighborlist,nnbrs)

    use config

    implicit none

    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ie,ik_core
    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(out) :: dyn(numprocs,norb,norb)

    real(kind=8),intent(in) :: cell(3,3),cell_sc(3,3),cell_all(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    real(kind=8),intent(in) :: mass_sc(:)
    real(kind=8),intent(in) :: pos_all(:,:),mass_all(:)
    integer(kind=4),intent(in) :: nxall,nyall,nzall
    integer(kind=4),intent(in) :: idx_sc2uc(:),idxall(:)
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4),intent(in) :: nx_sc,ny_sc,nz_sc
    integer(kind=4),intent(in) :: neighborlist(:,:),nnbrs(:)


    mm = myid+1
    nn = k_start_core(mm)+ik_core-1
    call gen_dyn_qe_k(k_grid(nn,:),dyn(mm,:,:),cell,cell_sc,cell_all,pos,pos_all,&
                idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nxall,nyall,nzall,&
                nx_sc,ny_sc,nz_sc,rws,rd,neighborlist,nnbrs)
 



    end subroutine gen_dyn_qe_mpi

    subroutine gen_dyn_qe_k(kgrid,dyn,cell,cell_sc,cell_all,pos,pos_all,&
             idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nx,ny,nz,&
                     nx_sc,ny_sc,nz_sc,rws,rd,neighborlist,nnbrs)
    use util

    implicit none


    real(kind=8),intent(in) :: kgrid(3)
    real(kind=8),intent(in) :: cell(3,3),cell_sc(3,3),cell_all(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    real(kind=8),intent(in) :: mass_sc(:)
    real(kind=8),intent(in) :: pos_all(:,:),mass_all(:)
    integer(kind=4),intent(in) :: nx,ny,nz
    integer(kind=4),intent(in) :: idx_sc2uc(:),idxall(:)
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4),intent(in) :: nx_sc,ny_sc,nz_sc

    integer(kind=4)            :: natm,natm_sc,natm_all
    real(kind=8)               :: pos_i(3),pos_j(3),dr(3)
    real(kind=8)               :: pos_i_uc(3),pos_j_uc(3)
    real(kind=8)               :: mi,mj
    integer(kind=4)            :: nscx,nscy,nscz
    integer(kind=4) :: nb 
    integer(kind=4) :: i,j,inb,a,b,c,nreq,jj,ir,m,n
    integer(kind=4) :: aa,bb,cc,t1,t2,t3
    real(kind=8)    :: ixyz(3),ivcell(3,3),df,weight
    complex(kind=8) :: dyn(:,:)
    complex(kind=8)    :: vec(size(dyn,1),size(dyn,2))
    real(kind=8)    :: eig(size(dyn,1))
    real(kind=8)    :: kpoint(3),cell_reci(3,3)
    real(kind=8)               :: cell_reci_sc(3,3)
    real(kind=8),allocatable   :: kpath_pt(:,:)
    real(kind=8),allocatable   :: kpath(:,:) 
    integer(kind=4),allocatable :: nkpath(:)
    integer(kind=4) :: npt,itemp,jtemp,npt_total,nline
    integer(kind=4) :: neighborlist(:,:),nnbrs(:)

    nb = size(pos_all,1)*3

    ivcell = inv_real(cell)
    cell_reci = 2*pi*transpose(ivcell)
    cell_reci_sc = 2*pi*transpose(inv_real(cell_all))
    natm = size(pos,1)
    natm_sc = size(pos,1)*nx*ny*nz
    natm_all = size(pos_all,1)

    nscx = ceiling(dble(nx)/dble(nx_sc))
    nscy = ceiling(dble(ny)/dble(ny_sc))
    nscz = ceiling(dble(nz)/dble(nz_sc))

    
        kpoint = matmul(kgrid,cell_reci) 
        dyn(:,:) = 0.0d0
        eig(:) = 0.0d0
        vec(:,:) = 0.0d0
        do i = 1,natm_all
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do inb = 1,nnbrs(i)
                j = neighborlist(i,inb)
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = -nscz,nscz
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        dyn((i-1)*3+m,(j-1)*3+n) = &
                                        dyn((i-1)*3+m,(j-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

!        call eigenH(dyn,eig,vec)

    end subroutine gen_dyn_qe_k

    ! quantum espresso ifc tridiagonal block matrix
    subroutine gen_dyn_qe_tdb_mpi(ie,ik_core,cell,cell_sc,cell_all,pos,pos_all,&
                 idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nxall,nyall,nzall,&
                 nx_sc,ny_sc,nz_sc,rws,rd,neighborlist,nnbrs,&
                 Hamd,Hamdv,Haml,Hamr,Hamlv,Hamrv,&
                 Hamld,Hamrd,nl,nr,nd,nlayer,norb_layer)

    use config

    implicit none

    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ie,ik_core
    integer(kind=4),intent(in) :: nl,nr,nd,nlayer,norb_layer
    complex(kind=8),intent(out) :: Hamd(:,:,:,:),Hamdv(:,:,:,:)
    complex(kind=8),intent(out) :: Haml(:,:,:),Hamr(:,:,:)
    complex(kind=8),intent(out) :: Hamlv(:,:,:),Hamrv(:,:,:)
    complex(kind=8),intent(out) :: Hamld(:,:,:),Hamrd(:,:,:)

    real(kind=8),intent(in) :: cell(3,3),cell_sc(3,3),cell_all(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    real(kind=8),intent(in) :: mass_sc(:)
    real(kind=8),intent(in) :: pos_all(:,:),mass_all(:)
    integer(kind=4),intent(in) :: nxall,nyall,nzall
    integer(kind=4),intent(in) :: idx_sc2uc(:),idxall(:)
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4),intent(in) :: nx_sc,ny_sc,nz_sc
    integer(kind=4),intent(in) :: neighborlist(:,:),nnbrs(:)


    mm = myid+1
    nn = k_start_core(mm)+ik_core-1
    call gen_dyn_qe_tdb_k(k_grid(nn,:),cell,cell_sc,cell_all,pos,pos_all,&
                idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nxall,nyall,nzall,&
                nx_sc,ny_sc,nz_sc,rws,rd,neighborlist,nnbrs,&
                Hamd(mm,:,:,:),Hamdv(mm,:,:,:),Haml(mm,:,:),Hamr(mm,:,:),&
                Hamlv(mm,:,:),Hamrv(mm,:,:),&
                Hamld(mm,:,:),Hamrd(mm,:,:),&
                nlayer,norb_layer)
 



    end subroutine gen_dyn_qe_tdb_mpi

    subroutine gen_dyn_qe_tdb_k(kgrid,cell,cell_sc,cell_all,pos,pos_all,&
             idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nx,ny,nz,&
                     nx_sc,ny_sc,nz_sc,rws,rd,neighborlist,nnbrs,&
                     Hamd,Hamdv,Haml,Hamr,Hamlv,Hamrv,&
                     Hamld,Hamrd,nlayer,norb_layer)
    use util
    use config ,  only : buffer_left, period_left,&
                         buffer_right, period_right

    implicit none


    real(kind=8),intent(in) :: kgrid(3)
    real(kind=8),intent(in) :: cell(3,3),cell_sc(3,3),cell_all(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    real(kind=8),intent(in) :: mass_sc(:)
    real(kind=8),intent(in) :: pos_all(:,:),mass_all(:)
    integer(kind=4),intent(in) :: nx,ny,nz
    integer(kind=4),intent(in) :: idx_sc2uc(:),idxall(:)
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4),intent(in) :: nx_sc,ny_sc,nz_sc

    integer(kind=4)            :: natm,natm_sc,natm_all
    real(kind=8)               :: pos_i(3),pos_j(3),dr(3)
    real(kind=8)               :: pos_i_uc(3),pos_j_uc(3)
    real(kind=8)               :: mi,mj
    integer(kind=4)            :: nscx,nscy,nscz
    integer(kind=4) :: nb 
    integer(kind=4) :: i,j,inb,a,b,c,nreq,jj,ir,m,n
    integer(kind=4) :: device_start,ilay
    integer(kind=4) :: iblk,jblk
    integer(kind=4) :: natm_layer
    integer(kind=4) :: aa,bb,cc,t1,t2,t3
    real(kind=8)    :: ixyz(3),ivcell(3,3),df,weight
    complex(kind=8) :: Hamd(:,:,:),Hamdv(:,:,:),Haml(:,:),Hamr(:,:)
    complex(kind=8) :: Hamlv(:,:),Hamrv(:,:)
    complex(kind=8) :: Hamld(:,:),Hamrd(:,:)
    integer(kind=4),intent(in) :: nlayer,norb_layer
    real(kind=8)    :: kpoint(3),cell_reci(3,3)
    real(kind=8)               :: cell_reci_sc(3,3)
    real(kind=8),allocatable   :: kpath_pt(:,:)
    real(kind=8),allocatable   :: kpath(:,:) 
    integer(kind=4),allocatable :: nkpath(:)
    integer(kind=4) :: npt,itemp,jtemp,npt_total,nline
    integer(kind=4) :: neighborlist(:,:),nnbrs(:)

    Hamd = 0.0d0
    Hamdv = 0.0d0
    Haml = 0.0d0
    Hamr = 0.0d0
    Hamlv = 0.0d0
    Hamrv = 0.0d0
    Hamld = 0.0d0
    Hamrd = 0.0d0

    nb = size(pos_all,1)*3

    ivcell = inv_real(cell)
    cell_reci = 2*pi*transpose(ivcell)
    cell_reci_sc = 2*pi*transpose(inv_real(cell_all))
    natm = size(pos,1)
    natm_sc = size(pos,1)*nx*ny*nz
    natm_all = size(pos_all,1)

    nscx = ceiling(dble(nx)/dble(nx_sc))
    nscy = ceiling(dble(ny)/dble(ny_sc))
    nscz = ceiling(dble(nz)/dble(nz_sc))

    
        kpoint = matmul(kgrid,cell_reci) 
        do i = buffer_left+1,buffer_left+period_left
            iblk = i-buffer_left
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = buffer_left+1,buffer_left+period_left
                jblk = j-buffer_left
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = 0,0
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        haml((iblk-1)*3+m,(jblk-1)*3+n) = &
                                        haml((iblk-1)*3+m,(jblk-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        do i = buffer_left+1,buffer_left+period_left
            iblk = i-buffer_left
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = buffer_left+period_left+1,buffer_left+2*period_left
                jblk = j-buffer_left-period_left
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = 0,0
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        hamlv((iblk-1)*3+m,(jblk-1)*3+n) = &
                                        hamlv((iblk-1)*3+m,(jblk-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        do i = natm_all-buffer_right-period_right+1,natm_all-buffer_right
            iblk = i-(natm_all-buffer_right-period_right)
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = natm_all-buffer_right-period_right+1,natm_all-buffer_right
                jblk = j-(natm_all-buffer_right-period_right)
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = 0,0
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        hamr((iblk-1)*3+m,(jblk-1)*3+n) = &
                                        hamr((iblk-1)*3+m,(jblk-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        do i = natm_all-buffer_right-period_right+1,natm_all-buffer_right
            iblk = i-(natm_all-buffer_right-period_right)
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = natm_all-buffer_right-2*period_right+1,natm_all-buffer_right-period_right
                jblk = j-(natm_all-buffer_right-2*period_right)
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = 0,0
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        hamrv((iblk-1)*3+m,(jblk-1)*3+n) = &
                                        hamrv((iblk-1)*3+m,(jblk-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        device_start = buffer_left+2*period_left
        natm_layer = norb_layer/3
        do ilay = 1,nlayer
            do i = device_start+(ilay-1)*natm_layer+1,&
                   device_start+ilay    *natm_layer
                iblk = i-(device_start+(ilay-1)*natm_layer)
                pos_i(:) = pos_all(i,:)
                pos_i_uc(:) = pos(idxall(i),:)
                mi = mass_all(i)
                do j = device_start+(ilay-1)*natm_layer+1,&
                       device_start+ilay    *natm_layer               
                    jblk = j-(device_start+(ilay-1)*natm_layer)
                    pos_j(:) = pos_all(j,:)
                    pos_j_uc(:) = pos(idxall(j),:)
                    mj = mass_all(j)
                    do a = -nscx,nscx
                        do b = -nscy,nscy
                            do c = 0,0
                                weight = 0.0d0
                                dr = -matmul((/dble(a),dble(b),dble(c)/),&
                                cell_all)+pos_j-pos_i
                                ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                                nreq = 1
                                jj = 0
                                do ir = 1,124
                                    df = dot_product(dr,rws(ir,:))-rd(ir)
                                    if (df.gt.1.0d-5) then
                                        jj = 1
                                        cycle
                                    end if
                                    if  (abs(df).le.1.0d-5) then
                                        nreq = nreq+1
                                    end if
                                end do
                                aa = -idnint(ixyz(1))
                                bb = -idnint(ixyz(2))
                                cc = -idnint(ixyz(3))
                                if (jj .eq. 0) then
                                    weight = 1.0d0/dble(nreq)
                                end if
                                if (weight .gt. 0.0d0) then
                                    t1 = mod(aa+1,nx)
                                    if(t1.le.0) then 
                                        t1=t1+nx
                                    end if
                                    t2 = mod(bb+1,ny)
                                    if(t2.le.0) then
                                        t2=t2+ny
                                    end if
                                    t3 = mod(cc+1,nz)
                                    if(t3.le.0) then 
                                         t3=t3+nz
                                    end if
                                    do m = 1,3
                                        do n = 1,3
                                            hamd(ilay,(iblk-1)*3+m,(jblk-1)*3+n) = &
                                            hamd(ilay,(iblk-1)*3+m,(jblk-1)*3+n) + &
                                            weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                            exp(i_imag*&
                                            dot_product(dr,kpoint))/sqrt(mi*mj)
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do    
                end do
            end do
        end do

        ! V_i_i+1
        do ilay = 1,nlayer-1
            do i = device_start+(ilay-1)*natm_layer+1,&
                   device_start+ilay    *natm_layer
                iblk = i-(device_start+(ilay-1)*natm_layer)
                pos_i(:) = pos_all(i,:)
                pos_i_uc(:) = pos(idxall(i),:)
                mi = mass_all(i)
                do j = device_start+ilay*natm_layer+1,&
                       device_start+(ilay+1)*natm_layer               
                    jblk = j-(device_start+ilay*natm_layer)
                    pos_j(:) = pos_all(j,:)
                    pos_j_uc(:) = pos(idxall(j),:)
                    mj = mass_all(j)
                    do a = -nscx,nscx
                        do b = -nscy,nscy
                            do c = 0,0
                                weight = 0.0d0
                                dr = -matmul((/dble(a),dble(b),dble(c)/),&
                                cell_all)+pos_j-pos_i
                                ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                                nreq = 1
                                jj = 0
                                do ir = 1,124
                                    df = dot_product(dr,rws(ir,:))-rd(ir)
                                    if (df.gt.1.0d-5) then
                                        jj = 1
                                        cycle
                                    end if
                                    if  (abs(df).le.1.0d-5) then
                                        nreq = nreq+1
                                    end if
                                end do
                                aa = -idnint(ixyz(1))
                                bb = -idnint(ixyz(2))
                                cc = -idnint(ixyz(3))
                                if (jj .eq. 0) then
                                    weight = 1.0d0/dble(nreq)
                                end if
                                if (weight .gt. 0.0d0) then
                                    t1 = mod(aa+1,nx)
                                    if(t1.le.0) then 
                                        t1=t1+nx
                                    end if
                                    t2 = mod(bb+1,ny)
                                    if(t2.le.0) then
                                        t2=t2+ny
                                    end if
                                    t3 = mod(cc+1,nz)
                                    if(t3.le.0) then 
                                         t3=t3+nz
                                    end if
                                    do m = 1,3
                                        do n = 1,3
                                            hamdv(ilay,(iblk-1)*3+m,(jblk-1)*3+n) = &
                                            hamdv(ilay,(iblk-1)*3+m,(jblk-1)*3+n) + &
                                            weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                            exp(i_imag*&
                                            dot_product(dr,kpoint))/sqrt(mi*mj)
                                        end do
                                    end do
                                end if
                            end do
                        end do
                    end do    
                end do
            end do
        end do

        ! V_l_d
        do i = buffer_left+  period_left+1,&
               buffer_left+2*period_left
            iblk = i-(buffer_left+period_left)
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = device_start+1,&
                   device_start+natm_layer               
                jblk = j-device_start
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                       do c = 0,0
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        hamld((iblk-1)*3+m,(jblk-1)*3+n) = &
                                        hamld((iblk-1)*3+m,(jblk-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        ! V_r_d
        do i = natm_all-buffer_right-2*period_right+1,natm_all-buffer_right-period_right
            iblk = i-(natm_all-buffer_right-2*period_right)
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = device_start+(nlayer-1)*natm_layer+1,&
                   device_start+nlayer*natm_layer               
                jblk = j-(device_start+(nlayer-1)*natm_layer)
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = 0,0
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        hamrd((iblk-1)*3+m,(jblk-1)*3+n) = &
                                        hamrd((iblk-1)*3+m,(jblk-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

    end subroutine gen_dyn_qe_tdb_k

    subroutine gen_dyn_qe(dyn,cell,cell_sc,cell_all,pos,pos_all,&
             idx_sc2uc,idxall,mass,mass_sc,mass_all,fc,nx,ny,nz,&
                     nx_sc,ny_sc,nz_sc,rws,rd)
    use util

    implicit none

    real(kind=8),intent(in) :: cell(3,3),cell_sc(3,3),cell_all(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    real(kind=8),intent(in) :: mass_sc(:)
    real(kind=8),intent(in) :: pos_all(:,:),mass_all(:)
    integer(kind=4),intent(in) :: nx,ny,nz
    integer(kind=4),intent(in) :: idx_sc2uc(:),idxall(:)
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4),intent(in) :: nx_sc,ny_sc,nz_sc
    integer(kind=4)            :: natm,natm_sc,natm_all
    real(kind=8)               :: pos_i(3),pos_j(3),dr(3)
    real(kind=8)               :: pos_i_uc(3),pos_j_uc(3)
    real(kind=8)               :: mi,mj
    integer(kind=4)            :: nscx,nscy,nscz
    integer(kind=4) :: nb 
    integer(kind=4) :: i,j,a,b,c,nreq,jj,ir,m,n
    integer(kind=4) :: aa,bb,cc,t1,t2,t3
    real(kind=8)    :: ixyz(3),ivcell(3,3),df,weight
    complex(kind=8) :: dyn(:,:)
    complex(kind=8)    :: vec(size(dyn,1),size(dyn,2))
    real(kind=8)    :: eig(size(dyn,1))
    real(kind=8)    :: kpoint(3),cell_reci(3,3)
    real(kind=8),allocatable   :: kpath_pt(:,:)
    real(kind=8),allocatable   :: kpath(:,:) 
    integer(kind=4),allocatable :: nkpath(:)
    integer(kind=4) :: npt,itemp,jtemp,npt_total,nline

    npt = 3
    nb = size(pos_all,1)*3
    allocate(kpath_pt(npt,3),nkpath(3))
    kpath_pt(1,:) = (/0.0d0,0.0d0,0.5d0/)
    kpath_pt(2,:) = (/0.0d0,0.0d0,0.0d0/)
    kpath_pt(3,:) = (/0.5d0,0.5d0,0.0d0/)
    nkpath(1) = 45
    nkpath(2) = 45
    nkpath(3) = 1
    itemp = 1
    npt_total = sum(nkpath)
    allocate(kpath(npt_total,3))
    do i = 1,npt-1
        do j = 1,nkpath(i)
            kpath(itemp,:) = (kpath_pt(i+1,:) - kpath_pt(i,:))*dble(j-1)/dble(nkpath(i)) +&
            kpath_pt(i,:)
            itemp = itemp + 1 
        end do
    end do
    kpath(itemp,:) = kpath_pt(npt,:) 
    open(unit=31,file="ph_qe_all.dat",status="UNKNOWN",action="write")
    write(31,*) nb,npt_total
    nline = ceiling(dble(nb)/6.0d0)

    ivcell = inv_real(cell)
    cell_reci = 2*pi*transpose(ivcell)
    natm = size(pos,1)
    natm_sc = size(pos,1)*nx*ny*nz
    natm_all = size(pos_all,1)

    nscx = ceiling(dble(nx)/dble(nx_sc))
    nscy = ceiling(dble(ny)/dble(ny_sc))
    nscz = ceiling(dble(nz)/dble(nz_sc))

   
    do itemp = 1,npt_total  
        kpoint = matmul(kpath(itemp,:),cell_reci) 
        dyn(:,:) = 0.0d0
        eig(:) = 0.0d0
        vec(:,:) = 0.0d0
        do i = 1,natm_all
            pos_i(:) = pos_all(i,:)
            pos_i_uc(:) = pos(idxall(i),:)
            mi = mass_all(i)
            do j = 1,natm_all
                pos_j(:) = pos_all(j,:)
                pos_j_uc(:) = pos(idxall(j),:)
                mj = mass_all(j)
                do a = -nscx,nscx
                    do b = -nscy,nscy
                        do c = -nscz,nscz
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell_all)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        dyn((i-1)*3+m,(j-1)*3+n) = &
                                        dyn((i-1)*3+m,(j-1)*3+n) + &
                                        weight*fc(m,n,idxall(i),idxall(j),t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        write(31,'(3F10.5)') kpath(itemp,:)
        call eigenH(dyn,eig,vec)
        eig = sqrt(abs(eig))
        if (nline .eq. 1) then
            write(31,'(6F10.5)') eig
        else
            do jtemp = 1,nline-1
                write(31,'(6F10.5)') eig((jtemp-1)*6+1:jtemp*6)
            end do
            write(31,'(6F10.5)') eig((nline-1)*6+1:size(eig))
        end if
    end do
    close(unit=31)

    end subroutine gen_dyn_qe

    subroutine gen_dyn_uc_qe(dyn,cell,cell_sc,pos,pos_sc,&
                        idx_sc2uc,mass,mass_sc,fc,nx,ny,nz,&
                            rws,rd)
    
    use util

    implicit none

    real(kind=8),intent(in) :: cell(3,3),cell_sc(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    real(kind=8),intent(in) :: pos_sc(:,:),mass_sc(:)
    integer(kind=4),intent(in) :: nx,ny,nz
    integer(kind=4),intent(in) :: idx_sc2uc(:)
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4)            :: natm,natm_sc
    real(kind=8)               :: pos_i(3),pos_j(3),dr(3)
    real(kind=8)               :: pos_i_uc(3),pos_j_uc(3)
    real(kind=8)               :: mi,mj
    integer(kind=4)            :: nscx,nscy,nscz
    integer(kind=4) :: nb 
    integer(kind=4) :: i,j,a,b,c,nreq,jj,ir,m,n
    integer(kind=4) :: aa,bb,cc,t1,t2,t3
    real(kind=8)    :: ixyz(3),ivcell(3,3),df,weight
    complex(kind=8) :: dyn(:,:)
    complex(kind=8)    :: vec(size(dyn,1),size(dyn,2))
    real(kind=8)    :: eig(size(dyn,1))
    real(kind=8)    :: kpoint(3),cell_reci(3,3)
    real(kind=8),allocatable   :: kpath_pt(:,:)
    real(kind=8),allocatable   :: kpath(:,:) 
    integer(kind=4),allocatable :: nkpath(:)
    integer(kind=4) :: npt,itemp,jtemp,npt_total,nline

    nb = size(pos,1)*3

    !3
    !0.5 0 0 45 !
    !0 0 0 45 ! G
    !0.5 0.5 0.0 1 !
    npt = 3
    allocate(kpath_pt(npt,3),nkpath(3))
    kpath_pt(1,:) = (/0.5d0,0.0d0,0.0d0/)
    kpath_pt(2,:) = (/0.0d0,0.0d0,0.0d0/)
    kpath_pt(3,:) = (/0.5d0,0.5d0,0.0d0/)
    nkpath(1) = 45
    nkpath(2) = 45
    nkpath(3) = 1
    itemp = 1
    npt_total = sum(nkpath)
    allocate(kpath(npt_total,3))
    do i = 1,npt-1
        do j = 1,nkpath(i)
            kpath(itemp,:) = (kpath_pt(i+1,:) - kpath_pt(i,:))*dble(j-1)/dble(nkpath(i)) +&
            kpath_pt(i,:)
            itemp = itemp + 1 
        end do
    end do
    kpath(itemp,:) = kpath_pt(npt,:) 
    open(unit=11,file="ph_qe.dat",status="UNKNOWN",action="write")
    write(11,*) nb,npt_total
    nline = ceiling(dble(nb)/6.0d0)

    ivcell = inv_real(cell)
    cell_reci = 2*pi*transpose(ivcell)
    natm = size(pos,1)
    natm_sc = size(pos,1)*nx*ny*nz

    do itemp = 1,npt_total  
        kpoint = matmul(kpath(itemp,:),cell_reci) 
        nscx = 1
        nscy = 1
        nscz = 1
        

        dyn(:,:) = 0.0d0
        eig(:) = 0.0d0
        vec(:,:) = 0.0d0
        do i = 1,natm
            pos_i(:) = pos(i,:)
            pos_i_uc(:) = pos(i,:) !pos(idx_sc2uc(i),:)
            mi = mass(i)
            do j = 1,natm
                pos_j(:) = pos(j,:)
                pos_j_uc(:) = pos(j,:) ! pos(idx_sc2uc(j),:)
                mj = mass(j)
                do a = -2*nx,2*nx
                    do b = -2*ny,2*ny
                        do c = -2*nz,2*nz
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        dyn((i-1)*3+m,(j-1)*3+n) = &
                                        dyn((i-1)*3+m,(j-1)*3+n) + &
                                        weight*fc(m,n,i,j,t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        write(11,'(3F10.5)') kpath(itemp,:)
        call eigenH(dyn,eig,vec)
        eig = sqrt(abs(eig))
        if (nline .eq. 1) then
            write(11,'(6F10.5)') eig
        else
            do jtemp = 1,nline-1
                write(11,'(6F10.5)') eig((jtemp-1)*6+1:jtemp*6)
            end do
            write(11,'(6F10.5)') eig((nline-1)*6+1:size(eig))
        end if
    end do
    close(unit=11)
    end subroutine

    subroutine gen_dyn_pc_qe(dyn,cell,pos,mass,fc,nx,ny,nz,&
                            rws,rd)
    
    use util
    use input

    implicit none

    real(kind=8),intent(in) :: cell(3,3)
    real(kind=8),intent(in) :: pos(:,:),mass(:),fc(:,:,:,:,:,:,:)
    integer(kind=4),intent(in) :: nx,ny,nz
    real(kind=8),intent(in)    :: rws(:,:),rd(:)
    integer(kind=4)            :: natm,natm_sc
    real(kind=8)               :: pos_i(3),pos_j(3),dr(3)
    real(kind=8)               :: pos_i_uc(3),pos_j_uc(3)
    real(kind=8)               :: mi,mj
    integer(kind=4)            :: nscx,nscy,nscz
    integer(kind=4) :: nb 
    integer(kind=4) :: i,j,a,b,c,nreq,jj,ir,m,n
    integer(kind=4) :: aa,bb,cc,t1,t2,t3
    real(kind=8)    :: ixyz(3),ivcell(3,3),df,weight
    complex(kind=8) :: dyn(:,:)
    complex(kind=8)    :: vec(size(dyn,1),size(dyn,2))
    real(kind=8)    :: eig(size(dyn,1))
    real(kind=8)    :: kpoint(3),cell_reci(3,3)
    real(kind=8),allocatable   :: kpath_pt(:,:)
    real(kind=8),allocatable   :: kpath(:,:) 
    integer(kind=4),allocatable :: nkpath(:)
    integer(kind=4) :: npt,itemp,jtemp,npt_total,nline

    nb = size(positions_pc_l,1)*3

    !3
    !0.5 0 0 45 !
    !0 0 0 45 ! G
    !0.5 0.5 0.0 1 !
    npt = 3
    allocate(kpath_pt(npt,3),nkpath(3))
    kpath_pt(1,:) = (/0.5d0,0.0d0,0.0d0/)
    kpath_pt(2,:) = (/0.0d0,0.0d0,0.0d0/)
    kpath_pt(3,:) = (/0.5d0,0.5d0,0.0d0/)
    nkpath(1) = 45
    nkpath(2) = 45
    nkpath(3) = 1
    itemp = 1
    npt_total = sum(nkpath)
    allocate(kpath(npt_total,3))
    do i = 1,npt-1
        do j = 1,nkpath(i)
            kpath(itemp,:) = (kpath_pt(i+1,:) - kpath_pt(i,:))*dble(j-1)/dble(nkpath(i)) +&
            kpath_pt(i,:)
            itemp = itemp + 1 
        end do
    end do
    kpath(itemp,:) = kpath_pt(npt,:) 
    open(unit=11,file="pc_ph_qe.dat",status="UNKNOWN",action="write")
    write(11,*) nb,npt_total
    nline = ceiling(dble(nb)/6.0d0)

    ivcell = inv_real(cell)
    cell_reci = recivec_pc_l
    natm = size(pos,1)
    natm_sc = size(pos,1)*nx*ny*nz
    do itemp = 1,npt_total  
        kpoint = matmul(kpath(itemp,:),cell_reci) 
        nscx = 1
        nscy = 1
        nscz = 1
        
        dyn(:,:) = 0.0d0
        eig(:) = 0.0d0
        vec(:,:) = 0.0d0
        do i = 1,natoms_pc_l
            pos_i(:) = pos(i,:)
            pos_i_uc(:) = pos(i,:) !pos(idx_sc2uc(i),:)
            mi = mass(i)
            do j = 1,natm
                pos_j(:) = pos(j,:)
                pos_j_uc(:) = pos(j,:) ! pos(idx_sc2uc(j),:)
                mj = mass(j)
                do a = -2*nx,2*nx
                    do b = -2*ny,2*ny
                        do c = -2*nz,2*nz
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            cell)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,rws(ir,:))-rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,nx)
                                if(t1.le.0) then 
                                    t1=t1+nx
                                end if
                                t2 = mod(bb+1,ny)
                                if(t2.le.0) then
                                    t2=t2+ny
                                end if
                                t3 = mod(cc+1,nz)
                                if(t3.le.0) then 
                                     t3=t3+nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        dyn((i-1)*3+m,mod(j-1,natoms_pc_l)*3+n) = &
                                        dyn((i-1)*3+m,mod(j-1,natoms_pc_l)*3+n) + &
                                        weight*fc(m,n,i,j,t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do

        write(11,'(3F10.5)') kpath(itemp,:)
        dyn = Hermitian(dyn)
        call eigenH(dyn,eig,vec)
        eig = sqrt(abs(eig))
        if (nline .eq. 1) then
            write(11,'(6F10.5)') eig
        else
            do jtemp = 1,nline-1
                write(11,'(6F10.5)') eig((jtemp-1)*6+1:jtemp*6)
            end do
            write(11,'(6F10.5)') eig((nline-1)*6+1:size(eig))
        end if
    end do
    close(unit=11)
    end subroutine

    subroutine gen_dyn_pc_qe_k(kpoint,dyn,mass_uc)
    
    use util
    use input
    use qe

    implicit none

    integer(kind=4)            :: natm,natm_sc
    real(kind=8)               :: pos_i(3),pos_j(3),dr(3)
    real(kind=8)               :: pos_i_uc(3),pos_j_uc(3)
    real(kind=8)               :: mi,mj
    integer(kind=4)            :: nscx,nscy,nscz
    integer(kind=4) :: nb 
    integer(kind=4) :: i,j,a,b,c,nreq,jj,ir,m,n
    integer(kind=4) :: aa,bb,cc,t1,t2,t3
    real(kind=8)    :: ixyz(3),ivcell(3,3),df,weight
    complex(kind=8) :: dyn(:,:)
    complex(kind=8)    :: vec(size(dyn,1),size(dyn,2))
    real(kind=8)    :: eig(size(dyn,1))
    real(kind=8)    :: kpoint(3),cell_reci(3,3)
    integer(kind=4) :: itemp,jtemp,npt_total,nline
    real(kind=8)    :: mass_uc(:)

    nb = size(positions_pc_l,1)*3

    !3
    !0.5 0 0 45 !
    !0 0 0 45 ! G
    !0.5 0.5 0.0 1 !
    itemp = 1

    ivcell = inv_real(qeph_cell)
    cell_reci = recivec_pc_l
    natm = size(qeph_pos,1)
    natm_sc = size(qeph_pos,1)*qeph_nx*qeph_ny*qeph_nz
        nscx = 1
        nscy = 1
        nscz = 1
        

        dyn(:,:) = 0.0d0
        eig(:) = 0.0d0
        vec(:,:) = 0.0d0
        do i = 1,natoms_pc_l
            pos_i(:) = qeph_pos(i,:)
            pos_i_uc(:) = qeph_pos(i,:) !pos(idx_sc2uc(i),:)
            mi = mass_uc(i)
            do j = 1,natm
                pos_j(:) = qeph_pos(j,:)
                pos_j_uc(:) = qeph_pos(j,:) ! pos(idx_sc2uc(j),:)
                mj = mass_uc(j)
                do a = -2*qeph_nx,2*qeph_nx
                    do b = -2*qeph_ny,2*qeph_ny
                        do c = -2*qeph_nz,2*qeph_nz
                            weight = 0.0d0
                            dr = -matmul((/dble(a),dble(b),dble(c)/),&
                            qeph_cell)+pos_j-pos_i
                            ixyz =  matmul(dr-(pos_j_uc-pos_i_uc),ivcell)
                            nreq = 1
                            jj = 0
                            do ir = 1,124
                                df = dot_product(dr,qeph_rws(ir,:))-qeph_rd(ir)
                                if (df.gt.1.0d-5) then
                                    jj = 1
                                    cycle
                                end if
                                if  (abs(df).le.1.0d-5) then
                                    nreq = nreq+1
                                end if
                            end do
                            aa = -idnint(ixyz(1))
                            bb = -idnint(ixyz(2))
                            cc = -idnint(ixyz(3))
                            if (jj .eq. 0) then
                                weight = 1.0d0/dble(nreq)
                            end if
                            if (weight .gt. 0.0d0) then
                                t1 = mod(aa+1,qeph_nx)
                                if(t1.le.0) then 
                                    t1=t1+qeph_nx
                                end if
                                t2 = mod(bb+1,qeph_ny)
                                if(t2.le.0) then
                                    t2=t2+qeph_ny
                                end if
                                t3 = mod(cc+1,qeph_nz)
                                if(t3.le.0) then 
                                     t3=t3+qeph_nz
                                end if
                                do m = 1,3
                                    do n = 1,3
                                        dyn((i-1)*3+m,mod(j-1,natoms_pc_l)*3+n) = &
                                        dyn((i-1)*3+m,mod(j-1,natoms_pc_l)*3+n) + &
                                        weight*qeph_fc(m,n,i,j,t1,t2,t3)* &
                                        exp(i_imag*&
                                        dot_product(dr,kpoint))/sqrt(mi*mj)
                                    end do
                                end do
                            end if
                        end do
                    end do
                end do    
            end do
        end do
        dyn = Hermitian(dyn)

    end subroutine

    subroutine find_neighbors(neighborlist,nnbrs,pos_all,cell_all,natom_all)

    use util
    use config
    
    implicit none

    integer(kind=4) :: neighborlist(:,:)
    integer(kind=4) :: nnbrs(:)
    integer(kind=4) :: natom_all
    real(kind=8) :: cell_all(3,3)   
    real(kind=8) :: pos_all(:,:)
    real(kind=8) :: mindist,rdiff(3),rdd
    integer(kind=4) :: i,j,l,m,n


    
    do i = 1,natom_all
        do j = 1,natom_all
            mindist = 1000.0d0
            do l =-1,1
                do m =-1,1
                    do n =-1,1
                        rdiff(1:3) = pos_all(j,1:3)-pos_all(i,1:3)+&
                        matmul((/dble(l),dble(m),dble(n)/),&
                        cell_all)
                        rdd = sqrt(dot_product(rdiff,rdiff))
                        if (rdd.lt. mindist) then
                            mindist = rdd
                        end if
                    end do
                end do
            end do
            if (mindist.lt.r_cutoff) then
                nnbrs(i) = nnbrs(i)+1
                neighborlist(i,nnbrs(i)) = j
            end if
        end do
    end do
 
    end subroutine find_neighbors



end module dyn
