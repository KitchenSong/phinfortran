module fc

    implicit none

    real(kind=8),allocatable :: force_constant_pc_l(:,:)
    real(kind=8),allocatable :: force_constant_pc_r(:,:)


contains

    subroutine read_fc(force_constant,nuc,nsc,filename)

    implicit none

    
    character(len=*),intent(in) :: filename
    integer(kind=4),intent(in) :: nuc, nsc
    real(kind=8),intent(out) :: force_constant(nuc,nsc)
    
    force_constant = read_binary_python(nuc,nsc,filename)


    end subroutine read_fc

    function read_binary_python(nr,nc,filename) result(D)

    use iso_c_binding
    use iso_fortran_env

    implicit none
    
    integer, parameter :: DPR = selected_real_kind(p=15)
    integer(kind=4),intent(in) :: nr,nc
    character(len=*),intent(in) :: filename
    real(DPR),dimension(nr,nc) :: dat
    real(kind=8),dimension(nr,nc) :: D


    open(20,file=trim(adjustl(filename)),status='old',access='stream',action='read')
    read(20) dat
    D(:,:) = dat(:,:)
    close(20)

    end function read_binary_python

    subroutine extract_fc_pc(fc,n_atm_pc,pos_pc,pos_sc_fc,&
                              latvec_pc,pos_fc_uc,fc_p,side)

    use input
    use config, only : period_left,period_right,&
                       buffer_left,buffer_right,&
                       n_fc_uc_x,n_fc_uc_y,n_fc_uc_z,&
                       n_bloch_x,n_bloch_y
    
    implicit none

    real(kind=8),intent(in) :: fc(:,:)
    real(kind=8),intent(in) :: pos_pc(:,:)
    real(kind=8),intent(out) :: pos_sc_fc(:,:)
    integer(kind=4),intent(in) :: n_atm_pc
    real(kind=8),intent(in) :: latvec_pc(3,3)
    real(kind=8),intent(out) :: fc_p(:,:)
    integer(kind=4) :: n_fc_sc,i,j,k
    integer(kind=4) :: n_fc_uc,icart,jcart
    real(kind=8) :: pos_fc_uc(:,:)
    real(kind=8) :: pos_diff(3)
    integer(kind=4) :: side
    integer(kind=4) :: istart,istart0

    n_fc_uc = size(fc,2)/3
    fc_p = 0.0d0
    call fc_sc_positions(pos_pc,1,1,1,&
         pos_sc_fc,latvec_pc)
    n_fc_sc = size(pos_sc_fc,1)
    fc_p = 0.0d0

    istart = natoms_uc*(n_fc_uc_x*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)&
           + n_fc_uc_y*(2*n_fc_uc_z+1)+n_fc_uc_z)

    if (side .eq. 1) then
         istart0 = (buffer_left)/(n_bloch_x*n_bloch_y)
    else 
         istart0 = natoms_uc -(buffer_right+period_right)/(n_bloch_x*n_bloch_y)
    end if
    do i = 1,n_atm_pc
        do j = 1,n_fc_uc
            do k = 1,n_fc_sc
                pos_diff(1:3) =  (pos_sc_fc(k,1:3)-pos_pc(i,1:3))-(pos_fc_uc(j,1:3)-pos_fc_uc(istart+istart0+i,1:3))
!                write(*,*) pos_pc(i,1:3)- pos_fc_uc(istart+i+istart0,1:3)

                if (sqrt(dot_product(pos_diff,pos_diff)).lt.1.0d-4) then
                    do icart = 1,3
                        do jcart = 1,3
                            fc_p(3*(i-1)+icart,3*(k-1)+jcart) = fc(3*(i+istart0-1)+icart,3*(j-1)+jcart)
                        end do
                    end do
                end if
            end do
        end do
    end do
    end subroutine extract_fc_pc

    subroutine read_qe_fc(cell_qe,cell_sc_qe,pos_qe,pos_sc_qe,&
               idx_sc2uc,idxall,mass_qe,mass_sc_qe,fc_qe,nx_qe,ny_qe,nz_qe,&
               rws,rd,cell,pos,mass_id)

    use config
    use util
    use qe
    use input, only : mass_no

    implicit none

    include "mpif.h"
    
    integer(kind=4)  :: ntype,nat,ibrav
    real(kind=8)     :: celldm(6),cell_qe(3,3),vol,fc
    real(kind=8),allocatable  :: fc_qe(:,:,:,:,:,:,:)
    integer(kind=4)  :: i,j,k,l,itemp,jtemp
    real(kind=8)               :: cell_sc_qe(3,3)
    real(kind=8),allocatable   :: mass_ele_qe(:),mass_qe(:)
    real(kind=8),allocatable   :: mass_sc_qe(:)
    real(kind=8),allocatable   :: pos_qe(:,:),pos_sc_qe(:,:)
    real(kind=8)               :: cell(3,3)
    real(kind=8),allocatable   :: pos(:,:)
    integer(kind=4),allocatable :: mass_id(:)
    integer(kind=4),allocatable :: idx_sc2uc(:),idxall(:)
    character(len=30) :: label
    real(kind=8)    :: bohr2ang = 0.529177
    integer(kind=4) :: nx_qe,ny_qe,nz_qe
    integer(kind=4) :: alpha,beta,iatm,jatm,idx,jdx,kdx
    real(kind=8)    :: rws(124,3),rd(124),dr(3)
    real(kind=8)    :: rmin,rtemp
    integer(kind=4) :: m1,m2,m3
    integer(kind=4) :: norb 


    
    open(23,file=trim(adjustl(flfrc)),status='old',action='read')
    read(23,*) ntype,nat,ibrav,celldm
    do i = 1,3
        read(23,*) cell_qe(i,:)
    end do
    ! number of orbitals
    norb = 3*nat

    allocate(mass_ele_qe(ntype),mass_qe(nat))
    allocate(pos_qe(nat,3))
    allocate(qeph_mass(nat),qeph_pos(nat,3))
    mass_ele_qe(:) = 0.0d0
    mass_qe(:) = 0.0d0
    
    do i = 1,ntype
        read(23,*) itemp, label, vol 
        mass_ele_qe(i) = ele_mass(label)
    end do
    do i = 1,nat
        read(23,*) itemp,jtemp, pos_qe(i,:)
        mass_qe(i) = mass_ele_qe(jtemp) 
    end do
    qeph_mass = mass_qe
    pos_qe(:,:) = pos_qe(:,:) * celldm(1) * bohr2ang
    qeph_pos = pos_qe
    cell_qe(:,:) = cell_qe(:,:) * celldm(1) * bohr2ang
    qeph_cell = cell_qe
    read(23,*) label
    read(23,*) nx_qe,ny_qe,nz_qe
    qeph_nx = nx_qe
    qeph_ny = ny_qe
    qeph_nz = nz_qe

    cell_sc_qe(1,:) = nx_qe*cell_qe(1,:)
    cell_sc_qe(2,:) = ny_qe*cell_qe(2,:)
    cell_sc_qe(3,:) = nz_qe*cell_qe(3,:)
 

    allocate(pos_sc_qe(nat*nx_qe*ny_qe*nz_qe,3))
    allocate(mass_sc_qe(nat*nx_qe*ny_qe*nz_qe))
    allocate(idx_sc2uc(nat*nx_qe*ny_qe*nz_qe))

    pos_sc_qe = 0.0d0
    mass_sc_qe = 0.0d0
    itemp = 1
    do i = 1,nx_qe
        do j = 1,ny_qe
            do k = 1,nz_qe
                do l = 1,nat
                    pos_sc_qe(itemp,:) = pos_qe(l,:)+&
                    matmul((/dble(i-1),dble(j-1),dble(k-1)/),&
                    cell_qe)
                    mass_sc_qe(itemp) = mass_qe(l)
                    idx_sc2uc(itemp) = l
                    itemp = itemp+1
                end do
            end do
        end do
    end do
    allocate(fc_qe(3,3,nat,nat,nx_qe,ny_qe,nz_qe))
    allocate(qeph_fc(3,3,nat,nat,nx_qe,ny_qe,nz_qe))
    do i = 1, (nat*3)**2
        read(23,*) alpha,beta,iatm,jatm
        do j = 1,nx_qe*ny_qe*nz_qe
            read(23,*) idx,jdx,kdx,fc
            fc_qe(alpha,beta,iatm,jatm,idx,jdx,kdx) = fc&
             *eleVolt/1.0d-20/mass_proton&
             /(1.0d12*2.0d0*pi)**2*(13.605662285137/0.529177249**2)

        end do
    end do
    close(23)
    if (update_nbr.eq.1) then
        open(unit=12,file="fc_qe_r.dat",status="UNKNOWN",action="write")

        ! cutoff radius
        do iatm = 1,nat
            do jatm = 1,nat
                do j = 1,nx_qe
                    do k = 1,ny_qe
                        do l = 1,nz_qe
                            rmin = 1.0d15
                            do m1 = -2,2
                                do m2 = -2,2
                                    do m3 = -2,2
                                        dr(:) = pos_qe(jatm,:) - pos_qe(iatm,:) &
                                        - matmul((/dble(j-1),dble(k-1),dble(l-1)/),&
                                        cell_qe)&
                                        + matmul((/dble(m1),dble(m2),dble(m3)/),&
                                        cell_sc_qe)
                                        rtemp = sqrt(dot_product(dr,dr))
                                        if ( rtemp.lt.rmin) then
                                            rmin = rtemp
                                        end if
                                    end do
                                end do
                            end do
                            write(12,'(10F20.5)') rmin,fc_qe(:,:,iatm,jatm,j,k,l)
                            if (rmin .gt. r_cutoff ) then
                                fc_qe(:,:,iatm,jatm,j,k,l) = 0.0d0
                            end if
                        end do
                    end do
                end do
            end do
        end do

        close(12)

    end if

    ! acoustic sum rule
    do i = 1,3
        do j = i,3
            do idx = 1,nat
                fc_qe(i,j,idx,idx,1,1,1) = -sum(fc_qe(i,j,idx,:,:,:,:))&
                +fc_qe(i,j,idx,idx,1,1,1)
            end do
        end do
    end do
    ! Invariance under the permutation of indices
    do i = 1,3
        do j = 1,i
            do idx = 1,nat
               fc_qe(i,j,idx,idx,1,1,1) = fc_qe(j,i,idx,idx,1,1,1)
            end do
        end do
    end do

    qeph_fc = fc_qe

    
    itemp = 1
    do m1 = -2,2
        do m2 = -2,2
            do m3 = -2,2
                if((m1.eq.0).and.(m2.eq.0).and.(m3.eq.0)) then
                    cycle
                end if 
                do i = 1,3
                    rws(itemp,i) = cell_sc_qe(1,i)*dble(m1)&
                                 + cell_sc_qe(2,i)*dble(m2)&
                                 + cell_sc_qe(3,i)*dble(m3)
                end do
                rd(itemp) = 0.5d0*dot_product(rws(itemp,:),rws(itemp,:))
                itemp = itemp+1
            end do
        end do
    end do
    qeph_rws = rws
    qeph_rd = rd
    ! construct super cell
    cell(1,:) = dble(n_bloch_x)*cell_qe(1,:)
    cell(2,:) = dble(n_bloch_y)*cell_qe(2,:)
    cell(3,:) = dble(nz)*cell_qe(3,:)
    allocate(pos(n_bloch_x*n_bloch_y*nz*nat,3))
    allocate(mass_id(n_bloch_x*n_bloch_y*nz*nat))
    allocate(idxall(n_bloch_x*n_bloch_y*nz*nat))
    pos(:,:) = 0.0d0

    ! set up mass profile
    open(24,file='mass_profile.dat',status='old',action='read')
    read(24,*) itemp
    do i = 1,itemp
        read(24,*) mass_id(i)
    end do
    close(24)


    itemp = 1
    do m1 = 1,nz
        do m2 = 1,n_bloch_y
            do m3 = 1,n_bloch_x
                do i = 1,nat
                    pos(itemp,:) = pos_qe(i,:) + &
                     matmul((/dble(m3-1),dble(m2-1),dble(m1-1)/),&
                     cell_qe)
        !            mass(itemp) = mass_qe(i)
                    idxall(itemp) = i
                    itemp = itemp + 1
                end do
            end do
        end do
    end do

 
    end subroutine read_qe_fc

    function ele_mass(ele) result(m)

    implicit none

    character(len=30),intent(in) :: ele
    real(kind=8)    :: m

    if (trim(adjustl(ele)) .eq. 'Si') then
        m = 28.0855
    end if

    end function ele_mass

end module fc
