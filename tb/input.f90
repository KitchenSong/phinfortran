module input

    implicit none

    integer(kind=4) :: natoms, nspecies
    real(kind=8)    :: lscale

    real(kind=8),allocatable :: positions(:,:)
    real(kind=8)             :: latvec(3,3)
    real(kind=8)             :: recivec(3,3)
    real(kind=8)             :: latvec_uc(3,3)
    real(kind=8) :: latvec_uc_l(3,3),latvec_uc_r(3,3)
    real(kind=8) :: reci_uc_l(3,3),reci_uc_r(3,3)
    integer(kind=4),allocatable  :: atm_no(:) 

    integer(kind=4) :: natoms_pc_l, natoms_uc_l,&
                       nspecies_pc_l
    real(kind=8)    :: lscale_pc_l

    real(kind=8),allocatable :: positions_pc_l(:,:),&
                                positions_uc_l(:,:)
    real(kind=8)             :: latvec_pc_l(3,3)
    real(kind=8)             :: recivec_pc_l(3,3)

    integer(kind=4) :: natoms_pc_r, natoms_uc_r,&
                       nspecies_pc_r
    real(kind=8)    :: lscale_pc_r

    real(kind=8),allocatable :: positions_pc_r(:,:),&
                                positions_uc_r(:,:)
    real(kind=8)             :: latvec_pc_r(3,3)
    real(kind=8)             :: recivec_pc_r(3,3)
    real(kind=8),allocatable :: pos_orb_pc_l(:,:)
    real(kind=8),allocatable :: pos_orb_uc_l(:,:)
    real(kind=8),allocatable :: pos_orb_pc_r(:,:)
    real(kind=8),allocatable :: pos_orb_uc_r(:,:)
    integer(kind=4),allocatable :: orb_uc_l_idx(:)
    integer(kind=4),allocatable :: orb_uc_r_idx(:)

    integer(kind=4),allocatable :: sc_to_pc_l(:)
    integer(kind=4),allocatable :: sc_to_pc_r(:)
    real(kind=8) :: rot_misori(3,3)
    real(kind=8) :: rot_z(3,3)
    real(kind=8) :: rot_all(3,3)
    real(kind=8) :: p_misori(1,3)


contains
    
    subroutine load_input()

    use config
    use util

    implicit none

    include "mpif.h"

    character(len=200) :: buffer, label
    integer(kind=4) :: ios, line, pos
    integer(kind=4) :: i, temp_i, temp
    real(kind=8) :: coord_crys_temp(3)
    real(kind=8),allocatable :: pos_temp(:,:) 
    integer(kind=4) :: i1,i2,i3

    ios = 0
    line = 0
        

    ! read tags for natoms, nspecies, lscale
    ! from input.fdf
    open(1,file=trim(adjustl(filename_input)),status="old")
       
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1:pos)   ! name of tag
            buffer = buffer(pos+1:) ! value of tag
            
            select case(label)
            case ("NumberOfAtoms")
                read(buffer, *, iostat=ios) natoms
            case ("NumberOfSpecies")
                read(buffer, *, iostat=ios) nspecies
            case ("LatticeConstant")
                read(buffer(1:scan(buffer,' ')), *, iostat=ios) lscale
            case default
            end select

        end if
    end do
    close(1)

    ! read lattice vectors, positions and 
    ! species of atoms from input.fdf

    allocate(positions(natoms,4))
    allocate(atm_no(nspecies))

    ios = 0

    open(1,file=filename_input,status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec(i,:)
                    latvec(i,:) = latvec(i,:) * lscale
                end do
            else if (index(buffer, "%block AtomicCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms
                    read(1,*,iostat=ios) positions(i,:)
                    positions(i,1:3) = positions(i,1:3) * lscale
                end do
            else if (index(buffer, "%block ChemicalSpeciesLabel")) then
                do i = 1,nspecies
                    read(1,*,iostat=ios) temp_i, temp
                    atm_no(temp_i) = temp ! the atomic number of each element
                                          ! species
                end do    
            end if    
        end if
    end do
    close(1)

    ! compute the reciprocal lattice vector
    ! in unit of ang-1
    recivec = 2*pi*transpose(inv_real(latvec))

    
    ! if doing unfolding, we need to know the 
    ! size of transverse unitcell
    latvec_uc = latvec
    latvec_uc(1,:) = latvec(1,:)/dble(n_bloch_x)
    latvec_uc(2,:) = latvec(2,:)/dble(n_bloch_y)

    latvec_uc_l = latvec_uc
    latvec_uc_r = latvec_uc
    
    ! read positions of atoms in primitive cell
    ! from left_primitive_cell_input.fdf
    
    ios = 0

    open(1,file=trim(adjustl(left_primitive_cell_input)),status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1:pos)   ! name of tag
            buffer = buffer(pos+1:) ! value of tag
            
            select case(label)
            case ("NumberOfAtoms")
                read(buffer, *, iostat=ios) natoms_uc_l
            case ("NumberOfSpecies")
                read(buffer, *, iostat=ios) nspecies_pc_l
            case ("LatticeConstant")
                read(buffer(1:scan(buffer,' ')), *, iostat=ios) lscale_pc_l
            case default
            end select

        end if
    end do

    close(1)

    allocate(positions_uc_l(natoms_uc_l,4))

    ios = 0
    theta = theta*pi/180
    rot_z(1,:) = (/sqrt(2.0)/2.0,-sqrt(2.0)/2.0, 0.0/)
    rot_z(2,:) = (/sqrt(2.0)/2.0, sqrt(2.0)/2.0, 0.0/)
    rot_z(3,:) = (/      0.0,       0.0, 1.0/)
    rot_misori(1,:) = (/1.0,0.0,0.0/)
    rot_misori(2,:) = (/0.0d0,dcos(theta),-dsin(theta)/)
    rot_misori(3,:) = (/0.0d0,dsin(theta),dcos(theta)/)
    rot_all = rot_z
    rot_all = matmul(rot_misori,rot_all)
    rot_all = matmul(transpose(rot_z),rot_all)

    open(1,file=trim(adjustl(left_primitive_cell_input)),status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec_pc_l(i,:)
                    latvec_pc_l(i,:) = latvec_pc_l(i,:) * lscale_pc_l
                    p_misori(1,:) = latvec_pc_l(i,:)
                    p_misori = transpose(matmul(rot_all,transpose(p_misori)))
                    latvec_pc_l(i,:) = p_misori(1,:)
                end do
            else if (index(buffer, "%block AtomicScaledCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc_l
                    read(1,*,iostat=ios) positions_uc_l(i,:)
                    positions_uc_l(i,1:3) = matmul(positions_uc_l(i,1:3),latvec_pc_l)
                end do
            end if    
        end if
    end do
    close(1)

    recivec_pc_l = inv_real(latvec_pc_l)
    natoms_pc_l = 0
    allocate(positions_pc_l(1,1),pos_temp(1,1))
    do i = 1, natoms_uc_l
        coord_crys_temp(:) = matmul(positions_uc_l(i,1:3),recivec_pc_l)
        if ((coord_crys_temp(1) .lt. 0.999).and.(coord_crys_temp(1) .gt. -0.001).and.&
            (coord_crys_temp(2) .lt. 0.999).and.(coord_crys_temp(2) .gt. -0.001).and.&
            (coord_crys_temp(3) .lt. 0.999).and.(coord_crys_temp(3) .gt. -0.001)) then
            natoms_pc_l = natoms_pc_l + 1
            deallocate(positions_pc_l)
            allocate(positions_pc_l(natoms_pc_l,4))
            if (natoms_pc_l.eq.1) then
                deallocate(pos_temp)
                allocate(pos_temp(1,4))
                positions_pc_l(natoms_pc_l,:) = positions_uc_l(i,:)
                pos_temp = positions_pc_l
            else
                positions_pc_l(1:natoms_pc_l-1,:) = pos_temp
                positions_pc_l(natoms_pc_l,:) = positions_uc_l(i,:)
                deallocate(pos_temp)
                allocate(pos_temp(natoms_pc_l,4))
                pos_temp = positions_pc_l
            end if
        end if
    end do

    recivec_pc_l = 2*pi*transpose(recivec_pc_l)

    ! read positions of atoms in primitive cell
    ! from right_primitive_cell_input.fdf
    
    ios = 0

    open(1,file=trim(adjustl(right_primitive_cell_input)),status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1:pos)   ! name of tag
            buffer = buffer(pos+1:) ! value of tag
            
            select case(label)
            case ("NumberOfAtoms")
                read(buffer, *, iostat=ios) natoms_uc_r
            case ("NumberOfSpecies")
                read(buffer, *, iostat=ios) nspecies_pc_r
            case ("LatticeConstant")
                read(buffer(1:scan(buffer,' ')), *, iostat=ios) lscale_pc_r
            case default
            end select

        end if
    end do

    close(1)

    allocate(positions_uc_r(natoms_uc_r,4))

    ios = 0

    rot_all = rot_z
    rot_all = matmul(transpose(rot_misori),rot_all)
    rot_all = matmul(transpose(rot_z),rot_all)

    open(1,file=trim(adjustl(right_primitive_cell_input)),status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec_pc_r(i,:)
                    latvec_pc_r(i,:) = latvec_pc_r(i,:) * lscale_pc_r
                    p_misori(1,:) = latvec_pc_r(i,:)
                    p_misori = transpose(matmul(rot_all,transpose(p_misori)))
                    latvec_pc_r(i,:) = p_misori(1,:)
                end do
            else if (index(buffer, "%block AtomicScaledCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc_r
                    read(1,*,iostat=ios) positions_uc_r(i,:)
                    positions_uc_r(i,1:3) = matmul(positions_uc_r(i,1:3),latvec_pc_r)
                end do
            end if    
        end if
    end do
    close(1)

    recivec_pc_r = inv_real(latvec_pc_r)
    natoms_pc_r = 0
    deallocate(pos_temp)
    allocate(positions_pc_r(1,1),pos_temp(1,1))
    do i = 1, natoms_uc_r
        coord_crys_temp(:) = matmul(positions_uc_r(i,1:3),recivec_pc_r)
        if ((coord_crys_temp(1) .lt. 0.999).and.(coord_crys_temp(1) .gt. -0.001).and.&
            (coord_crys_temp(2) .lt. 0.999).and.(coord_crys_temp(2) .gt. -0.001).and.&
            (coord_crys_temp(3) .lt. 0.999).and.(coord_crys_temp(3) .gt. -0.001)) then
            natoms_pc_r = natoms_pc_r + 1
            deallocate(positions_pc_r)
            allocate(positions_pc_r(natoms_pc_r,4))
            if (natoms_pc_r.eq.1) then
                deallocate(pos_temp)
                allocate(pos_temp(1,4))
                positions_pc_r(natoms_pc_r,:) = positions_uc_r(i,:)
                pos_temp = positions_pc_r
            else
                positions_pc_r(1:natoms_pc_r-1,:) = pos_temp
                positions_pc_r(natoms_pc_r,:) = positions_uc_r(i,:)
                deallocate(pos_temp)
                allocate(pos_temp(natoms_pc_r,4))
                pos_temp = positions_pc_r
            end if
        end if
    end do

    recivec_pc_r = 2*pi*transpose(recivec_pc_r)

    ! Since we now have the reciprocal lattice vector of primitive cell,
    ! we can translate the reciprocal lattice site of unit cell inside 
    ! the Wigner-Seitz brillouin zone.
!        do i1 =-3,3
!            do i2 =-3,3
!                do i3 =-3,3
!                    k3 = k_grid(i,:) + np.dot([i1,i2,i3],reci_sc)
!                    temp = in_ws(a_vectors,k_temp)
!                    if temp == 1:
!                        print i, np.dot(k_temp,PC),np.linalg.norm(k_temp)
!                end do
!            end do
!        end do

    end subroutine load_input 

    subroutine k_grid_in_ws()

    use config
    use util

    implicit none

    real(kind=8) :: a_vectors(3,3), k_pos(3)
    real(kind=8) :: k_all(25,3)
    real(kind=8) :: dist(25)
    integer(kind=4) :: i_min
    integer(kind=4) :: icnt,i,i1,i2
    real(kind=8) :: ndiff(3)
    integer(kind=4) :: output

    dist = 0.0d0
    k_all = 0.0d0
    a_vectors = recivec
    a_vectors(3,:) = 0.0d0
    do i = 1,nk 
        k_pos = matmul(k_grid(i,:),a_vectors)
        dist = 0
        icnt = 0
        do i1 =-2,2
            do i2 =-2,2
                ndiff(1) = dble(i1)
                ndiff(2) = dble(i2)
                ndiff(3) = 0.0d0

                icnt = icnt+ 1
                
                k_all(icnt,:) = ndiff 
                dist(icnt) = dot_product(k_pos-matmul(ndiff,a_vectors),&
                                         k_pos-matmul(ndiff,a_vectors)) 
            end do
        end do
    
        i_min = minloc(dist,1)
        k_grid(i,:) = k_grid(i,:) - k_all(i_min,:) ! shift it inside FBZ
    
    end do
    
    end subroutine k_grid_in_ws

end module input
