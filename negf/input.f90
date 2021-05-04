module input

    implicit none

    integer(kind=4) :: natoms, nspecies, natoms_sc
    integer(kind=4) :: natoms_uc
    real(kind=8)    :: lscale

    real(kind=8),allocatable :: positions(:,:)
    real(kind=8),allocatable :: positions_sc(:,:)
    real(kind=8),allocatable :: positions_uc(:,:)
    real(kind=8)             :: latvec(3,3)
    real(kind=8)             :: latvec_left(3,3)
    real(kind=8)             :: latvec_right(3,3)
    real(kind=8)             :: latvec_sc(3,3)
    real(kind=8)             :: recivec(3,3)
    real(kind=8)             :: latvec_uc(3,3)
    real(kind=8) :: latvec_uc_l(3,3),latvec_uc_r(3,3)
    real(kind=8) :: latvec_uc_l_sc(3,3),latvec_uc_r_sc(3,3)
    real(kind=8) :: reci_uc_l(3,3),reci_uc_r(3,3)
    integer(kind=4),allocatable  :: atm_no(:)
    real(kind=8),allocatable ::     mass_no(:) 

    integer(kind=4) :: natoms_pc_l, natoms_uc_l,&
                       nspecies_pc_l
    real(kind=8)    :: lscale_pc_l

    real(kind=8),allocatable :: positions_pc_l(:,:),&
                                positions_uc_l(:,:),&
                                pos_pc_l_fc_sc(:,:)
    real(kind=8)             :: latvec_pc_l(3,3)
    real(kind=8)             :: recivec_pc_l(3,3)

    integer(kind=4) :: natoms_pc_r, natoms_uc_r,&
                       nspecies_pc_r
    real(kind=8)    :: lscale_pc_r

    real(kind=8),allocatable :: positions_pc_r(:,:),&
                                positions_uc_r(:,:),&
                                pos_pc_r_fc_sc(:,:)
    real(kind=8)             :: latvec_pc_r(3,3)
    real(kind=8)             :: recivec_pc_r(3,3)

    real(kind=8),allocatable :: pos_orb_pc_l(:,:)
    real(kind=8),allocatable :: pos_orb_uc_l(:,:)
    real(kind=8),allocatable :: pos_orb_pc_r(:,:)
    real(kind=8),allocatable :: pos_orb_uc_r(:,:)
    integer(kind=4),allocatable :: orb_uc_l_idx(:)
    integer(kind=4),allocatable :: orb_uc_r_idx(:)

    integer(kind=4) :: norb_pc_l,norb_pc_r
    real(kind=8),allocatable :: mass_uc_l(:),mass_uc_r(:)


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
    allocate(mass_no(nspecies))

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
                latvec(:,1:3) = rotate_z(latvec(:,1:3),pi/4.0d0)
            else if (index(buffer, "%block AtomicCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms
                    read(1,*,iostat=ios) positions(i,:)
                    positions(i,1:3) = positions(i,1:3) * lscale
                end do
                ! rotate around z axis
                positions(:,1:3) = rotate_z(positions(:,1:3),pi/4.0d0)

            else if (index(buffer, "%block ChemicalSpeciesLabel")) then
                do i = 1,nspecies
                    read(1,*,iostat=ios) temp_i, temp
                    atm_no(temp_i) = temp ! the atomic number of each element
                                          ! species
                    mass_no(temp_i) = element_mass(temp)
                end do    
            end if    
        end if
    end do
    close(1)

    ! compute the reciprocal lattice vector
    ! in unit of ang-1
    recivec = transpose(2*pi*inv_real(latvec))
    latvec_left = latvec
    latvec_right = latvec


    ! read from unit cell
    ! read tags for natoms, nspecies, lscale
    ! from input_uc.fdf
    ios = 0

    open(1,file=trim(adjustl(filename_input_uc)),status="old")
       
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1:pos)   ! name of tag
            buffer = buffer(pos+1:) ! value of tag
            
            select case(label)
            case ("NumberOfAtoms")
                read(buffer, *, iostat=ios) natoms_uc
            case ("LatticeConstant")
                read(buffer(1:scan(buffer,' ')), *, iostat=ios) lscale
            case default
            end select

        end if
    end do
    close(1)

   
    ! read lattice vectors, positions and 
    ! species of atoms from input.fdf

    allocate(positions_uc(natoms_uc,4))

    ios = 0

    open(1,file=filename_input_uc,status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec_uc(i,:)
                    latvec_uc(i,:) = latvec_uc(i,:) * lscale
                end do
                latvec_uc(:,1:3) = rotate_z(latvec_uc(:,1:3),pi/4.0d0)
            else if (index(buffer, "%block AtomicCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc
                    read(1,*,iostat=ios) positions_uc(i,:)
                    positions_uc(i,1:3) = positions_uc(i,1:3) * lscale
                end do
                ! rotate around z axis
                positions_uc(:,1:3) = rotate_z(positions_uc(:,1:3),pi/4.0d0)
            end if    
        end if
    end do
    close(1)

    ! dimensions for the supercell in molecular 
    ! dynamics calculation
    latvec_sc(1,:) = latvec_uc(1,:) * dble(2*n_fc_uc_x+1)
    latvec_sc(2,:) = latvec_uc(2,:) * dble(2*n_fc_uc_y+1)
    latvec_sc(3,:) = latvec_uc(3,:) * dble(2*n_fc_uc_z+1)
    natoms_sc = (2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)*natoms_uc
    allocate(positions_sc(natoms_sc,4))
    call fc_sc_positions(positions_uc,n_fc_uc_x,n_fc_uc_y,n_fc_uc_z,positions_sc,latvec_uc)

    ! dimensions for lead
    latvec_uc_l = latvec_uc
    latvec_uc_r = latvec_uc

    latvec_uc_l_sc(1,:) = latvec_uc_l(1,:) * dble(2*n_fc_uc_x+1)
    latvec_uc_l_sc(2,:) = latvec_uc_l(2,:) * dble(2*n_fc_uc_y+1)
    latvec_uc_l_sc(3,:) = 0.0d0
    latvec_uc_r_sc(1,:) = latvec_uc_r(1,:) * dble(2*n_fc_uc_x+1)
    latvec_uc_r_sc(2,:) = latvec_uc_r(2,:) * dble(2*n_fc_uc_y+1)
    latvec_uc_r_sc(3,:) = 0.0d0


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
    allocate(mass_uc_l(natoms_uc_l))

    ios = 0

    open(1,file=trim(adjustl(left_primitive_cell_input)),status="old")
    do while (ios == 0)
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec_pc_l(i,:)
                    latvec_pc_l(i,:) = latvec_pc_l(i,:) * lscale_pc_l
                end do
            else if (index(buffer, "%block AtomicScaledCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc_l
                    read(1,*,iostat=ios) positions_uc_l(i,:)
                    positions_uc_l(i,1:3) = matmul(positions_uc_l(i,1:3),latvec_pc_l)
                    mass_uc_l(i) = mass_no(int(positions_uc_l(i,4)))
                end do
            end if
        end if
    end do
    close(1)

    latvec_pc_l(:,1:3) = rotate_z(latvec_pc_l(:,1:3),pi/4.0d0)
    positions_uc_l(:,1:3) = rotate_z(positions_uc_l(:,1:3),pi/4.0d0)

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
    norb_pc_l = 3* natoms_pc_l
    allocate(pos_pc_l_fc_sc(27*natoms_pc_l,4))
    pos_pc_l_fc_sc(:,:) = 0.0d0

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
    allocate(mass_uc_r(natoms_uc_r))

    ios = 0

    open(1,file=trim(adjustl(right_primitive_cell_input)),status="old")
    do while (ios == 0)
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec_pc_r(i,:)
                    latvec_pc_r(i,:) = latvec_pc_r(i,:) * lscale_pc_r
                end do
            else if (index(buffer, "%block AtomicScaledCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc_r
                    read(1,*,iostat=ios) positions_uc_r(i,:)
                    positions_uc_r(i,1:3) = matmul(positions_uc_r(i,1:3),latvec_pc_r)
                    mass_uc_r(i) = mass_no(int(positions_uc_r(i,4)))
                end do
            end if
        end if
    end do
    close(1)

    latvec_pc_r(:,1:3) = rotate_z(latvec_pc_r(:,1:3),pi/4.0d0)
    positions_uc_r(:,1:3) = rotate_z(positions_uc_r(:,1:3),pi/4.0d0)

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


    norb_pc_r = 3* natoms_pc_r
    allocate(pos_pc_r_fc_sc(27*natoms_pc_r,4))
    pos_pc_r_fc_sc(:,:) = 0.0d0


    recivec_pc_r = 2*pi*transpose(recivec_pc_r)

    end subroutine load_input

    function element_mass(no) result(mass)

    use config

    real(kind=8) :: mass
    integer(kind=4) :: no

    if (mass_only .eq. 0) then
        select case(no)
        case (14) ! Si
            mass = 28.0855
        case (32) ! Ge
            mass = 72.64 ! 33.1 !72.64
        case(33)
            mass = 32.0971
        case(34)
            mass = 31.0942
        case(35)
            mass = 30.0913
        case(36)
            mass = 29.0884
        end select
    else
        select case(no)
        case (14) ! Si
        mass = 28.0855
        case (32) ! Ge
        mass = 28.0855*mass_ratio
        end select
    end if
    end function element_mass

    subroutine fc_sc_positions(pos,nx,ny,nz,pos_sc,cell)

    implicit none

    integer(kind=4) :: nx,ny,nz,natoms,natoms_sc
    integer(kind=4) :: nucs
    integer(kind=4) :: i,j,k,n,idx,l
    real(kind=8),intent(in) :: pos(:,:)
    real(kind=8) :: cell(3,3)
    real(kind=8),intent(out) :: pos_sc(:,:) 

    pos_sc = 0.0d0

    natoms = size(pos,1)
    natoms_sc = natoms*(2*nx+1)*(2*ny+1)*(2*nz+1)

    do i = -nx,nx
        do j = -ny,ny
            do k = -nz,nz
                do n = 1,natoms
                    idx = ((i+nx)*(2*ny+1)*(2*nz+1)&
                          +(j+ny)*(2*nz+1)+nz+k)*natoms+n
                    pos_sc(idx,1:3) = pos(n,1:3)+&
                    matmul((/dble(i),dble(j),dble(k)/),&
                           cell)
                    pos_sc(idx,4) = pos(n,4)
                end do
            end do
        end do
    end do

    end subroutine fc_sc_positions
    
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

    subroutine read_pc_qe()

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
    integer(kind=4) :: mtemp

    ios = 0

    open(1,file=trim(adjustl(filename_input)),status="old")
       
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            line = line + 1
            pos = scan(buffer, ' ')
            label = buffer(1:pos)   ! name of tag
            buffer = buffer(pos+1:) ! value of tag
            
            select case(label)
            case ("NumberOfSpecies")
                read(buffer, *, iostat=ios) nspecies
            case default
            end select

        end if
    end do
    close(1)

   
    ! read lattice vectors, positions and 
    ! species of atoms from input.fdf
    allocate(mass_no(nspecies))

    ios = 0

    open(1,file=filename_input,status="old")
    do while (ios == 0) 
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block ChemicalSpeciesLabel")) then
                do i = 1,nspecies
                    read(1,*,iostat=ios) temp_i, temp
                    mass_no(temp_i) = element_mass(temp)
                end do    
            end if    
        end if
    end do
    close(1)

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
   
    allocate(positions_uc_l(natoms_uc_l,3))
    allocate(mass_uc_l(natoms_uc_l))

    ios = 0

    open(1,file=trim(adjustl(left_primitive_cell_input)),status="old")
    do while (ios == 0)
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios) latvec_pc_l(i,:)  
                    latvec_pc_l(i,:) = latvec_pc_l(i,:) * lscale_pc_l
                end do
                latvec_pc_l(:,1:3) = rotate_z(latvec_pc_l(:,1:3),pi/4.0d0)
            else if (index(buffer, "%block AtomicScaledCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc_l
                    read(1,*,iostat=ios) positions_uc_l(i,:),mtemp
                    positions_uc_l(i,1:3) = matmul(positions_uc_l(i,1:3),latvec_pc_l)
                    mass_uc_l(i) = mass_no(mtemp) 
                end do

            end if
        end if
    end do
    close(1)

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


    allocate(positions_uc_r(natoms_uc_r,3))
    allocate(mass_uc_r(natoms_uc_r))

    ios = 0

    open(1,file=trim(adjustl(right_primitive_cell_input)),status="old")
    do while (ios == 0)
        read(1, '(A)', iostat=ios) buffer
        if (ios == 0) then
            if (index(buffer, "%block LatticeVectors"))  then
                do i = 1, 3
                    read(1,*,iostat=ios)latvec_pc_r(i,:)
                    latvec_pc_r(i,:) = latvec_pc_r(i,:) * lscale_pc_r
                end do
                latvec_pc_r(:,1:3) = rotate_z(latvec_pc_r(:,1:3),pi/4.0d0)
            else if (index(buffer, "%block AtomicScaledCoordinatesAndAtomicSpecies")) then
                do i = 1,natoms_uc_r
                    read(1,*,iostat=ios) positions_uc_r(i,:),mtemp
                    positions_uc_r(i,1:3) = matmul(positions_uc_r(i,1:3),latvec_pc_r)
                    mass_uc_r(i) = mass_no(mtemp) 
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
            allocate(positions_pc_l(natoms_pc_l,3))
            if (natoms_pc_l.eq.1) then
                deallocate(pos_temp)
                allocate(pos_temp(1,3))
                positions_pc_l(natoms_pc_l,:) = positions_uc_l(i,:)
                pos_temp = positions_pc_l
            else
                positions_pc_l(1:natoms_pc_l-1,:) = pos_temp
                positions_pc_l(natoms_pc_l,:) = positions_uc_l(i,:)
                deallocate(pos_temp)
                allocate(pos_temp(natoms_pc_l,3))
                pos_temp = positions_pc_l
            end if
        end if
    end do
    recivec_pc_l = 2*pi*transpose(recivec_pc_l)
    norb_pc_l = 3* natoms_pc_l

    
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
            allocate(positions_pc_r(natoms_pc_r,3))
            if (natoms_pc_r.eq.1) then
                deallocate(pos_temp)
                allocate(pos_temp(1,3))
                positions_pc_r(natoms_pc_r,:) = positions_uc_r(i,:)
                pos_temp = positions_pc_r
            else
                positions_pc_r(1:natoms_pc_r-1,:) = pos_temp
                positions_pc_r(natoms_pc_r,:) = positions_uc_r(i,:)
                deallocate(pos_temp)
                allocate(pos_temp(natoms_pc_r,3))
                pos_temp = positions_pc_r
            end if
        end if
    end do
    recivec_pc_r = 2*pi*transpose(recivec_pc_r)
    norb_pc_r = 3* natoms_pc_r

    end subroutine read_pc_qe
end module input
