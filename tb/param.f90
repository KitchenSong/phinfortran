module param
    implicit none
    character(len=*), parameter :: filename_param="param"
    real(kind=8),allocatable :: V(:,:,:,:,:)  ! Harrison parameters
    real(kind=8),allocatable :: Vn(:,:,:,:,:)  ! Harrison parameters' power 
    integer(kind=4),allocatable :: orb_num(:) ! total number of orbitals 
    integer(kind=4),allocatable :: orb_n(:)   ! number of principal orbitals
    integer(kind=4),allocatable :: hop_num(:,:)  ! number of hopping parameters
    real(kind=8),allocatable :: r_d(:,:)
    character(len=8),allocatable :: atom_orb_list(:,:)
    integer(kind=4),allocatable :: atom_orb_num(:)
    integer(kind=4),allocatable :: atom_orb_num_pc_l(:)
    integer(kind=4),allocatable :: atom_orb_num_uc_l(:)
    integer(kind=4),allocatable :: atom_orb_num_pc_r(:)
    integer(kind=4),allocatable :: atom_orb_num_uc_r(:)

contains
    subroutine load_param()

    use input 
    
    implicit none

    include "mpif.h"

    integer(kind=4) :: ios, temp(3), tmp(2)
    integer(kind=4) :: i, j, line, n
    real(kind=8)    :: tempp(5)
    character(len=300) :: buffer

   
    ! read the number of orbitals for each atom
    ios = 0
    allocate(orb_num(nspecies))
    allocate(orb_n(nspecies)) 
    orb_num = 0
    orb_n = 0
    open(1,file=filename_param,status="old")
    read(1, '(A)', iostat=ios) buffer
    do while ((ios == 0).and.(.not.(index(buffer, "%endonsite"))))
        read(1, '(A)', iostat=ios) buffer
        if (index(buffer, "%block")) then 
            read(1,*, iostat=ios) temp 
            orb_n(temp(1)) = temp(2)   ! number of principal orbitals
            orb_num(temp(1)) = temp(3) ! number of orbitals that an atom has
        end if
    end do
    close(1)
    
    ! read the orbitals for each atom
    ios = 0
    allocate(atom_orb_list(nspecies,int(maxval(orb_num))))
    atom_orb_list(:,:) = 'none'
    open(1,file=filename_param,status="old")
    read(1, '(A)', iostat=ios) buffer
    n = 1
    do while ((ios == 0).and.(.not.(index(buffer, "%endonsite"))))
        read(1, '(A)', iostat=ios) buffer
        if (index(buffer, "%block")) then 
            read(1, '(A)', iostat=ios) buffer
            read(1,*,iostat=ios) atom_orb_list(n,:)
            n = n+1
        end if
    end do
    close(1)

    
    ! read the on site parameters
    ios = 0 
    allocate(V(nspecies,nspecies,4,4,4)) ! consider up to d orbitals.
                                         ! consider on-site, sigma, pi and delta bond
                                         ! This can be improved
    V = 0.0d0
    allocate(Vn(nspecies,nspecies,4,4,4)) ! the power law of bond length
    Vn = 0.0d0

    open(1,file=filename_param,status="old")
    read(1, '(A)', iostat=ios) buffer
    read(1, '(A)', iostat=ios) buffer
    if (index(buffer, "%block")) then
        do i = 1, nspecies
            read(1, '(A)', iostat=ios) buffer
            read(1, '(A)', iostat=ios) buffer
            do j = 1, orb_n(i)
                read(1,*, iostat=ios) tempp
                if (int(tempp(1)).eq.5) then
                    V(i,i,4,4,int(tempp(3)+1)) =&
                tempp(4)
                     Vn(i,i,4,4,int(tempp(3)+1)) =&
                tempp(5)
                else
                    V(i,i,int(tempp(1)+1),int(tempp(2)+1),int(tempp(3)+1)) =&
                tempp(4)
                    Vn(i,i,int(tempp(1)+1),int(tempp(2)+1),int(tempp(3)+1)) =&
                tempp(5)
                end if
            end do 
            read(1, '(A)', iostat=ios) buffer
            read(1, '(A)', iostat=ios) buffer
        end do
    end if
    close(1)
        
    ! determine the number of hopping parameters
    ios = 0
    allocate(hop_num(nspecies,nspecies))
    hop_num = 0
    open(1,file=filename_param,status="old")
    read(1, '(A)', iostat=ios) buffer
    do while ((ios == 0).and.(.not.(index(buffer,"%endonsite"))))
        read(1, '(A)', iostat=ios) buffer
    end do
    read(1, '(A)', iostat=ios) buffer
    do while ((ios == 0).and.(.not.(index(buffer, "%endhopping"))))
        read(1, '(A)', iostat=ios) buffer
        if (index(buffer, "%block")) then 
            read(1, *, iostat=ios) tmp
            line = 0
            do while ((ios == 0).and.(.not.(index(buffer, "%endblock"))))
                read(1, '(A)', iostat=ios) buffer
                line = line + 1
            end do
            hop_num(tmp(1)+1,tmp(2)+1) = line-1
        end if
    end do
    close(1)
   
    ! read hopping parameters and r_d
    ios = 0
    open(1,file=filename_param,status="old")
    do while ((ios == 0).and.(.not.(index(buffer,"%endonsite"))))
        read(1, '(A)', iostat=ios) buffer
    end do
    read(1, '(A)', iostat=ios) buffer
    do i = 1, nspecies**2
        read(1, '(A)', iostat=ios) buffer
        read(1, *, iostat=ios) tmp
        do j = 1, hop_num(tmp(1)+1,tmp(2)+1)
            read(1, *, iostat=ios) tempp
            if ((int(tempp(1)).eq.5).and.(int(tempp(2)).ne.5)) then
                V(tmp(1)+1,tmp(2)+1,4,int(tempp(2)+1),int(tempp(3)+2))&
            = tempp(4)
            else if ((int(tempp(1)).ne.5).and.(int(tempp(2)).eq.5)) then
                V(tmp(1)+1,tmp(2)+1,int(tempp(1)+1),4,int(tempp(3)+2))&
            = tempp(4)
            else if ((int(tempp(1)).eq.5).and.(int(tempp(2)).eq.5)) then
                V(tmp(1)+1,tmp(2)+1,4,4,int(tempp(3)+2))&
            = tempp(4)
            else
                V(tmp(1)+1,tmp(2)+1,int(tempp(1)+1),int(tempp(2)+1),int(tempp(3)+2))&
            = tempp(4)               
            end if
        end do
        read(1, '(A)', iostat=ios) buffer
    end do
     
    ! read r_d 
    allocate(r_d(nspecies,nspecies))
    read(1, '(A)', iostat=ios) buffer
    read(1, '(A)', iostat=ios) buffer
    if (index(buffer,"%rd")) then
    do j = 1,nspecies
        read(1,*,iostat=ios) r_d(j,:)
    end do
    end if

    close(1)
   
    ! save the number of orbital for each atom 
    allocate(atom_orb_num(natoms))
    do i = 1, natoms
        atom_orb_num(i) = orb_num(int(positions(i,4)))
    end do
    allocate(atom_orb_num_pc_l(natoms_pc_l))
    do i = 1, natoms_pc_l
        atom_orb_num_pc_l(i) = orb_num(int(positions_pc_l(i,4)))
    end do
    allocate(atom_orb_num_uc_l(natoms_uc_l))
    do i = 1, natoms_uc_l
        atom_orb_num_uc_l(i) = orb_num(int(positions_uc_l(i,4)))
    end do
    allocate(atom_orb_num_pc_r(natoms_pc_r))
    do i = 1, natoms_pc_r
        atom_orb_num_pc_r(i) = orb_num(int(positions_pc_r(i,4)))
    end do
    allocate(atom_orb_num_uc_r(natoms_uc_r))
    do i = 1, natoms_uc_r
        atom_orb_num_uc_r(i) = orb_num(int(positions_uc_r(i,4)))
    end do
 
    ! apparently this part is not general and can be modified
    ! for different orbitals combinations

    end subroutine load_param



    subroutine get_pos_orb(norb_pc_l,norb_uc_l,&
                           norb_pc_r,norb_uc_r)

    use input
    
    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: norb_pc_l,norb_uc_l,&
                                  norb_pc_r,norb_uc_r
    integer(kind=4) :: i,j,cont
    
    ! get the positions of each orbitals
    cont = 1
    do i = 1,natoms_uc_l
        do j = 1,atom_orb_num_uc_l(i)
            pos_orb_uc_l(cont,:) = positions_uc_l(i,:)
            cont = cont + 1
        end do
    end do
 
    cont = 1
    do i = 1,natoms_pc_l
        do j = 1,atom_orb_num_pc_l(i)
            pos_orb_pc_l(cont,:) = positions_pc_l(i,:)
            cont = cont + 1
        end do
    end do
 
    cont = 1
    do i = 1,natoms_uc_r
        do j = 1,atom_orb_num_uc_r(i)
            pos_orb_uc_r(cont,:) = positions_uc_r(i,:)
            cont = cont + 1
        end do
    end do
 
    cont = 1
    do i = 1,natoms_pc_r
        do j = 1,atom_orb_num_pc_r(i)
            pos_orb_pc_r(cont,:) = positions_pc_r(i,:)
            cont = cont + 1
        end do
    end do
    
    end subroutine get_pos_orb

end module param
