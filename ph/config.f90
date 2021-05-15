module config

implicit none

integer(kind=4) :: natm,ntype,ibrav
real(kind=8)    :: celldm(6)
real(kind=8),allocatable :: masss(:),pos(:,:)
real(kind=8) :: cell(3,3),cell_sc_dfpt(3,3)
real(kind=8) :: reci_cell(3,3),ivcell(3,3)
integer(kind=4) :: nx,ny,nz,nk
integer(kind=4) :: nx_sc,ny_sc,nz_sc
real(kind=8),allocatable ::fc(:,:,:,:,:,:,:)
real(kind=8) :: epsil(3,3)
real(kind=8),allocatable :: born(:,:,:)
character(len=30) :: filename_input
real(kind=8)      :: sigma ! width for delta function in cm-1
real(kind=8)      :: az ! the period length in z direction
integer(kind=4),allocatable   :: slz(:)
integer(kind=4),allocatable   :: nmixlist(:)
integer(kind=4)   :: sl(10000)
integer(kind=4)   :: mix(10000) 
integer(kind=4)   :: nsl
integer(kind=4)   :: nxy(2)
integer(kind=4)   :: ne
real(kind=8) :: emin,emax
real(kind=8),allocatable :: egrid(:)
integer(kind=4)   :: ipolar
integer(kind=4)   :: randsd,verbose
integer(kind=4)   :: dos,ndosx,ndosy,ndosz
real(kind=8),allocatable  :: kmesh 
real(kind=8)     :: mass_species(100)
integer(kind=4)  :: nspecies
real(kind=8),allocatable  :: mass_type(:)
real(kind=8)     :: mass_defect(100)
real(kind=8),allocatable  :: mass_d(:)
integer(kind=4)  :: perioda(100),periodb(100),defecta(100),defectb(100)
integer(kind=4),allocatable  :: prda(:),prdb(:)
integer(kind=4)  :: mass_read,unfold
real(kind=8)     :: tau



namelist/configlist/sigma,filename_input,az,sl,nxy,mix,emin,emax,ne,ipolar,nk,randsd,verbose,dos,ndosx,ndosy,ndosz,mass_species,nspecies,perioda,periodb,defecta,defectb,mass_defect,mass_read,unfold,tau
contains
    
subroutine load_configure()

use util
    
implicit none

integer(kind=4) :: i,j,itemp,jtemp,alpha,beta,iatm,jatm
integer(kind=4) :: idx,jdx,kdx
character(len=30) :: label
real(kind=8)   :: vol,fctemp
integer(kind=4) :: counts

unfold = 0
ne = 10
emax = 10
emin = 0
ipolar = 0 ! default: no long-range force constant
nk = 10
verbose = 0
dos = 0
ndosx = 1
ndosy = 1
ndosz = 1
nspecies = 1
mass_species = 0.0d0
perioda=0
periodb=0
defecta=0
defectb=0
mass_defect=0
mass_read=0
tau = 0.1 ! Linewidth in unit of THz = 2*hbar/lifetime

open(1,file="config",status="old")
read(1,nml=configlist)
counts = 1
do while (mass_species(counts).gt.0)
    counts = counts + 1
end do
allocate(mass_type(counts-1))
if (counts-1 .ne. nspecies) then
    write(*,*) 'nspecies incorrect!'
    stop
end if
mass_type(1:counts-1) = mass_species(1:counts-1)
counts = 1
do while (perioda(counts).gt.0)
    counts = counts + 1
end do
allocate(prda(counts-1))
prda(1:counts-1)=perioda(1:counts-1)
counts = 1
do while (periodb(counts).gt.0)
    counts = counts + 1
end do
allocate(prdb(counts-1))
prdb(1:counts-1)=periodb(1:counts-1)
counts = 1
do while (mass_defect(counts).gt.0)
    counts = counts + 1
end do
allocate(mass_d(counts-1))
mass_d(1:counts-1)=mass_defect(1:counts-1)

! energy grid for spectral function
allocate(egrid(ne))
do i = 1,ne
    egrid(i) = dble(i-1)/dble(ne-1)*(emax-emin) + emin
end do

counts = 1
do while (sl(counts).gt.0)
    counts = counts + 1
end do
allocate(slz(counts-1))
allocate(nmixlist(counts-1))
slz(1:counts-1) = sl(1:counts-1)
nmixlist(1:counts-1) = mix(1:counts-1)
nx_sc = nxy(1)
ny_sc = nxy(2)
nz_sc = sum(slz)

close(1)
open(23,file=trim(adjustl(filename_input)),status='old',action='read')
read(23,*) ntype,natm,ibrav,celldm
do i = 1,3
    read(23,*) cell(i,:)
end do
! number of orbitals

allocate(masss(ntype))
allocate(pos(natm,3))
masss(:) = 0.0d0

do i = 1,ntype
    read(23,*) itemp, label, vol 
    masss(i) = ele_mass(label)
end do
do i = 1,natm
    read(23,*) itemp,jtemp, pos(i,:)
end do
pos(:,:) = pos(:,:) * celldm(1) * bohr2ang
cell(:,:) = cell(:,:) * celldm(1) * bohr2ang
reci_cell= 2*pi*transpose(inv_real(cell))
read(23,*) label
allocate(born(natm,3,3))
if (trim(adjustl(label)) .eq. 'T') then
    do i = 1,3
        read(23,*) epsil(i,:)
    end do
    do i = 1,natm
        read(23,*) itemp
        do j = 1,3
            read(23,*) born(itemp,j,:)
        end do
    end do
end if
read(23,*) nx,ny,nz
cell_sc_dfpt(1,:) = nx*cell(1,:)
cell_sc_dfpt(2,:) = ny*cell(2,:)
cell_sc_dfpt(3,:) = nz*cell(3,:)
allocate(fc(3,3,natm,natm,nx,ny,nz))
do i = 1, (natm*3)**2
    read(23,*) alpha,beta,iatm,jatm
    do j = 1,nx*ny*nz
        read(23,*) idx,jdx,kdx,fctemp
        fc(alpha,beta,iatm,jatm,idx,jdx,kdx) = fctemp
    end do
end do
close(23)
! acoustic sum rule
do i = 1,3
    do j = 1,3
        do idx = 1,natm
            fc(i,j,idx,idx,1,1,1) = -sum(fc(i,j,idx,:,:,:,:))&
            +fc(i,j,idx,idx,1,1,1)
        end do
    end do
end do

ivcell = inv_real(cell)



end subroutine

function ele_mass(ele) result(m)

implicit none

character(len=30),intent(in) :: ele
real(kind=8)    :: m

if (trim(adjustl(ele)) .eq. 'Si') then
    m = 28.0855
else if (trim(adjustl(ele)) .eq. 'Ga') then
    m = 69.723
else if (trim(adjustl(ele)) .eq. 'As') then
    m = 74.921595
end if

end function ele_mass



end module
