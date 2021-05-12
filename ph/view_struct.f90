module view_struct

implicit none

integer(kind=4) :: natm_sc
real(kind=8),allocatable :: pos_sc(:,:)
real(kind=8),allocatable :: mass_sc(:)
real(kind=8),allocatable :: mass_ev(:)
integer(kind=4),allocatable :: sl_full(:)
integer(kind=4),allocatable :: sl_full_swap(:)
real(kind=8),allocatable :: interfaces_loc(:)
real(kind=8)   :: cell_sc(3,3),reci_cell_sc(3,3)
integer(kind=4),allocatable :: idx_scpc(:) 
real(kind=8)   :: volume
integer(kind=4),allocatable :: layerlist(:,:),layern(:)
integer(kind=4) :: nlay

contains

subroutine gen_struct

use config
use ifport
use util

implicit none

integer(kind=4) :: i,j,k,counts,iat,ii,iii,ilay,idx
real(kind=8) :: distz,randval 
integer(kind=4) :: nmix


CALL SEED(randsd)

natm_sc = nx_sc*ny_sc*nz_sc*natm
nlay = natm*nz_sc
allocate(layerlist(nlay,nx_sc*ny_sc))
allocate(layern(nlay))
layerlist = 0
layern = 0
allocate(pos_sc(natm_sc,3))
allocate(mass_sc(natm_sc))
allocate(mass_ev(3*natm_sc))
allocate(idx_scpc(natm_sc))

counts = 1
allocate(sl_full(4*sum(slz)))
allocate(sl_full_swap(4*sum(slz)))
allocate(interfaces_loc(size(slz,1)+1))
interfaces_loc(1) = -az/8.0d0
do i = 1,size(slz,1)
    do j = 1,slz(i)
        do k = 1,natm
            !if ((k .eq. 2) .or. (k .eq. 4)) then
            !    sl_full(counts) = 3 ! % 3 means As atom
            !else
                ii = mod(i-1,2)+1
                if (ii .eq. 1 ) then
                    sl_full(counts) = prda(k) ! 1 means Si
                    sl_full_swap(counts) = prdb(k) ! 1 means Si
                else
                    sl_full(counts) = prdb(k) ! 2 means Ge
                    sl_full_swap(counts) = prda(k) ! 2 means Ge
                end if
            !end if
            counts = counts + 1
        end do
    end do
    interfaces_loc(i+1) = az * dble(counts-1)/4.0d0-az/8.0d0 
end do

counts = 1
if (mass_read.ne.0) then
    open(23,file='mass_profile.dat',status='old',action='read')
    read(23,*) idx
end if


do i = 1,nx_sc
    do j = 1,ny_sc
        do k = 1,nz_sc
            do iat = 1,natm
                pos_sc(counts+iat-1,:) = pos(iat,:) + matmul((/dble(i-1),dble(j-1),dble(k-1)/),cell)
                idx_scpc(counts+iat-1) = iat
                ii = nint(pos_sc(counts+iat-1,3)/(az/4.0d0))+1
                if ( mass_read.eq.0 ) then
                    mass_sc(counts+iat-1) =  mass_type(sl_full(ii)) ! Si
                    do iii = 1,size(interfaces_loc,1)
                       if (nmixlist(iii) .gt. 0) then
                           nmix = nmixlist(iii)
                           if (abs(pos_sc(counts+iat-1,3) - interfaces_loc(iii)).lt.az/2.0d0*dble(nmix)) then
                                distz = abs(pos_sc(counts+iat-1,3) - interfaces_loc(iii))/(az/4.0d0*dble(nmix/2-0.5))
                                if (random(0) .lt. 0.8*exp(-distz**2)) then
                                    mass_sc(counts+iat-1) = &
                                    mass_type(sl_full_swap(ii)) ! Ge
                                end if
                                if (random(0) .lt. 0.8*exp(-distz**2)) then
                                    if (defecta(iat).gt.0) then
                                        mass_sc(counts+iat-1) = &
                                        mass_d(defecta(iat)) ! Ge
                                    end if
                                end if
                            end if
                        end if
                    end do 
               else
                   read(23,*) idx
                   mass_sc(counts+iat-1) = mass_type(idx)
               end if
               ilay = nint(pos_sc(counts+iat-1,3)/(az/dble(natm)))+1
               layern(ilay) = layern(ilay) + 1
               layerlist(ilay,layern(ilay)) = counts+iat-1
            end do
            counts = counts + natm
        end do
    end do
end do
close(23)
open(unit=1,file="pos_sc.dat",status="UNKNOWN",action="write")
do i = 1,size(pos_sc,1)
    WRITE(1,1000) pos_sc(i,:),mass_sc(i)
end do

mass_ev(1:3*natm_sc-2:3) = mass_sc(:)
mass_ev(2:3*natm_sc-1:3) = mass_sc(:)
mass_ev(3:3*natm_sc  :3) = mass_sc(:)

1000 format(4f20.10)
close(1)
cell_sc(1,:) = dble(nx_sc)*cell(1,:)
cell_sc(2,:) = dble(ny_sc)*cell(2,:)
cell_sc(3,:) = dble(nz_sc)*cell(3,:)
reci_cell_sc = 2*pi*transpose(inv_real(cell_sc))

! volume of the unit cell 
volume = dot_product(cell_sc(1,:),cross(cell_sc(2,:),cell_sc(3,:)))
end subroutine
end module
