module view_struct

implicit none

integer(kind=4) :: natm_sc
real(kind=8),allocatable :: pos_sc(:,:)
real(kind=8),allocatable :: mass_sc(:)
integer(kind=4),allocatable :: sl_full(:)
real(kind=8),allocatable :: interfaces_loc(:)
real(kind=8)   :: cell_sc(3,3),reci_cell_sc(3,3)
integer(kind=4),allocatable :: idx_scpc(:) 
real(kind=8)   :: volume

contains

subroutine gen_struct

use config
use ifport
use util

implicit none

integer(kind=4) :: i,j,k,counts,iat,ii,iii
real(kind=8) :: distz,randval 


CALL SEED (100)


      

natm_sc = nx_sc*ny_sc*nz_sc*natm
allocate(pos_sc(natm_sc,3))
allocate(mass_sc(natm_sc))
allocate(idx_scpc(natm_sc))

counts = 1
allocate(sl_full(4*sum(slz)))
allocate(interfaces_loc(size(slz,1)+1))
interfaces_loc(1) = -az/4.0d0
do i = 1,size(slz,1)
    do j = 1,slz(i)
        do k = 1,natm
            if ((k .eq. 2) .or. (k .eq. 4)) then
                sl_full(counts) = 3 ! % 3 means As atom
            else
                ii = mod(i-1,2)+1
                if (ii .eq. 1 ) then
                    sl_full(counts) = 1 ! 1 means Ga
                else
                    sl_full(counts) = 2 ! 2 means Al
                end if
            end if
            counts = counts + 1
        end do
    end do
    interfaces_loc(i+1) = az * dble(counts-2)/4.0d0 
end do

counts = 1

do i = 1,nx_sc
    do j = 1,ny_sc
        do k = 1,nz_sc
            do iat = 1,natm
                pos_sc(counts+iat-1,:) = pos(iat,:) + matmul((/dble(i-1),dble(j-1),dble(k-1)/),cell)
                idx_scpc(counts+iat-1) = iat
                ii = nint(pos_sc(counts+iat-1,3)/(az/4.0d0))+1
                if (sl_full(ii) .eq. 1) then
                    mass_sc(counts+iat-1) =  masss(1) ! Ga
                    if (nmix .gt. 0) then
                        do iii = 1,size(interfaces_loc,1)
                            if (abs(pos_sc(counts+iat-1,3) - interfaces_loc(iii)).lt.az/2.0d0*dble(nmix)) then
                                distz = abs(pos_sc(counts+iat-1,3) - interfaces_loc(iii))/(az/4.0d0*dble(nmix/2-0.5))
                                if (random(0) .lt. 0.8*exp(-distz**2)) then
                                    mass_sc(counts+iat-1) =  26.981539 ! Al
                                end if
                            end if
                        end do
                    end if
               else if (sl_full(ii) .eq. 2) then
                   mass_sc(counts+iat-1) =  26.981539 ! Al
                   if (nmix .gt. 0) then
                       do iii = 1,size(interfaces_loc,1)
                           if (abs(pos_sc(counts+iat-1,3) - interfaces_loc(iii))<az/2.0d0*dble(nmix)) then
                               distz = abs(pos_sc(counts+iat-1,3) - interfaces_loc(iii))/(az/4.0d0*dble(nmix/2-0.5))
                               if (random(0) .lt. 0.8*exp(-distz**2)) then
                                   mass_sc(counts+iat-1) =  masss(1) ! Ga
                               end if
                           end if
                       end do
                   end if
               else
                   mass_sc(counts+iat-1) =  masss(2) ! As
               end if             
            end do
            counts = counts + natm
        end do
    end do
end do
open(unit=1,file="pos_sc.dat",status="UNKNOWN",action="write")
do i = 1,size(pos_sc,1)
    WRITE(1,1000) pos_sc(i,:),mass_sc(i)
end do
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
