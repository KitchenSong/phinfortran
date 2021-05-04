module write_data
contains 
subroutine write_specular(E,mm,r,k_f,k_b,kp_f,kp_b,v_f,v_b,direction)

use config
use surface
use input

implicit none


real(kind=8),intent(in) :: E
integer(kind=4),intent(in) :: mm
complex(kind=8),intent(in) :: r(:,:)
real(kind=8),intent(in) :: k_f(:,:), k_b(:,:)
real(kind=8),intent(in) :: kp_f(:,:), kp_b(:,:)
real(kind=8),intent(in) :: v_f(:,:), v_b(:,:)
character(len=*) :: direction
real(kind=8) :: k_diff(3)
real(kind=8),allocatable :: r_s(:) ! specular part
real(kind=8),allocatable :: r_t(:) ! specular part+ non-specular part
integer(kind=4) :: i,j, n_f, n_b, ios
character(len=8) :: fmt,x1
integer(kind=4) :: count1

fmt = '(I6.6)'

n_f = size(k_f,1) ! number of forward-propagating modes
n_b = size(k_b,1) ! number of reflected modes

allocate(r_s(n_b),r_t(n_b))

r_s = 0.0d0
r_t = 0.0d0

write (x1,fmt) mm

open(unit=1, file=trim(adjustl(output_dir))//"/"//trim(direction)//"_specular_param_"//trim(x1)//".dat",status="UNKNOWN",action="write",iostat =ios,POSITION='APPEND')

count1 = 0
do i = 1,n_b
    if (k_b(i,3) .lt. 99.9) then
        count1 = count1 + 1
        do j = 1,n_f
            if (k_f(j,3) .lt. 99.9) then
                r_t(i) = r_t(i) + abs(r(j,i))**2
                k_diff = k_f(j,:) - k_b(i,:)
                k_diff(3) =  0.0d0 ! k_f(i,3) - k_b(j,3)
                if (dot_product(k_diff,k_diff)&
                    .lt.1.0d-5) then
                    r_s(i)  = r_s(i) + abs(r(j,i))**2
                end if
            end if
        end do
        if ((r_t(i).lt.1.0d-60).and.(r_t(i).gt.0)) then
            r_t(i) = 1.0d-60
        end if
        if ((r_s(i).lt.1.0d-60).and.(r_s(i).gt.0)) then
            r_s(i) = 1.0d-60
        end if 
        if (r_t(i).ne.0.0) then
            WRITE(1,1000,iostat=ios) count1,E,k_b(i,1),k_b(i,2),k_b(i,3),&
             kp_b(i,1),kp_b(i,2),kp_b(i,3),&
             v_b(i,1),v_b(i,2),v_b(i,3),r_s(i),r_t(i),r_s(i)/r_t(i)
        else
            WRITE(1,1000,iostat=ios) count1,E,k_b(i,1),k_b(i,2),k_b(i,3),&
             kp_b(i,1),kp_b(i,2),kp_b(i,3),&
             v_b(i,1),v_b(i,2),v_b(i,3),r_s(i),r_t(i),0.0d0
        end if

        1000 format(1I5,13e20.10)
    end if
end do

close(unit=1)
end subroutine

subroutine write_md(E,mm,t,r,k_ft,k_f,k_b,kp_ft,kp_f,kp_b,v_ft,v_f,v_b,direction)

use config
use surface
use input

implicit none


real(kind=8),intent(in) :: E
integer(kind=4),intent(in) :: mm
complex(kind=8),intent(in) :: t(:,:),r(:,:)
real(kind=8),intent(in) :: k_f(:,:), k_ft(:,:)
real(kind=8),intent(in) :: kp_f(:,:), kp_ft(:,:)
real(kind=8),intent(in) :: v_f(:,:), v_ft(:,:)
real(kind=8),intent(in) :: k_b(:,:)
real(kind=8),intent(in) :: kp_b(:,:)
real(kind=8),intent(in) :: v_b(:,:)
character(len=*) :: direction
real(kind=8) :: k_diff(3)
real(kind=8),allocatable :: r_s(:) ! specular part
real(kind=8),allocatable :: r_t(:) ! specular part+ non-specular part
real(kind=8),allocatable :: t_s(:) ! specular part
real(kind=8),allocatable :: t_t(:) ! specular part+ non-specular part
integer(kind=4) :: i,j, n_f,n_ft, n_b, ios
character(len=8) :: fmt,x1
integer(kind=4) :: count1

fmt = '(I6.6)'

n_f = size(k_f,1) ! number of forward-propagating modes
n_ft = size(k_ft,1)
n_b = size(k_b,1) ! number of reflected modes

allocate(r_s(n_b),r_t(n_b),t_s(n_b),t_t(n_b))

r_s = 0.0d0
r_t = 0.0d0
t_s = 0.0d0
t_t = 0.0d0

write (x1,fmt) mm

open(unit=1, file=trim(adjustl(output_dir))//"/"//trim(direction)//"_specular_param_"//trim(x1)//".dat",status="UNKNOWN",action="write",iostat =ios,POSITION='APPEND')

count1 = 0
do i = 1,n_b
    if (k_b(i,3) .lt. 99.9) then
        count1 = count1 + 1
        do j = 1,n_f
            if (k_f(j,3) .lt. 99.9) then
                r_t(i) = r_t(i) + abs(r(j,i))**2
                k_diff = k_f(j,:) - k_b(i,:)
                k_diff(3) =  0.0d0 ! k_f(i,3) - k_b(j,3)
                if (dot_product(k_diff,k_diff)&
                    .lt.1.0d-5) then
                    r_s(i)  = r_s(i) + abs(r(j,i))**2
                end if
            end if
        end do
        if ((r_t(i).lt.1.0d-60).and.(r_t(i).gt.0)) then
            r_t(i) = 1.0d-60
        end if
        if ((r_s(i).lt.1.0d-60).and.(r_s(i).gt.0)) then
            r_s(i) = 1.0d-60
        end if 

        do j = 1,n_ft
            if (k_ft(j,3) .lt. 99.9) then
                t_t(i) = t_t(i) + abs(t(j,i))**2
                k_diff = k_ft(j,:) - k_b(i,:)
                k_diff(3) =  0.0d0 ! k_f(i,3) - k_b(j,3)
                if (dot_product(k_diff,k_diff)&
                    .lt.1.0d-5) then
                    t_s(i)  = t_s(i) + abs(t(j,i))**2
                end if
            end if
        end do
        if ((t_t(i).lt.1.0d-60).and.(t_t(i).gt.0)) then
            t_t(i) = 1.0d-60
        end if
        if ((t_s(i).lt.1.0d-60).and.(t_s(i).gt.0)) then
            t_s(i) = 1.0d-60
        end if 
        
        WRITE(1,1000,iostat=ios) count1,E,k_b(i,1),k_b(i,2),k_b(i,3),&
             kp_b(i,1),kp_b(i,2),kp_b(i,3),&
             v_b(i,1),v_b(i,2),v_b(i,3),r_s(i),r_t(i),t_s(i),t_t(i)


        1000 format(1I5,15e20.10)
    end if
end do

close(unit=1)
end subroutine



subroutine write_md2md(E,mm,r,k_f,k_b,kpc_f,kpc_b,v_f,v_b,direction)
! mode to mode transmission
use config
use surface
use input

implicit none


real(kind=8),intent(in) :: E
integer(kind=4),intent(in) :: mm
complex(kind=8),intent(in) :: r(:,:)
real(kind=8),intent(in) :: k_f(:,:), k_b(:,:)
real(kind=8),intent(in) :: kpc_f(:,:), kpc_b(:,:)
real(kind=8),intent(in) :: v_f(:,:), v_b(:,:)
character(len=*) :: direction
real(kind=8) :: k_diff(3)
real(kind=8) :: ts,tns ! specular & non-specular part
integer(kind=4) :: i,j, n_f, n_b, ios
integer(kind=4) :: count1,count2
character(len=8) :: fmt,x1

fmt = '(I6.6)'

n_f = size(k_f,1) ! number of forward-propagating modes
n_b = size(k_b,1) ! number of reflected modes

write (x1,fmt) mm
open(unit=1, file=trim(adjustl(output_dir))//"/"//trim(direction)//"_md2md_"//trim(x1)//".dat",status="UNKNOWN",action="write",iostat =ios,POSITION='APPEND')

count1 = 0
count2 = 0

do j = 1,n_b
    if (k_b(j,3) .lt. 99.9) then
        count2 = count2 + 1
        do i = 1,n_f
            if (k_f(i,3) .lt. 99.9) then
                count1 = count1 + 1
                k_diff = k_f(i,:) - k_b(j,:)
                k_diff(3) =  0.0d0 ! k_f(i,3) - k_b(j,3)
                ts = 0.0d0
                tns = 0.0d0
                if (dot_product(k_diff,k_diff)&
                    .lt.1.0d-5) then
                    ts = abs(r(i,j))**2
                else
                    tns = abs(r(i,j))**2
                end if
      
                WRITE(1,1000,iostat=ios) count1,&
                count2,E,&
                k_f(i,1),k_f(i,2),k_f(i,3),&
                k_b(j,1),k_b(j,2),k_b(j,3),&
                kpc_f(i,1),kpc_f(i,2),kpc_f(i,3),&
                kpc_b(j,1),kpc_b(j,2),kpc_b(j,3),&
                v_f(i,1),v_f(i,2),v_f(i,3),&
                v_b(j,1),v_b(j,2),v_b(j,3),&
                ts,tns,i,j
            end if   
        end do

        1000 format(2I5,21e20.10,2I5)
    end if
end do

close(unit=1)
end subroutine

!subroutine write_t_vs_E(E,mm,ts,t)
!end subroutine

subroutine write_ldos(E,mm,ldos)

use config

implicit none

real(kind=8),intent(in) :: E
integer(kind=4),intent(in) :: mm
real(kind=8),intent(in) :: ldos(:,:)
character(len=8) :: fmt,x1
integer(kind=4) :: i,j,nl,no,ios

nl = size(ldos,1)
no = size(ldos,2)

fmt = '(I6.6)'
write (x1,fmt) mm

open(unit=1, file=trim(adjustl(output_dir))//"/ldos_"//trim(x1)//".dat",status="UNKNOWN",action="write",iostat =ios,POSITION='APPEND')

do i = 1,nl
    WRITE(1,1000,iostat=ios,advance='no') E
    do j = 1,no
        WRITE(1,1000,iostat=ios,advance='no') ldos(i,j)
        1000 format(1f20.10)
    end do
    WRITE(1,"(A)",advance="yes") " "
end do

close(unit=1)
end subroutine write_ldos

subroutine write_surface_ldos(sdosl,sdosr)

use config

implicit none

real(kind=8),intent(in) :: sdosl(:,:),sdosr(:,:)
character(len=8) :: fmt,x1
integer(kind=4) :: i,j,nl,no,ios

nl = size(sdosl,1)
no = size(sdosl,2)

fmt = '(I6.6)'

open(unit=1, file=trim(adjustl(output_dir))//"/sdosl.dat",status="UNKNOWN",action="write",iostat =ios,POSITION='APPEND')
open(unit=2, file=trim(adjustl(output_dir))//"/sdosr.dat",status="UNKNOWN",action="write",iostat =ios,POSITION='APPEND')

do i = 1,nl
    do j = 1,no
        WRITE(1,1000,iostat=ios,advance='no') sdosl(i,j)
        WRITE(2,1000,iostat=ios,advance='no') sdosr(i,j)
        1000 format(1f20.10)
    end do
    WRITE(1,"(A)",advance="yes") " "
    WRITE(2,"(A)",advance="yes") " "
end do

close(unit=1)
close(unit=2)
end subroutine write_surface_ldos



subroutine write_positions(positions_all)

use config

implicit none

real(kind=8),intent(in) :: positions_all(:,:)
integer(kind=4) :: i,j,m,n

open(unit=1, file=trim(adjustl(output_dir))//"/nkpoints.dat",status="REPLACE",action="write")
WRITE(1,1001) dble(nk)
1001 format(1f20.10)

close(unit=1)


open(unit=1, file=trim(adjustl(output_dir))//"/positions_all.dat",status="REPLACE",action="write")
m = size(positions_all,1)
do i = 1,m
    WRITE(1,1000) positions_all(i,:)
    1000 format(5f20.10)
end do
close(unit=1)
end subroutine

subroutine write_trans(E,transl_s,transl_ns,transr_s,transr_ns,transl,transr,&
                      refll_s,refll_ns,reflr_s,reflr_ns)

use config

implicit none

real(kind=8) :: transl_s(:,:),transl_ns(:,:),&
                transr_s(:,:),transr_ns(:,:),&
                transl(:,:), transr(:,:),&
                refll_s(:,:),refll_ns(:,:),&
                reflr_s(:,:),reflr_ns(:,:)

real(kind=8),intent(in) :: E(:,:)
integer(kind=4) :: i,j,m,n,m1,n1

m = size(transl_s,1)
n = size(transl_s,2)
m1 = size(transr_s,1)
n1 = size(transr_s,2)

open(unit=1, file=trim(adjustl(output_dir))//"/transl_s.dat",status="UNKNOWN",action="write")
open(unit=2, file=trim(adjustl(output_dir))//"/transl_ns.dat",status="UNKNOWN",action="write")

do i = 1,m
    WRITE(1,1000,advance='no') E(i,1)
    WRITE(2,1000,advance='no') E(i,1)
    do j = 1,n
        if ( transl_s(i,j) .lt. 1.0d-60) then
             transl_s(i,j) = 0.0d0
        end if
        if ( transl_ns(i,j) .lt. 1.0d-60) then
             transl_ns(i,j) = 0.0d0
        end if


        WRITE(1,1000,advance='no') transl_s(i,j)
        WRITE(2,1000,advance='no') transl_ns(i,j)
    end do
    WRITE(1,"(A)",advance="yes") " "
    WRITE(2,"(A)",advance="yes") " "
end do

close(unit=1)
close(unit=2)

open(unit=1, file=trim(adjustl(output_dir))//"/refll_s.dat",status="UNKNOWN",action="write")
open(unit=2, file=trim(adjustl(output_dir))//"/refll_ns.dat",status="UNKNOWN",action="write")

do i = 1,m
    WRITE(1,1000,advance='no') E(i,1)
    WRITE(2,1000,advance='no') E(i,1)
    do j = 1,n
        if ( refll_s(i,j) .lt. 1.0d-60) then
             refll_s(i,j) = 0.0d0
        end if
        if ( refll_ns(i,j) .lt. 1.0d-60) then
             refll_ns(i,j) = 0.0d0
        end if


        WRITE(1,1000,advance='no') refll_s(i,j)
        WRITE(2,1000,advance='no') refll_ns(i,j)
!        if  (refll_ns(i,j) .gt. 1.5) then
!            write(*,*) i,j
           ! stop
!        end if
    end do
    WRITE(1,"(A)",advance="yes") " "
    WRITE(2,"(A)",advance="yes") " "
end do

close(unit=1)
close(unit=2)

open(unit=1, file=trim(adjustl(output_dir))//"/transr_s.dat",status="UNKNOWN",action="write")
open(unit=2, file=trim(adjustl(output_dir))//"/transr_ns.dat",status="UNKNOWN",action="write")

do i = 1,m1
    WRITE(1,1000,advance='no') E(i,1)
    WRITE(2,1000,advance='no') E(i,1)
    do j = 1,n1
        if ( transr_s(i,j) .lt. 1.0d-60) then
            transr_s(i,j) = 0.0d0
        end if
        if ( transr_ns(i,j) .lt. 1.0d-60) then
            transr_ns(i,j) = 0.0d0
        end if

        WRITE(1,1000,advance='no') transr_s(i,j)
        WRITE(2,1000,advance='no') transr_ns(i,j)
    !    if  (transr_ns(i,j) .gt. 1.5) then
    !        write(*,*) i,j
            !stop
    !    end if
    end do
    WRITE(1,"(A)",advance="yes") ""
    WRITE(2,"(A)",advance="yes") ""
end do

close(unit=1)
close(unit=2)

open(unit=1, file=trim(adjustl(output_dir))//"/reflr_s.dat",status="UNKNOWN",action="write")
open(unit=2, file=trim(adjustl(output_dir))//"/reflr_ns.dat",status="UNKNOWN",action="write")

do i = 1,m1
    WRITE(1,1000,advance='no') E(i,1)
    WRITE(2,1000,advance='no') E(i,1)
    do j = 1,n1
!        if  (reflr_ns(i,j) .gt. 1.5) then
!            write(*,*) i,j
            !stop
!        end if
        if ( reflr_s(i,j) .lt. 1.0d-60) then
            reflr_s(i,j) = 0.0d0
        end if
        if ( reflr_ns(i,j) .lt. 1.0d-60) then
            reflr_ns(i,j) = 0.0d0
        end if

        WRITE(1,1000,advance='no') reflr_s(i,j)
        WRITE(2,1000,advance='no') reflr_ns(i,j)
    end do
    WRITE(1,"(A)",advance="yes") ""
    WRITE(2,"(A)",advance="yes") ""
end do

close(unit=1)
close(unit=2)

open(unit=1, file=trim(adjustl(output_dir))//"/transl_tot.dat",status="UNKNOWN",action="write")
open(unit=2, file=trim(adjustl(output_dir))//"/transr_tot.dat",status="UNKNOWN",action="write")

do i = 1,m1
    WRITE(1,1000,advance='no') E(i,1)
    WRITE(2,1000,advance='no') E(i,1)
    do j = 1,n1
        if (transl(i,j) .lt. 1.0d-60) then
            transl(i,j) = 0.0d0
        end if
        if (transr(i,j) .lt. 1.0d-60) then
            transr(i,j) = 0.0d0
        end if


        WRITE(1,1000,advance='no') transl(i,j)
        WRITE(2,1000,advance='no') transr(i,j)
    end do
    WRITE(1,"(A)",advance="yes") ""
    WRITE(2,"(A)",advance="yes") ""
end do
1000 format(1e24.16)

close(unit=1)
close(unit=2)

end subroutine

subroutine write_eig(fname,eigs)

use config

implicit none

character(len=*) :: fname
real(kind=8) :: eigs(:)
integer(kind=4) :: n,i

n = size(eigs,1)

open(unit=1, file=trim(adjustl(output_dir))//"/"//trim(adjustl(fname))//".dat",status="unknown",action="write",access="append")
do i = 1,n
    WRITE(1,1000,advance='no') eigs(i)
end do
WRITE(1,"(A)",advance="yes") ""
1000 format(1e20.10)

close(unit=1)


end subroutine write_eig

subroutine write_vel(fname,vels)

use config

implicit none

character(len=*) :: fname
real(kind=8) :: vels(:,:)
integer(kind=4) :: n,n1,i,j

n = size(vels,1)
n1 = size(vels,2)


open(unit=1, file=trim(adjustl(output_dir))//"/"//trim(adjustl(fname))//".dat",status="unknown",action="write",access="append")
do i = 1,n
    do j = 1,n1
        WRITE(1,1000,advance='no') vels(i,j)
    end do
end do
WRITE(1,"(A)",advance="yes") ""
1000 format(1e20.10)

close(unit=1)


end subroutine write_vel



end module write_data
