module dyn

real(kind=8) :: cell_g(3,3),cell_g_len(3)
integer(kind=4) :: ncell_g(3)

contains

subroutine gen_dyn()

use util
use config
use view_struct

implicit none

integer(kind=4) :: nscx,nscy,nscz
integer(kind=4) :: a,b,c
integer(kind=4) :: i,ii,j,jj,k,kk,ir,counts,m1,m2,m3,ie
real(kind=8) :: kpoint1(3),kpoint2(3),kpoint3(3)
real(kind=8),allocatable :: kps(:,:)
real(kind=8),allocatable :: kabs(:)
complex(kind=8),allocatable :: dynmat(:,:,:),dynmat_g(:,:,:)
complex(kind=8),allocatable :: matmat(:,:) 
real(kind=8),allocatable    :: omega(:,:)
complex(kind=8),allocatable :: evs(:,:)
real(kind=8)    :: rws(124,3),rd(124)
real(kind=8)    :: gmax,alpha,geg,exp_g
real(kind=8)    :: g(1,3),g_old(3),zig(1,3),zjg(1,3),auxi(3),gr,gtemp(1,3)
real(kind=8)    :: pos_i(3),pos_j(3),pos_i_pc(3),pos_j_pc(3),total_weight,weight
real(kind=8)    :: dr(3),ixyz(3),df
integer(kind=4) :: ipol,jpol,iat,jat,iidim,jdim,ik,nreq
integer(kind=4) :: aa,bb,cc,t1,t2,t3,m,n,nn,ll,ikz
real(kind=8)    :: kpoint(3),ktemp(3),mi,mj
real(kind=8),allocatable  :: Ggrid(:,:,:,:),weightk(:,:,:)
real(kind=8),allocatable :: spectral(:,:,:,:,:)
real(kind=8),allocatable :: spectral_proj(:,:,:,:,:,:,:)
real(kind=8),allocatable :: kmeshx(:,:,:,:,:),kmeshy(:,:,:,:,:),kmeshz(:,:,:,:,:)
real(kind=8),allocatable :: emesh(:,:,:,:,:)
real(kind=8),allocatable :: datamesh(:,:)
real(kind=8),allocatable :: datamesh_proj(:,:,:,:)
integer(kind=4) :: ib

kpoint1 = (/0.0,0.0,0.0/)
kpoint2 = (/0.0,0.0,-0.5/)
kpoint3 = (/0.0,0.0,0.5/)

nscx = ceiling(dble(nx)/dble(nx_sc))
nscy = ceiling(dble(ny)/dble(ny_sc))
nscz = ceiling(dble(nz)/dble(nz_sc))

if (dos.ne.0) then
    nk = ndosx*ndosy*ndosz
end if

allocate(kps(nk,3))
allocate(kabs(nk))
kps = 0.0d0
kabs = 0.0d0
if (dos.eq.0) then
    do i = 1,3
        kps(1:(nk)/2+1,i) = linspace(kpoint2(i),kpoint1(i),(nk)/2+1,1)    
    end do
    do i = 1,(nk)/2+1
        kabs(i) = sqrt(dot_product(kps(i,:)-kpoint2,kps(i,:)-kpoint2))
    end do
    do i = 1,3
        kps((nk)/2+1:nk,i) = linspace(kpoint1(i),kpoint3(i),(nk)/2,0)    
    end do
    do i = 1,(nk)/2
       kabs(i+(nk)/2) = kabs((nk)/2+1)+sqrt(dot_product(kps(i+(nk)/2,:)-kpoint1,kps(i+(nk)/2,:)-kpoint1))
    end do
else
    counts = 1
    do i= 1,ndosx
        do j = 1,ndosy
            do k=1,ndosz
                kps(counts,:) =(/dble(i-1)/dble(ndosx),dble(j-1)/dble(ndosy),dble(k-1)/dble(ndosz)/)
                counts = counts + 1
            end do
        end do
    end do
end if

allocate(dynmat(nk,natm_sc*3,natm_sc*3))
allocate(dynmat_g(nk,natm_sc*3,natm_sc*3))
allocate(omega(nk,natm_sc*3))
allocate(evs(natm_sc*3,natm_sc*3))


rws = 0.0d0
rd = 0.0d0
j = 1
do m1=-2,2
    do m2=-2,2
        do m3=-2,2
             if ((m1.eq.0) .and. (m2.eq.0) .and. (m3.eq.0)) then
                 cycle
             end if
             do i=1,3
                rws(j,i)=cell_sc_dfpt(1,i)*m1+cell_sc_dfpt(2,i)*m2+cell_sc_dfpt(3,i)*m3
             end do 
             rd(j)=0.5*dot_product(rws(j,1:3),rws(j,1:3))
             j=j+1
         end do
    end do
end do

cell_g = reci_cell_sc/angbohr
gmax = 14.0d0
alpha= (2.0*pi/sqrt(dot_product(cell_sc(1,:),cell_sc(1,:)))/angbohr)**2
geg=gmax*4.*alpha
do i = 1,3
    cell_g_len(i) = sqrt(dot_product(cell_g(i,:),cell_g(i,:)))
end do
do i = 1,3
    ncell_g(i)=nint(sqrt(geg)/cell_g_len(i))+1
end do
do m1=-ncell_g(1),ncell_g(1)
  do m2=-ncell_g(2),ncell_g(2)
     do m3=-ncell_g(3),ncell_g(3)
        g(1,:) =dble(m1)*cell_g(1,1:3)+dble(m2)*cell_g(2,1:3)+dble(m3)*cell_g(3,1:3)
        gtemp = transpose(matmul(epsil,transpose(g)))
        geg=dot_product(g(1,:),gtemp(1,:))
        if ((geg.gt.0.0d0) .and. (geg/alpha/4.0d0.lt.gmax)) then
            exp_g=exp(-geg/alpha/4.0d0)/geg
            do iat=1,natm_sc
               zig = matmul(g,born(idx_scpc(iat),1:3,1:3))
               auxi =  0.0
               do jat=1,natm_sc
                   gr=dot_product(g(1,:),pos_sc(jat,:)-pos_sc(iat,:))
                   zjg = matmul(g,born(idx_scpc(jat),1:3,1:3))
                   auxi(1:3)=auxi(1:3)+zjg(1,1:3)*real(exp(-i_imag*gr*angbohr))
                end do
                do ipol=1,3
                 iidim=(iat-1)*3+ipol
                 do jpol=1,3
                    jdim=(iat-1)*3+jpol
                    dynmat_g(1:nk,iidim,jdim)=dynmat_g(1:nk,iidim,jdim)-exp_g*zig(1,ipol)*auxi(jpol)
                 end do
                end do
             end do
        end if
        g_old(1:3)=g(1,1:3)
        do ik=1,nk
           g(1,1:3)=g_old(1:3)+matmul(kps(ik,1:3),reci_cell_sc)/angbohr
           gtemp = transpose(matmul(epsil,transpose(g)))
           geg=dot_product(g(1,:),gtemp(1,:))
           if ((geg.gt.0.0d0).and.(geg/alpha/4.0d0.lt.gmax)) then
              exp_g=exp(-geg/alpha/4.0d0)/geg
              do iat=1,natm_sc
                 zig = matmul(g,born(idx_scpc(iat),1:3,1:3))
                 do jat=1,natm_sc
                    gr=dot_product(g(1,:),pos_sc(jat,:)-pos_sc(iat,:))
                    zjg =matmul(g,born(idx_scpc(jat),1:3,1:3))
                    do ipol=1,3
                       iidim=(iat-1)*3+ipol
                       do jpol=1,3
                          jdim=(jat-1)*3+jpol
                          dynmat_g(ik,iidim,jdim)=dynmat_g(ik,iidim,jdim)+exp_g*zig(1,ipol)*zjg(1,jpol)*exp(-i_imag*gr*angbohr)
                       end do
                    end do
                  end do
               end do
             end if 
        end do
     end do
   end do
end do
dynmat_g = 8*pi*dynmat_g/volume/angbohr**3

do ik = 1,nk
    kpoint = kps(ik,:)
    do a = -nscx,nscx
        do b = -nscy,nscy
            do c = -nscz,nscz
                do i = 1,natm_sc
                   pos_i = pos_sc(i,:)
                   pos_i_pc = pos(idx_scpc(i),:)
                   do j = 1,natm_sc
                       pos_j = pos_sc(j,:)
                       pos_j_pc = pos(idx_scpc(j),:)
                       total_weight = 0.0
                       g(1,:) = matmul(kpoint,reci_cell_sc)
                       weight = 0.0d0
                       dr = -matmul((/dble(a),dble(b),dble(c)/),cell_sc) &
                       +pos_j-pos_i
                       ixyz = matmul(dr-(pos_j_pc-pos_i_pc),ivcell)
                       nreq = 1
                       jj = 0
                       do ir = 1,124
                           df = dot_product(dr,rws(ir,1:3))-rd(ir)
                           if (df.gt.1.0d-5) then
                               jj = 1
                               cycle
                           end if
                           if (abs(df).lt.1.0d-5) then
                               nreq = nreq+1
                           end if
                       end do
                       aa = -nint(ixyz(1))
                       bb = -nint(ixyz(2))
                       cc = -nint(ixyz(3))
                       if (jj.eq.0) then
                           weight = 1.0/dble(nreq)
                       end if
                       if (weight.gt.0) then
                           t1 = mod(aa+1,nx)
                           if (t1.le.0) then
                               t1 = t1+nx
                           end if
                           t2 = mod(bb+1,ny)
                           if (t2.le.0) then
                               t2 = t2+ny
                           end if
                           t3 = mod(cc+1,nz)
                           if (t3.le.0) then
                               t3 = t3+nz
                           end if
                           do m = 1,3
                               do n = 1,3
                                    dynmat(ik,(i-1)*3+m,(j-1)*3+n) = dynmat(ik,(i-1)*3+m,(j-1)*3+n)&
                                +weight*fc(m,n,idx_scpc(i),idx_scpc(j),t1,t2,t3)&
                                *exp(i_imag*dot_product(dr-(pos_j-pos_i),matmul(kpoint,reci_cell_sc)))
                               end do
                            end do



                       end if
                   end do
                end do
            end do
        end do
    end do
end do
do ik = 1,nk
    do iat = 1,natm_sc
        mi = mass_sc(iat)
        do jat = 1,natm_sc
            mj = mass_sc(jat)
            do ipol = 1,3
                do jpol = 1,3
                    dynmat(ik,(iat-1)*3+ipol,(jat-1)*3+jpol) = &
                    dynmat(ik,(iat-1)*3+ipol,(jat-1)*3+jpol)/sqrt(mi*mj)
                end do
            end do
        end do
    end do
end do

! possible reciprocal lattice vectors G
allocate(Ggrid(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,3))
Ggrid = 0.0d0
do i = 1,2*nx_sc+1
    do j = 1,2*ny_sc+1
        do k = 1,2*nz_sc+1
            Ggrid(i,j,k,:) = matmul((/dble(i-nx_sc-1),dble(j-ny_sc-1),dble(k-nz_sc-1)/),reci_cell_sc)
        end do
    end do
end do
allocate(weightk(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1))
allocate(spectral(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,nk,ne))
allocate(spectral_proj(3,nlay,2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,nk,ne))
allocate(kmeshx(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,nk,ne))
allocate(kmeshy(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,nk,ne))
allocate(kmeshz(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,nk,ne))
allocate(emesh(2*nx_sc+1,2*ny_sc+1,2*nz_sc+1,nk,ne))
allocate(matmat(3*natm_sc,3*natm_sc))
spectral = 0.0d0
spectral_proj = 0.0d0

open(unit=4,file="eigen.dat",status="UNKNOWN",action="write")


do ik = 1,nk
    matmat = (dynmat(ik,:,:)+dble(ipolar)*dynmat_g(ik,:,:)) &
    *eleVolt/1.0d-20/mass_proton/(1.0d12*2.0d0*pi)**2*33.35641**2*(13.605662285137/0.529177249**2)
    call eigenH(matmat,omega(ik,:),evs)
    if ((dos.ne.0).and.(verbose.ne.0)) then
        call write_evs(ik,evs,kps(ik,:))
    else if ((dos.ne.0).and.(verbose.eq.0)) then
        call write_evs_avg(ik,evs,kps(ik,:))
    end if
    omega(ik,:) = real(sqrt(abs(omega(ik,:))))
    
    do ib = 1,natm_sc*3
        weightk(:,:,:) = 0.0
        do ii = 1,2*nx_sc+1
            do jj = 1,2*ny_sc+1
                do kk = 1,2*nz_sc+1
                    do ie=1,size(egrid,1)
                        ktemp = matmul(kps(ik,:),reci_cell_sc) + Ggrid(ii,jj,kk,:)
                        kmeshx(ii,jj,kk,ik,ie) = ktemp(1)
                        kmeshy(ii,jj,kk,ik,ie) = ktemp(2)
                        kmeshz(ii,jj,kk,ik,ie) = ktemp(3)
                        emesh(ii,jj,kk,ik,ie) = egrid(ie) 
                    end do
                    weightk(ii,jj,kk) = weight_k(matmul(kps(ik,:),reci_cell_sc),Ggrid(ii,jj,kk,:),evs(:,ib))
                    spectral(ii,jj,kk,ik,:) = spectral(ii,jj,kk,ik,:)+ weightk(ii,jj,kk)*exp(-0.5*(omega(ik,ib)-egrid(:))**2/sigma**2)&
                        /sigma/sqrt(2*pi)
                    do nn= 1,nlay 
                        do ll = 1,layern(nn)
                            do alpha = 1,3
                                spectral_proj(alpha,nn,ii,jj,kk,ik,:) =spectral_proj(alpha,nn,ii,jj,kk,ik,:)&
                +abs(evs((layerlist(nn,ll)-1)*3+alpha,ib)/sqrt(mass_sc(layerlist(nn,ll))))**2&
                 *weightk(ii,jj,kk)*exp(-0.5*(omega(ik,ib)-egrid(:))**2/sigma**2)&
                /sigma/sqrt(2*pi)
                            end do
                        end do
                    end do 
                end do
            end do
        end do
        write(4,2000,advance="no") omega(ik,ib)
    end do
    WRITE(4,"(A)",advance="yes") " "
    write(*,140) dble(ik)/dble(nk)*100.0d0,'%'
end do
2000 format(1f10.5)
140 format(F8.2,A)
close(4)


allocate(datamesh((nk)*nz_sc+1,ne))
allocate(datamesh_proj(3,nlay,(nk)*nz_sc+1,ne))
datamesh = 0.0d0
datamesh_proj = 0.0d0
!! formatted
do ik = 1,nk
!    do ib = 1,natm_sc*3
        do ii = 1,2*nx_sc+1
            do jj = 1,2*ny_sc+1
                do kk = 1,2*nz_sc+1
                    if ((kmeshz(ii,jj,kk,ik,1).ge. (-1.0*pi/az)) .and.(kmeshz(ii,jj,kk,ik,1).le. (1.0*pi/az)).and.(abs(kmeshx(ii,jj,kk,ik,1)).lt.1.0d-4).and.(abs(kmeshy(ii,jj,kk,ik,1)).lt.1.0d-4)) then 
!                        write(*,*)(kmeshz(ii,jj,kk,ik,1)+1.0*pi/az)/(2*pi/az/nz_sc/(nk-1))-nint((kmeshz(ii,jj,kk,ik,1)+1.0*pi/az)/(2*pi/az/nz_sc/(nk-1)))
!                    do ie=1,size(egrid,1)
                         ikz = nint((kmeshz(ii,jj,kk,ik,1)+1.0*pi/az)/(2*pi/az/nz_sc/(nk)))+1
                         datamesh(ikz,:) = datamesh(ikz,:) +&
                         spectral(ii,jj,kk,ik,:)
                         do nn= 1,nlay 
                             do alpha=1,3
                                 datamesh_proj(alpha,nn,ikz,:) = datamesh_proj(alpha,nn,ikz,:) +&
                                 spectral_proj(alpha,nn,ii,jj,kk,ik,:)
                             end do                         
                         end do
                        
!                        WRITE(1,2000,advance='no') spectral(ii,jj,kk,ik,ie)
!                        WRITE(2,2000,advance='no') emesh(ii,jj,kk,ik,ie)
!                        WRITE(3,2000,advance='no') kmeshz(ii,jj,kk,ik,ie)
!                    end do
!                     WRITE(1,"(A)",advance="yes") " "
!                     WRITE(2,"(A)",advance="yes") " "
!                     WRITE(3,"(A)",advance="yes") " "
                     end if
                end do
            end do
        end do
!    end do
!    write(*,*) dble(ik)/dble(nk)
end do
!2000 format(1f10.5)
if (verbose .ne. 0) then
    ! write the spectral function
    open(unit=1,file="spectral.dat",status="UNKNOWN",action="write",form="unformatted")
    open(unit=2,file="emesh.dat",status="UNKNOWN",action="write",form="unformatted")
    open(unit=3,file="kmeshz.dat",status="UNKNOWN",action="write",form="unformatted")
    write(1) spectral
    write(2) emesh
    write(3) kmeshz
    close(1)
    close(2)
    close(3)
end if

open(unit=1,file="datamesh.dat",status="UNKNOWN",action="write",form="unformatted")
write(1) datamesh
close(1)
open(unit=1,file="datamesh_proj.dat",status="UNKNOWN",action="write",form="unformatted")
write(1) datamesh_proj
close(1)

end subroutine

function weight_k(kk,gg,vv) result(y)

use view_struct
use config
use util

implicit none

real(kind=8) :: kk(3),gg(3)
complex(kind=8) :: vv(:)
integer(kind=4) :: ii,jj,i,j,k
complex(kind=8) :: y3(natm*3)
real(kind=8) :: y

y3 = 0.0d0
do ii=1,natm
    do jj = 1,3
        do i = 1,nx_sc
            do j = 1,ny_sc
                do k = 1,nz_sc
                    y3(3*(ii-1)+jj) = y3(3*(ii-1)+jj) + &
                    vv((i-1)*ny_sc*nz_sc*natm*3 &
                           +(j-1)*nz_sc*natm*3 &
                           +(k-1)*natm*3 &
                           +(ii-1)*3+jj) &
                           *exp(-i_imag*dot_product( &
                           pos_sc((i-1)*ny_sc*nz_sc*natm &
                           +(j-1)*nz_sc*natm &
                           +(k-1)*natm+ii,:)&
                            -pos(ii,:),gg+kk))
                end do
            end do
        end do
    end do
end do
y = dot_product(y3,y3)

end function

subroutine write_evs(i,evs,kpts)

use view_struct

implicit none

integer(kind=4) :: i,m,n,nmode,nindex
complex(kind=8) :: evs(:,:)
real(kind=8)    :: kpts(3)
character(len=8) :: fmt,x1

nmode = size(evs,2)
nindex = size(evs,1)


fmt = '(I6.6)'
write (x1,fmt) i
open(unit=1,file="displacements/"//trim(x1)//".dat",status="UNKNOWN",action="write")
do m = 1,nmode
    WRITE(1,*) real(evs(:,m))/sqrt(mass_ev(:))
end do
close(unit=1)

end subroutine
subroutine write_evs_avg(i,evs,kpts)

use view_struct

implicit none

integer(kind=4) :: i,m,n,nmode,nindex
complex(kind=8) :: evs(:,:)
real(kind=8)    :: kpts(3)
real(kind=8),allocatable   :: disp(:)
character(len=8) :: fmt,x1

nmode = size(evs,2)
nindex = size(evs,1)
allocate(disp(nmode))
disp = 0.0d0

do m = 1,nmode
    disp = disp+real(evs(:,m))/sqrt(mass_ev(:))/dble(nmode)
end do

fmt = '(I6.6)'
write (x1,fmt) i
open(unit=1,file="displacements/"//trim(x1)//".dat",status="UNKNOWN",action="write")
    WRITE(1,*) disp
close(unit=1)

end subroutine


end module
