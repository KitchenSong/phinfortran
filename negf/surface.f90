module surface

implicit none

integer(kind=4) :: left_buffer_orb, right_buffer_orb
integer(kind=4) :: nl, nr, n_buffer_l, n_buffer_r
real(kind=8)    :: a_z_L,a_z_R

contains

    subroutine get_periodic_lead_orb_num()

    use config ,  only : buffer_left, period_left,&
                         buffer_right, period_right,&
                         n_fc_uc_z
    use input  ,  only : natoms,positions,&
                         latvec_left,&
                         latvec_right,&
                         latvec_uc_l,&
                         latvec_uc_r,&
                         latvec_uc_l_sc,&
                         latvec_uc_r_sc,&
                         reci_uc_l,reci_uc_r

    use util

    implicit none

    include "mpif.h"

    integer(kind=4) :: i
    real(kind=8) :: pos1,pos2

    pos1 = 1.0d8
    pos2 = 1.0d8

    n_buffer_l = 3*buffer_left
    n_buffer_r = 3*buffer_right

    nl = 3*period_left
    nr = 3*period_right

    ! get the periodic length in z direction for left lead
    do i = buffer_left+1, buffer_left+period_left
        if (pos1 .gt. positions(i,3)) then
            pos1 = positions(i,3)
        end if
        if (pos2 .gt. positions(period_left+i,3)) then
            pos2 = positions(period_left+i,3)
        end if
    end do
    a_z_L = pos2 - pos1

    latvec_uc_l(3,3) = a_z_L
    latvec_left(3,3) = a_z_L
    latvec_uc_l_sc(3,3) = a_z_L * dble(2*n_fc_uc_z+1)

    reci_uc_l = transpose(2.0*pi*inv_real(latvec_uc_l))

    ! get the periodic length in z direction for right lead
    pos1 = -1.0d8
    pos2 = -1.0d8
    do i = natoms - buffer_right-2*period_right+1,&
        natoms-buffer_right-period_right
        if (pos1 .lt. positions(i,3)) then
            pos1 = positions(i,3)
        end if
        if (pos2 .lt. positions(period_left+i,3)) then
            pos2 = positions(period_left+i,3)
        end if
    end do
    a_z_R = pos2 - pos1

    latvec_uc_r(3,3) = a_z_R
    latvec_right(3,3) = a_z_R
    latvec_uc_r_sc(3,3) = a_z_R * dble(2*n_fc_uc_z+1)
    reci_uc_r = transpose(2.0*pi*inv_real(latvec_uc_r))

    end subroutine get_periodic_lead_orb_num

    subroutine surface_green_mpi(ie,ik_core,Ham,norb,&
                    gl_adv,gr_ret,gl_plus,gr_minus,&
                    gml_adv,gmr_ret,&
                    sl_adv,sr_ret,E_k)
    
    use config

    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: ie,ik_core
    complex(kind=8),intent(in) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(out) :: gl_adv(numprocs,nl,nl)
    complex(kind=8),intent(out) :: gr_ret(numprocs,nr,nr)
    complex(kind=8),intent(out) :: gl_plus(numprocs,nl,nl)
    complex(kind=8),intent(out) :: gr_minus(numprocs,nr,nr)
    complex(kind=8),intent(out) :: gml_adv(numprocs,nl,nl)
    complex(kind=8),intent(out) :: gmr_ret(numprocs,nr,nr)
    complex(kind=8),intent(out) :: sl_adv(numprocs,nl,nl)
    complex(kind=8),intent(out) :: sr_ret(numprocs,nr,nr)

    integer(kind=4),intent(in) :: norb
    real(kind=8),intent(in) :: E_k(ne,nk)
    integer(kind=4) :: mm, nn

    mm = myid + 1
    nn = k_start_core(mm)+ik_core-1
 
    call surface_green(mm,Ham,norb,&
               gl_adv(mm,:,:),gr_ret(mm,:,:),&
               gl_plus(mm,:,:),gr_minus(mm,:,:),&
               gml_adv(mm,:,:),gmr_ret(mm,:,:),&
               sl_adv(mm,:,:),sr_ret(mm,:,:),E_k(ie,nn))

    end subroutine surface_green_mpi

    subroutine surface_green(mm,Ham,norb,&
                        gl_adv,gr_ret,gl_plus,gr_minus,&
                        gml_adv,gmr_ret,sl_adv,sr_ret,E)
    use input , only : reci_uc_l,latvec_uc_l,&
                       reci_uc_r,latvec_uc_r
    use config
    use util
    use circular

    implicit none

    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(in) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(out) :: gl_adv(nl,nl)
    complex(kind=8),intent(out) :: gr_ret(nr,nr)
    complex(kind=8),intent(out) :: gl_plus(nl,nl)
    complex(kind=8),intent(out) :: gr_minus(nr,nr)
    complex(kind=8),intent(out) :: gml_adv(nl,nl)
    complex(kind=8),intent(out) :: gmr_ret(nr,nr)
    complex(kind=8),intent(out) :: sl_adv(nl,nl)
    complex(kind=8),intent(out) :: sr_ret(nr,nr)
    integer(kind=4),intent(in) :: mm
    real(kind=8),intent(in) :: E
    integer(kind=4) :: i,j

    complex(kind=8) :: Ham_l_i_i(nl,nl)
    complex(kind=8) :: Ham_l_o_i(nl,nl)
    complex(kind=8) :: Ham_l_i_o(nl,nl)
    complex(kind=8) :: Ham_l_o_o(nl,nl)
    complex(kind=8) :: Ham_r_i_i(nr,nr)
    complex(kind=8) :: Ham_r_o_i(nr,nr)
    complex(kind=8) :: Ham_r_i_o(nr,nr)
    complex(kind=8) :: Ham_r_o_o(nr,nr)
    complex(kind=8) :: gloo_adv(nl,nl),groo_ret(nr,nr)
    complex(kind=8) :: gloo_plus(nl,nl),groo_minus(nr,nr)
    complex(kind=8) :: sigmal_adv(nl,nl)
    complex(kind=8) :: sigmar_ret(nr,nr)
    complex(kind=8) :: gammal_adv(nl,nl)
    complex(kind=8) :: gammar_ret(nr,nr)


    real(kind=8) :: maxdiff
    complex(kind=8),allocatable :: Es0(:,:),e1(:,:),es(:,:)
    complex(kind=8),allocatable :: alpha(:,:),beta(:,:),ive(:,:)

    integer(kind=4) :: na, nac

    na = nl/n_bloch_y
    nac = na/n_bloch_x

    ! J. Phys. F: Met. Phys. 15 (1985) 851-858. Eq. 11

    Ham_l_i_i = 0.0d0
    Ham_l_o_i = 0.0d0
    Ham_l_i_o = 0.0d0
    Ham_l_o_o = 0.0d0
    Ham_r_i_i = 0.0d0
    Ham_r_o_i = 0.0d0
    Ham_r_i_o = 0.0d0
    Ham_r_o_o = 0.0d0

    Ham_l_i_i(:,:) = Ham(mm,n_buffer_l+1:n_buffer_l+nl,&
    n_buffer_l+1:n_buffer_l+nl)
    Ham_l_o_i(:,:) = Ham(mm,n_buffer_l+1:n_buffer_l+nl,&
    n_buffer_l+nl+1:n_buffer_l+2*nl)
    Ham_l_i_o(:,:) = Ham(mm,n_buffer_l+nl+1:n_buffer_l+2*nl,&
    n_buffer_l+1:n_buffer_l+nl)
   
    Ham_l_o_o(:,:) = Ham(mm,n_buffer_l+nl+1:n_buffer_l+2*nl,&
    n_buffer_l+nl+1:n_buffer_l+2*nl)
    Ham_r_i_i(:,:) = Ham(mm,norb-n_buffer_r-nr+1:&
    norb-n_buffer_r,norb-n_buffer_r-nr+1:&
    norb-n_buffer_r)
    Ham_r_o_i(:,:) = Ham(mm,norb-n_buffer_r-2*nr+1:&
    norb-n_buffer_r-nr,norb-n_buffer_r-nr+1:&
    norb-n_buffer_r)
    Ham_r_i_o(:,:) = Ham(mm,norb-n_buffer_r-nr+1:&
    norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
    norb-n_buffer_r-nr)
    Ham_r_o_o(:,:) = Ham(mm,norb-n_buffer_r-2*nr+1:&
    norb-n_buffer_r-nr,norb-n_buffer_r-2*nr+1:&
    norb-n_buffer_r-nr)
!    write(*,*)'a'
!    write(*,*) Ham_l_o_i - transpose(dconjg(Ham_l_i_o))
!    write(*,*)'b'

!             do i = 1,12
!                 write(*,*) Ham_l_o_o(i,:)
!             end do
!             write(*,*) 'a'
    ! Left surface green function
    allocate(Es0(nl,nl),e1(nl,nl),es(nl,nl),&
    alpha(nl,nl),beta(nl,nl),ive(nl,nl))
    Es0 = (E-(E*eta+eta0)*i_imag)*eyemat(nl)
    e1(:,:) = Ham_l_o_o(:,:)
    es(:,:) = Ham_l_o_o(:,:)
    !write(*,*) 'df', maxval(abs(Ham_l_o_o(:,:)-Ham_l_i_i(:,:)))
    !write(*,*) 'max',minval(abs(Ham_l_o_o(:,:)))
    alpha(:,:) = Ham_l_i_o(:,:)
    beta(:,:) =  dconjg(transpose(Ham_l_i_o(:,:))) ! Ham_l_o_i(:,:)
    maxdiff = maxval(abs(alpha))

    if (circular_flag.eq.0) then
 !       if (maxdiff .gt. 15.0) then
 !       write(*,*) 'maxdiff'
 !       end if

        do while( maxdiff .gt. convergence)
            ive = inv(Es0-e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                    matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
        gloo_adv = inv(Es0-es)
    else
        gloo_adv =circular_g(Es0,e1,alpha,&
                 beta,convergence,&
                 nl,n_bloch_y,n_bloch_x,reci_uc_l(2,:),latvec_uc_l(2,:),&
                        reci_uc_l(1,:),latvec_uc_l(1,:))
    end if
    sigmal_adv =matmul(matmul(Ham_l_i_o,gloo_adv),&
                                    Ham_l_o_i)
    gammal_adv = i_imag*(sigmal_adv-&
    transpose(dconjg(sigmal_adv))) ! gammal_ret = -gammal_adv defined here

    ! Left surface green function, reflected
    Es0 = (E+(E*eta+eta0)*i_imag)*eyemat(nl)
    e1(:,:) = Ham_l_o_o(:,:)
    es(:,:) = Ham_l_o_o(:,:)
    alpha(:,:) = Ham_l_o_i(:,:)
    beta(:,:) = dconjg(transpose(Ham_l_o_i(:,:)))! Ham_l_i_o(:,:)

    if (circular_flag.eq.0) then
        maxdiff = maxval(abs(alpha))
        do while( maxdiff .gt. convergence)
            ive = inv(Es0-e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                    matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
!        write(*,*) maxdiff
        gloo_plus = inv(Es0-es)
    else
        gloo_plus =circular_g(Es0,e1,alpha,&
                 beta,convergence,&
                 nl,n_bloch_y,n_bloch_x,reci_uc_l(2,:),latvec_uc_l(2,:),&
                        reci_uc_l(1,:),latvec_uc_l(1,:))
    end if
    ! Right surface green function
    deallocate(Es0,e1,es,alpha,beta,ive)
    allocate(Es0(nr,nr),e1(nr,nr),es(nr,nr),&
    alpha(nr,nr),beta(nr,nr),ive(nr,nr))

    Es0 = (E+(E*eta+eta0)*i_imag)*eyemat(nr)
    alpha = Ham_r_o_i
    beta =  dconjg(transpose(Ham_r_o_i))! Ham_r_i_o
    e1 = Ham_r_o_o
    es = Ham_r_o_o
    maxdiff = maxval(abs(alpha))

    if (circular_flag.eq.0) then
        do while( maxdiff .gt. convergence)
            ive = inv(Es0 - e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                     matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
        groo_ret = inv(Es0-es)
    else
        groo_ret =circular_g(Es0,e1,alpha,&
             beta,convergence,&
             nr,n_bloch_y,n_bloch_x,reci_uc_r(2,:),latvec_uc_r(2,:),&
             reci_uc_r(1,:),latvec_uc_r(1,:))
    end if

    sigmar_ret =matmul(matmul(Ham_r_o_i,groo_ret),&
                                 Ham_r_i_o)

    gammar_ret = i_imag*(sigmar_ret - transpose(dconjg(sigmar_ret)))

    ! Right surface green function, reflected
    Es0 = (E-(E*eta+eta0)*i_imag)*eyemat(nr)
    alpha = Ham_r_i_o
    beta =  dconjg(transpose(Ham_r_i_o)) ! Ham_r_o_i
    e1 = Ham_r_o_o
    es = Ham_r_o_o
    maxdiff = maxval(abs(alpha))

    if (circular_flag.eq.0) then
        do while( maxdiff .gt. convergence)
            ive = inv(Es0 - e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                     matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
        groo_minus = inv(Es0-es)
    else
        groo_minus =circular_g(Es0,e1,alpha,&
                 beta,convergence,&
                 nr,n_bloch_y,n_bloch_x,reci_uc_r(2,:),latvec_uc_r(2,:),&
                 reci_uc_r(1,:),latvec_uc_r(1,:))
    end if

    gl_adv(:,:) = gloo_adv(:,:)
    gl_plus(:,:) = gloo_plus(:,:)
    gr_ret(:,:) = groo_ret(:,:)
    gr_minus(:,:) = groo_minus(:,:)
    gml_adv(:,:) = gammal_adv(:,:)
    gmr_ret(:,:) = gammar_ret(:,:)
    sl_adv(:,:) = sigmal_adv(:,:)
    sr_ret(:,:) = sigmar_ret(:,:)

    end subroutine surface_green

    subroutine surface_green_tdb_mpi(ie,ik_core,&
                    gl_adv,gr_ret,gl_plus,gr_minus,&
                    gml_adv,gmr_ret,&
                    sl_adv,sr_ret,E_k,&
                    Haml,Hamr,&
                    Hamlv,Hamrv)
    
    use config

    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: ie,ik_core
    complex(kind=8),intent(out) :: gl_adv(numprocs,nl,nl)
    complex(kind=8),intent(out) :: gr_ret(numprocs,nr,nr)
    complex(kind=8),intent(out) :: gl_plus(numprocs,nl,nl)
    complex(kind=8),intent(out) :: gr_minus(numprocs,nr,nr)
    complex(kind=8),intent(out) :: gml_adv(numprocs,nl,nl)
    complex(kind=8),intent(out) :: gmr_ret(numprocs,nr,nr)
    complex(kind=8),intent(out) :: sl_adv(numprocs,nl,nl)
    complex(kind=8),intent(out) :: sr_ret(numprocs,nr,nr)
    complex(kind=8),intent(in) :: Haml(:,:,:),Hamlv(:,:,:)
    complex(kind=8),intent(in) :: Hamr(:,:,:),Hamrv(:,:,:)

    real(kind=8),intent(in) :: E_k(ne,nk)
    integer(kind=4) :: mm, nn

    mm = myid + 1
    nn = k_start_core(mm)+ik_core-1
 
    call surface_green_tdb(mm,&
               gl_adv(mm,:,:),gr_ret(mm,:,:),&
               gl_plus(mm,:,:),gr_minus(mm,:,:),&
               gml_adv(mm,:,:),gmr_ret(mm,:,:),&
               sl_adv(mm,:,:),sr_ret(mm,:,:),E_k(ie,nn),&
               Haml(mm,:,:),Hamr(mm,:,:),&
               Hamlv(mm,:,:),Hamrv(mm,:,:))

    end subroutine surface_green_tdb_mpi

    subroutine surface_green_tdb(mm,&
                        gl_adv,gr_ret,gl_plus,gr_minus,&
                        gml_adv,gmr_ret,sl_adv,sr_ret,E,&
                        Haml,Hamr,&
                        Hamlv,Hamrv)
    use input , only : reci_uc_l,latvec_uc_l,&
                       reci_uc_r,latvec_uc_r
    use config
    use util
    use circular

    implicit none

    complex(kind=8),intent(in) :: Haml(:,:),Hamlv(:,:)
    complex(kind=8),intent(in) :: Hamr(:,:),Hamrv(:,:)
    complex(kind=8),intent(out) :: gl_adv(nl,nl)
    complex(kind=8),intent(out) :: gr_ret(nr,nr)
    complex(kind=8),intent(out) :: gl_plus(nl,nl)
    complex(kind=8),intent(out) :: gr_minus(nr,nr)
    complex(kind=8),intent(out) :: gml_adv(nl,nl)
    complex(kind=8),intent(out) :: gmr_ret(nr,nr)
    complex(kind=8),intent(out) :: sl_adv(nl,nl)
    complex(kind=8),intent(out) :: sr_ret(nr,nr)
    integer(kind=4),intent(in) :: mm
    real(kind=8),intent(in) :: E
    integer(kind=4) :: i,j

    complex(kind=8) :: Ham_l_i_i(nl,nl)
    complex(kind=8) :: Ham_l_o_i(nl,nl)
    complex(kind=8) :: Ham_l_i_o(nl,nl)
    complex(kind=8) :: Ham_l_o_o(nl,nl)
    complex(kind=8) :: Ham_r_i_i(nr,nr)
    complex(kind=8) :: Ham_r_o_i(nr,nr)
    complex(kind=8) :: Ham_r_i_o(nr,nr)
    complex(kind=8) :: Ham_r_o_o(nr,nr)
    complex(kind=8) :: gloo_adv(nl,nl),groo_ret(nr,nr)
    complex(kind=8) :: gloo_plus(nl,nl),groo_minus(nr,nr)
    complex(kind=8) :: sigmal_adv(nl,nl)
    complex(kind=8) :: sigmar_ret(nr,nr)
    complex(kind=8) :: gammal_adv(nl,nl)
    complex(kind=8) :: gammar_ret(nr,nr)


    real(kind=8) :: maxdiff
    complex(kind=8),allocatable :: Es0(:,:),e1(:,:),es(:,:)
    complex(kind=8),allocatable :: alpha(:,:),beta(:,:),ive(:,:)

    integer(kind=4) :: na, nac

    na = nl/n_bloch_y
    nac = na/n_bloch_x

    ! J. Phys. F: Met. Phys. 15 (1985) 851-858. Eq. 11

    Ham_l_i_i = 0.0d0
    Ham_l_o_i = 0.0d0
    Ham_l_i_o = 0.0d0
    Ham_l_o_o = 0.0d0
    Ham_r_i_i = 0.0d0
    Ham_r_o_i = 0.0d0
    Ham_r_i_o = 0.0d0
    Ham_r_o_o = 0.0d0

    Ham_l_i_i(:,:) = Haml
    Ham_l_o_i(:,:) = Hamlv
    !Ham(mm,n_buffer_l+1:n_buffer_l+nl,&
    !n_buffer_l+nl+1:n_buffer_l+2*nl)
    Ham_l_i_o(:,:) = transpose(dconjg(Hamlv))
    !Ham(mm,n_buffer_l+nl+1:n_buffer_l+2*nl,&
    !n_buffer_l+1:n_buffer_l+nl)
   
    Ham_l_o_o(:,:) = Haml
    !Ham(mm,n_buffer_l+nl+1:n_buffer_l+2*nl,&
    !n_buffer_l+nl+1:n_buffer_l+2*nl)
    Ham_r_i_i(:,:) = Hamr
    !Ham(mm,norb-n_buffer_r-nr+1:&
    !norb-n_buffer_r,norb-n_buffer_r-nr+1:&
    !norb-n_buffer_r)
    Ham_r_o_i(:,:) = transpose(dconjg(Hamrv))
    !Ham(mm,norb-n_buffer_r-2*nr+1:&
    !norb-n_buffer_r-nr,norb-n_buffer_r-nr+1:&
    !norb-n_buffer_r)
    Ham_r_i_o(:,:) = Hamrv 
    !Ham(mm,norb-n_buffer_r-nr+1:&
    !norb-n_buffer_r,norb-n_buffer_r-2*nr+1:&
    !norb-n_buffer_r-nr)
    Ham_r_o_o(:,:) = Hamr
    !Ham(mm,norb-n_buffer_r-2*nr+1:&
    !norb-n_buffer_r-nr,norb-n_buffer_r-2*nr+1:&
    !norb-n_buffer_r-nr)

    ! Left surface green function
    allocate(Es0(nl,nl),e1(nl,nl),es(nl,nl),&
    alpha(nl,nl),beta(nl,nl),ive(nl,nl))
    Es0 = (E-(E*eta+eta0)*i_imag)*eyemat(nl)
    e1(:,:) = Ham_l_o_o(:,:)
    es(:,:) = Ham_l_o_o(:,:)
    alpha(:,:) = Ham_l_i_o(:,:)
    beta(:,:) =  dconjg(transpose(Ham_l_i_o(:,:))) ! Ham_l_o_i(:,:)
    maxdiff = maxval(abs(alpha))

    if (circular_flag.eq.0) then

        do while( maxdiff .gt. convergence)
            ive = inv(Es0-e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                    matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
        gloo_adv = inv(Es0-es)
    else
        gloo_adv =circular_g(Es0,e1,alpha,&
                 beta,convergence,&
                 nl,n_bloch_y,n_bloch_x,reci_uc_l(2,:),latvec_uc_l(2,:),&
                        reci_uc_l(1,:),latvec_uc_l(1,:))
    end if
    sigmal_adv =matmul(matmul(Ham_l_i_o,gloo_adv),&
                                    Ham_l_o_i)
    gammal_adv = i_imag*(sigmal_adv-&
    transpose(dconjg(sigmal_adv))) ! gammal_ret = -gammal_adv defined here

    ! Left surface green function, reflected
    Es0 = (E+(E*eta+eta0)*i_imag)*eyemat(nl)
    e1(:,:) = Ham_l_o_o(:,:)
    es(:,:) = Ham_l_o_o(:,:)
    alpha(:,:) = Ham_l_o_i(:,:)
    beta(:,:) = dconjg(transpose(Ham_l_o_i(:,:)))! Ham_l_i_o(:,:)

    if (circular_flag.eq.0) then
        maxdiff = maxval(abs(alpha))
        do while( maxdiff .gt. convergence)
            ive = inv(Es0-e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                    matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
!        write(*,*) maxdiff
        gloo_plus = inv(Es0-es)
    else
        gloo_plus =circular_g(Es0,e1,alpha,&
                 beta,convergence,&
                 nl,n_bloch_y,n_bloch_x,reci_uc_l(2,:),latvec_uc_l(2,:),&
                        reci_uc_l(1,:),latvec_uc_l(1,:))
    end if
    ! Right surface green function
    deallocate(Es0,e1,es,alpha,beta,ive)
    allocate(Es0(nr,nr),e1(nr,nr),es(nr,nr),&
    alpha(nr,nr),beta(nr,nr),ive(nr,nr))

    Es0 = (E+(E*eta+eta0)*i_imag)*eyemat(nr)
    alpha = Ham_r_o_i
    beta =  dconjg(transpose(Ham_r_o_i))! Ham_r_i_o
    e1 = Ham_r_o_o
    es = Ham_r_o_o
    maxdiff = maxval(abs(alpha))

    if (circular_flag.eq.0) then
        do while( maxdiff .gt. convergence)
            ive = inv(Es0 - e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                     matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
        groo_ret = inv(Es0-es)
    else
        groo_ret =circular_g(Es0,e1,alpha,&
             beta,convergence,&
             nr,n_bloch_y,n_bloch_x,reci_uc_r(2,:),latvec_uc_r(2,:),&
             reci_uc_r(1,:),latvec_uc_r(1,:))
    end if

    sigmar_ret =matmul(matmul(Ham_r_o_i,groo_ret),&
                                 Ham_r_i_o)

    gammar_ret = i_imag*(sigmar_ret - transpose(dconjg(sigmar_ret)))

    ! Right surface green function, reflected
    Es0 = (E-(E*eta+eta0)*i_imag)*eyemat(nr)
    alpha = Ham_r_i_o
    beta =  dconjg(transpose(Ham_r_i_o)) ! Ham_r_o_i
    e1 = Ham_r_o_o
    es = Ham_r_o_o
    maxdiff = maxval(abs(alpha))

    if (circular_flag.eq.0) then
        do while( maxdiff .gt. convergence)
            ive = inv(Es0 - e1)
            e1 = e1 + matmul(matmul(alpha,ive),beta)+&
                     matmul(matmul(beta,ive),alpha)
            es = es + matmul(matmul(alpha,ive),beta)
            alpha = matmul(matmul(alpha,ive),alpha)
            beta = matmul(matmul(beta,ive),beta)
            maxdiff = maxval(abs(alpha))
        end do
        groo_minus = inv(Es0-es)
    else
        groo_minus =circular_g(Es0,e1,alpha,&
                 beta,convergence,&
                 nr,n_bloch_y,n_bloch_x,reci_uc_r(2,:),latvec_uc_r(2,:),&
                 reci_uc_r(1,:),latvec_uc_r(1,:))
    end if

    gl_adv(:,:) = gloo_adv(:,:)
    gl_plus(:,:) = gloo_plus(:,:)
    gr_ret(:,:) = groo_ret(:,:)
    gr_minus(:,:) = groo_minus(:,:)
    gml_adv(:,:) = gammal_adv(:,:)
    gmr_ret(:,:) = gammar_ret(:,:)
    sl_adv(:,:) = sigmal_adv(:,:)
    sr_ret(:,:) = sigmar_ret(:,:)

    end subroutine surface_green_tdb


    subroutine surface_ldos_mpi(ie,ik_core,&
                                gl_adv,gr_ret,&
                                sdos_l,sdos_r,E_k)

    use config
    use util

    implicit none

    include "mpif.h"

    integer(kind=4),intent(in) :: ie,ik_core
    complex(kind=8),intent(in) :: gl_adv(numprocs,nl,nl)
    complex(kind=8),intent(in) :: gr_ret(numprocs,nr,nr)
    real(kind=8) :: sdos_l(:,:),sdos_r(:,:)
    real(kind=8),intent(in) :: E_k(ne,nk)
    integer(kind=4) :: mm, nn
    
    mm = myid + 1
    nn = k_start_core(mm)+ik_core-1
    sdos_l(ie,nn) =  2.0*sqrt(E_k(ie,nn))*real(aimag(trace(gl_adv(mm,:,:)))/pi)
    sdos_r(ie,nn) = -2.0*sqrt(E_k(ie,nn))*real(aimag(trace(gr_ret(mm,:,:)))/pi)
    end subroutine surface_ldos_mpi

end module surface
