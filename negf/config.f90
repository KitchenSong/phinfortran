module config

    implicit none

    integer(kind=4) :: n_fc_uc_x,n_fc_uc_y,n_fc_uc_z
    integer(kind=8) :: nfc
    integer(kind=4) :: nz
    integer(kind=4) :: k_trans_mesh(3)
    integer(kind=4) :: ne, nk, nkcore
    integer(kind=4),allocatable :: k_start_core(:)
    real(kind=8) :: emin, emax
    real(kind=8),allocatable   :: Egrid(:,:)
    real(kind=8) :: domega
    character(len=30) :: filename_input
    character(len=30) :: filename_input_uc
    character(len=30) :: filename_sandwich
    character(len=30) :: filename_left_lead
    character(len=30) :: filename_right_lead
    character(len=30) :: left_primitive_cell_input,&
                         right_primitive_cell_input
    real(kind=8) ::  r_cutoff
    real(kind=8),allocatable :: k_grid(:,:)
    real(kind=8),allocatable :: k_shift(:,:)
    integer(kind=4) :: nthick
    integer(kind=4) :: qe_fc,update_nbr
    character(len=30) :: flfrc

    ! setup buttiker probes
    integer(kind=4) :: buttiker,lifetime
    real(kind=8)    :: t_left,t_right,B_probe,T_probe

    
    namelist /configlist/k_trans_mesh,filename_input,&
                         filename_input_uc,&
                         filename_sandwich,&
                         filename_left_lead,&
                         filename_right_lead,&
                         n_fc_uc_x,n_fc_uc_y,n_fc_uc_z,&
                         nz,&
                         r_cutoff,emin,emax,ne,&
                         left_primitive_cell_input,&
                         right_primitive_cell_input,&
                         nthick,qe_fc,flfrc,buttiker,&
                         t_left,t_right,lifetime,&
                         update_nbr,B_probe,T_probe


    integer(kind=4) :: period_left, period_right
    integer(kind=4) :: buffer_left, buffer_right
    integer(kind=4) :: period_left_uc, period_right_uc
    integer(kind=4) :: buffer_left_uc, buffer_right_uc
    integer(kind=4) :: mass_only
    real(kind=8)    :: mass_ratio

    namelist /setup/ period_left, period_right, buffer_left,&
                     buffer_right, mass_only, mass_ratio

    real(kind=8)  :: eta,convergence, lambda_eta,  eta0
    namelist /convergelist/ eta, convergence, lambda_eta, eta0

    integer(kind=4) :: n_bloch_x, n_bloch_y, circular_flag
    namelist /unfolding/ n_bloch_x, n_bloch_y, circular_flag

    integer(kind=4) :: path_mode,n_kpath_points,crystal
    character(len=30) :: kpath_file,side
    real(kind=8),allocatable :: kpath_points(:,:)
    real(kind=8) :: az1
    namelist /kpath/ kpath_file, path_mode, side, az1, crystal

    character(len=20) :: output_dir
    integer(kind=4) :: verbosity
    namelist /saving/ output_dir,verbosity

    integer(kind=4) :: myid, numprocs, error_unit, ierr

contains
    
    subroutine load_configure()
    
    implicit none

    include "mpif.h"

    integer(kind=4) :: i,j,n

    mass_only = 0
    mass_ratio = 1.0d0
    k_trans_mesh=0
    eta0 = 1.0d-6
    eta = 1.0d-6
    buttiker=0
    lifetime=0
    update_nbr=0
    t_left = 300
    t_right = 300
    qe_fc = 1
    B_probe = 5.0d-20*(2.0d0*3.1415926d0*1.0d12)
    T_probe = 430.0d0

    open(1,file="config",status="old")
    read(1,nml=configlist)
    read(1,nml=setup)
    read(1,nml=convergelist)
    read(1,nml=unfolding)
    read(1,nml=kpath)
    read(1,nml=saving)
    close(1)
   
    ! path_mode = 1 : eigenvalue mode
    ! path_mode = 2 : transmission on full mesh
    ! path_mode = 3 : transmission along k path
    ! path_mode = 4 : slab band structure
    ! path_mode = 5 : eigenvalue of uc
    ! path_mode = 6 : sdos

    ! number of unit cells in supercells in fc calculation
    nfc = (2*n_fc_uc_x+1)*(2*n_fc_uc_y+1)*(2*n_fc_uc_z+1)
    if (path_mode .ne. 2) then
        open(1,file=trim(adjustl(kpath_file)),status="old")
        read(1,*) n_kpath_points
        allocate(kpath_points(n_kpath_points,3))
        kpath_points = 0.0d0
        do i = 1,n_kpath_points
            read(1,*) kpath_points(i,:)
        end do
        close(1)
        nk = n_kpath_points
    else
        ! generate k point grid for transport calculation
        nk = k_trans_mesh(1) * k_trans_mesh(2)
    end if

    ! number of k points on the same CPU core
    nkcore = nk/numprocs

    if ((modulo(nk,numprocs).ne.0).or.(nk.lt.numprocs)) then
        write(error_unit,*) "Error: numprocs is not integer times of nk"
        stop
    end if

    if (nkcore .eq. 1) then
        allocate(k_start_core(numprocs))
        do i = 1,numprocs
            k_start_core(i) = i
        end do
    else
        allocate(k_start_core(nkcore))
        do i = 1,nkcore
            k_start_core(i) = nkcore*(i-1)+1
        end do
    end if

    allocate(k_grid(nk,3))
    allocate(k_shift(n_bloch_x*n_bloch_y,3))

    if (path_mode.eq.2) then
        k_grid = 0.0d0
        n = 1
        do i = 1, k_trans_mesh(1)
            do j = 1, k_trans_mesh(2)
                k_grid(n,1) = dble(i-1)/dble(k_trans_mesh(1))
                k_grid(n,2) = dble(j-1)/dble(k_trans_mesh(2))
                n = n+1
            end do
        end do

        k_grid(:,1) = k_grid(:,1) - 0.5 * &
             dble(k_trans_mesh(1)-1)/dble(k_trans_mesh(1))
        k_grid(:,2) = k_grid(:,2) - 0.5 * &
             dble(k_trans_mesh(2)-1)/dble(k_trans_mesh(2))

        if (nk .eq. 1) then
            k_grid(1,1) = 0.0
            k_grid(1,2) = 0.0
        end if
        k_shift = 0.0d0
        n = 1
        do i = 1, n_bloch_x
            do j = 1, n_bloch_y
                k_shift(n,1) = dble(i-1)
                k_shift(n,2) = dble(j-1)
                n = n + 1
            end do
        end do
    else
        k_grid = kpath_points
        k_grid(:,3) = 0.0d0 ! Transport direction
    end if

    ! load the energy grid
    if (emax.le.emin) then
        write(error_unit,*) "Error: wrong order of emin and emax!"
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
    else if (emin .lt. 0.0d0) then
        write(error_unit,*) "Error: emin has to be larger than zero!"
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
    end if

        allocate(Egrid(ne,nk)) ! Omega^2
        if (ne .eq. 1) then
            Egrid(1,:) = emax**2
        else
            do i = 1, ne
                Egrid(i,:) = (dble(i-1)/dble(ne-1)&
                           * (emax))**2
            end do
            ! replace 0 with emin**2
            Egrid(1,:) = emin**2
        end if
        if (ne .eq. 1) then
            domega = emax
        else
            domega = emax/dble(ne-1)
        end if
        

    end subroutine load_configure

end module
