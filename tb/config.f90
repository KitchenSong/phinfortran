module config
     
    implicit none

    integer(kind=4) :: k_trans_mesh(3)
    integer(kind=4) :: ne, nk, nkcore
    integer(kind=4),allocatable :: k_start_core(:)
    real(kind=8)    :: emin, emax
    real(kind=8),allocatable   :: Egrid(:,:)
    real(kind=8)    :: r_cutoff
    real(kind=8),allocatable :: k_grid(:,:)
    real(kind=8),allocatable :: k_grid_uc(:,:)
    real(kind=8),allocatable :: k_shift(:,:)
    integer(kind=4),allocatable :: k_uc_bz_idx(:)
    character(len=30) :: filename_input,&
                         left_primitive_cell_input,&
                         right_primitive_cell_input

    namelist /configlist/ k_trans_mesh, emin, emax, ne,&
                          r_cutoff,filename_input,&
                          left_primitive_cell_input,&
                          right_primitive_cell_input

    integer(kind=4) :: period_left, period_right
    integer(kind=4) :: buffer_left, buffer_right

    namelist /setup/ period_left, period_right, buffer_left,&
                     buffer_right

    real(kind=8)  :: eta,convergence
    namelist /convergelist/ eta, convergence

    integer(kind=4) :: n_bloch_x, n_bloch_y, circular_flag
    namelist /unfolding/ n_bloch_x, n_bloch_y, circular_flag

    integer(kind=4) :: path_mode,n_kpath_points 
    character(len=30) :: kpath_file,side
    real(kind=8),allocatable :: kpath_points(:,:)
    real(kind=8) :: az1
    namelist /kpath/ kpath_file, path_mode, side, az1

    character(len=20) :: output_dir
    integer(kind=4) :: verbosity
    namelist /saving/ output_dir,verbosity

    real(kind=8) :: uc_misori_x(3)
    real(kind=8) :: uc_misori_y(3)
    real(kind=8) :: uc_misori_z(3)
    real(kind=8) :: uc_misori(3,3)
    integer(kind=4) :: n_uc_misori1,n_uc_misori2,n_uc_misori3
    real(kind=8) :: theta
    namelist /misori/uc_misori_x,uc_misori_y,uc_misori_z,&
                    n_uc_misori1,n_uc_misori2,n_uc_misori3,&
                    theta


    integer(kind=4) :: myid, numprocs, error_unit, ierr


contains
    
    subroutine load_configure()
    
    implicit none

    include "mpif.h"

    integer(kind=4) :: i,j,n

    k_trans_mesh=0

    open(1,file="confige",status="old")
    read(1,nml=configlist)
    read(1,nml=setup)
    read(1,nml=convergelist)
    read(1,nml=unfolding)
    read(1,nml=kpath)
    read(1,nml=saving)
    read(1,nml=misori)
    close(1)

    ! path_mode = 1 : eigenvalue mode
    ! path_mode = 2 : transmission on full mesh
    ! path_mode = 3 : transmission along k path
    ! path_mode = 4 : surface band structure
    ! path_mode = 5 : eigenvalue mode for uc
    ! path_mode = 6 : surface dos mode

    ! generate k point grid for transport calculation
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
        ! number of k points on the same CPU core
        nk = k_trans_mesh(1) * k_trans_mesh(2)
    end if
    nkcore = nk/numprocs

    if ((modulo(nk,numprocs).ne.0).or.(nk.lt.numprocs)) then
        write(error_unit,*) "Error: numprocs is not integer times of nk"
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
    allocate(k_grid_uc(nk*n_bloch_x*n_bloch_y,3))
    allocate(k_shift(n_bloch_x*n_bloch_y,3))
 
    ! TO-DO: write a function that generates Wigner-Seitz
    ! cell in reciprocal space. This is required by folding
    ! algorithm. Those k points at boundary of Wigner-Seitz
    ! will be shared by neighboring Wigner-Seitz cells.
    k_grid = 0.0d0

    if (path_mode .eq. 2) then
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
         !do i = 1,k_trans_mesh(1)*k_trans_mesh(2)
         !    write(*,*) k_grid(i,1:2)
         !end do
         !stop

        if (nk .eq. 1) then
            k_grid(1,1) = 0.0
            k_grid(1,2) = 0.0
        end if

        k_grid_uc = 0.0d0
        allocate(k_uc_bz_idx(n_bloch_x*n_bloch_y*&
        k_trans_mesh(1)*k_trans_mesh(2)))
        k_uc_bz_idx = 0  
        n = 1
        do i = 1, k_trans_mesh(1)*n_bloch_x
            do j = 1, k_trans_mesh(2)*n_bloch_y
                k_grid_uc(n,1) = dble(i-1)/dble(k_trans_mesh(1))
                k_grid_uc(n,2) = dble(j-1)/dble(k_trans_mesh(2))
                k_uc_bz_idx(n) = (i-1)/n_bloch_x*&
                                  n_bloch_y+&
                                  (j-1)/n_bloch_y+1
                n = n + 1
            end do
        end do
        k_grid_uc(:,1) = k_grid_uc(:,1) - 0.5 * &
             dble(k_trans_mesh(1)*n_bloch_x-1)/dble(k_trans_mesh(1))
        k_grid_uc(:,2) = k_grid_uc(:,2) - 0.5 * &
             dble(k_trans_mesh(2)*n_bloch_y-1)/dble(k_trans_mesh(2))    
    elseif((path_mode.eq.3).or.(path_mode.eq.6)) then

        k_grid(:,:) = kpath_points(:,:) 
        k_grid(:,3) = 0.0d0

        allocate(k_uc_bz_idx(n_bloch_x*n_bloch_y*nk))
        k_uc_bz_idx = 0  
        n = 1
        do i = 1, nk*n_bloch_x*n_bloch_y
            k_uc_bz_idx(n) = n
            n = n + 1
        end do
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

    ! unit cell for misorientation
    uc_misori(1,:) = uc_misori_x(:)
    uc_misori(2,:) = uc_misori_y(:)
    uc_misori(3,:) = uc_misori_z(:) 
    
    ! load the energy grid
    if (emax.le.emin) then
        write(error_unit,*) "Error: wrong order of emin and emax!"
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_FINALIZE(ierr)
    end if

    allocate(Egrid(ne,nk))
    do i = 1, ne
        if (ne .eq. 1) then
            Egrid(i,:) = emax
        else
            Egrid(i,:) = emin + dble(i-1)/dble(ne-1)&
                       * (emax-emin)
        end if
    end do

    ! circular flag only works for supercell calculation
    if (circular_flag .eq. 1 .and. n_bloch_x.eq.1 .and. n_bloch_y.eq.1) then
        write(*,*) "circular flag only works for supercell calculation!"
        stop
    end if

    end subroutine load_configure

end module config
