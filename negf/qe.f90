module qe
    implicit none

    real(kind=8) :: qeph_cell(3,3)
    real(kind=8),allocatable :: qeph_pos(:,:)
    real(kind=8),allocatable :: qeph_mass(:)
    real(kind=8),allocatable :: qeph_fc(:,:,:,:,:,:,:)
    integer(kind=4)          :: qeph_nx,qeph_ny,qeph_nz
    real(kind=8)             :: qeph_rws(124,3),qeph_rd(124)


end module qe 
