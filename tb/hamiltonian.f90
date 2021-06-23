module hamiltonian

implicit none

contains
    function get_principal_number(orb)
    implicit none

    character(len=*),intent(in) :: orb
    integer(kind=4) :: get_principal_number
    character(len=8)   :: temp

    temp = trim(adjustl(orb))

    if (temp(1:1) .eq. 's') then
        get_principal_number = 1
    else if (temp(1:1) .eq. 'p') then
        get_principal_number = 2
    else if (temp(1:1) .eq. 'd') then
        get_principal_number = 3
    end if
        
    if (temp(1:2) .eq. 's*') then
        get_principal_number = 4
    end if
    end function get_principal_number
 
    function V_term(i,j,n1,n2,vtype,dist_ij,V,Vn,r_d,nspecies)

    implicit none
    ! physical constants
    integer(kind=4),intent(in) :: i,j,vtype
    integer(kind=4),intent(in) :: n1,n2
    integer(kind=4)            :: nspecies 
    real(kind=8),intent(in) :: r_d(nspecies,nspecies)
    real(kind=8),intent(in) :: V(nspecies,nspecies,4,4,4)
    real(kind=8),intent(in) :: Vn(nspecies,nspecies,4,4,4)
    real(kind=8),intent(in)    :: dist_ij
    real(kind=8),parameter :: hbar = 1.0545718d-34
    real(kind=8),parameter :: eV = 1.60217662d-19
    real(kind=8),parameter :: me = 9.10938356d-31
    real(kind=8),parameter :: ang = 1.0d-10
    real(kind=8)           :: V_term
    real(kind=8) :: d0,np

    ! the power law expression can be found in PhysRevB.79.245201

    V_term = 0.0d0
   
    d0 = r_d(i,j) ! equilibrium distance
    np = Vn(i,j,n1,n2,vtype)

    V_term = V(i,j,n1,n2,vtype)*(d0/dist_ij/ang)**np

    end function

    function onsite(i,j,orb1,orb2,V,nspecies)
    integer(kind=4),intent(in) :: nspecies, i, j
    character(len=*) :: orb1,orb2
    real(kind=8),intent(in) :: V(nspecies,nspecies,4,4,4)
    integer(kind=4) :: n1,n2
    real(kind=8)  :: onsite

    if (orb1 .eq. orb2) then
        n1 = get_principal_number(orb1)
        n2 = get_principal_number(orb2)
        onsite = V(i,j,n1,n2,1)
    else
        onsite = 0.0d0
    end if

    end function onsite

    function V_sk(i,j,orb1,orb2,&
                 r_ij,V,Vn,r_d,&
                 nspecies)
    
    ! the Slater-Koster transformation
    implicit none

    integer(kind=4),intent(in) :: i,j
    character(len=*),intent(in) :: orb1,orb2
    integer(kind=4) :: n1,n2
    real(kind=8),intent(in) :: r_ij(3)
    integer(kind=4),intent(in) :: nspecies
    real(kind=8),intent(in) :: r_d(nspecies,nspecies)
    real(kind=8),intent(in) :: V(nspecies,nspecies,4,4,4)
    real(kind=8),intent(in) :: Vn(nspecies,nspecies,4,4,4)
    real(kind=8)  :: dist_ij
    real(kind=8) :: V_sk    
    real(kind=8) :: l, m, n
    character(len=:),allocatable :: int_type
    integer(kind=4) :: clen

    dist_ij = sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)
    l = r_ij(1)/dist_ij
    m = r_ij(2)/dist_ij 
    n = r_ij(3)/dist_ij

    V_sk = 0.0d0

    clen = len( trim(adjustl(orb1))//trim(adjustl(orb2)))
    allocate(character(len=clen) :: int_type)
    int_type = trim(adjustl(orb1))//trim(adjustl(orb2))
    select case(int_type)
    case("ss")  
                V_sk=V_term(i,j,1,1,2,dist_ij,V,Vn,r_d,nspecies)
       ! l = 0 with l = 1
    case("spx") 
                V_sk=l * V_term(i,j,1,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("pxs") 
                V_sk=-l * V_term(j,i,1,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("spy") 
                V_sk=  m * V_term(i,j,1,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("pys") 
                V_sk= -m * V_term(j,i,1,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("spz") 
                V_sk=  n * V_term(i,j,1,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("pzs") 
                V_sk=-n * V_term(j,i,1,2,2,dist_ij,V,Vn,r_d,nspecies)
       ! l = 1 with l = 1
    case("pxpx")
                V_sk= (l**2) * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) + &
                (1 - l**2) * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pypy")
                V_sk= (m**2) * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) + &
                (1 - m**2) * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzpz")
                V_sk= (n**2) * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) + &
                (1 - n**2) * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pxpy")
                V_sk=l * m * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * m * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pypx")
                V_sk= l * m * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * m * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pypz")
                V_sk= m * n * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * n * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzpy")
                V_sk=m * n * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * n * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pxpz")
                V_sk=l * n * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * n * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzpx")
                V_sk=l * n * V_term(i,j,2,2,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * n * V_term(i,j,2,2,3,dist_ij,V,Vn,r_d,nspecies)
       ! l = 0 with l = 2
    case("sdxy") 
                V_sk=  sqrt(3.0) * (l*m) * V_term(i,j,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dxys") 
                V_sk= sqrt(3.0) * (l*m) * V_term(j,i,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("sdyz") 
                V_sk= sqrt(3.0) * (m*n) * V_term(i,j,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dyzs") 
                V_sk= sqrt(3.0) * (m*n) * V_term(j,i,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("sdxz") 
                V_sk= sqrt(3.0) * (l*n) * V_term(i,j,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dxzs") 
                V_sk= sqrt(3.0) * (l*n) * V_term(j,i,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("sdx2-y2")
                V_sk= (sqrt(3.0)/2) * (l**2 - m**2) * V_term(i,j,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2s")
                V_sk=(sqrt(3.0)/2) * (l**2 - m**2) * V_term(j,i,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("sdz2") 
                V_sk= (n**2 - (l**2 + m**2)/2) * V_term(i,j,1,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dz2s") 
                V_sk= (n**2 - (l**2 + m**2)/2) * V_term(j,i,1,3,2,dist_ij,V,Vn,r_d,nspecies)
       ! l = 1 with l = 2
    case("pxdxy")
                V_sk= sqrt(3.0) * ((l**2)*m) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * (1 - 2*(l**2)) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dxypx")
                V_sk=  -sqrt(3.0) * ((l**2)*m) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * (1 - 2*(l**2)) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pxdyz")
                V_sk= sqrt(3.0) * (l*m*n) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                2 * (l*m*n) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dyzpx")
                V_sk= -sqrt(3.0) * (l*m*n) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                2 * (l*m*n) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pxdxz")
                V_sk= sqrt(3.0) * ((l**2)*n) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                n * (1 - 2*(l**2)) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dxzpx")
                V_sk=-sqrt(3.0) * ((l**2)*n) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                n * (1 - 2*(l**2)) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)

    case("pydxy")
                V_sk= sqrt(3.0) * ((m**2)*l) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * (1 - 2*(m**2)) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dxypy")
                V_sk=  -sqrt(3.0) * ((m**2)*l) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * (1 - 2*(m**2)) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pydyz")
                V_sk=  sqrt(3.0) * ((m**2)*n) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                n * (1 - 2*(m**2)) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dyzpy")
                V_sk=  -sqrt(3.0) * ((m**2)*n) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                n * (1 - 2*(m**2)) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pydxz")
                V_sk=  sqrt(3.0) * (l*m*n) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                2 * (l*m*n) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dxzpy")
                V_sk=  -sqrt(3.0) * (l*m*n) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                2 * (l*m*n) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)

    case("pzdxy")
                V_sk=   sqrt(3.0) * (l*m*n) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                2 * (l*m*n) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dxypz")
                V_sk= -sqrt(3.0) * (l*m*n) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                2 * (l*m*n) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzdyz")
                V_sk= sqrt(3.0) * ((n**2)*m) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * (1 - 2*(n**2)) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dyzpz")
                V_sk= -sqrt(3.0) * ((n**2)*m) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * (1 - 2*(n**2)) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzdxz")
                V_sk= sqrt(3.0) * ((n**2)*l) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * (1 - 2*(n**2)) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dxzpz")
                V_sk= -sqrt(3.0) * ((n**2)*l) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * (1 - 2*(n**2)) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)

    case("pxdx2-y2")
                V_sk= (sqrt(3.0)/2) * l * (l**2 - m**2) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * (1 - l**2 + m**2) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2px")
                V_sk=-(sqrt(3.0)/2) * l * (l**2 - m**2) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                l * (1 - l**2 + m**2) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pydx2-y2")
                V_sk=  (sqrt(3.0)/2) * m * (l**2 - m**2) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * (1 + l**2 - m**2) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2py")
                V_sk= -(sqrt(3.0)/2) * m * (l**2 - m**2) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * (1 + l**2 - m**2) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzdx2-y2")
                V_sk= (sqrt(3.0)/2) * n * (l**2 - m**2) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                n * (l**2 - m**2) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2pz")
                V_sk= -(sqrt(3.0)/2) * n * (l**2 - m**2) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                n * (l**2 - m**2) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)

    case("pxdz2")
                V_sk= l * (n**2 - (l**2 + m**2)/2) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                sqrt(3.0) * l * (n**2) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dz2px")
                V_sk= -l * (n**2 - (l**2 + m**2)/2) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * l * (n**2) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pydz2")
                V_sk= m * (n**2 - (l**2 + m**2)/2) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                sqrt(3.0) * m * (n**2) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dz2py")
                V_sk= -m * (n**2 - (l**2 + m**2)/2) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * m * (n**2) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("pzdz2")
                V_sk= n * (n**2 - (l**2 + m**2)/2) * V_term(i,j,2,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * n * (l**2 + m**2) * V_term(i,j,2,3,3,dist_ij,V,Vn,r_d,nspecies)
    case("dz2pz")
                V_sk= -n * (n**2 - (l**2 + m**2)/2) * V_term(j,i,2,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                sqrt(3.0) * n * (l**2 + m**2) * V_term(j,i,2,3,3,dist_ij,V,Vn,r_d,nspecies)

       ! l = 2 with l = 2
    case("dxydxy")
                V_sk=  3 * (l**2) * (m**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                (l**2 + m**2 - 4*(l**2)*(m**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (n**2 + (l**2)*(m**2)) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dyzdyz")
                V_sk=  3 * (m**2) * (n**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                (m**2 + n**2 - 4*(m**2)*(n**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (l**2 + (m**2)*(n**2)) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dxzdxz")
                V_sk=  3 * (l**2) * (n**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                (l**2 + n**2 - 4*(l**2)*(n**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (m**2 + (l**2)*(n**2)) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)

    case("dxydyz")
                V_sk=   3 * l * (m**2) * n * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * n * (1 - 4*(m**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                l * n * (m**2 - 1) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dyzdxy")
                V_sk=  3 * l * (m**2) * n * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * n * (1 - 4*(m**2)) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                l * n * (m**2 - 1) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dxydxz")
                V_sk=  3 * m * (l**2) * n * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * n * (1 - 4*(l**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                m * n * (l**2 - 1) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dxzdxy")
                V_sk=  3 * m * (l**2) * n * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * n * (1 - 4*(l**2)) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                m * n * (l**2 - 1) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dyzdxz")
                V_sk= 3 * m * (n**2) * l * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * l * (1 - 4*(n**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                m * l * (n**2 - 1) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dxzdyz")
                V_sk= 3 * m * (n**2) * l * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                m * l * (1 - 4*(n**2)) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                m * l * (n**2 - 1) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)

    case("dxydx2-y2")
                V_sk= (3.0/2) * l * m * (l**2 - m**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                2 * l * m * (m**2 - l**2) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (1.0/2) * l * m * (l**2 - m**2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2dxy")
                V_sk=(3.0/2) * l * m * (l**2 - m**2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                2 * l * m * (m**2 - l**2) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (1.0/2) * l * m * (l**2 - m**2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dyzdx2-y2")
                V_sk= (3.0/2) * m * n * (l**2 - m**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * n * (1 + 2*(l**2 - m**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                m * n * (1 + (l**2 - m**2)/2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2dyz")
                V_sk=(3.0/2) * m * n * (l**2 - m**2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                m * n * (1 + 2*(l**2 - m**2)) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                m * n * (1 + (l**2 - m**2)/2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dxzdx2-y2")
                V_sk= (3.0/2) * l * n * (l**2 - m**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * n * (1 - 2*(l**2 - m**2)) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) - &
                l * n * (1 - (l**2 - m**2)/2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2dxz")
                V_sk= (3.0/2) * l * n * (l**2 - m**2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                l * n * (1 - 2*(l**2 - m**2)) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) - &
                l * n * (1 - (l**2 - m**2)/2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)

    case("dxydz2")
                V_sk= sqrt(3.0) * l * m * (n**2 - (l**2 + m**2)/2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                (2*sqrt(3.0)) * l * m * (n**2) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (sqrt(3.0)/2) * l * m * (1 + n**2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dz2dxy")
                V_sk= sqrt(3.0) * l * m * (n**2 - (l**2 + m**2)/2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) - &
                (2*sqrt(3.0)) * l * m * (n**2) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (sqrt(3.0)/2) * l * m * (1 + n**2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dyzdz2")
                V_sk=  sqrt(3.0) * m * n * (n**2 - (l**2 + m**2)/2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * m * n * (l**2 + m**2 - n**2) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) - &
                (sqrt(3.0)/2) * m * n * (l**2 + m**2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dz2dyz")
                V_sk=  sqrt(3.0) * m * n * (n**2 - (l**2 + m**2)/2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * m * n * (l**2 + m**2 - n**2) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) - &
                (sqrt(3.0)/2) * m * n * (l**2 + m**2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dxzdz2")
                V_sk= sqrt(3.0) * l * n * (n**2 - (l**2 + m**2)/2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * l * n * (l**2 + m**2 - n**2) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) - &
                (sqrt(3.0)/2) * l * n * (l**2 + m**2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dz2dxz")
                V_sk=  sqrt(3.0) * l * n * (n**2 - (l**2 + m**2)/2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * l * n * (l**2 + m**2 - n**2) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) - &
                (sqrt(3.0)/2) * l * n * (l**2 + m**2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)

    case("dx2-y2dx2-y2")
                V_sk=  (3.0/4) * ((l**2 - m**2)**2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                (l**2 + m**2 - (l**2 - m**2)**2) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (n**2 + ((l**2 - m**2)**2)/4) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2dz2")
                V_sk= (sqrt(3.0)/2) * (l**2 - m**2) * (n**2 - (l**2 + m**2)/2) * V_term(i,j,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * (n**2) * (m**2 - l**2) * V_term(i,j,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (sqrt(3.0)/4) * (1 + n**2) * (l**2 - m**2) * V_term(i,j,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dz2dx2-y2")
                V_sk=  (sqrt(3.0)/2) * (l**2 - m**2) * (n**2 - (l**2 + m**2)/2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                sqrt(3.0) * (n**2) * (m**2 - l**2) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (sqrt(3.0)/4) * (1 + n**2) * (l**2 - m**2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)
    case("dz2dz2")
                V_sk=  ((n**2 - (l**2 + m**2)/2)**2) * V_term(j,i,3,3,2,dist_ij,V,Vn,r_d,nspecies) + &
                3 * (n**2) * (l**2 + m**2) * V_term(j,i,3,3,3,dist_ij,V,Vn,r_d,nspecies) + &
                (3.0/4) * ((l**2 + m**2)**2) * V_term(j,i,3,3,4,dist_ij,V,Vn,r_d,nspecies)

       ! level s*
    case("s*s*")
                V_sk=    V_term(i,j,4,4,2,dist_ij,V,Vn,r_d,nspecies)
    case("ss*") 
                V_sk=   V_term(i,j,1,4,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*s") 
                V_sk=  V_term(j,i,1,4,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*px")
                V_sk=l * V_term(i,j,4,2,2,dist_ij,V,Vn,r_d,nspecies)
    case( "pxs*")
                V_sk=  -l * V_term(j,i,4,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*py")
                V_sk= m * V_term(i,j,4,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("pys*")
                V_sk=  -m * V_term(j,i,4,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*pz")
                V_sk=    n * V_term(i,j,4,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("pzs*")
                V_sk= -n * V_term(j,i,4,2,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*dxy")
                V_sk= sqrt(3.0) * (l*m) * V_term(i,j,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dxys*")
                V_sk= sqrt(3.0) * (l*m) * V_term(j,i,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*dyz")
                V_sk=  sqrt(3.0) * (m*n) * V_term(i,j,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dyzs*")
                V_sk=  sqrt(3.0) * (m*n) * V_term(j,i,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*dxz")
                V_sk= sqrt(3.0) * (l*n) * V_term(i,j,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dxzs*")
                V_sk=  sqrt(3.0) * (l*n) * V_term(j,i,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*dx2-y2")
                V_sk=  (sqrt(3.0)/2) * (l**2 - m**2) * V_term(i,j,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dx2-y2s*")
                V_sk=  (sqrt(3.0)/2) * (l**2 - m**2) * V_term(j,i,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("s*dz2")
                V_sk=   (n**2 - (l**2 + m**2)/2) * V_term(i,j,4,3,2,dist_ij,V,Vn,r_d,nspecies)
    case("dz2s*")
                V_sk= (n**2 - (l**2 + m**2)/2) * V_term(j,i,4,3,2,dist_ij,V,Vn,r_d,nspecies)

    case default
                V_sk = 0.0
    end select
    
    deallocate(int_type) 
    end function V_sk

    subroutine gen_ham_mpi(ie,ik_core,Ham,eigvec,eig,norb)
 
    use config
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ie,ik_core
    integer(kind=4),intent(in) :: norb
    complex(kind=8),intent(out) :: Ham(numprocs,norb,norb)
    complex(kind=8),intent(out) :: eigvec(numprocs,norb,norb)
    real(kind=8),intent(out) :: eig(numprocs,norb)

    mm = myid+1
    nn = k_start_core(mm)+ik_core-1
    call gen_ham(ie,k_grid(nn,:),&
         Ham(mm,:,:),eigvec(mm,:,:),eig(mm,:),norb)
     
    !call MPI_ALLREDUCE(Ham_red,Ham,&
    !         nk,MPI_COMPLEX,MPI_SUM,MPI_COMM_WORLD,mm)

    end subroutine gen_ham_mpi

    subroutine gen_ham(ie,kpoint,Ham,vec,eig,norb)
    
    use config
    use input
    use param
    use util
    
    implicit none
    
    integer(kind=4),intent(in) :: ie
    integer(kind=4) :: i,j, iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3 
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele, matrix_ele0
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=1
    real(kind=8),intent(in)  :: kpoint(3)
    real(kind=8)  :: kpoint_cart(3)
    integer(kind=4),intent(in) :: norb
    complex(kind=8) :: Ham(norb,norb)
    complex(kind=8) :: Ham1(norb,norb)
    complex(kind=8),intent(out) :: vec(norb,norb)
    real(kind=8),intent(out) :: eig(norb)
    Ham1 = 0.0d0
    eig = 0.0d0
    vec = 0.0d0
    Ham = 0.0d0 

    kpoint_cart = matmul(kpoint,recivec) ! new to convert

    i = 1
    do iatm = 1, natoms
        itype = int(positions(iatm,4))
        do iorb_idx = 1,atom_orb_num(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms
                jtype = int(positions(jatm,4))
                do jorb_idx =1,atom_orb_num(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    matrix_ele = 0.0d0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = positions(jatm,1:3)-positions(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec(1,1)+dble(i2)*latvec(2,1)+dble(i3)*latvec(3,1)
                               r_shift(2)=dble(i1)*latvec(1,2)+dble(i2)*latvec(2,2)+dble(i3)*latvec(3,2)
                               r_shift(3)=dble(i1)*latvec(1,3)+dble(i2)*latvec(2,3)+dble(i3)*latvec(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   matrix_ele = matrix_ele + &
                                   exp(i_imag*dot_product(r_ij,kpoint_cart))*V_ij
                               end if
                           end do
                        end do
                    end do 
                    Ham1(i,j) = matrix_ele
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do
    

    ! copy the Hamiltonian
    Ham = Ham1
    Ham = hermitian(Ham)
    end subroutine gen_ham

    subroutine gen_ham_Pmatrix(ik_core,Ham,eigvec,eig,norb,nl,nr)
 
    use config
    use input
    use param
    use util
  
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: mm,nn
    integer(kind=4),intent(in) :: ik_core
    integer(kind=4),intent(in) :: norb,nl,nr
    complex(kind=8) :: Ham(numprocs,norb,norb)
    complex(kind=8) :: Ham0(numprocs,norb,norb)
    complex(kind=8),intent(out) :: eigvec(numprocs,norb,norb)
    real(kind=8),intent(out) :: eig(numprocs,norb)
    complex(kind=8) :: Haml(nl,nl)
    complex(kind=8) :: Hamr(nr,nr)
    integer(kind=4) :: i,j, iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3 
    integer(kind=4) :: itxl, ityl, itxr, ityr
    integer(kind=4) :: jtxl, jtyl, jtxr, jtyr
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: mat_ele_core0(numprocs)
    complex(kind=8) :: mat_ele_core(numprocs)
    real(kind=8)    :: r_ij(3),r_shift(3)
    real(kind=8) :: kpoint_cart(3)
    integer(kind=4) :: l=1
    complex(kind=8) :: Px(norb,norb)
    complex(kind=8) :: Py(norb,norb)
    integer(kind=4) :: nuc

    nuc = n_bloch_x * n_bloch_y

    Ham0 = 0.0d0
    eigvec = 0.0d0
    eig = 0.0d0

    mm = myid+1
    nn = k_start_core(mm)+ik_core-1
    mat_ele_core = 0.0d0
    mat_ele_core0 = 0.0d0
  
    Ham = 0.0d0 
    Px = 0.0d0
    Py = 0.0d0

    kpoint_cart = matmul(k_grid(nn,:),recivec) ! new to convert
    ! crystal coordinates to cartesian coordinates

    ! for each atom, the number of orbitals is saved
 
    ! order of transverse uc positions: 
    ! (x1,y1), (x2,y1), (x3,y1)
    ! (x1,y2), (x2,y2), (x3,y2)
    ! (x1,y3), (x2,y3), (x3,y3)

    i = 1
    do iatm = 1, natoms
        itype = int(positions(iatm,4))
        ityl = (iatm-period_left-buffer_left)/(nuc*natoms_uc_l)
        itxl = modulo(iatm-period_left-buffer_left,&
                    nuc*natoms_uc_l)&
                    /(n_bloch_x*natoms_uc_l)
        ityr = (iatm-(natoms-period_right-buffer_right))&
                   /(nuc*natoms_uc_r)
        itxr = modulo(iatm-(natoms-period_right-buffer_right),&
                    nuc*natoms_uc_r)&
                    /(n_bloch_x*natoms_uc_r)
        do iorb_idx = 1,atom_orb_num(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms
                jtype = int(positions(jatm,4))
                jtyl = (jatm-period_left-buffer_left)/(nuc*natoms_uc_l)
                jtxl = modulo(jatm-period_left-buffer_left,&
                    nuc*natoms_uc_l)&
                    /(n_bloch_x*natoms_uc_l)
                jtyr = (jatm-(natoms-period_right-buffer_right))&
                   /(nuc*natoms_uc_r)
                jtxr = modulo(jatm-(natoms-period_right-buffer_right),&
                    nuc*natoms_uc_r)&
                    /(n_bloch_x*natoms_uc_r)
                do jorb_idx =1,atom_orb_num(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    mat_ele_core(mm) = 0.0d0
                    mat_ele_core0(mm) = 0.0d0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = positions(jatm,1:3)-positions(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec(1,1)+dble(i2)*latvec(2,1)+dble(i3)*latvec(3,1)
                               r_shift(2)=dble(i1)*latvec(1,2)+dble(i2)*latvec(2,2)+dble(i3)*latvec(3,2)
                               r_shift(3)=dble(i1)*latvec(1,3)+dble(i2)*latvec(2,3)+dble(i3)*latvec(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   mat_ele_core(mm) = mat_ele_core(mm) + &
                                   exp(i_imag*dot_product(r_shift,kpoint_cart))*V_ij
                                   mat_ele_core0(mm) = mat_ele_core0(mm) + V_ij
                               end if
                           end do
                        end do
                    end do 
                    Ham(mm,i,j) = mat_ele_core(mm)
                    if ((ityl.ge.0) .and. (ityl.le.nuc) &
                       .and. (jtyl.ge.0) .and. (jtyl.le.nuc)) then
                        Ham0(mm,i,j) = mat_ele_core(mm) *&
                        exp(-i_imag*dot_product(&
                        matmul((/float(jtxl-itxl),float(jtyl-ityl),0.0/),&
                        latvec_uc_l), kpoint_cart))
                    end if
                    if ((ityr.ge.0) .and. (ityr.le.nuc) &
                       .and. (jtyr.ge.0) .and. (jtyr.le.nuc)) then
                        Ham0(mm,i,j) = mat_ele_core(mm) *&
                        exp(-i_imag*dot_product(&
                        matmul((/float(jtxr-itxr),float(jtyr-ityr),0.0/),&
                        latvec_uc_r), kpoint_cart))
                    end if


                    !Px(i,j) = exp(i_imag*dot_product(&
                    !          matmul((/float(itx-jtx),0.0,0.0/),&
                    !          latvec_uc_l),kpoint_cart))
                    !Py(i,j) = exp(i_imag*dot_product(&
                    !          matmul((/0.0,float(ity-jty),0.0/),&
                    !          latvec_uc_l),kpoint_cart))
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do
    
    end subroutine gen_ham_Pmatrix


    subroutine find_k_pc_l(E,norb_pc_l,&
                          k_sc_cart,&
                          k_uc_cart,&
                          norb_uc_l,nl,u,pos_l,&
                          ksc,kpc_unfold,&
                          i_deg_pc_tot,vxyz_unfold)

    use config
    use input
    use util

    implicit none

    integer(kind=4),intent(in) :: norb_pc_l 
    real(kind=8),intent(in) :: E
    integer(kind=4) :: l,i
    real(kind=8),intent(in) :: k_sc_cart(3)
    real(kind=8),intent(in) :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_uc_l
    integer(kind=4),intent(in) :: nl
    real(kind=8),intent(in) :: pos_l(nl,3)
    complex(kind=8),intent(in) :: u(nl)
    complex(kind=8) :: u1(norb_uc_l,1)
    complex(kind=8) :: u2(norb_uc_l,1)
    complex(kind=8) :: u3(norb_uc_l,1)
    real(kind=8),intent(in) :: ksc(3)
    integer(kind=4) :: i1,i2,i3,ibz
    real(kind=8) :: kpc(3),weight
    integer(kind=4) :: i_deg_pc = 0
    integer(kind=4),intent(out) :: i_deg_pc_tot
    real(kind=8),intent(out) :: kpc_unfold(125,3) ! the first dimension is just a
                                    ! choice, as most materials 4-fold
                                    ! degeneracy hardly exist
    real(kind=8),intent(out) :: vxyz_unfold(125,3)
    real(kind=8) :: vxyz(3)
!    integer(kind=8),intent(in) :: nob
!    complex(kind=8) :: Ham_all(nob,nob)
    complex(kind=8) :: Ham_uc_l(norb_uc_l,norb_uc_l)
    real(kind=8) :: eigs(norb_uc_l)
    complex(kind=8) :: evs(norb_uc_l,norb_uc_l)
    real(kind=8) :: eigs_pc(norb_pc_l) 
    complex(kind=8) :: evs_pc(norb_pc_l,norb_pc_l) 
    real(kind=8)    :: v_pc(norb_pc_l,3)
    !call gen_ham_uc_l(norb_uc_l,&
    !                  k_uc_cart,Ham_uc_l,eigs,evs)
    !write(*,*) 'uc'
    !write(*,*) eigs
    !stop
    !u1(:,1) = u(:) ! from Ham for left lead
    !write(*,*) nl,norb_uc_l,size(evs,1)
    !u2 = 0.0d0
    !do i = 1,nl
    !   l = modulo(i-1,norb_uc_l)+1
    !    u2(l,1) = u2(l,1)+evs(i,7)*&
        !exp(i_imag*dot_product(positions_uc_l((i-1)/10+1,1:3),k_uc_cart))
    !    exp(i_imag*dot_product(pos_l(i,1:3),k_uc_cart)) 

    !end do
    !u1 = 0.0d0
    !do i = 1,nl
    !   l = modulo(i-1,norb_uc_l)+1
    !    u1(l,1) = u1(l,1) + u(i) *exp(i_imag*dot_product(pos_l(i,1:3),k_sc_cart))
    !end do
     ! u1(1:40,1) = u1(1:40,1)/sqrt(dot_product(u1(1:40,1),u1(1:40,1)))
   !u2(1:40,1) = u2(1:40,1)/sqrt(dot_product(u2(1:40,1),u2(1:40,1)))
!  write(*,*) dot_product(u1(1:40,1),u1(1:40,1))
!    u1 = u1/sqrt(dot_product(u1(:,1),u1(:,1)))
!    u2 = u2/sqrt(dot_product(u2(:,1),u2(:,1)))
!    write(*,*) abs(dot_product(u2(:,1),u1(:,1)))**2
!    stop

    !write(*,*) abs(dot_product(u1(:,1),u2(:,1)))**2
    !write(*,*) eigs(:)
   ! write(*,*) 'eg',maxval(abs(matmul(Ham_uc_l,u2(1:40,1))-eigs(7)*u2(1:40,1)))
    !    stop
!    do i =1 , 40
    ! write(*,*) (matmul(Ham_uc_l,u2(1:40,i))-eigs(7)*u2(1:40,1))
!    write(*,*) 'eg',maxval(abs(matmul(Ham_uc_l,evs(:,i))-eigs(i)*evs(:,i)))
!    end do
!    stop
 !   u3(:,1) = evs(:,7)
!    write(*,*)'u2', maxval(abs(matmul(Ham_all,u2(:,1))-E*u2(:,1)))
!    do i = 1,40
!       write(*,*) abs(u1(i,1))-abs(u2(i,1))!dot_product(u1(1:40,1),u2(1:40,1))
!   end do
    !log(u1(1:40,1)/u2(1:40,1))

    !write(*,*) dot_product(u1(1:40,1),u2(1:40,1))
    !do i = 1,40
    !  write(*,*) u1(i,1)/u1(80+i,1)-u2(i,1)/u2(80+i,1)
      !write(*,*) u2(1,1)/u2(41,1)
    !end do
    !u1(:,1)/u1(1,1)-u2(:,1)/u2(1,1) !/u1(42,1) 
!     write(*,*) u2(:,1) !/u2(42,1) 
!     stop

!    write(*,*) abs(dot_product(u1(:,1),u2(:,1)))
!    stop
    !do l = 1,size(eigs,1)
    ! if (abs(eigs(l)-E).lt.0.01) then
    !  write(*,*) l
  !end if
   ! end do
   ! stop 

    l = max(n_bloch_x,n_bloch_y)
    i_deg_pc_tot = 0
    kpc_unfold(:,:) = 0.0d0
    weight = 0.0
    vxyz = 0.0d0
    vxyz_unfold = 0.0d0
    
    ibz = in_ws(reci_uc_l,k_uc_cart)
!    write(*,*)'ibz111',ibz,matmul(k_uc_cart,inv_real(reci_uc_l))

    if (ibz.eq.2)then
        i_deg_pc = 1
        i_deg_pc_tot = i_deg_pc +i_deg_pc_tot
        kpc_unfold(i_deg_pc_tot,:) = k_uc_cart

        call group_velocity_pc_l(norb_pc_l,&
             k_uc_cart,eigs_pc,evs_pc,v_pc)
        do i1 = 1,norb_pc_l
            if (abs(E-eigs_pc(i1)).lt.1.0d-3) then
                vxyz = v_pc(i1,:)
                vxyz_unfold(1,:) = vxyz(:)
            end if
        end do
    else
        do i1 = -l,l
            do i2 = -l,l
                do i3 = -l,l
                    kpc(:) = k_uc_cart(:) &
                           + matmul((/dble(i1),dble(i2),dble(i3)/),reci_uc_l)
                    ibz = in_ws(recivec_pc_l,kpc)
                    if (ibz.ne.0) then
                        call compare_evs_l(E,norb_pc_l,nl,&
                              kpc,norb_uc_l,u,pos_l,&
                              k_sc_cart,k_uc_cart,ksc,i_deg_pc,&
                              weight,vxyz)
                        if (i_deg_pc .gt. 0) then
                            i_deg_pc_tot = i_deg_pc + i_deg_pc_tot
                            kpc_unfold(i_deg_pc_tot,:) = kpc
                            vxyz_unfold(i_deg_pc_tot,:) = vxyz
                        end if
                    end if
                end do
            end do
        end do
    end if
    
    if (i_deg_pc_tot .eq.0 )then
        write(*,*) 'one K point does not find the corresponding k in FBZ of left PC'
        write(*,*) vxyz,kpc
    end if

    end subroutine find_k_pc_l

    subroutine find_k_pc_r(E,norb_pc_r,&
                          k_sc_cart,&
                          k_uc_cart,&
                          norb_uc_r,nr,u,pos_r,&
                          ksc,kpc_unfold,&
                          i_deg_pc_tot,vxyz_unfold)

    use config
    use input
    use util

    implicit none

    integer(kind=4),intent(in) :: norb_pc_r
    real(kind=8),intent(in) :: E
    integer(kind=4) :: l
    real(kind=8),intent(in) :: k_sc_cart(3)
    real(kind=8),intent(in) :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_uc_r
    integer(kind=4),intent(in) :: nr
    real(kind=8),intent(in) :: pos_r(nr,3)
    complex(kind=8),intent(in) :: u(nr)
    real(kind=8),intent(in) :: ksc(3)
    integer(kind=4) :: i1,i2,i3,ibz
    real(kind=8) :: kpc(3),weight
    integer(kind=4) :: i_deg_pc = 0
    integer(kind=4),intent(out) :: i_deg_pc_tot
    real(kind=8),intent(out) :: kpc_unfold(125,3) ! the first dimension is just a
                                    ! choice, as most materials 4-fold
                                    ! degeneracy hardly exist
    real(kind=8),intent(out) :: vxyz_unfold(125,3)
    real(kind=8) :: vxyz(3)
    real(kind=8) :: eigs_pc(norb_pc_r) 
    complex(kind=8) :: evs_pc(norb_pc_r,norb_pc_r) 
    real(kind=8)    :: v_pc(norb_pc_r,3)

    l = max(n_bloch_x,n_bloch_y)
    i_deg_pc_tot = 0
    kpc_unfold(:,:) = 0.0d0
    weight = 0.0
    vxyz_unfold = 0.d0
    vxyz = 0.0d0

    ibz = in_ws(reci_uc_r,k_uc_cart)
    if (ibz.eq.2)then
        i_deg_pc = 1
        i_deg_pc_tot = i_deg_pc +i_deg_pc_tot
        kpc_unfold(i_deg_pc_tot,:) = k_uc_cart

        call group_velocity_pc_r(norb_pc_r,&
             k_uc_cart,eigs_pc,evs_pc,v_pc)    
        do i1 = 1,norb_pc_r
            if (abs(E-eigs_pc(i1)).lt.1.0d-3) then
                vxyz = v_pc(i1,:)
                vxyz_unfold(1,:) = vxyz(:)
            end if
        end do

    else
        do i1 = -l,l
            do i2 = -l,l
                do i3 = -l,l
                    kpc(:) = k_uc_cart(:) &
                           + matmul((/dble(i1),dble(i2),dble(i3)/),reci_uc_r)
                    ibz = in_ws(recivec_pc_r,kpc)

                    if (ibz.ne.0) then
                        call compare_evs_r(E,norb_pc_r,nr,&
                              kpc,norb_uc_r,u,pos_r,&
                              k_sc_cart,k_uc_cart,ksc,i_deg_pc,&
                              weight,vxyz)

                    !write(*,*) 'ibz',ibz,'i_deg_pc',i_deg_pc
                    !write(*,*) 'kpc', kpc
                    !write(*,*)'compare_evs_r_out',i1, vxyz
!                   ! write(*,*) 'ibz',i1,i2,i3,i_deg_pc
                        if (i_deg_pc .gt. 0) then
                            i_deg_pc_tot = i_deg_pc + i_deg_pc_tot
                            kpc_unfold(i_deg_pc_tot,:) = kpc
                            vxyz_unfold(i_deg_pc_tot,:) = vxyz
                        end if
                    end if
                end do
            end do
        end do
    end if
                        
    if (i_deg_pc_tot .eq.0 )then
        write(*,*) 'one K point does not find the corresponding k in FBZ of right PC'
        write(*,*) k_uc_cart,ksc
    end if

    end subroutine find_k_pc_r

    subroutine gen_ham_pc_l(norb_pc_l,&
                           kpoint_cart,Ham_pc_l)
    
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=2
    real(kind=8),intent(in)  :: kpoint_cart(3)
    integer(kind=4),intent(in) :: norb_pc_l
    complex(kind=8),intent(out) :: Ham_pc_l(norb_pc_l,norb_pc_l)
   ! primitive cell of left side
    Ham_pc_l = 0.0d0
 
    i = 1
    do iatm = 1, natoms_pc_l
        itype = int(positions_pc_l(iatm,4))
        do iorb_idx = 1,atom_orb_num_pc_l(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms_pc_l
                jtype = int(positions_pc_l(jatm,4))
                do jorb_idx =1,atom_orb_num_pc_l(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    matrix_ele = 0.0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = positions_pc_l(jatm,1:3)-positions_pc_l(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec_pc_l(1,1)+dble(i2)*latvec_pc_l(2,1)&
                                         +dble(i3)*latvec_pc_l(3,1)
                               r_shift(2)=dble(i1)*latvec_pc_l(1,2)+dble(i2)*latvec_pc_l(2,2)&
                                         +dble(i3)*latvec_pc_l(3,2)
                               r_shift(3)=dble(i1)*latvec_pc_l(1,3)+dble(i2)*latvec_pc_l(2,3)&
                                         +dble(i3)*latvec_pc_l(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   matrix_ele = matrix_ele + &
                                   exp(i_imag*dot_product(r_ij,kpoint_cart))*V_ij
                               end if
                           end do
                        end do
                    end do
                    
                    Ham_pc_l(i,j) = matrix_ele
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do
    Ham_pc_l = hermitian(Ham_pc_l)
    
    end subroutine gen_ham_pc_l

    ! Rectangular unit cell in lead in NEGF 
    subroutine gen_ham_uc_l(norb_uc_l,&
                           k_uc_cart,Ham_uc_l,&
                           eigs,evs)
    
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=1
    real(kind=8),intent(in)  :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_uc_l
    complex(kind=8),intent(out) :: Ham_uc_l(norb_uc_l,norb_uc_l)
    real(kind=8),intent(out) :: eigs(norb_uc_l)
    complex(kind=8),intent(out) :: evs(norb_uc_l,norb_uc_l)

   ! unit cell of left side
    Ham_uc_l = 0.0d0
 
    i = 1
    do iatm = 1, natoms_uc_l
        itype = int(positions_uc_l(iatm,4))
        do iorb_idx = 1,atom_orb_num_uc_l(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms_uc_l
                jtype = int(positions_uc_l(jatm,4))
                do jorb_idx =1,atom_orb_num_uc_l(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    matrix_ele = 0.0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = positions_uc_l(jatm,1:3)-positions_uc_l(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec_uc_l(1,1)+dble(i2)*latvec_uc_l(2,1)&
                                         +dble(i3)*latvec_uc_l(3,1)
                               r_shift(2)=dble(i1)*latvec_uc_l(1,2)+dble(i2)*latvec_uc_l(2,2)&
                                         +dble(i3)*latvec_uc_l(3,2)
                               r_shift(3)=dble(i1)*latvec_uc_l(1,3)+dble(i2)*latvec_uc_l(2,3)&
                                         +dble(i3)*latvec_uc_l(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   matrix_ele = matrix_ele + &
                                   exp(i_imag*dot_product(r_ij,k_uc_cart))*V_ij
                               end if
                           end do
                        end do
                    end do
                    
                    Ham_uc_l(i,j) = matrix_ele
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do

    Ham_uc_l = hermitian(Ham_uc_l)

    call eigenH(Ham_uc_l,eigs,evs)

    end subroutine gen_ham_uc_l

    subroutine gen_ham_uc_r(norb_uc_r,&
                           k_uc_cart,Ham_uc_r,&
                           eigs,evs)
    
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=1
    real(kind=8),intent(in)  :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_uc_r
    complex(kind=8),intent(out) :: Ham_uc_r(norb_uc_r,norb_uc_r)
    real(kind=8),intent(out) :: eigs(norb_uc_r)
    complex(kind=8),intent(out) :: evs(norb_uc_r,norb_uc_r)

   ! unit cell of left side
    Ham_uc_r = 0.0d0
 
    i = 1
    do iatm = 1, natoms_uc_r
        itype = int(positions_uc_r(iatm,4))
        do iorb_idx = 1,atom_orb_num_uc_r(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms_uc_r
                jtype = int(positions_uc_r(jatm,4))
                do jorb_idx =1,atom_orb_num_uc_r(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    matrix_ele = 0.0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = positions_uc_r(jatm,1:3)-positions_uc_r(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec_uc_r(1,1)+dble(i2)*latvec_uc_r(2,1)&
                                         +dble(i3)*latvec_uc_r(3,1)
                               r_shift(2)=dble(i1)*latvec_uc_r(1,2)+dble(i2)*latvec_uc_r(2,2)&
                                         +dble(i3)*latvec_uc_r(3,2)
                               r_shift(3)=dble(i1)*latvec_uc_r(1,3)+dble(i2)*latvec_uc_r(2,3)&
                                         +dble(i3)*latvec_uc_r(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   matrix_ele = matrix_ele + &
                                   exp(i_imag*dot_product(r_ij,k_uc_cart))*V_ij
                               end if
                           end do
                        end do
                    end do
                    
                    Ham_uc_r(i,j) = matrix_ele
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do

    Ham_uc_r = hermitian(Ham_uc_r)

    call eigenH(Ham_uc_r,eigs,evs)

    end subroutine gen_ham_uc_r

    subroutine gen_ham_slab(norb_slab,&
                           k_uc_cart,Ham_slab,&
                           eigs,evs,pos_slab,&
                           latvec_slab,atom_orb_num_slab)
    
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=1
    real(kind=8),intent(in)  :: k_uc_cart(3)
    integer(kind=4),intent(in) :: norb_slab
    complex(kind=8),intent(out) :: Ham_slab(norb_slab,norb_slab)
    real(kind=8),intent(out) :: eigs(norb_slab)
    complex(kind=8),intent(out) :: evs(norb_slab,norb_slab)
    real(kind=8),intent(in) :: pos_slab(:,:)
    real(kind=8),intent(in) :: latvec_slab(3,3)
    integer(kind=4),intent(in) :: atom_orb_num_slab(:)
    integer(kind=4) :: natoms_slab

    natoms_slab = size(pos_slab,1)
    Ham_slab = 0.0d0
 
    i = 1
    do iatm = 1, natoms_slab
        itype = int(pos_slab(iatm,4))
        do iorb_idx = 1,atom_orb_num_slab(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms_slab
                jtype = int(pos_slab(jatm,4))
                do jorb_idx =1,atom_orb_num_slab(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    matrix_ele = 0.0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = pos_slab(jatm,1:3)-pos_slab(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec_slab(1,1)+dble(i2)*latvec_slab(2,1)&
                                         +dble(i3)*latvec_slab(3,1)
                               r_shift(2)=dble(i1)*latvec_slab(1,2)+dble(i2)*latvec_slab(2,2)&
                                         +dble(i3)*latvec_slab(3,2)
                               r_shift(3)=dble(i1)*latvec_slab(1,3)+dble(i2)*latvec_slab(2,3)&
                                         +dble(i3)*latvec_slab(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   matrix_ele = matrix_ele + &
                                   exp(i_imag*dot_product(r_ij,k_uc_cart))*V_ij
                               end if
                           end do
                        end do
                    end do
                    
                    Ham_slab(i,j) = matrix_ele
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do

    Ham_slab = hermitian(Ham_slab)

    call eigenH(Ham_slab,eigs,evs)

    end subroutine gen_ham_slab


    subroutine gen_hf_pc_l(norb_pc_l,&
                           kpoint_cart,HF_pc_l)
    
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx,&
                       ncart
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=2
    real(kind=8),intent(in)  :: kpoint_cart(3)
    integer(kind=4),intent(in) :: norb_pc_l
    complex(kind=8),intent(out) :: HF_pc_l(3,norb_pc_l,norb_pc_l)
   ! primitive cell of left side
    HF_pc_l = 0.0d0
 
    do ncart = 1,3
        i = 1
        do iatm = 1, natoms_pc_l
            itype = int(positions_pc_l(iatm,4))
            do iorb_idx = 1,atom_orb_num_pc_l(iatm)
                iorbit = atom_orb_list(itype,iorb_idx)
                j = 1
                do jatm = 1, natoms_pc_l
                    jtype = int(positions_pc_l(jatm,4))
                    do jorb_idx =1,atom_orb_num_pc_l(jatm)
                        jorbit = atom_orb_list(jtype,jorb_idx)
                        matrix_ele = 0.0
                        do i1 = -l,l
                            do i2 = -l,l
                               do i3 = -l,l
                                   r_ij(:) = positions_pc_l(jatm,1:3)-positions_pc_l(iatm,1:3)
                                   r_shift(1)=dble(i1)*latvec_pc_l(1,1)+dble(i2)*latvec_pc_l(2,1)&
                                             +dble(i3)*latvec_pc_l(3,1)
                                   r_shift(2)=dble(i1)*latvec_pc_l(1,2)+dble(i2)*latvec_pc_l(2,2)&
                                             +dble(i3)*latvec_pc_l(3,2)
                                   r_shift(3)=dble(i1)*latvec_pc_l(1,3)+dble(i2)*latvec_pc_l(2,3)&
                                             +dble(i3)*latvec_pc_l(3,3)
                                   r_ij(:) = r_ij(:) + r_shift(:)
                                   if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                       if (iatm.eq.jatm) then
                                           V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                       else    
                                           V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                       end if
                                       matrix_ele = matrix_ele + &
                                       exp(i_imag*dot_product(r_ij,kpoint_cart))*V_ij&
                                       *r_ij(ncart)
                                   end if
                               end do
                            end do
                        end do
                        
                        HF_pc_l(ncart,i,j) = matrix_ele
                        j = j+1
                    end do
                end do
                i = i +1
            end do
        end do
        HF_pc_l(ncart,:,:) = i_imag*antihermitian(HF_pc_l(ncart,:,:))
    end do
        
    end subroutine gen_HF_pc_l


    subroutine group_velocity_pc_l(norb_pc,&
               kpoint_cart,eigs,evs,vels)

    use util

    implicit none

    integer(kind=4),intent(in) :: norb_pc
    real(kind=8),intent(in)  :: kpoint_cart(3)
    real(kind=8),intent(out) :: eigs(norb_pc)
    complex(kind=8),intent(out) :: evs(norb_pc,norb_pc)
    real(kind=8),intent(out) :: vels(norb_pc,3)
    complex(kind=8)  :: Ham_pc(norb_pc,norb_pc)
    complex(kind=8)  :: Ham_pc_p(norb_pc,norb_pc)
    complex(kind=8)  :: Ham_pc_m(norb_pc,norb_pc)
    complex(kind=8) :: HF_pc(3,norb_pc,norb_pc)
    complex(kind=8) :: dHdk(norb_pc,norb_pc)
    complex(kind=8) :: Vmat(norb_pc,norb_pc)
    real(kind=8) :: devi = 1.0d-5
    real(kind=8) :: ktemp1(3),ktemp2(3),delta_k
    integer(kind=4) :: i,j
    complex(kind=8) :: temp(norb_pc,norb_pc)


    call gen_ham_pc_l(norb_pc,kpoint_cart,Ham_pc)
    call eigenH(Ham_pc,eigs,evs)

    ! finite difference method

    !do i = 1,3
    !    delta_k = devi
        
    !    ktemp1 = kpoint_cart
    !    ktemp1(i) = ktemp1(i)-delta_k 
    !    ktemp2 = kpoint_cart
    !    ktemp2(i) = ktemp2(i)+delta_k

    !    call gen_ham_pc_l(norb_pc,&
    !                 ktemp1,Ham_pc_m) 
    !    call gen_ham_pc_l(norb_pc,&
    !                 ktemp2,Ham_pc_p) 
        
    !    dHdk = (Ham_pc_p-Ham_pc_m)/(2.0*delta_k)
    !    Vmat = matmul(matmul(&
    !           transpose(dconjg(evs)),dHdk),evs)
    !    do j = 1,norb_pc
    !        vels(j,i) = real(Vmat(j,j))  ! v_x/y/z
    !        write(*,*) vels(j,i)
    !    end do 
    !end do 

    
    call  gen_hf_pc_l(norb_pc,&
                kpoint_cart,HF_pc)

    do i = 1,3
        Vmat = matmul(matmul(transpose(dconjg(evs)),HF_pc(i,:,:)),evs)
        do j = 1,norb_pc
            vels(j,i) = real(Vmat(j,j))  ! v_x/y/z
        end do
    end do

    end subroutine group_velocity_pc_l

    subroutine gen_ham_pc_r(norb_pc_r,&
                          kpoint_cart,Ham_pc_r)
 
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=2
    real(kind=8),intent(in)  :: kpoint_cart(3)
    integer(kind=4),intent(in) :: norb_pc_r
    complex(kind=8),intent(out) :: Ham_pc_r(norb_pc_r,norb_pc_r)

    ! primitive cell of left side
    Ham_pc_r = 0.0d0
 
    i = 1
    do iatm = 1, natoms_pc_r
        itype = int(positions_pc_r(iatm,4))
        do iorb_idx = 1,atom_orb_num_pc_r(iatm)
            iorbit = atom_orb_list(itype,iorb_idx)
            j = 1
            do jatm = 1, natoms_pc_r
                jtype = int(positions_pc_r(jatm,4))
                do jorb_idx =1,atom_orb_num_pc_r(jatm)
                    jorbit = atom_orb_list(jtype,jorb_idx)
                    matrix_ele = 0.0
                    do i1 = -l,l
                        do i2 = -l,l
                           do i3 = -l,l
                               r_ij(:) = positions_pc_r(jatm,1:3)-positions_pc_r(iatm,1:3)
                               r_shift(1)=dble(i1)*latvec_pc_r(1,1)+dble(i2)*latvec_pc_r(2,1)&
                                         +dble(i3)*latvec_pc_r(3,1)
                               r_shift(2)=dble(i1)*latvec_pc_r(1,2)+dble(i2)*latvec_pc_r(2,2)&
                                         +dble(i3)*latvec_pc_r(3,2)
                               r_shift(3)=dble(i1)*latvec_pc_r(1,3)+dble(i2)*latvec_pc_r(2,3)&
                                         +dble(i3)*latvec_pc_r(3,3)
                               r_ij(:) = r_ij(:) + r_shift(:)
                               if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                   if (iatm.eq.jatm) then
                                       V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                   else    
                                       V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                   end if
                                   matrix_ele = matrix_ele + &
                                   exp(i_imag*dot_product(r_ij,kpoint_cart))*V_ij
                               end if
                           end do
                        end do
                    end do
                    
                    Ham_pc_r(i,j) = matrix_ele
                    j = j+1
                end do
            end do
            i = i +1
        end do
    end do    

    Ham_pc_r = hermitian(Ham_pc_r)

    end subroutine gen_ham_pc_r

    subroutine gen_hf_pc_r(norb_pc_r,&
                          kpoint_cart,HF_pc_r)

    ! hellman-feynman theorem
 
    use config
    use input
    use param
    use util
   
    implicit none
    
    include "mpif.h"

    integer(kind=4) :: i, j, k, id,&
                       iatm, jatm, iorb_idx, jorb_idx,&
                       ncart
    integer(kind=4) :: itype, jtype, i1, i2, i3
    character(len=8) :: iorbit, jorbit
    real(kind=8) :: V_ij, V_on
    real(kind=8) :: r(3)
    complex(kind=8) :: matrix_ele
    real(kind=8)    :: r_ij(3),r_shift(3)
    integer(kind=4) :: l=2
    real(kind=8),intent(in)  :: kpoint_cart(3)
    integer(kind=4),intent(in) :: norb_pc_r
    complex(kind=8),intent(out) :: HF_pc_r(3,norb_pc_r,norb_pc_r)

    ! primitive cell of left side
    HF_pc_r = 0.0d0
    do ncart = 1,3 
        i = 1
        do iatm = 1, natoms_pc_r
            itype = int(positions_pc_r(iatm,4))
            do iorb_idx = 1,atom_orb_num_pc_r(iatm)
                iorbit = atom_orb_list(itype,iorb_idx)
                j = 1
                do jatm = 1, natoms_pc_r
                    jtype = int(positions_pc_r(jatm,4))
                    do jorb_idx =1,atom_orb_num_pc_r(jatm)
                        jorbit = atom_orb_list(jtype,jorb_idx)
                        matrix_ele = 0.0
                        do i1 = -l,l
                            do i2 = -l,l
                               do i3 = -l,l
                                   r_ij(:) = positions_pc_r(jatm,1:3)-positions_pc_r(iatm,1:3)
                                   r_shift(1)=dble(i1)*latvec_pc_r(1,1)+dble(i2)*latvec_pc_r(2,1)&
                                             +dble(i3)*latvec_pc_r(3,1)
                                   r_shift(2)=dble(i1)*latvec_pc_r(1,2)+dble(i2)*latvec_pc_r(2,2)&
                                             +dble(i3)*latvec_pc_r(3,2)
                                   r_shift(3)=dble(i1)*latvec_pc_r(1,3)+dble(i2)*latvec_pc_r(2,3)&
                                             +dble(i3)*latvec_pc_r(3,3)
                                   r_ij(:) = r_ij(:) + r_shift(:)
                                   if (sqrt(r_ij(1)**2+r_ij(2)**2+r_ij(3)**2)<r_cutoff) then
                                       if (iatm.eq.jatm) then
                                           V_ij = onsite(itype,jtype,iorbit,jorbit,V,nspecies)
                                       else    
                                           V_ij = V_sk(itype,jtype,iorbit,jorbit,r_ij,V,Vn,r_d,nspecies)
                                       end if
                                       matrix_ele = matrix_ele + &
                                       exp(i_imag*dot_product(r_ij,kpoint_cart))*V_ij&
                                       *r_ij(ncart)
                                   end if
                               end do
                            end do
                        end do
                        
                        HF_pc_r(ncart,i,j) = matrix_ele
                        j = j+1
                    end do
                end do
                i = i +1
            end do
        end do
        HF_pc_r(ncart,:,:) = i_imag*antihermitian(HF_pc_r(ncart,:,:)) 
    end do

    end subroutine gen_hf_pc_r

    subroutine group_velocity_pc_r(norb_pc,&
               kpoint_cart,eigs,evs,vels)

    use util

    implicit none

    integer(kind=4),intent(in) :: norb_pc
    real(kind=8),intent(in)  :: kpoint_cart(3)
    real(kind=8),intent(out) :: eigs(norb_pc)
    complex(kind=8),intent(out) :: evs(norb_pc,norb_pc)
    real(kind=8),intent(out) :: vels(norb_pc,3)
    complex(kind=8)  :: Ham_pc(norb_pc,norb_pc)
    complex(kind=8)  :: Ham_pc_p(norb_pc,norb_pc)
    complex(kind=8)  :: Ham_pc_m(norb_pc,norb_pc)
    complex(kind=8) :: HF_pc(3,norb_pc,norb_pc)
    complex(kind=8) :: dHdk(norb_pc,norb_pc)
    complex(kind=8) :: Vmat(norb_pc,norb_pc)
    real(kind=8) :: devi = 1.0d-5
    real(kind=8) :: ktemp1(3),ktemp2(3),delta_k
    integer(kind=4) :: i,j


    call gen_ham_pc_r(norb_pc, kpoint_cart,Ham_pc)
    call eigenH(Ham_pc,eigs,evs)

    ! finite difference method to calculate group velocity

    !do i = 1,3
    !    delta_k = devi
        
    !    ktemp1 = kpoint_cart
    !    ktemp1(i) = ktemp1(i)-delta_k 
    !    ktemp2 = kpoint_cart
    !    ktemp2(i) = ktemp2(i)+delta_k

    !    call gen_ham_pc_r(norb_pc,&
    !                 ktemp1,Ham_pc_m) 
    !    call gen_ham_pc_r(norb_pc,&
    !                 ktemp2,Ham_pc_p) 
        
    !    dHdk = (Ham_pc_p-Ham_pc_m)/(2.0*delta_k)
    !    Vmat = matmul(matmul(&
    !           transpose(dconjg(evs)),dHdk),evs)
    !    do j = 1,norb_pc
    !        vels(j,i) = real(Vmat(j,j))  ! v_x/y/z
    !        write(*,*) vels(j,i)
    !    end do 
    !end do 
    call  gen_hf_pc_r(norb_pc,&
                kpoint_cart,HF_pc)
    do i = 1,3
        Vmat = matmul(matmul(transpose(dconjg(evs)),HF_pc(i,:,:)),evs)
        do j = 1,norb_pc
            vels(j,i) = real(Vmat(j,j))  ! v_x/y/z
        end do
    end do


    end subroutine group_velocity_pc_r


    subroutine compare_evs_l(E,norb_pc_l,nl,&
                          kpoint_cart,norb_uc_l,u,posl,&
                          k_sc_cart,k_uc_cart,ksc,i_deg_pc,&
                          weight,vxyz)
    use config
    use input
    use param
    use util

    implicit none

    real(kind=8),intent(in)  :: E
    real(kind=8),intent(in)  :: k_sc_cart(3)
    real(kind=8),intent(in)  :: k_uc_cart(3)
    real(kind=8),intent(in)  :: kpoint_cart(3)
    real(kind=8),intent(in)  :: ksc(3)
    integer(kind=4),intent(in) :: nl
    real(kind=8),intent(in) :: posl(nl,4)
    real(kind=8) :: pos_l_shiftz(nl,4)
    integer(kind=4),intent(in) :: norb_pc_l,norb_uc_l
    complex(kind=8),intent(in):: u(nl)
    real(kind=8),intent(out) :: weight
    real(kind=8),intent(out) :: vxyz(3)
    real(kind=8) :: eigs(norb_pc_l)
    complex(kind=8) :: evs(norb_pc_l,norb_pc_l)
    complex(kind=8) :: u_pc(norb_pc_l)
    complex(kind=8) :: u_sc(nl)
    real(kind=8) :: vels(norb_pc_l,3)
    complex(kind=8),allocatable :: evs_uc(:)
    complex(kind=8),allocatable :: evs_l(:)
    complex(kind=8) :: temp_sum
    integer(kind=4),intent(out) :: i_deg_pc
    integer(kind=4) :: npc ! number of 
                            ! primitive cells in unitcell
    integer(kind=4) :: i,j,k,id,idx,l,j1


    call group_velocity_pc_l(norb_pc_l,&
               kpoint_cart,eigs,evs,vels)

    npc = natoms_uc_l/natoms_pc_l
        
    allocate(evs_uc(norb_pc_l*npc)) 
    allocate(evs_l(norb_uc_l*n_bloch_x*n_bloch_y))  
    evs_l = 0.0d0
    evs_uc = 0.0d0
    u_pc = 0.0d0
    u_sc = 0.0d0
    
    pos_l_shiftz  = posl
    pos_l_shiftz(:,3) = pos_l_shiftz(:,3) - minval(posl(:,3))
   
    i_deg_pc = 0
    weight = 0.0d0
    vxyz = 0.0d0
!    write(*,*) 'lst' 
    do i =1,norb_pc_l
!        write(*,*) E-eigs(i)
        if (abs(E-eigs(i)).lt.1.0d-4)then
            vxyz = vels(i,:)

            i_deg_pc = 1
            do j = 1,n_bloch_y
                do j1 = 1,n_bloch_x
                    do l = 1,npc 
                        do k = 1,norb_pc_l
                            id = (((j-1)*n_bloch_x+j1-1)*npc+l-1)*norb_pc_l+k
                            idx = int(pos_l_shiftz(id,4)) 
                            evs_l(id) =&
                            evs(k,i)*exp(i_imag*(dot_product(pos_l_shiftz(id,1:3), &
                            kpoint_cart))) 

                            u_sc(id) = u(id) *& 
                            exp(-i_imag*dot_product(pos_l_shiftz(id,1:3),(/0.0d0,0.0d0,k_sc_cart(3)/))) * & 
                            exp(i_imag*(dot_product(pos_l_shiftz(id,1:3),k_sc_cart)))

                        end do
                    end do
                 end do
            end do
            evs_l = evs_l / sqrt(dble(n_bloch_x*n_bloch_y*npc))
            temp_sum = 0.0
            temp_sum = dot_product(u_sc(:),evs_l(:))
            weight = abs(temp_sum)**2
!            if (abs(temp_sum) < 0.50) then          
!                i_deg_pc = i_deg_pc -1
!            end if
        end if
    end do
   
    end subroutine compare_evs_l

    subroutine compare_evs_r(E,norb_pc_r,nr,&
                          kpoint_cart,norb_uc_r,u,posr,&
                          k_sc_cart,&
                          k_uc_cart,ksc,i_deg_pc,&
                          weight,vxyz)
    use config
    use input
    use param
    use util

    implicit none

    real(kind=8),intent(in) :: E
    real(kind=8),intent(in) :: k_sc_cart(3)
    real(kind=8),intent(in)  :: k_uc_cart(3)
    real(kind=8),intent(in)  :: kpoint_cart(3)
    real(kind=8),intent(in)  :: ksc(3)
    integer(kind=4),intent(in) :: nr
    real(kind=8),intent(in) :: posr(nr,4)
    real(kind=8) :: pos_r_shiftz(nr,4)
    integer(kind=4),intent(in) :: norb_pc_r,norb_uc_r
    complex(kind=8),intent(in):: u(nr)
    real(kind=8),intent(out) :: weight
    real(kind=8),intent(out) :: vxyz(3)
    real(kind=8) :: eigs(norb_pc_r)
    complex(kind=8) :: evs(norb_pc_r,norb_pc_r)
    real(kind=8) :: vels(norb_pc_r,3)
    complex(kind=8),allocatable :: evs_uc(:)
    complex(kind=8),allocatable :: evs_r(:)
    complex(kind=8)  :: u_sc(nr)
    complex(kind=8) :: temp_sum
    integer(kind=4) :: i_deg_pc
    integer(kind=4) :: npc ! number of 
                            ! primitive cells in unitcell
    integer(kind=4) :: i,j,k,id,idx,l,j1

    call group_velocity_pc_r(norb_pc_r,&
            kpoint_cart,eigs,evs,vels)

    npc = natoms_uc_r/natoms_pc_r
        
    allocate(evs_uc(norb_pc_r*npc)) 
    allocate(evs_r(norb_uc_r*n_bloch_x*n_bloch_y))  
    
    pos_r_shiftz  = posr
    pos_r_shiftz(:,3) = pos_r_shiftz(:,3) - minval(posr(:,3))
  
    i_deg_pc = 0
    weight = 0.0d0
    u_sc = 0.0d0
    vxyz = 0.0d0

    do i =1,norb_pc_r
        if (abs(E-eigs(i)).lt.1.0d-4)then
            vxyz = vels(i,:)
            i_deg_pc = 1

            do j = 1,n_bloch_y
                do j1 = 1,n_bloch_y
                    do l = 1,npc 
                        do k = 1,norb_pc_r
                            id = (((j-1)*n_bloch_x+j1-1)*npc+l-1)*norb_pc_r+k
                            idx = int(pos_r_shiftz(id,4)) 
                            evs_r(id) =&
                            evs(k,i)*exp(i_imag*(dot_product(pos_r_shiftz(id,1:3), &
                            kpoint_cart))) 

                            u_sc(id) = u(id) *& 
                            exp(-i_imag*dot_product(pos_r_shiftz(id,1:3),(/0.0d0,0.0d0,k_sc_cart(3)/))) * & 
                            exp(i_imag*(dot_product(pos_r_shiftz(id,1:3),k_sc_cart)))

                        end do
                    end do
                end do
            end do
            evs_r = evs_r / sqrt(dble(n_bloch_x*n_bloch_y*npc))
            temp_sum = 0.0
            temp_sum = dot_product(u_sc(:),evs_r(:))
            weight = abs(temp_sum)**2

!            if (abs(temp_sum) < 0.5) then          
!                i_deg_pc = i_deg_pc -1
!            end if
         end if
    end do

    end subroutine compare_evs_r

end module hamiltonian
