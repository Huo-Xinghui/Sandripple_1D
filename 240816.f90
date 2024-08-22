module public_val
    ! constants
    implicit none
    integer, parameter :: dbPc = 8 !selected_real_kind(15, 307)
    real(kind=dbPc), parameter :: pi = 3.14159265358979323846
    ! computational domain
    real(kind=dbPc), parameter :: xMax = 0.5 ! x size
    real(kind=dbPc), parameter :: yMax = 0.01 ! y sizw
    real(kind=dbPc), parameter :: zMax = 0.3 ! z size
    real(kind=dbPc), parameter :: area = xMax*yMax
    integer, parameter :: mx = 502 ! x grid num +2
    integer, parameter :: my = 12 ! y grid num +2
    integer, parameter :: mz = 150 ! z grid num
    integer, parameter :: mzUni = 30 ! z grid number above which zDiff becomes uniform
    real(kind=dbPc), parameter :: xDiff = xMax/(mx - 2)
    real(kind=dbPc), parameter :: yDiff = yMax/(my - 2)
    integer, parameter :: nNodes = 5 ! num of subdomain
    integer, parameter :: mxNode = (mx - 2)/nNodes + 2 ! x direction grid num for every subdomain
    ! particle
    logical, parameter :: midairCollision = .false.
    ! particle diameter: whichDiameterDist = 0: normal distribution, 1: uniform diameter, 2: Bernoulli distribution
    ! whichDiameterDist=0: npdf must >= 3, mu=dpa, sigma=dpStddDev, range:mu-3*sigma ~ mu+3*sigma
    ! whichDiameterDist=1: npdf must = 1, d=dpa
    ! whichDiameterDist=2: npdf must = 2, p1=prob1, p2=1-prob1, d1=dpa-dpStddDev, d2=dpa+dpStddDev
    integer, parameter :: whichDiameterDist = 1
    integer, parameter :: npdf = 1 ! bin num of particle distribution
    integer, parameter :: pNumInit = 10 ! initial particle num
    integer, parameter :: maxEjectNum = 10000 ! max eject particle num in one time step
    integer, parameter :: maxNum = 100000 ! max particle num in one subdomain
    integer, parameter :: pNumInGridMax = maxNum/mxNode !/(my)
    integer, parameter :: pNumExchMax = maxNum/10
    real(kind=dbPc), parameter :: dpa = 2.5e-4 ! average particle diameter
    real(kind=dbPc), parameter :: dpStddDev = 1.0e-4 ! particle diameter standard deviation
    real(kind=dbPc), parameter :: prob1 = 0.5 ! probability one of Bernoulli distribution
    real(kind=dbPc), parameter :: binWidth = 4.0*dpStddDev/npdf
    real(kind=dbPc), parameter :: binStart = dpa - 2.0*dpStddDev
    real(kind=dbPc), parameter :: binEnd = dpa + 2.0*dpStddDev
    real(kind=dbPc), parameter :: resN = 0.9 ! normal restitution coefficient
    real(kind=dbPc), parameter :: resT = 0.0 ! tangential restitution coefficient
    real(kind=dbPc), parameter :: rhoP = 2650.0 ! particle density
    real(kind=dbPc), parameter :: por = 0.6 ! bedform porosity
    ! bed surface
    logical, parameter :: predefineSurface = .false.
    real(kind=dbPc), parameter :: initSurfElevation = 5.0e-2 ! initial bed height
    real(kind=dbPc), parameter :: amp = 0.005 ! amplitude of the predefined surface
    real(kind=dbPc), parameter :: omg = 32.0*pi ! wave number (friquence) of the predefined surface
    real(kind=dbPc), parameter :: z0 = dpa/30.0 ! roughness height
    real(kind=dbPc), parameter :: repostAngle = 35.0 ! repost angle
    real(kind=dbPc), parameter :: tanRepostAngle = tan(repostAngle/180.0*pi)
    real(kind=dbPc), parameter :: blockHeight = dpa*5.0
    real(kind=dbPc), parameter :: initBlock = xDiff*yDiff*blockHeight*por
    ! fluid
    real(kind=dbPc), parameter :: uStar = 0.45 ! fractional velocity
    real(kind=dbPc), parameter :: rho = 1.263 ! fluid density
    real(kind=dbPc), parameter :: nu = 1.51e-5 ! kinetic viscosity
    real(kind=dbPc), parameter :: kapa = 0.4 ! von Kaman's constant
    real(kind=dbPc), parameter :: zDiffMin = nu/uStar ! smallest z grid size
    ! iteration
    integer, parameter :: parStart = 1 ! Iteration num when the particle calculation starts
    integer, parameter :: oneSecond = 20000
    integer, parameter :: intervalField = oneSecond*120
    integer, parameter :: intervalProfile = oneSecond*60
    integer, parameter :: intervalMonitor = oneSecond
    integer, parameter :: intervalParticle = oneSecond*60
    integer, parameter :: intervalSurface = oneSecond*30
    integer, parameter :: intervalStatistics = oneSecond*60
    integer, parameter :: intervalCreateFile = oneSecond*240
    real(kind=dbPc), parameter :: dt = 5.0e-5 ! time step
    real(kind=dbPc), parameter :: endTime = 3600.0 ! The time that the simulation lasts

    ! variables
    ! MPI
    integer :: realtype
    integer :: inttype
    integer :: procs
    integer :: comm
    integer :: myID
    integer :: sliceType
    integer :: edgeType1, edgeType2
    integer, dimension(2) :: neighbor
    ! Fluid field
    real(kind=dbPc), dimension(mxNode) :: x
    real(kind=dbPc), dimension(my) :: y
    real(kind=dbPc), dimension(mz) :: zDiff
    real(kind=dbPc), dimension(mz) :: z
    real(kind=dbPc), dimension(mxNode, my, mz) :: zReal
    real(kind=dbPc), dimension(mxNode) :: xu
    real(kind=dbPc), dimension(my) :: yv
    real(kind=dbPc), dimension(mz + 1) :: zw
    real(kind=dbPc), dimension(mxNode, my, mz + 1) :: zwReal
    real(kind=dbPc), dimension(mxNode, my, mz) :: pfrac
    real(kind=dbPc), dimension(mz) :: gridVolume
    real(kind=dbPc), dimension(mz) :: u
    real(kind=dbPc), dimension(mz) :: tau_p
    real(kind=dbPc), dimension(mz) :: tau_f
    real(kind=dbPc), dimension(mz) :: phi_p
    real(kind=dbPc), dimension(mz) :: F_p
    real(kind=dbPc), dimension(mz) :: F_pNode
    ! Particle bed
    integer, dimension(mxNode, my) :: rollTo, rollFrom
    real(kind=dbPc), dimension(mxNode + 1) :: xsf
    real(kind=dbPc), dimension(my) :: ysf
    real(kind=dbPc), dimension(mxNode + 1, my) :: zsf
    real(kind=dbPc), dimension(mxNode + 1, my) :: dsf
    real(kind=dbPc), dimension(npdf) :: initDiameterDist
    real(kind=dbPc), dimension(mxNode + 1, my, npdf) :: hist
    real(kind=dbPc), dimension(mxNode + 1, my) :: zChange
    real(kind=dbPc), dimension(mxNode + 1, my, npdf) :: binChange
    ! particle
    integer :: pNum
    integer, dimension(maxNum, 3) :: pIndex
    integer, dimension(maxNum) :: collCount
    real(kind=dbPc) :: xflux, zflux
    real(kind=dbPc), dimension(maxNum) :: xp, yp, zp
    real(kind=dbPc), dimension(maxNum) :: up, vp, wp
    real(kind=dbPc), dimension(maxNum) :: dp
    real(kind=dbPc), dimension(maxNum) :: hp
    real(kind=dbPc), dimension(maxNum) :: maxhp
    real(kind=dbPc), dimension(maxNum) :: survTime
    real(kind=dbPc), dimension(maxNum) :: survLength
    real(kind=dbPc), dimension(mz) :: xfluxPf
    real(kind=dbPc), dimension(mz) :: zfluxPf
    ! iteration
    integer :: iter
    real(kind=dbPc) :: time
end module public_val

module vector_cal
    use public_val
    implicit none
contains
    function dot_prod(vec1, vec2)
        real(kind=dbPc) :: dot_prod
        real(kind=dbPc), intent(in), dimension(3) :: vec1, vec2
        !
        dot_prod = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
    end function
    !
    function cross_prod(vec1, vec2)
        real(kind=dbPc), dimension(3) :: cross_prod
        real(kind=dbPc), intent(in), dimension(3) :: vec1, vec2
        !
        cross_prod(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        cross_prod(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        cross_prod(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    end function
    !
    function norm_2(vec)
        real(kind=dbPc) :: norm_2
        real(kind=dbPc), intent(in), dimension(3) :: vec
        !
        norm_2 = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
    end function
    !
    function unit_vec(vec)
        real(kind=dbPc), dimension(3) :: unit_vec
        real(kind=dbPc), intent(in), dimension(3) :: vec
        real(kind=dbPc) :: normv
        !
        normv = norm_2(vec)
        if (normv > 0.0) then
            unit_vec = vec/normv
        else
            unit_vec = 0.0
        end if
    end function
    !
    function dist_p(vec1, vec2)
        real(kind=dbPc) :: dist_p
        real(kind=dbPc), intent(in), dimension(3) :: vec1, vec2
        !
        dist_p = sqrt((vec1(1) - vec2(1))**2 + (vec1(2) - vec2(2))**2 + (vec1(3) - vec2(3))**2)
    end function
end module vector_cal

module math_operations
    use public_val
    implicit none
    private
    logical :: flag = .true.
    public :: valObeyCertainPDF, generateNormalDistHistogram, valObeyNormalDist

contains

    function valObeyCertainPDF(histogram)
        implicit none
        real(kind=dbPc) :: valObeyCertainPDF
        real(kind=dbPc), dimension(npdf), intent(in) :: histogram
        integer :: i
        real(kind=dbPc) :: rand, val, cumulativeProb
        !
        call random_number(rand)
        cumulativeProb = 0.0
        do i = 1, npdf
            cumulativeProb = cumulativeProb + histogram(i)
            if (rand <= cumulativeProb) then
                val = binStart + (i - 0.5)*binWidth
                exit
            end if
        end do
        valObeyCertainPDF = val
    end function valObeyCertainPDF

    subroutine generateNormalDistHistogram(histogram)
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram
        real(kind=dbPc) :: total

        ! Initialize the histogram array to zero
        histogram = 0.0

        ! Fill the histogram array
        do i = 1, npdf
            histogram(i) = exp(-0.5*((binStart + (i - 0.5)*binWidth - dpa)/dpStddDev)**2)/ &
                           (sqrt(2.0*pi)*dpStddDev)
            histogram(i) = histogram(i)*binWidth
        end do
        total = sum(histogram)
        histogram = histogram/total
    end subroutine generateNormalDistHistogram

    function valObeyNormalDist(mean, stdd)
        implicit none
        real(kind=dbPc), intent(in) :: mean, stdd
        real(kind=dbPc) :: valObeyNormalDist
        real(kind=dbPc) :: rand1, rand2
        real(kind=dbPc) :: r1, r2

        call random_number(rand1)
        call random_number(rand2)

        if (flag) then
            r1 = sqrt(-2.0*log(rand1))*cos(2.0*pi*rand2)
            valObeyNormalDist = mean + stdd*r1
            flag = .false.
        else
            r2 = sqrt(-2.0*log(rand1))*sin(2.0*pi*rand2)
            valObeyNormalDist = mean + stdd*r2
            flag = .true.
        end if
    end function valObeyNormalDist

end module math_operations

module gather_xyz
    implicit none
    private

    integer, parameter :: dbPc = 8 !selected_real_kind(15, 307)

    interface gatherxyz
        module procedure gxyz_real
        module procedure gxyz_int
    end interface

    public :: gatherx, gatherxy, gatherxyz

contains

    subroutine gxyz_real(comm, mxNode, mx, my, mz, f, tf)
        include "mpif.h"
        ! public
        integer, intent(in) :: mxNode, mx, my, mz
        integer, intent(in) :: comm
        real(kind=dbPc), intent(in), dimension(mxNode, my, mz) :: f
        real(kind=dbPc), dimension(mx, my, mz) :: tf
        ! local
        integer :: valtype
        integer :: arrySize, arrySizeTot
        integer :: ierr
        integer :: i, j, k
        integer :: ijk
        integer :: blk
        real(kind=dbPc), allocatable, dimension(:, :, :) :: ff
        real(kind=dbPc), allocatable, dimension(:) :: a
        real(kind=dbPc), allocatable, dimension(:) :: aa
        !
        allocate (ff(mxNode - 2, my, mz))
        arrySize = (mxNode - 2)*my*mz
        arrySizeTot = (mx - 2)*my*mz
        allocate (a(arrySize))
        allocate (aa(arrySizeTot))
        !
        valtype = MPI_DOUBLE
        a = 0.
        aa = 0.
        tf = 0.
        ff = f(2:mxNode - 1, 1:my, 1:mz)
        do k = 1, mz
        do j = 1, my
        do i = 1, mxNode - 2
            ijk = i + (j - 1)*(mxNode - 2) + (k - 1)*(mxNode - 2)*my
            a(ijk) = ff(i, j, k)
        end do
        end do
        end do
        call MPI_ALLGATHER(a, arrySize, valtype, aa, arrySize, valtype, comm, ierr)
        do k = 1, mz
        do j = 1, my
        do i = 2, mx - 1
            blk = (i - 2)/(mxNode - 2)
            ijk = (i - 1) - blk*(mxNode - 2) + (j - 1)*(mxNode - 2) + (k - 1)*(mxNode - 2)*my + blk*(mxNode - 2)*my*mz
            tf(i, j, k) = aa(ijk)
        end do
        end do
        end do
        deallocate (ff)
        deallocate (a)
        deallocate (aa)
    end subroutine gxyz_real
    !
    subroutine gxyz_int(comm, mxNode, mx, my, mz, f, tf)
        include "mpif.h"
        ! public
        integer, intent(in) :: mxNode, mx, my, mz
        integer, intent(in) :: comm
        integer, intent(in), dimension(mxNode, my, mz) :: f
        integer, dimension(mx, my, mz) :: tf
        ! local
        integer :: valtype
        integer :: arrySize, arrySizeTot
        integer :: ierr
        integer :: i, j, k
        integer :: ijk
        integer :: blk
        integer, allocatable, dimension(:, :, :) :: ff
        integer, allocatable, dimension(:) :: a
        integer, allocatable, dimension(:) :: aa

        allocate (ff(mxNode - 2, my, mz))
        arrySize = (mxNode - 2)*my*mz
        arrySizeTot = (mx - 2)*my*mz
        allocate (a(arrySize))
        allocate (aa(arrySizeTot))

        valtype = MPI_INTEGER
        a = 0.
        aa = 0.
        tf = 0.
        ff = f(2:mxNode - 1, 1:my, 1:mz)
        do k = 1, mz
        do j = 1, my
        do i = 1, mxNode - 2
            ijk = i + (j - 1)*(mxNode - 2) + (k - 1)*(mxNode - 2)*my
            a(ijk) = ff(i, j, k)
        end do
        end do
        end do
        call MPI_ALLGATHER(a, arrySize, valtype, aa, arrySize, valtype, comm, ierr)
        do k = 1, mz
        do j = 1, my
        do i = 2, mx - 1
            blk = (i - 2)/(mxNode - 2)
            ijk = (i - 1) - blk*(mxNode - 2) + (j - 1)*(mxNode - 2) + (k - 1)*(mxNode - 2)*my + blk*(mxNode - 2)*my*mz
            tf(i, j, k) = aa(ijk)
        end do
        end do
        end do
        deallocate (ff)
        deallocate (a)
        deallocate (aa)
    end subroutine gxyz_int

    subroutine gatherx(comm, mxNode, mx, f, tf)
        implicit none
        include "mpif.h"
        ! public
        integer, intent(in) :: mxNode, mx
        integer, intent(in) :: comm
        real(kind=dbPc), intent(in), dimension(mxNode) :: f
        real(kind=dbPc), dimension(mx) :: tf
        ! local
        integer :: valtype
        integer :: arrySize, arrySizeTot
        integer :: ierr
        real(kind=dbPc), allocatable, dimension(:) :: ff
        real(kind=dbPc), allocatable, dimension(:) :: aa
        !
        arrySize = mxNode - 2
        arrySizeTot = mx - 2
        allocate (ff(arrySize))
        allocate (aa(arrySizeTot))
        !
        valtype = MPI_DOUBLE
        aa = 0.
        tf = 0.
        ff = f(2:mxNode - 1)
        call MPI_ALLGATHER(ff, arrySize, valtype, aa, arrySize, valtype, comm, ierr)
        tf(2:mx - 1) = aa
        deallocate (ff)
        deallocate (aa)
    end subroutine gatherx

    subroutine gatherxy(comm, mxNode, mx, my, f, tf)
        include "mpif.h"
        ! public
        integer, intent(in) :: mxNode, mx, my
        integer, intent(in) :: comm
        real(kind=dbPc), intent(in), dimension(mxNode, my) :: f
        real(kind=dbPc), dimension(mx, my) :: tf
        ! local
        integer :: valtype
        integer :: arrySize, arrySizeTot
        integer :: ierr
        integer :: i, j
        integer :: ij
        integer :: blk
        real(kind=dbPc), allocatable, dimension(:, :) :: ff
        real(kind=dbPc), allocatable, dimension(:) :: a
        real(kind=dbPc), allocatable, dimension(:) :: aa
        !
        allocate (ff(mxNode - 2, my))
        arrySize = (mxNode - 2)*my
        arrySizeTot = (mx - 2)*my
        allocate (a(arrySize))
        allocate (aa(arrySizeTot))
        !
        valtype = MPI_DOUBLE
        a = 0.0
        aa = 0.0
        tf = 0.0
        ff = f(2:mxNode - 1, 1:my)
        do j = 1, my
            do i = 1, mxNode - 2
                ij = i + (j - 1)*(mxNode - 2)
                a(ij) = ff(i, j)
            end do
        end do
        call MPI_ALLGATHER(a, arrySize, valtype, aa, arrySize, valtype, comm, ierr)
        do j = 1, my
            do i = 2, mx - 1
                blk = (i - 2)/(mxNode - 2)
                ij = (i - 1) - blk*(mxNode - 2) + (j - 1)*(mxNode - 2) + blk*(mxNode - 2)*my
                tf(i, j) = aa(ij)
            end do
        end do
        deallocate (ff)
        deallocate (a)
        deallocate (aa)
    end subroutine gatherxy

end module gather_xyz

program main
    use public_val
    implicit none
    include "mpif.h"
    integer :: ierr
    logical :: periods
    integer :: coords
    integer :: nbrleft, nbrright
    !
    call MPI_INIT(ierr)
    call random_seed()
    realtype = MPI_DOUBLE
    inttype = MPI_INTEGER
    ! create MPI Cartesian topology
    procs = nNodes
    periods = .true.
    call MPI_CART_CREATE(MPI_COMM_WORLD, 1, procs, periods, .true., comm, ierr)
    call MPI_COMM_RANK(comm, myID, ierr)
    call MPI_CART_GET(comm, 1, procs, periods, coords, ierr)
    ! find the neighbors
    !
    !       |           |
    !       |           |
    !     -----------------
    !       |           |
    !       |           |
    !     1 |   myID    | 2
    !       |           |
    !       |           |
    !     -----------------
    !       |           |
    !       |           |
    !
    call MPI_CART_SHIFT(comm, 0, 1, nbrleft, nbrright, ierr)
    neighbor(1) = nbrleft
    neighbor(2) = nbrright
    call MPI_BARRIER(comm, ierr)
    ! find indices of subdomain and check that dimensions of arrays are sufficient
    if (mod(mx - 2, nNodes) /= 0) then
        print *, 'mx-2 cannot diveded by nNodes'
        stop
    end if
    if (mod(mx - 2, nNodes) /= 0) then
        print *, 'mx-2 cannot diveded by nNodes'
        stop
    end if
    call generateSurfGrid
    call initializeSurface
    call generateFieldGrid
    call initializeField
    call initializeParticle
    call determineParIJK
    call defineExchangeType
    call generateOutputFile
    iter = 0
    time = 0.0
    do while (time <= endTime)
        time = time + dt
        iter = iter + 1
        if (iter >= parStart) then
            call determineParRollDir
            call calculateSplash
            call calculateParMov
            call reallocateParticle
            call determineParIJK
            if (midairCollision) then
                call calculateMidAirColl
            end if
            call updateSurfGrid
            call updateFieldGrid
        end if
        call calculateFluidField
        call outputAll
    end do
    stop
    call MPI_FINALIZE(ierr)
end program main

subroutine generateSurfGrid
    use public_val
    implicit none
    include "mpif.h"
    integer :: i, j
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)

    zsf = initSurfElevation
    if (myID == 0) then
        xsf(1) = -xDiff
        do i = 2, mxNode + 1
            xsf(i) = xsf(i - 1) + xDiff
        end do
        call MPI_SEND(xsf(mxNode - 1), 1, realtype, neighbor(2), 3, comm, ierr)
    else
        call MPI_RECV(xsf(1), 1, realtype, neighbor(1), 3, comm, status, ierr)
        do i = 2, mxNode + 1
            xsf(i) = xsf(i - 1) + xDiff
        end do
        call MPI_SEND(xsf(mxNode - 1), 1, realtype, neighbor(2), 3, comm, ierr)
    end if
    ysf(1) = -yDiff
    do j = 2, my
        ysf(j) = ysf(j - 1) + yDiff
    end do
end subroutine generateSurfGrid

subroutine initializeSurface
    use public_val
    use math_operations
    implicit none
    integer :: i, j, n

    ! Initialize the diameter distribution
    select case (whichDiameterDist)
    case (0)
        if (npdf >= 3) then
            call generateNormalDistHistogram(initDiameterDist)
        else
            print *, 'Error: npdf must be greater than or equal to 3 for normal distribution'
            stop
        end if
    case (1)
        if (npdf == 1) then
            initDiameterDist = 1.0
        else
            print *, 'Error: npdf must be equal to 1 for uniform diameter distribution'
            stop
        end if
    case (2)
        if (npdf == 2) then
            initDiameterDist(1) = prob1
            initDiameterDist(2) = 1.0 - prob1
        else
            print *, 'Error: npdf must be equal to 2 for Bernoulli distribution'
            stop
        end if
    end select
    ! Initialize the particle bed
    do j = 1, my
        do i = 1, mxNode + 1
            if (predefineSurface) then
                zsf(i, j) = amp*sin(omg*xsf(i)) + initSurfElevation
            end if
            dsf(i, j) = dpa
            do n = 1, npdf
                hist(i, j, n) = initDiameterDist(n)
            end do
        end do
    end do
end subroutine initializeSurface

subroutine generateFieldGrid
    use public_val
    implicit none
    integer :: i, j, k, n
    real(kind=dbPc) :: zDiffMax, refineRatio

    zw(mz) = 0.0
    zDiff(mz) = 0.0
    n = 0
    do while (zMax - zw(mz) > zDiff(mz)*2.0)
        zDiffMax = zMax/(mz - n)
        refineRatio = (zDiffMax/zDiffMin)**(1.0/(mzUni - 1))
        zDiff(1) = zDiffMin
        zw(1) = 0.0
        z(1) = zDiff(1)*0.5
        do k = 2, mz
            if (k <= mzUni) then
                zDiff(k) = zDiff(k - 1)*refineRatio
            else
                zDiff(k) = zDiffMax
            end if
            zw(k) = zw(k - 1) + zDiff(k - 1)
            z(k) = zw(k) + zDiff(k)*0.5
        end do
        n = n + 1
    end do
    zDiff(mz) = zMax - zw(mz)
    z(mz) = zMax - zDiff(mz)*0.5
    zw(mz + 1) = zMax

    do i = 1, mxNode
        xu(i) = xsf(i)
        x(i) = xu(i) + xDiff*0.5
    end do
    do j = 1, my
        yv(j) = ysf(j)
        y(j) = yv(j) + yDiff*0.5
    end do

    do k = 1, mz
        do j = 1, my
            do i = 1, mxNode
                zwReal(i, j, k) = zw(k) + zsf(i, j)
                zReal(i, j, k) = z(k) + zsf(i, j)
            end do
        end do
    end do
    do j = 1, my
        do i = 1, mxNode
            zwReal(i, j, mz + 1) = zMax + zsf(i, j)
        end do
    end do
end subroutine generateFieldGrid

subroutine initializeField
    use public_val
    implicit none
    integer :: k

    pfrac = 0.0
    gridVolume = xDiff*yDiff*zDiff
    tau_p = 0.0
    tau_f = rho*uStar**2
    F_pNode = 0.0
    do k = 1, mz
        if (z(k) > z0) then
            u(k) = uStar/kapa*log(z(k)/z0)
        else
            u(k) = 0.0
        end if
    end do
end subroutine initializeField

subroutine initializeParticle
    use public_val
    use math_operations
    implicit none
    integer :: n
    real(kind=dbPc) :: rand1, rand2, rand3

    if (binStart + 0.5*binWidth < 0.0) then
        print *, 'Error: particle diameter must be greater than or equal to 0'
        stop
    end if
    pNum = pNumInit
    do n = 1, pNum
        call random_number(rand1)
        call random_number(rand2)
        call random_number(rand3)
        xp(n) = (xu(mxNode) - xu(2))*rand1 + xu(2)
        yp(n) = (yv(my) - yv(2))*rand2 + yv(2)
        hp(n) = 0.05*rand3
        maxhp(n) = hp(n)
        zp(n) = hp(n) + initSurfElevation
        up(n) = 0.0
        vp(n) = 0.0
        wp(n) = 0.0
        dp(n) = valObeyCertainPDF(initDiameterDist)
        collCount(n) = 0
        survTime(n) = 0.0
        survLength(n) = 0.0
    end do
end subroutine initializeParticle

subroutine determineParIJK
    use public_val
    implicit none

    integer :: n, ip, jp, kp
    integer :: tempKp
    real(kind=dbPc) :: currentZ

    do n = 1, pNum
        ip = floor((xp(n) - xu(2))/xDiff) + 2
        jp = floor((yp(n) - yv(2))/yDiff) + 2
        if (ip < 2) ip = 2
        if (jp < 2) jp = 2
        if (ip > mxNode - 1) ip = mxNode - 1
        if (jp > my - 1) jp = my - 1
        do tempKp = 1, mz
            kp = tempKp
            currentZ = zwReal(ip, jp, tempKp + 1)
            if (zp(n) < currentZ) then
                exit
            end if
        end do
        pIndex(n, 1) = ip
        pIndex(n, 2) = jp
        pIndex(n, 3) = kp
    end do
end subroutine determineParIJK

subroutine defineExchangeType
    use public_val
    implicit none
    include "mpif.h"
    integer dataLength, tempType1, tempType2, ierr
    ! sliceType: i=const planes
    ! datatype for i=const,k=const line
    call MPI_TYPE_VECTOR(my, 1, mxNode, realtype, tempType1, ierr)
    call MPI_TYPE_COMMIT(tempType1, ierr)
    ! datatype for i=const plane
    call MPI_TYPE_EXTENT(realtype, dataLength, ierr)
    call MPI_TYPE_HVECTOR(mz, 1, mxNode*my*dataLength, tempType1, sliceType, ierr)
    call MPI_TYPE_COMMIT(sliceType, ierr)
    ! datatype for i=const line
    call MPI_TYPE_VECTOR(my, 1, mxNode + 1, realtype, edgeType1, ierr)
    call MPI_TYPE_COMMIT(edgeType1, ierr)

    call MPI_TYPE_VECTOR(my, 1, mxNode + 1, realtype, tempType2, ierr)
    call MPI_TYPE_COMMIT(tempType2, ierr)
    call MPI_TYPE_HVECTOR(npdf, 1, (mxNode + 1)*my*dataLength, tempType2, edgeType2, ierr)
    call MPI_TYPE_COMMIT(edgeType2, ierr)

end subroutine defineExchangeType

subroutine generateOutputFile
    use public_val
    implicit none
    character(len=32) bashCmd
    external :: system

    if (myID == 0) then
        bashCmd = 'rm -rf Particle'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'rm -rf Field'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'rm -rf Surface'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'rm -rf Statistics'
        call system(trim(adjustl(bashCmd)))

        bashCmd = 'mkdir Particle'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Field'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Surface'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Statistics'
        call system(trim(adjustl(bashCmd)))

        open (unit=10, file='./Particle/ParticleNum.plt')
        write (10, *) 'variables = "t", "Number", "iteration"'
        close (10)

        open (unit=11, file='./Particle/ParticleData_0.plt')
        write (11, "(A100)") 'variables = "x", "y", "z", "u", "v", "w", "u_mag", "h", "h_max", "d", "L", "t", "n_coll"'
        close (11)

        open (unit=12, file='./Field/Profile.plt')
        write (12, *) 'variables = "z", "u", "tau_p", "tau_f", "F_p", "phi_p"'
        close (12)

        open (unit=13, file='./Field/FieldData_0.plt')
        write (13, *) 'variables = "x", "y", "z", "u", "Phi_p"'
        close (13)

        open (unit=14, file='./Surface/SurfaceData_0.plt')
        write (14, *) 'variables = "x", "y", "z", "d"'
        close (14)

        open (unit=15, file='./Statistics/vsTime.plt')
        write (15, *) 'variables = "t", "Q", "Psi", "HopLength"'
        close (15)

        open (unit=16, file='./Statistics/vsHeight.plt')
        write (16, *) 'variables = "z", "Q", "Psi", "Phi_p", "betterQ"'
        close (16)
    end if
end subroutine generateOutputFile

subroutine determineParRollDir
    use public_val
    implicit none
    integer :: i, j, k
    real(kind=dbPc) :: centerZ, westZ, eastZ, northZ, southZ
    real(kind=dbPc), dimension(4) :: slopes
    real(kind=dbPc) :: maxSlope

    ! Determine the direction of particle rolling
    ! avalanchTo and avalanchFrom are defined as follows:
    ! avalanchTo: 1 -> east, 2 -> west, 3 -> north, 4 -> south
    !       3 ^
    !         |
    !  2 <----0----> 1
    !         |
    !         v 4

    ! avalanchFrom: 1 -> east, 2 -> west, 3 -> north, 4 -> south
    !       3 |
    !         v
    !  2 ---->0<---- 1
    !         ^
    !         | 4

    rollTo = 0
    rollFrom = 0
    do j = 2, my
        do i = 2, mxNode
            centerZ = zsf(i, j)
            eastZ = zsf(i + 1, j)
            westZ = zsf(i - 1, j)
            if (j < my) then
                northZ = zsf(i, j + 1)
            else
                northZ = zsf(i, 3)
            end if
            southZ = zsf(i, j - 1)

            slopes(1) = (centerZ - eastZ)/xDiff
            slopes(2) = (centerZ - westZ)/xDiff
            slopes(3) = (centerZ - northZ)/yDiff
            slopes(4) = (centerZ - southZ)/yDiff

            maxSlope = tanRepostAngle
            do k = 1, 4
                if (slopes(k) > maxSlope) then
                    maxSlope = slopes(k)
                    rollTo(i, j) = k
                end if
            end do

            slopes(1) = (eastZ - centerZ)/xDiff
            slopes(2) = (westZ - centerZ)/xDiff
            slopes(3) = (northZ - centerZ)/yDiff
            slopes(4) = (southZ - centerZ)/yDiff

            maxSlope = tanRepostAngle
            do k = 1, 4
                if (slopes(k) > maxSlope) then
                    maxSlope = slopes(k)
                    rollFrom(i, j) = k
                end if
            end do
        end do
    end do
end subroutine determineParRollDir

subroutine calculateSplash
    use public_val
    use vector_cal
    use math_operations
    implicit none

    integer :: tempNum
    integer :: nAddGlobal
    integer :: n
    integer :: ip, jp
    integer :: whichTri
    integer :: whichVer
    integer :: ipp, jpp
    integer :: iterAng
    logical :: rebound
    integer :: ii, jj
    integer :: iBin
    integer :: ejectNum
    integer :: nadd
    real(kind=dbPc) :: estimateAltitude
    real(kind=dbPc) :: localXP, localYP
    real(kind=dbPc) :: d1, d2
    real(kind=dbPc) :: dd1, dd2
    real(kind=dbPc) :: m1, m2
    real(kind=dbPc) :: v1, v2
    real(kind=dbPc) :: eta, alpha, beta, gama, lambda, sigma, mu
    real(kind=dbPc) :: eBar
    real(kind=dbPc) :: angin1
    real(kind=dbPc) :: angout1, angout1Min, angout1Max, angout1Mid
    real(kind=dbPc) :: maxPDF
    real(kind=dbPc) :: rr1, rr2, rr3
    real(kind=dbPc) :: xPDF, yPDF, pdf
    real(kind=dbPc) :: E1, E2, Ed1, Ed2, Eeff
    real(kind=dbPc) :: E2Bar
    real(kind=dbPc) :: merfc
    real(kind=dbPc) :: ejectVol, rollVol
    real(kind=dbPc) :: vch
    real(kind=dbPc) :: tau_s
    real(kind=dbPc), dimension(npdf) :: currentHist
    real(kind=dbPc), dimension(maxEjectNum) :: tempx, tempy, tempz
    real(kind=dbPc), dimension(maxEjectNum) :: tempw
    real(kind=dbPc), dimension(maxEjectNum) :: tempd
    real(kind=dbPc), dimension(maxEjectNum) :: tempi, tempj
    real(kind=dbPc), dimension(4) :: adjacentZ
    real(kind=dbPc), dimension(3) :: point0, point1, point2, point3
    real(kind=dbPc), dimension(3) :: vector12, vector13
    real(kind=dbPc), dimension(3) :: nVec
    real(kind=dbPc), dimension(3) :: pVol1, pVol2
    real(kind=dbPc), dimension(3) :: localCoordX, localCoordY, localCoordZ
    real(kind=dbPc), dimension(3) :: distance
    real(kind=dbPc), dimension(3) :: vin
    real(kind=dbPc), dimension(3) :: vout
    real(kind=dbPc), dimension(3) :: vec1, vec2, vec3

    zChange = 0.0
    binChange = 0.0
    tempNum = 0
    nAddGlobal = 0
    xflux = 0.0
    zflux = 0.0
    do n = 1, pNum
        ip = pIndex(n, 1)
        jp = pIndex(n, 2)
        adjacentZ(1) = zsf(ip, jp)
        adjacentZ(2) = zsf(ip + 1, jp)
        adjacentZ(3) = zsf(ip, jp + 1)
        adjacentZ(4) = zsf(ip + 1, jp + 1)
        estimateAltitude = zp(n) - maxval(adjacentZ)
        if (estimateAltitude > 0.0) then
            tempNum = tempNum + 1
            xp(tempNum) = xp(n)
            yp(tempNum) = yp(n)
            zp(tempNum) = zp(n)
            hp(tempNum) = estimateAltitude
            up(tempNum) = up(n)
            vp(tempNum) = vp(n)
            wp(tempNum) = wp(n)
            dp(tempNum) = dp(n)
            collCount(tempNum) = collCount(n)
            survTime(tempNum) = survTime(n)
            survLength(tempNum) = survLength(n)
            maxhp(tempNum) = maxhp(n)
            pIndex(tempNum, :) = pIndex(n, :)
            cycle
        end if
        localXP = xp(n) - xsf(ip)
        localYP = yp(n) - ysf(jp)
        if (mod(ip + jp, 2) == 0) then
            ! ip+jp = Even num
            !   3 ------- 2
            !     |    /|
            !     |  /  |
            !     |/    |
            !   1 ------- 3
            ! Left triangle: whichTriangle=1
            ! Right triangle: whichTriangle=2
            ! Current surfGrid%location is at point 1, it is the local origin
            point1 = [xsf(ip), ysf(jp), zsf(ip, jp)]
            point2 = [xsf(ip + 1), ysf(jp + 1), zsf(ip + 1, jp + 1)]
            if (localXP/localYP <= xDiff/yDiff) then
                whichTri = 1
                point3 = [xsf(ip), ysf(jp + 1), zsf(ip, jp + 1)]
                vector12 = point2 - point1
                vector13 = point3 - point1
                nVec = unit_vec(cross_prod(vector12, vector13))
            else
                whichTri = 2
                point3 = [xsf(ip + 1), ysf(jp), zsf(ip + 1, jp)]
                vector12 = point2 - point1
                vector13 = point3 - point1
                nVec = unit_vec(cross_prod(vector13, vector12))
            end if
        else
            ! ip+jp = Odd num
            !   2 ------- 3
            !     |\    |
            !     |  \  |
            !     |    \|
            !   3 ------- 1
            ! Left triangle: whichTriangle=3
            ! Right triangle: whichTriangle=4
            ! Current surfGrid%location is at point 3, it is the local origin
            point1 = [xsf(ip + 1), ysf(jp), zsf(ip + 1, jp)]
            point2 = [xsf(ip), ysf(jp + 1), zsf(ip, jp + 1)]
            if (localXP/(yDiff - localYP) <= xDiff/yDiff) then
                whichTri = 3
                point3 = [xsf(ip), ysf(jp), zsf(ip, jp)]
                vector12 = point2 - point1
                vector13 = point3 - point1
                nVec = unit_vec(cross_prod(vector12, vector13))
            else
                whichTri = 4
                point3 = [xsf(ip + 1), ysf(jp + 1), zsf(ip + 1, jp + 1)]
                vector12 = point2 - point1
                vector13 = point3 - point1
                nVec = unit_vec(cross_prod(vector13, vector12))
            end if
        end if
        point0(1) = xp(n)
        point0(2) = yp(n)
        point0(3) = point3(3) - (nVec(1)*(point0(1) - point3(1)) + nVec(2)*(point0(2) - point3(2)))/nVec(3)
        hp(n) = zp(n) - point0(3)
        if (hp(n) <= 0.0 .and. wp(n) < 0.0) then
            pVol1(1) = up(n)
            pVol1(2) = vp(n)
            pVol1(3) = wp(n)
            localCoordZ = nVec
            localCoordY = unit_vec(cross_prod(nVec, pVol1))
            localCoordX = unit_vec(cross_prod(localCoordY, nVec))
            vin(1) = dot_prod(pVol1, localCoordX)
            vin(2) = dot_prod(pVol1, localCoordY)
            vin(3) = dot_prod(pVol1, localCoordZ)
            distance(1) = dist_p(point0, point1)
            distance(2) = dist_p(point0, point2)
            distance(3) = dist_p(point0, point3)
            whichVer = minloc(distance, 1)
            select case (whichTri)
            case (1)
                select case (whichVer)
                case (1)
                    ipp = ip
                    jpp = jp
                case (2)
                    ipp = ip + 1
                    jpp = jp + 1
                case (3)
                    ipp = ip
                    jpp = jp + 1
                end select
            case (2)
                select case (whichVer)
                case (1)
                    ipp = ip
                    jpp = jp
                case (2)
                    ipp = ip + 1
                    jpp = jp + 1
                case (3)
                    ipp = ip + 1
                    jpp = jp
                end select
            case (3)
                select case (whichVer)
                case (1)
                    ipp = ip + 1
                    jpp = jp
                case (2)
                    ipp = ip
                    jpp = jp + 1
                case (3)
                    ipp = ip
                    jpp = jp
                end select
            case (4)
                select case (whichVer)
                case (1)
                    ipp = ip + 1
                    jpp = jp
                case (2)
                    ipp = ip
                    jpp = jp + 1
                case (3)
                    ipp = ip + 1
                    jpp = jp + 1
                end select
            end select
            d1 = dp(n)
            d2 = dsf(ipp, jpp)
            dd1 = d1/(0.5*(d1 + d2))
            dd2 = d2/(0.5*(d1 + d2))
            m1 = (pi*d1**3)/6.0*rhoP
            m2 = (pi*d2**3)/6.0*rhoP
            v1 = norm_2(pVol1)
            eta = resN*dd1**3/(dd1**3 + resN*dd2**3)
            alpha = (1.0 + resN)/(1.0 + eta) - 1.0
            beta = 1.0 - (2.0/7.0)*(1.0 - resT)/(1.0 + eta)
            angin1 = atan(abs(vin(3)/vin(1)))
            eBar = beta - (beta**2 - alpha**2)*dd2*angin1/(2.0*beta)
            v2 = eBar*v1
            gama = 4.0/9.0*beta**2/(alpha + beta)**2/dd2
            angout1Min = -angin1
            angout1Max = sqrt(angin1/gama)*2.0 - angin1
            if (angout1Max > pi) angout1Max = pi
            angout1Mid = 0.5*(angout1Min + angout1Max)
            maxPDF = gama*(angin1 + angout1Mid)/angin1*log(2.0*angin1/gama/(angin1 + angout1Mid)**2)*1.2
            iterAng = 0
            do
                iterAng = iterAng + 1
                call random_number(rr1)
                call random_number(rr2)
                xPDF = rr1*(angout1Max - angout1Min) + angout1Min
                yPDF = rr2*min(maxPDF, 1.0)
                pdf = gama*(xPDF + angin1)/angin1*log(2.0*angin1/gama/(xPDF + angin1)**2)
                if (yPDF <= pdf .or. iterAng > 10000) then
                    angout1 = xPDF + angin1
                    exit
                end if
            end do
            vout(3) = v2*sin(angout1)
            E2 = 0.5*m1*vout(3)**2
            Ed1 = m1*9.8*0.5*(d1 + d2)
            if (E2 < Ed1 .or. eBar <= 0.0 .or. angout1 <= 0.0) then
                rebound = .false.
                !eBar = 0.0
            else
                rebound = .true.
            end if
            if (rebound) then
                tempNum = tempNum + 1
                vout(1) = v2*cos(angout1)
                vout(2) = 0.0
                vec1 = vout(1)*localCoordX
                vec2 = vout(2)*localCoordY
                vec3 = vout(3)*localCoordZ
                pVol2 = vec1 + vec2 + vec3
                xp(tempNum) = point0(1)
                yp(tempNum) = point0(2)
                zp(tempNum) = point0(3)
                hp(tempNum) = 0.0
                up(tempNum) = pVol2(1)
                vp(tempNum) = pVol2(2)
                wp(tempNum) = pVol2(3)
                dp(tempNum) = d1
                survLength(tempNum) = 0.0
                survTime(tempNum) = 0.0
                collCount(tempNum) = 0
                maxhp(tempNum) = 0.0
                pIndex(tempNum, 1) = ip
                pIndex(tempNum, 2) = jp
                pIndex(tempNum, 3) = 1
                zflux = zflux + m1/area/dt
            else
                select case (rollTo(ipp, jpp))
                case (0)
                    ii = ipp
                    jj = jpp
                case (1)
                    ii = ipp + 1
                    jj = jpp
                case (2)
                    ii = ipp - 1
                    jj = jpp
                case (3)
                    ii = ipp
                    jj = jpp + 1
                case (4)
                    ii = ipp
                    jj = jpp - 1
                end select
                if (jj > my) jj = 3
                if (whichDiameterDist /= 1) then
                    iBin = floor((d1 - binStart)/binWidth) + 1
                    iBin = max(iBin, 1)
                    iBin = min(iBin, npdf)
                else
                    iBin = 1
                end if
                vch = (pi*d1**3)/6.0
                zChange(ii, jj) = zChange(ii, jj) + vch
                binChange(ii, jj, iBin) = binChange(ii, jj, iBin) + vch
            end if
            E1 = 0.5*m1*v1**2
            Ed2 = m2*9.8*d2
            tau_s = rho*0.0123*(rhoP/rho*9.8*d2 + 3.0e-4/(rho*d2))
            Eeff = Ed2*(1.0 - tau_f(1)/tau_s)
            Eeff = max(Eeff, 0.1*Ed2)
            lambda = 2.0*log((1.0 - eBar**2)*E1/Eeff)
            sigma = sqrt(lambda)*log(2.0)
            mu = log((1.0 - eBar**2)*E1) - lambda*log(2.0)
            E2Bar = Eeff*((1.0 - eBar**2)*E1/Eeff)**(1.0 - (2.0 - log(2.0))*log(2.0))
            merfc = erfc((log(Eeff) - mu)/(sqrt(2.0)*sigma))
            ejectNum = int(0.06*((1.0 - eBar**2)*E1/(2.0*E2Bar))*merfc)
            if (ejectNum > 0) then
                ejectVol = 0.0
                do nadd = 1, ejectNum
                    nAddGlobal = nAddGlobal + 1
                    E2 = exp(valObeyNormalDist(mu, sigma))
                    currentHist = hist(ipp, jpp, :)
                    d2 = valObeyCertainPDF(currentHist)
                    m2 = (pi*d2**3)/6.0*rhoP
                    v2 = sqrt(2.0*E2/m2)
                    rr3 = -1.0
                    do while (rr3 < 0.0)
                        call random_number(rr1)
                        call random_number(rr2)
                        rr3 = 1.0 - rr1 - rr2
                    end do
                    tempx(nAddGlobal) = rr1*point1(1) + rr2*point2(1) + rr3*point3(1)
                    tempy(nAddGlobal) = rr1*point1(2) + rr2*point2(2) + rr3*point3(2)
                    tempz(nAddGlobal) = rr1*point1(3) + rr2*point2(3) + rr3*point3(3)
                    tempw(nAddGlobal) = v2
                    tempd(nAddGlobal) = d2
                    tempi(nAddGlobal) = ip
                    tempj(nAddGlobal) = jp
                    zflux = zflux + m2/area/dt
                    vch = (pi*d2**3)/6.0
                    ejectVol = ejectVol + vch
                    if (whichDiameterDist /= 1) then
                        iBin = floor((d2 - binStart)/binWidth) + 1
                        iBin = max(iBin, 1)
                        iBin = min(iBin, npdf)
                    else
                        iBin = 1
                    end if
                    binChange(ipp, jpp, iBin) = binChange(ipp, jpp, iBin) - vch
                    zChange(ipp, jpp) = zChange(ipp, jpp) - vch
                end do
                if (rollFrom(ipp, jpp) /= 0) then
                    select case (rollFrom(ipp, jpp))
                    case (1)
                        ii = ipp + 1
                        jj = jpp
                    case (2)
                        ii = ipp - 1
                        jj = jpp
                    case (3)
                        ii = ipp
                        jj = jpp + 1
                    case (4)
                        ii = ipp
                        jj = jpp - 1
                    end select
                    if (jj > my) jj = 3
                    rollVol = 0.0
                    do while (rollVol < ejectVol)
                        currentHist = hist(ii, jj, :)
                        d2 = valObeyCertainPDF(currentHist)
                        vch = (pi*d2**3)/6.0
                        rollVol = rollVol + vch
                        if (whichDiameterDist /= 1) then
                            iBin = floor((d2 - binStart)/binWidth) + 1
                            iBin = max(iBin, 1)
                            iBin = min(iBin, npdf)
                        else
                            iBin = 1
                        end if
                        zChange(ii, jj) = zChange(ii, jj) - vch
                        zChange(ipp, jpp) = zChange(ipp, jpp) + vch
                        binChange(ii, jj, iBin) = binChange(ii, jj, iBin) - vch
                        binChange(ipp, jpp, iBin) = binChange(ipp, jpp, iBin) + vch
                    end do
                end if
            end if
        else
            tempNum = tempNum + 1
            xp(tempNum) = xp(n)
            yp(tempNum) = yp(n)
            zp(tempNum) = zp(n)
            hp(tempNum) = hp(n)
            up(tempNum) = up(n)
            vp(tempNum) = vp(n)
            wp(tempNum) = wp(n)
            dp(tempNum) = dp(n)
            collCount(tempNum) = collCount(n)
            survLength(tempNum) = survLength(n)
            survTime(tempNum) = survTime(n)
            maxhp(tempNum) = maxhp(n)
            pIndex(tempNum, :) = pIndex(n, :)
        end if
    end do
    if (nAddGlobal > 0) then
        pNum = tempNum + nAddGlobal
        if (pNum > maxNum) then
            print *, pNum, maxNum
            print *, "particle number reach the limit"
            stop
        else
            xp(tempNum + 1:tempNum + nAddGlobal) = tempx(1:nAddGlobal)
            yp(tempNum + 1:tempNum + nAddGlobal) = tempy(1:nAddGlobal)
            zp(tempNum + 1:tempNum + nAddGlobal) = tempz(1:nAddGlobal)
            hp(tempNum + 1:tempNum + nAddGlobal) = 0.0
            up(tempNum + 1:tempNum + nAddGlobal) = 0.0
            vp(tempNum + 1:tempNum + nAddGlobal) = 0.0
            wp(tempNum + 1:tempNum + nAddGlobal) = tempw(1:nAddGlobal)
            dp(tempNum + 1:tempNum + nAddGlobal) = tempd(1:nAddGlobal)
            survLength(tempNum + 1:tempNum + nAddGlobal) = 0.0
            survTime(tempNum + 1:tempNum + nAddGlobal) = 0.0
            collCount(tempNum + 1:tempNum + nAddGlobal) = 0
            maxhp(tempNum + 1:tempNum + nAddGlobal) = 0.0
            pIndex(tempNum + 1:tempNum + nAddGlobal, 1) = tempi(1:nAddGlobal)
            pIndex(tempNum + 1:tempNum + nAddGlobal, 2) = tempj(1:nAddGlobal)
            pIndex(tempNum + 1:tempNum + nAddGlobal, 3) = 1
        end if
    else
        pNum = tempNum
    end if
end subroutine calculateSplash

subroutine calculateParMov
    use public_val
    implicit none
    integer :: i, j, k, n
    real(kind=dbPc) :: diaP, volP, mP
    real(kind=dbPc), dimension(3) :: ufp, fDrag
    real(kind=dbPc), dimension(3) :: bulkForce
    real(kind=dbPc), dimension(3) :: tempU1, tempU2, tempU3, tempU4
    real(kind=dbPc), dimension(3) :: tempA1, tempA2, tempA3, tempA4

    pfrac = 0.0
    F_pNode = 0.0
    do n = 1, pNum
        tempU1(1) = up(n)
        tempU1(2) = vp(n)
        tempU1(3) = wp(n)
        diaP = dp(n)
        volP = (pi*diaP**3)/6.0
        mP = rhoP*volP
        i = pIndex(n, 1)
        j = pIndex(n, 2)
        k = pIndex(n, 3)
        ufp(1) = u(k)
        ufp(2) = 0.0
        ufp(3) = 0.0
        bulkForce(1) = 0.0
        bulkForce(2) = 0.0
        bulkForce(3) = -9.8*(1.0 - rho/rhoP)
        fDrag = dragForce(tempU1)
        tempA1 = fDrag/mP + bulkForce
        tempU2 = tempU1 + 0.5*tempA1*dt
        tempA2 = dragForce(tempU2)/mP + bulkForce
        tempU3 = tempU1 + 0.5*tempA2*dt
        tempA3 = dragForce(tempU3)/mP + bulkForce
        tempU4 = tempU1 + tempA3*dt
        tempA4 = dragForce(tempU4)/mP + bulkForce
        xp(n) = xp(n) + (tempU1(1) + 2.0*tempU2(1) + 2.0*tempU3(1) + tempU4(1))/6.0*dt
        yp(n) = yp(n) + (tempU1(2) + 2.0*tempU2(2) + 2.0*tempU3(2) + tempU4(2))/6.0*dt
        zp(n) = zp(n) + (tempU1(3) + 2.0*tempU2(3) + 2.0*tempU3(3) + tempU4(3))/6.0*dt
        hp(n) = hp(n) + (tempU1(3) + 2.0*tempU2(3) + 2.0*tempU3(3) + tempU4(3))/6.0*dt
        up(n) = up(n) + (tempA1(1) + 2.0*tempA2(1) + 2.0*tempA3(1) + tempA4(1))/6.0*dt
        vp(n) = vp(n) + (tempA1(2) + 2.0*tempA2(2) + 2.0*tempA3(2) + tempA4(2))/6.0*dt
        wp(n) = wp(n) + (tempA1(3) + 2.0*tempA2(3) + 2.0*tempA3(3) + tempA4(3))/6.0*dt
        maxhp(n) = max(maxhp(n), hp(n))
        survTime(n) = survTime(n) + dt
        survLength(n) = survLength(n) + (tempU1(1) + 2.0*tempU2(1) + 2.0*tempU3(1) + tempU4(1))/6.0*dt
        F_pNode(k) = F_pNode(k) + fDrag(1)
        pfrac(i, j, k) = pfrac(i, j, k) + volP/gridVolume(k)
        xflux = xflux + up(n)*mP/area
        xfluxPf(k) = xfluxPf(k) + up(n)*mP/(area*zDiff(k))
        zfluxPf(k) = zfluxPf(k) + wp(n)*mP/(area*zDiff(k))
    end do

contains

    function dragForce(upp) result(ffd)
        implicit none
        ! public
        real(kind=dbPc), dimension(3), intent(in) :: upp
        real(kind=dbPc), dimension(3) :: ffd
        ! local
        integer :: ii
        real(kind=dbPc), dimension(3) :: relativeU, C_d, Re_p
        !
        do ii = 1, 3
            relativeU(ii) = ufp(ii) - upp(ii)
            if (relativeU(ii) /= 0.0) then
                Re_p(ii) = abs(relativeU(ii))*diaP/nu
                C_d(ii) = (sqrt(0.5) + sqrt(24.0/Re_p(ii)))**2
                ffd(ii) = pi/8.0*rho*diaP**2*C_d(ii)*abs(relativeU(ii))*relativeU(ii)
            else
                ffd(ii) = 0.0
            end if
        end do
    end function dragForce
end subroutine calculateParMov

subroutine reallocateParticle
    use public_val
    implicit none
    include 'mpif.h'

    integer :: n, ierr
    integer :: tempNum
    integer :: nESend, nWSend
    integer :: nERecv, nWRecv
    integer, dimension(pNumExchMax) :: ccE
    integer, dimension(pNumExchMax) :: ccW
    integer :: status(MPI_STATUS_SIZE)
    real(kind=dbPc), dimension(pNumExchMax) :: xpE, ypE, zpE
    real(kind=dbPc), dimension(pNumExchMax) :: upE, vpE, wpE
    real(kind=dbPc), dimension(pNumExchMax) :: dpE, hpE, mhE
    real(kind=dbPc), dimension(pNumExchMax) :: slE, stE
    real(kind=dbPc), dimension(pNumExchMax) :: xpW, ypW, zpW
    real(kind=dbPc), dimension(pNumExchMax) :: upW, vpW, wpW
    real(kind=dbPc), dimension(pNumExchMax) :: dpW, hpW, mhW
    real(kind=dbPc), dimension(pNumExchMax) :: slW, stW
    real(kind=dbPc), allocatable, dimension(:) :: exchESend, exchWSend
    real(kind=dbPc), allocatable, dimension(:) :: exchERecv, exchWRecv
    real(kind=dbPc), allocatable, dimension(:) :: exchESendi, exchWSendi
    real(kind=dbPc), allocatable, dimension(:) :: exchERecvi, exchWRecvi

    tempNum = 0
    nESend = 0
    nWSend = 0
    do n = 1, pNum
        if (xp(n) >= xu(mxNode)) then
            nESend = nESend + 1
            if (xp(n) >= xMax) then
                xpE(nESend) = xp(n) - xMax
            else
                xpE(nESend) = xp(n)
            end if
            ypE(nESend) = yp(n)
            zpE(nESend) = zp(n)
            upE(nESend) = up(n)
            vpE(nESend) = vp(n)
            wpE(nESend) = wp(n)
            dpE(nESend) = dp(n)
            hpE(nESend) = hp(n)
            mhE(nESend) = maxhp(n)
            slE(nESend) = survLength(n)
            stE(nESend) = survTime(n)
            ccE(nESend) = collCount(n)
        else if (xp(n) < xu(2)) then
            nWSend = nWSend + 1
            if (xp(n) < 0.0) then
                xpW(nWSend) = xp(n) + xMax
            else
                xpW(nWSend) = xp(n)
            end if
            ypW(nWSend) = yp(n)
            zpW(nWSend) = zp(n)
            upW(nWSend) = up(n)
            vpW(nWSend) = vp(n)
            wpW(nWSend) = wp(n)
            dpW(nWSend) = dp(n)
            hpW(nWSend) = hp(n)
            mhW(nWSend) = maxhp(n)
            slW(nWSend) = survLength(n)
            stW(nWSend) = survTime(n)
            ccW(nWSend) = collCount(n)
        else
            tempNum = tempNum + 1
            xp(tempNum) = xp(n)
            yp(tempNum) = yp(n)
            zp(tempNum) = zp(n)
            up(tempNum) = up(n)
            vp(tempNum) = vp(n)
            wp(tempNum) = wp(n)
            dp(tempNum) = dp(n)
            hp(tempNum) = hp(n)
            maxhp(tempNum) = maxhp(n)
            survLength(tempNum) = survLength(n)
            survTime(tempNum) = survTime(n)
            collCount(tempNum) = collCount(n)
            pIndex(tempNum, :) = pIndex(n, :)
        end if
    end do
    pNum = tempNum

    call MPI_SENDRECV(nESend, 1, inttype, neighbor(2), 22, &
                      nERecv, 1, inttype, neighbor(1), 22, comm, status, ierr)
    call MPI_SENDRECV(nWSend, 1, inttype, neighbor(1), 23, &
                      nWRecv, 1, inttype, neighbor(2), 23, comm, status, ierr)
    allocate (exchESend(11*nESend))
    allocate (exchERecv(11*nERecv))
    allocate (exchESendi(nESend))
    allocate (exchERecvi(nERecv))
    exchESend(1:nESend) = xpE(1:nESend)
    exchESend(1*nESend + 1:2*nESend) = ypE(1:nESend)
    exchESend(2*nESend + 1:3*nESend) = zpE(1:nESend)
    exchESend(3*nESend + 1:4*nESend) = upE(1:nESend)
    exchESend(4*nESend + 1:5*nESend) = vpE(1:nESend)
    exchESend(5*nESend + 1:6*nESend) = wpE(1:nESend)
    exchESend(6*nESend + 1:7*nESend) = dpE(1:nESend)
    exchESend(7*nESend + 1:8*nESend) = hpE(1:nESend)
    exchESend(8*nESend + 1:9*nESend) = mhE(1:nESend)
    exchESend(9*nESend + 1:10*nESend) = slE(1:nESend)
    exchESend(10*nESend + 1:11*nESend) = stE(1:nESend)
    exchESendi(1:nESend) = ccE(1:nESend)
    call MPI_SENDRECV(exchESend, nESend*11, realtype, neighbor(2), 24, &
                      exchERecv, nERecv*11, realtype, neighbor(1), 24, &
                      comm, status, ierr)
    call MPI_SENDRECV(exchESendi, nESend, inttype, neighbor(2), 25, &
                      exchERecvi, nERecv, inttype, neighbor(1), 25, &
                      comm, status, ierr)
    if (nERecv > 0) then
        xp(pNum + 1:pNum + nERecv) = exchERecv(1:nERecv)
        yp(pNum + 1:pNum + nERecv) = exchERecv(1*nERecv + 1:2*nERecv)
        zp(pNum + 1:pNum + nERecv) = exchERecv(2*nERecv + 1:3*nERecv)
        up(pNum + 1:pNum + nERecv) = exchERecv(3*nERecv + 1:4*nERecv)
        vp(pNum + 1:pNum + nERecv) = exchERecv(4*nERecv + 1:5*nERecv)
        wp(pNum + 1:pNum + nERecv) = exchERecv(5*nERecv + 1:6*nERecv)
        dp(pNum + 1:pNum + nERecv) = exchERecv(6*nERecv + 1:7*nERecv)
        hp(pNum + 1:pNum + nERecv) = exchERecv(7*nERecv + 1:8*nERecv)
        maxhp(pNum + 1:pNum + nERecv) = exchERecv(8*nERecv + 1:9*nERecv)
        survLength(pNum + 1:pNum + nERecv) = exchERecv(9*nERecv + 1:10*nERecv)
        survTime(pNum + 1:pNum + nERecv) = exchERecv(10*nERecv + 1:11*nERecv)
        collCount(pNum + 1:pNum + nERecv) = exchERecvi(1:nERecv)
        pNum = pNum + nERecv
    end if
    deallocate (exchESend)
    deallocate (exchERecv)
    deallocate (exchESendi)
    deallocate (exchERecvi)

    allocate (exchWSend(11*nWSend))
    allocate (exchWRecv(11*nWRecv))
    allocate (exchWSendi(nWSend))
    allocate (exchWRecvi(nWRecv))
    exchWSend(1:nWSend) = xpW(1:nWSend)
    exchWSend(1*nWSend + 1:2*nWSend) = ypW(1:nWSend)
    exchWSend(2*nWSend + 1:3*nWSend) = zpW(1:nWSend)
    exchWSend(3*nWSend + 1:4*nWSend) = upW(1:nWSend)
    exchWSend(4*nWSend + 1:5*nWSend) = vpW(1:nWSend)
    exchWSend(5*nWSend + 1:6*nWSend) = wpW(1:nWSend)
    exchWSend(6*nWSend + 1:7*nWSend) = dpW(1:nWSend)
    exchWSend(7*nWSend + 1:8*nWSend) = hpW(1:nWSend)
    exchWSend(8*nWSend + 1:9*nWSend) = mhW(1:nWSend)
    exchWSend(9*nWSend + 1:10*nWSend) = slW(1:nWSend)
    exchWSend(10*nWSend + 1:11*nWSend) = stW(1:nWSend)
    exchWSendi(1:nWSend) = ccW(1:nWSend)
    call MPI_SENDRECV(exchWSend, nWSend*11, realtype, neighbor(1), 26, &
                      exchWRecv, nWRecv*11, realtype, neighbor(2), 26, &
                      comm, status, ierr)
    call MPI_SENDRECV(exchWSendi, nWSend, inttype, neighbor(1), 27, &
                      exchWRecvi, nWRecv, inttype, neighbor(2), 27, &
                      comm, status, ierr)
    if (nWRecv > 0) then
        xp(pNum + 1:pNum + nWRecv) = exchWRecv(1:nWRecv)
        yp(pNum + 1:pNum + nWRecv) = exchWRecv(1*nWRecv + 1:2*nWRecv)
        zp(pNum + 1:pNum + nWRecv) = exchWRecv(2*nWRecv + 1:3*nWRecv)
        up(pNum + 1:pNum + nWRecv) = exchWRecv(3*nWRecv + 1:4*nWRecv)
        vp(pNum + 1:pNum + nWRecv) = exchWRecv(4*nWRecv + 1:5*nWRecv)
        wp(pNum + 1:pNum + nWRecv) = exchWRecv(5*nWRecv + 1:6*nWRecv)
        dp(pNum + 1:pNum + nWRecv) = exchWRecv(6*nWRecv + 1:7*nWRecv)
        hp(pNum + 1:pNum + nWRecv) = exchWRecv(7*nWRecv + 1:8*nWRecv)
        maxhp(pNum + 1:pNum + nWRecv) = exchWRecv(8*nWRecv + 1:9*nWRecv)
        survLength(pNum + 1:pNum + nWRecv) = exchWRecv(9*nWRecv + 1:10*nWRecv)
        survTime(pNum + 1:pNum + nWRecv) = exchWRecv(10*nWRecv + 1:11*nWRecv)
        collCount(pNum + 1:pNum + nWRecv) = exchWRecvi(1:nWRecv)
        pNum = pNum + nWRecv
    end if
    deallocate (exchWSend)
    deallocate (exchWRecv)
    deallocate (exchWSendi)
    deallocate (exchWRecvi)

    do n = 1, pNum
        if (yp(n) >= yMax) then
            yp(n) = yp(n) - yMax
        else if (yp(n) < 0.0) then
            yp(n) = yp(n) + yMax
        end if
        if (zp(n) > zMax) then
            zp(n) = zMax
            wp(n) = -abs(wp(n))
        end if
    end do
end subroutine reallocateParticle

subroutine calculateMidAirColl
    use public_val
    use vector_cal
    implicit none

    integer :: globalN1, globalN2
    integer :: i, j
    integer :: localN1, localN2, localNMax
    integer, dimension(mxNode, my) :: pNumInGrid
    integer, dimension(mxNode, my, pNumInGridMax) :: globalN
    real(kind=dbPc) :: distance12
    real(kind=dbPc) :: contactDist
    real(kind=dbPc) :: rv12N, rv12T
    real(kind=dbPc) :: eta1, eta2
    real(kind=dbPc) :: alpha1, beta1
    real(kind=dbPc) :: alpha2, beta2
    real(kind=dbPc) :: resTMidAir
    real(kind=dbPc) :: d1, d2
    real(kind=dbPc), dimension(3) :: pLoc1, pLoc2
    real(kind=dbPc), dimension(3) :: pVol1, pVol2
    real(kind=dbPc), dimension(3) :: rv12
    real(kind=dbPc), dimension(3) :: nVec

    pNumInGrid = 0
    do globalN1 = 1, pNum
        i = pIndex(globalN1, 1)
        j = pIndex(globalN1, 2)
        pNumInGrid(i, j) = pNumInGrid(i, j) + 1
        localN1 = pNumInGrid(i, j)
        globalN(i, j, localN1) = globalN1
    end do
    do j = 1, my
        do i = 1, mxNode
            localNMax = pNumInGrid(i, j)
            if (localNMax < 2) cycle
            do localN1 = 1, localNMax - 1
                globalN1 = globalN(i, j, localN1)
                do localN2 = localN1 + 1, localNMax
                    globalN2 = globalN(i, j, localN2)
                    pLoc1(1) = xp(globalN1)
                    pLoc1(2) = yp(globalN1)
                    pLoc1(3) = zp(globalN1)
                    pLoc2(1) = xp(globalN2)
                    pLoc2(2) = yp(globalN2)
                    pLoc2(3) = zp(globalN2)
                    distance12 = dist_p(pLoc1, pLoc2)
                    d1 = dp(globalN1)
                    d2 = dp(globalN2)
                    contactDist = 0.5*(d1 + d2)
                    if (distance12 > contactDist .or. distance12 <= 0.0) cycle
                    pVol1(1) = up(globalN1)
                    pVol1(2) = vp(globalN1)
                    pVol1(3) = wp(globalN1)
                    pVol2(1) = up(globalN2)
                    pVol2(2) = vp(globalN2)
                    pVol2(3) = wp(globalN2)
                    ! nVec = n (1->2)
                    nVec = unit_vec(pLoc2 - pLoc1)
                    ! rv12=v1-v2
                    rv12 = pVol1 - pVol2
                    ! rv12N=n.v12
                    rv12N = dot_prod(nVec, rv12)
                    if (rv12N <= 0.0) cycle
                    rv12T = sqrt(norm_2(rv12)**2 - rv12N**2)
                    eta1 = d1**3/d2**3
                    eta2 = 1.0/eta1
                    alpha1 = (1.0 + resN)/(1.0 + eta1)
                    alpha2 = (1.0 + resN)/(1.0 + eta2)
                    resTMidAir = 1.0 - 0.4*(1.0 + resN)/(2.0/7.0)*rv12N/rv12T
                    resTMidAir = max(0.0, resTMidAir)
                    beta1 = (2.0/7.0)*(1.0 - resTMidAir)/(1.0 + eta1)
                    beta2 = (2.0/7.0)*(1.0 - resTMidAir)/(1.0 + eta2)
                    up(globalN1) = pVol1(1) - alpha1*rv12N*nVec(1) - beta1*(rv12(1) - rv12N*nVec(1))
                    vp(globalN1) = pVol1(2) - alpha1*rv12N*nVec(2) - beta1*(rv12(2) - rv12N*nVec(2))
                    wp(globalN1) = pVol1(3) - alpha1*rv12N*nVec(3) - beta1*(rv12(3) - rv12N*nVec(3))
                    up(globalN2) = pVol2(1) + alpha2*rv12N*nVec(1) + beta2*(rv12(1) - rv12N*nVec(1))
                    vp(globalN2) = pVol2(2) + alpha2*rv12N*nVec(2) + beta2*(rv12(2) - rv12N*nVec(2))
                    wp(globalN2) = pVol2(3) + alpha2*rv12N*nVec(3) + beta2*(rv12(3) - rv12N*nVec(3))
                    collCount(globalN1) = collCount(globalN1) + 1
                    collCount(globalN2) = collCount(globalN2) + 1
                    xp(globalN2) = xp(globalN2) + (contactDist - distance12)*nVec(1)
                    yp(globalN2) = yp(globalN2) + (contactDist - distance12)*nVec(2)
                    zp(globalN2) = zp(globalN2) + (contactDist - distance12)*nVec(3)
                    !exit
                end do
            end do
        end do
    end do
end subroutine calculateMidAirColl

subroutine addGhostData
    use public_val
    implicit none
    include "mpif.h"

    integer :: i, j, k, jk
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    real(kind=dbPc), dimension(my) :: zChangeESend, zChangeERecv
    real(kind=dbPc), dimension(my) :: zChangeWSend, zChangeWRecv
    real(kind=dbPc), dimension(my) :: zChangeEESend, zChangeEERecv
    real(kind=dbPc), dimension(my*npdf) :: binChangeESend, binChangeERecv
    real(kind=dbPc), dimension(my*npdf) :: binChangeWSend, binChangeWRecv
    real(kind=dbPc), dimension(my*npdf) :: binChangeEESend, binChangeEERecv

    ! because the value of ghost cell has changed
    ! need to add ghost value back to real domain before exchange
    ! zChange, binChange add back
    ! x=mxNode+1 add to x=3: send to 2 and receive from 1
    ! x=mxNode add to x=2: send to 2 and receive from 1
    ! x=1 add to x=mxNode-1: send to 1 and receive from 2
    do j = 1, my
        zChangeESend(j) = zChange(mxNode, j)
        zChangeEESend(j) = zChange(mxNode + 1, j)
        zChangeWSend(j) = zChange(1, j)
        do k = 1, npdf
            jk = k + (j - 1)*npdf
            binChangeESend(jk) = binChange(mxNode, j, k)
            binChangeEESend(jk) = binChange(mxNode + 1, j, k)
            binChangeWSend(jk) = binChange(1, j, k)
        end do
    end do
    call MPI_SENDRECV(zChangeESend, my, realtype, neighbor(2), 14, &
                      zChangeERecv, my, realtype, neighbor(1), 14, comm, status, ierr)
    call MPI_SENDRECV(binChangeESend, my*npdf, realtype, neighbor(2), 15, &
                      binChangeERecv, my*npdf, realtype, neighbor(1), 15, comm, status, ierr)
    call MPI_SENDRECV(zChangeEESend, my, realtype, neighbor(2), 12, &
                      zChangeEERecv, my, realtype, neighbor(1), 12, comm, status, ierr)
    call MPI_SENDRECV(binChangeEESend, my*npdf, realtype, neighbor(2), 13, &
                      binChangeEERecv, my*npdf, realtype, neighbor(1), 13, comm, status, ierr)
    call MPI_SENDRECV(zChangeWSend, my, realtype, neighbor(1), 16, &
                      zChangeWRecv, my, realtype, neighbor(2), 16, comm, status, ierr)
    call MPI_SENDRECV(binChangeWSend, my*npdf, realtype, neighbor(1), 17, &
                      binChangeWRecv, my*npdf, realtype, neighbor(2), 17, comm, status, ierr)
    do j = 1, my
        zChange(2, j) = zChange(2, j) + zChangeERecv(j)
        zChange(3, j) = zChange(3, j) + zChangeEERecv(j)
        zChange(mxNode - 1, j) = zChange(mxNode - 1, j) + zChangeWRecv(j)
        do k = 1, npdf
            jk = k + (j - 1)*npdf
            binChange(2, j, k) = binChange(2, j, k) + binChangeERecv(jk)
            binChange(3, j, k) = binChange(3, j, k) + binChangeEERecv(jk)
            binChange(mxNode - 1, j, k) = binChange(mxNode - 1, j, k) + binChangeWRecv(jk)
        end do
    end do
    ! y=1 add to y=my-1, y=my add to y=2
    do i = 1, mxNode + 1
        zChange(i, 2) = zChange(i, 2) + zChange(i, my)
        zChange(i, my - 1) = zChange(i, my - 1) + zChange(i, 1)
        do k = 1, npdf
            binChange(i, 2, k) = binChange(i, 2, k) + binChange(i, my, k)
            binChange(i, my - 1, k) = binChange(i, my - 1, k) + binChange(i, 1, k)
        end do
    end do
end subroutine addGhostData

subroutine updateSurfGrid
    use public_val
    implicit none

    integer :: i, j, n
    real(kind=dbPc) :: currentBlock, patchBlock
    real(kind=dbPc), dimension(npdf) :: patchBin
    real(kind=dbPc), dimension(mxNode + 1, my, npdf) :: bin

    call addGhostData
    zsf = zsf + zChange/(xDiff*yDiff*por)
    bin = hist*initBlock
    bin = bin + binChange
    do j = 1, my
        do i = 1, mxNode + 1
            if (zsf(i, j) < 0. .or. zsf(i, j) > zMax) then
                print *, 'Error: zsf reach the lower/upper boundary', iter
                stop
            end if
        end do
    end do

    if (whichDiameterDist /= 1) then
        do j = 2, my - 1
            do i = 2, mxNode - 1
                currentBlock = sum(bin(i, j, :))
                patchBlock = initBlock - currentBlock
                if (patchBlock < 0.0) then
                    do n = 1, npdf
                        patchBin(n) = patchBlock*hist(i, j, n)
                        bin(i, j, n) = bin(i, j, n) + patchBin(n)
                    end do
                else
                    do n = 1, npdf
                        patchBin(n) = patchBlock*initDiameterDist(n)
                        bin(i, j, n) = bin(i, j, n) + patchBin(n)
                    end do
                end if
            end do
        end do
        do j = 2, my - 1
            do i = 2, mxNode - 1
                dsf(i, j) = 0.0
                currentBlock = sum(bin(i, j, :))
                do n = 1, npdf
                    hist(i, j, n) = bin(i, j, n)/currentBlock
                    dsf(i, j) = dsf(i, j) + hist(i, j, n)*(binStart + (n - 0.5)*binWidth)
                end do
            end do
        end do
    else
        do n = 1, npdf
            do j = 2, my - 1
                do i = 2, mxNode - 1
                    hist(i, j, n) = initDiameterDist(n)
                end do
            end do
        end do
        dsf = dpa
    end if
    call edgeExch1(zsf)
    call edgeExch1(dsf)
    call edgeExch2(hist)
end subroutine updateSurfGrid

subroutine edgeExch1(data1)
    use public_val
    implicit none
    include "mpif.h"
    ! public
    real(kind=dbPc), dimension(mxNode + 1, my) :: data1
    ! local
    integer :: i, ierr
    integer :: status(MPI_STATUS_SIZE)
    !
    ! send to 2 and receive from 1
    call MPI_SENDRECV(data1(mxNode - 1, 1), 1, edgeType1, neighbor(2), 18, &
                      data1(1, 1), 1, edgeType1, neighbor(1), 18, comm, status, ierr)
    ! send to 1 and receive from 2
    call MPI_SENDRECV(data1(2, 1), 1, edgeType1, neighbor(1), 20, &
                      data1(mxNode, 1), 1, edgeType1, neighbor(2), 20, comm, status, ierr)
    call MPI_SENDRECV(data1(3, 1), 1, edgeType1, neighbor(1), 28, &
                      data1(mxNode + 1, 1), 1, edgeType1, neighbor(2), 28, comm, status, ierr)
    do i = 1, mxNode + 1
        data1(i, my) = data1(i, 2)
        data1(i, 1) = data1(i, my - 1)
    end do
end subroutine edgeExch1

subroutine edgeExch2(data2)
    use public_val
    implicit none
    include "mpif.h"
    ! public
    real(kind=dbPc), dimension(mxNode + 1, my, npdf) :: data2
    ! local
    integer :: i, k, ierr
    integer :: status(MPI_STATUS_SIZE)
    !
    ! send to 2 and receive from 1
    call MPI_SENDRECV(data2(mxNode - 1, 1, 1), 1, edgeType2, neighbor(2), 19, &
                      data2(1, 1, 1), 1, edgeType2, neighbor(1), 19, comm, status, ierr)
    ! send to 1 and receive from 2
    call MPI_SENDRECV(data2(2, 1, 1), 1, edgeType2, neighbor(1), 21, &
                      data2(mxNode, 1, 1), 1, edgeType2, neighbor(2), 21, comm, status, ierr)
    call MPI_SENDRECV(data2(3, 1, 1), 1, edgeType2, neighbor(1), 29, &
                      data2(mxNode + 1, 1, 1), 1, edgeType2, neighbor(2), 29, comm, status, ierr)
    do k = 1, npdf
        do i = 1, mxNode + 1
            data2(i, my, k) = data2(i, 2, k)
            data2(i, 1, k) = data2(i, my - 1, k)
        end do
    end do
end subroutine edgeExch2

subroutine updateFieldGrid
    use public_val
    implicit none
    integer :: i, j, k
    real(kind=dbPc) :: bedElevation

    do i = 1, mxNode
        do j = 1, my
            bedElevation = zsf(i, j)
            do k = 1, mz
                zwReal(i, j, k) = bedElevation + zw(k)
                zReal(i, j, k) = bedElevation + z(k)
            end do
            zwReal(i, j, mz + 1) = bedElevation + zMax
        end do
    end do
end subroutine updateFieldGrid

subroutine calculateFluidField
    use public_val
    implicit none
    include "mpif.h"

    integer :: k
    integer :: ierr
    real(kind=dbPc) :: tau_t
    real(kind=dbPc) :: uT, uB
    real(kind=dbPc) :: dudz, mixl
    real(kind=dbPc) :: dzP, dzT, dzB
    real(kind=dbPc) :: wtt, wtb, wbb, wbt
    real(kind=dbPc) :: nut, nutot
    real(kind=dbPc), dimension(mz) :: phi_pNode

    do k = 1, mz
        phi_pNode(k) = sum(pfrac(2:mxNode - 1, 2:my - 1, k))/(mxNode - 2)/(my - 2)
    end do
    call MPI_ALLREDUCE(phi_pNode, phi_p, mz, realtype, MPI_SUM, comm, ierr)
    call MPI_ALLREDUCE(F_pNode, F_p, mz, realtype, MPI_SUM, comm, ierr)
    tau_t = rho*uStar**2
    tau_p(mz) = F_p(mz)/area
    tau_f(mz) = tau_t - tau_p(mz)
    do k = mz - 1, 1, -1
        tau_p(k) = tau_p(k + 1) + F_p(k)/area
        tau_f(k) = tau_t - tau_p(k)
    end do

    dzP = zDiff(1)
    dzT = 0.5*(zDiff(1) + zDiff(2))
    wtt = 0.5*zDiff(1)/dzT
    wtb = 0.5*zDiff(2)/dzT
    uT = u(2)*wtt + u(1)*wtb
    uB = 0.0
    mixl = kapa*z(1)*(1.0 - exp(-1.0/26.0*z(1)*uStar/nu))
    dudz = (uT - uB)/dzP
    nut = mixl**2*abs(dudz)
    nutot = nu + nut
    uT = tau_f(1)/(rho*nutot)*dzP + uB
    u(1) = 0.5*(uT + uB)
    do k = 2, mz - 1
        dzP = zDiff(k)
        dzT = 0.5*(zDiff(k) + zDiff(k + 1))
        dzB = 0.5*(zDiff(k) + zDiff(k - 1))
        wtt = 0.5*zDiff(k)/dzT
        wtb = 0.5*zDiff(k + 1)/dzT
        wbb = 0.5*zDiff(k)/dzB
        wbt = 0.5*zDiff(k - 1)/dzB
        uT = u(k + 1)*wtt + u(k)*wtb
        uB = u(k - 1)*wbb + u(k)*wbt
        mixl = kapa*z(k)*(1.0 - exp(-1.0/26.0*z(k)*uStar/nu))
        dudz = (uT - uB)/dzP
        nut = mixl**2*abs(dudz)
        nutot = nu + nut
        uT = tau_f(k)/(rho*nutot)*dzP + uB
        u(k) = 0.5*(uT + uB)
    end do
    dzP = zDiff(mz)
    dzB = 0.5*(zDiff(mz) + zDiff(mz - 1))
    wbb = 0.5*zDiff(mz)/dzB
    wbt = 0.5*zDiff(mz - 1)/dzB
    uT = u(mz)
    uB = u(mz - 1)*wbb + u(mz)*wbt
    mixl = kapa*z(mz)*(1.0 - exp(-1.0/26.0*z(mz)*uStar/nu))
    dudz = (uT - uB)/dzP
    nut = mixl**2*abs(dudz)
    nutot = nu + nut
    uT = tau_f(mz)/(rho*nutot)*dzP + uB
    u(mz) = 0.5*(uT + uB)
end subroutine calculateFluidField

subroutine outputAll
    use public_val
    use gather_xyz
    implicit none
    include "mpif.h"
    character(len=300) :: filename, line
    character(len=3) :: str
    integer :: pNumTot
    integer :: serialNum
    integer :: amode
    integer :: fh
    integer :: ierr
    integer :: n
    integer :: status(MPI_STATUS_SIZE)
    integer :: i, j, k
    real(kind=dbPc) :: uMag
    real(kind=dbPc) :: xfluxTot, zfluxTot, hopLength
    real(kind=dbPc), dimension(mx) :: xsfTot, xTot
    real(kind=dbPc), dimension(mxNode) :: xsfTemp
    real(kind=dbPc), dimension(mz) :: xfluxPfTot, zfluxPfTot
    real(kind=dbPc), dimension(mx, my) :: zsfTot, dsfTot
    real(kind=dbPc), dimension(mxNode, my) :: zsfTemp, dsfTemp
    real(kind=dbPc), dimension(mx, my, mz) :: pfracTot, zRealTot

    if (mod(iter, intervalMonitor) == 0) then
        call MPI_ALLREDUCE(pNum, pNumTot, 1, inttype, MPI_SUM, comm, ierr)
        if (myID == 0) then
            open (unit=10, position='append', file='./Particle/ParticleNum.plt')
            write (10, *) time, pNumTot, iter
            close (10)
        end if
    end if

    serialNum = iter/intervalCreateFile
    write (filename, '("./Particle/ParticleData_", I0, ".plt")') serialNum
    write (str, '(I3)') serialNum
    if (mod(iter, intervalCreateFile) == 0) then
        if (myID == 0) then
            open (unit=11, file=trim(adjustl(filename)))
            write (11, "(A100)") 'variables = "x", "y", "z", "u", "v", "w", "u_mag", "h", "h_max", "d", "L", "t", "n_coll"'
            close (11)

            open (unit=13, file='./Field/FieldData_'//trim(adjustl(str))//'.plt')
            write (13, *) 'variables = "x", "y", "z", "u", "phi_p"'
            close (13)

            open (unit=14, file='./Surface/SurfaceData_'//trim(adjustl(str))//'.plt')
            write (14, *) 'variables = "x", "y", "z", "d"'
            close (14)
        end if
    end if

    if (mod(iter, intervalParticle) == 0) then
        if (myID == 0) then
            open (unit=11, position='append', file=trim(adjustl(filename)))
            write (11, *) 'zone', ' T = "', time, '"'
            close (11)
        end if
        amode = MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_APPEND
        call MPI_FILE_OPEN(comm, trim(adjustl(filename)), amode, MPI_INFO_NULL, fh, ierr)
        do n = 1, pNum
            uMag = sqrt(up(n)**2 + vp(n)**2 + wp(n)**2)
            line = ''
            write (line, "(12E15.4, I)") xp(n), yp(n), zp(n), up(n), vp(n), wp(n), uMag, hp(n), maxhp(n), &
                dp(n), survLength(n), survTime(n), collCount(n)
            line = trim(adjustl(line))//char(10)
            call MPI_FILE_WRITE_SHARED(fh, trim(line), len(trim(line)), MPI_CHARACTER, status, ierr)
        end do
        call MPI_FILE_CLOSE(fh, ierr)
    end if

    if (mod(iter, intervalProfile) == 0) then
        if (myID == 0) then
            open (unit=12, position='append', file='./Field/Profile.plt', action='write')
            write (12, *) 'zone', ' T = "', time, '"'
            do k = 1, mz
                write (12, "(6E15.4)") z(k), u(k), tau_p(k), tau_f(k), F_p(k), phi_p(k)
            end do
            close (12)
        end if
    end if

    if (mod(iter, intervalField) == 0) then
        call gatherx(comm, mxNode, mx, x, xTot)
        call gatherxyz(comm, mxNode, mx, my, mz, pfrac, pfracTot)
        call gatherxyz(comm, mxNode, mx, my, mz, zReal, zRealTot)
        if (myID == 0) then
            open (unit=13, position='append', file='./Field/FieldData_'//trim(adjustl(str))//'.plt', action='write')
            write (13, *) 'zone', ' T = "', time, '"'
            write (13, *) 'i=', mx - 2, ' j=', my - 2, ' k=', mz, ' datapacking=point'
            do k = 1, mz
                do j = 2, my - 1
                    do i = 2, mx - 1
                        write (13, "(5E15.4)") xTot(i), y(j), zRealTot(i, j, k), u(k), pfracTot(i, j, k)
                    end do
                end do
            end do
            close (13)
        end if
    end if

    if (mod(iter, intervalSurface) == 0) then
        xsfTemp(1:mxNode) = xsf(1:mxNode)
        zsfTemp(1:mxNode, 1:my) = zsf(1:mxNode, 1:my)
        dsfTemp(1:mxNode, 1:my) = dsf(1:mxNode, 1:my)
        call gatherx(comm, mxNode, mx, xsfTemp, xsfTot)
        call gatherxy(comm, mxNode, mx, my, zsfTemp, zsfTot)
        call gatherxy(comm, mxNode, mx, my, dsfTemp, dsfTot)
        if (myID == 0) then
            open (unit=14, position='append', file='./Surface/SurfaceData_'//trim(adjustl(str))//'.plt', action='write')
            write (14, *) 'zone', ' T = "', time, '"'
            write (14, *) 'i=', mx - 2, ' j=', my - 2, ' datapacking=point'
            do j = 2, my - 1
                do i = 2, mx - 1
                    write (14, "(4E15.4)") xsfTot(i), ysf(j), zsfTot(i, j), dsfTot(i, j)
                end do
            end do
            close (14)
        end if
    end if

    if (mod(iter, intervalStatistics) == 0) then
        call MPI_ALLREDUCE(xflux, xfluxTot, 1, realtype, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(zflux, zfluxTot, 1, realtype, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(xfluxPf, xfluxPfTot, mz, realtype, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(zfluxPf, zfluxPfTot, mz, realtype, MPI_SUM, comm, ierr)
        if (myID == 0) then
            xfluxTot = xfluxTot/nNodes
            zfluxTot = zfluxTot/nNodes
            xfluxPfTot = xfluxPfTot/nNodes
            zfluxPfTot = zfluxPfTot/nNodes
            if (zfluxTot /= 0.0) then
                hopLength = xfluxTot/zfluxTot
            else
                hopLength = 0.0
            end if
            open (unit=15, position='append', file='./Statistics/vsTime.plt', action='write')
            write (15, "(4E15.7)") time, xfluxTot, zfluxTot, hopLength
            close (15)

            open (unit=16, position='append', file='./Statistics/vsHeight.plt', action='write')
            write (16, *) 'zone', ' T = "', time, '"'
            do k = 1, mz
                write (16, "(5E15.7)") z(k), xfluxPfTot(k), zfluxPfTot(k), phi_p(k), phi_p(k)*u(k)*rhoP
            end do
            close (16)
        end if
    end if

end subroutine outputAll
