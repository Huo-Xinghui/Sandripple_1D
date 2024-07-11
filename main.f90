module public_parameter
    ! constants
    implicit none
    integer, parameter :: dbPc = selected_real_kind(14, 307)
    real(kind=dbPc), parameter :: pi = 3.14159265358979323846

    ! computational domain

    real(kind=dbPc), parameter :: xMax = 1.0 ! x size
    real(kind=dbPc), parameter :: yMax = 0.2 ! y size
    real(kind=dbPc), parameter :: zMax = 0.5 ! z size
    integer, parameter :: nx = 500 ! x grid num
    integer, parameter :: ny = 100 ! y grid num
    integer, parameter :: nz = 250 ! z grid num
    integer, parameter :: nzUni = 100 ! z grid number above which zDiff becomes uniform
    real(kind=dbPc), parameter :: xDiff = xMax/nx
    real(kind=dbPc), parameter :: yDiff = yMax/nx
    integer, parameter :: nNodes = 5 ! num of subdomain

    ! time

    real(kind=dbPc), parameter :: dt = 1.0e-4 ! time step
    real(kind=dbPc), parameter :: tla = 600.0 ! time last

    ! fluid

    real(kind=dbPc), parameter :: uStar = 0.5 ! fractional velocity
    real(kind=dbPc), parameter :: rho = 1.263 ! fluid density
    real(kind=dbPc), parameter :: nu = 1.51e-5 ! kinetic viscosity
    real(kind=dbPc), parameter :: kapa = 0.4 ! von Kaman's constant
    real(kind=dbPc), parameter :: zDiffMin = nu/uStar ! smallest z grid size

    ! particle

    integer, parameter :: ipar = 1 ! calculating particles: ipar = 0: no, 1: yes
    ! particle diameter: ipd = 0: normal distribution, 1: uniform diameter, 2: Bernoulli distribution
    ! ipd=0: npdf must >= 3, mu=dpa, sigma=dpStddDev, range:mu-3*sigma ~ mu+3*sigma
    ! ipd=1: npdf must = 1, d=dpa
    ! ipd=2: npdf must = 2, p1=prob1, p2=1-prob1, d1=dpa-dpStddDev, d2=dpa+dpStddDev
    integer, parameter :: ipd = 1
    integer, parameter :: npdf = 3 ! bin num of particle distribution
    integer, parameter :: pNumInit = 100 ! initial particle num
    real(kind=dbPc), parameter :: dpa = 2.5e-4 ! average particle diameter
    real(kind=dbPc), parameter :: dpStddDev = 2.0e-4 ! particle diameter standard deviation
    real(kind=dbPc), parameter :: prob1 = 0.5 ! probability one of Bernoulli distribution
    real(kind=dbPc), parameter :: repostAngle = 30.0 ! repose angle
    real(kind=dbPc), parameter :: resCoeffN = 0.9 ! normal restitution coefficient
    real(kind=dbPc), parameter :: resCoeffT = 0.0 ! tangential restitution coefficient
    real(kind=dbPc), parameter :: els1 = 0.9 ! normal restitution coefficient (mid-air collision)
    real(kind=dbPc), parameter :: fric1 = 0.0 ! tangential restitution coefficient (mid-air collision)
    real(kind=dbPc), parameter :: bedCellTknessInit = dpa*5.0 ! thickness of the thin surface layer on particle bed
    real(kind=dbPc), parameter :: rhoP = 2650.0 ! particle density
    real(kind=dbPc), parameter :: nkl = 1.0 ! one particle stands for x particles
    real(kind=dbPc), parameter :: por = 0.6 ! bedform porosity
    integer, parameter :: pNumMax = 1000000 ! max particle num in one subdomain
    integer, parameter :: nspmax = 10000 ! max eject particle num in one time step

    ! bed surface

    integer, parameter :: isf = 0 ! initial surface isf=0: flat, isf=1: prerippled
    real(kind=dbPc), parameter :: initSurfElevation = 0.05 ! initial average bed height
    real(kind=dbPc), parameter :: initAmp = 8.0*dpa ! amplitude of prerippled surface
    real(kind=dbPc), parameter :: initOmg = 4.0*pi ! wave number of prerippled surface
    real(kind=dbPc), parameter :: wavl = 2.0*pi/initOmg ! wavelength of prerippled surface

    ! method

    integer, parameter :: isp = 0 ! splash function: isp=0:lammel, isp=1:kok.
    ! boundary condition of particles: 0 periodic, 1 IO, 2 wall
    integer, parameter :: ikbx = 0 ! x direction
    integer, parameter :: ikby = 0 ! y direction
    integer, parameter :: ikbz = 1 ! uppper surface (never be 0)
    integer, parameter :: icol = 1 ! mid-air collision: icol=1: on, icol=0: off
    integer, parameter :: irsf = 0 ! surface irsf=0: erodable, irsf=1: rigid
    ! output per x steps
    integer, parameter :: nnf = 1e5 ! field
    integer, parameter :: nns = 1e4 ! tau_a & tau_p
    integer, parameter :: nnc = 1e4 ! num of moving particles
    integer, parameter :: nnkl = 1e5 ! particle information
    integer, parameter :: nnsf = 1e4 ! surface, surfaced
    integer, parameter :: nnfx = 1e4 ! sand flux
    ! the initial step
    integer, parameter :: sstart = 1 ! the initial step of surface calculation
    integer, parameter :: pistart = 1 ! the initial step of particle info output
    ! file
    integer, parameter :: nnfi = 1e6 ! iter num contained in a file
    ! others
    integer, parameter :: mx = nx + 2 ! x grid num +2
    integer, parameter :: my = ny + 2 ! y grid num +2
end module public_parameter

module vector_operations
    use public_parameter
    implicit none
    private
    public :: dotProduct, crossProduct, vectorMagnitude, unitVector

    interface dotProduct
        module procedure dotProductNd
    end interface

    interface crossProduct
        module procedure crossProduct3d  ! Only for 3D vectors
    end interface

    interface vectorMagnitude
        module procedure vectorMagnitudeNd
    end interface

    interface unitVector
        module procedure unitVectorNd
    end interface

contains

    ! Compute the dot product of two arbitrary-dimensional vectors
    pure function dotProductNd(a, b)
        real(kind=dbPc), dimension(:), intent(in) :: a, b
        real(kind=dbPc) :: dotProductNd
        integer :: i
        dotProductNd = 0.0
        do i = 1, size(a)
            dotProductNd = dotProductNd + a(i)*b(i)
        end do
    end function dotProductNd

    ! Compute the cross product of two 3D vectors
    pure function crossProduct3d(a, b)
        real(kind=dbPc), dimension(3), intent(in) :: a, b
        real(kind=dbPc), dimension(3) :: crossProduct3d
        crossProduct3d(1) = a(2)*b(3) - a(3)*b(2)
        crossProduct3d(2) = a(3)*b(1) - a(1)*b(3)
        crossProduct3d(3) = a(1)*b(2) - a(2)*b(1)
    end function crossProduct3d

    ! Compute the magnitude of an arbitrary-dimensional vector
    pure function vectorMagnitudeNd(v)
        real(kind=dbPc), dimension(:), intent(in) :: v
        real(kind=dbPc) :: vectorMagnitudeNd
        integer :: i
        vectorMagnitudeNd = 0.0
        do i = 1, size(v)
            vectorMagnitudeNd = vectorMagnitudeNd + v(i)**2
        end do
        vectorMagnitudeNd = sqrt(vectorMagnitudeNd)
    end function vectorMagnitudeNd

    ! Compute the unit vector of an arbitrary-dimensional vector
    pure function unitVectorNd(v)
        real(kind=dbPc), dimension(:), intent(in) :: v
        real(kind=dbPc), dimension(size(v)) :: unitVectorNd
        real(kind=dbPc) :: magnitude
        integer :: i
        magnitude = vectorMagnitudeNd(v)
        if (magnitude /= 0.0) then
            do i = 1, size(v)
                unitVectorNd(i) = v(i)/magnitude
            end do
        else
            unitVectorNd = 0.0  ! Return the zero vector if the input vector has zero magnitude
        end if
    end function unitVectorNd
end module vector_operations

module surface_operations
    use public_parameter
    implicit none
    private

    type surfaceGridType
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc) :: averageDiameter
        real(kind=dbPc), dimension(npdf + 2) :: diameterDistribution
        integer :: avalanchTo, avalanchFrom
    end type surfaceGridType

    type(surfaceGridType), dimension(mx, my) :: surfGrid
    real(kind=dbPc), dimension(npdf + 2) :: initDiameterDist

    public :: generateSurfaceGrid, surfaceInitiation, determineParticleRollDirection
    public :: surfGrid, initDiameterDist

contains

    subroutine generateSurfaceGrid
        implicit none
        integer :: i, j

        ! In fact the location of surfGrid(i, j) here is the location of the south-west node of the grid
        surfGrid(1, :)%location(1) = -xDiff
        surfGrid(:, 1)%location(2) = -yDiff
        surfGrid%location(3) = initSurfElevation

        do j = 2, my
            do i = 2, mx
                surfGrid(i, j)%location(1) = surfGrid(i - 1, j)%location(1) + xDiff
                surfGrid(i, j)%location(2) = surfGrid(i, j - 1)%location(2) + yDiff
            end do
        end do

        ! Justify the location to the center of the surface grid
        surfGrid%location(1) = surfGrid%location(1) + 0.5*xDiff
        surfGrid%location(2) = surfGrid%location(2) + 0.5*yDiff
    end subroutine generateSurfaceGrid

    subroutine surfaceInitiation
        implicit none
        integer :: i, j, k

        if (ipd == 0) then
            if (npdf >= 3) then
                call normalDistHistogram(initDiameterDist)
            else
                print *, 'Bin number (npdf) must >= 3 for normal distribution'
                stop
            end if
        else if (ipd == 1) then
            if (npdf == 1) then
                initDiameterDist(1) = 1.0
                initDiameterDist(2) = 0.0 ! Bin width
                initDiameterDist(3) = dpa ! Bin start
            else
                print *, 'Bin number (npdf) must = 1 for uniform particle diameter'
                stop
            end if
        else
            if (npdf == 2) then
                initDiameterDist(1) = prob1
                initDiameterDist(2) = 1.0 - initDiameterDist(1)
                initDiameterDist(3) = 2.0*dpStddDev ! Bin width
                initDiameterDist(4) = dpa - 2.0*dpStddDev ! Bin start
            else
                print *, 'Bin number (npdf) must = 2 for Bernoulli distribution'
                stop
            end if
        end if
        do j = 1, my
            do i = 1, mx
                if (isf == 1) then
                    surfGrid(i, j)%location(3) = initAmp*sin(initOmg*surfGrid(i, j)%location(1)) + initSurfElevation
                end if
                surfGrid(i, j)%averageDiameter = dpa
                do k = 1, npdf + 2
                    surfGrid(i, j)%diameterDistribution(k) = initDiameterDist(k)
                end do
            end do
        end do
    end subroutine surfaceInitiation

    subroutine normalDistHistogram(histogram)
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf + 2) :: histogram
        real(kind=dbPc) :: binWidth, binStart, binEnd

        ! Calculate the range and width of each bin
        binWidth = 6.0*dpStddDev/npdf
        binStart = dpa - 3.0*dpStddDev
        binEnd = dpa + 3.0*dpStddDev

        ! Initialize the histogram array to zero
        histogram = 0.0

        ! Fill the histogram array
        do i = 1, npdf
            histogram(i) = exp(-0.5*((binStart + (i - 0.5)*binWidth - dpa)/dpStddDev)**2)/(sqrt(2.0*pi)*dpStddDev)
            histogram(i) = histogram(i)*binWidth
        end do
        histogram = histogram/sum(histogram)
        histogram(npdf + 1) = binWidth
        histogram(npdf + 2) = binStart
    end subroutine normalDistHistogram

    subroutine determineParticleRollDirection
        use omp_lib
        implicit none
        integer :: i, j, k, threadID, numThreads
        real(kind=dbPc) :: centerZ, westZ, eastZ, northZ, southZ
        real(kind=dbPc), dimension(4) :: slopes
        real(kind=dbPc) :: tanRepostAngle, maxSlope

        ! avalanchTo:
        !       3 ^
        !         |
        !  2 <----0----> 1
        !         |
        !         v 4

        ! avalanchFrom:
        !       3 |
        !         v
        !  2 ---->0<---- 1
        !         ^
        !         | 4

        tanRepostAngle = tan(repostAngle/180.0*pi)
        ! Initialize OpenMP
        numThreads = omp_get_max_threads()
        call omp_set_num_threads(numThreads)

        !$omp parallel do private(i, j, k, centerZ, westZ, eastZ, northZ, southZ, slopes, maxSlope) &
        !$omp& shared(surfGrid, xDiff, yDiff, repostAngle, pi, nx, ny)
        do j = 2, ny
            do i = 2, nx
                centerZ = surfGrid(i, j)%location(3)
                westZ = surfGrid(i - 1, j)%location(3)
                eastZ = surfGrid(i + 1, j)%location(3)
                northZ = surfGrid(i, j + 1)%location(3)
                southZ = surfGrid(i, j - 1)%location(3)

                slopes(1) = (centerZ - eastZ)/xDiff
                slopes(2) = (centerZ - westZ)/xDiff
                slopes(3) = (centerZ - northZ)/yDiff
                slopes(4) = (centerZ - southZ)/yDiff

                surfGrid(i, j)%avalanchTo = 0
                maxSlope = tanRepostAngle
                do k = 1, 4
                    if (slopes(k) > maxSlope) then
                        maxSlope = slopes(k)
                        surfGrid(i, j)%avalanchTo = k
                    end if
                end do

                slopes(1) = (eastZ - centerZ)/xDiff
                slopes(2) = (westZ - centerZ)/xDiff
                slopes(3) = (northZ - centerZ)/yDiff
                slopes(4) = (southZ - centerZ)/yDiff

                surfGrid(i, j)%avalanchFrom = 0
                maxSlope = tanRepostAngle
                do k = 1, 4
                    if (slopes(k) > maxSlope) then
                        maxSlope = slopes(k)
                        surfGrid(i, j)%avalanchFrom = k
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine determineParticleRollDirection
end module surface_operations

module field_operations
    use public_parameter
    implicit none
    private

    type gridType
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc) :: zDiff
    end type gridType

    type(gridType), dimension(mx + 1, my + 1, nz + 1) :: vectorGrid
    type(gridType), dimension(mx, my, nz) :: scalarGrid
    real(kind=dbPc) :: zDiffMax, refineRatio

    public :: generateGrid, scalarGrid

contains

    subroutine generateGrid
        use surface_operations
        implicit none

        integer :: i, j, k

        do i = 1, mx
            do j = 1, my
                if (i == 1) then
                    vectorGrid(i, j, 1)%location(1) = -xDiff
                else
                    vectorGrid(i, j, 1)%location(1) = vectorGrid(i - 1, j, 1)%location(1) + xDiff
                end if
                if (j == 1) then
                    vectorGrid(i, j, 1)%location(2) = -yDiff
                else
                    vectorGrid(i, j, 1)%location(2) = vectorGrid(i, j - 1, 1)%location(2) + yDiff
                end if
                vectorGrid(i, j, 1)%location(3) = surfGrid(i, j)%location(3)
                vectorGrid(i, j, 1)%zDiff = zDiffMin
                zDiffMax = (zMax - surfGrid(i, j)%location(3))/nz
                refineRatio = (zDiffMax/zDiffMin)**(1.0/(nzUni - 1))
                do k = 2, nz
                    if (i == 1) then
                        vectorGrid(i, j, k)%location(1) = -xDiff
                    else
                        vectorGrid(i, j, k)%location(1) = vectorGrid(i - 1, j, k)%location(1) + xDiff
                    end if
                    if (j == 1) then
                        vectorGrid(i, j, k)%location(2) = -yDiff
                    else
                        vectorGrid(i, j, k)%location(2) = vectorGrid(i, j - 1, k)%location(2) + yDiff
                    end if
                    vectorGrid(i, j, k)%location(3) = vectorGrid(i, j, k - 1)%location(3) + vectorGrid(i, j, k - 1)%zDiff
                    if (k <= nzUni) then
                        vectorGrid(i, j, k)%zDiff = vectorGrid(i, j, k - 1)%zDiff*refineRatio
                    else
                        vectorGrid(i, j, k)%zDiff = zDiffMax
                    end if
                end do
                if (i == 1) then
                    vectorGrid(i, j, nz + 1)%location(1) = -xDiff
                else
                    vectorGrid(i, j, nz + 1)%location(1) = vectorGrid(i - 1, j, nz + 1)%location(1) + xDiff
                end if
                if (j == 1) then
                    vectorGrid(i, j, nz + 1)%location(2) = -yDiff
                else
                    vectorGrid(i, j, nz + 1)%location(2) = vectorGrid(i, j - 1, nz + 1)%location(2) + yDiff
                end if
                vectorGrid(i, j, nz + 1)%location(3) = zMax
                vectorGrid(i, j, nz)%zDiff = zMax - vectorGrid(i, j, nz)%location(3)
                vectorGrid(i, j, nz + 1)%zDiff = vectorGrid(i, j, nz)%zDiff
            end do
        end do
        vectorGrid(mx + 1, :, :) = vectorGrid(mx, :, :)
        vectorGrid(:, my + 1, :) = vectorGrid(:, my, :)

        do k = 1, nz
            do j = 1, my
                do i = 1, mx
                    scalarGrid(i, j, k)%location(1) = (vectorGrid(i + 1, j, k)%location(1) - &
                                                       vectorGrid(i, j, k)%location(1))/2.0 + &
                                                      vectorGrid(i, j, k)%location(1)
                    scalarGrid(i, j, k)%location(2) = (vectorGrid(i, j + 1, k)%location(2) - &
                                                       vectorGrid(i, j, k)%location(2))/2.0 + &
                                                      vectorGrid(i, j, k)%location(2)
                    scalarGrid(i, j, k)%location(3) = (vectorGrid(i, j, k + 1)%location(3) - &
                                                       vectorGrid(i, j, k)%location(3))/2.0 + &
                                                      vectorGrid(i, j, k)%location(3)
                    scalarGrid(i, j, k)%zDiff = vectorGrid(i, j, k)%zDiff
                end do
            end do
        end do
    end subroutine generateGrid
end module field_operations

module particle_operations
    use public_parameter
    implicit none
    private

    type particleList
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc), dimension(3) :: velocity
        real(kind=dbPc) :: diameter
        real(kind=dbPc) :: altitude
        integer, dimension(3) :: indices
        type(particleList), pointer :: next => null()
    end type particleList

    type(particleList), pointer :: pListHead => null(), pListTail => null(), particle => null()
    integer pNum

    public :: particleInitiation

contains

    subroutine particleInitiation
        use surface_operations
        implicit none

        integer :: n
        real(kind=dbPc) :: rand1, rand2, rand3

        pNum = pNumInit
        do n = 1, pNum
            call random_number(rand1)
            call random_number(rand2)
            call random_number(rand3)
            allocate (particle)
            particle%location(1) = xMax*rand1
            particle%location(2) = yMax*rand2
            particle%location(3) = zMax*rand3
            particle%velocity = 0.0
            particle%diameter = valObeyCertainPDF(initDiameterDist)
            call determineParticleIndices
            if (.not. associated(pListHead)) then
                pListHead => particle
                pListTail => particle
            else
                pListTail%next => particle
                pListTail => particle
            end if
        end do
    end subroutine particleInitiation

    subroutine determineParticleIndices
        use surface_operations
        use field_operations
        implicit none

        integer :: ip, jp, kp
        real(kind=dbPc) :: currentZ

        ip = int(mod(particle%location(1), xMax)/xDiff) + 2
        jp = int(mod(particle%location(2), yMax)/yDiff) + 2
        currentZ = surfGrid(ip, jp)%location(3)
        kp = 0
        if (particle%location(3) <= currentZ) then
            kp = -1
        else if (particle%location(3) > zMax) then
            kp = 0
        else
            do kp = 1, nz
                currentZ = currentZ + scalarGrid(ip, jp, kp)%zDiff
                if (particle%location(3) < currentZ) then
                    exit
                end if
            end do
        end if
        particle%indices(1) = ip
        particle%indices(2) = jp
        particle%indices(3) = kp
    end subroutine determineParticleIndices

    function valObeyCertainPDF(histogram)
        implicit none

        real(kind=dbPc) :: valObeyCertainPDF
        real(kind=dbPc), dimension(npdf + 2), intent(in) :: histogram

        integer :: i
        real(kind=dbPc) :: rand
        real(kind=dbPc) :: cumulativeProb

        call random_number(rand)
        cumulativeProb = 0.0
        do i = 1, npdf
            cumulativeProb = cumulativeProb + histogram(i)
            if (cumulativeProb >= rand) then
                valObeyCertainPDF = histogram(npdf + 2) + (i - 0.5)*histogram(npdf + 1)
                exit
            end if
        end do
    end function valObeyCertainPDF

    subroutine calculateImpactSplash
        use surface_operations
        use vector_operations
        implicit none

        integer ip, jp
        integer :: iterateNum
        real(kind=dbPc), dimension(2) :: LocalLoc ! particle local location
        real(kind=dbPc), dimension(3) :: node1, node2, node3
        real(kind=dbPc), dimension(3) :: vector12, vector13
        real(kind=dbPc), dimension(3) :: surfaceNormalVector
        real(kind=dbPc), dimension(3) :: particleProjection
        real(kind=dbPc), dimension(3) :: impactVelocity
        real(kind=dbPc), dimension(3) :: impactCoordinateX, impactCoordinateY, impactCoordinateZ
        real(kind=dbPc) :: eta, alpha, beta, gama, sigma, lambda, mu
        real(kind=dbPc) :: d1, d2, averageD1D2, nonDimD1, nonDimD2
        real(kind=dbPc) :: m1, m2
        real(kind=dbPc) :: v1, v2
        real(kind=dbPc) :: theta1, theta2, theta2Min, theta2Max, theta2Mid
        real(kind=dbPc) :: pT2T1, pMax
        real(kind=dbPc) :: rand1, rand2
        real(kind=dbPc) :: pointX, pointY
        real(kind=dbPc) :: E1, E2
        real(kind=dbPc) :: eBar

        particle => pListHead
        do while (associated(particle))
            ip = particle%indices(1)
            jp = particle%indices(2)
            LocalLoc(1) = particle%location(1) - surfGrid(ip, jp)%location(1)
            LocalLoc(2) = particle%location(2) - surfGrid(ip, jp)%location(2)

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
                node1 = surfGrid(ip, jp)%location
                node2 = surfGrid(ip + 1, jp + 1)%location
                if (LocalLoc(1)/LocalLoc(2) <= xDiff/yDiff) then
                    whichTriangle = 1
                    node3 = surfGrid(ip, jp + 1)%location
                    vector12 = node2 - node1
                    vector13 = node3 - node1
                    surfaceNormalVector = unitVector(crossProduct(vector12, vector13))
                else
                    whichTriangle = 2
                    node3 = surfGrid(ip + 1, jp)%location
                    vector12 = node2 - node1
                    vector13 = node3 - node1
                    surfaceNormalVector = unitVector(crossProduct(vector13, vector12))
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
                ! Current surfGrid%location is at the lower left point 3, it is the local origin
                node1 = surfGrid(ip + 1, jp)%location
                node2 = surfGrid(ip, jp + 1)%location
                if (LocalLoc(1)/(yDiff - LocalLoc(2)) <= xDiff/yDiff) then
                    whichTriangle = 3
                    node3 = surfGrid(ip, jp)%location
                    vector12 = node2 - node1
                    vector13 = node3 - node1
                    surfaceNormalVector = unitVector(crossProduct(vector12, vector13))
                else
                    whichTriangle = 4
                    node3 = surfGrid(ip + 1, jp + 1)%location
                    vector12 = node2 - node1
                    vector13 = node3 - node1
                    surfaceNormalVector = unitVector(crossProduct(vector13, vector12))
                end if
            end if
            particleProjection(1) = particle%location(1)
            particleProjection(2) = particle%location(2)
            particleProjection(3) = node3(3) + &
                                    (surfaceNormalVector(1)*(particle%location(1) - node3(1)) + &
                                     surfaceNormalVector(2)*(particle%location(2) - node3(2)))/ &
                                    (-surfaceNormalVector(3))
            particle%altitude = particle%location(3) - particleProjection(3)
            if (particle%altitude <= 0.0) then
                impactCoordinateZ = surfaceNormalVector
                impactVelocity(3) = dotProduct(particle%velocity, impactCoordinateZ)
                if (impactVelocity(3) < 0.0) then
                    impactCoordinateY = crossProduct(impactCoordinateZ, particle%velocity)
                    impactCoordinateX = crossProduct(impactCoordinateY, impactCoordinateZ)
                    impactVelocity(1) = dotProduct(particle%velocity, impactCoordinateX)
                    impactVelocity(2) = dotProduct(particle%velocity, impactCoordinateY)
                    ! Splash function of Lammel et al. 2017
                    d1 = particle%diameter
                    d2 = surfGrid(ip, jp)%averageDiameter
                    averageD1D2 = (d1 + d2)/2.0
                    nonDimD1 = d1/averageD1D2
                    nonDimD2 = d2/averageD1D2
                    eta = resCoeffN*nonDimD1**3/(nonDimD1**3 + resCoeffN*nonDimD2**3)
                    alpha = (1.0 + resCoeffN)/(1.0 + eta) - 1.0
                    beta = 1.0 - (2.0/7.0)*(1.0 - resCoeffT)/(1.0 + eta)
                    theta1 = atan(abs(impactVelocity(3)/impactVelocity(1)))
                    eBar = beta - (beta**2 - alpha**2)*nonDimD2*theta1/(2.0*beta)
                    v1 = vectorMagnitude(impactVelocity)
                    v2 = v1*eBar

                    gama = 4.0/9.0/nonDimD2*(beta/(alpha + beta))**2
                    theta2Min = -theta1
                    theta2Max = 2.0*sqrt(theta1/gama) - theta1
                    theta2Mid = (theta2Max + theta2Min)/2.0
                    pMax = gama*(theta1 + theta2Mid)/theta1*log(2.0*theta1/gama/(theta1 + theta2Mid)**2)*1.2
                    iterateNum = 0
                    do
                        iterateNum = iterateNum + 1
                        call random_number(rand1)
                        call random_number(rand2)
                        pointX = (theta2Max - theta2Min)*rand1 + theta2Min
                        pointY = rand2*min(pMax, 1)
                        pointY = max(pointY, 0.0)
                        pT2T1 = gama*(theta1 + pointX)/theta1*log(2.0*theta1/gama/(theta1 + pointX)**2)
                        if (pointY <= pT2T1 .or. iterateNum > 10000) then
                            theta2 = pointX
                            exit
                        end if
                    end do

                    norm_vin = norm_2(pVelVec)
                    ee1 = 0.5*mm1*norm_vin**2
                    ! particle rebound
                    angout1 = arebound(alpha, beta, angin1, dd2)
                    norm_vout = pp*norm_vin
                    ee2 = mm1/2.0*norm_vout**2
                    vout(3) = norm_vout*sin(angout1)
                    if (vout(3) < sqrt(2.0*gg3*0.5*(d1 + d2)) .or. pp <= 0.0 .or. pp > 1.0 .or. angout1 <= 0.0) then
                        nne = 0
                        pp = 0.0
                    else
                        nne = 1
                    end if
                end if
            end if
        end do
        point0 = closestPoint2Triangle(point0, point1, point2, point3)
        unitSurfNormal = unit_vec(surfNormal)
        !
        pLocVec(1) = xp(n)
        pLocVec(2) = yp(n)
        pLocVec(3) = zp(n)
        vec1p = pLocVec - point1
        dotP1 = dotProduct(vec12, vec1p)
        dotP2 = dotProduct(vec13, vec1p)
        if (dotP1 <= 0.0 .and. dotP2 <= 0.0) then
            point0 = point1
        else
            vec2p = pLocVec - point2
            dotP3 = dotProduct(vec12, vec2p)
            dotP4 = dotProduct(vec13, vec2p)
            if (dotP3 >= 0.0 .and. dotP4 <= dotP3) then
                point0 = point2
            else
                lam23 = dotP1*dotP4 - dotP3*dotP2
                if (lam23 <= 0.0 .and. dotP1 >= 0.0 .and. dotP3 <= 0.0) then
                end if
            end if
            lambda3 = (vec1(1)*vec3(2) - vec1(2)*vec3(1))/(vec1(1)*vec2(2) - vec1(2)*vec2(1))
            lambda2 = (vec3(1) - vec2(1)*lambda3)/vec1(1)
            lambda1 = 1.0 - lambda2 - lambda3
            point0(1) = pLocVec(1)
            point0(2) = pLocVec(2)
            point0(3) = lambda1*point1(3) + lambda2*point2(3) + lambda3*point3(3)
            !
            fz(n) = pLocVec(3) - point0(3)
            !
            vin(3) = dotProduct(pVelVec, unitSurfNormal)
            ifimpact: if (fz(n) <= 0.0 .and. vin(3) < 0.0) then
                ! information of inject particle
                pVelVec(1) = up(n)
                pVelVec(2) = vp(n)
                pVelVec(3) = wp(n)
                ttvec = crossProduct(surfNormal, pVelVec)
                tnvec = crossProduct(ttvec, surfNormal)
                utnvec = unit_vec(tnvec)
                uttvec = unit_vec(ttvec)
                vin(1) = dotProduct(pVelVec, utnvec)
                vin(2) = dotProduct(pVelVec, uttvec)
                gg(1) = 0.0
                gg(2) = 0.0
                gg(3) = 9.8
                gg3 = 9.8 !abs(dotProduct(gg, unitSurfNormal))
                norm_vin = norm_2(pVelVec)
                ! nearest point to impact point
                d01 = dist_p(point0, point1)
                d02 = dist_p(point0, point2)
                d03 = dist_p(point0, point3)
                select case (whichSurfTri)
                case (1)
                    if (d01 <= d02 .and. d01 <= d03) then
                        ipp = ik
                        jpp = jk
                    else if (d02 <= d01 .and. d02 <= d03) then
                        ipp = ik + 1
                        jpp = jk + 1
                    else if (d03 <= d01 .and. d03 <= d02) then
                        ipp = ik
                        jpp = jk + 1
                    end if
                case (2)
                    if (d01 <= d02 .and. d01 <= d03) then
                        ipp = ik
                        jpp = jk
                    else if (d02 <= d01 .and. d02 <= d03) then
                        ipp = ik + 1
                        jpp = jk + 1
                    else if (d03 <= d01 .and. d03 <= d02) then
                        ipp = ik + 1
                        jpp = jk
                    end if
                case (3)
                    if (d01 <= d02 .and. d01 <= d03) then
                        ipp = ik + 1
                        jpp = jk
                    else if (d02 <= d01 .and. d02 <= d03) then
                        ipp = ik
                        jpp = jk + 1
                    else if (d03 <= d01 .and. d03 <= d02) then
                        ipp = ik
                        jpp = jk
                    end if
                case (4)
                    if (d01 <= d02 .and. d01 <= d03) then
                        ipp = ik + 1
                        jpp = jk
                    else if (d02 <= d01 .and. d02 <= d03) then
                        ipp = ik
                        jpp = jk + 1
                    else if (d03 <= d01 .and. d03 <= d02) then
                        ipp = ik + 1
                        jpp = jk + 1
                    end if
                case default
                    print *, 'whichSurfTri error1'
                    stop
                end select
                if (nne == 0) then
                    ! influence of repose angle
                    select case (rollDirBump(ipp, jpp))
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
                    case default
                        print *, rollDirBump(ipp, jpp), 'rollDirBump error'
                        stop
                    end select
                    if (ipd /= 2) then
                        iii = int((d1 - dpa + dSigma*3.0)/dSigma/6.0*dfloat(npdf)) + 1
                        if (iii > npdf) iii = npdf
                        if (iii <= 0) iii = 1
                    else
                        if (d1 < dpa) then
                            iii = 1
                        else if (d1 > dpa) then
                            iii = 2
                        else
                            print *, 'error on d1=', d1, '/=', dpa, 'or', dpa + dSigma
                            stop
                        end if
                    end if
                    vch = nkl*(pi*d1**3)/6.0/por
                    if (jj >= mky + 1) jj = 3
                    if (ii <= mkxNode) then
                        Dkz(ii, jj) = Dkz(ii, jj) + vch/kArea
                        DbedPDist(ii, jj, iii) = DbedPDist(ii, jj, iii) + vch
                    else
                        eepnch(jj) = eepnch(jj) + vch
                        jjkk = iii + (jj - 1)*npdf
                        eepdfch(jjkk) = eepdfch(jjkk) + vch
                    end if
                else
                    ammu2 = 0.0
                    aSigma2 = 10.0/180.0*pi
                    !angout2 = normal(ammu2, aSigma2)
                    vout(1) = norm_vout*cos(angout1) !*cos(angout2)
                    vout(2) = 0.0 !norm_vout*cos(angout1)*sin(angout2)
                    vout(3) = norm_vout*sin(angout1)
                    vector1 = vout(1)*utnvec
                    vector2 = vout(2)*uttvec
                    vector3 = vout(3)*unitSurfNormal
                    upvec2 = vector1 + vector2 + vector3
                    !upvec2(1) = upvec2(1) + ucreep(ipp, jpp)
                    !upvec2(2) = upvec2(2) + vcreep(ipp, jpp)
                    pNumTemp = pNumTemp + 1
                    xp(pNumTemp) = point0(1) !+ upvec2(1)*dt
                    yp(pNumTemp) = point0(2) !+ upvec2(2)*dt
                    zp(pNumTemp) = point0(3) !+ upvec2(3)*dt
                    up(pNumTemp) = upvec2(1)
                    vp(pNumTemp) = upvec2(2)
                    wp(pNumTemp) = upvec2(3)
                    dp(pNumTemp) = d1
                    fk(pNumTemp) = 0.0 !fk(n) !upvec2(1)*dt
                    fh(pNumTemp) = 0.0
                    fg(pNumTemp) = zp(pNumTemp)
                    ft(pNumTemp) = 0.0
                    fz(pNumTemp) = 0.0 !upvec2(3)*dt
                    wflx = wflx + nkl*mm1/xMax/yMax/dt
                    vpout = vpout + vout
                    norm_vpout = norm_vpout + norm_vout
                    vvpout = vvpout + norm_vout**2
                    mpout = mpout + mm1
                    mupout(ik, jk) = mupout(ik, jk) + mm2*vout(1)
                    mvpout(ik, jk) = mvpout(ik, jk) + mm2*vout(2)
                    npout = npout + 1.0
                end if
                ! particle splash
                if (isp == 0) then ! lammel
                    utaot = sqrt(0.0123*(rhos/rho*9.8*bedPD(ipp, jpp) + 3.0e-4/(rho*bedPD(ipp, jpp))))
                    taot = rho*utaot**2
                    eed2x = mm2*gg3*d2
                    eed2 = eed2x*(1.0 - htao(1)/taot)
                    if (eed2/eed2x <= 0.1) then
                        eed2 = eed2x*0.1
                    end if
                    lambda = 2.0*log((1.0 - pp**2)*ee1/eed2)
                    if (lambda <= 0.0) then
                        ne = 0
                    else
                        sigma = sqrt(lambda)*log(2.0)
                        mmu = log((1.0 - pp**2)*ee1) - lambda*log(2.0)
                        merfc = myerfc((log(eed2) - mmu)/(sqrt(2.0)*sigma))
                        ee2bar = eed2*((1.0 - pp**2)*ee1/eed2)**(1.0 - (2.0 - log(2.0))*log(2.0))
                        ne = int(0.06*((1.0 - pp**2)*ee1/(2.0*ee2bar))*merfc)
                    end if
                else ! kok
                    mm2 = (pi*bedPD(ipp, jpp)**3)/6.0*rhos
                    nee = 0.03*norm_vin/sqrt(9.8*bedPD(ipp, jpp))
                    ne = int(nee)
                end if
                if (ne >= 1) then
                    ! influence of repose angle
                    select case (rollDirSink(ipp, jpp))
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
                    case default
                        print *, rollDirSink(ipp, jpp), 'rollDirSink error'
                        stop
                    end select
                    splp: do kd = 1, ne
                        ammu1 = 60.0/180.0*pi
                        ammu2 = 0.0
                        aSigma1 = 15.0/180.0*pi
                        aSigma2 = 10.0/180.0*pi
                        if (ipd == 0) then
                            ppdf = bedPDist(ipp, jpp, :)
                            dpd = valObeyCertainProbDist(ppdf, dpa, dSigma, npdf, rpdf)
                        else if (ipd == 1) then
                            dpd = dpa
                        else if (ipd == 2) then
                            ppdf = bedPDist(ipp, jpp, 1)
                            nbi = biDist(ppdf)
                            if (nbi == 1) then
                                dpd = dpa - dSigma
                            else
                                dpd = dpa + dSigma
                            end if
                        end if
                        mm2 = (pi*dpd**3)/6.0*rhos
                        if (isp == 0) then ! lammel
                            ee2 = exp(normal(mmu, sigma))
                            norm_vout = sqrt(2.0*ee2/mm2)
                        else ! kok
                            mmu = 0.08*norm_vin
                            lambda = 1.0/mmu
                            norm_vout = expdev(lambda)
                            ee2 = 0.5*mm2*norm_vout**2
                        end if
                        !lambda = 1.0/ammu1
                        !angout1 = expdev(lambda)
                        !angout1 = abs(normal(ammu1, aSigma1)) !Dupont
                        !angout2 = normal(ammu2, aSigma2)
                        vout(1) = 0.0 !norm_vout*cos(angout1) !*cos(angout2)
                        vout(2) = 0.0 !norm_vout*cos(angout1)*sin(angout2)
                        vout(3) = norm_vout !*sin(angout1)
                        vector1 = vout(1)*utnvec
                        vector2 = vout(2)*uttvec
                        vector3 = vout(3)*unitSurfNormal
                        upvec2 = vector1 + vector2 + vector3
                        !upvec2(1) = upvec2(1) + ucreep(ipp, jpp)
                        !upvec2(2) = upvec2(2) + vcreep(ipp, jpp)
                        call random_number(rr1)
                        call random_number(rr2)
                        call random_number(rr3)
                        pNumAdd = pNumAdd + 1
                        tempx(pNumAdd) = point0(1) + upvec2(1)*dt*rr1
                        tempy(pNumAdd) = point0(2) + upvec2(2)*dt*rr2
                        tempz(pNumAdd) = point0(3) + upvec2(3)*dt*rr3
                        tempu(pNumAdd) = 0.0 !upvec2(1)
                        tempv(pNumAdd) = 0.0 !upvec2(2)
                        tempw(pNumAdd) = norm_vout !upvec2(3)
                        tempd(pNumAdd) = dpd
                        tempfk(pNumAdd) = 0.0 !upvec2(1)*dt
                        tempfh(pNumAdd) = 0.0
                        tempfg(pNumAdd) = tempz(pNumAdd)
                        tempft(pNumAdd) = 0.0
                        tempfz(pNumAdd) = 0.0 !upvec2(3)*dt
                        wflx = wflx + nkl*mm2/xMax/yMax/dt
                        vpout = vpout + vout
                        norm_vpout = norm_vpout + norm_vout
                        vvpout = vvpout + norm_vout**2
                        mpout = mpout + mm2
                        mupout(ik, jk) = mupout(ik, jk) + mm2*vout(1)
                        mvpout(ik, jk) = mvpout(ik, jk) + mm2*vout(2)
                        npout = npout + 1.0
                        if (ipd /= 2) then
                            iii = int((dpd - dpa + dSigma*3.0)/dSigma/6.0*dfloat(npdf)) + 1
                            if (iii > npdf) iii = npdf
                            if (iii <= 0) iii = 1
                        else
                            if (dpd < dpa) then
                                iii = 1
                            else if (dpd > dpa) then
                                iii = 2
                            else
                                print *, 'error on dpd=', dpd, '/=', dpa, 'or', dpa + dSigma
                                stop
                            end if
                        end if
                        vch = nkl*(pi*dpd**3)/6.0/por
                        if (jj == mky + 1) jj = 3
                        if (ii <= mkxNode) then
                            Dkz(ii, jj) = Dkz(ii, jj) - vch/kArea
                            DbedPDist(ii, jj, iii) = DbedPDist(ii, jj, iii) - vch
                        else
                            eepnch(jj) = eepnch(jj) - vch
                            jjkk = iii + (jj - 1)*npdf
                            eepdfch(jjkk) = eepdfch(jjkk) - vch
                        end if
                    end do splp
                end if
            else
                pNumTemp = pNumTemp + 1
                xp(pNumTemp) = xp(n)
                yp(pNumTemp) = yp(n)
                zp(pNumTemp) = zp(n)
                up(pNumTemp) = up(n)
                vp(pNumTemp) = vp(n)
                wp(pNumTemp) = wp(n)
                dp(pNumTemp) = dp(n)
                fk(pNumTemp) = fk(n)
                fh(pNumTemp) = fh(n)
                fg(pNumTemp) = fg(n)
                ft(pNumTemp) = ft(n)
                fz(pNumTemp) = fz(n)
            end if ifimpact

            end subroutine calculateImpactSplash
            end module particle_operations

            module output_operations
                use public_parameter
                implicit none

                private

                public :: generateOutPutFile

            contains

                subroutine generateOutPutFile
                    implicit none

                    character(len=32) bashCmd

                    bashCmd = 'mkdir particle_loc'
                    call system(trim(adjustl(bashCmd)))
                    bashCmd = 'mkdir surface'
                    call system(trim(adjustl(bashCmd)))
                    bashCmd = 'mkdir field3D'
                    call system(trim(adjustl(bashCmd)))

                    open (unit=10, file='./field/field0.plt')
                    write (10, *) 'variables = "X", "Y", "Z", "U"'
                    close (10)

                    open (unit=11, file='./surface/surface0.plt')
                    write (11, *) 'variables = "X", "Y", "Z", "D"'
                    close (11)

                    open (unit=12, file='./particle_loc/particle_loc0.plt')
                    write (12, "(A82)") 'variables = "X", "Y", "Z", "U", "V", "W", "D"'
                    close (12)

                    !open (unit=31, file='particle_num.plt')
                    !write (31, *) 'variables = "T", "Num"'
                    !close (31)

                    !open (unit=35, file='average_flux.plt')
                    !write (35, *) 'variables = "T", "uFlux", "wFlux", "salength"'
                    !close (35)

                    !open (unit=36, file='flux_vs_height.plt')
                    !write (36, *) 'variables = "Z", "uFlux", "wFlux"'
                    !close (36)

                    !open (unit=43, file='htao.plt')
                    !write (43, *) 'variables = "Z", "taoa", "taop", "vfrac", "u", "fptx"'
                    !close (43)

                    !open (unit=39, file='vin.plt')
                    !write (39, *) 'variables = "T", "upin", "vpin", "wpin", "norm_vpin"'
                    !close (39)

                    !open (unit=46, file='vout.plt')
                    !write (46, *) 'variables = "T", "upout", "vpout", "wpout", "norm_vpout"'
                    !close (46)

                    !open (unit=44, file='eminout.plt')
                    !write (44, *) 'variables = "T", "vvpin", "vvpout", "mpin", "mpout"'
                    !close (44)

                    !open (unit=45, file='numinout.plt')
                    !write (45, *) 'variables = "T", "numin", "numout"'
                    !close (45)
                end subroutine generateOutPutFile
            end module output_operations

            program main
                use public_parameter
                use field_operations
                use surface_operations
                use output_operations
                implicit none

                integer :: last

                call random_seed()
                ! generate surfGrid and initial bed
                call generateSurfaceGrid
                ! initiate surface
                call surfaceInitiation
                ! generate grid
                call generateGrid
                ! initiate particle
                call particleInitiation
                ! creat output file
                call generateOutPutFile
                last = 1

                ! start iteration loop
                do
                    if (ipar == 1) then
                        call determineParticleRollDirection
                        call pickOutImpactParticles
                        call calculateParticleSplash
                        if (last < sstart) then
                            Dkz = 0.0
                            do i = 1, mkxNode
                                do j = 1, mky
                                    bedCellTkness(i, j) = bedCellTknessInit
                                    if (irsf == 0) then
                                        do k = 1, npdf
                                            bedPDist(i, j, k) = prob(k)
                                        end do
                                        bedPD(i, j) = dpa
                                    else
                                        bedPDist(i, j, 2) = 0.5*(0.5*sin(initOmg*kx(i)) + 0.5)
                                        bedPDist(i, j, 1) = 1.0 - bedPDist(i, j, 2)
                                        bedPD(i, j) = bedPDist(i, j, 1)*(dpa - dpStddDev) + bedPDist(i, j, 2)*(dpa + dpStddDev)
                                    end if
                                end do
                            end do
                        end if
                    end if
                    ! calculate fluid field
                    call fluidField
                    ! generate boundary key point
                    call imgd
                    phirho = 1.0
                    ! output result
                    call output
                    ! time advance
                    time = time + dt
                    last = last + 1
                    if (time > tla) exit
                end do
            end program main

