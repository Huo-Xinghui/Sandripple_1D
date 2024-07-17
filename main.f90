module public_parameter
    ! constants
    implicit none
    integer, parameter :: dbPc = selected_real_kind(15, 307)
    real(kind=dbPc), parameter :: pi = 3.14159265358979323846

    ! computational domain

    real(kind=dbPc), parameter :: xMax = 1.0 ! x size
    real(kind=dbPc), parameter :: yMax = 0.2 ! y size
    real(kind=dbPc), parameter :: zMax = 0.3 ! z size
    integer, parameter :: nx = 500 ! x grid num
    integer, parameter :: ny = 100 ! y grid num
    integer, parameter :: nz = 150 ! z grid num
    integer, parameter :: nzUni = 100 ! z grid number above which zDiff becomes uniform
    real(kind=dbPc), parameter :: xDiff = xMax/nx
    real(kind=dbPc), parameter :: yDiff = yMax/nx
    integer, parameter :: nNodes = 5 ! num of subdomain
    integer, parameter :: mx = nx + 2 ! x grid num +2
    integer, parameter :: my = ny + 2 ! y grid num +2

    ! time

    real(kind=dbPc), parameter :: dt = 1.0e-4 ! time step
    real(kind=dbPc), parameter :: EndTime = 600.0 ! the time program ends

    ! fluid

    real(kind=dbPc), parameter :: uStar = 0.5 ! fractional velocity
    real(kind=dbPc), parameter :: rho = 1.263 ! fluid density
    real(kind=dbPc), parameter :: nu = 1.51e-5 ! kinetic viscosity
    real(kind=dbPc), parameter :: kapa = 0.4 ! von Kaman's constant
    real(kind=dbPc), parameter :: zDiffMin = nu/uStar ! smallest z grid size

    ! particle

    logical, parameter :: ifCalculateParticle = .true.
    logical, parameter :: ifMidairCollision = .true.
    ! particle diameter: whichDiameterDist = 0: normal distribution, 1: uniform diameter, 2: Bernoulli distribution
    ! whichDiameterDist=0: npdf must >= 3, mu=dpa, sigma=dpStddDev, range:mu-3*sigma ~ mu+3*sigma
    ! whichDiameterDist=1: npdf must = 1, d=dpa
    ! whichDiameterDist=2: npdf must = 2, p1=prob1, p2=1-prob1, d1=dpa-dpStddDev, d2=dpa+dpStddDev
    integer, parameter :: whichDiameterDist = 1
    integer, parameter :: npdf = 3 ! bin num of particle distribution
    integer, parameter :: pNumInit = 100 ! initial particle num
    real(kind=dbPc), parameter :: dpa = 2.5e-4 ! average particle diameter
    real(kind=dbPc), parameter :: dpStddDev = 2.0e-4 ! particle diameter standard deviation
    real(kind=dbPc), parameter :: prob1 = 0.5 ! probability one of Bernoulli distribution
    ! The range and width of each bin
    real(kind=dbPc), parameter :: binWidth = 6.0*dpStddDev/npdf
    real(kind=dbPc), parameter :: binStart = dpa - 3.0*dpStddDev
    real(kind=dbPc), parameter :: binEnd = dpa + 3.0*dpStddDev
    real(kind=dbPc), parameter :: repostAngle = 30.0 ! repose angle
    real(kind=dbPc), parameter :: resCoeffN = 0.9 ! normal restitution coefficient
    real(kind=dbPc), parameter :: resCoeffT = 0.0 ! tangential restitution coefficient
    real(kind=dbPc), parameter :: els1 = 0.9 ! normal restitution coefficient (mid-air collision)
    real(kind=dbPc), parameter :: fric1 = 0.0 ! tangential restitution coefficient (mid-air collision)
    real(kind=dbPc), parameter :: bedCellTknessInit = dpa*5.0 ! thickness of the thin surface layer on particle bed
    real(kind=dbPc), parameter :: rhoP = 2650.0 ! particle density
    real(kind=dbPc), parameter :: nkl = 1.0 ! one particle stands for x particles
    real(kind=dbPc), parameter :: por = 0.6 ! bedform porosity
    integer, parameter :: maxEjectNum = 10000 ! max eject particle num in one time step

    ! bed surface

    logical, parameter :: ifPreformedSurface = .false.
    real(kind=dbPc), parameter :: initSurfElevation = 0.05 ! initial average bed height
    real(kind=dbPc), parameter :: initAmp = 8.0*dpa ! amplitude of prerippled surface
    real(kind=dbPc), parameter :: initOmg = 4.0*pi ! wave number of prerippled surface
    real(kind=dbPc), parameter :: wavl = 2.0*pi/initOmg ! wavelength of prerippled surface
    real(kind=dbPc), parameter :: z0 = dpa/30.0 ! initial roughness height

    ! method

    ! boundary condition of particles: 0 periodic, 1 IO, 2 wall
    integer, parameter :: ikbx = 0 ! x direction
    integer, parameter :: ikby = 0 ! y direction
    integer, parameter :: ikbz = 1 ! uppper surface (never be 0)
    integer, parameter :: irsf = 0 ! surface irsf=0: erodable, irsf=1: rigid
    ! output per x steps
    integer, parameter :: nnf = 1e5 ! field
    integer, parameter :: nns = 1e4 ! tau_f & tau_p
    integer, parameter :: nnc = 1e4 ! num of moving particles
    integer, parameter :: nnkl = 1e5 ! particle information
    integer, parameter :: nnsf = 1e4 ! surface, surfaced
    integer, parameter :: nnfx = 1e4 ! sand flux
    ! the initial step
    integer, parameter :: sstart = 1 ! the initial step of surface calculation
    integer, parameter :: pistart = 1 ! the initial step of particle info output
    ! file
    integer, parameter :: nnfi = 1e6 ! iter num contained in a file

end module public_parameter

module parallel_operations
    use public_parameter
    implicit none
    private

    type parallelType
        integer :: ID
        integer :: nx, mx
        integer :: i1, in
        integer :: i0, im
        integer, dimension(2) :: neighbor
        real(kind=dbPc) :: sx, ex
    end type parallelType

    type(parallelType) :: currentNode
    integer :: comm, ierr, coords, nbrLeft, nbrRight, MPI_I, MPI_D
    integer :: MPI_SG_TYPE, MPI_G_TYPE

    public :: initiateParallel, createMpiStructure, freeMpiStructure
    public :: comm, currentNode, MPI_I, MPI_D, MPI_SG_TYPE, MPI_G_TYPE

contains

    subroutine initiateParallel(comm3d)
        implicit none
        integer, intent(in) :: comm3d

        ! create MPI Cartesian topology
        comm = comm3d
        call MPI_COMM_RANK(comm, currentNode%ID, ierr)
        call MPI_CART_GET(comm, 1, nNodes, .true., coords, ierr)
        ! find neighbors
        !
        !       |           |
        !       |           |
        !     -----------------
        !       |           |
        !       |           |
        !     1 |   myid    | 2
        !       |           |
        !       |           |
        !     -----------------
        !       |           |
        !       |           |
        !
        call MPI_CART_SHIFT(comm, 0, 1, nbrLeft, nbrRight, ierr)
        currentNode%neighbor(1) = nbrLeft
        currentNode%neighbor(2) = nbrRight
        call MPI_BARRIER(comm, ierr)

        !             sx         ex
        ! 1  |  2 3 4  |  5 6 7  |  8 9 10  |  11
        !   sx         ex       sx          ex
        ! for the 3rd node shown above:
        ! currentNode%nx = 3 (8, 9, 10)
        ! currentNode%mx = 5 (7, 8, 9, 10, 11)
        ! currentNode%i0 = 7
        ! currentNode%i1 = 8
        ! currentNode%in = 10
        ! currentNode%im = 11
        currentNode%nx = nx/nNodes
        currentNode%mx = nx + 2
        currentNode%i1 = currentNode%ID*currentNode%nx + 2
        currentNode%i0 = currentNode%i1 - 1
        currentNode%in = (currentNode%ID + 1)*currentNode%nx + 1
        currentNode%im = currentNode%in + 1
    end subroutine initiateParallel

    ! need to check in the end
    subroutine createMpiStructure(mpiInteger, mpiDouble)
        implicit none
        integer, intent(in) :: mpiInteger, mpiDouble
        integer, allocatable, dimension(:):: blockLen
        integer, allocatable, dimension(:):: disp
        integer, allocatable, dimension(:):: oldType

        allocate (blockLen(7))
        allocate (disp(7))
        allocate (oldType(7))
        MPI_I = mpiInteger
        MPI_D = mpiDouble
        blockLen = [1, 1, 1, 1, 3, npdf, npdf]
        disp = [0, kind(0), 2*kind(0), 2*kind(0) + kind(0.0_dbPc), &
                2*kind(0) + 2*kind(0.0_dbPc), 2*kind(0) + 5*kind(0.0_dbPc), &
                2*kind(0) + 5*kind(0.0_dbPc) + npdf*kind(0.0_dbPc)]
        oldType = [MPI_I, MPI_I, MPI_D, MPI_D, MPI_D, MPI_D, MPI_D]
        call MPI_TYPE_CREATE_STRUCT(7, blockLen, disp, oldType, MPI_SG_TYPE, ierr)
        call MPI_TYPE_COMMIT(MPI_SG_TYPE, ierr)
        deallocate (blockLen)
        deallocate (disp)
        deallocate (oldType)

        allocate (blockLen(3))
        allocate (disp(3))
        allocate (oldType(3))
        blockLen = [1, 3, 3]
        disp = [0, kind(0.0_dbPc), 4*kind(0.0_dbPc)]
        oldType = [MPI_D, MPI_D, MPI_D]
        call MPI_TYPE_CREATE_STRUCT(3, blockLen, disp, oldType, MPI_G_TYPE, ierr)
        call MPI_TYPE_COMMIT(MPI_G_TYPE, ierr)
        deallocate (blockLen)
        deallocate (disp)
        deallocate (oldType)
    end subroutine createMpiStructure

    subroutine freeMpiStructure
        implicit none

        call MPI_TYPE_FREE(MPI_SG_TYPE, ierr)
        call MPI_TYPE_FREE(MPI_G_TYPE, ierr)
    end subroutine freeMpiStructure
!
!    !subroutine parallelExchangeParticle
!    !    implicit none
!    !    integer :: pNumSend, pNumRecv
!
!    !end subroutine parallelExchangeParticle
!
!    !subroutine parallelGatherParticle
!    !    !use particle_operations
!    !    implicit none
!    !    integer, dimension(nNodes) :: count, displacement
!
!    !    !call MPI_BARRIER(comm3d, ierr)
!    !    !call MPI_ALLREDUCE(pNum, pNumTotal, 1, MPI_INTEGER, MPI_SUM, comm3d, ierr)
!    !    !allocate (allParticle(pNumTotal))
!    !    !displs(1) = 0
!    !    !call MPI_GATHER(nnp, 1, inttype, cnt, 1, inttype, 0, comm3d, ierr)
!    !    !do i = 2, dims
!    !    !    displs(i) = displs(i - 1) + cnt(i - 1)
!    !    !end do
!    !    !call MPI_GATHERV(xp, nnp, realtype, txp, cnt, displs, realtype, 0, comm3d, ierr)
!
!    !end subroutine parallelGatherParticle
!
end module parallel_operations

module vector_operations
    use public_parameter
    implicit none
    private
    public :: dotProduct, crossProduct, vectorMagnitude, unitVector, distance2d, distance3d

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

    ! Compute the distance between two 2D points
    pure function distance2d(p1, p2)
        real(kind=dbPc), dimension(2), intent(in) :: p1, p2
        real(kind=dbPc) :: distance2d

        distance2d = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2)
    end function distance2d

    ! Compute the distance between two 3D points
    pure function distance3d(p1, p2)
        real(kind=dbPc), dimension(3), intent(in) :: p1, p2
        real(kind=dbPc) :: distance3d

        distance3d = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
    end function distance3d

end module vector_operations

module math_operations
    use public_parameter
    implicit none
    private
    public :: valObeyCertainPDF, erfinv, generateNormalDistHistogram !normalDistRand, uniformRand, bernoulliRand,

contains

    function valObeyCertainPDF(histogram)
        implicit none

        ! generate a random number obeying a certain PDF
        real(kind=dbPc) :: valObeyCertainPDF
        real(kind=dbPc), dimension(npdf), intent(in) :: histogram

        integer :: i
        real(kind=dbPc) :: rand
        real(kind=dbPc) :: cumulativeProb

        call random_number(rand)
        cumulativeProb = 0.0
        do i = 1, npdf
            cumulativeProb = cumulativeProb + histogram(i)
            if (cumulativeProb >= rand) then
                valObeyCertainPDF = binStart + (i - 0.5)*binWidth
                exit
            end if
        end do
    end function valObeyCertainPDF

    function erfinv(x)
        implicit none

        ! Inverse error function
        real(kind=dbPc) :: erfinv
        real(kind=dbPc), intent(in) :: x

        real(kind=dbPc) :: t, w, y
        t = 2.0/(pi*0.147) + 0.5*log(1.0 - x**2)
        w = sqrt(t**2 - log(1.0 - x**2)/0.147)
        y = t - w
        erfinv = sign(sqrt(y), x)
    end function erfinv

    subroutine generateNormalDistHistogram(histogram)
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram

        ! Initialize the histogram array to zero
        histogram = 0.0

        ! Fill the histogram array
        do i = 1, npdf
            histogram(i) = exp(-0.5*((binStart + (i - 0.5)*binWidth - dpa)/dpStddDev)**2)/(sqrt(2.0*pi)*dpStddDev)
            histogram(i) = histogram(i)*binWidth
        end do
        histogram = histogram/sum(histogram)
    end subroutine generateNormalDistHistogram

end module math_operations

module surface_operations
    use public_parameter
    implicit none
    private

    type surfaceGridType
        integer :: avalanchTo
        integer :: avalanchFrom
        real(kind=dbPc) :: averageDiameter
        real(kind=dbPc) :: zChange
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc), dimension(npdf) :: diameterDistribution
        real(kind=dbPc), dimension(npdf) :: binVolumeChange
    end type surfaceGridType

    type(surfaceGridType), dimension(mx, my) :: surfGrid
    real(kind=dbPc), dimension(npdf) :: initDiameterDist

    public :: generateSurfaceGrid, initiateSurface, determineParticleRollDirection
    public :: surfGrid, initDiameterDist

contains

    subroutine generateSurfaceGrid
        use parallel_operations
        implicit none
        integer :: i, j, ierr

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
        call MPI_BARRIER(comm, ierr)
        call MPI_BCAST(surfGrid, mx*my, MPI_SG_TYPE, 0, comm, ierr)
    end subroutine generateSurfaceGrid

    subroutine initiateSurface
        use math_operations
        implicit none
        integer :: i, j, k

        if (whichDiameterDist == 0) then
            if (npdf >= 3) then
                call generateNormalDistHistogram(initDiameterDist)
            else
                print *, 'Bin number (npdf) must >= 3 for normal distribution'
                stop
            end if
        else if (whichDiameterDist == 1) then
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
                if (ifPreformedSurface) then
                    surfGrid(i, j)%location(3) = initAmp*sin(initOmg*surfGrid(i, j)%location(1)) + initSurfElevation
                end if
                surfGrid(i, j)%averageDiameter = dpa
                do k = 1, npdf
                    surfGrid(i, j)%diameterDistribution(k) = initDiameterDist(k)
                end do
            end do
        end do
    end subroutine initiateSurface

    subroutine determineParticleRollDirection
        implicit none
        integer :: i, j, k
        real(kind=dbPc) :: centerZ, westZ, eastZ, northZ, southZ
        real(kind=dbPc), dimension(4) :: slopes
        real(kind=dbPc) :: tanRepostAngle, maxSlope

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

        tanRepostAngle = tan(repostAngle/180.0*pi)
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
    end subroutine determineParticleRollDirection

end module surface_operations

module field_operations
    use public_parameter
    implicit none
    private

    type gridType
        real(kind=dbPc) :: zDiff
        ! The location of the bottom-south-west node of the grid
        real(kind=dbPc), dimension(3) :: vLocation
        ! The location of the grid center
        real(kind=dbPc), dimension(3) :: cLocation
        ! for vectorGrid, data(1)=x velocity, data(2)=particle shear stress tau_p, data(3)=total shear stress tau_t
        ! for scalarGrid, data(1)= grid volume, data(2)=particle volume fraction, data(3)=forcing term
        real(kind=dbPc), dimension(3) :: data
    end type gridType

    ! ekalhxh
    type profileType
        real(kind=dbPc) :: zDiff
        real(kind=dbPc), dimension(:), allocatable :: data
    end type profileType

    type(gridType), dimension(mx, my, nz + 1) :: vectorGrid
    type(gridType), dimension(mx, my, nz) :: scalarGrid
    real(kind=dbPc) :: zDiffMax, refineRatio

    public :: generateGrid, initiateField
    public :: vectorGrid, scalarGrid

contains

    subroutine generateGrid
        use parallel_operations
        use surface_operations
        implicit none

        integer :: i, j, k
        integer :: ierr

        do i = 1, mx
            do j = 1, my
                vectorGrid(i, j, 1)%location = surfGrid(i, j)%location
                vectorGrid(i, j, 1)%zDiff = zDiffMin
                zDiffMax = (zMax - surfGrid(i, j)%location(3))/nz
                refineRatio = (zDiffMax/zDiffMin)**(1.0/(nzUni - 1))
                do k = 2, nz
                    vectorGrid(i, j, k)%location(1) = surfGrid(i, j)%location(1)
                    vectorGrid(i, j, k)%location(2) = surfGrid(i, j)%location(2)
                    vectorGrid(i, j, k)%location(3) = vectorGrid(i, j, k - 1)%location(3) + vectorGrid(i, j, k - 1)%zDiff
                    if (k <= nzUni) then
                        vectorGrid(i, j, k)%zDiff = vectorGrid(i, j, k - 1)%zDiff*refineRatio
                    else
                        vectorGrid(i, j, k)%zDiff = zDiffMax
                    end if
                end do
                vectorGrid(i, j, nz + 1)%location(1) = surfGrid(i, j)%location(1)
                vectorGrid(i, j, nz + 1)%location(2) = surfGrid(i, j)%location(2)
                vectorGrid(i, j, nz + 1)%location(3) = zMax
                vectorGrid(i, j, nz)%zDiff = zMax - vectorGrid(i, j, nz)%location(3)
                vectorGrid(i, j, nz + 1)%zDiff = vectorGrid(i, j, nz)%zDiff
            end do
        end do
        currentNode%sx = vectorGrid(currentNode%i1, ny, nz)%location(1)
        currentNode%ex = vectorGrid(currentNode%im, ny, nz)%location(1)

        do k = 1, nz
            do j = 1, my - 1
                do i = 1, mx - 1
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
                scalarGrid(mx, j, k)%location(1) = scalarGrid(mx - 1, j, k)%location(1) + xDiff
                scalarGrid(mx, j, k)%location(2) = scalarGrid(mx - 1, j, k)%location(2)
                scalarGrid(mx, j, k)%location(3) = (vectorGrid(mx, j, k + 1)%location(3) - &
                                                    vectorGrid(mx, j, k)%location(3))/2.0 + &
                                                   vectorGrid(mx, j, k)%location(3)
                scalarGrid(mx, j, k)%zDiff = vectorGrid(mx, j, k)%zDiff
            end do
            do i = 1, mx
                scalarGrid(i, my, k)%location(1) = scalarGrid(i, my - 1, k)%location(1)
                scalarGrid(i, my, k)%location(2) = scalarGrid(i, my - 1, k)%location(2) + yDiff
                scalarGrid(i, my, k)%location(3) = (vectorGrid(i, my, k + 1)%location(3) - &
                                                    vectorGrid(i, my, k)%location(3))/2.0 + &
                                                   vectorGrid(i, my, k)%location(3)
                scalarGrid(i, my, k)%zDiff = vectorGrid(i, my, k)%zDiff
            end do
        end do
        call MPI_BARRIER(comm, ierr)
        call MPI_BCAST(vectorGrid, mx*my*(nz + 1), MPI_G_TYPE, 0, comm, ierr)
        call MPI_BCAST(scalarGrid, mx*my*nz, MPI_G_TYPE, 0, comm, ierr)
    end subroutine generateGrid

    subroutine initiateField
        implicit none
        integer :: i, j, k
        real(kind=dbPc) :: velocityAltitude

        do i = 1, mx
            do j = 1, my
                do k = 1, nz + 1
                    velocityAltitude = scalarGrid(i, j, k)%location(3)
                    vectorGrid(i, j, k)%data(1) = uStar/kapa*log(velocityAltitude/z0)
                    vectorGrid(i, j, k)%data(2) = 0.0
                    vectorGrid(i, j, k)%data(3) = rho*uStar**2
                end do
            end do
        end do
        do i = 1, mx
            do j = 1, my
                do k = 1, nz
                    scalarGrid(i, j, k)%data(1) = xDiff*yDiff*scalarGrid(i, j, k)%zDiff
                    scalarGrid(i, j, k)%data(2) = 0.0
                end do
            end do
        end do
    end subroutine initiateField

end module field_operations

module particle_operations
    use public_parameter
    implicit none
    private

    type particleType
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc), dimension(3) :: velocity
        real(kind=dbPc) :: diameter
        real(kind=dbPc) :: altitude
        integer, dimension(3) :: indices
    end type particleType

    type(particleType), allocatable, dimension(:) :: particle
    type(particleType), allocatable, dimension(:) :: allParticle
    integer pNum, pNumTotal

    public :: particleInitiation, calculateParticleCollisions, calculateParticleMovement
    public :: particle, allParticle, pNum, pNumTotal

contains

    subroutine particleInitiation
        use parallel_operations
        use math_operations
        use surface_operations
        use field_operations
        implicit none

        integer :: n
        real(kind=dbPc) :: rand1, rand2, rand3
        type(particleType) :: currentParticle

        pNum = pNumInit
        allocate (particle(pNum))
        do n = 1, pNum
            currentParticle = particle(n)
            call random_number(rand1)
            call random_number(rand2)
            call random_number(rand3)
            currentParticle%location(1) = currentNode%sx + (currentNode%ex - currentNode%sx)*rand1
            currentParticle%location(2) = yMax*rand2
            currentParticle%location(3) = zMax*rand3
            currentParticle%velocity = 0.0
            currentParticle%diameter = valObeyCertainPDF(initDiameterDist)
            call determineParticleIndices(currentParticle)
            particle(n) = currentParticle
        end do

    end subroutine particleInitiation

    subroutine determineParticleIndices(currentP)
        use surface_operations
        use field_operations
        implicit none

        integer :: ip, jp, kp, tempKp
        real(kind=dbPc) :: currentZ
        type(particleType) :: currentP

        ip = floor(mod(currentP%location(1), xMax)/xDiff) + 2
        jp = floor(mod(currentP%location(2), yMax)/yDiff) + 2
        if (ip == 1) then
            ip = mx - 1
        else if (ip <= 0) then
            ip = nx + ip
        end if
        if (jp == 1) then
            jp = my - 1
        else if (jp <= 0) then
            jp = ny + jp
        end if
        currentZ = surfGrid(ip, jp)%location(3)
        kp = 0
        do tempKp = 1, nz
            currentZ = currentZ + scalarGrid(ip, jp, tempKp)%zDiff
            if (currentP%location(3) < currentZ) then
                kp = tempKp
                exit
            end if
        end do
        currentP%indices(1) = ip
        currentP%indices(2) = jp
        currentP%indices(3) = kp
        scalarGrid(ip, jp, kp)%data(2) = scalarGrid(ip, jp, kp)%data(2) + &
                                         pi*currentP%diameter**3/6.0/scalarGrid(ip, jp, kp)%data(1)
    end subroutine determineParticleIndices

    subroutine calculateParticleCollisions
        use math_operations
        use field_operations
        use surface_operations
        use vector_operations
        implicit none

        integer :: n, nadd, n2
        integer :: totalEjectNum
        integer :: i, j, k
        integer :: ip, jp
        integer :: whichTriangle
        integer :: iterateNum
        integer :: ejectNum
        integer :: currentTotalNum
        integer :: whichVertex
        integer :: closestIP, closestJP
        integer :: whichBin
        integer :: changedIP, changedJP
        integer :: localN, globalN2
        integer, dimension(mx, my) :: pNumInGrid
        integer, allocatable, dimension(:, :, :) :: globalN
        real(kind=dbPc) :: estimateAltitude
        real(kind=dbPc) :: eta, alpha, beta, gama, sigma, lambda, mu
        real(kind=dbPc) :: d1, d2, averageD1D2, nonDimD1, nonDimD2
        real(kind=dbPc) :: v1, v2
        real(kind=dbPc) :: theta1, theta2, theta2Min, theta2Max, theta2Mid
        real(kind=dbPc) :: pT2T1, pMax
        real(kind=dbPc) :: rand1, rand2, rand3
        real(kind=dbPc) :: pointX, pointY
        real(kind=dbPc) :: eBar
        real(kind=dbPc) :: m1, m2
        real(kind=dbPc) :: E1, E2, Ed2, Eeff, E2Bar
        real(kind=dbPc) :: tau_s, tau_tw
        real(kind=dbPc) :: minDistance
        real(kind=dbPc) :: ejectVolume, rollVolume
        real(kind=dbPc) :: distance12, contactDistance
        real(kind=dbPc) :: eta1, eta2, alpha1, alpha2, beta1, beta2
        real(kind=dbPc) :: relativeV12Normal, relativeV12Tangent
        real(kind=dbPc) :: resCoeffTMidAir
        real(kind=dbPc), dimension(4) :: adjacentSurfGridZ
        real(kind=dbPc), dimension(2) :: LocalLoc ! particle local location
        real(kind=dbPc), dimension(3) :: vertex1, vertex2, vertex3
        real(kind=dbPc), dimension(3) :: vector12, vector13
        real(kind=dbPc), dimension(3) :: surfaceNormalVector
        real(kind=dbPc), dimension(3) :: particleProjection
        real(kind=dbPc), dimension(3) :: impactVelocity, reboundVelocity
        real(kind=dbPc), dimension(3) :: impactCoordinateX, impactCoordinateY, impactCoordinateZ
        real(kind=dbPc), dimension(3) :: particleVertexDistance
        real(kind=dbPc), dimension(3) :: relativeV12
        real(kind=dbPc), dimension(npdf) :: initialBinVolumeChange
        type(particleType) :: currentParticle, currentParticle2
        type(particleType), dimension(maxEjectNum) :: addParticle
        type(particleType), allocatable, dimension(:) :: tempParticle

        allocate (tempParticle(pNum + maxEjectNum))
        surfGrid%zChange = 0.0
        initialBinVolumeChange = 0.0
        do ip = 1, mx
            do jp = 1, my
                surfGrid(ip, jp)%binVolumeChange = initialBinVolumeChange
            end do
        end do
        currentTotalNum = 0
        totalEjectNum = 0
        do n = 1, pNum
            currentParticle = particle(n)
            ip = currentParticle%indices(1)
            jp = currentParticle%indices(2)
            adjacentSurfGridZ(1) = surfGrid(ip, jp)%location(3)
            adjacentSurfGridZ(2) = surfGrid(ip + 1, jp)%location(3)
            adjacentSurfGridZ(3) = surfGrid(ip, jp + 1)%location(3)
            adjacentSurfGridZ(4) = surfGrid(ip + 1, jp + 1)%location(3)
            estimateAltitude = currentParticle%location(3) - maxval(adjacentSurfGridZ)
            if (estimateAltitude > 0.5*currentParticle%diameter) then
                currentTotalNum = currentTotalNum + 1
                tempParticle(currentTotalNum) = currentParticle
                cycle
            end if

            LocalLoc(1) = currentParticle%location(1) - surfGrid(ip, jp)%location(1)
            LocalLoc(2) = currentParticle%location(2) - surfGrid(ip, jp)%location(2)
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
                vertex1 = surfGrid(ip, jp)%location
                vertex2 = surfGrid(ip + 1, jp + 1)%location
                if (LocalLoc(1)/LocalLoc(2) <= xDiff/yDiff) then
                    whichTriangle = 1
                    vertex3 = surfGrid(ip, jp + 1)%location
                    vector12 = vertex2 - vertex1
                    vector13 = vertex3 - vertex1
                    surfaceNormalVector = unitVector(crossProduct(vector12, vector13))
                else
                    whichTriangle = 2
                    vertex3 = surfGrid(ip + 1, jp)%location
                    vector12 = vertex2 - vertex1
                    vector13 = vertex3 - vertex1
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
                ! Current surfGrid%location is at point 3, it is the local origin
                vertex1 = surfGrid(ip + 1, jp)%location
                vertex2 = surfGrid(ip, jp + 1)%location
                if (LocalLoc(1)/(yDiff - LocalLoc(2)) <= xDiff/yDiff) then
                    whichTriangle = 3
                    vertex3 = surfGrid(ip, jp)%location
                    vector12 = vertex2 - vertex1
                    vector13 = vertex3 - vertex1
                    surfaceNormalVector = unitVector(crossProduct(vector12, vector13))
                else
                    whichTriangle = 4
                    vertex3 = surfGrid(ip + 1, jp + 1)%location
                    vector12 = vertex2 - vertex1
                    vector13 = vertex3 - vertex1
                    surfaceNormalVector = unitVector(crossProduct(vector13, vector12))
                end if
            end if

            particleProjection(1) = currentParticle%location(1)
            particleProjection(2) = currentParticle%location(2)
            particleProjection(3) = vertex3(3) + &
                                    (surfaceNormalVector(1)*(currentParticle%location(1) - vertex3(1)) + &
                                     surfaceNormalVector(2)*(currentParticle%location(2) - vertex3(2)))/ &
                                    (-surfaceNormalVector(3))
            currentParticle%altitude = currentParticle%location(3) - particleProjection(3)

            particleVertexDistance(1) = distance3d(particleProjection, vertex1)
            particleVertexDistance(2) = distance3d(particleProjection, vertex2)
            particleVertexDistance(3) = distance3d(particleProjection, vertex3)
            minDistance = particleVertexDistance(1)
            whichVertex = 1
            do k = 2, 3
                if (particleVertexDistance(k) < minDistance) then
                    minDistance = particleVertexDistance(k)
                    whichVertex = k
                end if
            end do
            select case (whichTriangle)
                !   3 ------- 2
                !     |    /|
                !     |  /  |
                !     |/    |
                !   1 ------- 3
                ! Left triangle: whichTriangle=1
                ! Right triangle: whichTriangle=2
            case (1)
                select case (whichVertex)
                case (1)
                    closestIP = ip
                    closestJP = jp
                case (2)
                    closestIP = ip + 1
                    closestJP = jp + 1
                case (3)
                    closestIP = ip
                    closestJP = jp + 1
                end select
            case (2)
                select case (whichVertex)
                case (1)
                    closestIP = ip
                    closestJP = jp
                case (2)
                    closestIP = ip + 1
                    closestJP = jp + 1
                case (3)
                    closestIP = ip + 1
                    closestJP = jp
                end select
                !   2 ------- 3
                !     |\    |
                !     |  \  |
                !     |    \|
                !   3 ------- 1
                ! Left triangle: whichTriangle=3
                ! Right triangle: whichTriangle=4
            case (3)
                select case (whichVertex)
                case (1)
                    closestIP = ip + 1
                    closestJP = jp
                case (2)
                    closestIP = ip
                    closestJP = jp + 1
                case (3)
                    closestIP = ip
                    closestJP = jp
                end select
            case (4)
                select case (whichVertex)
                case (1)
                    closestIP = ip + 1
                    closestJP = jp
                case (2)
                    closestIP = ip
                    closestJP = jp + 1
                case (3)
                    closestIP = ip + 1
                    closestJP = jp + 1
                end select
            end select
            if (closestIP >= mx) then
                closestIP = closestIP - mx + 2
            end if
            if (closestJP >= my) then
                closestJP = closestJP - my + 2
            end if

            impactCoordinateZ = surfaceNormalVector
            impactVelocity(3) = dotProduct(currentParticle%velocity, impactCoordinateZ)
            if (currentParticle%altitude < 0.5*currentParticle%diameter .and. impactVelocity(3) < 0.0) then
                impactCoordinateY = unitVector(crossProduct(impactCoordinateZ, currentParticle%velocity))
                impactCoordinateX = unitVector(crossProduct(impactCoordinateY, impactCoordinateZ))
                impactVelocity(1) = dotProduct(currentParticle%velocity, impactCoordinateX)
                impactVelocity(2) = dotProduct(currentParticle%velocity, impactCoordinateY)

                ! Splash function of Lammel et al. 2017
                d1 = currentParticle%diameter
                d2 = surfGrid(closestIP, closestJP)%averageDiameter
                averageD1D2 = (d1 + d2)/2.0
                nonDimD1 = d1/averageD1D2
                nonDimD2 = d2/averageD1D2
                eta = resCoeffN*nonDimD1**3/(nonDimD1**3 + resCoeffN*nonDimD2**3)
                alpha = (1.0 + resCoeffN)/(1.0 + eta) - 1.0
                beta = 1.0 - (2.0/7.0)*(1.0 - resCoeffT)/(1.0 + eta)
                theta1 = atan(abs(impactVelocity(3)/impactVelocity(1)))
                eBar = beta - (beta**2 - alpha**2)*nonDimD2*theta1/(2.0*beta)

                gama = 4.0/9.0/nonDimD2*(beta/(alpha + beta))**2
                theta2Min = -theta1
                theta2Max = 2.0*sqrt(theta1/gama) - theta1
                theta2Max = min(theta2Max, pi)
                theta2Mid = (theta2Max + theta2Min)/2.0
                pMax = gama*(theta1 + theta2Mid)/theta1*log(2.0*theta1/gama/(theta1 + theta2Mid)**2)*1.2
                iterateNum = 0
                do
                    iterateNum = iterateNum + 1
                    call random_number(rand1)
                    call random_number(rand2)
                    pointX = (theta2Max - theta2Min)*rand1 + theta2Min
                    pointY = rand2*min(pMax, 1.0)
                    pointY = max(pointY, 0.0)
                    pT2T1 = gama*(theta1 + pointX)/theta1*log(2.0*theta1/gama/(theta1 + pointX)**2)
                    if (pointY <= pT2T1 .or. iterateNum > 10000) then
                        theta2 = pointX
                        exit
                    end if
                end do

                m1 = rhoP*(pi*d1**3)/6.0
                m2 = rhoP*(pi*d2**3)/6.0
                v1 = vectorMagnitude(impactVelocity)
                E1 = (m1*v1**2)/2.0
                Ed2 = m2*9.8*d2
                ! Shao et al. 2000
                tau_s = rho*0.0123*(rhoP/rho*9.8*d2 + 3.0e-4/(rho*d2))
                tau_tw = vectorGrid(closestIP, closestJP, 1)%data(3)
                Eeff = Ed2*(1.0 - tau_tw/tau_s)
                lambda = 2.0*log((1.0 - eBar**2)*E1/Ed2)
                sigma = sqrt(lambda)*log(2.0)
                mu = log((1.0 - eBar**2)*E1) - lambda*log(2.0)
                E2Bar = Ed2*((1.0 - eBar**2)*E1/Ed2)**(1.0 - (2.0 - log(2.0))*log(2.0))
                ejectNum = nint(0.06*((1.0 - eBar**2)*E1/(2.0*E2Bar))*erfc((log(Eeff) - mu)/(sqrt(2.0)*sigma)))

                v2 = v1*eBar
                E2 = (m2*v2**2)/2.0
                ! whether rebound
                if (theta2 <= 0.0 .or. E2 < Ed2 .or. eBar < 0.0) then
                    select case (surfGrid(closestIP, closestJP)%avalanchTo)
                    case (1)
                        changedIP = closestIP + 1
                        changedJP = closestJP
                    case (2)
                        changedIP = closestIP - 1
                        changedJP = closestJP
                    case (3)
                        changedIP = closestIP
                        changedJP = closestJP + 1
                    case (4)
                        changedIP = closestIP
                        changedJP = closestJP - 1
                    case default
                        changedIP = closestIP
                        changedJP = closestJP
                    end select
                    if (changedIP >= mx) then
                        changedIP = changedIP - mx + 2
                    else if (changedIP <= 1) then
                        changedIP = changedIP + mx - 2
                    end if
                    if (changedJP >= my) then
                        changedJP = changedJP - my + 2
                    else if (changedJP <= 1) then
                        changedJP = changedJP + my - 2
                    end if

                    surfGrid(changedIP, changedJP)%zChange = surfGrid(changedIP, changedJP)%zChange &
                                                             + (pi*d1**3)/6.0/xDiff/yDiff

                    if (whichDiameterDist /= 1) then
                        whichBin = floor((currentParticle%diameter - binStart)/binWidth) + 1
                        whichBin = max(whichBin, 1)
                        whichBin = min(whichBin, npdf)
                        surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) = &
                            surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) + (pi*d1**3)/6.0/xDiff/yDiff
                    end if
                else
                    currentTotalNum = currentTotalNum + 1
                    currentParticle%altitude = 0.5*d1
                    currentParticle%location(3) = particleProjection(3) + 0.5*d1
                    reboundVelocity(1) = v2*cos(theta2)
                    reboundVelocity(2) = 0.0
                    reboundVelocity(3) = v2*sin(theta2)
                    currentParticle%velocity(1) = dotProduct(reboundVelocity, impactCoordinateX)
                    currentParticle%velocity(2) = dotProduct(reboundVelocity, impactCoordinateY)
                    currentParticle%velocity(3) = dotProduct(reboundVelocity, impactCoordinateZ)
                    call determineParticleIndices(currentParticle)
                    tempParticle(currentTotalNum) = currentParticle
                end if

                ejectVolume = 0.0
                do nadd = 1, ejectNum
                    call random_number(rand1)
                    ! E2 = random number obeying log normal distribution with mean mu and standard deviation sigma
                    E2 = exp(mu + sigma*sqrt(2.0)*erfinv(2.0*rand1 - 1.0))
                    d2 = valObeyCertainPDF(surfGrid(closestIP, closestJP)%diameterDistribution)
                    m2 = rhoP*(pi*d2**3)/6.0
                    rand3 = -1.0
                    do while (rand3 < 0.0)
                        call random_number(rand1)
                        call random_number(rand2)
                        rand3 = 1.0 - rand1 - rand2
                    end do
                    totalEjectNum = totalEjectNum + 1
                    addParticle(totalEjectNum)%location(1) = rand1*vertex1(1) + rand2*vertex2(1) + rand3*vertex3(1)
                    addParticle(totalEjectNum)%location(2) = rand1*vertex1(2) + rand2*vertex2(2) + rand3*vertex3(2)
                    addParticle(totalEjectNum)%location(3) = rand1*vertex1(3) + rand2*vertex2(3) + rand3*vertex3(3) + 0.5*d2
                    addParticle(totalEjectNum)%altitude = 0.5*d2
                    addParticle(totalEjectNum)%velocity(1) = 0.0
                    addParticle(totalEjectNum)%velocity(2) = 0.0
                    addParticle(totalEjectNum)%velocity(3) = sqrt(2.0*E2/m2)
                    addParticle(totalEjectNum)%diameter = d2
                    addParticle(totalEjectNum)%indices(1) = ip
                    addParticle(totalEjectNum)%indices(2) = jp
                    addParticle(totalEjectNum)%indices(3) = 1

                    ejectVolume = ejectVolume + (pi*d2**3)/6.0

                    if (whichDiameterDist /= 1) then
                        whichBin = floor((d2 - binStart)/binWidth) + 1
                        whichBin = max(whichBin, 1)
                        whichBin = min(whichBin, npdf)
                        surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) = &
                            surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) - (pi*d2**3)/6.0/xDiff/yDiff
                    end if
                end do
                surfGrid(closestIP, closestJP)%zChange = surfGrid(closestIP, closestJP)%zChange &
                                                         - ejectVolume/xDiff/yDiff

                select case (surfGrid(closestIP, closestJP)%avalanchFrom)
                case (1)
                    changedIP = closestIP + 1
                    changedJP = closestJP
                case (2)
                    changedIP = closestIP - 1
                    changedJP = closestJP
                case (3)
                    changedIP = closestIP
                    changedJP = closestJP + 1
                case (4)
                    changedIP = closestIP
                    changedJP = closestJP - 1
                case default
                    changedIP = closestIP
                    changedJP = closestJP
                end select
                if (changedIP >= mx) then
                    changedIP = changedIP - mx + 2
                else if (changedIP <= 1) then
                    changedIP = changedIP + mx - 2
                end if
                if (changedJP >= my) then
                    changedJP = changedJP - my + 2
                else if (changedJP <= 1) then
                    changedJP = changedJP + my - 2
                end if

                rollVolume = 0.0
                do while (rollVolume < ejectVolume .and. ejectNum > 0)
                    d2 = valObeyCertainPDF(surfGrid(changedIP, changedJP)%diameterDistribution)
                    rollVolume = rollVolume + (pi*d2**3)/6.0
                    if (whichDiameterDist /= 1) then
                        whichBin = floor((d2 - binStart)/binWidth) + 1
                        whichBin = max(whichBin, 1)
                        whichBin = min(whichBin, npdf)
                        surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) = &
                            surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) - (pi*d2**3)/6.0/xDiff/yDiff
                        surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) = &
                            surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) + (pi*d2**3)/6.0/xDiff/yDiff
                    end if
                end do
                surfGrid(changedIP, changedJP)%zChange = surfGrid(changedIP, changedJP)%zChange &
                                                         - rollVolume/xDiff/yDiff
                surfGrid(closestIP, closestJP)%zChange = surfGrid(closestIP, closestJP)%zChange &
                                                         + rollVolume/xDiff/yDiff
            else
                currentTotalNum = currentTotalNum + 1
                tempParticle(currentTotalNum) = currentParticle
            end if
        end do
        pNum = currentTotalNum
        deallocate (particle)

        if (ifMidairCollision) then
            allocate (globalN(mx, my, pNum/nx/ny*2))
            pNumInGrid = 0
            do n = 1, pNum
                currentParticle = tempParticle(n)
                i = currentParticle%indices(1)
                j = currentParticle%indices(2)
                pNumInGrid(i, j) = pNumInGrid(i, j) + 1
                localN = pNumInGrid(i, j)
                globalN(i, j, localN) = n
                if (LocalN < 2) cycle
                do n2 = 1, localN - 1
                    globalN2 = globalN(i, j, n2)
                    currentParticle2 = tempParticle(globalN2)
                    vector12 = currentParticle2%location - currentParticle%location
                    distance12 = vectorMagnitude(vector12)
                    contactDistance = 0.5*(currentParticle%diameter + currentParticle2%diameter)
                    if (distance12 >= contactDistance) cycle
                    vector12 = unitVector(vector12)
                    relativeV12 = currentParticle%velocity - currentParticle2%velocity
                    relativeV12Normal = dotProduct(vector12, relativeV12)
                    if (relativeV12Normal <= 0.0) cycle
                    relativeV12Tangent = sqrt(vectorMagnitude(relativeV12)**2 - relativeV12Normal**2)
                    eta1 = currentParticle%diameter**3/currentParticle2%diameter**3
                    eta2 = 1.0/eta1
                    alpha1 = (1.0 + resCoeffN)/(1.0 + eta1)
                    alpha2 = (1.0 + resCoeffN)/(1.0 + eta2)
                    resCoeffTMidAir = max(0.0, 1.0 - 0.4*(1.0 + resCoeffN)/(2.0/7.0)*relativeV12Normal/relativeV12Tangent)
                    beta1 = (2.0/7.0)*(1.0 - resCoeffTMidAir)/(1.0 + eta1)
                    beta2 = (2.0/7.0)*(1.0 - resCoeffTMidAir)/(1.0 + eta2)
                    tempParticle(globalN2)%velocity = currentParticle2%velocity + alpha2*relativeV12Normal*vector12 &
                                                      + beta2*(vector12 - relativeV12Normal*vector12)
                    tempParticle(globalN2)%location = currentParticle2%location + 0.5*(contactDistance - distance12)*vector12
                    exit
                end do
                tempParticle(n)%velocity = currentParticle%velocity - alpha1*relativeV12Normal*vector12 &
                                           - beta1*(vector12 - relativeV12Normal*vector12)
                tempParticle(n)%location = currentParticle%location - 0.5*(contactDistance - distance12)*vector12
            end do
            deallocate (globalN)
        end if

        if (totalEjectNum > 0) then
            pNum = pNum + totalEjectNum
            tempParticle(currentTotalNum + 1:pNum) = addParticle(1:totalEjectNum)
        else if (totalEjectNum > maxEjectNum) then
            print *, "The eject number in one time step", totalEjectNum, ">", maxEjectNum
            stop
        end if
        allocate (particle(pNum))
        particle = tempParticle(1:pNum)
        deallocate (tempParticle)

    end subroutine calculateParticleCollisions

    subroutine calculateParticleMovement
        use field_operations
        implicit none
        integer :: n, ip, jp, kp, i
        real(kind=dbPc) :: dp, mp
        real(kind=dbPc) :: C_d, Re_p
        real(kind=dbPc), dimension(3) :: up, uf, relativeU
        real(kind=dbPc), dimension(3) :: bulkForce, totalForce
        real(kind=dbPc), dimension(3) :: u1, u2, u3, u4
        real(kind=dbPc), dimension(3) :: a1, a2, a3, a4
        type(particleType) :: currentParticle

        scalarGrid%data(3) = 0.0
        do n = 1, pNum
            currentParticle = particle(n)
            up = currentParticle%velocity
            dp = currentParticle%diameter
            ip = currentParticle%indices(1)
            jp = currentParticle%indices(2)
            kp = currentParticle%indices(3)
            uf(1) = vectorGrid(ip, jp, kp)%data(1)
            uf(2) = 0.0
            uf(3) = 0.0
            mp = rhoP*(pi*dp**3)/6.0
            bulkForce(1:2) = 0.0
            bulkForce(3) = -9.8*(1.0 - rho/rhoP)*mp

            ! use 4th order Runge-Kutta method to solve the ODE
            u1 = up
            totalForce = dragForce() + bulkForce
            a1 = totalForce/mp
            up = u1 + 0.5*dt*a1
            u2 = up
            totalForce = dragForce() + bulkForce
            a2 = totalForce/mp
            up = u1 + 0.5*dt*a2
            u3 = up
            totalForce = dragForce() + bulkForce
            a3 = totalForce/mp
            up = u1 + dt*a3
            u4 = up
            totalForce = dragForce() + bulkForce
            a4 = totalForce/mp
            currentParticle%location = currentParticle%location + (u1 + 2.0*u2 + 2.0*u3 + u4)/6.0*dt
            currentParticle%velocity = currentParticle%velocity + (a1 + 2.0*a2 + 2.0*a3 + a4)/6.0*dt
            call determineParticleIndices(currentParticle)
            particle(n) = currentParticle
            scalarGrid(ip, jp, kp)%data(3) = scalarGrid(ip, jp, kp)%data(3) + a1(1)
        end do

    contains

        function dragForce() result(Fd)
            implicit none
            real(kind=dbPc), dimension(3) :: Fd

            relativeU = uf - up
            Re_p = abs(relativeU(1))*dp/nu
            C_d = (sqrt(0.5) + sqrt(24.0/Re_p))**2
            do i = 1, 3
                Fd(i) = -pi/8.0*rho*dp**2*C_d*abs(relativeU(i))*relativeU(i)
            end do
        end function dragForce

    end subroutine calculateParticleMovement

    subroutine reallocateParticle
        use parallel_operations
        implicit none
        integer :: n, ip
        type(particleType) :: currentParticle

        do n = 1, pNum
            currentParticle = particle(n)
        end do

    end subroutine reallocateParticle

end module particle_operations

module output_operations
    use public_parameter
    implicit none

    private

    public :: generateOutPutFile

contains

    subroutine generateOutPutFile
        implicit none
        character(len=32) :: bashCmd
        external :: system

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
    use parallel_operations
    use surface_operations
    use field_operations
    use particle_operations
    use output_operations
    implicit none
    include "mpif.h"

    integer :: ierr, comm3d, mpiInteger, mpiDouble
    integer :: iteration
    real(kind=dbPc) :: time

    ! find indices of subdomain and check that dimensions of arrays are sufficient
    if (mod(nx, nNodes) /= 0) then
        print *, 'nx cannot diveded by nNodes'
        stop
    end if
    call MPI_INIT(ierr)
    call MPI_CART_CREATE(MPI_COMM_WORLD, 1, nNodes, .true., .true., comm3d, ierr)
    call initiateParallel(comm3d)
    call random_seed()
    mpiInteger = MPI_INTEGER
    mpiDouble = MPI_DOUBLE
    call createMpiStructure(mpiInteger, mpiDouble)
    ! generate surfGrid and initial bed
    call generateSurfaceGrid
    ! initiate surface
    call initiateSurface
    ! generate grid
    call generateGrid
    ! initiate fluid field
    call initiateField
    ! initiate particle
    call particleInitiation
    ! create output file
    call generateOutPutFile
    iteration = 1
    time = 0.0

    ! start iteration loop
    do while (time < Endtime)
        time = time + dt
        iteration = iteration + 1
        if (ifCalculateParticle) then
            call determineParticleRollDirection
            call calculateParticleCollisions
            call calculateParticleMovement
            call reallocateParticle
            !        if (iteration < sstart) then
            !            Dkz = 0.0
            !            do i = 1, mkxNode
            !                do j = 1, mky
            !                    bedCellTkness(i, j) = bedCellTknessInit
            !                    if (irsf == 0) then
            !                        do k = 1, npdf
            !                            bedPDist(i, j, k) = prob(k)
            !                        end do
            !                        bedPD(i, j) = dpa
            !                    else
            !                        bedPDist(i, j, 2) = 0.5*(0.5*sin(initOmg*kx(i)) + 0.5)
            !                        bedPDist(i, j, 1) = 1.0 - bedPDist(i, j, 2)
            !                        bedPD(i, j) = bedPDist(i, j, 1)*(dpa - dpStddDev) + bedPDist(i, j, 2)*(dpa + dpStddDev)
            !                    end if
            !                end do
            !            end do
            !        end if
        end if
        !    ! calculate fluid field
        !    call fluidField
        !    ! generate boundary key point
        !    call imgd
        !    phirho = 1.0
        !    ! output result
        !    call output
        !    ! time advance
    end do
    call freeMpiStructure

    call MPI_FINALIZE(ierr)
end program main

