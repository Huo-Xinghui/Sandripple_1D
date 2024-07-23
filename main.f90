! This is a Fortran 90 code for simulating the sandripple formation in a 1D wind field
! version 2.0
! Author: Xianghui Huo
! Date: 2024-07-24

! *****************************************

! The module contains the public parameters

! *****************************************
module public_parameter
    ! constants
    implicit none
    integer, parameter :: dbPc = selected_real_kind(15, 307)
    real(kind=dbPc), parameter :: pi = 3.14159265358979323846

    ! computational domain

    real(kind=dbPc), parameter :: xMax = 1.0 ! x size
    real(kind=dbPc), parameter :: yMax = 0.1 ! y size
    real(kind=dbPc), parameter :: zMax = 0.3 ! z size
    real(kind=dbPc), parameter :: area = xMax*yMax ! computational area
    integer, parameter :: nx = 500 ! x grid num
    integer, parameter :: ny = 50 ! y grid num
    integer, parameter :: nz = 50 ! z grid num
    integer, parameter :: nzUni = 20 ! z grid number above which zDiff becomes uniform
    real(kind=dbPc), parameter :: xDiff = xMax/nx
    real(kind=dbPc), parameter :: yDiff = yMax/ny
    integer, parameter :: nNodes = 2 ! num of subdomain
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
    integer, parameter :: whichDiameterDist = 0
    integer, parameter :: npdf = 11 ! bin num of particle distribution
    integer, parameter :: pNumInit = 10 ! initial particle num
    integer, parameter :: maxEjectNum = 10000 ! max eject particle num in one time step
    integer, parameter :: maxSendRecvNum = 1000 ! max particle num to send or receive in one time step
    integer, parameter :: maxParticleNum = 10000000 ! max particle num
    real(kind=dbPc), parameter :: dpa = 2.5e-4 ! average particle diameter
    real(kind=dbPc), parameter :: dpStddDev = 8.0e-5 ! particle diameter standard deviation
    real(kind=dbPc), parameter :: prob1 = 0.5 ! probability one of Bernoulli distribution
    ! The range and width of each bin
    real(kind=dbPc), parameter :: binWidth = 6.0*dpStddDev/npdf
    real(kind=dbPc), parameter :: binStart = dpa - 3.0*dpStddDev
    real(kind=dbPc), parameter :: binEnd = dpa + 3.0*dpStddDev
    real(kind=dbPc), parameter :: repostAngle = 30.0 ! repose angle
    real(kind=dbPc), parameter :: tanRepostAngle = tan(repostAngle*pi/180.0)
    real(kind=dbPc), parameter :: resCoeffN = 0.9 ! normal restitution coefficient
    real(kind=dbPc), parameter :: resCoeffT = 0.0 ! tangential restitution coefficient
    real(kind=dbPc), parameter :: rhoP = 2650.0 ! particle density
    real(kind=dbPc), parameter :: nkl = 1.0 ! one particle stands for x particles
    real(kind=dbPc), parameter :: por = 0.6 ! bedform porosity

    ! bed surface

    logical, parameter :: ifPreformedSurface = .true.
    real(kind=dbPc), parameter :: initSurfElevation = 0.05 ! initial average bed height
    real(kind=dbPc), parameter :: initAmp = 0.01 !8.0*dpa ! amplitude of prerippled surface
    real(kind=dbPc), parameter :: initOmg = 4.0*pi ! wave number of prerippled surface
    real(kind=dbPc), parameter :: wavl = 2.0*pi/initOmg ! wavelength of prerippled surface
    real(kind=dbPc), parameter :: z0 = dpa/30.0 ! roughness height
    real(kind=dbPc), parameter :: blockHeight = 6.0*dpa ! bin height
    real(kind=dbPc), parameter :: initBlock = xDiff*yDiff*blockHeight ! initial block volume

    ! output after every x steps

    integer, parameter :: intervalField = 1
    integer, parameter :: intervalProfile = 1
    integer, parameter :: intervalMonitor = 1
    integer, parameter :: intervalParticle = 1
    integer, parameter :: intervalSurface = 1
    integer, parameter :: intervalStatistic = 1

    ! file

    integer, parameter :: intervalCreateFile = 1e6 ! iter num contained in a file

end module public_parameter

! **************************************

! The module contains the MPI operations

! **************************************
module parallel_operations
    use public_parameter
    implicit none
    include "mpif.h"
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
    integer :: comm, ierr
    integer :: MPI_I, MPI_D
    integer :: MPI_SG_TYPE, MPI_G_TYPE, MPI_PF_TYPE

    public :: comm, currentNode
    public :: initializeParallel, handleError, createMpiStructure, freeMpiStructure
    public :: MPI_I, MPI_D, MPI_SG_TYPE, MPI_G_TYPE, MPI_PF_TYPE

contains

    ! *********************************************************************
    ! Create the MPI Cartesian topology and find the neighbors of each node
    ! *********************************************************************
    subroutine initializeParallel
        implicit none
        integer :: nbrLeft, nbrRight

        ! create MPI Cartesian topology
        call MPI_CART_CREATE(MPI_COMM_WORLD, 1, nNodes, .true., .true., comm, ierr)
        call handleError(ierr)
        call MPI_COMM_RANK(comm, currentNode%ID, ierr)
        call handleError(ierr)
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
        call handleError(ierr)
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
    end subroutine initializeParallel

    ! *********************
    ! Handle the MPI errors
    ! *********************
    subroutine handleError(ir)
        implicit none
        integer, intent(in) :: ir
        character(len=MPI_MAX_ERROR_STRING) :: errorString
        integer :: resultLen

        if (ir /= 0) then
            print *, 'Error code: ', ir
            call MPI_ERROR_STRING(ir, errorString, resultLen, ierr)
            print *, 'Error message: ', trim(errorString)
            stop
        end if
    end subroutine handleError

    ! ******************************************************************************************
    ! Create the MPI data types for the surfaceGridType, gridType, profileType, and particleType
    ! ******************************************************************************************
    subroutine createMpiStructure
        implicit none
        integer, allocatable, dimension(:):: blockLen
        integer, allocatable, dimension(:):: disp
        integer, allocatable, dimension(:):: oldType

        ! create MPI data type for surfaceGridType
        allocate (blockLen(7))
        allocate (disp(7))
        allocate (oldType(7))
        MPI_I = MPI_INTEGER
        MPI_D = MPI_DOUBLE
        blockLen = [1, 1, 1, 1, 3, npdf, npdf]
        disp = [0, kind(0), 2*kind(0), 2*kind(0) + kind(0.0_dbPc), &
                2*kind(0) + 2*kind(0.0_dbPc), 2*kind(0) + 5*kind(0.0_dbPc), &
                2*kind(0) + 5*kind(0.0_dbPc) + npdf*kind(0.0_dbPc)]
        oldType = [MPI_I, MPI_I, MPI_D, MPI_D, MPI_D, MPI_D, MPI_D]
        call MPI_TYPE_CREATE_STRUCT(7, blockLen, disp, oldType, MPI_SG_TYPE, ierr)
        call handleError(ierr)
        call MPI_TYPE_COMMIT(MPI_SG_TYPE, ierr)
        deallocate (blockLen)
        deallocate (disp)
        deallocate (oldType)

        ! create MPI data type for gridType
        allocate (blockLen(5))
        allocate (disp(5))
        allocate (oldType(5))
        blockLen = [1, 1, 1, 3, 3]
        disp = [0, kind(0.0_dbPc), 2*kind(0.0_dbPc), 3*kind(0.0_dbPc), 6*kind(0.0_dbPc)]
        oldType = [MPI_D, MPI_D, MPI_D, MPI_D, MPI_D]
        call MPI_TYPE_CREATE_STRUCT(5, blockLen, disp, oldType, MPI_G_TYPE, ierr)
        call handleError(ierr)
        call MPI_TYPE_COMMIT(MPI_G_TYPE, ierr)
        deallocate (blockLen)
        deallocate (disp)
        deallocate (oldType)

        ! create MPI data type for profileType
        allocate (blockLen(7))
        allocate (disp(7))
        allocate (oldType(7))
        blockLen = [1, 1, 1, 1, 1, 1, 1]
        disp = [0, kind(0.0_dbPc), 2*kind(0.0_dbPc), 3*kind(0.0_dbPc), 4*kind(0.0_dbPc), 5*kind(0.0_dbPc), 6*kind(0.0_dbPc)]
        oldType = [MPI_D, MPI_D, MPI_D, MPI_D, MPI_D, MPI_D, MPI_D]
        call MPI_TYPE_CREATE_STRUCT(7, blockLen, disp, oldType, MPI_PF_TYPE, ierr)
        call handleError(ierr)
        call MPI_TYPE_COMMIT(MPI_PF_TYPE, ierr)
        deallocate (blockLen)
        deallocate (disp)
        deallocate (oldType)
    end subroutine createMpiStructure

    ! ***********************
    ! Free the MPI data types
    ! ***********************
    subroutine freeMpiStructure
        implicit none

        call MPI_TYPE_FREE(MPI_SG_TYPE, ierr)
        call MPI_TYPE_FREE(MPI_G_TYPE, ierr)
        call MPI_TYPE_FREE(MPI_PF_TYPE, ierr)
    end subroutine freeMpiStructure

end module parallel_operations

! *****************************************

! The module contains the vector operations

! *****************************************
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

    ! ************************************************************
    ! Compute the dot product of two arbitrary-dimensional vectors
    ! ************************************************************
    pure function dotProductNd(a, b)
        real(kind=dbPc), dimension(:), intent(in) :: a, b
        real(kind=dbPc) :: dotProductNd
        integer :: i
        dotProductNd = 0.0
        do i = 1, size(a)
            dotProductNd = dotProductNd + a(i)*b(i)
        end do
    end function dotProductNd

    ! *******************************************
    ! Compute the cross product of two 3D vectors
    ! *******************************************
    pure function crossProduct3d(a, b)
        real(kind=dbPc), dimension(3), intent(in) :: a, b
        real(kind=dbPc), dimension(3) :: crossProduct3d
        crossProduct3d(1) = a(2)*b(3) - a(3)*b(2)
        crossProduct3d(2) = a(3)*b(1) - a(1)*b(3)
        crossProduct3d(3) = a(1)*b(2) - a(2)*b(1)
    end function crossProduct3d

    ! ********************************************************
    ! Compute the magnitude of an arbitrary-dimensional vector
    ! ********************************************************
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

    ! **********************************************************
    ! Compute the unit vector of an arbitrary-dimensional vector
    ! **********************************************************
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

    ! ******************************************
    ! Compute the distance between two 2D points
    ! ******************************************
    pure function distance2d(p1, p2)
        real(kind=dbPc), dimension(2), intent(in) :: p1, p2
        real(kind=dbPc) :: distance2d

        distance2d = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2)
    end function distance2d

    ! ******************************************
    ! Compute the distance between two 3D points
    ! ******************************************
    pure function distance3d(p1, p2)
        real(kind=dbPc), dimension(3), intent(in) :: p1, p2
        real(kind=dbPc) :: distance3d

        distance3d = sqrt((p1(1) - p2(1))**2 + (p1(2) - p2(2))**2 + (p1(3) - p2(3))**2)
    end function distance3d

end module vector_operations

! ***************************************

! The module contains the math operations

! ***************************************
module math_operations
    use public_parameter
    implicit none
    private
    logical :: flag = .true.
    public :: valObeyCertainPDF, generateNormalDistHistogram, valObeyNormalDist

contains

    ! ********************************************************
    ! Generate a random number obeying a certain PDF histogram
    ! ********************************************************
    function valObeyCertainPDF(histogram)
        implicit none
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
        if (cumulativeProb < rand) then
            !valObeyCertainPDF = binEnd
            print *, 'Error: cumulativeProb', cumulativeProb, '< random number', rand
            stop
        end if
    end function valObeyCertainPDF

    ! **************************************************************
    ! Generate a normal distribution histogram for particle diameter
    ! **************************************************************
    subroutine generateNormalDistHistogram(histogram)
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram
        real(kind=dbPc) :: total

        ! Initialize the histogram array to zero
        histogram = 0.0

        ! Fill the histogram array
        do i = 1, npdf
            histogram(i) = exp(-0.5*((binStart + (i - 0.5)*binWidth - dpa)/dpStddDev)**2)/(sqrt(2.0*pi)*dpStddDev)
            histogram(i) = histogram(i)*binWidth
        end do
        total = sum(histogram)
        histogram = histogram/total
    end subroutine generateNormalDistHistogram

    ! ******************************************************
    ! Generate a random number obeying a normal distribution
    ! ******************************************************
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

! ******************************************

! The module contains the surface operations

! ******************************************
module surface_operations
    use public_parameter
    implicit none
    include "mpif.h"

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
    real(kind=dbPc), dimension(mx, my) :: avgD, dz
    real(kind=dbPc), dimension(mx, my, 3) :: loc
    real(kind=dbPc), dimension(mx, my, npdf) :: hist
    real(kind=dbPc), dimension(mx, my, npdf) :: dBin

    public :: generateSurfaceGrid, initializeSurface, determineParticleRollDirection, &
              updateSurfaceGrid, outputSurface
    public :: surfGrid, initDiameterDist

contains

    ! *************************************
    ! Generate the grid for the bed surface
    ! *************************************
    subroutine generateSurfaceGrid
        use parallel_operations
        implicit none
        integer :: i, j

        loc(1, :, 1) = -xDiff
        loc(:, 1, 2) = -yDiff
        loc(:, :, 3) = initSurfElevation
        do j = 1, my
            do i = 2, mx
                loc(i, j, 1) = loc(i - 1, j, 1) + xDiff
            end do
        end do
        do j = 2, my
            do i = 1, mx
                loc(i, j, 2) = loc(i, j - 1, 2) + yDiff
            end do
        end do

        ! In fact the location of surfGrid(i, j) here is the location of the south-west node of the grid
        do i = 1, mx
            do j = 1, my
                surfGrid(i, j)%location(1) = loc(i, j, 1)
                surfGrid(i, j)%location(2) = loc(i, j, 2)
                surfGrid(i, j)%location(3) = loc(i, j, 3)
            end do
        end do
    end subroutine generateSurfaceGrid

    ! ************************
    ! Initialize the bed surface
    ! ************************
    subroutine initializeSurface
        use math_operations
        implicit none
        integer :: i, j, n

        ! Initialize the diameter distribution
        select case (whichDiameterDist)
        case (0)
            if (npdf >= 3) then
                call generateNormalDistHistogram(initDiameterDist)
            else
                print *, 'Bin number (npdf) must >= 3 for normal distribution'
                stop
            end if
        case (1)
            if (npdf == 1) then
                initDiameterDist(1) = 1.0
            else
                print *, 'Bin number (npdf) must = 1 for uniform particle diameter'
                stop
            end if
        case (2)
            if (npdf == 2) then
                initDiameterDist(1) = prob1
                initDiameterDist(2) = 1.0 - initDiameterDist(1)
            else
                print *, 'Bin number (npdf) must = 2 for Bernoulli distribution'
                stop
            end if
        end select

        ! Initialize the bed surface
        do i = 1, mx
            do j = 1, my
                if (ifPreformedSurface) then
                    loc(i, j, 3) = initAmp*sin(initOmg*loc(i, j, 1)) + initSurfElevation
                    surfGrid(i, j)%location(3) = loc(i, j, 3)
                end if
                avgD(i, j) = dpa
                surfGrid(i, j)%averageDiameter = avgD(i, j)
                do n = 1, npdf
                    hist(i, j, n) = initDiameterDist(n)
                    surfGrid(i, j)%diameterDistribution(n) = hist(i, j, n)
                end do
            end do
        end do
    end subroutine initializeSurface

    ! *******************************************
    ! Determine the direction of particle rolling
    ! *******************************************
    subroutine determineParticleRollDirection
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

        do j = 1, my
            do i = 1, mx
                surfGrid(i, j)%avalanchTo = 0
                surfGrid(i, j)%avalanchFrom = 0
            end do
        end do
        do j = 2, my - 1
            do i = 2, mx - 1
                centerZ = surfGrid(i, j)%location(3)
                westZ = surfGrid(i - 1, j)%location(3)
                eastZ = surfGrid(i + 1, j)%location(3)
                northZ = surfGrid(i, j + 1)%location(3)
                southZ = surfGrid(i, j - 1)%location(3)

                slopes(1) = (centerZ - eastZ)/xDiff
                slopes(2) = (centerZ - westZ)/xDiff
                slopes(3) = (centerZ - northZ)/yDiff
                slopes(4) = (centerZ - southZ)/yDiff

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

    ! ******************************************
    ! Update the bed surface after particle move
    ! ******************************************
    subroutine updateSurfaceGrid
        use parallel_operations
        implicit none
        integer :: i, j, n, ierr
        real(kind=dbPc) :: currentBlock, patchBlock
        real(kind=dbPc), dimension(npdf) :: patchBin
        real(kind=dbPc), dimension(mx, my) :: tempDz
        real(kind=dbPc), dimension(mx, my, npdf) :: tempDBin, bin

        surfGrid%zChange = surfGrid%zChange/(xDiff*yDiff*por)
        do j = 1, my
            do i = 1, mx
                tempDz(i, j) = surfGrid(i, j)%zChange
                do n = 1, npdf
                    tempDBin(i, j, n) = surfGrid(i, j)%binVolumeChange(n)
                    bin(i, j, n) = hist(i, j, n)*initBlock
                end do
            end do
        end do
        call MPI_BARRIER(comm, ierr)
        call MPI_ALLREDUCE(tempDz, dz, mx*my, MPI_D, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(tempDBin, dBin, mx*my*npdf, MPI_D, MPI_SUM, comm, ierr)
        do j = 2, my - 1
            do i = 2, mx - 1
                loc(i, j, 3) = loc(i, j, 3) + dz(i, j)
                bin(i, j, :) = bin(i, j, :) + dBin(i, j, :)
            end do
        end do
        do i = 1, mx
            loc(i, 1, 3) = loc(i, my - 1, 3)
            loc(i, my, 3) = loc(i, 2, 3)
            bin(i, 1, :) = bin(i, my - 1, :)
            bin(i, my, :) = bin(i, 2, :)
        end do
        do j = 1, my
            loc(1, j, 3) = loc(mx - 1, j, 3)
            loc(mx, j, 3) = loc(2, j, 3)
            bin(1, j, :) = bin(mx - 1, j, :)
            bin(mx, j, :) = bin(2, j, :)
        end do

        do j = 1, my
            do i = 1, mx
                if (whichDiameterDist > 1) then
                    currentBlock = sum(bin(i, j, :))
                    patchBlock = initBlock - currentBlock
                    avgD(i, j) = 0.0
                    do n = 1, npdf
                        patchBin(n) = patchBlock*initDiameterDist(n)
                        bin(i, j, n) = (bin(i, j, n) + patchBin(n))/initBlock
                    end do
                end if
            end do
        end do

        do i = 1, mx
            do j = 1, my
                surfGrid(i, j)%location(3) = loc(i, j, 3)
                do n = 1, npdf
                    hist(i, j, n) = bin(i, j, n)/sum(bin(i, j, :))
                    surfGrid(i, j)%diameterDistribution(n) = hist(i, j, n)
                    avgD(i, j) = avgD(i, j) + hist(i, j, n)*(binStart + (n - 0.5)*binWidth)
                end do
                surfGrid(i, j)%averageDiameter = avgD(i, j)
            end do
        end do
        call MPI_BCAST(loc, mx*my*3, MPI_D, 0, comm, ierr)
        call MPI_BCAST(hist, mx*my*npdf, MPI_D, 0, comm, ierr)

    end subroutine updateSurfaceGrid

    ! **********************
    ! Output the bed surface
    ! **********************
    subroutine outputSurface(iter, t)
        use parallel_operations
        implicit none
        integer, intent(in) :: iter
        real(kind=dbPc), intent(in) :: t
        character(len=3) :: ctemp
        integer :: i, j, nameNum
        integer :: ierr

        nameNum = iter/intervalCreateFile
        write (ctemp, '(i3)') nameNum
        if (mod(iter, intervalCreateFile) == 0) then
            open (unit=14, file='./Surface/SurfaceData_'//trim(adjustl(ctemp))//'.plt')
            write (14, *) 'variables = "x", "y", "z", "dz", "d"'
            close (14)
        end if
        if (mod(iter, intervalSurface) == 0) then
            if (currentNode%ID == 0) then
                open (unit=14, position='append', file='./Surface/SurfaceData_'//trim(adjustl(ctemp))//'.plt')
                write (14, *) 'zone', ' T = "', t, '"'
                write (14, *) 'i=', mx, ' j=', my, ' datapacking=point'
                do j = 1, my
                    do i = 1, mx
                        write (14, "(5E15.4)") loc(i, j, 1), loc(i, j, 2), loc(i, j, 3), dz(i, j), avgD(i, j)
                    end do
                end do
                close (14)
            end if
        end if
    end subroutine outputSurface

end module surface_operations

! **********************************************

! The module contains the fluid field operations

! **********************************************
module field_operations
    use public_parameter
    implicit none
    include "mpif.h"
    private

    type gridType
        real(kind=dbPc) :: zDiff
        real(kind=dbPc) :: volume
        real(kind=dbPc) :: particleVolumeFraction
        ! z location starts from the elevation of bed surface
        ! The location of the bottom-south-west node of the grid
        real(kind=dbPc), dimension(3) :: vLocation
        ! The location of the grid center
        real(kind=dbPc), dimension(3) :: cLocation
    end type gridType

    type profileType
        real(kind=dbPc) :: zDiff
        ! z location starts from 0
        real(kind=dbPc) :: vlocation
        real(kind=dbPc) :: clocation
        real(kind=dbPc) :: xVelocity
        real(kind=dbPc) :: particleShearStress
        real(kind=dbPc) :: fluidShearStress
        real(kind=dbPc) :: forcingTerm
    end type profileType

    type(gridType), dimension(mx, my, nz + 1) :: grid
    type(profileType), dimension(nz) :: profile
    real(kind=dbPc), dimension(nz) :: zd, zc, zcReal, u, tau_p, tau_f, F, phi
    real(kind=dbPc), dimension(mx) :: xc, xv
    real(kind=dbPc), dimension(my) :: yc, yv
    real(kind=dbPc), dimension(nz + 1) :: zv, zvReal
    real(kind=dbPc), dimension(mx, my, nz) :: pfrac

    public :: generateGrid, initializeFluidField, calculateFluidField, &
              updateFieldGrid, outputField
    public :: grid, profile

contains

    ! ***********************
    ! Generate the fluid grid
    ! ***********************
    subroutine generateGrid
        use parallel_operations
        use surface_operations
        implicit none

        integer :: i, j, k, n
        real(kind=dbPc) :: zDiffMax, refineRatio
        real(kind=dbPc) :: bedElevation

        zv(nz) = 0.0
        zd(nz) = 0.0
        n = 0
        do while (zMax - zv(nz) > zd(nz)*2.0)
            zDiffMax = zMax/(nz - n)
            refineRatio = (zDiffMax/zDiffMin)**(1.0/(nzUni - 1))
            zd(1) = zDiffMin
            zv(1) = 0.0
            do k = 2, nz
                if (k <= nzUni) then
                    zd(k) = zd(k - 1)*refineRatio
                else
                    zd(k) = zDiffMax
                end if
                zv(k) = zv(k - 1) + zd(k - 1)
                zc(k) = zv(k) + zd(k)*0.5
            end do
            n = n + 1
        end do
        zd(nz) = zMax - zv(nz)
        zc(nz) = zMax - zd(nz)*0.5
        zv(nz + 1) = zMax

        do i = 1, mx
            xv(i) = surfGrid(i, 2)%location(1)
            xc(i) = xv(i) + xDiff*0.5
        end do
        currentNode%sx = xv(currentNode%i1)
        currentNode%ex = xv(currentNode%im) + xDiff
        do j = 1, my
            yv(j) = surfGrid(2, j)%location(2)
            yc(j) = yv(j) + yDiff*0.5
        end do

        do i = 1, mx
            do j = 1, my
                bedElevation = surfGrid(i, j)%location(3)
                do k = 1, nz
                    zvReal(k) = bedElevation + zv(k)
                    zcReal(k) = bedElevation + zc(k)
                    profile(k)%zDiff = zd(k)
                    profile(k)%vlocation = zv(k)
                    profile(k)%clocation = zc(k)
                    grid(i, j, k)%zDiff = zd(k)
                    grid(i, j, k)%vLocation(1) = xv(i)
                    grid(i, j, k)%vLocation(2) = yv(j)
                    grid(i, j, k)%vLocation(3) = zvReal(k)
                    grid(i, j, k)%clocation(1) = xc(i)
                    grid(i, j, k)%clocation(2) = yc(j)
                    grid(i, j, k)%clocation(3) = zcReal(k)
                end do
            end do
        end do
        do i = 1, mx
            do j = 1, my
                bedElevation = surfGrid(i, j)%location(3)
                zvReal(nz + 1) = bedElevation + zMax
                grid(i, j, nz + 1)%zDiff = zd(nz)
                grid(i, j, nz + 1)%vLocation(1) = xv(i)
                grid(i, j, nz + 1)%vLocation(2) = yv(j)
                grid(i, j, nz + 1)%vLocation(3) = zvReal(nz + 1)
                grid(i, j, nz + 1)%clocation(1) = xc(i)
                grid(i, j, nz + 1)%clocation(2) = yc(j)
                grid(i, j, nz + 1)%clocation(3) = zvReal(nz + 1) + zd(nz)*0.5
            end do
        end do
    end subroutine generateGrid

    ! **************************
    ! Initialize the fluid field
    ! **************************
    subroutine initializeFluidField
        implicit none
        integer :: i, j, k

        do i = 1, mx
            do j = 1, my
                do k = 1, nz
                    pfrac(i, j, k) = 0.0
                    grid(i, j, k)%volume = xDiff*yDiff*zd(k)
                    grid(i, j, k)%particleVolumeFraction = pfrac(i, j, k)
                end do
            end do
        end do
        do i = 1, mx
            do j = 1, my
                grid(i, j, nz + 1)%volume = xDiff*yDiff*zd(nz)
                grid(i, j, nz + 1)%particleVolumeFraction = 0.0
            end do
        end do
        do k = 1, nz
            if (zc(k) > z0) then
                u(k) = uStar/kapa*log(zc(k)/z0)
            else
                u(k) = 0.0
            end if
            tau_p(k) = 0.0
            tau_f(k) = rho*uStar**2
            F(k) = 0.0
            profile(k)%xVelocity = u(k)
            profile(k)%particleShearStress = tau_p(k)
            profile(k)%fluidShearStress = tau_f(k)
            profile(k)%forcingTerm = F(k)
        end do
    end subroutine initializeFluidField

    ! ****************************************************
    ! Update the 1D fluid field using the Thomas algorithm
    ! ****************************************************
    subroutine calculateFluidField
        use parallel_operations
        implicit none
        integer :: i, j, k, ierr
        real(kind=dbPc) :: ap, at, ab, b
        real(kind=dbPc) :: dzP, dzT, dzB, uT, uB
        real(kind=dbPc) :: nut, nutot, dudz, mixl
        real(kind=dbPc) :: wtt, wtb, wbb, wbt
        real(kind=dbPc), dimension(nz) :: p, q
        real(kind=dbPc), dimension(nz) :: tempF
        real(kind=dbPc), dimension(mx, my, nz) :: tempPfrac
        real(kind=dbPc), dimension(3*nz) :: buffer

        do k = 1, nz
            tempF(k) = profile(k)%forcingTerm
        end do
        do j = 1, my
            do i = 1, mx
                do k = 1, nz
                    tempPfrac(i, j, k) = grid(i, j, k)%particleVolumeFraction
                end do
            end do
        end do
        call MPI_BARRIER(comm, ierr)
        call MPI_ALLREDUCE(tempF, F, nz, MPI_D, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(tempPfrac, pfrac, mx*my*nz, MPI_D, MPI_SUM, comm, ierr)

        ! Calculate the shear stress
        phi(nz) = sum(pfrac(2:mx - 1, 2:my - 1, nz))/(nx*ny)
        tau_p(nz) = F(nz)/(area*(1.0 - phi(nz)))
        tau_f(nz) = rho*uStar**2 + tau_p(nz)
        do k = nz - 1, 1, -1
            phi(k) = sum(pfrac(2:mx - 1, 2:my - 1, k))/(nx*ny)
            tau_p(k) = tau_p(k + 1) + F(k)/(area*(1.0 - phi(k)))
            tau_f(k) = rho*uStar**2 + tau_p(k)
        end do

        ! Calculate the velocity profile
        ! Forward elimination
        dzP = zd(1)
        dzT = (zd(1) + zd(2))/2.0
        wtt = 0.5*zd(1)/dzT
        wtb = 0.5*zd(2)/dzT
        uT = u(2)*wtt + u(1)*wtb
        uB = 0.0
        mixl = kapa*zc(1)*(1.0 - exp(-1.0/26.0*zc(1)*uStar/nu))
        dudz = (uT - uB)/dzP
        nut = mixl**2*abs(dudz)
        nutot = nu + nut
        ap = nutot/(dzP*dzT) + 1.0/dt
        at = -nutot/(dzP*dzT)
        b = F(1)/(area*dzP*rho*(1.0 - phi(1))) + u(1)/dt - tau_f(1)/dzP
        p(1) = -at/ap
        q(1) = b/ap
        do k = 2, nz - 1
            dzP = zd(k)
            dzT = (zd(k) + zd(k + 1))/2.0
            dzB = (zd(k) + zd(k - 1))/2.0
            wtt = 0.5*zd(k)/dzT
            wtb = 0.5*zd(k + 1)/dzT
            wbb = 0.5*zd(k)/dzB
            wbt = 0.5*zd(k - 1)/dzB
            uT = u(k + 1)*wtt + u(k)*wtb
            uB = u(k - 1)*wbb + u(k)*wbt
            mixl = kapa*zc(k)*(1.0 - exp(-1.0/26.0*zc(k)*uStar/nu))
            dudz = (uT - uB)/dzP
            nut = mixl**2*abs(dudz)
            nutot = nu + nut
            ap = nutot/(dzP*dzT) + nutot/(dzP*dzB) + 1.0/dt
            at = -nutot/(dzP*dzT)
            ab = -nutot/(dzP*dzB)
            b = F(k)/(area*dzP*rho*(1.0 - phi(k))) + u(k)/dt
            p(k) = -at/(ap + ab*p(k - 1))
            q(k) = (b - ab*q(k - 1))/(ap + ab*p(k - 1))
        end do
        dzP = zd(nz)
        dzB = (zd(nz) + zd(nz - 1))/2.0
        wbb = 0.5*zd(nz)/dzB
        wbt = 0.5*zd(nz - 1)/dzB
        mixl = kapa*zc(nz)*(1.0 - exp(-1.0/26.0*zc(nz)*uStar/nu))
        nutot = nu
        ap = nutot/(dzP*dzB) + 1.0/dt
        ab = -nutot/(dzP*dzB)
        b = F(nz)/(area*dzP*rho*(1.0 - phi(nz))) + u(nz)/dt + rho*uStar**2/dzP
        q(nz) = (b - ab*q(nz - 1))/(ap + ab*p(nz - 1))

        ! Backward substitution
        u(nz) = q(nz)
        do k = nz - 1, 1, -1
            u(k) = p(k)*u(k + 1) + q(k)
        end do

        do k = 1, nz
            buffer(k) = u(k)
            buffer(k + nz) = tau_p(k)
            buffer(k + 2*nz) = tau_f(k)
        end do

        call MPI_BARRIER(comm, ierr)
        call MPI_BCAST(buffer, 3*nz, MPI_D, 0, comm, ierr)
        call handleError(ierr)
        do k = 1, nz
            u(k) = buffer(k)
            tau_p(k) = buffer(k + nz)
            tau_f(k) = buffer(k + 2*nz)
            profile(k)%xVelocity = u(k)
            profile(k)%particleShearStress = tau_p(k)
            profile(k)%fluidShearStress = tau_f(k)
        end do
    end subroutine calculateFluidField

    ! ***************************************************
    ! Update the field grid after the bed surface changes
    ! ***************************************************
    subroutine updateFieldGrid
        use surface_operations
        implicit none
        integer :: i, j, k
        real(kind=dbPc) :: bedElevation

        do i = 1, mx
            do j = 1, my
                bedElevation = surfGrid(i, j)%location(3)
                do k = 1, nz
                    zvReal(k) = bedElevation + zv(k)
                    zcReal(k) = bedElevation + zc(k)
                    grid(i, j, k)%vLocation(3) = zvReal(k)
                    grid(i, j, k)%clocation(3) = zcReal(k)
                end do
                zvReal(nz + 1) = bedElevation + zMax
                grid(i, j, nz + 1)%vLocation(3) = zvReal(nz + 1)
                grid(i, j, nz + 1)%clocation(3) = zvReal(nz + 1) + zd(nz)*0.5
            end do
        end do
    end subroutine updateFieldGrid

    ! **********************
    ! Output the fluid field
    ! **********************
    subroutine outputField(iter, t)
        use parallel_operations
        implicit none
        integer, intent(in) :: iter
        real(kind=dbPc), intent(in) :: t
        character(len=3) :: ctemp
        integer :: i, j, k
        integer :: nameNum

        ! Output the profile data
        if (mod(iter, intervalProfile) == 0) then
            if (currentNode%ID == 0) then
                open (unit=12, position='append', file='./Field/Profile.plt')
                write (12, *) 'zone', ' T = "', t, '"'
                do k = 1, nz
                    write (12, "(6E15.4)") zc(k), u(k), -tau_p(k), tau_f(k), F(k), phi(k)
                end do
                close (12)
            end if
        end if

        ! Output the field data
        nameNum = iter/intervalCreateFile
        write (ctemp, '(i3)') nameNum
        if (mod(iter, intervalCreateFile) == 0) then
            open (unit=13, file='./Field/FieldData_'//trim(adjustl(ctemp))//'.plt')
            write (13, *) 'variables = "x", "y", "z", "Phi_p"'
            close (13)
        end if
        if (mod(iter, intervalField) == 0) then
            if (currentNode%ID == 0) then
                open (unit=13, position='append', file='./Field/FieldData_'//trim(adjustl(ctemp))//'.plt')
                write (13, *) 'zone', ' T = "', t, '"'
                write (13, *) 'i=', mx, ' j=', my, ' k=', nz, ' datapacking=point'
                do k = 1, nz
                    do j = 1, my
                        do i = 1, mx
                            write (13, "(4E15.4)") xc(i), yc(j), grid(i, j, k)%cLocation(3), pfrac(i, j, k)
                        end do
                    end do
                end do
                close (13)
            end if
        end if
    end subroutine outputField

end module field_operations

! *******************************************

! The module contains the particle operations

! *******************************************
module particle_operations
    use public_parameter
    implicit none
    include "mpif.h"

    private

    type particleType
        integer, dimension(3) :: indices
        real(kind=dbPc) :: diameter
        real(kind=dbPc) :: altitude
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc), dimension(3) :: velocity
    end type particleType

    type(particleType), allocatable, dimension(:) :: particle
    integer pNum

    public :: initializeParticle, calculateParColl, calculateParticleMovement, &
              reallocateParticle, outputParticle
    public :: particle, pNum

contains

    ! ************************
    ! Initialize the particles
    ! ************************
    subroutine initializeParticle
        use parallel_operations
        use math_operations
        use surface_operations
        implicit none
        integer :: n
        real(kind=dbPc) :: rand1, rand2, rand3
        type(particleType) :: currentParticle

        if (binStart < 0.0) then
            print *, 'The starting point of the bin must be >= 0.0'
            stop
        end if
        pNum = pNumInit
        allocate (particle(pNum))
        do n = 1, pNum
            currentParticle = particle(n)
            call random_number(rand1)
            call random_number(rand2)
            call random_number(rand3)
            currentParticle%location(1) = xMax*rand1
            currentParticle%location(2) = yMax*rand2
            currentParticle%location(3) = 0.1*dpa + initSurfElevation !zMax*rand3*0.1 + initSurfElevation
            currentParticle%velocity = [1.0, 0.0, -0.1]
            currentParticle%diameter = valObeyCertainPDF(initDiameterDist)
            call determineParIJK(currentParticle)
            particle(n) = currentParticle
        end do

    end subroutine initializeParticle

    ! *************************************
    ! Determine the indices of the particle
    ! *************************************
    subroutine determineParIJK(currentP)
        use field_operations
        implicit none

        integer :: ip, jp, kp, tempKp
        real(kind=dbPc) :: currentZ
        type(particleType) :: currentP

        ip = floor(currentP%location(1)/xDiff) + 2
        jp = floor(currentP%location(2)/yDiff) + 2
        if (ip > mx - 1) then
            ip = ip - (mx - 1) + 1
            currentP%location(1) = currentP%location(1) - xMax
        elseif (ip < 2) then
            ip = mx - (2 - ip)
            currentP%location(1) = currentP%location(1) + xMax
        end if
        if (jp > my - 1) then
            jp = jp - (my - 1) + 1
            currentP%location(2) = currentP%location(2) - yMax
        elseif (jp < 2) then
            jp = my - (2 - jp)
            currentP%location(2) = currentP%location(2) + yMax
        end if
        do tempKp = 1, nz
            kp = tempKp
            currentZ = grid(ip, jp, tempKp + 1)%vLocation(3)
            if (currentP%location(3) < currentZ) then
                exit
            end if
        end do
        currentP%indices(1) = ip
        currentP%indices(2) = jp
        currentP%indices(3) = kp
    end subroutine determineParIJK

    ! *******************************************************
    ! Calculate particle-bed and particle-particle collisions
    ! *******************************************************
    subroutine calculateParColl
        use math_operations
        use parallel_operations
        use field_operations
        use surface_operations
        use vector_operations
        implicit none

        integer :: n, nadd, n2
        integer :: totalEjectNum
        integer :: i, j, k
        integer :: ip, jp
        integer :: whichTriangle
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
        real(kind=dbPc) :: sigma, lambda, mu
        real(kind=dbPc) :: d1, d2
        real(kind=dbPc) :: v1, v2
        real(kind=dbPc) :: theta2
        real(kind=dbPc) :: rand1, rand2, rand3
        real(kind=dbPc) :: eBar
        real(kind=dbPc) :: m1, m2
        real(kind=dbPc) :: E1, E2, Ed2, Eeff, E2Bar
        real(kind=dbPc) :: tau_s, tau_fw
        real(kind=dbPc) :: ejectVolume, rollVolume
        real(kind=dbPc) :: distance12, contactDistance
        real(kind=dbPc) :: eta1, eta2, alpha1, alpha2, beta1, beta2
        real(kind=dbPc) :: relativeV12Normal, relativeV12Tangent
        real(kind=dbPc) :: resCoeffTMidAir
        real(kind=dbPc), dimension(4) :: adjacentSurfGridZ
        real(kind=dbPc), dimension(3) :: vertex1, vertex2, vertex3
        real(kind=dbPc), dimension(3) :: vector12
        real(kind=dbPc), dimension(3) :: surfaceNormalVector
        real(kind=dbPc), dimension(3) :: particleProjection
        real(kind=dbPc), dimension(3) :: impactVelocity, reboundVelocity
        real(kind=dbPc), dimension(3) :: impactCoordinateX, impactCoordinateY, impactCoordinateZ
        real(kind=dbPc), dimension(3) :: relativeV12
        type(particleType) :: currentParticle, currentParticle2
        type(particleType), dimension(maxEjectNum) :: addParticle
        type(particleType), allocatable, dimension(:) :: tempParticle

        ! initialization
        allocate (tempParticle(pNum + maxEjectNum))
        do i = 1, mx
            do j = 1, my
                surfGrid(i, j)%zChange = 0.0_dbPc
                do k = 1, npdf
                    surfGrid(i, j)%binVolumeChange(k) = 0.0_dbPc
                end do
            end do
        end do
        currentTotalNum = 0
        totalEjectNum = 0

        ! particle-bed collision
        do n = 1, pNum
            currentParticle = particle(n)
            ip = currentParticle%indices(1)
            jp = currentParticle%indices(2)

            ! Estimate the altitude of the particle
            adjacentSurfGridZ(1) = surfGrid(ip, jp)%location(3)
            adjacentSurfGridZ(2) = surfGrid(ip + 1, jp)%location(3)
            adjacentSurfGridZ(3) = surfGrid(ip, jp + 1)%location(3)
            adjacentSurfGridZ(4) = surfGrid(ip + 1, jp + 1)%location(3)
            estimateAltitude = currentParticle%location(3) - maxval(adjacentSurfGridZ)
            if (estimateAltitude > 0.5*currentParticle%diameter) then
                currentTotalNum = currentTotalNum + 1
                currentParticle%altitude = estimateAltitude
                tempParticle(currentTotalNum) = currentParticle
                cycle
            end if

            ! Determine which triangle the particle will collide with
            call determineSurfTri(currentParticle, ip, jp, whichTriangle, vertex1, vertex2, vertex3, &
                                  surfaceNormalVector)

            ! Calculate the altitude of the particle
            particleProjection(1) = currentParticle%location(1)
            particleProjection(2) = currentParticle%location(2)
            particleProjection(3) = vertex3(3) + &
                                    (surfaceNormalVector(1)*(currentParticle%location(1) - vertex3(1)) + &
                                     surfaceNormalVector(2)*(currentParticle%location(2) - vertex3(2)))/ &
                                    (-surfaceNormalVector(3))
            currentParticle%altitude = currentParticle%location(3) - particleProjection(3)

            ! Process the impact event
            impactCoordinateZ = surfaceNormalVector
            impactVelocity(3) = dotProduct(currentParticle%velocity, impactCoordinateZ)
            if (currentParticle%altitude < 0.5*currentParticle%diameter .and. impactVelocity(3) < 0.0) then
                ! Change to local coordinate system
                impactCoordinateY = unitVector(crossProduct(impactCoordinateZ, currentParticle%velocity))
                impactCoordinateX = unitVector(crossProduct(impactCoordinateY, impactCoordinateZ))
                impactVelocity(1) = dotProduct(currentParticle%velocity, impactCoordinateX)
                impactVelocity(2) = dotProduct(currentParticle%velocity, impactCoordinateY)

                ! Find the closest vertex to the particle
                call findClosestVertex(particleProjection, ip, jp, vertex1, vertex2, vertex3, &
                                       whichTriangle, whichVertex, closestIP, closestJP)

                ! Particle and bed information
                d1 = currentParticle%diameter
                d2 = surfGrid(closestIP, closestJP)%averageDiameter
                m1 = rhoP*(pi*d1**3)/6.0
                m2 = rhoP*(pi*d2**3)/6.0
                v1 = vectorMagnitude(impactVelocity)

                ! Calculate rebound
                call calculateReboundPara(impactVelocity, d1, d2, eBar, theta2)
                v2 = v1*eBar
                reboundVelocity(1) = v2*cos(theta2)
                reboundVelocity(2) = 0.0
                reboundVelocity(3) = v2*sin(theta2)
                E2 = (m2*reboundVelocity(3)**2)*0.5
                Ed2 = m2*9.8*d2
                if (theta2 <= 0.0 .or. E2 < Ed2 .or. eBar < 0.0) then ! no rebound
                    call findChangedGridNode(1, closestIP, closestJP, changedIP, changedJP)
                    surfGrid(changedIP, changedJP)%zChange = surfGrid(changedIP, changedJP)%zChange &
                                                             + (pi*d1**3)/6.0
                    if (whichDiameterDist /= 1) then
                        whichBin = floor((currentParticle%diameter - binStart)/binWidth) + 1
                        whichBin = max(whichBin, 1)
                        whichBin = min(whichBin, npdf)
                        surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) = &
                            surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) + (pi*d1**3)/6.0
                    end if
                else ! rebound
                    currentTotalNum = currentTotalNum + 1
                    currentParticle%altitude = 0.5*d1
                    currentParticle%location(3) = particleProjection(3) + 0.5*d1
                    currentParticle%velocity(1) = dotProduct(reboundVelocity, impactCoordinateX)
                    currentParticle%velocity(2) = dotProduct(reboundVelocity, impactCoordinateY)
                    currentParticle%velocity(3) = dotProduct(reboundVelocity, impactCoordinateZ)
                    call determineParIJK(currentParticle)
                    tempParticle(currentTotalNum) = currentParticle
                end if

                ! Calculate eject
                ! Calculate eject number according to Lammel et al. 2017
                tau_s = rho*0.0123*(rhoP/rho*9.8*d2 + 3.0e-4/(rho*d2)) ! Shao and Lu 2000
                tau_fw = profile(1)%fluidShearStress
                if (tau_fw > 0.9*tau_s) then
                    Eeff = 0.1*Ed2
                else
                    Eeff = Ed2*(1.0 - tau_fw/tau_s)
                end if
                E1 = (m1*v1**2)/2.0
                lambda = 2.0*log((1.0 - eBar**2)*E1/Ed2)
                sigma = sqrt(lambda)*log(2.0)
                mu = log((1.0 - eBar**2)*E1) - lambda*log(2.0)
                E2Bar = Ed2*((1.0 - eBar**2)*E1/Ed2)**(1.0 - (2.0 - log(2.0))*log(2.0))
                ejectNum = floor(0.06*((1.0 - eBar**2)*E1/(2.0*E2Bar))*erfc((log(Eeff) - mu)/(sqrt(2.0)*sigma)))
                ! process eject particles and surface grid change
                if (ejectNum > 0) then
                    ejectVolume = 0.0
                    do nadd = 1, ejectNum
                        call random_number(rand1)
                        ! E2 = random number obeying log normal distribution with mean mu and standard deviation sigma
                        E2 = exp(valObeyNormalDist(mu, sigma))
                        v2 = sqrt(2.0*E2/m2)
                        d2 = valObeyCertainPDF(surfGrid(closestIP, closestJP)%diameterDistribution)
                        m2 = rhoP*(pi*d2**3)/6.0
                        rand3 = -1.0
                        do while (rand3 < 0.0)
                            call random_number(rand1)
                            call random_number(rand2)
                            rand3 = 1.0 - rand1 - rand2
                        end do
                        totalEjectNum = totalEjectNum + 1
                        addParticle(totalEjectNum)%indices(1) = ip
                        addParticle(totalEjectNum)%indices(2) = jp
                        addParticle(totalEjectNum)%indices(3) = 1
                        addParticle(totalEjectNum)%diameter = d2
                        addParticle(totalEjectNum)%altitude = v2*dt
                        addParticle(totalEjectNum)%location(1) = rand1*vertex1(1) + rand2*vertex2(1) + rand3*vertex3(1)
                        addParticle(totalEjectNum)%location(2) = rand1*vertex1(2) + rand2*vertex2(2) + rand3*vertex3(2)
                        addParticle(totalEjectNum)%location(3) = rand1*vertex1(3) + rand2*vertex2(3) + rand3*vertex3(3) + v2*dt
                        addParticle(totalEjectNum)%velocity(1) = 0.0
                        addParticle(totalEjectNum)%velocity(2) = 0.0
                        addParticle(totalEjectNum)%velocity(3) = v2
                        ejectVolume = ejectVolume + (pi*d2**3)/6.0
                        if (whichDiameterDist /= 1) then
                            whichBin = floor((d2 - binStart)/binWidth) + 1
                            whichBin = max(whichBin, 1)
                            whichBin = min(whichBin, npdf)
                            surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) = &
                                surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) - (pi*d2**3)/6.0
                        end if
                    end do
                    surfGrid(closestIP, closestJP)%zChange = surfGrid(closestIP, closestJP)%zChange &
                                                             - ejectVolume
                    ! Process particle roll due to avalanche
                    call findChangedGridNode(2, closestIP, closestJP, changedIP, changedJP)
                    rollVolume = 0.0
                    if (closestIP /= changedIP .or. closestJP /= changedJP) then
                        do while (rollVolume < ejectVolume)
                            d2 = valObeyCertainPDF(surfGrid(changedIP, changedJP)%diameterDistribution)
                            rollVolume = rollVolume + (pi*d2**3)/6.0
                            if (whichDiameterDist /= 1) then
                                whichBin = floor((d2 - binStart)/binWidth) + 1
                                whichBin = max(whichBin, 1)
                                whichBin = min(whichBin, npdf)
                                surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) = &
                                    surfGrid(changedIP, changedJP)%binVolumeChange(whichBin) - (pi*d2**3)/6.0
                                surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) = &
                                    surfGrid(closestIP, closestJP)%binVolumeChange(whichBin) + (pi*d2**3)/6.0
                            end if
                        end do
                        surfGrid(changedIP, changedJP)%zChange = surfGrid(changedIP, changedJP)%zChange &
                                                                 - rollVolume
                        surfGrid(closestIP, closestJP)%zChange = surfGrid(closestIP, closestJP)%zChange &
                                                                 + rollVolume
                    end if
                end if
            else
                currentTotalNum = currentTotalNum + 1
                tempParticle(currentTotalNum) = currentParticle
            end if
        end do
        pNum = currentTotalNum

        ! Calculate particle-particle collision
        if (ifMidairCollision) then
            allocate (globalN(mx, my, pNum/(currentNode%nx*ny)*2))
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

        ! Reallocate the particle array
        if (totalEjectNum > maxEjectNum) then
            print *, "The eject number in one time step", totalEjectNum, ">", maxEjectNum
            stop
        end if
        if (totalEjectNum > 0) then
            pNum = pNum + totalEjectNum
            tempParticle(currentTotalNum + 1:pNum) = addParticle(1:totalEjectNum)
        end if
        deallocate (particle)
        allocate (particle(pNum))
        particle = tempParticle(1:pNum)
        deallocate (tempParticle)
    end subroutine calculateParColl

    ! ***************************************************************
    ! Determine which surface triangle the particle will collide with
    ! ***************************************************************
    subroutine determineSurfTri(currentP, i, j, whichTri, v1, v2, v3, nVec)
        use surface_operations
        use vector_operations
        implicit none
        type(particleType), intent(in) :: currentP
        integer, intent(in) :: i, j
        integer, intent(out) :: whichTri
        real(kind=dbPc), dimension(3), intent(out) :: v1, v2, v3, nVec
        real(kind=dbPc), dimension(2) :: LocalLoc ! particle local location
        real(kind=dbPc), dimension(3) :: vector12, vector13

        LocalLoc(1) = currentP%location(1) - surfGrid(i, j)%location(1)
        LocalLoc(2) = currentP%location(2) - surfGrid(i, j)%location(2)
        if (mod(i + j, 2) == 0) then
            ! ip+jp = Even num
            !   3 ------- 2
            !     |    /|
            !     |  /  |
            !     |/    |
            !   1 ------- 3
            ! Left triangle: whichTriangle=1
            ! Right triangle: whichTriangle=2
            ! Current surfGrid%location is at point 1, it is the local origin
            v1 = surfGrid(i, j)%location
            v2 = surfGrid(i + 1, j + 1)%location
            if (LocalLoc(1)/LocalLoc(2) <= xDiff/yDiff) then
                whichTri = 1
                v3 = surfGrid(i, j + 1)%location
                vector12 = v2 - v1
                vector13 = v3 - v1
                nVec = unitVector(crossProduct(vector12, vector13))
            else
                whichTri = 2
                v3 = surfGrid(i + 1, j)%location
                vector12 = v2 - v1
                vector13 = v3 - v1
                nVec = unitVector(crossProduct(vector13, vector12))
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
            v1 = surfGrid(i + 1, j)%location
            v2 = surfGrid(i, j + 1)%location
            if (LocalLoc(1)/(yDiff - LocalLoc(2)) <= xDiff/yDiff) then
                whichTri = 3
                v3 = surfGrid(i, j)%location
                vector12 = v2 - v1
                vector13 = v3 - v1
                nVec = unitVector(crossProduct(vector12, vector13))
            else
                whichTri = 4
                v3 = surfGrid(i + 1, j + 1)%location
                vector12 = v2 - v1
                vector13 = v3 - v1
                nVec = unitVector(crossProduct(vector13, vector12))
            end if
        end if
    end subroutine determineSurfTri

    ! ***************************************
    ! Find the closest vertex to the particle
    ! ***************************************
    subroutine findClosestVertex(particleP, i, j, v1, v2, v3, whichTri, whichV, closestI, closestJ)
        use vector_operations
        implicit none
        real(kind=dbPc), dimension(3), intent(in) :: particleP
        real(kind=dbPc), dimension(3), intent(in) :: v1, v2, v3
        integer, intent(in) :: i, j, whichTri
        integer, intent(out) :: whichV, closestI, closestJ
        integer :: k
        real(kind=dbPc), dimension(3) :: particleVertexDistance
        real(kind=dbPc) :: minDistance

        particleVertexDistance(1) = distance3d(particleP, v1)
        particleVertexDistance(2) = distance3d(particleP, v2)
        particleVertexDistance(3) = distance3d(particleP, v3)
        minDistance = particleVertexDistance(1)
        whichV = 1
        do k = 2, 3
            if (particleVertexDistance(k) < minDistance) then
                minDistance = particleVertexDistance(k)
                whichV = k
            end if
        end do
        select case (whichTri)
            !   3 ------- 2
            !     |    /|
            !     |  /  |
            !     |/    |
            !   1 ------- 3
            ! Left triangle: whichTriangle=1
            ! Right triangle: whichTriangle=2
        case (1)
            select case (whichV)
            case (1)
                closestI = i
                closestJ = j
            case (2)
                closestI = i + 1
                closestJ = j + 1
            case (3)
                closestI = i
                closestJ = j + 1
            end select
        case (2)
            select case (whichV)
            case (1)
                closestI = i
                closestJ = j
            case (2)
                closestI = i + 1
                closestJ = j + 1
            case (3)
                closestI = i + 1
                closestJ = j
            end select
            !   2 ------- 3
            !     |\    |
            !     |  \  |
            !     |    \|
            !   3 ------- 1
            ! Left triangle: whichTriangle=3
            ! Right triangle: whichTriangle=4
        case (3)
            select case (whichV)
            case (1)
                closestI = i + 1
                closestJ = j
            case (2)
                closestI = i
                closestJ = j + 1
            case (3)
                closestI = i
                closestJ = j
            end select
        case (4)
            select case (whichV)
            case (1)
                closestI = i + 1
                closestJ = j
            case (2)
                closestI = i
                closestJ = j + 1
            case (3)
                closestI = i + 1
                closestJ = j + 1
            end select
        end select
        if (closestI >= mx) then
            closestI = closestI - mx + 2
        end if
        if (closestJ >= my) then
            closestJ = closestJ - my + 2
        end if
    end subroutine findClosestVertex

    ! **********************************************************
    ! Calculate the rebound parameters of particle-bed collision
    ! **********************************************************
    subroutine calculateReboundPara(impactV, d1, d2, eBar, theta2)
        use vector_operations
        implicit none
        real(kind=dbPc), intent(in) :: d1, d2
        real(kind=dbPc), dimension(3), intent(in) :: impactV
        real(kind=dbPc), intent(out) :: eBar, theta2
        integer :: iterateNum
        real(kind=dbPc) :: theta1, theta2Min, theta2Max, theta2Mid
        real(kind=dbPc) :: pT2T1, pMax
        real(kind=dbPc) :: rand1, rand2
        real(kind=dbPc) :: pointX, pointY
        real(kind=dbPc) :: nonDimD1, nonDimD2
        real(kind=dbPc) :: eta, alpha, beta, gama

        ! Splash function of Lammel et al. 2017
        nonDimD1 = 2.0*d1/(d1 + d2)
        nonDimD2 = 2.0*d2/(d1 + d2)
        eta = resCoeffN*nonDimD1**3/(nonDimD1**3 + resCoeffN*nonDimD2**3)
        alpha = (1.0 + resCoeffN)/(1.0 + eta) - 1.0
        beta = 1.0 - (2.0/7.0)*(1.0 - resCoeffT)/(1.0 + eta)
        theta1 = atan(abs(impactV(3)/impactV(1)))
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
    end subroutine calculateReboundPara

    ! *******************************************
    ! Find the changed grid node due to avalanche
    ! *******************************************
    subroutine findChangedGridNode(tag, i, j, ci, cj)
        use surface_operations
        implicit none
        integer, intent(in) :: i, j, tag
        integer, intent(out) :: ci, cj
        integer :: direction

        if (tag == 1) then
            direction = surfGrid(i, j)%avalanchTo
        else
            direction = surfGrid(i, j)%avalanchFrom
        end if

        select case (direction)
        case (1)
            ci = i + 1
            cj = j
        case (2)
            ci = i - 1
            cj = j
        case (3)
            ci = i
            cj = j + 1
        case (4)
            ci = i
            cj = j - 1
        case default
            ci = i
            cj = j
        end select
        if (ci >= mx) then
            ci = ci - mx + 2
        elseif (ci <= 1) then
            ci = ci + mx - 2
        end if
        if (cj >= my) then
            cj = cj - my + 2
        elseif (cj <= 1) then
            cj = cj + my - 2
        end if
    end subroutine findChangedGridNode

    ! ***********************************
    ! Calculate the movement of particles
    ! ***********************************
    subroutine calculateParticleMovement
        use field_operations
        implicit none
        integer :: n, ip, jp, kp, i, j, k
        real(kind=dbPc) :: dp, mp
        real(kind=dbPc) :: fDrag
        real(kind=dbPc), dimension(3) :: up, uf
        real(kind=dbPc), dimension(3) :: bulkForce, totalForce
        real(kind=dbPc), dimension(3) :: u1, u2, u3, u4
        real(kind=dbPc), dimension(3) :: a1, a2, a3, a4
        type(particleType) :: currentParticle

        do k = 1, nz
            profile(k)%forcingTerm = 0.0
            do j = 1, my
                do i = 1, mx
                    grid(i, j, k)%particleVolumeFraction = 0.0
                end do
            end do
        end do
        do n = 1, pNum
            currentParticle = particle(n)
            up = currentParticle%velocity
            dp = currentParticle%diameter
            kp = currentParticle%indices(3)
            uf(1) = profile(kp)%xVelocity
            uf(2) = 0.0
            uf(3) = 0.0
            mp = rhoP*(pi*dp**3)/6.0
            bulkForce(1) = 0.0
            bulkForce(2) = 0.0
            bulkForce(3) = -9.8*mp*(1.0 - rho/rhoP)

            ! use 4th order Runge-Kutta method to solve the ODE
            u1 = up
            totalForce = dragForce()
            fDrag = -totalForce(1)
            totalForce = totalForce + bulkForce
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
            call determineParIJK(currentParticle)
            ip = currentParticle%indices(1)
            jp = currentParticle%indices(2)
            kp = currentParticle%indices(3)
            particle(n) = currentParticle
            profile(kp)%forcingTerm = profile(kp)%forcingTerm + fDrag
            grid(ip, jp, kp)%particleVolumeFraction = grid(ip, jp, kp)%particleVolumeFraction + &
                                                      pi*dp**3/6.0/grid(ip, jp, kp)%volume
        end do

    contains

        ! ***********************************************
        ! Calculate the drag force acting on the particle
        ! ***********************************************
        function dragForce() result(Fd)
            implicit none
            real(kind=dbPc), dimension(3) :: Fd
            real(kind=dbPc), dimension(3) :: C_d, Re_p, relativeU

            do i = 1, 3
                relativeU(i) = uf(i) - up(i)
                if (abs(relativeU(i)) > 1.0e-8) then
                    Re_p(i) = abs(relativeU(i))*dp/nu
                    C_d(i) = (sqrt(0.5) + sqrt(24.0/Re_p(i)))**2
                    Fd(i) = pi/8.0*rho*dp**2*C_d(i)*abs(relativeU(i))*relativeU(i)
                else
                    Fd(i) = 0.0
                end if
            end do
        end function dragForce
    end subroutine calculateParticleMovement

    ! ********************************************
    ! Exchange particles between neighboring nodes
    ! ********************************************
    subroutine reallocateParticle
        use parallel_operations
        implicit none
        integer :: n, ip, jp, kp, currentN, ierr
        integer :: sendRightNum, recvLeftNum
        integer :: status(MPI_STATUS_SIZE)
        integer, allocatable, dimension(:) :: sendArryInt
        integer, allocatable, dimension(:) :: recvArryInt
        real(kind=dbPc), allocatable, dimension(:) :: sendArryDouble
        real(kind=dbPc), allocatable, dimension(:) :: recvArryDouble
        type(particleType) :: currentParticle
        type(particleType), dimension(maxSendRecvNum) :: tempSendRight
        type(particleType), allocatable, dimension(:) :: tempP

        allocate (tempP(pNum))
        sendRightNum = 0
        currentN = 0
        do n = 1, pNum
            currentParticle = particle(n)
            ip = currentParticle%indices(1)
            jp = currentParticle%indices(2)
            kp = currentParticle%indices(3)
            if (ip > currentNode%in) then
                sendRightNum = sendRightNum + 1
                tempSendRight(sendRightNum) = currentParticle
            elseif (currentParticle%altitude <= zMax) then
                currentN = currentN + 1
                tempP(currentN) = currentParticle
            end if
        end do

        call MPI_SENDRECV(sendRightNum, 1, MPI_I, currentNode%neighbor(2), 10, &
                          recvLeftNum, 1, MPI_I, currentNode%neighbor(1), 10, comm, status, ierr)
        if (sendRightNum > 0) then
            allocate (sendArryInt(3*sendRightNum))
            allocate (sendArryDouble(8*sendRightNum))
            do n = 1, sendRightNum
                currentParticle = tempSendRight(n)
                sendArryInt(n) = currentParticle%indices(1)
                sendArryInt(n + sendRightNum) = currentParticle%indices(2)
                sendArryInt(n + 2*sendRightNum) = currentParticle%indices(3)
                sendArryDouble(n) = currentParticle%location(1)
                sendArryDouble(n + sendRightNum) = currentParticle%location(2)
                sendArryDouble(n + 2*sendRightNum) = currentParticle%location(3)
                sendArryDouble(n + 3*sendRightNum) = currentParticle%velocity(1)
                sendArryDouble(n + 4*sendRightNum) = currentParticle%velocity(2)
                sendArryDouble(n + 5*sendRightNum) = currentParticle%velocity(3)
                sendArryDouble(n + 6*sendRightNum) = currentParticle%altitude
                sendArryDouble(n + 7*sendRightNum) = currentParticle%diameter
            end do
            call MPI_SEND(sendArryInt, sendRightNum*3, MPI_I, currentNode%neighbor(2), 20, comm, ierr)
            call MPI_SEND(sendArryDouble, sendRightNum*8, MPI_D, currentNode%neighbor(2), 20, comm, ierr)
            deallocate (sendArryInt)
            deallocate (sendArryDouble)
        end if
        call MPI_BARRIER(comm, ierr)
        if (recvLeftNum > 0) then
            pNum = currentN + recvLeftNum
            deallocate (particle)
            allocate (particle(pNum))
            particle(1:currentN) = tempP(1:currentN)
            deallocate (tempP)

            allocate (recvArryInt(3*recvLeftNum))
            allocate (recvArryDouble(8*recvLeftNum))
            call MPI_RECV(recvArryInt, recvLeftNum*3, MPI_I, currentNode%neighbor(1), 20, comm, status, ierr)
            call MPI_RECV(recvArryDouble, recvLeftNum*8, MPI_D, currentNode%neighbor(1), 20, comm, status, ierr)
            do n = 1, recvLeftNum
                currentN = currentN + 1
                currentParticle = particle(currentN)
                currentParticle%indices(1) = recvArryInt(n)
                currentParticle%indices(2) = recvArryInt(n + recvLeftNum)
                currentParticle%indices(3) = recvArryInt(n + 2*recvLeftNum)
                currentParticle%location(1) = recvArryDouble(n)
                currentParticle%location(2) = recvArryDouble(n + recvLeftNum)
                currentParticle%location(3) = recvArryDouble(n + 2*recvLeftNum)
                currentParticle%velocity(1) = recvArryDouble(n + 3*recvLeftNum)
                currentParticle%velocity(2) = recvArryDouble(n + 4*recvLeftNum)
                currentParticle%velocity(3) = recvArryDouble(n + 5*recvLeftNum)
                currentParticle%altitude = recvArryDouble(n + 6*recvLeftNum)
                currentParticle%diameter = recvArryDouble(n + 7*recvLeftNum)
                call determineParIJK(currentParticle)
                particle(currentN) = currentParticle
            end do
            deallocate (recvArryInt)
            deallocate (recvArryDouble)
        else
            pNum = currentN
            deallocate (particle)
            allocate (particle(pNum))
            particle = tempP(1:pNum)
            deallocate (tempP)
        end if
    end subroutine reallocateParticle

    ! *****************************
    ! Output particle data to files
    ! *****************************
    subroutine outputParticle(iter, t)
        use parallel_operations
        use vector_operations
        implicit none
        integer, intent(in) :: iter
        real(kind=dbPc), intent(in) :: t
        character(len=200) :: filename, line
        integer :: ierr, nameNum, n, amode, fh
        integer :: pNumTotal
        integer :: status(MPI_STATUS_SIZE)
        real(kind=dbPc) :: h, d, vMag
        real(kind=dbPc), dimension(3) :: loc
        real(kind=dbPc), dimension(3) :: vel

        if (mod(iter, intervalMonitor) == 0) then
            call MPI_ALLREDUCE(pNum, pNumTotal, 1, MPI_I, MPI_SUM, comm, ierr)
            call handleError(ierr)
            if (currentNode%ID == 0) then
                open (unit=10, position='append', file='./Particle/ParticleNum.plt')
                write (10, "(E15.2, I5)") t, pNumTotal
                close (10)
            end if
        end if

        nameNum = iter/intervalCreateFile
        write (filename, '("./Particle/ParticleData_", I0, ".plt")') nameNum
        if (mod(iter, intervalCreateFile) == 0) then
            if (currentNode%ID == 0) then
                open (unit=11, file=trim(adjustl(filename)))
                write (11, "(A80)") 'variables = "x", "y", "z", "h", "u", "v", "w", "velocity_mag", "d"'
                close (11)
            end if
        end if
        if (mod(iter, intervalParticle) == 0) then
            if (currentNode%ID == 0) then
                open (unit=11, position='append', file=filename)
                write (11, *) 'zone', ' T = "', t, '"'
                close (11)
            end if
            amode = MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_APPEND
            call MPI_FILE_OPEN(comm, trim(adjustl(filename)), amode, MPI_INFO_NULL, fh, ierr)
            do n = 1, pNum
                loc = particle(n)%location
                vel = particle(n)%velocity
                h = particle(n)%altitude
                d = particle(n)%diameter
                vMag = vectorMagnitude(vel)
                write (line, '(9E15.4)') loc(1), loc(2), loc(3), h, vel(1), vel(2), vel(3), vMag, d
                line = trim(adjustl(line))//char(10)
                call MPI_FILE_WRITE_SHARED(fh, trim(line), len(trim(line)), MPI_CHARACTER, status, ierr)
            end do
            call MPI_FILE_CLOSE(fh, ierr)
        end if
    end subroutine outputParticle

end module particle_operations

module output_operations
    use public_parameter
    implicit none
    include "mpif.h"

    private

    public :: generateOutPutFile

contains

    subroutine generateOutPutFile
        use parallel_operations
        implicit none
        character(len=32) :: bashCmd
        external :: system
        integer :: ierr

        if (currentNode%ID == 0) then
            bashCmd = 'rm -rf Particle'
            call system(trim(adjustl(bashCmd)))
            bashCmd = 'rm -rf Field'
            call system(trim(adjustl(bashCmd)))
            bashCmd = 'rm -rf Surface'
            call system(trim(adjustl(bashCmd)))

            bashCmd = 'mkdir Particle'
            call system(trim(adjustl(bashCmd)))
            bashCmd = 'mkdir Field'
            call system(trim(adjustl(bashCmd)))
            bashCmd = 'mkdir Surface'
            call system(trim(adjustl(bashCmd)))

            open (unit=10, file='./Particle/ParticleNum.plt')
            write (10, *) 'variables = "t", "Number"'
            close (10)

            open (unit=11, file='./Particle/ParticleData_0.plt')
            write (11, "(A80)") 'variables = "x", "y", "z", "h", "u", "v", "w", "velocity_mag", "d"'
            close (11)

            open (unit=12, file='./Field/Profile.plt')
            write (12, *) 'variables = "z", "u", "tau_p", "tau_f", "F_p", "phi_p"'
            close (12)

            open (unit=13, file='./Field/FieldData_0.plt')
            write (13, *) 'variables = "x", "y", "z", "Phi_p"'
            close (13)

            open (unit=14, file='./Surface/SurfaceData_0.plt')
            write (14, *) 'variables = "x", "y", "z", "dz", "d"'
            close (14)
        end if
        call MPI_BARRIER(comm, ierr)
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

    integer :: ierr
    integer :: iteration
    real(kind=dbPc) :: time

    ! find indices of subdomain and check that dimensions of arrays are sufficient
    if (mod(nx, nNodes) /= 0) then
        print *, 'nx cannot diveded by nNodes'
        stop
    end if
    call MPI_INIT(ierr)
    call initializeParallel
    call random_seed()
    call createMpiStructure
    ! generate surfGrid and initial bed
    call generateSurfaceGrid
    ! initialize surface
    call initializeSurface
    ! generate grid
    call generateGrid
    ! initialize fluid field
    call initializeFluidField
    ! initialize particle
    call initializeParticle
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
            call calculateParColl
            call calculateParticleMovement
            call reallocateParticle
            call calculateFluidField
            call updateSurfaceGrid
            call updateFieldGrid
            call outputParticle(iteration, time)
            call outputField(iteration, time)
            call outputSurface(iteration, time)
            ! **********************************Check1***********************************
            call MPI_BARRIER(comm, ierr)
            print *, currentNode%ID
            stop
            ! **************************************************************************
        end if
    end do
    call freeMpiStructure

    call MPI_FINALIZE(ierr)
end program main
