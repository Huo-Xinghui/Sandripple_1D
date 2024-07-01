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
    real(kind=dbPc), parameter :: repoAng = 30.0 ! repose angle
    real(kind=dbPc), parameter :: els = 0.9 ! normal restitution coefficient
    real(kind=dbPc), parameter :: fric = 0.0 ! tangential restitution coefficient
    real(kind=dbPc), parameter :: els1 = 0.9 ! normal restitution coefficient (mid-air collision)
    real(kind=dbPc), parameter :: fric1 = 0.0 ! tangential restitution coefficient (mid-air collision)
    real(kind=dbPc), parameter :: bedCellTknessInit = dpa*5.0 ! thickness of the thin surface layer on particle bed
    real(kind=dbPc), parameter :: rhos = 2650.0 ! particle density
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

module surface_operations
    use public_parameter
    implicit none
    private

    type surfaceGridType
        real(kind=dbPc), dimension(3) :: location
        real(kind=dbPc) :: averageDiameter
        real(kind=dbPc), dimension(npdf + 2) :: diameterDistribution
    end type surfaceGridType

    type(surfaceGridType), dimension(mx, my) :: surfGrid
    real(kind=dbPc), dimension(npdf + 2) :: initDiameterDist

    public :: generateSurfGrid, surfaceInitiation, initDiameterDist, surfGrid

contains

    subroutine generateSurfGrid
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
    end subroutine generateSurfGrid

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

        ip = int((particle%location(1) + xDiff)/xDiff) + 1
        jp = int((particle%location(2) + yDiff)/yDiff) + 1
        if (ip < 2) then
            ip = 1
        else if (ip > mx - 1) then
            ip = mx
        end if
        if (jp < 2) then
            jp = 1
        else if (jp > my - 1) then
            jp = my
        end if
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
    call generateSurfGrid
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
            call particleCal
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

