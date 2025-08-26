module public_val
    ! constants
    implicit none
    integer, parameter :: dbPc = 8 !selected_real_kind(15, 307)
    real(kind=dbPc), parameter :: pi = 3.14159265358979323846
    ! computational domain
    real(kind=dbPc), parameter :: xMax = 0.1 ! x size
    real(kind=dbPc), parameter :: yMax = 0.01 ! y size
    real(kind=dbPc), parameter :: zMax = 0.3 ! z size
    real(kind=dbPc), parameter :: area = xMax*yMax ! bottom area of the computational domain
    integer, parameter :: mx = 102 ! x grid num +2
    integer, parameter :: my = 12 ! y grid num +2
    integer, parameter :: mz = 150 ! z grid num
    integer, parameter :: mzUni = 30 ! z grid number above which zDiff becomes uniform
    real(kind=dbPc), parameter :: xDiff = xMax/(mx - 2) ! x grid size
    real(kind=dbPc), parameter :: yDiff = yMax/(my - 2) ! y grid size
    integer, parameter :: nNodes = 5 ! num of subdomain
    integer, parameter :: mxNode = (mx - 2)/nNodes + 2 ! x direction grid num for every subdomain
    ! particle
    logical, parameter :: midairCollision = .true. ! whether to consider midair collision
    ! particle diameter: whichDiameterDist = 0: normal distribution, 1: uniform diameter, 2: Bernoulli distribution, 3: lognormal distribution
    ! whichDiameterDist=0: npdf must >= 3, mu=dpa, sigma=dpStddDev
    ! whichDiameterDist=1: npdf must = 1, d=dpa
    ! whichDiameterDist=2: npdf must = 2, p1=prob1, p2=1-prob1, d1=dpa-dpStddDev, d2=dpa+dpStddDev
    ! whichDiameterDist=3: npdf must >= 3, mu=logMu, sigma=logSigma
    integer, parameter :: whichDiameterDist = 3 ! particle diameter distribution type
    integer, parameter :: npdf = 15 ! bin num of particle distribution
    integer, parameter :: pNumInit = 10 ! initial particle num
    integer, parameter :: maxEjectNum = 10000 ! max eject particle num in one time step
    integer, parameter :: maxNum = 100000 ! max particle num in one subdomain
    integer, parameter :: pNumInGridMax = maxNum/mxNode !/(my) ! max particle num in one x-y grid
    integer, parameter :: pNumExchMax = maxNum/10 ! max particle num for exchange between processors
    real(kind=dbPc), parameter :: dpa = 3.0e-4 ! average particle diameter
    real(kind=dbPc), parameter :: rhoP = 2650.0 ! particle density
    real(kind=dbPc), parameter :: dpStddDev = 2.0e-4 ! particle diameter standard deviation
    real(kind=dbPc), parameter :: logMu = -8.3221 ! mu of lognormal distribution
    real(kind=dbPc), parameter :: logSigma = 0.6167 ! sigma of lognormal distribution
    real(kind=dbPc), parameter :: prob1 = 0.5 ! probability one of Bernoulli distribution
    real(kind=dbPc), parameter :: binStart = 1.0e-4 ! start of the particle diameter distribution
    real(kind=dbPc), parameter :: binEnd = 10.0e-4 ! end of the particle diameter distribution
    real(kind=dbPc), parameter :: binWidth = (binEnd - binStart)/npdf ! bin width of the particle diameter distribution
    real(kind=dbPc), parameter :: resN = 0.7251 ! normal restitution coefficient
    real(kind=dbPc), parameter :: resT = -0.9288 ! tangential restitution coefficient
    real(kind=dbPc), parameter :: resN2 = 0.7 ! normal restitution coefficient for midair collision
    real(kind=dbPc), parameter :: por = 0.6 ! bedform porosity
    ! bed surface
    logical, parameter :: predefineSurface = .false. ! whether to predefine the bed surface
    real(kind=dbPc), parameter :: initSurfElevation = 4.0e-2 ! initial bed height
    real(kind=dbPc), parameter :: gammaE = 0.042 ! energy fraction that returns to bed surface
    real(kind=dbPc), parameter :: amp = 0.005 ! amplitude of the predefined surface
    real(kind=dbPc), parameter :: omg = 32.0*pi ! wave number (friquence) of the predefined surface
    real(kind=dbPc), parameter :: z0 = dpa/30.0 ! roughness height
    real(kind=dbPc), parameter :: repostAngle = 35.0 ! repost angle
    real(kind=dbPc), parameter :: tanRepostAngle = tan(repostAngle/180.0*pi) ! tan of the repost angle
    real(kind=dbPc), parameter :: chunkHeight = dpa ! height of particle diameter grid
    integer, parameter :: zChkNum = floor(initSurfElevation/chunkHeight)*2 ! num of sand bed chunks in z direction
    real(kind=dbPc), parameter :: chunkVol = xDiff*yDiff*chunkHeight ! volume of a particle diameter grid
    real(kind=dbPc), parameter :: blockHeight = dpa*3.0 ! height of the block, must>chunkHeight. block is a temporary grid on surface
    real(kind=dbPc), parameter :: blockVol = xDiff*yDiff*blockHeight ! volume of a block
    ! fluid
    real(kind=dbPc), parameter :: uStar = 0.30 ! fractional velocity
    real(kind=dbPc), parameter :: rho = 1.263 ! fluid density
    real(kind=dbPc), parameter :: nu = 1.51e-5 ! kinetic viscosity
    real(kind=dbPc), parameter :: kapa = 0.42 ! von Kaman's constant
    real(kind=dbPc), parameter :: zDiffMin = nu/uStar ! smallest z grid size
    real(kind=dbPc), parameter :: gHat = 9.81*(1.0 - rho/rhoP) ! reduced gravity
    ! iteration
    integer, parameter :: parStart = 1 ! iteration num when the particle calculation starts
    integer, parameter :: oneSecond = 20000 ! iteration num of one second
    integer, parameter :: intervalField = oneSecond*120 ! interval between fluid field outputs
    integer, parameter :: intervalProfile = oneSecond*60 ! interval between profile outputs
    integer, parameter :: intervalMonitor = oneSecond ! interval between monitor outputs
    integer, parameter :: intervalParticle = oneSecond*60 ! interval between particle outputs
    integer, parameter :: intervalSurface = oneSecond*30 ! interval between surface outputs
    integer, parameter :: intervalInside = oneSecond*120 ! interval between inside outputs
    integer, parameter :: intervalStatistics = oneSecond*60 ! interval between statistics outputs
    integer, parameter :: intervalCreateFile = oneSecond*240 ! interval between creating output files
    real(kind=dbPc), parameter :: dt = 5.0e-5 ! time step
    real(kind=dbPc), parameter :: endTime = 300.0 ! The time that the simulation lasts

    ! variables
    ! MPI
    integer :: realtype ! MPI_DOUBLE
    integer :: inttype ! MPI_INTEGER
    integer :: procs ! number of processors
    integer :: comm ! MPI communicator
    integer :: myID ! current processor ID
    integer :: sliceType ! MPI data type for a x-y slice
    integer :: edgeType1, edgeType2 ! MPI data type for a edgeline of x-y plan, edgeType1: for data(i, j), edgeType2: for data(i, j, k)
    integer, dimension(2) :: neighbor ! the index of the left and right neighbors
    ! Fluid field
    real(kind=dbPc), dimension(mxNode) :: x ! x coordinate
    real(kind=dbPc), dimension(my) :: y ! y coordinate
    real(kind=dbPc), dimension(mz) :: zDiff ! z grid size
    real(kind=dbPc), dimension(mz) :: z ! z coordinate
    real(kind=dbPc), dimension(mxNode, my, mz) :: zReal ! real z coordinate (z + zsf)
    real(kind=dbPc), dimension(mxNode) :: xu ! x coordinate of u
    real(kind=dbPc), dimension(my) :: yv ! y coordinate of v
    real(kind=dbPc), dimension(mz + 1) :: zw ! z coordinate of w
    real(kind=dbPc), dimension(mxNode, my, mz + 1) :: zwReal ! real z coordinate of w (zw + zsf)
    real(kind=dbPc), dimension(mxNode, my, mz) :: pfrac ! particle volume fraction
    real(kind=dbPc), dimension(mz) :: gridVolume ! grid volume
    real(kind=dbPc), dimension(mz) :: u ! u velocity
    real(kind=dbPc), dimension(mz) :: tau_p ! particle shear stress
    real(kind=dbPc), dimension(mz) :: tau_f ! fluid shear stress
    real(kind=dbPc), dimension(mz) :: phi_p ! particle volume fraction
    real(kind=dbPc), dimension(mz) :: F_p ! total particle drag force
    real(kind=dbPc), dimension(mz) :: F_pNode ! particle drag force per grid
    ! Particle bed
    integer, dimension(mxNode, my) :: rollTo, rollFrom ! the direction of particle rolling
    real(kind=dbPc), dimension(mxNode + 1) :: xsf ! x coordinate of the bed surface
    real(kind=dbPc), dimension(my) :: ysf ! y coordinate of the bed surface
    real(kind=dbPc), dimension(mxNode + 1, my) :: zsf ! z coordinate of the bed surface
    real(kind=dbPc), dimension(mxNode + 1, my) :: dsf ! particle diameter on the bed surface
    real(kind=dbPc), dimension(npdf) :: initDiameterDist ! initial particle diameter distribution
    real(kind=dbPc), dimension(mxNode + 1, my, npdf) :: hist ! particle diameter distribution
    real(kind=dbPc), dimension(mxNode + 1, my, npdf) :: binChange ! volume change of the particle diameter distribution bin
    real(kind=dbPc), dimension(zChkNum) :: zChk ! the z coordinate of the particle diameter grid center
    real(kind=dbPc), dimension(mxNode, my, zChkNum) :: chunkDia ! average diameter of a particle diameter grid
    real(kind=dbPc), dimension(mxNode, my, zChkNum, npdf) :: chunkHist ! diameter distribution of a particle diameter grid
    ! particle
    logical :: flag1
    logical :: flag2
    integer :: pNum ! particle number
    integer, dimension(maxNum, 3) :: pIndex ! the indices of the grid that the particle is in
    real(kind=dbPc) :: xflux, zflux ! the flux of the particle in x and z direction
    real(kind=dbPc), dimension(maxNum) :: xp, yp, zp ! particle position
    real(kind=dbPc), dimension(maxNum) :: up, vp, wp ! particle velocity
    real(kind=dbPc), dimension(maxNum) :: dp ! particle diameter
    real(kind=dbPc), dimension(maxNum) :: hp ! particle height
    real(kind=dbPc), dimension(maxNum) :: maxhp ! max particle height
    real(kind=dbPc), dimension(maxNum) :: survTime ! particle survival time
    real(kind=dbPc), dimension(maxNum) :: survLength ! particle survival length
    real(kind=dbPc), dimension(maxNum) :: collCount ! particle collision count
    real(kind=dbPc), dimension(mz) :: xfluxPf ! the x flux density of the particle
    real(kind=dbPc), dimension(mz) :: zfluxPf ! the z flux density of the particle
    ! iteration
    integer :: iter ! iteration number
    real(kind=dbPc) :: time ! current time
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
    public :: valObeyCertainPDF, generateNormalDistHistogram, &
              generateLogNormalDistHistogram, valObeyNormalDist, &
              valObeyPoissonDist, calculateD50, calculateD90

contains

    function valObeyCertainPDF(histogram) ! this function will normalize the input histogram
        implicit none
        real(kind=dbPc) :: valObeyCertainPDF
        real(kind=dbPc), dimension(npdf) :: histogram
        integer :: i
        real(kind=dbPc) :: rand, val, cumulativeProb
        !
        do i = 1, npdf
            histogram(i) = max(0.0, histogram(i))
        end do
        histogram = histogram/sum(histogram)
        cumulativeProb = 0.0
        val = dpa
        call random_number(rand)

        if (whichDiameterDist == 2) then
            if (rand <= histogram(1)) then
                val = binStart + 0.5*binWidth
            else
                val = binEnd - 0.5*binWidth
            end if
        else if (whichDiameterDist /= 1) then
            if (flag1) then
                do i = 1, npdf
                    cumulativeProb = cumulativeProb + histogram(i)
                    if (rand <= cumulativeProb) then
                        cumulativeProb = cumulativeProb - histogram(i)
                        val = binStart + real(i - 1, kind=dbPc)*binWidth + (rand - cumulativeProb)*binWidth/histogram(i)
                        exit
                    end if
                end do
                flag1 = .false.
            else
                do i = npdf, 1, -1
                    cumulativeProb = cumulativeProb + histogram(i)
                    if (rand <= cumulativeProb) then
                        cumulativeProb = cumulativeProb - histogram(i)
                        val = binStart + real(i, kind=dbPc)*binWidth - (rand - cumulativeProb)*binWidth/histogram(i)
                        exit
                    end if
                end do
                flag1 = .true.
            end if
        end if
        valObeyCertainPDF = val
    end function valObeyCertainPDF

    subroutine generateNormalDistHistogram(histogram) ! generate normal distribution histogram
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram
        real(kind=dbPc) :: total

        ! Initialize the histogram array to zero
        histogram = 0.0

        ! Fill the histogram array
        do i = 1, npdf
            histogram(i) = exp(-0.5*((binStart + (real(i, kind=dbPc) - 0.5)*binWidth - dpa)/dpStddDev)**2)/ &
                           (sqrt(2.0*pi)*dpStddDev)
            histogram(i) = histogram(i)*binWidth
        end do
        total = sum(histogram)
        histogram = histogram/total
    end subroutine generateNormalDistHistogram

    subroutine generateLogNormalDistHistogram(histogram) ! generate lognormal distribution histogram
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram
        real(kind=dbPc) :: total

        ! Initialize the histogram array to zero
        histogram = 0.0

        ! Fill the histogram array
        do i = 1, npdf
            histogram(i) = exp(-0.5*((log(binStart + (real(i, kind=dbPc) - 0.5)*binWidth) - logMu)/logSigma)**2)/ &
                           ((binStart + (real(i, kind=dbPc) - 0.5)*binWidth)*logSigma*sqrt(2.0*pi))
            histogram(i) = histogram(i)*binWidth
        end do
        total = sum(histogram)
        histogram = histogram/total
    end subroutine generateLogNormalDistHistogram

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

    function valObeyPoissonDist(lambda) result(k)
        implicit none
        real(kind=dbPc), intent(in) :: lambda
        integer :: k
        real(kind=dbPc) :: prob, rand, cumProb

        call random_number(rand)
        prob = exp(-lambda)
        cumProb = prob
        k = 0
        do while (rand > cumProb)
            k = k + 1
            prob = lambda*prob/real(k, kind=dbPc)
            cumProb = cumProb + prob
        end do
    end function valObeyPoissonDist

    function calculateD90(histogram) result(d90)
        ! calculate d90 of a certain diameter histogram
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram
        real(kind=dbPc) :: cumulativeProb, total, d90, res

        cumulativeProb = 0.0
        total = sum(histogram)
        do i = 1, npdf
            cumulativeProb = cumulativeProb + histogram(i)
            if (cumulativeProb >= 0.9*total) then
                res = cumulativeProb - 0.9*total
                d90 = binStart + real(i - 1, kind=dbPc)*binWidth + (res/histogram(i))*binWidth
                exit
            end if
        end do
    end function calculateD90

    function calculateD50(histogram) result(d50)
        ! calculate d50 of a certain diameter histogram
        implicit none
        integer :: i
        real(kind=dbPc), dimension(npdf) :: histogram
        real(kind=dbPc) :: cumulativeProb, total, d50, res

        cumulativeProb = 0.0
        total = sum(histogram)
        do i = 1, npdf
            cumulativeProb = cumulativeProb + histogram(i)
            if (cumulativeProb >= 0.5*total) then
                res = cumulativeProb - 0.5*total
                d50 = binStart + real(i - 1, kind=dbPc)*binWidth + (res/histogram(i))*binWidth
                exit
            end if
        end do

    end function calculateD50

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
    do while (time <= endTime .and. pNum > 0)
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
    integer :: i, j, k
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)

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
    do k = 1, zChkNum
        zChk(k) = real(k, kind=dbPc)*0.5*chunkHeight
    end do
end subroutine generateSurfGrid

subroutine initializeSurface
    use public_val
    use math_operations
    implicit none
    integer :: i, j, k, kk, n

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
    case (3)
        if (npdf >= 3) then
            call generateLogNormalDistHistogram(initDiameterDist)
        else
            print *, 'Error: npdf must be greater than or equal to 3 for lognormal distribution'
            stop
        end if
        !if (myID == 0) then
        !    print*, initDiameterDist
        !    stop
        !end if
    end select
    ! Initialize the particle bed
    zsf = initSurfElevation
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
    chunkDia = 0.0
    chunkHist = 0.0
    do j = 1, my
        do i = 1, mxNode
            kk = floor(zsf(i, j)/chunkHeight) + 1
            do k = 1, kk
                chunkDia(i, j, k) = dsf(i, j)
                do n = 1, npdf
                    chunkHist(i, j, k, n) = hist(i, j, n)
                end do
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
    real(kind=dbPc), dimension(npdf) :: currentHist

    if (binStart + 0.5*binWidth < 0.0) then
        print *, 'Error: particle diameter must be greater than or equal to 0'
        stop
    end if
    flag1 = .true.
    flag2 = .true.
    pNum = pNumInit
    currentHist = initDiameterDist
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
        if (whichDiameterDist == 1) then
            dp(n) = dpa
        else
            dp(n) = valObeyCertainPDF(currentHist)
        end if
        collCount(n) = 0.0
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
        bashCmd = 'rm -rf Inside'
        call system(trim(adjustl(bashCmd)))

        bashCmd = 'mkdir Particle'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Field'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Surface'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Statistics'
        call system(trim(adjustl(bashCmd)))
        bashCmd = 'mkdir Inside'
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

        open (unit=17, file='./Inside/InsideData_0.plt')
        write (17, *) 'variables = "x", "y", "z", "d"'
        close (17)
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
    integer :: n, nn, nnn
    integer :: ip, jp
    integer :: whichTri
    integer :: whichVer
    integer :: ipp, jpp
    logical :: rebound
    integer :: ii, jj
    integer :: iBin
    integer :: ejectNum, ejectNumMax
    integer :: nadd
    integer :: nxMax
    integer :: nAng
    real(kind=dbPc) :: estimateAltitude
    real(kind=dbPc) :: localXP, localYP
    real(kind=dbPc) :: d1, d2, d3, dTemp, zc, dCube, dRoll
    real(kind=dbPc) :: binLeft, binRight
    real(kind=dbPc) :: dd3, d13, d23, d13Temp, d23Temp
    real(kind=dbPc) :: sZeta, sZetaMax, inSq, cosz
    real(kind=dbPc) :: cosPsi, psi, angleC, xc
    real(kind=dbPc) :: m1, m2
    real(kind=dbPc) :: v1, v2, v2z, v2x
    real(kind=dbPc) :: eta, alpha, beta, lambda, sigma, mu
    real(kind=dbPc) :: e, eVx, eVz
    real(kind=dbPc) :: angin1
    real(kind=dbPc) :: cscAngin1, cotAngin1
    real(kind=dbPc) :: rhs, dxin, xinMaxMax, xinMaxMin, xinMaxTemp, xsin, xcos
    real(kind=dbPc) :: eq
    real(kind=dbPc) :: xinMin, xinMax, xin
    real(kind=dbPc) :: remx, hh, ht, wt
    real(kind=dbPc) :: rr1, rr2, rr3
    real(kind=dbPc) :: E1, E2, Ed1, Ed2, Eeff, Ec
    real(kind=dbPc) :: tau_s
    real(kind=dbPc) :: E2Bar
    real(kind=dbPc) :: merfc, erfc1, erfc2
    real(kind=dbPc) :: ejectVol, rollVol
    real(kind=dbPc) :: vch
    real(kind=dbPc), dimension(npdf) :: currentHist, currentDia
    real(kind=dbPc), dimension(maxEjectNum) :: tempx, tempy, tempz
    real(kind=dbPc), dimension(maxEjectNum) :: tempu, tempv, tempw
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
    real(kind=dbPc), dimension(60) :: angoutHist

    do n = 1, npdf
        currentDia(n) = binStart + (real(n, kind=dbPc) - 0.5)*binWidth
    end do
    nAng = size(angoutHist)
    nxMax = 100
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
            if (whichDiameterDist == 1) then
                d2 = dpa
                d3 = dpa
            else
                currentHist = hist(ipp, jpp, :)
                d2 = valObeyCertainPDF(currentHist)
                d3 = valObeyCertainPDF(currentHist)
            end if
            dd3 = d3/(0.5*(d1 + d2))
            d13 = 0.5*(d1 + d3)/(0.5*(d1 + d2))
            d23 = 0.5*(d2 + d3)/(0.5*(d1 + d2))
            !sZetaMax = dd3/(2.0*d23)
            !if (abs(sZetaMax) > 1.0) then
            !    sZetaMax = 1.0
            !end if
            !call random_number(rr1)
            !sZeta = sZetaMax * rr1
            !inSq = 1.0 - (sZeta*d23/d13)**2
            !if (inSq < 0.0) then
            !    inSq = 1.0
            !end if
            !d13Temp = d13*sqrt(inSq)
            !cosz = sqrt(1.0 - sZeta**2)
            !d23Temp = d23*cosz
            !do while (d23Temp - d13Temp >= 1.0)
            !    call random_number(rr1)
            !    sZeta = sZetaMax * rr1
            !    inSq = 1.0 - (sZeta*d23/d13)**2
            !    if (inSq < 0.0) then
            !        inSq = 1.0
            !    end if
            !    d13Temp = d13*sqrt(inSq)
            !    cosz = sqrt(1.0 - sZeta**2)
            !    d23Temp = d23*cosz
            !end do
            !d13 = d13Temp
            !d23 = d23Temp
            m1 = (pi*d1**3)/6.0*rhoP
            v1 = norm_2(pVol1)
            eta = resN*d1**3/(d1**3 + resN*d2**3)
            alpha = (1.0 + resN)/(1.0 + eta) - 1.0
            beta = 1.0 - (2.0/7.0)*(1.0 - resT)/(1.0 + eta)
            angin1 = atan(abs(vin(3)/vin(1)))
            cosPsi = (d13**2 + d23**2 - 1.0)/(2.0*d13*d23)
            psi = acos(cosPsi)
            angleC = 0.5*pi - psi
            xc = (1.0 + d23**2 - d13**2)/(2.0*d23)
            cscAngin1 = 1.0/sin(angin1)
            cotAngin1 = 1.0/tan(angin1)
            if (angin1 < angleC) then
                xinMin = d13*cscAngin1 - d23
            else
                xinMin = cotAngin1*sqrt(1.0 - xc**2) - xc
            end if
            rhs = 1.0/(1.0 + beta/alpha)
            xinMaxMax = 1.0/sin(angin1)
            xinMaxMin = 1.0/tan(angin1)
            dxin = (xinMaxMax - xinMaxMin)/real(nxMax, kind=dbPc)
            xinMax = xinMaxMin
            do nn = 1, nxMax
                xinMaxTemp = xinMaxMin + dxin*real(nn, kind=dbPc)
                xsin = xinMaxTemp*sin(angin1)
                xcos = xinMaxTemp*cos(angin1)
                eq = xsin**2 - xcos*sqrt(1.0 - xsin**2)
                if (eq > rhs) exit
                xinMax = xinMaxTemp
            end do
            if (xinMax < xinMin) then
                dTemp = d2
                d2 = d3
                d3 = dTemp
                d13 = 0.5*(d1 + d3)/(0.5*(d1 + d2))
                d23 = 0.5*(d2 + d3)/(0.5*(d1 + d2))
                !sZetaMax = dd3/(2.0*d23)
                !if (abs(sZetaMax) > 1.0) then
                !    sZetaMax = 1.0
                !end if
                !call random_number(rr1)
                !sZeta = sZetaMax * rr1
                !inSq = 1.0 - (sZeta*d23/d13)**2
                !if (inSq < 0.0) then
                !    inSq = 1.0
                !end if
                !d13Temp = d13*sqrt(inSq)
                !cosz = sqrt(1.0 - sZeta**2)
                !d23Temp = d23*cosz
                !do while (d23Temp - d13Temp >= 1.0)
                !    call random_number(rr1)
                !    sZeta = sZetaMax * rr1
                !    inSq = 1.0 - (sZeta*d23/d13)**2
                !    if (inSq < 0.0) then
                !        inSq = 1.0
                !    end if
                !    d13Temp = d13*sqrt(inSq)
                !    cosz = sqrt(1.0 - sZeta**2)
                !    d23Temp = d23*cosz
                !end do
                !d13 = d13Temp
                !d23 = d23Temp
                m1 = (pi*d1**3)/6.0*rhoP
                v1 = norm_2(pVol1)
                eta = resN*d1**3/(d1**3 + resN*d2**3)
                alpha = (1.0 + resN)/(1.0 + eta) - 1.0
                beta = 1.0 - (2.0/7.0)*(1.0 - resT)/(1.0 + eta)
                angin1 = atan(abs(vin(3)/vin(1)))
                cosPsi = (d13**2 + d23**2 - 1.0)/(2.0*d13*d23)
                psi = acos(cosPsi)
                angleC = 0.5*pi - psi
                xc = (1.0 + d23**2 - d13**2)/(2.0*d23)
                cscAngin1 = 1.0/sin(angin1)
                cotAngin1 = 1.0/tan(angin1)
                if (angin1 < angleC) then
                    xinMin = d13*cscAngin1 - d23
                else
                    xinMin = cotAngin1*sqrt(1.0 - xc**2) - xc
                end if
                rhs = 1.0/(1.0 + beta/alpha)
                xinMaxMax = 1.0/sin(angin1)
                xinMaxMin = 1.0/tan(angin1)
                dxin = (xinMaxMax - xinMaxMin)/real(nxMax, kind=dbPc)
                xinMax = xinMaxMin
            end if
            if (xinMax < xinMin) then
                xin = xinMax
            else
                call random_number(rr1)
                xin = xinMin + (xinMax - xinMin)*rr1
            end if
            xsin = xin*sin(angin1)
            xcos = xin*cos(angin1)
            eVx = -alpha*cos(angin1) + (alpha + beta)*xsin**2*cos(angin1) + &
                  (alpha + beta)*xsin*sin(angin1)*sqrt(1.0 - xsin**2)
            eVz = alpha*sin(angin1) - (alpha + beta)*xsin**2*sin(angin1) + &
                  (alpha + beta)*xsin*cos(angin1)*sqrt(1.0 - xsin**2)
            e = sqrt(eVx**2 + eVz**2)
            v2 = e*v1
            v2x = eVx*v1
            v2z = eVz*v1
            hh = (d1 + d3)*sin(psi)
            if (v2x > 0.0) then
                !ht = 0.5*(d1 + d2) - hh
                !wt = 0.5*(d1 + d2)*xc
                ht = 0.5*(d1 + d2) - 0.5*(d1+d3)*cos(angleC)
                wt = 0.5*(d1 + d3)*xc
                if (v2z**2/(2.0*gHat) > ht) then
                    remx = (sqrt((v2z/gHat)**2 - (2.0*ht)/gHat) + v2z/gHat)*v2x - wt
                else
                    remx = -1.0
                end if
            else if (v2x < 0.0) then
                !ht = 0.5*(d1 + d3) - hh
                !wt = 0.5*(d1 + d3)*cosPsi
                ht = 0.5*(d1 + d3)*(1.0 - cos(angleC))
                wt = 0.5*(d1 + d3)*sin(angleC)
                if (v2z**2/(2.0*gHat) > ht) then
                    remx = (sqrt((v2z/gHat)**2 - (2.0*ht)/gHat) + v2z/gHat)*(-v2x) - wt
                else
                    remx = -1.0
                end if
            else
                remx = -1.0
            end if
            if (remx <= 0.0 .or. e <= 0.0) then
                rebound = .false.
            else
                rebound = .true.
            end if
            if (rebound) then
                tempNum = tempNum + 1
                vout(1) = v2x
                vout(2) = 0.0
                vout(3) = v2z
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
                collCount(tempNum) = 0.0
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
                binChange(ii, jj, iBin) = binChange(ii, jj, iBin) + vch
            end if
            E1 = 0.5*m1*v1**2
            Ed1 = gammaE*(1.0 - e**2)*E1
            if (whichDiameterDist == 1) then
                dCube = dpa**3
                zc = dpa
            else if (whichDiameterDist == 0 .or. whichDiameterDist ==3) then
                dCube = 0.0
                !do nnn = 1, npdf
                !    binLeft = binStart + real(nnn - 1, kind=dbPc)*binWidth
                !    binRight = binStart + real(nnn, kind=dbPc)*binWidth
                !    dCube = dCube + (binLeft**2 + binRight**2)*(binLeft + binRight)*currentHist(nnn)/4.0
                !end do
                d2 = valObeyCertainPDF(currentHist)
                dCube = d2**3
                zc = calculateD90(currentHist)
            else if (whichDiameterDist == 2) then
                d2 = valObeyCertainPDF(currentHist)
                dCube = d2**3
                zc = binEnd - 0.5*binWidth
            end if
            m2 = (pi*dCube)/6.0*rhoP
            Ed2 = m2*gHat*zc
            !tau_s = rho*0.0123*(rhoP/rho*gHat*dpa + 3.0e-4/(rho*dpa))
            !Ec = Ed2*(1.0 - tau_f(1)/tau_s)
            !Eeff = max(Ec, 0.1*Ed2)
            Eeff = Ed2
            if (Ed1 > Eeff) then
                lambda = 2.0*log((1.0 - e**2)*E1/Eeff)
                sigma = sqrt(lambda)*log(2.0)
                mu = log((1.0 - e**2)*E1) - lambda*log(2.0)
                E2Bar = exp(mu + 0.5*sigma**2)
                merfc = erfc((log(Eeff) - mu)/(sqrt(2.0)*sigma))
                ejectNum = nint(Ed1/(2.0*E2Bar)*merfc)
            else
                ejectNum = 0
            end if
            !vch = (pi*d2**3)/6.0
            !ejectNumMax = blockVol/vch
            !ejectNum = min(ejectNum, ejectNumMax)
            if (ejectNum > 0) then
                erfc1 = erfc((log(Eeff) - mu - sigma**2)/(log(2.0)*sigma))
                erfc2 = erfc((log(Eeff) - mu)/(log(2.0)*sigma))
                E2 = E2Bar * erfc1/erfc2
                ejectVol = 0.0
                do nadd = 1, ejectNum
                    if (whichDiameterDist /= 1) then
                        d2 = valObeyCertainPDF(currentHist)
                    else
                        d2 = dpa
                    end if
                    vch = (pi*d2**3)/6.0
                    m2 = vch*rhoP
                    v2 = sqrt(2.0*E2/m2)
                    if (E2/(m2*gHat) < zc) then
                        cycle
                    end if
                    rr3 = -1.0
                    do while (rr3 < 0.0)
                        call random_number(rr1)
                        call random_number(rr2)
                        rr3 = 1.0 - rr1 - rr2
                    end do
                    nAddGlobal = nAddGlobal + 1
                    tempx(nAddGlobal) = rr1*point1(1) + rr2*point2(1) + rr3*point3(1)
                    tempy(nAddGlobal) = rr1*point1(2) + rr2*point2(2) + rr3*point3(2)
                    tempz(nAddGlobal) = rr1*point1(3) + rr2*point2(3) + rr3*point3(3)
                    tempu(nAddGlobal) = 0.0
                    tempv(nAddGlobal) = 0.0
                    tempw(nAddGlobal) = v2
                    tempd(nAddGlobal) = d2
                    tempi(nAddGlobal) = ip
                    tempj(nAddGlobal) = jp
                    zflux = zflux + m2/area/dt
                    ejectVol = ejectVol + vch
                    if (whichDiameterDist /= 1) then
                        iBin = floor((d2 - binStart)/binWidth) + 1
                        iBin = max(iBin, 1)
                        iBin = min(iBin, npdf)
                    else
                        iBin = 1
                    end if
                    binChange(ipp, jpp, iBin) = binChange(ipp, jpp, iBin) - vch
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
                    currentHist = hist(ii, jj, :)
                    do while (rollVol < ejectVol)
                        if (whichDiameterDist /= 1) then
                            dRoll = valObeyCertainPDF(currentHist)
                            iBin = floor((dRoll - binStart)/binWidth) + 1
                            iBin = max(iBin, 1)
                            iBin = min(iBin, npdf)
                        else
                            dRoll = dpa
                            iBin = 1
                        end if
                        vch = (pi*dRoll**3)/6.0
                        rollVol = rollVol + vch
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
    binChange = binChange/por
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
            up(tempNum + 1:tempNum + nAddGlobal) = tempu(1:nAddGlobal)
            vp(tempNum + 1:tempNum + nAddGlobal) = tempv(1:nAddGlobal)
            wp(tempNum + 1:tempNum + nAddGlobal) = tempw(1:nAddGlobal)
            dp(tempNum + 1:tempNum + nAddGlobal) = tempd(1:nAddGlobal)
            survLength(tempNum + 1:tempNum + nAddGlobal) = 0.0
            survTime(tempNum + 1:tempNum + nAddGlobal) = 0.0
            collCount(tempNum + 1:tempNum + nAddGlobal) = 0.0
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
        bulkForce(3) = -gHat
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
    integer :: status(MPI_STATUS_SIZE)
    real(kind=dbPc), dimension(pNumExchMax) :: xpE, ypE, zpE
    real(kind=dbPc), dimension(pNumExchMax) :: upE, vpE, wpE
    real(kind=dbPc), dimension(pNumExchMax) :: dpE, hpE, mhE
    real(kind=dbPc), dimension(pNumExchMax) :: slE, stE, ccE
    real(kind=dbPc), dimension(pNumExchMax) :: xpW, ypW, zpW
    real(kind=dbPc), dimension(pNumExchMax) :: upW, vpW, wpW
    real(kind=dbPc), dimension(pNumExchMax) :: dpW, hpW, mhW
    real(kind=dbPc), dimension(pNumExchMax) :: slW, stW, ccW
    real(kind=dbPc), allocatable, dimension(:) :: exchESend, exchWSend
    real(kind=dbPc), allocatable, dimension(:) :: exchERecv, exchWRecv

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
    allocate (exchESend(12*nESend))
    allocate (exchERecv(12*nERecv))
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
    exchESend(11*nESend + 1:12*nESend) = ccE(1:nESend)
    call MPI_SENDRECV(exchESend, nESend*12, realtype, neighbor(2), 24, &
                      exchERecv, nERecv*12, realtype, neighbor(1), 24, &
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
        collCount(pNum + 1:pNum + nERecv) = exchERecv(11*nERecv + 1:12*nERecv)
        pNum = pNum + nERecv
    end if
    deallocate (exchESend)
    deallocate (exchERecv)

    allocate (exchWSend(12*nWSend))
    allocate (exchWRecv(12*nWRecv))
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
    exchWSend(11*nWSend + 1:12*nWSend) = ccW(1:nWSend)
    call MPI_SENDRECV(exchWSend, nWSend*12, realtype, neighbor(1), 26, &
                      exchWRecv, nWRecv*12, realtype, neighbor(2), 26, &
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
        collCount(pNum + 1:pNum + nWRecv) = exchWRecv(11*nWRecv + 1:12*nWRecv)
        pNum = pNum + nWRecv
    end if
    deallocate (exchWSend)
    deallocate (exchWRecv)

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
        if (pNumInGrid(i, j) < pNumInGridMax) then
            pNumInGrid(i, j) = pNumInGrid(i, j) + 1
        else
            pNumInGrid(i, j) = pNumInGridMax
        end if
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
                    alpha1 = (1.0 + resN2)/(1.0 + eta1)
                    alpha2 = (1.0 + resN2)/(1.0 + eta2)
                    resTMidAir = 1.0 - 0.4*(1.0 + resN2)/(2.0/7.0)*rv12N/rv12T
                    resTMidAir = max(0.0, resTMidAir)
                    beta1 = (2.0/7.0)*(1.0 - resTMidAir)/(1.0 + eta1)
                    beta2 = (2.0/7.0)*(1.0 - resTMidAir)/(1.0 + eta2)
                    up(globalN1) = pVol1(1) - alpha1*rv12N*nVec(1) - beta1*(rv12(1) - rv12N*nVec(1))
                    vp(globalN1) = pVol1(2) - alpha1*rv12N*nVec(2) - beta1*(rv12(2) - rv12N*nVec(2))
                    wp(globalN1) = pVol1(3) - alpha1*rv12N*nVec(3) - beta1*(rv12(3) - rv12N*nVec(3))
                    up(globalN2) = pVol2(1) + alpha2*rv12N*nVec(1) + beta2*(rv12(1) - rv12N*nVec(1))
                    vp(globalN2) = pVol2(2) + alpha2*rv12N*nVec(2) + beta2*(rv12(2) - rv12N*nVec(2))
                    wp(globalN2) = pVol2(3) + alpha2*rv12N*nVec(3) + beta2*(rv12(3) - rv12N*nVec(3))
                    xp(globalN2) = xp(globalN2) + (contactDist - distance12)*nVec(1)
                    yp(globalN2) = yp(globalN2) + (contactDist - distance12)*nVec(2)
                    zp(globalN2) = zp(globalN2) + (contactDist - distance12)*nVec(3)
                    collCount(globalN1) = collCount(globalN1) + 1
                    collCount(globalN2) = collCount(globalN2) + 1
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
    real(kind=dbPc), dimension(my*npdf) :: binChangeESend, binChangeERecv
    real(kind=dbPc), dimension(my*npdf) :: binChangeWSend, binChangeWRecv
    real(kind=dbPc), dimension(my*npdf) :: binChangeEESend, binChangeEERecv

    ! because the value of ghost cell has changed
    ! need to add ghost value back to real domain before exchange
    ! binChange add back
    ! x=mxNode+1 add to x=3: send to 2 and receive from 1
    ! x=mxNode add to x=2: send to 2 and receive from 1
    ! x=1 add to x=mxNode-1: send to 1 and receive from 2
    do j = 1, my
        do k = 1, npdf
            jk = k + (j - 1)*npdf
            binChangeESend(jk) = binChange(mxNode, j, k)
            binChangeEESend(jk) = binChange(mxNode + 1, j, k)
            binChangeWSend(jk) = binChange(1, j, k)
        end do
    end do
    call MPI_SENDRECV(binChangeESend, my*npdf, realtype, neighbor(2), 15, &
                      binChangeERecv, my*npdf, realtype, neighbor(1), 15, comm, status, ierr)
    call MPI_SENDRECV(binChangeEESend, my*npdf, realtype, neighbor(2), 13, &
                      binChangeEERecv, my*npdf, realtype, neighbor(1), 13, comm, status, ierr)
    call MPI_SENDRECV(binChangeWSend, my*npdf, realtype, neighbor(1), 17, &
                      binChangeWRecv, my*npdf, realtype, neighbor(2), 17, comm, status, ierr)
    do j = 1, my
        do k = 1, npdf
            jk = k + (j - 1)*npdf
            binChange(2, j, k) = binChange(2, j, k) + binChangeERecv(jk)
            binChange(3, j, k) = binChange(3, j, k) + binChangeEERecv(jk)
            binChange(mxNode - 1, j, k) = binChange(mxNode - 1, j, k) + binChangeWRecv(jk)
        end do
    end do
    ! y=1 add to y=my-1, y=my add to y=2
    do i = 1, mxNode + 1
        do k = 1, npdf
            binChange(i, 2, k) = binChange(i, 2, k) + binChange(i, my, k)
            binChange(i, my - 1, k) = binChange(i, my - 1, k) + binChange(i, 1, k)
        end do
    end do
end subroutine addGhostData

subroutine updateSurfGrid
    use public_val
    implicit none

    integer :: i, j, k, ks, ks0, kb, n
    !real(kind=dbPc) :: floatKs
    real(kind=dbPc) :: blkBottom
    real(kind=dbPc) :: vChange, zChange
    real(kind=dbPc) :: bottomChkVol
    real(kind=dbPc) :: surfChkVolOld, surfChkVolNew
    real(kind=dbPc) :: zsfOld, zsfNew
    real(kind=dbPc) :: binOld
    !real(kind=dbPc) :: fixHist
    real(kind=dbPc), dimension(npdf) :: currentDia
    real(kind=dbPc), dimension(npdf) :: bin
    real(kind=dbPc), dimension(npdf) :: tempHist
    logical :: badflag
    !real(kind=dbPc) :: check, check1, check2

    call addGhostData
    if (whichDiameterDist /= 1) then
        do n = 1, npdf
            currentDia(n) = binStart + (real(n, kind=dbPc) - 0.5)*binWidth
        end do
        do j = 2, my - 1
            do i = 2, mxNode - 1
                vChange = sum(binChange(i, j, :))
                if (vChange == 0.0) cycle
                zChange = vChange/(xDiff*yDiff)
                zsfOld = zsf(i, j)
                ks0 = floor(zsfOld/chunkHeight) + 1
                !floatKs = real(ks0 - 1, kind=dbPc)
                surfChkVolOld = (zsfOld - (ks0 - 1)*chunkHeight)*xDiff*yDiff
                zsfNew = zsfOld + zChange
                badflag = .false.
                if (zsfNew - blockHeight < 0.0 .or. zsf(i, j) > zMax .or. abs(zChange) > 0.2*zsfOld) then
                    print *, 'Error: zsf reach the lower/upper boundary', zsfOld, zsfNew, zChange, dsf(i, j)
                    badflag = .true.
                    !stop
                end if
                !___________________________________________
                !zsf(i, j) = zsfNew
                !cycle ! ekalhxh: just change bed elevation!!
                !-------------------------------------------
                ks = floor(zsfNew/chunkHeight) + 1
                !floatKs = real(ks - 1, kind=dbPc)
                surfChkVolNew = (zsfNew - (ks - 1)*chunkHeight)*xDiff*yDiff

                if (zChange < 0.0) then
                    chunkDia(i, j, ks) = 0.0
                    do n = 1, npdf
                        bin(n) = chunkHist(i, j, ks0, n)*surfChkVolOld
                        do k = ks - 1, ks0 - 1
                            bin(n) = bin(n) + chunkHist(i, j, k, n)*chunkVol
                        end do
                        binOld = bin(n)
                        bin(n) = bin(n) + binChange(i, j, n)
                        bin(n) = max(bin(n), 0.1*binOld)
                        chunkHist(i, j, ks, n) = bin(n)/(surfChkVolNew + chunkVol)
                        chunkHist(i, j, ks - 1, n) = chunkHist(i, j, ks, n)
                        chunkDia(i, j, ks) = chunkDia(i, j, ks) + chunkHist(i, j, ks, n)*currentDia(n)
                    end do
                    chunkDia(i, j, ks - 1) = chunkDia(i, j, ks)
                    !fixHist = sum(chunkHist(i, j, ks, :)) - 1.0
                    !chunkHist(i, j, ks, :) = chunkHist(i, j, ks, :) - fixHist/npdf
                    !chunkHist(i, j, ks - 1, :) = chunkHist(i, j, ks - 1, :) - fixHist/npdf
                    if (ks0 > ks) then
                        do k = ks + 1, ks0
                            chunkDia(i, j, k) = 0.0
                            do n = 1, npdf
                                chunkHist(i, j, k, n) = 0.0
                            end do
                        end do
                    end if
                    !floatKs = real(ks - 2, kind=dbPc)
                    !zsf(i, j) = (ks - 2)*chunkHeight + sum(bin)/(xDiff*yDiff)
                    if (badflag) then
                        zsf(i, j) = zsfOld
                    else
                        zsf(i, j) = zsfNew
                    end if
                    !zsf(i, j) = zsfNew
                    !check = sum(bin) - surfChkVolNew - chunkVol
                    !if (abs(check) > 1.0e-10 .or. zsf(i, j) - blockHeight < 0.0) then
                    !    print *, 'Error1', zsf(i, j), zsfOld, zsfNew, zChange, check, fixHist/npdf, chunkDia(i, j, ks), vchange, chunkVol, surfChkVolNew
                    !    stop
                    !end if
                else
                    chunkDia(i, j, ks) = 0.0
                    do n = 1, npdf
                        bin(n) = chunkHist(i, j, ks0, n)*surfChkVolOld + chunkHist(i, j, ks0 - 1, n)*chunkVol
                        bin(n) = bin(n) + binChange(i, j, n)
                        chunkHist(i, j, ks, n) = bin(n)/(surfChkVolOld + chunkVol + vChange)
                        chunkDia(i, j, ks) = chunkDia(i, j, ks) + chunkHist(i, j, ks, n)*currentDia(n)
                    end do
                    !fixHist = sum(chunkHist(i, j, ks, :)) - 1.0
                    do k = ks0 - 1, ks - 1
                        do n = 1, npdf
                            chunkHist(i, j, k, n) = chunkHist(i, j, ks, n)
                        end do
                        chunkDia(i, j, k) = chunkDia(i, j, ks)
                    end do
                    !floatKs = real(ks0 - 2, kind=dbPc)
                    !zsf(i, j) = (ks0 - 2)*chunkHeight + sum(bin)/(xDiff*yDiff)
                    if (badflag) then
                        zsf(i, j) = zsfOld
                    else
                        zsf(i, j) = zsfNew
                    end if
                    !zsf(i, j) = zsfNew
                    !check = sum(bin) - surfChkVolOld - chunkVol - vChange
                    !if (abs(check) > 1.0e-10 .or. zsf(i, j) - blockHeight < 0.0) then
                    !    print *, 'Error2', zsf(i, j), zsfOld, zsfNew, zChange, sum(bin)/(xDiff*yDiff), check, fixHist/npdf, chunkDia(i, j, ks)
                    !    stop
                    !end if
                end if

                blkBottom = zsf(i, j) - blockHeight
                kb = floor(blkBottom/chunkHeight) + 1
                bottomChkVol = (kb*chunkHeight - blkBottom)*xDiff*yDiff
                do n = 1, npdf
                    bin(n) = surfChkVolNew*chunkHist(i, j, ks, n) + bottomChkVol*chunkHist(i, j, kb, n)
                    if (kb + 1 < ks) then
                        do k = kb + 1, ks - 1
                            bin(n) = bin(n) + chunkHist(i, j, k, n)*chunkVol
                        end do
                    end if
                    hist(i, j, n) = bin(n)/blockVol
                    tempHist(n) = max(hist(i, j, n), 0.0)
                end do
                tempHist = tempHist/sum(tempHist)
                dsf(i, j) = 0.0
                do n = 1, npdf
                    dsf(i, j) = dsf(i, j) + tempHist(n)*currentDia(n)
                end do
                if (dsf(i, j) > binEnd .or. dsf(i, j) < binStart) then
                    print *, 'Error: dsf out of range', dsf(i, j), hist(i, j, :)
                    stop
                end if
            end do
        end do
    else
        do j = 2, my - 1
            do i = 2, mxNode - 1
                vChange = sum(binChange(i, j, :))
                if (vChange == 0.0) cycle
                zsf(i, j) = zsf(i, j) + vChange/(xDiff*yDiff)
                if (zsf(i, j) < 0.0 .or. zsf(i, j) > zMax) then
                    print *, 'Error: zsf reach the lower/upper boundary', iter
                    stop
                end if
                do n = 1, npdf
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
    real(kind=dbPc), dimension(mx, my, zChkNum) :: chunkDiaTot

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

            open (unit=17, file='./Inside/InsideData_'//trim(adjustl(str))//'.plt')
            write (17, *) 'variables = "x", "y", "z", "d"'
            close (17)
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
            write (line, "(13E15.4)") xp(n), yp(n), zp(n), up(n), vp(n), wp(n), uMag, hp(n), maxhp(n), &
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

    if (mod(iter, intervalInside) == 0) then
        call gatherx(comm, mxNode, mx, x, xTot)
        call gatherxyz(comm, mxNode, mx, my, zChkNum, chunkDia, chunkDiaTot)
        if (myID == 0) then
            open (unit=17, position='append', file='./Inside/InsideData_'//trim(adjustl(str))//'.plt', action='write')
            write (17, *) 'zone', ' T = "', time, '"'
            write (17, *) 'i=', mx - 2, ' j=', my - 2, ' k=', zChkNum, ' datapacking=point'
            do k = 1, zChkNum
                do j = 2, my - 1
                    do i = 2, mx - 1
                        write (17, "(4E15.4)") xTot(i), y(j), zChk(k), chunkDiaTot(i, j, k)
                    end do
                end do
            end do
            close (17)
        end if
    end if

end subroutine outputAll
