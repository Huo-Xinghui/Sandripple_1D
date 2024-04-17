module public_val
  ! constants
  implicit none
  ! computational domain
  double precision, parameter :: xMax = 1.0 ! x size
  double precision, parameter :: yMax = 0.2 ! y size
  double precision, parameter :: zMax = 0.5 ! z size
  double precision, parameter :: zu = 0.1 ! the altitude that grid becomes parallel
  integer, parameter :: nx = 500 ! x grid num
  integer, parameter :: ny = 100 ! y grid num
  integer, parameter :: nz = 250 ! z grid num
  double precision, parameter :: zGridCtrl = 1.06 ! the control parameter of z grid
  integer, parameter :: nkx = 1000 ! x key point num
  integer, parameter :: nky = 200 ! y key point num
  integer, parameter :: nNodes = 5 ! num of subdomain
  ! time
  double precision, parameter :: dt = 1.0e-4 ! time step
  double precision, parameter :: tla = 600.0 ! time last
  ! fluid
  double precision, parameter :: wind = 0.5 ! fractional velocity
  double precision, parameter :: rho = 1.263 ! fluid density
  double precision, parameter :: nu = 1.49e-5 ! kinetic viscosity
  double precision, parameter :: kapa = 0.4 ! von Kaman's constant
  ! particle
  integer, parameter :: ikl = 1 ! calculating particles: ikl = 0: no, 1: yes
  integer, parameter :: nnps = 100 ! initial particle num
  double precision, parameter :: els = 0.9 ! normal restitution coefficient
  double precision, parameter :: fric = 0.0 ! tangential restitution coefficient
  double precision, parameter :: els1 = 0.9 ! normal restitution coefficient (mid-air collision)
  double precision, parameter :: fric1 = 0.0 ! tangential restitution coefficient (mid-air collision)
  double precision, parameter :: dpa = 2.5e-4 ! average particle diameter
  double precision, parameter :: dcgma = 2.0e-4 ! particle diameter standard deviation x 2
  integer, parameter :: npdf = 3 ! bin num of particle distribution. must be odd number when iud=0
  double precision, parameter :: szpdf = dpa*5.0 ! thickness of the thin surface layer on particle bed
  double precision, parameter :: rhos = 2650.0 ! particle density
  double precision, parameter :: nkl = 1.0 ! one particle stands for x particles
  double precision, parameter :: por = 0.6 ! bedform porosity
  integer, parameter :: nnpmax = 100000 ! max particle num in one subdomain
  integer, parameter :: nspmax = 10000 ! max eject particle num in one time step
  ! pre-formed surface
  double precision, parameter :: initAveKz = 0.05 ! initial average bed height
  double precision, parameter :: initAmp = 8.0*dpa ! amplitude of prerippled surface
  double precision, parameter :: initOmg = 4.0*pi ! wave number of prerippled surface
  double precision, parameter :: wavl = 2.0*pi/initOmg ! wavelength of prerippled surface
  ! method
  integer, parameter :: isp = 0 ! splash function: isp=0:lammel, isp=1:kok.
  ! boundary condition of particles: 0 periodic, 1 IO, 2 wall
  integer, parameter :: ikbx = 0 ! x direction
  integer, parameter :: ikby = 0 ! y direction
  integer, parameter :: ikbz = 1 ! uppper surface (never be 0)
  ! particle diameter: iud = 0: polydisperse, 1: monodisperse, 2: bidisperse
  ! if iud=0, normal distribution, mu=dpa, sigma=dcgma, dpa-3dcgma~dpa+3dcgma
  ! if iud=1, d=dpa.
  ! if iud=2, d1=dpa-dcgma, d2=dpa+dcgma, npdf must equal to 2.
  integer, parameter :: iud = 1
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
  integer, parameter :: pstart = 1 ! the initial step of particle calculation
  integer, parameter :: sstart = 1 ! the initial step of surface calculation
  integer, parameter :: pistart = 1 ! the initial step of particle info output
  ! file
  integer, parameter :: nnfi = 1e6 ! iter num contained in a file
  ! others
  integer, parameter :: mx = nx+2 ! x grid num +2
  integer, parameter :: my = ny+2 ! y grid num +2
  integer, parameter :: mz = nz+2 ! z grid num +2
  integer, parameter :: mkx = nkx+2 ! x key point num +2
  integer, parameter :: mky = nky+2 ! y key point num +2
  integer, parameter :: mxNode = nx/nNodes + 2 ! x grid num for every node
  integer, parameter :: mkxNode = nkx/nNodes + 2 ! x key point num for every node
  double precision, parameter :: pi = acos(-1.0) ! define Pi

  ! variables
  ! define in main
  integer :: realType
  integer :: intType
  integer :: comm1D
  integer :: myID
  integer, dimension(2) :: neighbor
  ! define in generateGrid
  double precision :: xDiff, yDiff
  double precision, dimension(mz) :: zDiff
  double precision, dimension(mxNode+1) :: xu
  double precision, dimension(my+1) :: yv
  double precision, dimension(mz+1) :: zw
  double precision, dimension(mxNode) :: x
  double precision, dimension(my) :: y
  double precision, dimension(mz) :: z
  ! define in dataExType
  integer :: fieldExType
  integer :: surfExType
  ! define in generateInitBed
  double precision :: kxDiff, kyDiff
  double precision, dimension(mkxNode) :: kx
  double precision, dimension(mky) :: ky
  double precision, dimension(mkxNode, mky) :: kz
  !
  double precision :: norm_vpin, norm_vpout, vvpin, vvpout, mpin, mpout
  double precision :: npin, npout
  double precision :: tot_nvpin, tot_nvpout, tot_vvpin, tot_vvpout, tot_mpin, tot_mpout
  double precision :: tot_npin, tot_npout
  double precision :: utaot, taot
  double precision :: aucreep
  double precision, dimension(3) :: vpin, vpout
  double precision, dimension(3) :: tot_vpin, tot_vpout
  double precision, dimension(mxNode, my, mz) :: phirho
  double precision, dimension(mz) :: ampu, ampd
  double precision, dimension(mz) :: fptx
  double precision, dimension(mz) :: totfptx, totvolpar
  integer :: nnp
  double precision :: uflx, wflx
  double precision, dimension(mz) :: uflxz
  double precision, dimension(mz) :: wflxz
  double precision, dimension(mkxNode, mky) :: dpaa
  double precision, dimension(mkxNode, mky) :: pnch
  double precision, dimension(mkxNode, mky, npdf) :: pdfch
  double precision, dimension(mky) :: eepnch, eepnchr
  double precision, dimension(mky*npdf) :: eepdfch, eepdfchr
  double precision :: time
  integer :: last
  integer :: gtypei
  double precision, dimension(mz) :: htaop, thtaop
  double precision, dimension(mz) :: htao, thtao
  double precision, dimension(mz) :: ahff, tahff
  double precision, dimension(mz) :: hru, thru
  double precision, dimension(mz) :: pcoll, ttpcoll
  double precision, dimension(npdf) :: spdf
  double precision, dimension(nnpmax) :: xp, yp, zp
  double precision, dimension(nnpmax) :: up, vp, wp
  double precision, dimension(nnpmax) :: dp, fk, fz
  double precision, dimension(nnpmax) :: fh, fg, ft
  double precision, dimension(mkxNode, mky, npdf) :: pdf
  double precision, dimension(mkxNode, mky) :: zpdf
  double precision, dimension(mkxNode, mky) :: ucreep, vcreep
end module public_val

! vector calculation
module vector_cal
  implicit none
contains
  function dot_prod(vec1, vec2)
    double precision :: dot_prod
    double precision, intent(in), dimension(3) :: vec1, vec2
    !
    dot_prod = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
  end function
  !
  function cross_prod(vec1, vec2)
    double precision, dimension(3) :: cross_prod
    double precision, intent(in), dimension(3) :: vec1, vec2
    !
    cross_prod(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
    cross_prod(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
    cross_prod(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  end function
  !
  function norm_2(vec)
    double precision :: norm_2
    double precision, intent(in), dimension(3) :: vec
    !
    norm_2 = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
  end function
  !
  function unit_vec(vec)
    double precision, dimension(3) :: unit_vec
    double precision, intent(in), dimension(3) :: vec
    double precision :: normv
    !
    normv = norm_2(vec)
    if (normv>0.0) then
      unit_vec = vec/normv 
    else
      unit_vec = 0.0
    end if
  end function
  !
  function dist_p(vec1, vec2)
    double precision :: dist_p
    double precision, intent(in), dimension(3) :: vec1, vec2
    !
    dist_p = sqrt((vec1(1)-vec2(1))**2+(vec1(2)-vec2(2))**2+(vec1(3)-vec2(3))**2)
  end function
end module vector_cal

! gather the matrix of different nodes into one single large matrix 
module gather_xyz
  implicit none
  interface gatherxyz
    module procedure gxyz_real
    module procedure gxyz_int
  end interface
contains
  subroutine gxyz_real(MPI_Comm, subMx, totMx, totMy, totMz, subField, totField)
    include "mpif.h"
    ! public
    integer, intent(in) :: subMx, totMx, totMy, totMz
    integer, intent(in) :: MPI_Comm
    double precision, intent(in), dimension(subMx, totMy, totMz) :: subField
    double precision, dimension(totMx, totMy, totMz) :: totField
    ! local
    integer :: subNx, totNx
    integer :: dataType
    integer :: subGridNum, totGridNum
    integer :: ierr
    integer :: i, j, k
    integer :: serialNum
    integer :: nodeID
    double precision, allocatable, dimension(:, :, :) :: subField_noGhost
    double precision, allocatable, dimension(:) :: subArray
    double precision, allocatable, dimension(:) :: totArray
    !
    subNx = subMx - 2
    totNx = totMx - 2
    subGridNum = subNx*totMy*totMz
    totGridNum = totNx*totMy*totMz
    allocate(subField_noGhost(subNx, totMy, totMz))
    allocate(subArray(subGridNum))
    allocate(totArray(totGridNum))
    !
    dataType = MPI_DOUBLE_PRECISION
    subArray = 0.0
    totArray = 0.0
    totField = 0.0
    subField_noGhost = subField(2:subMx-1, 1:totMy, 1:totMz)
    do k = 1, totMz
    do j = 1, totMy
    do i = 1, subNx
    serialNum = i + (j-1)*subNx + (k-1)*subNx*totMy
    subArray(serialNum) = subField_noGhost(i, j, k)
    end do
    end do
    end do
    call MPI_ALLGATHER(subArray,subGridNum,dataType,totArray,subGridNum,dataType,MPI_Comm,ierr)
    do k = 1, totMz
    do j = 1, totMy
    do i = 2, totMx-1
    nodeID = (i-2)/subNx
    serialNum = (i-1) - nodeID*subNx + (j-1)*subNx + (k-1)*subNx*totMy + nodeID*subNx*totMy*totMz
    totField(i, j, k) = totArray(serialNum)
    end do
    end do
    end do
    deallocate(subField_noGhost)
    deallocate(subArray)
    deallocate(totArray)
  end subroutine gxyz_real
  !
  subroutine gxyz_int(MPI_Comm, subMx, totMx, totMy, totMz, subField, totField)
    include "mpif.h"
    ! public
    integer, intent(in) :: subMx, totMx, totMy, totMz
    integer, intent(in) :: MPI_Comm
    integer, intent(in), dimension(subMx, totMy, totMz) :: subField
    integer, dimension(totMx, totMy, totMz) :: totField
    ! local
    integer :: subNx, totNx
    integer :: dataType
    integer :: subGridNum, totGridNum
    integer :: ierr
    integer :: i, j, k
    integer :: serialNum
    integer :: nodeID
    integer, allocatable, dimension(:, :, :) :: subField_noGhost
    integer, allocatable, dimension(:) :: subArray
    integer, allocatable, dimension(:) :: totArray
    !
    subNx = subMx - 2
    totNx = subMx - 2
    subGridNum = subNx*totMy*totMz
    totGridNum = totNx*totMy*totMz
    allocate(subField_noGhost(subNx, totMy, totMz))
    allocate(subArray(subGridNum))
    allocate(totArray(totGridNum))
    !
    dataType = MPI_INTEGER
    subArray = 0.0
    totArray = 0.0
    totField = 0.0
    subField_noGhost = subField(2:subMx-1, 1:totMy, 1:totMz)
    do k = 1, totMz
    do j = 1, totMy
    do i = 1, subNx
    serialNum = i + (j-1)*subNx + (k-1)*subNx*totMy
    subArray(serialNum) = subField_noGhost(i, j, k)
    end do
    end do
    end do
    call MPI_ALLGATHER(subArray,subGridNum,dataType,totArray,subGridNum,dataType,MPI_Comm,ierr)
    do k = 1, totMz
    do j = 1, totMy
    do i = 2, totMx-1
    nodeID = (i-2)/subNx
    serialNum = (i-1) - nodeID*subNx + (j-1)*subNx + (k-1)*subNx*totMy + nodeID*subNx*totMy*totMz
    totField(i, j, k) = totArray(serialNum)
    end do
    end do
    end do
    deallocate(subField_noGhost)
    deallocate(subArray)
    deallocate(totArray)
  end subroutine gxyz_int
  !
  subroutine gatherx(MPI_Comm, subMx, totMx, subField, totField)
    implicit none
    include "mpif.h"
    ! public
    integer, intent(in) :: subMx, totMx
    integer, intent(in) :: MPI_Comm
    double precision, intent(in), dimension(subMx) :: subField
    double precision, dimension(totMx) :: totField
    ! local
    integer :: dataType
    integer :: subGridNum, totGridNum
    integer :: ierr
    double precision, allocatable, dimension(:) :: subArray
    double precision, allocatable, dimension(:) :: totArray
    !
    subGridNum = subMx - 2
    totGridNum = totMx - 2
    allocate(subArray(subGridNum))
    allocate(totArray(totGridNum))
    !
    dataType = MPI_DOUBLE_PRECISION
    totArray = 0.0
    totField = 0.0
    subArray = subField(2:subMx-1)
    call MPI_ALLGATHER(subArray,subGridNum,dataType,totArray,subGridNum,dataType,MPI_Comm,ierr)
    totField(2:totMx-1) = totArray
    deallocate(subArray)
    deallocate(totArray)
  end subroutine gatherx
end module gather_xyz

program main
  use public_val
  implicit none
  include "mpif.h"
  integer :: ierr
  !integer :: nodes
  logical :: periodic
  integer :: nbrLeft, nbrRight
  integer :: i, j, k
  !
  call MPI_INIT(ierr)
  call random_seed()
  realType = MPI_DOUBLE_PRECISION
  intType = MPI_INTEGER
  ! create MPI Cartesian topology
  !nodes = nNodes ! May be this line can be deleted. By ekalhxh, in Apr 2nd 2024
  periodic = .true.
  call MPI_CART_CREATE(MPI_COMM_WORLD,1,nNodes,periodic,.true.,comm1D,ierr)
  call MPI_COMM_RANK(comm1D,myID,ierr)
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
  call MPI_CART_SHIFT(comm1D,0,1,nbrLeft,nbrRight,ierr)
  neighbor(1) = nbrLeft
  neighbor(2) = nbrRight
  call MPI_BARRIER(comm1D,ierr)
  ! check the number of nodes
  if (mod(nx, nNodes)/=0) then
    print*, 'x grid num cannot exactly be diveded by the num of nodes'
    stop
  end if
  if (mod(nkx, nNodes)/=0) then
    print*, 'x grid num of bed surface cannot exactly be diveded by the num of nodes'
    stop
  end if
  ! generate grid
  call generateGrid
  ! generate initial bed
  call generateInitBed
  ! define data exchange type for MPI
  call dataExType
  ! creat output file
  call outPutFile
  ! initialization
  fptx = 0.0
  pcoll = 0.0
  thtaop = 0.0
  thtao = 0.0
  tahff = 0.0
  thru = 0.0
  ttpcoll = 0.0
  nnp = 0
  uflx = 0.0
  wflx = 0.0
  uflxz = 0.0
  wflxz = 0.0
  dpaa = dpa
  time = 0.0
  last = 1
  phirho = 1.0
  tot_vpin = 0.0
  tot_vpout = 0.0
  tot_nvpin = 0.0
  tot_nvpout = 0.0
  tot_vvpin = 0.0
  tot_vvpout = 0.0
  tot_mpin = 0.0
  tot_mpout = 0.0
  tot_npin = 0.0
  tot_npout = 0.0
  totfptx = 0.0
  totvolpar = 0.0
  ucreep = 0.0
  vcreep = 0.0
  !
  ! start iteration loop
  !
  do
  ! calculate fluid field
  call fluidField
  ! generate boundary key point
  call imgd
  phirho = 1.0
  ! calculate particles
  if (ikl==1) then
    if (last==pstart) then
      call parstart
    end if
    if (last>=pstart) then
      call parcalculate
    end if
    if (last<sstart) then
      pnch = 0.0
      do i = 1, mkxNode
      do j = 1, mky
      zpdf(i, j) = szpdf
      if (irsf==0) then
        do k = 1, npdf
        pdf(i, j, k) = spdf(k)
        end do
        dpaa(i, j) = dpa
      else
        pdf(i, j, 2) = 0.5*(0.5*sin(initOmg*kx(i))+0.5)
        pdf(i, j, 1) = 1.0 - pdf(i, j, 2)
        dpaa(i, j) = pdf(i, j, 1)*(dpa-dcgma) + pdf(i, j, 2)*(dpa+dcgma)
      end if
      end do
      end do
    end if
  end if
  ! output result
  call output
  ! time advance
  time = time + dt
  last = last + 1
  if (time>tla) exit
  end do
  stop
  call MPI_FINALIZE(ierr)
end program main

subroutine generateGrid
  use public_val
  implicit none
  include "mpif.h"
  integer :: i, j, k
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  !
  ! define delta_x, delta_y, delta_z
  xDiff = xMax/dfloat(nx)
  yDiff = yMax/dfloat(ny)
  do k = 1, mz
  zDiff(k) = zGridCtrl**(k-1)*zMax/((zGridCtrl**mz-1.0)/(zGridCtrl-1.0))
  end do
  ! define the grid of vector
  xu = 0.0
  yv = 0.0
  zw = 0.0
  if (myID==0) then
    xu(1) = -xDiff
    do i = 2, mxNode+1
    xu(i) = xu(i-1) + xDiff
    end do
    call MPI_SEND(xu(mxNode-1),1,realType,neighbor(2),0,comm1D,ierr)
  else
    call MPI_RECV(xu(1),1,realType,neighbor(1),0,comm1D,status,ierr)
    do i = 2, mxNode+1
    xu(i) = xu(i-1) + xDiff
    end do
    call MPI_SEND(xu(mxNode-1),1,realType,neighbor(2),0,comm1D,ierr)
  end if
  yv(1) = -yDiff
  do j = 2, my+1
  yv(j) = yv(j-1) + yDiff
  end do
  zw(1) = 0.0
  do k = 2, mz+1
  zw(k) = zw(k-1) + zDiff(k-1)
  end do
  ! define the grid of scalar
  x = 0.0
  y = 0.0
  z = 0.0
  do i = 1, mxNode
  x(i) = (xu(i+1)-xu(i))/2.0 + xu(i)
  end do
  do j = 1, my
  y(j) = (yv(j+1)-yv(j))/2.0 + yv(j)
  end do
  do k = 1, mz
  z(k) = (zw(k+1)-zw(k))/2.0 + zw(k)
  end do
end subroutine generateGrid

subroutine generateInitBed
  use public_val
  implicit none
  include "mpif.h"
  ! local
  integer :: i, j
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  !
  ! define delta_kx, delta_ky, delta_kz
  kxDiff = xMax/dfloat(nkx)
  kyDiff = yMax/dfloat(nky)
  ! define the grid of surface key points
  kx = 0.0
  ky = 0.0
  kz = 0.0
  if (myID==0) then
      kx(1) = -kxDiff
      do i = 2, mkxNode
      kx(i) = kx(i-1) + kxDiff
      end do
      call MPI_SEND(kx(mkxNode-1),1,realType,neighbor(2),3,comm1D,ierr)
  else
      call MPI_RECV(kx(1),1,realType,neighbor(1),3,comm1D,status,ierr)
      do i = 2, mkxNode
      kx(i) = kx(i-1) + kxDiff
      end do
      call MPI_SEND(kx(mkxNode-1),1,realType,neighbor(2),3,comm1D,ierr)
  end if
  ky(1) = -kyDiff
  do j = 2, mky
  ky(j) = ky(j-1) + kyDiff
  end do
  do j = 1, mky
  do i = 1, mkxNode
  kz(i, j) = initAmp*sin(initOmg*kx(i)) + initAveKz
  end do
  end do
end subroutine generateInitBed

subroutine dataExType
  use public_val
  implicit none
  include "mpif.h"
  integer i, tempExType, ierr
  !
  ! fieldExType: i=const planes
  ! tempExType: i=const,k=const line 
  call MPI_TYPE_VECTOR(my,1,mxNode,realType,tempExType,ierr)
  call MPI_TYPE_COMMIT(tempExType,ierr)
  call MPI_TYPE_EXTENT(realType,i,ierr)
  call MPI_TYPE_HVECTOR(mz,1,mxNode*my*i,tempExType,fieldExType,ierr)
  call MPI_TYPE_COMMIT(fieldExType,ierr)
  ! surfExType: i=const line
  call MPI_TYPE_VECTOR(mky,1,mkxNode,realType,surfExType,ierr)
  call MPI_TYPE_COMMIT(surfExType,ierr)
end subroutine dataExType

subroutine outPutFile
  use public_val
  implicit none
  include "mpif.h"
  character(len=32) bashCmd
  !
  if (myID==0) then
    bashCmd = 'mkdir particle_loc'
    call system(trim(adjustl(bashCmd)))
    bashCmd = 'mkdir surface'
    call system(trim(adjustl(bashCmd)))
    bashCmd = 'mkdir field3D'
    call system(trim(adjustl(bashCmd)))
    !
    open (unit=32, file='./particle_loc/particle_loc0.plt')
    write (32, "(A82)") 'variables = "XP", "YP", "ZP", "DP", "UP", "VP", "WP", "FK", "FZ", "FH", "FG", "FT"'
    close(32)
    !
    open (unit=33, file='./surface/surface0.plt')
    write (33, *) 'variables = "X", "Y", "Z", "DP"'
    close(33)
    !
    open (unit=42, file='./field3D/field3D0.plt')
    write (42, *) 'variables = "X", "Y", "Z", "concentration"'
    close(42)
    !
    open (unit=31, file='particle_num.plt')
    write (31, *) 'variables = "T", "Num"'
    close(31)
    !
    open (unit=35, file='average_flux.plt')
    write (35, *) 'variables = "T", "uFlux", "wFlux", "salength"'
    close(35)
    !
    open (unit=36, file='flux_vs_height.plt')
    write (36, *) 'variables = "Z", "uFlux", "wFlux"'
    close(36)
    !
    open (unit=43, file='htao.plt')
    write (43, *) 'variables = "Z", "taoa", "taop", "vfrac", "u", "fptx"'
    close(43)
    !
    open (unit=39, file='vin.plt')
    write (39, *) 'variables = "T", "upin", "vpin", "wpin", "norm_vpin"'
    close(39)
    !
    open (unit=46, file='vout.plt')
    write (46, *) 'variables = "T", "upout", "vpout", "wpout", "norm_vpout"'
    close(46)
    !
    open (unit=44, file='eminout.plt')
    write (44, *) 'variables = "T", "vvpin", "vvpout", "mpin", "mpout"'
    close(44)
    !
    open (unit=45, file='numinout.plt')
    write (45, *) 'variables = "T", "numin", "numout"'
    close(45)
  end if
end subroutine outPutFile

subroutine fluidField
  use public_val
  implicit none
  include "mpif.h"
  integer :: i, j, k, n
  integer :: kk
  integer :: h
  integer :: ierr
  double precision :: mixl
  double precision :: dudz
  double precision :: oo
  double precision :: lmd
  double precision :: shru
  double precision :: ddz
  double precision :: chru
  double precision :: relax
  double precision, dimension(mz) :: volpar, numpar, tvolpar
  double precision, dimension(mz) :: afptx
  double precision, dimension(mz) :: pfptx
  double precision, dimension(mz) :: tfptx
  double precision, dimension(mz) :: ataop
  double precision, dimension(mz) :: ptaop
  double precision, dimension(mz) :: tampd, tampu, aampd, aampu
  ! function
  double precision :: ffd, ntmixl
  !
  !volpar = 0.0
  !do k = 1, mz
  !do j = 2, my-1
  !do i = 2, mxNode-1
  !volpar(k) = volpar(k) + (1.0-phirho(i, j, k))*xDiff*yDiff*zDiff(k)
  !end do
  !end do
  !end do
  !call MPI_ALLREDUCE(volpar,tvolpar,mz,realType,MPI_SUM,comm1D,ierr)
  !call MPI_ALLREDUCE(ampu,tampu,mz,realType,MPI_SUM,comm1D,ierr)
  !afptx = tampu/xMax/yMax
  !totfptx = totfptx + afptx/dfloat(nna)
  !totvolpar = totvolpar + tvolpar/dfloat(nna)
    htao = 0.0
    ahff = 1.0
    do k = 1, mz
    ahff(k) = 1.0 - totvolpar(k)/(xMax*yMax*zDiff(k))
    end do
    !
    fptx = totfptx/ahff
    !
    relax = 0.01
    do k = 1, mz
    tfptx(k) = sum(fptx(k:mz))
    htaop(k) = htaop(k)*(1.0-relax) + relax*tfptx(k)
    if (ahff(k)<=(1.0-por)) then
      htaop(k) = rho*wind**2
    end if
    htao(k) = rho*wind**2 - htaop(k)
    if (htao(k)<0.0) htao(k) = 0.0
    !if (htao(k)>rho*wind**2) htao(k) = rho*wind**2
    end do
    !
    mixl = kapa*z(1)*(1.0-exp(-1.0/26.0*z(1)*wind/nu))
    dudz = (sqrt(nu**2+4.0*mixl**2*htao(1)/rho)-nu)/2.0/mixl**2
    hru(1) = dudz*z(1)
    do k = 2, mz
    ddz = z(k) - z(k-1)
    mixl = kapa*z(k)*(1.0-exp(-1.0/26.0*z(k)*wind/nu))
    dudz = (sqrt(nu**2+4.0*mixl**2*htao(k)/rho)-nu)/2.0/mixl**2
    hru(k) = dudz*ddz + hru(k-1)
    end do
    2001 continue
    !
    totfptx = 0.0
    totvolpar = 0.0
end subroutine fluidField

subroutine imgd
  use public_val
  implicit none
  include "mpif.h"
  ! local
  integer :: i, j
  integer :: n
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  integer :: wn
  double precision :: aaa
  double precision :: totpz, avepz, tavepz
  double precision :: posit
  !
  if (last==1 .or. ikl==0) then
    kx = 0.
    ky = 0.
    kz = 0.
    pnch = 0.
    kxDiff = xMax/dfloat(mkx-2)
    kyDiff = yMax/dfloat(mky-2)
    if (myID==0) then
      kx(1) = -kxDiff
      do i = 2, mkxNode
      kx(i) = kx(i-1) + kxDiff
      end do
      call MPI_SEND(kx(mkxNode-1),1,realType,neighbor(2),3,comm1D,ierr)
    else
      call MPI_RECV(kx(1),1,realType,neighbor(1),3,comm1D,status,ierr)
      do i = 2, mkxNode
      kx(i) = kx(i-1) + kxDiff
      end do
      call MPI_SEND(kx(mkxNode-1),1,realType,neighbor(2),3,comm1D,ierr)
    end if
    ky(1) = -kyDiff
    do j = 2, mky
    ky(j) = ky(j-1) + kyDiff
    end do
    do j = 1, mky
    do i = 1, mkxNode
    kz(i, j) = initAmp*sin(dfloat(int(time/20.0)+2)*8.0*pi*kx(i)) + initAveKz
    !wn = int(kx(i)/wavl)
    !posit = kx(i)/wavl - dfloat(wn) - 0.5
    !if (posit>=0.0) then
    !  kz(i, j) = initAmp*0.5 - 2.0*initAmp*posit + initAveKz
    !else
    !  kz(i, j) = initAmp*0.5 + 2.0*initAmp*posit + initAveKz
    !end if
    end do
    end do
  else
    aaa = kxDiff*kyDiff
    !do j = 1, mky
    !do i = 1, mkxNode
    !kz(i, j) = kz(i, j) + pnch(i, j)/aaa/por
    !end do
    !end do
    do j = 1, mky
    do i = 1, mkxNode
    kz(i, j) = initAmp*sin(dfloat(int(time/20.0)+2)*8.0*pi*kx(i)) + initAveKz
    end do
    end do
    !do j = 1, mky
    !do i = 1, mkxNode
    !if (kz(i, j)<0. .or. kz(i, j)>zMax .or. abs(pnch(i, j)/aaa/por)>=0.01) then
    ! print*, 'error: kz out of lower/upper boundary'
    !  print*, 'z=', kz(i, j), '    z change=', pnch(i, j)/aaa/por
    !  print*, 'i=', i, '  j=', j, '  myID=', myID
    !  kz(i, j) = kz(i, j) - pnch(i, j)/aaa/por
    !end if
    !end do
    !end do
    call pxch(mkxNode, mky, kz, surfExType, neighbor, comm1D)
    do i = 1, mkxNode
    kz(i, 1) = kz(i, mky-1)
    kz(i, mky) = kz(i, 2)
    end do
    !
    n = 0
    totpz = 0.0
    do j = 2, mky-1
    do i = 2, mkxNode-1
    n = n + 1
    totpz = totpz + kz(i, j)
    end do
    end do
    avepz = totpz/dfloat(n)
    call MPI_ALLREDUCE(avepz,tavepz,1,realType,MPI_SUM,comm1D,ierr)
  end if
end subroutine imgd

subroutine pxch(mkxNode, mky, kz, surfExType, neighbor, comm1D)
  implicit none
  include "mpif.h"
  ! public
  integer :: mkxNode, mky
  integer :: comm1D
  integer :: surfExType
  integer, dimension(2) :: neighbor
  double precision, dimension(mkxNode, mky) :: kz
  ! local
  integer :: i, j
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  !
  ! send to 2 and receive from 1
  call MPI_SENDRECV(kz(mkxNode-1, 1),1,surfExType,neighbor(2),1,kz(1, 1),1,surfExType,neighbor(1),1,comm1D,status,ierr)
  ! send to 1 and receive from 2
  call MPI_SENDRECV(kz(2, 1),1,surfExType,neighbor(1),2,kz(mkxNode, 1),1,surfExType,neighbor(2),2,comm1D,status,ierr)
end subroutine pxch

function ntmixl(lm, k, nu, ux, dz)
  implicit none
  ! public
  double precision :: ntmixl
  double precision, intent(in) :: lm, k, nu, ux, dz
  ! local
  integer :: n
  double precision :: xr, x1
  double precision :: fx, gx, fxr
  x1 = 1.0
  n = 0
  do
  n = n + 1
  if (n>10000) then
    print*, 'ntmixl', lm, k, nu, ux, dz
    stop
  end if
  fx = k*(1.0-exp(-sqrt(ux*x1/7.0/nu))) - (x1-lm)/dz
  gx = exp(-sqrt(ux*x1/7.0/nu))*k*ux/(2.0*nu*7.0*sqrt(ux*x1/7.0/nu)) - 1.0/dz
  xr = x1 - fx/gx
  fxr = k*(1.0-exp(-sqrt(ux*xr/7.0/nu))) - (xr-lm)/dz
  if (abs(fxr)<1.0e-6) exit
  x1 = xr
  end do
  ntmixl = xr
end function ntmixl

subroutine parstart
  use public_val
  implicit none
  include "mpif.h"
  integer :: n
  integer :: bnldev
  integer :: i, j, k
  integer :: nbi
  double precision :: rr1, rr2, rr3
  double precision :: parsiz
  double precision :: xx
  double precision :: tspdf
  double precision :: ppdf
  double precision :: cgma
  double precision :: mu
  !
  if (iud==0) then
    cgma = dcgma*1.0e4
    mu = dpa*1.0e4
    do i = 1, npdf
    xx = (dfloat(i-1)+0.5)*cgma*6.0/dfloat(npdf) + (mu-cgma*3.0)
    spdf(i) = 1.0/(sqrt(2.0*pi)*cgma)*exp(-(xx-mu)**2/(2.0*cgma**2))
    tspdf = tspdf + spdf(i)
    end do
    spdf = spdf/tspdf
  else if (iud==2) then
    spdf(1) = 0.5
    spdf(2) = 0.5
  end if
  do n = 1, nnps
  call random_number(rr1)
  call random_number(rr2)
  call random_number(rr3)
  xp(nnp+n) = (xu(mxNode)-xu(2))*rr1 + xu(2)
  yp(nnp+n) = (yv(my)-yv(2))*rr2 + yv(2)
  zp(nnp+n) = 0.05*rr3 + initAveKz
  up(nnp+n) = 0.0
  vp(nnp+n) = 0.0
  wp(nnp+n) = 0.0
  fk(nnp+n) = 0.0
  fh(nnp+n) = 0.0
  fg(nnp+n) = 0.0
  ft(nnp+n) = 0.0
  fz(nnp+n) = 0.0
  if (iud==0) then
    dp(nnp+n) = parsiz(spdf, dpa, dcgma, npdf)
  else if (iud==1) then
    dp(nnp+n) = dpa
  else if (iud==2) then
    ppdf = spdf(1)
    nbi = bnldev(ppdf, 1)
    if (nbi==1) then
      dp(nnp+n) = dpa - dcgma
    else
      dp(nnp+n) = dpa + dcgma
    end if
  end if
  end do
  nnp = nnp + nnps
  do j = 1, mky
  do i = 1, mkxNode
  zpdf(i, j) = szpdf
  dpaa(i, j) = dpa
  do k = 1, npdf
  pdf(i, j, k) = spdf(k)
  end do
  end do
  end do
end subroutine parstart

subroutine parcalculate
  use public_val
  use vector_cal
  implicit none
  include "mpif.h"
  integer :: nn, nni
  integer :: n
  integer :: i, j, k
  integer :: ierr
  integer :: ip, jp
  integer :: ne
  integer :: nrol
  integer :: kd
  integer :: nks, tnk, nnks
  integer :: bnldev
  integer :: nnn, nnni
  integer :: nnnp
  integer :: addnnp
  integer :: nksmax, chmax
  integer :: ii, jj, kk
  integer :: ipp, jpp
  integer :: jk
  integer :: iii, jjj
  integer :: nne
  integer :: nbi
  integer :: kkk
  integer :: hpl
  integer :: status(MPI_STATUS_SIZE)
  integer :: nkk(mxNode, my)
  double precision :: mp, volp
  double precision :: lmt1, lmt2, lmt3
  double precision :: norm_vin
  double precision :: mmu
  double precision :: cgma
  double precision :: alpha, beta
  double precision :: alpha1, beta1
  double precision :: alpha2, beta2
  double precision :: eta1, eta2
  double precision :: rr1, rr2, rr3
  double precision :: normal
  double precision :: rrr, rrrx, rrrd
  double precision :: nv1, nv2, tv1, tv2
  double precision :: pp
  double precision :: lambda
  double precision :: ammu1, ammu2
  double precision :: gg1, gg2, gg3
  double precision :: dzz
  double precision :: acgma1, acgma2
  double precision :: angout1, angout2
  double precision :: norm_vout
  double precision :: expdev
  double precision :: ppz0, ppz1, ppz2, ppz3, ppz4
  double precision :: tam
  double precision :: dd1, dd2
  double precision :: angin1
  double precision :: ee1, ee2, eed2, eed2x
  double precision :: mm1, mm2
  double precision :: d1, d2
  double precision :: myerfc
  double precision :: merfc
  double precision :: ee2bar
  double precision :: parsiz
  double precision :: prebound
  double precision :: arebound
  double precision :: nee
  double precision :: vch
  double precision :: dpd
  double precision :: svlayer, xvlayer
  double precision :: norm_nnvec, norm_ttvec, norm_tnvec
  double precision :: d01, d02, d03
  double precision :: mmp, nump
  double precision :: hrlx1, dhcld
  double precision :: xp0
  integer, dimension(mkxNode, mky) :: rolldir1, rolldir2
  double precision, dimension(npdf) :: ppdf
  double precision, dimension(mky) :: eepz, eepzr
  double precision, dimension(mky*npdf) :: ppdfch, ppdfchr
  double precision, dimension(mky*npdf) :: ipdfch, ipdfchr
  double precision, dimension(mky) :: pzpdfch, pzpdfchr
  double precision, dimension(mky) :: izpdfch, izpdfchr
  double precision, dimension(mkxNode, mky) :: vlayer
  double precision, dimension(mkxNode, mky, npdf) :: vbin
  double precision, dimension(mkxNode, mky) :: tpdf
  double precision, dimension(mxNode, my, mz) :: vpar
  double precision, dimension(mz) :: ncoll, ntotc
  double precision, dimension(nspmax) :: tempx, tempy, tempz, tempu, tempv, tempw, tempd, tempfk, tempfz, &
    tempfh, tempfg, tempft
  double precision, dimension(mkxNode, mky) :: mupin, mvpin, mupout, mvpout, qcreepx, qcreepy
  double precision, dimension(4) :: ta
  double precision, dimension(3) :: point0, point1, point2, point3
  double precision, dimension(3) :: vec1, vec2, vec3
  double precision, dimension(3) :: nnvec, ttvec, tnvec
  double precision, dimension(3) :: unnvec, uttvec, utnvec
  double precision, dimension(3) :: upvec1, upvec2, xpvec1, xpvec2
  double precision, dimension(3) :: vin
  double precision, dimension(3) :: vout
  double precision, dimension(3) :: gg
  integer, allocatable, dimension(:, :, :) :: scp
  double precision, allocatable, dimension(:) :: chxp, chyp, chzp
  double precision, allocatable, dimension(:) :: chup, chvp, chwp
  double precision, allocatable, dimension(:) :: chdp, chfk, chfz
  double precision, allocatable, dimension(:) :: chfh, chfg, chft
  double precision, allocatable, dimension(:) :: chxpi, chypi, chzpi
  double precision, allocatable, dimension(:) :: chupi, chvpi, chwpi
  double precision, allocatable, dimension(:) :: chdpi, chfki, chfzi
  double precision, allocatable, dimension(:) :: chfhi, chfgi, chfti
  double precision, allocatable, dimension(:) :: exch, exchi
  double precision, allocatable, dimension(:) :: exchr, exchir
  nksmax = nnpmax !/(mxNode) !/(my)
  chmax = nnpmax/10
  allocate(scp(mxNode, my, nksmax))
  allocate(chxp(chmax), chyp(chmax), chzp(chmax))
  allocate(chup(chmax), chvp(chmax), chwp(chmax))
  allocate(chdp(chmax), chfk(chmax), chfz(chmax))
  allocate(chfh(chmax), chfg(chmax), chft(chmax))
  allocate(chxpi(chmax), chypi(chmax), chzpi(chmax))
  allocate(chupi(chmax), chvpi(chmax), chwpi(chmax))
  allocate(chdpi(chmax), chfki(chmax), chfzi(chmax))
  allocate(chfhi(chmax), chfgi(chmax), chfti(chmax))
  !
  ! influence of repose angle and particle roll direction
  eepz(:) = kz(3, :)
  call MPI_SENDRECV(eepz,mky,realType,neighbor(1),105, &
    eepzr,mky,realType,neighbor(2),105,comm1D,status,ierr)
  !
  rolldir1 = 5
  rolldir2 = 5
  do jp = 2, mky
  do ip = 2, mkxNode
  ppz0 = kz(ip, jp)
  if (ip<mkxNode) then
    ppz1 = kz(ip+1, jp)
  else
    ppz1 = eepzr(jp)
  end if
  ppz2 = kz(ip-1, jp)
  if (jp<mky) then
    ppz3 = kz(ip, jp+1)
  else
    ppz3 = kz(ip, 3)
  end if
  ppz4 = kz(ip, jp-1)
  !
  ta(1) = (ppz0-ppz1)/kxDiff
  ta(2) = (ppz0-ppz2)/kxDiff
  ta(3) = (ppz0-ppz3)/kyDiff
  ta(4) = (ppz0-ppz4)/kyDiff
  tam = ta(1)
  rolldir1(ip, jp) = 1
  do k = 2, 4
  if (ta(k)>tam) then
    tam = ta(k)
    rolldir1(ip, jp) = k
  end if
  end do
  if (tam<=tan(30.0/180.0*pi)) rolldir1(ip, jp) = 0
  !
  ta(1) = -(ppz0-ppz1)/kxDiff
  ta(2) = -(ppz0-ppz2)/kxDiff
  ta(3) = -(ppz0-ppz3)/kyDiff
  ta(4) = -(ppz0-ppz4)/kyDiff
  tam = ta(1)
  rolldir2(ip, jp) = 1
  do k = 2, 4
  if (ta(k)>tam) then
    tam = ta(k)
    rolldir2(ip, jp) = k
  end if
  end do
  if (tam<=tan(30.0/180.0*pi)) rolldir2(ip, jp) = 0
  end do
  end do
  ! rebound and splash
  nnnp = 0
  addnnp = 0
  pnch = 0.0
  pdfch = 0.0
  eepnch = 0.0
  eepdfch = 0.0
  vpar = 0.0
  ampu = 0.0
  ampd = 0.0
  vpin = 0.0
  vpout = 0.0
  norm_vpin = 0.0
  norm_vpout = 0.0
  vvpin = 0.0
  vvpout = 0.0
  mpin = 0.0
  mpout = 0.0
  mupin = 0.0
  mvpin = 0.0
  mupout = 0.0
  mvpout = 0.0
  npin = 0.0
  npout = 0.0
  dospl: do n = 1, nnp
  ip = int((xp(n)-kx(1))/kxDiff) + 1
  jp = int((yp(n)-ky(1))/kyDiff) + 1
  if (ip<2) ip = 2
  if (jp<2) jp = 2
  if (ip>mkxNode-1) ip = mkxNode-1
  if (jp>mky-1) jp = mky-1
  if (mod(ip+jp, 2)==0) then
    if (abs(yp(n)-ky(jp))>=kyDiff/kxDiff*abs(xp(n)-kx(ip))) then
      kk = 1
      point1(1) = kx(ip)
      point1(2) = ky(jp)
      point1(3) = kz(ip, jp)
      point2(1) = kx(ip+1)
      point2(2) = ky(jp+1)
      point2(3) = kz(ip+1, jp+1)
      point3(1) = kx(ip)
      point3(2) = ky(jp+1)
      point3(3) = kz(ip, jp+1)
      vec1 = point2 - point1
      vec2 = point3 - point1
      nnvec = cross_prod(vec1, vec2)
    else 
      kk = 2
      point1(1) = kx(ip)
      point1(2) = ky(jp)
      point1(3) = kz(ip, jp)
      point2(1) = kx(ip+1)
      point2(2) = ky(jp+1)
      point2(3) = kz(ip+1, jp+1)
      point3(1) = kx(ip+1)
      point3(2) = ky(jp)
      point3(3) = kz(ip+1, jp)
      vec1 = point2 - point1
      vec2 = point3 - point1
      nnvec = cross_prod(vec2, vec1)
    end if
  else
    if (abs(yp(n)-ky(jp))<=-kyDiff/kxDiff*abs(xp(n)-kx(ip))+kyDiff) then
      kk = 3
      point1(1) = kx(ip+1)
      point1(2) = ky(jp)
      point1(3) = kz(ip+1, jp)
      point2(1) = kx(ip)
      point2(2) = ky(jp+1)
      point2(3) = kz(ip, jp+1)
      point3(1) = kx(ip)
      point3(2) = ky(jp)
      point3(3) = kz(ip, jp)
      vec1 = point2 - point1
      vec2 = point3 - point1
      nnvec = cross_prod(vec1, vec2)
    else 
      kk = 4
      point1(1) = kx(ip+1)
      point1(2) = ky(jp)
      point1(3) = kz(ip+1, jp)
      point2(1) = kx(ip)
      point2(2) = ky(jp+1)
      point2(3) = kz(ip, jp+1)
      point3(1) = kx(ip+1)
      point3(2) = ky(jp+1)
      point3(3) = kz(ip+1, jp+1)
      vec1 = point2 - point1
      vec2 = point3 - point1
      nnvec = cross_prod(vec2, vec1)
    end if
  end if
  xpvec1(1) = xp(n)
  xpvec1(2) = yp(n)
  xpvec1(3) = zp(n)
  upvec1(1) = up(n)
  upvec1(2) = vp(n)
  upvec1(3) = wp(n)
  ttvec = cross_prod(nnvec, upvec1)
  tnvec = cross_prod(ttvec, nnvec)
  !
  vec3 = xpvec1 - point1
  lmt3 = (vec1(1)*vec3(2)-vec1(2)*vec3(1))/(vec1(1)*vec2(2)-vec1(2)*vec2(1))
  lmt2 = (vec3(1)-vec2(1)*lmt3)/vec1(1)
  lmt1 = 1.0 - lmt2 - lmt3
  point0(1) = xpvec1(1)
  point0(2) = xpvec1(2)
  point0(3) = lmt1*point1(3) + lmt2*point2(3) + lmt3*point3(3)
  !
  fz(n) = xpvec1(3) - point0(3)
  !
  unnvec = unit_vec(nnvec)
  vin(3) = dot_prod(upvec1, unnvec)
  ifimpact: if (fz(n)<=0.0 .and. vin(3)<0.0) then
    ! information of inject particle
    utnvec = unit_vec(tnvec)
    uttvec = unit_vec(ttvec)
    vin(1) = dot_prod(upvec1, utnvec)
    vin(2) = dot_prod(upvec1, uttvec)
    gg(1) = 0.0
    gg(2) = 0.0
    gg(3) = 9.8
    gg3 = 9.8 !abs(dot_prod(gg, unnvec))
    norm_vin = norm_2(upvec1)
    ! nearest point to impact point
    d01 = dist_p(point0, point1)
    d02 = dist_p(point0, point2)
    d03 = dist_p(point0, point3)
    select case(kk)
    case(1)
      if (d01<=d02 .and. d01<=d03) then
        ipp = ip
        jpp = jp
      else if (d02<=d01 .and. d02<=d03) then
        ipp = ip + 1
        jpp = jp + 1
      else if (d03<=d01 .and. d03<=d02) then
        ipp = ip
        jpp = jp + 1
      end if
    case(2)
      if (d01<=d02 .and. d01<=d03) then
        ipp = ip
        jpp = jp
      else if (d02<=d01 .and. d02<=d03) then
        ipp = ip + 1
        jpp = jp + 1
      else if (d03<=d01 .and. d03<=d02) then
        ipp = ip + 1
        jpp = jp
      end if
    case(3)
      if (d01<=d02 .and. d01<=d03) then
        ipp = ip + 1
        jpp = jp
      else if (d02<=d01 .and. d02<=d03) then
        ipp = ip
        jpp = jp + 1
      else if (d03<=d01 .and. d03<=d02) then
        ipp = ip
        jpp = jp
      end if
    case(4)
      if (d01<=d02 .and. d01<=d03) then
        ipp = ip + 1
        jpp = jp
      else if (d02<=d01 .and. d02<=d03) then
        ipp = ip
        jpp = jp + 1
      else if (d03<=d01 .and. d03<=d02) then
        ipp = ip + 1
        jpp = jp + 1
      end if
    case default
      print*, 'kk error1'
      stop
    end select
    ! properties of injector and particles in the bed
    d1 = dp(n)
    d2 = dpaa(ipp, jpp)
    mm1 = (pi*d1**3)/6.0*rhos
    mm2 = (pi*d2**3)/6.0*rhos
    if (abs(vin(1))>0.0) then
      angin1 = atan(abs(vin(3)/vin(1)))
    else
      angin1 = pi/2.0
    end if
    vpin = vpin + vin
    norm_vpin = norm_vpin + norm_vin
    vvpin = vvpin + norm_vin**2
    mpin = mpin + mm1
    npin = npin + 1.0
    mupin(ip, jp) = mupin(ip, jp) + mm1*vin(1)
    mvpin(ip, jp) = mvpin(ip, jp) + mm1*vin(3)
    !upvec1(1) = upvec1(1) - ucreep(ipp, jpp)
    !upvec1(2) = upvec1(2) - vcreep(ipp, jpp)
    norm_vin = norm_2(upvec1)
    ee1 = 0.5*mm1*norm_vin**2
    ! particle rebound
    if (isp==0) then ! lammel
      dd1 = d1/(0.5*(d1+d2))
      dd2 = d2/(0.5*(d1+d2))
      mmu = els*dd1**3/(dd1**3+els*dd2**3)
      alpha = (1.0+els)/(1.0+mmu) - 1.0
      beta = 1.0 - (2.0/7.0)*(1.0-fric)/(1.0+mmu)
      ! particle rebound
      pp = beta - (beta**2-alpha**2)*dd2*angin1/(2.0*beta)
      angout1 = arebound(alpha, beta, angin1, dd2)
      norm_vout = pp*norm_vin
      ee2 = mm1/2.0*norm_vout**2
      vout(3) = norm_vout*sin(angout1)
      if (vout(3)<sqrt(2.0*gg3*0.5*(d1+d2)) .or. pp<=0.0 .or. pp>1.0 .or. angout1<=0.0) then
        nne = 0
        pp = 0.0
      else
        nne = 1
      end if
    else ! kok
      pp = 0.95*(1.0-exp(-2.0*norm_vin))
      nne = bnldev(pp, 1)
      if (nne>0) then
        mmu = 0.6*norm_vin
        cgma = 0.25*norm_vin
        norm_vout = normal(mmu, cgma)
        ee2 = 0.5*mm2*norm_vout**2
        ammu1 = 30.0/180.0*pi
        acgma1 = 15.0/180.0*pi
        angout1 = normal(ammu1, acgma1)
        if (norm_vout<=0. .or. angout1<=0.) nne = 0
      end if
    end if
    if (nne==0) then
      ! influence of repose angle
      select case(rolldir1(ipp, jpp))
      case(0)
        ii = ipp
        jj = jpp
      case(1)
        ii = ipp+1
        jj = jpp
      case(2)
        ii = ipp-1
        jj = jpp
      case(3)
        ii = ipp
        jj = jpp+1
      case(4)
        ii = ipp
        jj = jpp-1
      case default
        print*, rolldir1(ipp, jpp), 'rolldir1 error'
        stop
      end select
      if (iud/=2) then
        iii = int((d1-dpa+dcgma*3.0)/dcgma/6.0*dfloat(npdf)) + 1
        if (iii>npdf) iii = npdf
        if (iii<=0) iii = 1
      else
        if (d1<dpa) then
          iii = 1
        else if (d1>dpa) then
          iii = 2
        else
          print*, 'error on d1=', d1, '/=', dpa, 'or', dpa+dcgma
          stop
        end if
      end if
      vch = nkl*(pi*d1**3)/6.0
      if (jj>=mky+1) jj = 3
      if (ii<=mkxNode) then
        pnch(ii, jj) = pnch(ii, jj) + vch
        pdfch(ii, jj, iii) = pdfch(ii, jj, iii) + vch
      else
        eepnch(jj) = eepnch(jj) + vch
        jk = iii + (jj-1)*npdf
        eepdfch(jk) = eepdfch(jk) + vch
      end if
    else
      ammu2 = 0.0
      acgma2 = 10.0/180.0*pi
      !angout2 = normal(ammu2, acgma2)
      vout(1) = norm_vout*cos(angout1) !*cos(angout2)
      vout(2) = 0.0 !norm_vout*cos(angout1)*sin(angout2)
      vout(3) = norm_vout*sin(angout1)
      vec1 = vout(1)*utnvec
      vec2 = vout(2)*uttvec
      vec3 = vout(3)*unnvec
      upvec2 = vec1 + vec2 + vec3
      !upvec2(1) = upvec2(1) + ucreep(ipp, jpp)
      !upvec2(2) = upvec2(2) + vcreep(ipp, jpp)
      nnnp = nnnp + 1
      xp(nnnp) = point0(1) !+ upvec2(1)*dt
      yp(nnnp) = point0(2) !+ upvec2(2)*dt
      zp(nnnp) = point0(3) !+ upvec2(3)*dt
      up(nnnp) = upvec2(1)
      vp(nnnp) = upvec2(2)
      wp(nnnp) = upvec2(3)
      dp(nnnp) = d1
      fk(nnnp) = 0.0 !fk(n) !upvec2(1)*dt
      fh(nnnp) = 0.0
      fg(nnnp) = zp(nnnp)
      ft(nnnp) = 0.0
      fz(nnnp) = 0.0 !upvec2(3)*dt
      wflx = wflx + nkl*mm1/xMax/yMax/dt
      vpout = vpout + vout
      norm_vpout = norm_vpout + norm_vout
      vvpout = vvpout + norm_vout**2
      mpout = mpout + mm1
      mupout(ip, jp) = mupout(ip, jp) + mm2*vout(1)
      mvpout(ip, jp) = mvpout(ip, jp) + mm2*vout(2)
      npout = npout + 1.0
    end if
    ! particle splash
    if (isp==0) then ! lammel
      utaot = sqrt(0.0123*(rhos/rho*9.8*dpaa(ipp, jpp)+3.0e-4/(rho*dpaa(ipp, jpp))))
      taot = rho*utaot**2
      eed2x = mm2*gg3*d2
      eed2 = eed2x*(1.0-htao(1)/taot)
      if (eed2/eed2x<=0.1) then
        eed2 = eed2x*0.1
      end if
      lambda = 2.0*log((1.0-pp**2)*ee1/eed2)
      if (lambda<=0.0) then
        ne = 0
      else
        cgma = sqrt(lambda)*log(2.0)
        mmu = log((1.0-pp**2)*ee1) - lambda*log(2.0)
        merfc = myerfc((log(eed2)-mmu)/(sqrt(2.0)*cgma))
        ee2bar = eed2*((1.0-pp**2)*ee1/eed2)**(1.0-(2.0-log(2.0))*log(2.0))
        ne = int(0.06*((1.0-pp**2)*ee1/(2.0*ee2bar))*merfc)
      end if
    else ! kok
      mm2 = (pi*dpaa(ipp, jpp)**3)/6.0*rhos
      nee = 0.03*norm_vin/sqrt(9.8*dpaa(ipp, jpp))
      ne = int(nee)
    end if
    if (ne>=1) then
      ! influence of repose angle
      select case(rolldir2(ipp, jpp))
      case(0)
        ii = ipp
        jj = jpp
      case(1)
        ii = ipp+1
        jj = jpp
      case(2)
        ii = ipp-1
        jj = jpp
      case(3)
        ii = ipp
        jj = jpp+1
      case(4)
        ii = ipp
        jj = jpp-1
      case default
        print*, rolldir2(ipp, jpp), 'rolldir2 error'
        stop
      end select
      splp: do kd = 1, ne
      ammu1 = 60.0/180.0*pi
      ammu2 = 0.0
      acgma1 = 15.0/180.0*pi
      acgma2 = 10.0/180.0*pi
      if (iud==0) then
        ppdf = pdf(ipp, jpp, :)
        dpd = parsiz(ppdf, dpa, dcgma, npdf)
      else if (iud==1) then
        dpd = dpa
      else if (iud==2) then
        ppdf = pdf(ipp, jpp, 1)
        nbi = bnldev(ppdf, 1)
        if (nbi==1) then
          dpd = dpa - dcgma
        else
          dpd = dpa + dcgma
        end if
      end if
      mm2 = (pi*dpd**3)/6.0*rhos
      if (isp==0) then ! lammel
        ee2 = exp(normal(mmu, cgma))
        norm_vout = sqrt(2.0*ee2/mm2)
      else ! kok
        mmu = 0.08*norm_vin
        lambda = 1.0/mmu
        norm_vout = expdev(lambda)
        ee2 = 0.5*mm2*norm_vout**2
      end if
      !lambda = 1.0/ammu1
      !angout1 = expdev(lambda)
      !angout1 = abs(normal(ammu1, acgma1)) !Dupont
      !angout2 = normal(ammu2, acgma2)
      vout(1) = 0.0 !norm_vout*cos(angout1) !*cos(angout2)
      vout(2) = 0.0 !norm_vout*cos(angout1)*sin(angout2)
      vout(3) = norm_vout !*sin(angout1)
      vec1 = vout(1)*utnvec
      vec2 = vout(2)*uttvec
      vec3 = vout(3)*unnvec
      upvec2 = vec1 + vec2 + vec3
      !upvec2(1) = upvec2(1) + ucreep(ipp, jpp)
      !upvec2(2) = upvec2(2) + vcreep(ipp, jpp)
      call random_number(rr1)
      call random_number(rr2)
      call random_number(rr3)
      addnnp = addnnp + 1
      tempx(addnnp) = point0(1) + upvec2(1)*dt*rr1
      tempy(addnnp) = point0(2) + upvec2(2)*dt*rr2
      tempz(addnnp) = point0(3) + upvec2(3)*dt*rr3
      tempu(addnnp) = 0.0 !upvec2(1)
      tempv(addnnp) = 0.0 !upvec2(2)
      tempw(addnnp) = norm_vout !upvec2(3)
      tempd(addnnp) = dpd
      tempfk(addnnp) = 0.0 !upvec2(1)*dt
      tempfh(addnnp) = 0.0
      tempfg(addnnp) = tempz(addnnp)
      tempft(addnnp) = 0.0
      tempfz(addnnp) = 0.0 !upvec2(3)*dt
      wflx = wflx + nkl*mm2/xMax/yMax/dt
      vpout = vpout + vout
      norm_vpout = norm_vpout + norm_vout
      vvpout = vvpout + norm_vout**2
      mpout = mpout + mm2
      mupout(ip, jp) = mupout(ip, jp) + mm2*vout(1)
      mvpout(ip, jp) = mvpout(ip, jp) + mm2*vout(2)
      npout = npout + 1.0
      if (iud/=2) then
        iii = int((dpd-dpa+dcgma*3.0)/dcgma/6.0*dfloat(npdf)) + 1
        if (iii>npdf) iii = npdf
        if (iii<=0) iii = 1
      else
        if (dpd<dpa) then
          iii = 1
        else if (dpd>dpa) then
          iii = 2
        else
          print*, 'error on dpd=', dpd, '/=', dpa, 'or', dpa+dcgma
          stop
        end if
      end if
      vch = nkl*(pi*dpd**3)/6.0
      if (jj==mky+1) jj = 3
      if (ii<=mkxNode) then
        pnch(ii, jj) = pnch(ii, jj) - vch
        pdfch(ii, jj, iii) = pdfch(ii, jj, iii) - vch
      else
        eepnch(jj) = eepnch(jj) - vch
        jk = iii + (jj-1)*npdf
        eepdfch(jk) = eepdfch(jk) - vch
      end if
      end do splp
    end if
  else
    nnnp = nnnp + 1
    xp(nnnp) = xp(n)
    yp(nnnp) = yp(n)
    zp(nnnp) = zp(n)
    up(nnnp) = up(n)
    vp(nnnp) = vp(n)
    wp(nnnp) = wp(n)
    dp(nnnp) = dp(n)
    fk(nnnp) = fk(n)
    fh(nnnp) = fh(n)
    fg(nnnp) = fg(n)
    ft(nnnp) = ft(n)
    fz(nnnp) = fz(n)
  end if ifimpact
  end do dospl
  if (addnnp>=1) then
    nnp = nnnp+addnnp
    if (nnp>nnpmax) then
      print*, nnp, nnpmax
      print*, "particle number reach the threshold"
      stop
    else
      xp(nnnp+1:nnnp+addnnp) = tempx(1:addnnp)
      yp(nnnp+1:nnnp+addnnp) = tempy(1:addnnp)
      zp(nnnp+1:nnnp+addnnp) = tempz(1:addnnp)
      up(nnnp+1:nnnp+addnnp) = tempu(1:addnnp)
      vp(nnnp+1:nnnp+addnnp) = tempv(1:addnnp)
      wp(nnnp+1:nnnp+addnnp) = tempw(1:addnnp)
      dp(nnnp+1:nnnp+addnnp) = tempd(1:addnnp)
      fk(nnnp+1:nnnp+addnnp) = tempfk(1:addnnp)
      fh(nnnp+1:nnnp+addnnp) = tempfh(1:addnnp)
      fg(nnnp+1:nnnp+addnnp) = tempfg(1:addnnp)
      ft(nnnp+1:nnnp+addnnp) = tempft(1:addnnp)
      fz(nnnp+1:nnnp+addnnp) = tempfz(1:addnnp)
    end if
  else
    nnp = nnnp
  end if
  !
  do n = 1, nnp
  dpd = dp(n)
  volp = (pi*dpd**3)/6.0
  mp = rhos*volp
  point0(3) = zp(n) - fz(n)
  if (zp(n)<zu) then
    hpl = 1
    hrlx1 = (zu - point0(3))/(zu-initAveKz)
  else
    hpl = 0
    hrlx1 = 1.0
  end if
  call parloc(i, j, k, kk, xp(n), yp(n), zp(n), fz(n), hrlx1, hpl, dhcld)
  ! calculate particle volume (vpar) within every cell
  vpar(i, j, kk) = vpar(i, j, kk) + volp
  ! calculate particle flux (uflx, uflxz, wflxz)
  uflx = uflx + nkl*mp*up(n)/xMax/yMax
  uflxz(kk) = uflxz(kk) + nkl*mp*up(n)/xMax/yMax/zDiff(kk)
  if (wp(n)>0.0) then
    wflxz(kk) = wflxz(kk) + nkl*mp*wp(n)/xMax/yMax/zDiff(kk)
  end if
  ! up, xp development
  call parvol(up(n), vp(n), wp(n), xp(n), yp(n), fz(n), dpd, mp, k, kk, hrlx1, dhcld, fk(n))
  zp(n) = fz(n) + point0(3)
  if (fz(n)>fh(n)) fh(n) = fz(n)
  if (zp(n)>fg(n)) fg(n) = zp(n)
  ft(n) = ft(n) + dt
  end do
  !
  !aucreep = 0.0
  !qcreepx = 0.0
  !qcreepy = 0.0
  !do jp = 2, mky-1
  !do ip = 2, mkxNode-1
  !dzz = kz(ip, jp) - kz(ip-1, jp)
  !gg1 = 9.8*dzz/sqrt(dzz**2+kxDiff**2)
  !gg2 = 9.8*kxDiff/sqrt(dzz**2+kxDiff**2)
  !qcreepx(ip, jp) = dt*(mupin(ip, jp) - dt*kxDiff*kyDiff*dpa*por*rhos*(gg1+gg2*0.5) - mupout(ip, jp))/kxDiff/rhos/por
  !if (qcreepx(ip, jp)<0.0) qcreepx(ip, jp) = 0.0
  !dzz = kz(ip, jp) - kz(ip, jp-1)
  !gg1 = 9.8*dzz/sqrt(dzz**2+kyDiff**2)
  !gg2 = 9.8*kyDiff/sqrt(dzz**2+kyDiff**2)
  !qcreepy(ip, jp) = dt*(mvpin(ip, jp) - dt*kxDiff*kyDiff*dpa*por*rhos*(gg1+gg2*0.5) - mvpout(ip, jp))/kyDiff/rhos/por
  !if (qcreepy(ip, jp)<0.0) qcreepy(ip, jp) = 0.0
  !ucreep(ip, jp) = qcreepx(ip, jp)/dt/kyDiff/dpa
  !vcreep(ip, jp) = qcreepy(ip, jp)/dt/kxDiff/dpa
  !aucreep = aucreep + ucreep(ip, jp)/dfloat((mkxNode-2)*(mky-2))
  !end do
  !end do
  !call pxch(mkxNode, mky, qcreepx, surfExType, neighbor, comm1D)
  !call pxch(mkxNode, mky, ucreep, surfExType, neighbor, comm1D)
  !call pxch(mkxNode, mky, vcreep, surfExType, neighbor, comm1D)
  !do i = 1, mkxNode
  !qcreepy(i, mky) = qcreepy(i, 2)
  !qcreepy(i, 1) = qcreepy(i, mky-1)
  !ucreep(i, mky) = ucreep(i, 2)
  !ucreep(i, 1) = ucreep(i, mky-1)
  !vcreep(i, mky) = vcreep(i, 2)
  !vcreep(i, 1) = vcreep(i, mky-1)
  !end do
  !mmp = dpa**3*pi/6.0*rhos
  !nump = xMax*yMax*dpa*por/(pi*dpa**3/6.0)
  !uflx = uflx + nump*mmp*aucreep/xMax/yMax
  !do jp = 2, mky-1
  !do ip = 2, mkxNode-1
  !vch = qcreepx(ip-1, jp) - qcreepx(ip, jp) + qcreepy(ip, jp-1) - qcreepy(ip, jp)
  !pnch(ip, jp) = pnch(ip, jp) + vch
  !end do
  !end do
  call surfexch
  ! calculate fluid volume friction phirho
  do k = 1, mz
  do j = 1, my
  do i = 1, mxNode
  phirho(i, j, k) = phirho(i, j, k) - vpar(i, j, k)/(xDiff*yDiff*zDiff(k))
  if (phirho(i, j, k)<1.0-por) phirho(i, j, k) = 1.0 - por + 1.0e-6
  end do
  end do
  end do
  ! boundary condition of phirho
  call gxch(mxNode, my, mz, phirho, comm1D, neighbor, fieldExType, 1)
  do i = 1, mxNode
  do k = 1, mz
  phirho(i, 1, k) = phirho(i, my-1, k)
  phirho(i, my, k) = phirho(i, 2, k)
  end do
  end do
  ! mid-air collision
  ncoll = 0.0
  ntotc = 0.0
  if (icol==1) then
    nkk = 0
    do n = 1, nnp
    i = int((xp(n)-xu(1))/xDiff) + 1
    j = int((yp(n)-yv(1))/yDiff) + 1
    if (i<1) i=1
    if (j<1) j=1
    if (i>mxNode) i=mxNode
    if (j>my) j=my
    nkk(i, j) = nkk(i, j) + 1
    nks = nkk(i, j)
    scp(i, j, nks) = n
    end do
    !
    do j = 1, my
    do i = 1, mxNode
    tnk = nkk(i, j)
    if (tnk<2) cycle
    do nks = 1, tnk-1
    n = scp(i, j, nks)
    kk = mz
    if (fz(n)<=0.0) then
      kk = 1
    else
      do k = 1, mz
      if (fz(n)>=zw(k) .and. fz(n)<zw(k+1)) then
        kk = k
        exit
      end if
      end do
    end if
    do nnks = nks+1, tnk
    ntotc(kk) = ntotc(kk) + 1.0
    nn = scp(i, j, nnks)
    xpvec1(1) = xp(n)
    xpvec1(2) = yp(n)
    xpvec1(3) = zp(n)
    xpvec2(1) = xp(nn)
    xpvec2(2) = yp(nn)
    xpvec2(3) = zp(nn)
    rrr = dist_p(xpvec1, xpvec2)
    rrrd = (dp(n)+dp(nn))/2.0
    if (rrr<=rrrd .and. rrr>0.0) then
      ncoll(kk) = ncoll(kk) + 1.0
      upvec1(1) = up(n)
      upvec1(2) = vp(n)
      upvec1(3) = wp(n)
      upvec2(1) = up(nn)
      upvec2(2) = vp(nn)
      upvec2(3) = wp(nn)
      point1 = xpvec1 + upvec1*dt
      point2 = xpvec2 + upvec2*dt
      nnvec = xpvec2 - xpvec1
      ! unnvec = n1 (1->2), uttvec=n2 (2->1)
      unnvec = unit_vec(nnvec)
      uttvec = -unnvec
      ! vec1=v12 (v1-v2), vec2=v21 (v2-v1)
      vec1 = upvec1 - upvec2
      vec2 = -vec1
      ! nv1=n.v12, nv2=n.v21
      nv1 = dot_prod(unnvec, vec1)
      nv2 = nv1
      tv1 = sqrt((vec1(1)-nv1*unnvec(1))**2+(vec1(2)-nv1*unnvec(2))**2+(vec1(3)-nv1*unnvec(3))**2)
      tv2 = tv1
      mm1 = (pi*dp(n)**3)/6.0*rhos
      mm2 = (pi*dp(nn)**3)/6.0*rhos
      eta1 = mm1/mm2
      eta2 = mm2/mm1
      alpha1 = (1.0+els1)/(1.0+eta1)
      alpha2 = (1.0+els1)/(1.0+eta2)
      !if (tv1/nv1<7.0/2.0*fric1*(1.0+els1)) then
        beta1 = (2.0/7.0)/(1.0+eta1)
        beta2 = (2.0/7.0)/(1.0+eta2)
      !else
      !  beta1 = fric1*(1.0+els1)*nv1/tv1/(1.0+eta1)
      !  beta2 = fric1*(1.0+els1)*nv2/tv2/(1.0+eta2)
      !end if
      !
      up(n) = upvec1(1) - alpha1*nv1*unnvec(1) - beta1*(vec1(1)-nv1*unnvec(1))
      vp(n) = upvec1(2) - alpha1*nv1*unnvec(2) - beta1*(vec1(2)-nv1*unnvec(2))
      wp(n) = upvec1(3) - alpha1*nv1*unnvec(3) - beta1*(vec1(3)-nv1*unnvec(3))
      up(nn) = upvec2(1) - alpha2*nv2*uttvec(1) - beta2*(vec2(1)-nv2*uttvec(1))
      vp(nn) = upvec2(2) - alpha2*nv2*uttvec(2) - beta2*(vec2(2)-nv2*uttvec(2))
      wp(nn) = upvec2(3) - alpha2*nv2*uttvec(3) - beta2*(vec2(3)-nv2*uttvec(3))
      upvec1(1) = up(n)
      upvec1(2) = vp(n)
      upvec1(3) = wp(n)
      upvec2(1) = up(nn)
      upvec2(2) = vp(nn)
      upvec2(3) = wp(nn)
      nv1 = dot_prod(unnvec, upvec1)
      nv2 = dot_prod(unnvec, upvec2)
      if (nv1>=nv2) then
        xpvec2 = xpvec1 - rrrd*unnvec
      else
        xpvec2 = xpvec1 + rrrd*unnvec
      end if
      !xp(n) = xpvec1(1) + dt*up(n)
      !yp(n) = xpvec1(2) + dt*vp(n)
      !zp(n) = xpvec1(3) + dt*wp(n)
      xp0 = xp(nn)
      xp(nn) = xpvec2(1) !+ dt*up(nn)
      yp(nn) = xpvec2(2) !+ dt*vp(nn)
      zp(nn) = xpvec2(3) !+ dt*wp(nn)
      fk(nn) = fk(nn) + xp(nn) - xp0
      !exit
    end if
    end do
    !xp(n) = xpvec1(1) + dt*up(n)
    !yp(n) = xpvec1(2) + dt*vp(n)
    !zp(n) = xpvec1(3) + dt*wp(n)
    end do
    end do
    end do
    do k = 1, mz
    if (ntotc(k)>0.9) then
      pcoll(k) = ncoll(k)/ntotc(k)
    else
      pcoll(k) = 0.0
    end if
    end do
  end if
  ! pick out particles out of boundary
  nnnp = 0
  nn = 0
  nni = 0
  pick: do n = 1, nnp
  if (xp(n)>=xu(mxNode)) then
    nn = nn + 1
    if (myID==nNodes-1) then
      chxp(nn) = xp(n) - xMax
    else
      chxp(nn) = xp(n)
    end if
    chyp(nn) = yp(n)
    chzp(nn) = zp(n)
    chup(nn) = up(n)
    chvp(nn) = vp(n)
    chwp(nn) = wp(n)
    chdp(nn) = dp(n)
    chfk(nn) = fk(n)
    chfz(nn) = fz(n)
    chfh(nn) = fh(n)
    chfg(nn) = fg(n)
    chft(nn) = ft(n)
  else if (xp(n)<xu(2)) then
    nni = nni + 1
    if (xp(n)<0.) then
      chxpi(nni) = xp(n) + xMax
    else
      chxpi(nni) = xp(n)
    end if
    chypi(nni) = yp(n)
    chzpi(nni) = zp(n)
    chupi(nni) = up(n)
    chvpi(nni) = vp(n)
    chwpi(nni) = wp(n)
    chdpi(nni) = dp(n)
    chfki(nni) = fk(n)
    chfzi(nni) = fz(n)
    chfhi(nni) = fh(n)
    chfgi(nni) = fg(n)
    chfti(nni) = ft(n)
  else
    nnnp = nnnp + 1
    xp(nnnp) = xp(n)
    yp(nnnp) = yp(n)
    zp(nnnp) = zp(n)
    up(nnnp) = up(n)
    vp(nnnp) = vp(n)
    wp(nnnp) = wp(n)
    dp(nnnp) = dp(n)
    fk(nnnp) = fk(n)
    fz(nnnp) = fz(n)
    fh(nnnp) = fh(n)
    fg(nnnp) = fg(n)
    ft(nnnp) = ft(n)
  end if
  end do pick
  nnp = nnnp
  ! particle exchange between processes
  call MPI_SENDRECV(nn,1,intType,neighbor(2),18,  &
    nnn,1,intType,neighbor(1),18,comm1D,status,ierr)
  call MPI_SENDRECV(nni,1,intType,neighbor(1),27, &
    nnni,1,intType,neighbor(2),27,comm1D,status,ierr)
  ! from 1 to 2
  allocate(exch(12*nn))
  allocate(exchr(12*nnn))
  exch(1:nn) = chxp(1:nn)
  exch(1*nn+1:2*nn) = chyp(1:nn)
  exch(2*nn+1:3*nn) = chzp(1:nn)
  exch(3*nn+1:4*nn) = chup(1:nn)
  exch(4*nn+1:5*nn) = chvp(1:nn)
  exch(5*nn+1:6*nn) = chwp(1:nn)
  exch(6*nn+1:7*nn) = chdp(1:nn)
  exch(7*nn+1:8*nn) = chfk(1:nn)
  exch(8*nn+1:9*nn) = chfz(1:nn)
  exch(9*nn+1:10*nn) = chfh(1:nn)
  exch(10*nn+1:11*nn) = chfg(1:nn)
  exch(11*nn+1:12*nn) = chft(1:nn)
  call MPI_SENDRECV(exch,nn*12,realType,neighbor(2),10, &
    exchr,nnn*12,realType,neighbor(1),10, &
    comm1D,status,ierr)
  if (nnn>0) then
    xp(nnp+1:nnp+nnn) = exchr(1:nnn)
    yp(nnp+1:nnp+nnn) = exchr(1*nnn+1:2*nnn)
    zp(nnp+1:nnp+nnn) = exchr(2*nnn+1:3*nnn)
    up(nnp+1:nnp+nnn) = exchr(3*nnn+1:4*nnn)
    vp(nnp+1:nnp+nnn) = exchr(4*nnn+1:5*nnn)
    wp(nnp+1:nnp+nnn) = exchr(5*nnn+1:6*nnn)
    dp(nnp+1:nnp+nnn) = exchr(6*nnn+1:7*nnn)
    fk(nnp+1:nnp+nnn) = exchr(7*nnn+1:8*nnn)
    fz(nnp+1:nnp+nnn) = exchr(8*nnn+1:9*nnn)
    fh(nnp+1:nnp+nnn) = exchr(9*nnn+1:10*nnn)
    fg(nnp+1:nnp+nnn) = exchr(10*nnn+1:11*nnn)
    ft(nnp+1:nnp+nnn) = exchr(11*nnn+1:12*nnn)
    if (myID==0) then
      if (ikbx==0) then
        nnp = nnp + nnn
      end if
    else
      nnp = nnp + nnn
    end if
  end if
  deallocate(exch)
  deallocate(exchr)
  ! from 2 to 1
  allocate(exchi(12*nni))
  allocate(exchir(12*nnni))
  exchi(1:nni) = chxpi(1:nni)
  exchi(1*nni+1:2*nni) = chypi(1:nni)
  exchi(2*nni+1:3*nni) = chzpi(1:nni)
  exchi(3*nni+1:4*nni) = chupi(1:nni)
  exchi(4*nni+1:5*nni) = chvpi(1:nni)
  exchi(5*nni+1:6*nni) = chwpi(1:nni)
  exchi(6*nni+1:7*nni) = chdpi(1:nni)
  exchi(7*nni+1:8*nni) = chfki(1:nni)
  exchi(8*nni+1:9*nni) = chfzi(1:nni)
  exchi(9*nni+1:10*nni) = chfhi(1:nni)
  exchi(10*nni+1:11*nni) = chfgi(1:nni)
  exchi(11*nni+1:12*nni) = chfti(1:nni)
  call MPI_SENDRECV(exchi,nni*12,realType,neighbor(1),19, &
    exchir,nnni*12,realType,neighbor(2),19, &
    comm1D,status,ierr)
  if (nnni>0) then
    xp(nnp+1:nnp+nnni) = exchir(1:nnni)
    yp(nnp+1:nnp+nnni) = exchir(1*nnni+1:2*nnni)
    zp(nnp+1:nnp+nnni) = exchir(2*nnni+1:3*nnni)
    up(nnp+1:nnp+nnni) = exchir(3*nnni+1:4*nnni)
    vp(nnp+1:nnp+nnni) = exchir(4*nnni+1:5*nnni)
    wp(nnp+1:nnp+nnni) = exchir(5*nnni+1:6*nnni)
    dp(nnp+1:nnp+nnni) = exchir(6*nnni+1:7*nnni)
    fk(nnp+1:nnp+nnni) = exchir(7*nnni+1:8*nnni)
    fz(nnp+1:nnp+nnni) = exchir(8*nnni+1:9*nnni)
    fh(nnp+1:nnp+nnni) = exchir(9*nnni+1:10*nnni)
    fg(nnp+1:nnp+nnni) = exchir(10*nnni+1:11*nnni)
    ft(nnp+1:nnp+nnni) = exchir(11*nnni+1:12*nnni)
    if (myID==nNodes-1) then
      if (ikbx==0) then
        nnp = nnp + nnni
      end if
    else
      nnp = nnp + nnni
    end if
  end if
  deallocate(exchi)
  deallocate(exchir)
  ! y z boundary condition of particles
  do n = 1, nnp
  if (yp(n)>=yMax) then
    yp(n) = yp(n) - yMax
  else if (yp(n)<0.0) then
    yp(n) = yp(n) + yMax
  end if
  if (zp(n)>zMax) then
    zp(n) = zMax
    wp(n) = -abs(wp(n))
  end if
  end do
  ! bed particle size distributation
  if (iud/=1 .and. irsf==0) then
    ! vlayer change
    do j = 2, mky-1
    do i = 2, mkxNode-1
    vlayer(i, j) = kxDiff*kyDiff*zpdf(i, j)*por
    do k = 1, npdf
    vbin(i, j, k) = pdf(i, j, k)*vlayer(i, j) + pdfch(i, j, k)
    vlayer(i, j) = vlayer(i, j) + pdfch(i, j, k)
    end do
    end do
    end do
    ! calculate pdf, zpdf
    do j = 2, mky-1
    do i = 2, mkxNode-1
    zpdf(i, j) = vlayer(i, j)/(kxDiff*kyDiff)/por
    if (zpdf(i, j)<=0.0) then
      print*, 'zpdf<0'
      zpdf(i, j) = szpdf
      do k = 1, npdf
      pdf(i, j, k) = spdf(k)
      end do
    else if (zpdf(i, j)>=szpdf) then
      do k = 1, npdf
      pdf(i, j, k) = vbin(i, j, k)/vlayer(i, j)
      end do
      if (zpdf(i, j)>2.0*szpdf) zpdf(i, j) = 2.0*szpdf
    else
      svlayer = kxDiff*kyDiff*szpdf*por
      xvlayer = kxDiff*kyDiff*(szpdf-zpdf(i, j))*por
      do k = 1, npdf
      pdf(i, j, k) = (vbin(i, j, k) + xvlayer*spdf(k))/svlayer
      end do
      zpdf(i, j) = szpdf
    end if
    end do
    end do
    ! pdf, zpdf exchange
    ! i=mkxNode-1 send to i=1: send to 2 and receive from 1
    do j = 1, mky
    pzpdfch(j) = zpdf(mkxNode-1, j)
    do k = 1, npdf
    jk = k + (j-1)*npdf
    ppdfch(jk) = pdf(mkxNode-1, j, k)
    end do
    end do
    call MPI_SENDRECV(ppdfch,mky*npdf,realType,neighbor(2),201,  &
      ppdfchr,mky*npdf,realType,neighbor(1),201,comm1D,status,ierr)
    call MPI_SENDRECV(pzpdfch,mky,realType,neighbor(2),203,  &
      pzpdfchr,mky,realType,neighbor(1),203,comm1D,status,ierr)
    do j = 1, mky
    zpdf(1, j) = pzpdfchr(j)
    do k = 1, npdf
    jk = k + (j-1)*npdf
    pdf(1, j, k) = ppdfchr(jk)
    end do
    end do
    ! i=2 send to i=mkxNode: send to 1 and receive from 2
    do j = 1, mky
    izpdfch(j) = zpdf(2, j)
    do k = 1, npdf
    jk = k + (j-1)*npdf
    ipdfch(jk) = pdf(2, j, k)
    end do
    end do
    call MPI_SENDRECV(ipdfch,mky*npdf,realType,neighbor(1),202,  &
      ipdfchr,mky*npdf,realType,neighbor(2),202,comm1D,status,ierr)
    call MPI_SENDRECV(izpdfch,mky,realType,neighbor(1),204,  &
      izpdfchr,mky,realType,neighbor(2),204,comm1D,status,ierr)
    do j = 1, mky
    zpdf(mkxNode, j) = izpdfchr(j)
    do k = 1, npdf
    jk = k + (j-1)*npdf
    pdf(mkxNode, j, k) = ipdfchr(jk)
    end do
    end do
    ! y=mky-1 send to y=1, y=2 send to y=mky
    do i = 2, mkxNode-1
    zpdf(i, 1) = zpdf(i, mky-1)
    zpdf(i, mky) = zpdf(i, 2)
    do k = 1, npdf
    pdf(i, 1, k) = pdf(i, mky-1, k)
    pdf(i, mky, k) = pdf(i, 2, k)
    end do
    end do
    !
    tpdf = 0.0
    do i = 1, mkxNode
    do j = 1, mky
    do k = 1, npdf
    tpdf(i, j) = tpdf(i, j) + pdf(i, j, k)
    end do
    if (tpdf(i, j)<=0.0) then
      print*, tpdf(i, j), 'tpdf<=0'
      stop
    end if
    end do
    end do
    !
    dpaa = 0.0
    do j = 1, mky
    do i = 1, mkxNode
    if (iud==0) then
      do k = 1, npdf
      pdf(i, j, k) = pdf(i, j, k)/tpdf(i, j)
      if (pdf(i, j, k)<0.0) then
        pdf(i, j, k) = 0.0
      end if
      dpaa(i, j) = dpaa(i, j) + pdf(i, j, k)*((dfloat(k-1)+0.5)*dcgma*6.0/dfloat(npdf) + (dpa-dcgma*3.0))
      end do
    else if (iud==2) then
      pdf(i, j, 1) = pdf(i, j, 1)/tpdf(i, j)
      pdf(i, j, 2) = pdf(i, j, 2)/tpdf(i, j)
      if (pdf(i, j, 1)<0.0) then
        pdf(i, j, 1) = 0.0
        pdf(i, j, 2) = 1.0
      else if (pdf(i, j, 2)<0.0) then
        pdf(i, j, 2) = 0.0
        pdf(i, j, 1) = 1.0
      end if
      dpaa(i, j) = pdf(i, j, 1)*(dpa-dcgma) + pdf(i, j, 2)*(dpa+dcgma)
      if (pdf(i, j, 1)+pdf(i, j, 2)>1.01 .or. dpaa(i, j)>dpa+dcgma .or. dpaa(i, j)<dpa-dcgma) then
        print*, 'dpaa error', dpaa(i, j), pdf(i, j, 1), pdf(i, j, 2),  zpdf(i, j), szpdf, vlayer(i, j), &
          vbin(i, j, 1), vbin(i, j, 2)
        stop
      end if
    end if
    if (dpaa(i, j)>dpa+dcgma*3.0) dpaa(i, j) = dpa + dcgma*3.0
    if (dpaa(i, j)<dpa-dcgma*3.0) dpaa(i, j) = dpa - dcgma*3.0
    end do
    end do
  else
    dpaa = dpa
  end if
  !
  deallocate(scp)
  deallocate(chxp, chyp, chzp)
  deallocate(chup, chvp, chwp)
  deallocate(chdp, chfk, chfz)
  deallocate(chfh, chfg, chft)
  deallocate(chxpi, chypi, chzpi)
  deallocate(chupi, chvpi, chwpi)
  deallocate(chdpi, chfki, chfzi)
  deallocate(chfhi, chfgi, chfti)
end subroutine parcalculate

subroutine parloc(iii, jjj, kkk1, kkk2, xxp, yyp, zzp, hhp, hr1, hl, dhcld)
  use public_val
  implicit none
  include "mpif.h"
  ! public
  integer :: iii, jjj, kkk1, kkk2, hl
  double precision, intent(in) :: xxp, yyp, zzp, hhp, hr1
  double precision :: dhcld
  ! local
  integer :: k
  double precision :: zzw, zzw0, zzw1, zzw2, rzp
  !
  iii = int((xxp-xu(1))/xDiff) + 1
  jjj = int((yyp-yv(1))/yDiff) + 1
  if (iii<1) iii=1
  if (jjj<1) jjj=1
  if (iii>mxNode) iii=mxNode
  if (jjj>my) jjj=my
  dhcld = 0.0
  if (hl==1) then
    if (hhp<zw(1)) then
      kkk2 = 1
      kkk1 = kkk2
      dhcld = 0.5*zDiff(1)*hr1
    else
      do k = 1, mz
      if (k==1) then
        zzw = 0.0
      else
        zzw = zzw + zDiff(k-1)*hr1
      end if
      zzw0 = zzw + zDiff(k)*hr1
      if (hhp<zzw0) then
        kkk2 = k
        zzw1 = zzw + 0.5*zDiff(kkk2)*hr1
        if (hhp<zzw1) then
          kkk1 = kkk2
          dhcld = abs(zzw1-hhp)
        else
          kkk1 = kkk2 + 1
          zzw2 = zzw0 + 0.5*zDiff(kkk1)*hr1
          dhcld = abs(zzw2-hhp)
        end if
        exit
      end if
      end do
    end if
  else if (hl==0) then
    rzp = zzp - initAveKz
    if (rzp>=z(mz)) then
      kkk2 = mz
      kkk1 = kkk2
      dhcld = 0.0
    else
      do k = 1, mz
      if (rzp<zw(k+1)) then
        kkk2 = k
        if (rzp<z(kkk2)) then
          kkk1 = kkk2
        else
          kkk1 = kkk2 + 1
        end if
        dhcld = abs(z(kkk1)-rzp)
        exit
      end if
      end do
    end if
  end if
end subroutine parloc

subroutine parvol(uup, vvp, wwp, xxp, yyp, hhp, ddp, mmp, k, kk, hrlx1, dhcld, ffk)
  use public_val
  implicit none
  include "mpif.h"
  ! public
  integer, intent(in) :: k, kk
  double precision :: uup, vvp, wwp
  double precision :: xxp, yyp, hhp, ffk
  double precision, intent(in) :: ddp, mmp, hrlx1, dhcld
  ! local
  double precision :: alpha, beta
  double precision :: ufp, vfp, wfp
  double precision :: d1x, d1y, d1z
  double precision :: d1u, d1v, d1w
  double precision :: d2x, d2y, d2z
  double precision :: d2u, d2v, d2w
  double precision :: d3x, d3y, d3z
  double precision :: d3u, d3v, d3w
  double precision :: d4x, d4y, d4z
  double precision :: d4u, d4v, d4w
  double precision :: xxp0
  ! function
  double precision :: ffd
  ! ufp, vfp, wfp
  if (k>mz) then
    ufp = hru(mz)
  else if (k==1) then
    if (hhp<=0.0) then
      ufp = 0.0
    else
      alpha = hhp/(z(1)*hrlx1)
      ufp = hru(1)*alpha
    end if
  else
    beta = dhcld/(zDiff(k)+zDiff(k-1))/hrlx1*2.0
    alpha = 1.0 - beta
    ufp = hru(k-1)*beta + hru(k)*alpha
  end if
  vfp = 0.0
  wfp = 0.0
  ! up, xp development
  d1x = uup
  d1y = vvp
  d1z = wwp
  d1u = ffd(d1x, ufp, ddp, nu, rho, rhos)/mmp
  d1v = ffd(d1y, vfp, ddp, nu, rho, rhos)/mmp
  d1w = ffd(d1z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0-rho/rhos)
  d2x = uup + dt/2.0*d1u
  d2y = vvp + dt/2.0*d1v
  d2z = wwp + dt/2.0*d1w
  d2u = ffd(d2x, ufp, ddp, nu, rho, rhos)/mmp
  d2v = ffd(d2y, vfp, ddp, nu, rho, rhos)/mmp
  d2w = ffd(d2z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0-rho/rhos)
  d3x = uup + dt/2.0*d2u
  d3y = vvp + dt/2.0*d2v
  d3z = wwp + dt/2.0*d2w
  d3u = ffd(d3x, ufp, ddp, nu, rho, rhos)/mmp
  d3v = ffd(d3y, vfp, ddp, nu, rho, rhos)/mmp
  d3w = ffd(d3z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0-rho/rhos)
  d4x = uup + dt*d3u
  d4y = vvp + dt*d3v
  d4z = wwp + dt*d3w
  d4u = ffd(d4x, ufp, ddp, nu, rho, rhos)/mmp
  d4v = ffd(d4y, vfp, ddp, nu, rho, rhos)/mmp
  d4w = ffd(d4z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0-rho/rhos)
  xxp0 = xxp
  xxp = xxp + (d1x+2.0*d2x+2.0*d3x+d4x)/6.0*dt
  yyp = yyp + (d1y+2.0*d2y+2.0*d3y+d4y)/6.0*dt
  hhp = hhp + (d1z+2.0*d2z+2.0*d3z+d4z)/6.0*dt
  uup = uup + (d1u+2.0*d2u+2.0*d3u+d4u)/6.0*dt
  vvp = vvp + (d1v+2.0*d2v+2.0*d3v+d4v)/6.0*dt
  wwp = wwp + (d1w+2.0*d2w+2.0*d3w+d4w)/6.0*dt
  ffk = ffk + xxp - xxp0
  ampd(kk) = (d1u+2.0*d2u+2.0*d3u+d4u)/6.0*mmp
  ampu(kk) = ampu(kk) + ampd(kk)
end subroutine parvol

subroutine surfexch
  use public_val
  implicit none
  include "mpif.h"
  ! local
  integer :: i, j, k, jk
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  double precision, dimension(mky) :: epnch, epnchr
  double precision, dimension(mky) :: spnch, spnchr
  double precision, dimension(mky*npdf) :: epdfch, epdfchr
  double precision, dimension(mky*npdf) :: spdfch, spdfchr
  ! because the value of ghost cell has changed
  ! need to add ghost value back to real domain before exchange
  ! pnch, pdfch add back
  ! x=mkxNode+1 add to x=3: send to 2 and receive from 1
  ! x=mkxNode add to x=2: send to 2 and receive from 1
  ! x=1 add to x=mkxNode-1: send to 1 and receive from 2
  do j = 1, mky
  epnch(j) = pnch(mkxNode, j)
  spnch(j) = pnch(1, j)
  do k = 1, npdf
  jk = k + (j-1)*npdf
  epdfch(jk) = pdfch(mkxNode, j, k)
  spdfch(jk) = pdfch(1, j, k)
  end do
  end do
  !
  call MPI_SENDRECV(eepnch,mky,realType,neighbor(2),107,  &
    eepnchr,mky,realType,neighbor(1),107,comm1D,status,ierr)
  call MPI_SENDRECV(eepdfch,mky*npdf,realType,neighbor(2),106,  &
    eepdfchr,mky*npdf,realType,neighbor(1),106,comm1D,status,ierr)
  !
  call MPI_SENDRECV(epnch,mky,realType,neighbor(2),103,  &
    epnchr,mky,realType,neighbor(1),103,comm1D,status,ierr)
  call MPI_SENDRECV(epdfch,mky*npdf,realType,neighbor(2),200,  &
    epdfchr,mky*npdf,realType,neighbor(1),200,comm1D,status,ierr)
  !
  call MPI_SENDRECV(spnch,mky,realType,neighbor(1),104,  &
    spnchr,mky,realType,neighbor(2),104,comm1D,status,ierr)
  call MPI_SENDRECV(spdfch,mky*npdf,realType,neighbor(1),108,  &
    spdfchr,mky*npdf,realType,neighbor(2),108,comm1D,status,ierr)
  !
  do j = 1, mky
  pnch(3, j) = pnch(3, j) + eepnchr(j)
  pnch(2, j) = pnch(2, j) + epnchr(j)
  pnch(mkxNode-1, j) = pnch(mkxNode-1, j) + spnchr(j)
  do k = 1, npdf
  jk = k + (j-1)*npdf
  pdfch(3, j, k) = pdfch(3, j, k) + eepdfchr(jk)
  pdfch(2, j, k) = pdfch(2, j, k) + epdfchr(jk)
  pdfch(mkxNode-1, j, k) = pdfch(mkxNode-1, j, k) + spdfchr(jk)
  end do
  end do
  ! y=1 add to y=mky-1, y=mky add to y=2
  do i = 1, mkxNode
  pnch(i, 2) = pnch(i, 2) + pnch(i, mky)
  pnch(i, mky-1) = pnch(i, mky-1) + pnch(i, 1)
  do k = 1, npdf
  pdfch(i, 2, k) = pdfch(i, 2, k) + pdfch(i, mky, k)
  pdfch(i, mky-1, k) = pdfch(i, mky-1, k) + pdfch(i, 1, k)
  end do
  end do
  ! pnch, pdfch exchange
  ! x=2 send to x=mkxNode, x=mkxNode-1 send to x=1
  call pxch(mkxNode, mky, pnch, surfExType, neighbor, comm1D)
  ! y=2 send to y=mky, y=mky-1 send to y=1
  do i = 1, mkxNode
  pnch(i, mky) = pnch(i, 2)
  pnch(i, 1) = pnch(i, mky-1)
  do k = 1, npdf
  pdfch(i, mky, k) = pdfch(i, 2, k)
  pdfch(i, 1, k) = pdfch(i, mky-1, k)
  end do
  end do
end subroutine surfexch

subroutine gxch(mxNode, my, mz, a, comm1D, neighbor, fieldExType, tag)
  implicit none
  include "mpif.h"
  ! public
  integer, intent(in) :: mxNode, my, mz
  integer, intent(in) :: comm1D
  integer, intent(in) :: fieldExType
  integer, intent(in) :: tag
  integer, intent(in), dimension(2) :: neighbor
  double precision, dimension(mxNode, my, mz) :: a
  ! local
  integer :: status(MPI_STATUS_SIZE)
  integer :: ierr
  !
  ! planes i=constant
  !
  ! neighbor:
  !       |           |
  !      ---------------                j
  !       |           |               ^ 
  !       |           |               |
  !      1|    myID   |2              |
  !       |           |              ------>
  !       |           |               |     i
  !      ---------------
  !       |           |
  !
  ! send to 2 and receive from 1
  call MPI_SENDRECV(a(mxNode-1, 1, 1),1,fieldExType,neighbor(2),tag,   &
    a(1, 1, 1),1,fieldExType,neighbor(1),tag,comm1D,status,ierr)
  ! send to 1 and receive from 2
  call MPI_SENDRECV(a(2, 1, 1),1,fieldExType,neighbor(1),tag+1,   &
    a(mxNode, 1, 1),1,fieldExType,neighbor(2),tag+1,comm1D,status,ierr)
end subroutine gxch

subroutine output
  use public_val
  use gather_xyz
  implicit none
  include "mpif.h"
  character(len=3) :: ctemp
  integer :: i, j, k
  integer :: ij
  integer :: blk
  integer :: nf, ns, nc, np, nfi, nfii, nfx, nsf
  integer :: ierr
  integer :: n
  integer :: tnnp
  integer :: numa
  integer, dimension(nNodes) :: cnt, displs
  double precision :: tuflx, twflx
  double precision :: tnorm_vpin, tnorm_vpout, tvvpin, tvvpout, tmpin, tmpout
  double precision :: tnpin, tnpout
  double precision, dimension(3) :: tvpin, tvpout
  double precision, dimension(mx, my, mz) :: tu, tv, tw, tp
  integer, dimension(mx, my, mz) :: tfp
  double precision, dimension(mx) :: tx
  double precision, dimension(mkx) :: tpx
  double precision, dimension(mkx, mky) :: tpz, tpz4
  double precision, dimension((mkxNode-2)*mky) :: apz
  double precision, dimension((mkx-2)*mky) :: tapz
  double precision, dimension((mkxNode-2)*mky) :: apz4
  double precision, dimension((mkx-2)*mky) :: tapz4
  double precision, dimension(mz) :: tuflxz, twflxz
  double precision, dimension(mz) :: tpcoll, apcoll
  double precision, allocatable, dimension(:) :: txp, typ, tzp, tdp, tup, tvp, twp, tfk, tfz, tfh, tfg, tft
  !
  nf = mod(last, nnf)
  ns = mod(last, nns)
  nc = mod(last, nnc)
  np = mod(last, nnkl)
  nsf = mod(last, nnsf)
  nfi = mod(last, nnfi)
  nfx = mod(last, nnfx)
  nfii = (last-nfi)/nnfi
  write(ctemp, '(i3)') nfii
  if (nfi==0) then
    if (myID==0) then
      open (unit=32, file='./particle_loc/particle_loc'//trim(adjustl(ctemp))//'.plt')
      write (32, "(A82)") 'variables = "XP", "YP", "ZP", "DP", "UP", "VP", "WP", "FK", "FZ", "FH", "FG", "FT"'
      close(32)
      !
      open (unit=33, file='./surface/surface'//trim(adjustl(ctemp))//'.plt')
      write (33, *) 'variables = "X", "Y", "Z", "DP"'
      close(33)
      !
      open (unit=42, file='./field3D/field3D'//trim(adjustl(ctemp))//'.plt')
      write (42, *) 'variables = "X", "Y", "Z", "concentration"'
      close(42)
    end if
  end if
  !if (nf==0) then
  !  ! Gather all phirho
  !  call gatherxyz(comm1D, mxNode, mx, my, mz, phirho, tp)
  !  ! Gather all x
  !  call gatherx(comm1D, mxNode, mx, x, tx)
  !  if (myID==0) then
  !    open (unit=42, position='append', file='./field3D/field3D'//trim(adjustl(ctemp))//'.plt')
  !    write (42, *) 'zone', ' T = "', time, '"'
  !    write (42, *) 'i=', mx-2, ' j=', my-2, ' k=', mz, ' datapacking=point'
  !    do k = 1, mz
  !    do j = 2, my-1
  !    do i = 2, mx-1
  !    write (42, "(5E15.7)") tx(i), y(j), z(k), tp(i, j, k)
  !    end do
  !    end do
  !    end do
  !    close(42)
  !  end if
  !end if
  !
  thtaop = thtaop + htaop/dfloat(nns)
  thtao = thtao + htao/dfloat(nns)
  tahff = tahff + ahff/dfloat(nns)
  thru = thru + hru/dfloat(nns)
  ttpcoll = ttpcoll + pcoll/dfloat(nns)
  if (ns==0) then
    call MPI_ALLREDUCE(ttpcoll,tpcoll,mz,realType,MPI_SUM,comm1D,ierr)
    apcoll = tpcoll/nNodes
    if (myID==0) then
      open (unit=43, position='append', file='htao.plt')
      write (43, *) 'zone', ' T = "', time, '"'
      do k = 1, mz
      write (43, "(5E15.7)") z(k), thtao(k), thtaop(k), tahff(k), thru(k), apcoll(k)
      end do
      close(43)
    end if
    thtaop = 0.0
    thtao = 0.0
    tahff = 0.0
    thru = 0.0
    ttpcoll = 0.0
  end if
  !
  if (ikl==1) then
    call MPI_ALLREDUCE(nnp,tnnp,1,intType,MPI_SUM,comm1D,ierr)
    if (nc==0) then
      if (myID==0) then
        open (unit=31, position='append', file='particle_num.plt')
        write (31, "(5E15.7)") time, real(tnnp)
        close(31)
      end if
    end if
  end if
  !
  allocate(txp(tnnp), typ(tnnp), tzp(tnnp), tdp(tnnp), tup(tnnp), tvp(tnnp), twp(tnnp), tfk(tnnp), tfz(tnnp), &
    tfh(tnnp), tfg(tnnp), tft(tnnp))
  if (ikl==1) then
    if (np==0) then
      if (last>=pistart) then
        displs(1) = 0
        call MPI_GATHER(nnp,1,intType,cnt,1,intType,0,comm1D,ierr)
        do i = 2, nNodes
        displs(i) = displs(i-1) + cnt(i-1)
        end do
        call MPI_GATHERV(xp,nnp,realType,txp,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(yp,nnp,realType,typ,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(zp,nnp,realType,tzp,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(dp,nnp,realType,tdp,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(up,nnp,realType,tup,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(vp,nnp,realType,tvp,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(wp,nnp,realType,twp,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(fk,nnp,realType,tfk,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(fz,nnp,realType,tfz,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(fh,nnp,realType,tfh,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(fg,nnp,realType,tfg,cnt,displs,realType,0,comm1D,ierr)
        call MPI_GATHERV(ft,nnp,realType,tft,cnt,displs,realType,0,comm1D,ierr)
        if (myID==0) then
          open (unit=32, position='append', file='./particle_loc/particle_loc'//trim(adjustl(ctemp))//'.plt')
          write (32, *) 'zone', ' T = "', time, '"'
          do n = 1, tnnp
          write (32, "(5E15.7)") txp(n), typ(n), tzp(n), tdp(n), tup(n), tvp(n), &
            twp(n), tfk(n), tfz(n), tfh(n), tfg(n), tft(n)
          end do
          close(32)
        end if
      end if
    end if
  end if
  !
  if (nsf==0) then
    ! Gather all kx to myID=0
    call gatherx(comm1D, mkxNode, mkx, kx, tpx)
    ! Gather all kz and dpaa to myID=0
    do j = 1, mky
    do i = 2, mkxNode-1
    ij = (i-1) + (j-1)*(mkxNode-2)
    apz(ij) = kz(i, j)
    apz4(ij) = dpaa(i, j)
    end do
    end do
    numa = (mkxNode-2)*mky
    call MPI_ALLGATHER(apz,numa,realType,tapz,numa,realType,comm1D,ierr)
    call MPI_ALLGATHER(apz4,numa,realType,tapz4,numa,realType,comm1D,ierr)
    do j = 1, mky
    do i = 2, mkx-1
    blk = (i-2)/(mkxNode-2)
    ij = (i-1) - blk*(mkxNode-2) + (j-1)*(mkxNode-2) + blk*(mkxNode-2)*mky
    tpz(i, j) = tapz(ij)
    tpz4(i, j) = tapz4(ij)
    end do
    end do
    if (myID==0) then
      open (unit=33, position='append', file='./surface/surface'//trim(adjustl(ctemp))//'.plt')
      write (33, *) 'zone', ' T = "', time, '"'
      write (33, *) 'i=', mkx-2, ' j=', mky-2, ' datapacking=point'
      do j = 2, mky-1
      do i = 2, mkx-1
      write (33, "(5E15.7)") tpx(i), ky(j), tpz(i, j), tpz4(i, j)
      end do
      end do
      close(33)
      !
      !open (unit=34, position='append', file='./surfaced/surfaced'//trim(adjustl(ctemp))//'.plt')
      !write (34, *) 'zone', ' T = "', time, '"'
      !write (34, *) 'i=', mkx-2, ' j=', mky-2, ' datapacking=point'
      !do j = 2, mky-1
      !do i = 2, mkx-1
      !write (34, "(5E15.7)") tpx(i), ky(j), tpz4(i, j)
      !end do
      !end do
      !close(34)
    end if
  end if
  call MPI_ALLREDUCE(vpin,tvpin,3,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(vpout,tvpout,3,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(norm_vpin,tnorm_vpin,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(norm_vpout,tnorm_vpout,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(vvpin,tvvpin,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(vvpout,tvvpout,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(mpin,tmpin,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(mpout,tmpout,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(npin,tnpin,1,realType,MPI_SUM,comm1D,ierr)
  call MPI_ALLREDUCE(npout,tnpout,1,realType,MPI_SUM,comm1D,ierr)
  if (tnpin<=0.0) then
    tot_vpin = tot_vpin + 0.0
    tot_nvpin = tot_nvpin + 0.0
    tot_vvpin = tot_vvpin + 0.0
    tot_mpin = tot_mpin + 0.0
    tot_npin = tot_npin + 0.0
  else
    tot_vpin = tot_vpin + tvpin/tnpin
    tot_nvpin = tot_nvpin + tnorm_vpin/tnpin
    tot_vvpin = tot_vvpin + tvvpin/tnpin
    tot_mpin = tot_mpin + tmpin/tnpin
    tot_npin = tot_npin + tnpin
  end if
  if (tnpout<=0.0) then
    tot_vpout = tot_vpout + 0.0
    tot_nvpout = tot_nvpout + 0.0
    tot_vvpout = tot_vvpout + 0.0
    tot_mpout = tot_mpout + 0.0
    tot_npout = tot_npout + 0.0
  else
    tot_vpout = tot_vpout + tvpout/tnpout
    tot_nvpout = tot_nvpout + tnorm_vpout/tnpout
    tot_vvpout = tot_vvpout + tvvpout/tnpout
    tot_mpout = tot_mpout + tmpout/tnpout
    tot_npout = tot_npout + tnpout
  end if
  if (nfx==0) then
    call MPI_ALLREDUCE(uflx,tuflx,1,realType,MPI_SUM,comm1D,ierr)
    call MPI_ALLREDUCE(wflx,twflx,1,realType,MPI_SUM,comm1D,ierr)
    call MPI_ALLREDUCE(uflxz,tuflxz,mz,realType,MPI_SUM,comm1D,ierr)
    call MPI_ALLREDUCE(wflxz,twflxz,mz,realType,MPI_SUM,comm1D,ierr)
    tuflx = tuflx/dfloat(nnfx)
    twflx = twflx/dfloat(nnfx)
    tuflxz = tuflxz/dfloat(nnfx)
    twflxz = twflxz/dfloat(nnfx)
    tot_vpin = tot_vpin/dfloat(nnfx)
    tot_vpout = tot_vpout/dfloat(nnfx)
    tot_nvpin = tot_nvpin/dfloat(nnfx)
    tot_nvpout = tot_nvpout/dfloat(nnfx)
    tot_vvpin = tot_vvpin/dfloat(nnfx)
    tot_vvpout = tot_vvpout/dfloat(nnfx)
    tot_mpin = tot_mpin/dfloat(nnfx)
    tot_mpout = tot_mpout/dfloat(nnfx)
    tot_npin = tot_npin/dfloat(nnfx)
    tot_npout = tot_npout/dfloat(nnfx)
    if (myID==0) then
      open (unit=35, position='append', file='average_flux.plt')
      if (twflx/=0.) then
        write (35, "(5E15.7)") time, tuflx, twflx, tuflx/twflx
      else
        write (35, "(5E15.7)") time, tuflx, twflx, 0.0
      end if
      close(35)
      !
      open (unit=36, position='append', file='flux_vs_height.plt')
      write (36, *) 'zone', ' T = "', time, '"'
      do k = 1, mz
      write (36, "(5E15.7)") z(k), tuflxz(k), twflxz(k)
      end do
      close(36)
      !
      open (unit=39, position='append', file='vin.plt')
      write (39, "(5E15.7)") time, tot_vpin(1), tot_vpin(2), tot_vpin(3), tot_nvpin
      close(39)
      !
      open (unit=46, position='append', file='vout.plt')
      write (46, "(5E15.7)") time, tot_vpout(1), tot_vpout(2), tot_vpout(3), tot_nvpout
      close(46)
      !
      open (unit=44, position='append', file='eminout.plt')
      write (44, "(5E15.7)") time, tot_vvpin, tot_vvpout, tot_mpin, tot_mpout
      close(44)
      !
      open (unit=45, position='append', file='numinout.plt')
      write (45, "(5E15.7)") time, tot_npin, tot_npout
      close(45)
    end if
    uflx = 0.0
    wflx = 0.0
    uflxz = 0.0
    wflxz = 0.0
    tot_vpin = 0.0
    tot_vpout = 0.0
    tot_nvpin = 0.0
    tot_nvpout = 0.0
    tot_vvpin = 0.0
    tot_vvpout = 0.0
    tot_mpin = 0.0
    tot_mpout = 0.0
    tot_npin = 0.0
    tot_npout = 0.0
  end if
  !
  deallocate(txp, typ, tzp, tdp, tup, tvp, twp, tfk, tfz, tfh, tfg, tft)
end subroutine output

function ffd(upp, ufp, ddp, nu, rho, rhos)
  implicit none
  ! public
  double precision :: ffd
  double precision, intent(in) :: upp, ufp
  double precision, intent(in) :: ddp
  double precision, intent(in) :: nu, rho, rhos
  ! local
  double precision :: cd
  double precision :: rep, frep
  double precision :: beta
  double precision :: mp
  double precision :: ttp
  double precision, parameter :: pi = acos(-1.0)
  !
  rep = abs(upp-ufp)*ddp/nu
  if (rep==0.) then
    ffd = 0.0
  else
    !if (rep<1.0) then
    !  frep = 1.0
    !else if (rep<1000.0) then
    !  frep = 1.0 + 0.15*rep**0.687
    !else
    !  frep = 0.0183*rep
    !end if
    !cd = 24./rep*frep
    cd = (sqrt(0.5)+sqrt(24.0/rep))**2
    ffd = -pi/8.*cd*rho*ddp**2*abs(upp-ufp)*(upp-ufp)
    !mp = rhos*pi*ddp**3/6.0
    !ttp = rhos*ddp**2/18.0/rho/nu
    !ffd = mp*(ufp-upp)/ttp*frep
  end if
end function ffd

function bnldev(pp, mm)
  implicit none
  ! public
  integer :: bnldev
  integer, intent(in) :: mm
  double precision, intent(in) :: pp
  ! local
  integer :: i
  double precision :: rr
  !
  bnldev = 0
  do i = 1, mm
  call random_number(rr)
  if (rr<=pp) then
    bnldev = bnldev + 1
  end if
  end do
end function bnldev

function normal(mmu, sigma)
  implicit none
  ! public
  double precision :: normal
  double precision, intent(in) :: mmu, sigma
  !local
  integer :: flg
  double precision :: pi, u1, u2, y1, y2
  save flg
  data flg /0/
  parameter(pi = acos(-1.))
  call random_number(u1)
  call random_number(u2)
  if (flg==0) then
    y1 = sqrt(-2.0*log(u1))*cos(2.0*pi*u2)
    normal = mmu + sigma*y1
    flg = 1
  else
    y2 = sqrt(-2.0*log(u1))*sin(2.0*pi*u2)
    normal = mmu + sigma*y2
    flg = 0
  end if
end function normal

function expdev(lambda)
  implicit none
  ! public
  double precision :: expdev
  double precision, intent(in) :: lambda
  ! local
  double precision :: pv
  !
  do while (.true.)
  call random_number(pv)
  if (pv<1.) then
    exit
  end if
  end do
  expdev = (-1./lambda)*log(1.-pv)
  !
  return
end function expdev

function myerfc(x)
  implicit none
  ! public
  double precision :: myerfc
  double precision, intent(in) :: x
  ! function
  double precision :: gammp, gammq
  ! local
  double precision :: a
  !
  a = 0.5
  if (x<0.0) then
    myerfc = 1.0 + gammp(a, x**2)
  else
    myerfc = gammq(a, x**2)
  end if
end function myerfc

function gammq(a, x)
  implicit none
  ! public
  double precision :: gammq
  double precision, intent(in) :: a
  double precision, intent(in) :: x
  ! local
  integer :: igser, igcf
  double precision :: gammcf, gamser
  !
  if (x<0. .or. a<=0.) return
  if (x<a+1.0) then
    call gser(gamser, a, x, igser)
    if (igser==0) then
      gammq = 1.0 - gamser
    else
      gammq = 1.0
    end if
  else
    call gcf(gammcf, a, x, igcf)
    if (igcf==0) then
      gammq = gammcf
    else
      gammq = 1.0
    end if
  end if
end function gammq

function gammp(a, x)
  implicit none
  ! public
  double precision :: gammp
  double precision, intent(in) :: a
  double precision, intent(in) :: x
  ! local
  integer :: igser, igcf
  double precision :: gammcf, gamser
  !
  if (x<0. .or. a<=0.) return
  if (x<a+1.0) then
    call gser(gamser, a, x, igser)
    if (igser==0) then
      gammp = gamser
    else
      gammp = 1.0
    end if
  else
    call gcf(gammcf, a, x, igcf)
    if (igcf==0) then
      gammp = 1.0 - gammcf
    else
      gammp = 1.0
    end if
  end if
end function gammp

subroutine gser(gamser, a, x, igser)
  implicit none
  ! public
  integer :: igser
  double precision :: gamser
  double precision, intent(in) :: a
  double precision, intent(in) :: x
  ! function
  double precision :: gammln
  ! local
  integer, parameter :: itmax = 100
  integer :: n
  double precision, parameter :: eps = 3.0e-7
  double precision :: gln
  double precision :: ap
  double precision :: del
  double precision :: asum
  !
  gln = gammln(a)
  if (x<=0.) then
    if (x<0.) return
    gamser=0.
    return
  end if
  ap = a
  asum = 1.0/a
  del = asum
  do n = 1, itmax
  ap = ap + 1.0
  del = del*x/ap
  asum = asum + del
  if (abs(del)<abs(asum)*eps) goto 4000
  end do
  print*, 'iter num > itmax in gser'
  igser = 1
  return
  4000 continue
  gamser = asum*exp(-x+a*log(x)-gln)
  igser = 0
end subroutine gser

subroutine gcf(gammcf, a, x, igcf)
  implicit none
  ! public
  integer :: igcf
  double precision :: gammcf
  double precision, intent(in) :: a
  double precision, intent(in) :: x
  ! function
  double precision :: gammln
  ! local
  integer, parameter :: itmax = 100
  integer :: i
  double precision :: gln
  double precision, parameter :: eps = 3.0e-7
  double precision, parameter :: fpmin = 1.0e-30
  double precision :: an, b, c, d
  double precision :: del
  double precision :: h
  !
  gln = gammln(a)
  b = x + 1.0 - a
  c = 1.0/fpmin
  d = 1.0/b
  h = d
  do i = 1, itmax
  an = -i*(i-a)
  b = b + 2.0
  d = an*d + b
  if (abs(d)<fpmin) d = fpmin
  c = b + an/c
  if (abs(c)<fpmin) c = fpmin
  d = 1.0/d
  del = d*c
  h = h*del
  if (abs(del-1.0)<eps) goto 4001
  end do
  print*, 'iter num > itmax in gcf'
  igcf = 1
  return
  4001 continue
  gammcf = exp(-x+a*log(x)-gln)*h
  igcf = 0
end subroutine gcf

function gammln(xx)
  implicit none
  ! public
  double precision :: gammln
  double precision, intent(in) :: xx
  ! local
  integer :: j
  double precision :: ser
  double precision :: tmp
  double precision :: x, y
  double precision, save :: cof(6)
  double precision, save :: stp
  data cof, stp/76.18009172947146d0,-86.50532032941677d0, &
    24.01409824083091d0,-1.231739572450155d0, &
    .1208650973866179d-2,-.5395239384953d-5,  &
    2.5066282746310005d0/
  !
  x = xx
  y = x
  tmp = x + 5.5d0
  tmp = (x+0.5d0)*log(tmp) - tmp
  ser = 1.000000000190015d0
  do j = 1, 6
  y = y + 1.d0
  ser = ser + cof(j)/y
  end do
  gammln = tmp + log(stp*ser/x)
end function gammln

function parsiz(ppp, dpa, dcgma, npdf)
  implicit none
  ! public
  double precision :: parsiz
  integer, intent(in) :: npdf
  double precision, intent(in), dimension(npdf) :: ppp
  double precision, intent(in) :: dpa, dcgma
  ! local
  integer :: n, ii
  double precision :: rr, pdf
  double precision :: x, y, r1, r2
  double precision :: cgma, mu
  cgma = dcgma*1.0e4
  mu = dpa*1.0e4
  y = 1.0
  pdf = 0.0
  rr = maxval(ppp)
  n = 0
  do while (y>pdf)
  n = n+1
  call random_number(r1)
  call random_number(r2)
  x = cgma*6.0*r1 + (mu-cgma*3.0)
  y = r2*rr
  ii = int((x-mu+cgma*3.0)/cgma/6.0*dfloat(npdf)) + 1
  if (ii>npdf .or. ii<1) cycle
  pdf = ppp(ii)
  if (n>10000) then
    print*, 'parsiz, n>10000', ppp
    x = dpa
    exit
  end if
  end do
  parsiz = x*1.0e-4
end function parsiz

function arebound(alpha, beta, angin, dp2)
  implicit none
  ! public
  double precision :: arebound
  double precision, intent(in) :: alpha, beta
  double precision, intent(in) :: angin, dp2
  ! local
  integer :: n, iii
  double precision :: x, y
  double precision :: r1, r2
  double precision :: pdf
  double precision :: gama
  double precision :: xMax, xmin, xmid
  double precision :: da
  double precision, parameter :: pi = acos(-1.0)
  !
  gama = 4.0/9.0*beta**2/(alpha+beta)**2/dp2
  xmin = -angin
  xMax = sqrt(angin/gama)*2.0 - angin
  if (xMax>pi) then
    xMax = pi
  end if
  da = (xMax-xmin)/500.0
  y = 1.0
  pdf = 0.0
  n = 0
  do
  call random_number(r1)
  call random_number(r2)
  n = n + 1
  if (n>10000) then
    print*, 'arebound, n>10000', dp2
    x = 60.0/180.0*pi
    exit
  end if
  iii = int(r1*(xMax-xmin)/da)
  x = dfloat(iii)*da + 0.5*da + xmin
  y = r2
  pdf = gama*(x+angin)/angin*log(2.0/gama*angin/(x+angin)**2)
  pdf = pdf*da
  if (y<=pdf) exit
  end do
  arebound = x + angin
end function arebound
