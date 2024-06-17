module public_parameter
   ! constants
   implicit none
   integer, parameter :: dbPc = selected_real_kind(14, 307)
   real(kind=dbPc), parameter :: pi = 3.14159265358979323846
   ! computational domain
   real(kind=dbPc), parameter :: xmax = 1.0 ! x size
   real(kind=dbPc), parameter :: ymax = 0.2 ! y size
   real(kind=dbPc), parameter :: zmax = 0.5 ! z size
   real(kind=dbPc), parameter :: zUni = 0.1 ! the altitude that grid z becomes uniform
   integer, parameter :: nx = 500 ! x grid num
   integer, parameter :: ny = 100 ! y grid num
   integer, parameter :: nz = 250 ! z grid num
   real(kind=dbPc), parameter :: zGridCtrl = 1.06 ! the control parameter of z grid
   integer, parameter :: nkx = 1000 ! x key point num
   integer, parameter :: nky = 200 ! y key point num
   integer, parameter :: nNodes = 5 ! num of subdomain
   ! time
   real(kind=dbPc), parameter :: dt = 1.0e-4 ! time step
   real(kind=dbPc), parameter :: tla = 600.0 ! time last
   ! fluid
   real(kind=dbPc), parameter :: wind = 0.5 ! fractional velocity
   real(kind=dbPc), parameter :: rho = 1.263 ! fluid density
   real(kind=dbPc), parameter :: nu = 1.51e-5 ! kinetic viscosity
   real(kind=dbPc), parameter :: kapa = 0.4 ! von Kaman's constant
   ! particle
   integer, parameter :: ipar = 1 ! calculating particles: ipar = 0: no, 1: yes
   integer, parameter :: npdf = 3 ! bin num of particle distribution. must be odd number when ipd=0
   integer, parameter :: pNumInit = 100 ! initial particle num
   real(kind=dbPc), parameter :: repoAng = 30.0 ! repose angle
   real(kind=dbPc), parameter :: els = 0.9 ! normal restitution coefficient
   real(kind=dbPc), parameter :: fric = 0.0 ! tangential restitution coefficient
   real(kind=dbPc), parameter :: els1 = 0.9 ! normal restitution coefficient (mid-air collision)
   real(kind=dbPc), parameter :: fric1 = 0.0 ! tangential restitution coefficient (mid-air collision)
   real(kind=dbPc), parameter :: dpa = 2.5e-4 ! average particle diameter
   real(kind=dbPc), parameter :: dSigma = 2.0e-4 ! particle diameter standard deviation x 2
   real(kind=dbPc), parameter :: rpdf = 3.0*dSigma ! the half range of particle distribution.
   real(kind=dbPc), parameter :: bedCellTknessInit = dpa*5.0 ! thickness of the thin surface layer on particle bed
   real(kind=dbPc), parameter :: rhos = 2650.0 ! particle density
   real(kind=dbPc), parameter :: nkl = 1.0 ! one particle stands for x particles
   real(kind=dbPc), parameter :: por = 0.6 ! bedform porosity
   integer, parameter :: pNumMax = 100000 ! max particle num in one subdomain
   integer, parameter :: nspmax = 10000 ! max eject particle num in one time step
   ! pre-formed surface
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
   ! particle diameter: ipd = 0: polydisperse, 1: monodisperse, 2: bidisperse
   ! if ipd=0, normal distribution, mu=dpa, sigma=dSigma, dpa-rpdf*dSigma~dpa+rpdf*dSigma
   ! if ipd=1, d=dpa.
   ! if ipd=2, d1=dpa-dSigma, d2=dpa+dSigma, npdf must equal to 2.
   integer, parameter :: ipd = 1
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
end module public_parameter

module grid_operations
   use public_parameter
   implicit none
   private
   public :: generateGrid

   type :: gridType
      real(kind=dbPc) :: xDiff
      real(kind=dbPc) :: yDiff
      real(kind=dbPc), dimension(mz) :: zDiff
      real(kind=dbPc), dimension(mx+1) :: xu
      real(kind=dbPc), dimension(my+1) :: yv
      real(kind=dbPc), dimension(mz+1) :: zw
      real(kind=dbPc), dimension(mx) :: x
      real(kind=dbPc), dimension(my) :: y
      real(kind=dbPc), dimension(mz) :: z
   end type gridType

   type(gridType) :: grid
   public :: grid

contains

   subroutine generateGrid
      integer :: i, j, k

      ! Define delta_x, delta_y, delta_z of grid
      grid%xDiff = xmax / nx
      grid%yDiff = ymax / ny
      do k = 1, mz
         grid%zDiff(k) = zGridCtrl**(k-1) * zmax / ((zGridCtrl**mz - 1.0) / (zGridCtrl - 1.0))
      end do

      ! Define the grid node coordinates for storing vectors in a staggered grid
      grid%xu(1) = -grid%xDiff
      do i = 2, mx + 1
         grid%xu(i) = grid%xu(i-1) + grid%xDiff
      end do

      grid%yv(1) = -grid%yDiff
      do j = 2, my + 1
         grid%yv(j) = grid%yv(j-1) + grid%yDiff
      end do

      grid%zw(1) = 0.0
      do k = 2, mz + 1
         grid%zw(k) = grid%zw(k-1) + grid%zDiff(k-1)
      end do

      ! Define the grid node coordinates for storing scalars in a staggered grid
      do i = 1, mx
         grid%x(i) = (grid%xu(i+1) - grid%xu(i)) / 2.0 + grid%xu(i)
      end do

      do j = 1, my
         grid%y(j) = (grid%yv(j+1) - grid%yv(j)) / 2.0 + grid%yv(j)
      end do

      do k = 1, mz
         grid%z(k) = (grid%zw(k+1) - grid%zw(k)) / 2.0 + grid%zw(k)
      end do
   end subroutine generateGrid

end module grid_operations

module surface_grid_operations
   use public_parameter
   implicit none
   private
   public :: generateSurfGrid

   type :: surfaceGridType
      real(kind=dbPc) :: xDiff
      real(kind=dbPc) :: yDiff
      real(kind=dbPc), dimension(mkx) :: x
      real(kind=dbPc), dimension(mky) :: y
      real(kind=dbPc), dimension(mkx, mky) :: z
   end type surfaceGridType

   type(surfaceGridType) :: surfGrid
   public :: surfGrid

contains

   subroutine generateSurfGrid
      integer :: i, j

      ! Define delta_x, delta_y of surface grid
      surfGrid%xDiff = xmax / nkx
      surfGrid%yDiff = ymax / nky

      ! Define the surface grid node coordinates and the corresponding elevation
      surfGrid%x(1) = -surfGrid%xDiff
      do i = 2, mkx
         surfGrid%x(i) = surfGrid%x(i-1) + surfGrid%xDiff
      end do

      surfGrid%y(1) = -surfGrid%yDiff
      do j = 2, mky
         surfGrid%y(j) = surfGrid%y(j-1) + surfGrid%yDiff
      end do

      surfGrid%z = initSurfElevation
   end subroutine generateSurfGrid

end module surface_grid_operations

module output_file_generation
   use public_parameter
   implicit none
   private
   public :: outPutFile

contains

   subroutine outPutFile
      character(len=32) bashCmd

      bashCmd = 'mkdir particle_loc'
      call system(trim(adjustl(bashCmd)))
      bashCmd = 'mkdir surface'
      call system(trim(adjustl(bashCmd)))
      bashCmd = 'mkdir field3D'
      call system(trim(adjustl(bashCmd)))

      open (unit=32, file='./particle_loc/particle_loc0.plt')
      write (32, "(A82)") 'variables = "XP", "YP", "ZP", "DP", "UP", "VP", "WP", "FK", "FZ", "FH", "FG", "FT"'
      close(32)

      open (unit=33, file='./surface/surface0.plt')
      write (33, *) 'variables = "X", "Y", "Z", "DP"'
      close(33)

      open (unit=42, file='./field3D/field3D0.plt')
      write (42, *) 'variables = "X", "Y", "Z", "concentration"'
      close(42)

      open (unit=31, file='particle_num.plt')
      write (31, *) 'variables = "T", "Num"'
      close(31)

      open (unit=35, file='average_flux.plt')
      write (35, *) 'variables = "T", "uFlux", "wFlux", "salength"'
      close(35)

      open (unit=36, file='flux_vs_height.plt')
      write (36, *) 'variables = "Z", "uFlux", "wFlux"'
      close(36)

      open (unit=43, file='htao.plt')
      write (43, *) 'variables = "Z", "taoa", "taop", "vfrac", "u", "fptx"'
      close(43)

      open (unit=39, file='vin.plt')
      write (39, *) 'variables = "T", "upin", "vpin", "wpin", "norm_vpin"'
      close(39)

      open (unit=46, file='vout.plt')
      write (46, *) 'variables = "T", "upout", "vpout", "wpout", "norm_vpout"'
      close(46)

      open (unit=44, file='eminout.plt')
      write (44, *) 'variables = "T", "vvpin", "vvpout", "mpin", "mpout"'
      close(44)

      open (unit=45, file='numinout.plt')
      write (45, *) 'variables = "T", "numin", "numout"'
      close(45)
   end subroutine outPutFile

end module output_file_generation

program main
   use public_parameter
   use grid_operations
   use surface_grid_operations
   use output_file_generation
   implicit none
   integer :: last

   ! generate grid
   call generateGrid
   ! generate initial bed
   call generateSurfGrid
   ! creat output file
   call outPutFile
   !
   ! start iteration loop
   !
   last = 1
   !do
   !   ! calculate particle movement
   !   if (ipar==1) then
   !      if (last>=pstart) then
   !         if (last==pstart) then
   !            call particleInit
   !         end if
   !         call particleCal
   !      end if
   !      if (last<sstart) then
   !         Dkz = 0.0
   !         do i = 1, mkxNode
   !            do j = 1, mky
   !               bedCellTkness(i, j) = bedCellTknessInit
   !               if (irsf==0) then
   !                  do k = 1, npdf
   !                     bedPDist(i, j, k) = prob(k)
   !                  end do
   !                  bedPD(i, j) = dpa
   !               else
   !                  bedPDist(i, j, 2) = 0.5*(0.5*sin(initOmg*kx(i))+0.5)
   !                  bedPDist(i, j, 1) = 1.0 - bedPDist(i, j, 2)
   !                  bedPD(i, j) = bedPDist(i, j, 1)*(dpa-dSigma) + bedPDist(i, j, 2)*(dpa+dSigma)
   !               end if
   !            end do
   !         end do
   !      end if
   !   end if
   !   ! calculate fluid field
   !   call fluidField
   !   ! generate boundary key point
   !   call imgd
   !   phirho = 1.0
   !   ! output result
   !   call output
   !   ! time advance
   !   time = time + dt
   !   last = last + 1
   !   if (time>tla) exit
   !end do
end program main

