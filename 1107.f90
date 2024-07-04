! 不用dfloat，改用隐式类型转换，即直接将整型变量赋值给另一个双精度类型变量
! function(arry)采用假设形状数组，可通过size命令得到数组中元素个数
! 修改变量及函数命名
! 将共用函数放置于module中，专用函数放置于program或subroutine的contains中
! 定义函数时使用result属性来指定一个变量作为函数的返回值：function square(x) result(y)
! 调用这个函数：result=square(x)
! 计算分布时不用乘1e4
! 使用派生类型传递全局变量
module public_val
    ! define in dataExType
    integer :: fieldExType
    integer :: surfExType
    ! define in main
    integer :: realType
    integer :: intType
    integer :: comm1D
    integer :: myID
    integer, dimension(2) :: neighbor
    integer :: last
    ! define in generateGrid
    ! define in generateInitBed
    real(kind=dbPc) :: kxDiff, kyDiff
    real(kind=dbPc), dimension(mkxNode) :: kx
    real(kind=dbPc), dimension(mky) :: ky
    real(kind=dbPc), dimension(mkxNode, mky) :: bedPD
    real(kind=dbPc), dimension(mkxNode, mky) :: bedCellTkness
    ! define in particleInit
    real(kind=dbPc), dimension(npdf) :: prob
    integer :: pNum
    integer, dimension(pNumMax, 5) :: pLoc
    integer, dimension(pNumMax, 2) :: pTag
    real(kind=dbPc), dimension(pNumMax) :: xp, yp, zp
    real(kind=dbPc), dimension(pNumMax) :: up, vp, wp
    real(kind=dbPc), dimension(pNumMax) :: dp, fk, fz
    real(kind=dbPc), dimension(pNumMax) :: fh, fg, ft
    real(kind=dbPc), dimension(mkxNode, mky, npdf) :: bedPDist
    ! define in particleCal
    real(kind=dbPc), dimension(mkxNode, mky) :: Dkz
    real(kind=dbPc), dimension(mxNode, my, mz) :: phirho
    real(kind=dbPc), dimension(mkxNode, mky, npdf) :: DbedPDist
    !
    real(kind=dbPc) :: norm_vpin, norm_vpout, vvpin, vvpout, mpin, mpout
    real(kind=dbPc) :: npin, npout
    real(kind=dbPc) :: tot_nvpin, tot_nvpout, tot_vvpin, tot_vvpout, tot_mpin, tot_mpout
    real(kind=dbPc) :: tot_npin, tot_npout
    real(kind=dbPc) :: utaot, taot
    real(kind=dbPc) :: aucreep
    real(kind=dbPc), dimension(3) :: vpin, vpout
    real(kind=dbPc), dimension(3) :: tot_vpin, tot_vpout
    real(kind=dbPc), dimension(mz) :: ampu, ampd
    real(kind=dbPc), dimension(mz) :: fptx
    real(kind=dbPc), dimension(mz) :: totfptx, totvolpar
    real(kind=dbPc) :: uflx, wflx
    real(kind=dbPc), dimension(mz) :: uflxz
    real(kind=dbPc), dimension(mz) :: wflxz
    real(kind=dbPc), dimension(mky) :: eepnch, eepnchr
    real(kind=dbPc), dimension(mky*npdf) :: eepdfch, eepdfchr
    real(kind=dbPc) :: time
    integer :: gtypei
    real(kind=dbPc), dimension(mz) :: htaop, thtaop
    real(kind=dbPc), dimension(mz) :: htao, thtao
    real(kind=dbPc), dimension(mz) :: ahff, tahff
    real(kind=dbPc), dimension(mz) :: hru, thru
    real(kind=dbPc), dimension(mz) :: pcoll, ttpcoll
    real(kind=dbPc), dimension(mkxNode, mky) :: ucreep, vcreep
end module public_val

module prob_distribution
    implicit none
contains
    ! need to check

    function biDist(pp)
        implicit none
        ! public
        integer :: biDist
        real(kind=dbPc), intent(in) :: pp
        ! local
        real(kind=dbPc) :: rand
        !
        biDist = 0
        call random_number(rand)
        if (rand >= pp) biDist = 1
    end function biDist
end module prob_distribution

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
        real(kind=dbPc), intent(in), dimension(subMx, totMy, totMz) :: subField
        real(kind=dbPc), dimension(totMx, totMy, totMz) :: totField
        ! local
        integer :: subNx, totNx
        integer :: dataType
        integer :: subGridNum, totGridNum
        integer :: ierr
        integer :: i, j, k
        integer :: serialNum
        integer :: nodeID
        real(kind=dbPc), allocatable, dimension(:, :, :) :: subField_noGhost
        real(kind=dbPc), allocatable, dimension(:) :: subArray
        real(kind=dbPc), allocatable, dimension(:) :: totArray
        !
        subNx = subMx - 2
        totNx = totMx - 2
        subGridNum = subNx*totMy*totMz
        totGridNum = totNx*totMy*totMz
        allocate (subField_noGhost(subNx, totMy, totMz))
        allocate (subArray(subGridNum))
        allocate (totArray(totGridNum))
        !
        dataType = MPI_DOUBLE_PRECISION
        subArray = 0.0
        totArray = 0.0
        totField = 0.0
        subField_noGhost = subField(2:subMx - 1, 1:totMy, 1:totMz)
        do k = 1, totMz
            do j = 1, totMy
                do i = 1, subNx
                    serialNum = i + (j - 1)*subNx + (k - 1)*subNx*totMy
                    subArray(serialNum) = subField_noGhost(i, j, k)
                end do
            end do
        end do
        call MPI_ALLGATHER(subArray, subGridNum, dataType, totArray, subGridNum, dataType, MPI_Comm, ierr)
        do k = 1, totMz
            do j = 1, totMy
                do i = 2, totMx - 1
                    nodeID = (i - 2)/subNx
                    serialNum = (i - 1) - nodeID*subNx + (j - 1)*subNx + (k - 1)*subNx*totMy + nodeID*subNx*totMy*totMz
                    totField(i, j, k) = totArray(serialNum)
                end do
            end do
        end do
        deallocate (subField_noGhost)
        deallocate (subArray)
        deallocate (totArray)
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
        allocate (subField_noGhost(subNx, totMy, totMz))
        allocate (subArray(subGridNum))
        allocate (totArray(totGridNum))
        !
        dataType = MPI_INTEGER
        subArray = 0.0
        totArray = 0.0
        totField = 0.0
        subField_noGhost = subField(2:subMx - 1, 1:totMy, 1:totMz)
        do k = 1, totMz
            do j = 1, totMy
                do i = 1, subNx
                    serialNum = i + (j - 1)*subNx + (k - 1)*subNx*totMy
                    subArray(serialNum) = subField_noGhost(i, j, k)
                end do
            end do
        end do
        call MPI_ALLGATHER(subArray, subGridNum, dataType, totArray, subGridNum, dataType, MPI_Comm, ierr)
        do k = 1, totMz
            do j = 1, totMy
                do i = 2, totMx - 1
                    nodeID = (i - 2)/subNx
                    serialNum = (i - 1) - nodeID*subNx + (j - 1)*subNx + (k - 1)*subNx*totMy + nodeID*subNx*totMy*totMz
                    totField(i, j, k) = totArray(serialNum)
                end do
            end do
        end do
        deallocate (subField_noGhost)
        deallocate (subArray)
        deallocate (totArray)
    end subroutine gxyz_int
    !
    subroutine gatherx(MPI_Comm, subMx, totMx, subField, totField)
        implicit none
        include "mpif.h"
        ! public
        integer, intent(in) :: subMx, totMx
        integer, intent(in) :: MPI_Comm
        real(kind=dbPc), intent(in), dimension(subMx) :: subField
        real(kind=dbPc), dimension(totMx) :: totField
        ! local
        integer :: dataType
        integer :: subGridNum, totGridNum
        integer :: ierr
        real(kind=dbPc), allocatable, dimension(:) :: subArray
        real(kind=dbPc), allocatable, dimension(:) :: totArray
        !
        subGridNum = subMx - 2
        totGridNum = totMx - 2
        allocate (subArray(subGridNum))
        allocate (totArray(totGridNum))
        !
        dataType = MPI_DOUBLE_PRECISION
        totArray = 0.0
        totField = 0.0
        subArray = subField(2:subMx - 1)
        call MPI_ALLGATHER(subArray, subGridNum, dataType, totArray, subGridNum, dataType, MPI_Comm, ierr)
        totField(2:totMx - 1) = totArray
        deallocate (subArray)
        deallocate (totArray)
    end subroutine gatherx
end module gather_xyz

subroutine dataExType
    use public_val
    implicit none
    include "mpif.h"
    integer i, tempExType, ierr
    !
    ! fieldExType: i=const planes
    ! tempExType: i=const,k=const line
    call MPI_TYPE_VECTOR(my, 1, mxNode, realType, tempExType, ierr)
    call MPI_TYPE_COMMIT(tempExType, ierr)
    call MPI_TYPE_EXTENT(realType, i, ierr)
    call MPI_TYPE_HVECTOR(mz, 1, mxNode*my*i, tempExType, fieldExType, ierr)
    call MPI_TYPE_COMMIT(fieldExType, ierr)
    ! surfExType: i=const line
    call MPI_TYPE_VECTOR(mky, 1, mkxNode, realType, surfExType, ierr)
    call MPI_TYPE_COMMIT(surfExType, ierr)
end subroutine dataExType

subroutine dataEx(mxNode, my, mz, a, comm1D, neighbor, fieldExType, tag)
    implicit none
    include "mpif.h"
    ! public
    integer, intent(in) :: mxNode, my, mz
    integer, intent(in) :: comm1D
    integer, intent(in) :: fieldExType
    integer, intent(in) :: tag
    integer, intent(in), dimension(2) :: neighbor
    real(kind=dbPc), dimension(mxNode, my, mz) :: a
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
    call MPI_SENDRECV(a(mxNode - 1, 1, 1), 1, fieldExType, neighbor(2), tag, &
                      a(1, 1, 1), 1, fieldExType, neighbor(1), tag, comm1D, status, ierr)
    ! send to 1 and receive from 2
    call MPI_SENDRECV(a(2, 1, 1), 1, fieldExType, neighbor(1), tag + 1, &
                      a(mxNode, 1, 1), 1, fieldExType, neighbor(2), tag + 1, comm1D, status, ierr)
end subroutine dataEx

subroutine particleCal
    use public_val
    use vector_cal
    use prob_distribution
    implicit none
    include "mpif.h"
    ! MPI
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierr
    ! roll direction
    real(kind=dbPc) :: kzC, kzE, kzW, kzN, kzS
    real(kind=dbPc) :: rampMax
    integer, dimension(mkxNode, mky) :: rollDirBump, rollDirSink
    real(kind=dbPc), dimension(4) :: ramp
    ! splash
    integer :: whichSurfTri
    real(kind=dbPc) :: dotP1, dotP2, dotP3, dotP4, dotP5, dotP6
    real(kind=dbPc) :: lam22, lam23, lam2, lam3
    real(kind=dbPc), dimension(3) :: point1, point2, point3
    real(kind=dbPc), dimension(3) :: vec12, vec13, vec1p, vec2p, vec3p
    real(kind=dbPc), dimension(3) :: surfNormal, unitSurfNormal
    real(kind=dbPc), dimension(3) :: pLocVec, pVelVec
    real(kind=dbPc), dimension(3) :: point0
    real(kind=dbPc), dimension(3) :: vector1, vector2, vector3
    !
    integer :: ik, jk
    integer :: i, j, k
    real(kind=dbPc), dimension(mky) :: kz3Send, kz3Recv
    !
    integer :: nn, nni
    integer :: n
    integer :: ne
    integer :: nrol
    integer :: kd
    integer :: nks, tnk, nnks
    integer :: biDist
    integer :: nnn, nnni
    integer :: pNumTemp
    integer :: pNumAdd
    integer :: nksmax, chmax
    integer :: ii, jj, kk
    integer :: ipp, jpp
    integer :: jjkk
    integer :: iii, jjj
    integer :: nne
    integer :: nbi
    integer :: kkk
    integer :: hpl
    integer :: nkk(mxNode, my)
    real(kind=dbPc) :: kArea
    real(kind=dbPc) :: mp, volp
    real(kind=dbPc) :: norm_vin
    real(kind=dbPc) :: mmu
    real(kind=dbPc) :: sigma
    real(kind=dbPc) :: alpha, beta
    real(kind=dbPc) :: alpha1, beta1
    real(kind=dbPc) :: alpha2, beta2
    real(kind=dbPc) :: eta1, eta2
    real(kind=dbPc) :: rr1, rr2, rr3
    real(kind=dbPc) :: normal
    real(kind=dbPc) :: rrr, rrrx, rrrd
    real(kind=dbPc) :: nv1, nv2, tv1, tv2
    real(kind=dbPc) :: pp
    real(kind=dbPc) :: lambda
    real(kind=dbPc) :: ammu1, ammu2
    real(kind=dbPc) :: gg1, gg2, gg3
    real(kind=dbPc) :: dzz
    real(kind=dbPc) :: aSigma1, aSigma2
    real(kind=dbPc) :: angout1, angout2
    real(kind=dbPc) :: norm_vout
    real(kind=dbPc) :: expdev
    real(kind=dbPc) :: dd1, dd2
    real(kind=dbPc) :: angin1
    real(kind=dbPc) :: ee1, ee2, eed2, eed2x
    real(kind=dbPc) :: mm1, mm2
    real(kind=dbPc) :: d1, d2
    real(kind=dbPc) :: myerfc
    real(kind=dbPc) :: merfc
    real(kind=dbPc) :: ee2bar
    real(kind=dbPc) :: prebound
    real(kind=dbPc) :: arebound
    real(kind=dbPc) :: nee
    real(kind=dbPc) :: vch
    real(kind=dbPc) :: dpd
    real(kind=dbPc) :: svlayer, xvlayer
    real(kind=dbPc) :: norm_nnvec, norm_ttvec, norm_tnvec
    real(kind=dbPc) :: d01, d02, d03
    real(kind=dbPc) :: mmp, nump
    real(kind=dbPc) :: hrlx1, dhcld
    real(kind=dbPc) :: xp0
    real(kind=dbPc), dimension(3) :: nnvec, ttvec, tnvec
    real(kind=dbPc), dimension(3) :: unnvec, uttvec, utnvec
    real(kind=dbPc), dimension(3) :: upvec1, upvec2, xpvec1, xpvec2
    real(kind=dbPc), dimension(3) :: vin
    real(kind=dbPc), dimension(3) :: vout
    real(kind=dbPc), dimension(3) :: gg
    real(kind=dbPc), dimension(npdf) :: ppdf
    real(kind=dbPc), dimension(mky*npdf) :: ppdfch, ppdfchr
    real(kind=dbPc), dimension(mky*npdf) :: ipdfch, ipdfchr
    real(kind=dbPc), dimension(mky) :: pzpdfch, pzpdfchr
    real(kind=dbPc), dimension(mky) :: izpdfch, izpdfchr
    real(kind=dbPc), dimension(mkxNode, mky) :: vlayer
    real(kind=dbPc), dimension(mkxNode, mky, npdf) :: vbin
    real(kind=dbPc), dimension(mkxNode, mky) :: tpdf
    real(kind=dbPc), dimension(mxNode, my, mz) :: vpar
    real(kind=dbPc), dimension(mz) :: ncoll, ntotc
    real(kind=dbPc), dimension(nspmax) :: tempx, tempy, tempz, tempu, tempv, tempw, tempd, tempfk, tempfz, &
                                          tempfh, tempfg, tempft
    real(kind=dbPc), dimension(mkxNode, mky) :: mupin, mvpin, mupout, mvpout, qcreepx, qcreepy
    integer, allocatable, dimension(:, :, :) :: scp
    real(kind=dbPc), allocatable, dimension(:) :: chxp, chyp, chzp
    real(kind=dbPc), allocatable, dimension(:) :: chup, chvp, chwp
    real(kind=dbPc), allocatable, dimension(:) :: chdp, chfk, chfz
    real(kind=dbPc), allocatable, dimension(:) :: chfh, chfg, chft
    real(kind=dbPc), allocatable, dimension(:) :: chxpi, chypi, chzpi
    real(kind=dbPc), allocatable, dimension(:) :: chupi, chvpi, chwpi
    real(kind=dbPc), allocatable, dimension(:) :: chdpi, chfki, chfzi
    real(kind=dbPc), allocatable, dimension(:) :: chfhi, chfgi, chfti
    real(kind=dbPc), allocatable, dimension(:) :: exch, exchi
    real(kind=dbPc), allocatable, dimension(:) :: exchr, exchir
    nksmax = pNumMax !/(mxNode) !/(my)
    chmax = pNumMax/10
    allocate (scp(mxNode, my, nksmax))
    allocate (chxp(chmax), chyp(chmax), chzp(chmax))
    allocate (chup(chmax), chvp(chmax), chwp(chmax))
    allocate (chdp(chmax), chfk(chmax), chfz(chmax))
    allocate (chfh(chmax), chfg(chmax), chft(chmax))
    allocate (chxpi(chmax), chypi(chmax), chzpi(chmax))
    allocate (chupi(chmax), chvpi(chmax), chwpi(chmax))
    allocate (chdpi(chmax), chfki(chmax), chfzi(chmax))
    allocate (chfhi(chmax), chfgi(chmax), chfti(chmax))
    !
    ! rebound and splash
    !
    pNumTemp = 0
    pNumAdd = 0
    Dkz = 0.0
    kArea = kxDiff*kyDiff
    DbedPDist = 0.0
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
    if (pNumAdd >= 1) then
        pNum = pNumTemp + pNumAdd
        if (pNum > pNumMax) then
            print *, pNum, pNumMax
            print *, "particle number reach the threshold"
            stop
        else
            xp(pNumTemp + 1:pNumTemp + pNumAdd) = tempx(1:pNumAdd)
            yp(pNumTemp + 1:pNumTemp + pNumAdd) = tempy(1:pNumAdd)
            zp(pNumTemp + 1:pNumTemp + pNumAdd) = tempz(1:pNumAdd)
            up(pNumTemp + 1:pNumTemp + pNumAdd) = tempu(1:pNumAdd)
            vp(pNumTemp + 1:pNumTemp + pNumAdd) = tempv(1:pNumAdd)
            wp(pNumTemp + 1:pNumTemp + pNumAdd) = tempw(1:pNumAdd)
            dp(pNumTemp + 1:pNumTemp + pNumAdd) = tempd(1:pNumAdd)
            fk(pNumTemp + 1:pNumTemp + pNumAdd) = tempfk(1:pNumAdd)
            fh(pNumTemp + 1:pNumTemp + pNumAdd) = tempfh(1:pNumAdd)
            fg(pNumTemp + 1:pNumTemp + pNumAdd) = tempfg(1:pNumAdd)
            ft(pNumTemp + 1:pNumTemp + pNumAdd) = tempft(1:pNumAdd)
            fz(pNumTemp + 1:pNumTemp + pNumAdd) = tempfz(1:pNumAdd)
        end if
    else
        pNum = pNumTemp
    end if
    !
    do n = 1, pNum
        dpd = dp(n)
        volp = (pi*dpd**3)/6.0
        mp = rhos*volp
        point0(3) = zp(n) - fz(n)
        if (zp(n) < zUni) then
            hpl = 1
            hrlx1 = (zUni - point0(3))/(zUni - initAveKz)
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
        if (wp(n) > 0.0) then
            wflxz(kk) = wflxz(kk) + nkl*mp*wp(n)/xMax/yMax/zDiff(kk)
        end if
        ! up, xp development
        call parvol(up(n), vp(n), wp(n), xp(n), yp(n), fz(n), dpd, mp, k, kk, hrlx1, dhcld, fk(n))
        zp(n) = fz(n) + point0(3)
        if (fz(n) > fh(n)) fh(n) = fz(n)
        if (zp(n) > fg(n)) fg(n) = zp(n)
        ft(n) = ft(n) + dt
    end do
    !
    !aucreep = 0.0
    !qcreepx = 0.0
    !qcreepy = 0.0
    !do jk = 2, mky-1
    !do ik = 2, mkxNode-1
    !dzz = kz(ik, jk) - kz(ik-1, jk)
    !gg1 = 9.8*dzz/sqrt(dzz**2+kxDiff**2)
    !gg2 = 9.8*kxDiff/sqrt(dzz**2+kxDiff**2)
    !qcreepx(ik, jk) = dt*(mupin(ik, jk) - dt*kxDiff*kyDiff*dpa*por*rhos*(gg1+gg2*0.5) - mupout(ik, jk))/kxDiff/rhos/por
    !if (qcreepx(ik, jk)<0.0) qcreepx(ik, jk) = 0.0
    !dzz = kz(ik, jk) - kz(ik, jk-1)
    !gg1 = 9.8*dzz/sqrt(dzz**2+kyDiff**2)
    !gg2 = 9.8*kyDiff/sqrt(dzz**2+kyDiff**2)
    !qcreepy(ik, jk) = dt*(mvpin(ik, jk) - dt*kxDiff*kyDiff*dpa*por*rhos*(gg1+gg2*0.5) - mvpout(ik, jk))/kyDiff/rhos/por
    !if (qcreepy(ik, jk)<0.0) qcreepy(ik, jk) = 0.0
    !ucreep(ik, jk) = qcreepx(ik, jk)/dt/kyDiff/dpa
    !vcreep(ik, jk) = qcreepy(ik, jk)/dt/kxDiff/dpa
    !aucreep = aucreep + ucreep(ik, jk)/dfloat((mkxNode-2)*(mky-2))
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
    !do jk = 2, mky-1
    !do ik = 2, mkxNode-1
    !vch = qcreepx(ik-1, jk) - qcreepx(ik, jk) + qcreepy(ik, jk-1) - qcreepy(ik, jk)
    !Dkz(ik, jk) = Dkz(ik, jk) + vch/kArea/por
    !end do
    !end do
    call surfexch
    ! calculate fluid volume friction phirho
    do k = 1, mz
        do j = 1, my
            do i = 1, mxNode
                phirho(i, j, k) = phirho(i, j, k) - vpar(i, j, k)/(xDiff*yDiff*zDiff(k))
                if (phirho(i, j, k) < 1.0 - por) phirho(i, j, k) = 1.0 - por + 1.0e-6
            end do
        end do
    end do
    ! boundary condition of phirho
    call dataEx(mxNode, my, mz, phirho, comm1D, neighbor, fieldExType, 1)
    do i = 1, mxNode
        do k = 1, mz
            phirho(i, 1, k) = phirho(i, my - 1, k)
            phirho(i, my, k) = phirho(i, 2, k)
        end do
    end do
    ! mid-air collision
    ncoll = 0.0
    ntotc = 0.0
    if (icol == 1) then
        nkk = 0
        do n = 1, pNum
            i = int((xp(n) - xu(1))/xDiff) + 1
            j = int((yp(n) - yv(1))/yDiff) + 1
            if (i < 1) i = 1
            if (j < 1) j = 1
            if (i > mxNode) i = mxNode
            if (j > my) j = my
            nkk(i, j) = nkk(i, j) + 1
            nks = nkk(i, j)
            scp(i, j, nks) = n
        end do
        !
        do j = 1, my
            do i = 1, mxNode
                tnk = nkk(i, j)
                if (tnk < 2) cycle
                do nks = 1, tnk - 1
                    n = scp(i, j, nks)
                    kk = mz
                    if (fz(n) <= 0.0) then
                        kk = 1
                    else
                        do k = 1, mz
                            if (fz(n) >= zw(k) .and. fz(n) < zw(k + 1)) then
                                kk = k
                                exit
                            end if
                        end do
                    end if
                    do nnks = nks + 1, tnk
                        ntotc(kk) = ntotc(kk) + 1.0
                        nn = scp(i, j, nnks)
                        xpvec1(1) = xp(n)
                        xpvec1(2) = yp(n)
                        xpvec1(3) = zp(n)
                        xpvec2(1) = xp(nn)
                        xpvec2(2) = yp(nn)
                        xpvec2(3) = zp(nn)
                        rrr = dist_p(xpvec1, xpvec2)
                        rrrd = (dp(n) + dp(nn))/2.0
                        if (rrr <= rrrd .and. rrr > 0.0) then
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
                            ! vector1=v12 (v1-v2), vector2=v21 (v2-v1)
                            vector1 = upvec1 - upvec2
                            vector2 = -vector1
                            ! nv1=n.v12, nv2=n.v21
                            nv1 = dotProduct(unnvec, vector1)
                            nv2 = nv1
                     tv1 = sqrt((vector1(1) - nv1*unnvec(1))**2 + (vector1(2) - nv1*unnvec(2))**2 + (vector1(3) - nv1*unnvec(3))**2)
                            tv2 = tv1
                            mm1 = (pi*dp(n)**3)/6.0*rhos
                            mm2 = (pi*dp(nn)**3)/6.0*rhos
                            eta1 = mm1/mm2
                            eta2 = mm2/mm1
                            alpha1 = (1.0 + els1)/(1.0 + eta1)
                            alpha2 = (1.0 + els1)/(1.0 + eta2)
                            !if (tv1/nv1<7.0/2.0*fric1*(1.0+els1)) then
                            beta1 = (2.0/7.0)/(1.0 + eta1)
                            beta2 = (2.0/7.0)/(1.0 + eta2)
                            !else
                            !  beta1 = fric1*(1.0+els1)*nv1/tv1/(1.0+eta1)
                            !  beta2 = fric1*(1.0+els1)*nv2/tv2/(1.0+eta2)
                            !end if
                            !
                            up(n) = upvec1(1) - alpha1*nv1*unnvec(1) - beta1*(vector1(1) - nv1*unnvec(1))
                            vp(n) = upvec1(2) - alpha1*nv1*unnvec(2) - beta1*(vector1(2) - nv1*unnvec(2))
                            wp(n) = upvec1(3) - alpha1*nv1*unnvec(3) - beta1*(vector1(3) - nv1*unnvec(3))
                            up(nn) = upvec2(1) - alpha2*nv2*uttvec(1) - beta2*(vector2(1) - nv2*uttvec(1))
                            vp(nn) = upvec2(2) - alpha2*nv2*uttvec(2) - beta2*(vector2(2) - nv2*uttvec(2))
                            wp(nn) = upvec2(3) - alpha2*nv2*uttvec(3) - beta2*(vector2(3) - nv2*uttvec(3))
                            upvec1(1) = up(n)
                            upvec1(2) = vp(n)
                            upvec1(3) = wp(n)
                            upvec2(1) = up(nn)
                            upvec2(2) = vp(nn)
                            upvec2(3) = wp(nn)
                            nv1 = dotProduct(unnvec, upvec1)
                            nv2 = dotProduct(unnvec, upvec2)
                            if (nv1 >= nv2) then
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
            if (ntotc(k) > 0.9) then
                pcoll(k) = ncoll(k)/ntotc(k)
            else
                pcoll(k) = 0.0
            end if
        end do
    end if
    ! pick out particles out of boundary
    pNumTemp = 0
    nn = 0
    nni = 0
    pick: do n = 1, pNum
        if (xp(n) >= xu(mxNode)) then
            nn = nn + 1
            if (myID == nNodes - 1) then
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
        else if (xp(n) < xu(2)) then
            nni = nni + 1
            if (xp(n) < 0.) then
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
            pNumTemp = pNumTemp + 1
            xp(pNumTemp) = xp(n)
            yp(pNumTemp) = yp(n)
            zp(pNumTemp) = zp(n)
            up(pNumTemp) = up(n)
            vp(pNumTemp) = vp(n)
            wp(pNumTemp) = wp(n)
            dp(pNumTemp) = dp(n)
            fk(pNumTemp) = fk(n)
            fz(pNumTemp) = fz(n)
            fh(pNumTemp) = fh(n)
            fg(pNumTemp) = fg(n)
            ft(pNumTemp) = ft(n)
        end if
    end do pick
    pNum = pNumTemp
    ! particle exchange between processes
    call MPI_SENDRECV(nn, 1, intType, neighbor(2), 18, &
                      nnn, 1, intType, neighbor(1), 18, comm1D, status, ierr)
    call MPI_SENDRECV(nni, 1, intType, neighbor(1), 27, &
                      nnni, 1, intType, neighbor(2), 27, comm1D, status, ierr)
    ! from 1 to 2
    allocate (exch(12*nn))
    allocate (exchr(12*nnn))
    exch(1:nn) = chxp(1:nn)
    exch(1*nn + 1:2*nn) = chyp(1:nn)
    exch(2*nn + 1:3*nn) = chzp(1:nn)
    exch(3*nn + 1:4*nn) = chup(1:nn)
    exch(4*nn + 1:5*nn) = chvp(1:nn)
    exch(5*nn + 1:6*nn) = chwp(1:nn)
    exch(6*nn + 1:7*nn) = chdp(1:nn)
    exch(7*nn + 1:8*nn) = chfk(1:nn)
    exch(8*nn + 1:9*nn) = chfz(1:nn)
    exch(9*nn + 1:10*nn) = chfh(1:nn)
    exch(10*nn + 1:11*nn) = chfg(1:nn)
    exch(11*nn + 1:12*nn) = chft(1:nn)
    call MPI_SENDRECV(exch, nn*12, realType, neighbor(2), 10, &
                      exchr, nnn*12, realType, neighbor(1), 10, &
                      comm1D, status, ierr)
    if (nnn > 0) then
        xp(pNum + 1:pNum + nnn) = exchr(1:nnn)
        yp(pNum + 1:pNum + nnn) = exchr(1*nnn + 1:2*nnn)
        zp(pNum + 1:pNum + nnn) = exchr(2*nnn + 1:3*nnn)
        up(pNum + 1:pNum + nnn) = exchr(3*nnn + 1:4*nnn)
        vp(pNum + 1:pNum + nnn) = exchr(4*nnn + 1:5*nnn)
        wp(pNum + 1:pNum + nnn) = exchr(5*nnn + 1:6*nnn)
        dp(pNum + 1:pNum + nnn) = exchr(6*nnn + 1:7*nnn)
        fk(pNum + 1:pNum + nnn) = exchr(7*nnn + 1:8*nnn)
        fz(pNum + 1:pNum + nnn) = exchr(8*nnn + 1:9*nnn)
        fh(pNum + 1:pNum + nnn) = exchr(9*nnn + 1:10*nnn)
        fg(pNum + 1:pNum + nnn) = exchr(10*nnn + 1:11*nnn)
        ft(pNum + 1:pNum + nnn) = exchr(11*nnn + 1:12*nnn)
        if (myID == 0) then
            if (ikbx == 0) then
                pNum = pNum + nnn
            end if
        else
            pNum = pNum + nnn
        end if
    end if
    deallocate (exch)
    deallocate (exchr)
    ! from 2 to 1
    allocate (exchi(12*nni))
    allocate (exchir(12*nnni))
    exchi(1:nni) = chxpi(1:nni)
    exchi(1*nni + 1:2*nni) = chypi(1:nni)
    exchi(2*nni + 1:3*nni) = chzpi(1:nni)
    exchi(3*nni + 1:4*nni) = chupi(1:nni)
    exchi(4*nni + 1:5*nni) = chvpi(1:nni)
    exchi(5*nni + 1:6*nni) = chwpi(1:nni)
    exchi(6*nni + 1:7*nni) = chdpi(1:nni)
    exchi(7*nni + 1:8*nni) = chfki(1:nni)
    exchi(8*nni + 1:9*nni) = chfzi(1:nni)
    exchi(9*nni + 1:10*nni) = chfhi(1:nni)
    exchi(10*nni + 1:11*nni) = chfgi(1:nni)
    exchi(11*nni + 1:12*nni) = chfti(1:nni)
    call MPI_SENDRECV(exchi, nni*12, realType, neighbor(1), 19, &
                      exchir, nnni*12, realType, neighbor(2), 19, &
                      comm1D, status, ierr)
    if (nnni > 0) then
        xp(pNum + 1:pNum + nnni) = exchir(1:nnni)
        yp(pNum + 1:pNum + nnni) = exchir(1*nnni + 1:2*nnni)
        zp(pNum + 1:pNum + nnni) = exchir(2*nnni + 1:3*nnni)
        up(pNum + 1:pNum + nnni) = exchir(3*nnni + 1:4*nnni)
        vp(pNum + 1:pNum + nnni) = exchir(4*nnni + 1:5*nnni)
        wp(pNum + 1:pNum + nnni) = exchir(5*nnni + 1:6*nnni)
        dp(pNum + 1:pNum + nnni) = exchir(6*nnni + 1:7*nnni)
        fk(pNum + 1:pNum + nnni) = exchir(7*nnni + 1:8*nnni)
        fz(pNum + 1:pNum + nnni) = exchir(8*nnni + 1:9*nnni)
        fh(pNum + 1:pNum + nnni) = exchir(9*nnni + 1:10*nnni)
        fg(pNum + 1:pNum + nnni) = exchir(10*nnni + 1:11*nnni)
        ft(pNum + 1:pNum + nnni) = exchir(11*nnni + 1:12*nnni)
        if (myID == nNodes - 1) then
            if (ikbx == 0) then
                pNum = pNum + nnni
            end if
        else
            pNum = pNum + nnni
        end if
    end if
    deallocate (exchi)
    deallocate (exchir)
    ! y z boundary condition of particles
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
    ! bed particle size distributation
    if (ipd /= 1 .and. irsf == 0) then
        ! vlayer change
        do j = 2, mky - 1
            do i = 2, mkxNode - 1
                vlayer(i, j) = kxDiff*kyDiff*bedCellTkness(i, j)*por
                do k = 1, npdf
                    vbin(i, j, k) = bedPDist(i, j, k)*vlayer(i, j) + DbedPDist(i, j, k)
                    vlayer(i, j) = vlayer(i, j) + DbedPDist(i, j, k)
                end do
            end do
        end do
        ! calculate bedPDist, bedCellTkness
        do j = 2, mky - 1
            do i = 2, mkxNode - 1
                bedCellTkness(i, j) = vlayer(i, j)/(kxDiff*kyDiff)/por
                if (bedCellTkness(i, j) <= 0.0) then
                    print *, 'bedCellTkness<0'
                    bedCellTkness(i, j) = bedCellTknessInit
                    do k = 1, npdf
                        bedPDist(i, j, k) = prob(k)
                    end do
                else if (bedCellTkness(i, j) >= bedCellTknessInit) then
                    do k = 1, npdf
                        bedPDist(i, j, k) = vbin(i, j, k)/vlayer(i, j)
                    end do
                    if (bedCellTkness(i, j) > 2.0*bedCellTknessInit) bedCellTkness(i, j) = 2.0*bedCellTknessInit
                else
                    svlayer = kxDiff*kyDiff*bedCellTknessInit*por
                    xvlayer = kxDiff*kyDiff*(bedCellTknessInit - bedCellTkness(i, j))*por
                    do k = 1, npdf
                        bedPDist(i, j, k) = (vbin(i, j, k) + xvlayer*prob(k))/svlayer
                    end do
                    bedCellTkness(i, j) = bedCellTknessInit
                end if
            end do
        end do
        ! bedPDist, bedCellTkness exchange
        ! i=mkxNode-1 send to i=1: send to 2 and receive from 1
        do j = 1, mky
            pzpdfch(j) = bedCellTkness(mkxNode - 1, j)
            do k = 1, npdf
                jjkk = k + (j - 1)*npdf
                ppdfch(jjkk) = bedPDist(mkxNode - 1, j, k)
            end do
        end do
        call MPI_SENDRECV(ppdfch, mky*npdf, realType, neighbor(2), 201, &
                          ppdfchr, mky*npdf, realType, neighbor(1), 201, comm1D, status, ierr)
        call MPI_SENDRECV(pzpdfch, mky, realType, neighbor(2), 203, &
                          pzpdfchr, mky, realType, neighbor(1), 203, comm1D, status, ierr)
        do j = 1, mky
            bedCellTkness(1, j) = pzpdfchr(j)
            do k = 1, npdf
                jjkk = k + (j - 1)*npdf
                bedPDist(1, j, k) = ppdfchr(jjkk)
            end do
        end do
        ! i=2 send to i=mkxNode: send to 1 and receive from 2
        do j = 1, mky
            izpdfch(j) = bedCellTkness(2, j)
            do k = 1, npdf
                jjkk = k + (j - 1)*npdf
                ipdfch(jjkk) = bedPDist(2, j, k)
            end do
        end do
        call MPI_SENDRECV(ipdfch, mky*npdf, realType, neighbor(1), 202, &
                          ipdfchr, mky*npdf, realType, neighbor(2), 202, comm1D, status, ierr)
        call MPI_SENDRECV(izpdfch, mky, realType, neighbor(1), 204, &
                          izpdfchr, mky, realType, neighbor(2), 204, comm1D, status, ierr)
        do j = 1, mky
            bedCellTkness(mkxNode, j) = izpdfchr(j)
            do k = 1, npdf
                jjkk = k + (j - 1)*npdf
                bedPDist(mkxNode, j, k) = ipdfchr(jjkk)
            end do
        end do
        ! y=mky-1 send to y=1, y=2 send to y=mky
        do i = 2, mkxNode - 1
            bedCellTkness(i, 1) = bedCellTkness(i, mky - 1)
            bedCellTkness(i, mky) = bedCellTkness(i, 2)
            do k = 1, npdf
                bedPDist(i, 1, k) = bedPDist(i, mky - 1, k)
                bedPDist(i, mky, k) = bedPDist(i, 2, k)
            end do
        end do
        !
        tpdf = 0.0
        do i = 1, mkxNode
            do j = 1, mky
                do k = 1, npdf
                    tpdf(i, j) = tpdf(i, j) + bedPDist(i, j, k)
                end do
                if (tpdf(i, j) <= 0.0) then
                    print *, tpdf(i, j), 'tpdf<=0'
                    stop
                end if
            end do
        end do
        !
        bedPD = 0.0
        do j = 1, mky
            do i = 1, mkxNode
                if (ipd == 0) then
                    do k = 1, npdf
                        bedPDist(i, j, k) = bedPDist(i, j, k)/tpdf(i, j)
                        if (bedPDist(i, j, k) < 0.0) then
                            bedPDist(i, j, k) = 0.0
                        end if
                  bedPD(i, j) = bedPD(i, j) + bedPDist(i, j, k)*((dfloat(k - 1) + 0.5)*dSigma*6.0/dfloat(npdf) + (dpa - dSigma*3.0))
                    end do
                else if (ipd == 2) then
                    bedPDist(i, j, 1) = bedPDist(i, j, 1)/tpdf(i, j)
                    bedPDist(i, j, 2) = bedPDist(i, j, 2)/tpdf(i, j)
                    if (bedPDist(i, j, 1) < 0.0) then
                        bedPDist(i, j, 1) = 0.0
                        bedPDist(i, j, 2) = 1.0
                    else if (bedPDist(i, j, 2) < 0.0) then
                        bedPDist(i, j, 2) = 0.0
                        bedPDist(i, j, 1) = 1.0
                    end if
                    bedPD(i, j) = bedPDist(i, j, 1)*(dpa - dSigma) + bedPDist(i, j, 2)*(dpa + dSigma)
              if (bedPDist(i, j, 1) + bedPDist(i, j, 2) > 1.01 .or. bedPD(i, j) > dpa + dSigma .or. bedPD(i, j) < dpa - dSigma) then
                        print *, 'bedPD error', bedPD(i, j), bedPDist(i, j, 1), bedPDist(i, j, 2), bedCellTkness(i, j), &
                            bedCellTknessInit, vlayer(i, j), vbin(i, j, 1), vbin(i, j, 2)
                        stop
                    end if
                end if
                if (bedPD(i, j) > dpa + dSigma*3.0) bedPD(i, j) = dpa + dSigma*3.0
                if (bedPD(i, j) < dpa - dSigma*3.0) bedPD(i, j) = dpa - dSigma*3.0
            end do
        end do
    else
        bedPD = dpa
    end if
    !
    deallocate (scp)
    deallocate (chxp, chyp, chzp)
    deallocate (chup, chvp, chwp)
    deallocate (chdp, chfk, chfz)
    deallocate (chfh, chfg, chft)
    deallocate (chxpi, chypi, chzpi)
    deallocate (chupi, chvpi, chwpi)
    deallocate (chdpi, chfki, chfzi)
    deallocate (chfhi, chfgi, chfti)
end subroutine particleCal

subroutine parloc(iLoc, jLoc, kkk1, kkk2, xxp, yyp, zzp, hhp, hr1, hl, dhcld)
    use public_val
    implicit none
    include "mpif.h"
    ! public
    integer :: iLoc, jLoc, kkk1, kkk2, hl
    real(kind=dbPc), intent(in) :: xxp, yyp, zzp, hhp, hr1
    real(kind=dbPc) :: dhcld
    ! local
    integer :: k
    real(kind=dbPc) :: zzw, zzw0, zzw1, zzw2, rzp
    !
    iLoc = int((xxp - xu(1))/xDiff) + 1
    jLoc = int((yyp - yv(1))/yDiff) + 1
    if (iLoc < 1) iLoc = 1
    if (jLoc < 1) jLoc = 1
    if (iLoc > mxNode - 1) iLoc = mxNode - 1
    if (jLoc > my - 1) jLoc = my - 1
    dhcld = 0.0
    if (hl == 1) then
        if (hhp < zw(1)) then
            kkk2 = 1
            kkk1 = kkk2
            dhcld = 0.5*zDiff(1)*hr1
        else
            do k = 1, mz
                if (k == 1) then
                    zzw = 0.0
                else
                    zzw = zzw + zDiff(k - 1)*hr1
                end if
                zzw0 = zzw + zDiff(k)*hr1
                if (hhp < zzw0) then
                    kkk2 = k
                    zzw1 = zzw + 0.5*zDiff(kkk2)*hr1
                    if (hhp < zzw1) then
                        kkk1 = kkk2
                        dhcld = abs(zzw1 - hhp)
                    else
                        kkk1 = kkk2 + 1
                        zzw2 = zzw0 + 0.5*zDiff(kkk1)*hr1
                        dhcld = abs(zzw2 - hhp)
                    end if
                    exit
                end if
            end do
        end if
    else if (hl == 0) then
        rzp = zzp - initAveKz
        if (rzp >= z(mz)) then
            kkk2 = mz
            kkk1 = kkk2
            dhcld = 0.0
        else
            do k = 1, mz
                if (rzp < zw(k + 1)) then
                    kkk2 = k
                    if (rzp < z(kkk2)) then
                        kkk1 = kkk2
                    else
                        kkk1 = kkk2 + 1
                    end if
                    dhcld = abs(z(kkk1) - rzp)
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
    real(kind=dbPc) :: uup, vvp, wwp
    real(kind=dbPc) :: xxp, yyp, hhp, ffk
    real(kind=dbPc), intent(in) :: ddp, mmp, hrlx1, dhcld
    ! local
    real(kind=dbPc) :: alpha, beta
    real(kind=dbPc) :: ufp, vfp, wfp
    real(kind=dbPc) :: d1x, d1y, d1z
    real(kind=dbPc) :: d1u, d1v, d1w
    real(kind=dbPc) :: d2x, d2y, d2z
    real(kind=dbPc) :: d2u, d2v, d2w
    real(kind=dbPc) :: d3x, d3y, d3z
    real(kind=dbPc) :: d3u, d3v, d3w
    real(kind=dbPc) :: d4x, d4y, d4z
    real(kind=dbPc) :: d4u, d4v, d4w
    real(kind=dbPc) :: xxp0
    ! function
    real(kind=dbPc) :: ffd
    ! ufp, vfp, wfp
    if (k > mz) then
        ufp = hru(mz)
    else if (k == 1) then
        if (hhp <= 0.0) then
            ufp = 0.0
        else
            alpha = hhp/(z(1)*hrlx1)
            ufp = hru(1)*alpha
        end if
    else
        beta = dhcld/(zDiff(k) + zDiff(k - 1))/hrlx1*2.0
        alpha = 1.0 - beta
        ufp = hru(k - 1)*beta + hru(k)*alpha
    end if
    vfp = 0.0
    wfp = 0.0
    ! up, xp development
    d1x = uup
    d1y = vvp
    d1z = wwp
    d1u = ffd(d1x, ufp, ddp, nu, rho, rhos)/mmp
    d1v = ffd(d1y, vfp, ddp, nu, rho, rhos)/mmp
    d1w = ffd(d1z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0 - rho/rhos)
    d2x = uup + dt/2.0*d1u
    d2y = vvp + dt/2.0*d1v
    d2z = wwp + dt/2.0*d1w
    d2u = ffd(d2x, ufp, ddp, nu, rho, rhos)/mmp
    d2v = ffd(d2y, vfp, ddp, nu, rho, rhos)/mmp
    d2w = ffd(d2z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0 - rho/rhos)
    d3x = uup + dt/2.0*d2u
    d3y = vvp + dt/2.0*d2v
    d3z = wwp + dt/2.0*d2w
    d3u = ffd(d3x, ufp, ddp, nu, rho, rhos)/mmp
    d3v = ffd(d3y, vfp, ddp, nu, rho, rhos)/mmp
    d3w = ffd(d3z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0 - rho/rhos)
    d4x = uup + dt*d3u
    d4y = vvp + dt*d3v
    d4z = wwp + dt*d3w
    d4u = ffd(d4x, ufp, ddp, nu, rho, rhos)/mmp
    d4v = ffd(d4y, vfp, ddp, nu, rho, rhos)/mmp
    d4w = ffd(d4z, wfp, ddp, nu, rho, rhos)/mmp - 9.8*(1.0 - rho/rhos)
    xxp0 = xxp
    xxp = xxp + (d1x + 2.0*d2x + 2.0*d3x + d4x)/6.0*dt
    yyp = yyp + (d1y + 2.0*d2y + 2.0*d3y + d4y)/6.0*dt
    hhp = hhp + (d1z + 2.0*d2z + 2.0*d3z + d4z)/6.0*dt
    uup = uup + (d1u + 2.0*d2u + 2.0*d3u + d4u)/6.0*dt
    vvp = vvp + (d1v + 2.0*d2v + 2.0*d3v + d4v)/6.0*dt
    wwp = wwp + (d1w + 2.0*d2w + 2.0*d3w + d4w)/6.0*dt
    ffk = ffk + xxp - xxp0
    ampd(kk) = (d1u + 2.0*d2u + 2.0*d3u + d4u)/6.0*mmp
    ampu(kk) = ampu(kk) + ampd(kk)
end subroutine parvol

subroutine fluidField
    use public_val
    implicit none
    include "mpif.h"
    integer :: i, j, k, n
    integer :: kk
    integer :: h
    integer :: ierr
    real(kind=dbPc) :: mixl
    real(kind=dbPc) :: dudz
    real(kind=dbPc) :: oo
    real(kind=dbPc) :: lmd
    real(kind=dbPc) :: shru
    real(kind=dbPc) :: ddz
    real(kind=dbPc) :: chru
    real(kind=dbPc) :: relax
    real(kind=dbPc), dimension(mz) :: volpar, numpar, tvolpar
    real(kind=dbPc), dimension(mz) :: afptx
    real(kind=dbPc), dimension(mz) :: pfptx
    real(kind=dbPc), dimension(mz) :: tfptx
    real(kind=dbPc), dimension(mz) :: ataop
    real(kind=dbPc), dimension(mz) :: ptaop
    real(kind=dbPc), dimension(mz) :: tampd, tampu, aampd, aampu
    ! function
    real(kind=dbPc) :: ffd, ntmixl
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
        htaop(k) = htaop(k)*(1.0 - relax) + relax*tfptx(k)
        if (ahff(k) <= (1.0 - por)) then
            htaop(k) = rho*wind**2
        end if
        htao(k) = rho*wind**2 - htaop(k)
        if (htao(k) < 0.0) htao(k) = 0.0
        !if (htao(k)>rho*wind**2) htao(k) = rho*wind**2
    end do
    !
    mixl = kapa*z(1)*(1.0 - exp(-1.0/26.0*z(1)*wind/nu))
    dudz = (sqrt(nu**2 + 4.0*mixl**2*htao(1)/rho) - nu)/2.0/mixl**2
    hru(1) = dudz*z(1)
    do k = 2, mz
        ddz = z(k) - z(k - 1)
        mixl = kapa*z(k)*(1.0 - exp(-1.0/26.0*z(k)*wind/nu))
        dudz = (sqrt(nu**2 + 4.0*mixl**2*htao(k)/rho) - nu)/2.0/mixl**2
        hru(k) = dudz*ddz + hru(k - 1)
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
    real(kind=dbPc) :: aaa
    real(kind=dbPc) :: totpz, avepz, tavepz
    real(kind=dbPc) :: posit
    !
    if (last == 1 .or. ipar == 0) then
        kx = 0.
        ky = 0.
        kz = 0.
        Dkz = 0.0
        kxDiff = xMax/dfloat(mkx - 2)
        kyDiff = yMax/dfloat(mky - 2)
        if (myID == 0) then
            kx(1) = -kxDiff
            do i = 2, mkxNode
                kx(i) = kx(i - 1) + kxDiff
            end do
            call MPI_SEND(kx(mkxNode - 1), 1, realType, neighbor(2), 3, comm1D, ierr)
        else
            call MPI_RECV(kx(1), 1, realType, neighbor(1), 3, comm1D, status, ierr)
            do i = 2, mkxNode
                kx(i) = kx(i - 1) + kxDiff
            end do
            call MPI_SEND(kx(mkxNode - 1), 1, realType, neighbor(2), 3, comm1D, ierr)
        end if
        ky(1) = -kyDiff
        do j = 2, mky
            ky(j) = ky(j - 1) + kyDiff
        end do
        do j = 1, mky
            do i = 1, mkxNode
                kz(i, j) = initAmp*sin(dfloat(int(time/20.0) + 2)*8.0*pi*kx(i)) + initAveKz
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
        !kz(i, j) = kz(i, j) + Dkz(i, j)
        !end do
        !end do
        do j = 1, mky
            do i = 1, mkxNode
                kz(i, j) = initAmp*sin(dfloat(int(time/20.0) + 2)*8.0*pi*kx(i)) + initAveKz
            end do
        end do
        !do j = 1, mky
        !do i = 1, mkxNode
        !if (kz(i, j)<0. .or. kz(i, j)>zMax .or. abs(Dkz(i, j))>=0.01) then
        ! print*, 'error: kz out of lower/upper boundary'
        !  print*, 'z=', kz(i, j), '    z change=', Dkz(i, j)
        !  print*, 'i=', i, '  j=', j, '  myID=', myID
        !  kz(i, j) = kz(i, j) - Dkz(i, j)
        !end if
        !end do
        !end do
        call pxch(mkxNode, mky, kz, surfExType, neighbor, comm1D)
        do i = 1, mkxNode
            kz(i, 1) = kz(i, mky - 1)
            kz(i, mky) = kz(i, 2)
        end do
        !
        n = 0
        totpz = 0.0
        do j = 2, mky - 1
            do i = 2, mkxNode - 1
                n = n + 1
                totpz = totpz + kz(i, j)
            end do
        end do
        avepz = totpz/dfloat(n)
        call MPI_ALLREDUCE(avepz, tavepz, 1, realType, MPI_SUM, comm1D, ierr)
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
    real(kind=dbPc), dimension(mkxNode, mky) :: kz
    ! local
    integer :: i, j
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    !
    ! send to 2 and receive from 1
 call MPI_SENDRECV(kz(mkxNode - 1, 1), 1, surfExType, neighbor(2), 1, kz(1, 1), 1, surfExType, neighbor(1), 1, comm1D, status, ierr)
    ! send to 1 and receive from 2
    call MPI_SENDRECV(kz(2, 1), 1, surfExType, neighbor(1), 2, kz(mkxNode, 1), 1, surfExType, neighbor(2), 2, comm1D, status, ierr)
end subroutine pxch

function ntmixl(lm, k, nu, ux, dz)
    implicit none
    ! public
    real(kind=dbPc) :: ntmixl
    real(kind=dbPc), intent(in) :: lm, k, nu, ux, dz
    ! local
    integer :: n
    real(kind=dbPc) :: xr, x1
    real(kind=dbPc) :: fx, gx, fxr
    x1 = 1.0
    n = 0
    do
        n = n + 1
        if (n > 10000) then
            print *, 'ntmixl', lm, k, nu, ux, dz
            stop
        end if
        fx = k*(1.0 - exp(-sqrt(ux*x1/7.0/nu))) - (x1 - lm)/dz
        gx = exp(-sqrt(ux*x1/7.0/nu))*k*ux/(2.0*nu*7.0*sqrt(ux*x1/7.0/nu)) - 1.0/dz
        xr = x1 - fx/gx
        fxr = k*(1.0 - exp(-sqrt(ux*xr/7.0/nu))) - (xr - lm)/dz
        if (abs(fxr) < 1.0e-6) exit
        x1 = xr
    end do
    ntmixl = xr
end function ntmixl

subroutine surfexch
    use public_val
    implicit none
    include "mpif.h"
    ! local
    integer :: i, j, k, jk
    integer :: ierr
    integer :: status(MPI_STATUS_SIZE)
    real(kind=dbPc), dimension(mky) :: epnch, epnchr
    real(kind=dbPc), dimension(mky) :: spnch, spnchr
    real(kind=dbPc), dimension(mky*npdf) :: epdfch, epdfchr
    real(kind=dbPc), dimension(mky*npdf) :: spdfch, spdfchr
    ! because the value of ghost cell has changed
    ! need to add ghost value back to real domain before exchange
    ! Dkz, DbedPDist add back
    ! x=mkxNode+1 add to x=3: send to 2 and receive from 1
    ! x=mkxNode add to x=2: send to 2 and receive from 1
    ! x=1 add to x=mkxNode-1: send to 1 and receive from 2
    do j = 1, mky
        epnch(j) = Dkz(mkxNode, j)
        spnch(j) = Dkz(1, j)
        do k = 1, npdf
            jk = k + (j - 1)*npdf
            epdfch(jk) = DbedPDist(mkxNode, j, k)
            spdfch(jk) = DbedPDist(1, j, k)
        end do
    end do
    !
    call MPI_SENDRECV(eepnch, mky, realType, neighbor(2), 107, &
                      eepnchr, mky, realType, neighbor(1), 107, comm1D, status, ierr)
    call MPI_SENDRECV(eepdfch, mky*npdf, realType, neighbor(2), 106, &
                      eepdfchr, mky*npdf, realType, neighbor(1), 106, comm1D, status, ierr)
    !
    call MPI_SENDRECV(epnch, mky, realType, neighbor(2), 103, &
                      epnchr, mky, realType, neighbor(1), 103, comm1D, status, ierr)
    call MPI_SENDRECV(epdfch, mky*npdf, realType, neighbor(2), 200, &
                      epdfchr, mky*npdf, realType, neighbor(1), 200, comm1D, status, ierr)
    !
    call MPI_SENDRECV(spnch, mky, realType, neighbor(1), 104, &
                      spnchr, mky, realType, neighbor(2), 104, comm1D, status, ierr)
    call MPI_SENDRECV(spdfch, mky*npdf, realType, neighbor(1), 108, &
                      spdfchr, mky*npdf, realType, neighbor(2), 108, comm1D, status, ierr)
    !
    do j = 1, mky
        Dkz(3, j) = Dkz(3, j) + eepnchr(j)
        Dkz(2, j) = Dkz(2, j) + epnchr(j)
        Dkz(mkxNode - 1, j) = Dkz(mkxNode - 1, j) + spnchr(j)
        do k = 1, npdf
            jk = k + (j - 1)*npdf
            DbedPDist(3, j, k) = DbedPDist(3, j, k) + eepdfchr(jk)
            DbedPDist(2, j, k) = DbedPDist(2, j, k) + epdfchr(jk)
            DbedPDist(mkxNode - 1, j, k) = DbedPDist(mkxNode - 1, j, k) + spdfchr(jk)
        end do
    end do
    ! y=1 add to y=mky-1, y=mky add to y=2
    do i = 1, mkxNode
        Dkz(i, 2) = Dkz(i, 2) + Dkz(i, mky)
        Dkz(i, mky - 1) = Dkz(i, mky - 1) + Dkz(i, 1)
        do k = 1, npdf
            DbedPDist(i, 2, k) = DbedPDist(i, 2, k) + DbedPDist(i, mky, k)
            DbedPDist(i, mky - 1, k) = DbedPDist(i, mky - 1, k) + DbedPDist(i, 1, k)
        end do
    end do
    ! Dkz, DbedPDist exchange
    ! x=2 send to x=mkxNode, x=mkxNode-1 send to x=1
    call pxch(mkxNode, mky, Dkz, surfExType, neighbor, comm1D)
    ! y=2 send to y=mky, y=mky-1 send to y=1
    do i = 1, mkxNode
        Dkz(i, mky) = Dkz(i, 2)
        Dkz(i, 1) = Dkz(i, mky - 1)
        do k = 1, npdf
            DbedPDist(i, mky, k) = DbedPDist(i, 2, k)
            DbedPDist(i, 1, k) = DbedPDist(i, mky - 1, k)
        end do
    end do
end subroutine surfexch

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
    real(kind=dbPc) :: tuflx, twflx
    real(kind=dbPc) :: tnorm_vpin, tnorm_vpout, tvvpin, tvvpout, tmpin, tmpout
    real(kind=dbPc) :: tnpin, tnpout
    real(kind=dbPc), dimension(3) :: tvpin, tvpout
    real(kind=dbPc), dimension(mx, my, mz) :: tu, tv, tw, tp
    integer, dimension(mx, my, mz) :: tfp
    real(kind=dbPc), dimension(mx) :: tx
    real(kind=dbPc), dimension(mkx) :: tpx
    real(kind=dbPc), dimension(mkx, mky) :: tpz, tpz4
    real(kind=dbPc), dimension((mkxNode - 2)*mky) :: apz
    real(kind=dbPc), dimension((mkx - 2)*mky) :: tapz
    real(kind=dbPc), dimension((mkxNode - 2)*mky) :: apz4
    real(kind=dbPc), dimension((mkx - 2)*mky) :: tapz4
    real(kind=dbPc), dimension(mz) :: tuflxz, twflxz
    real(kind=dbPc), dimension(mz) :: tpcoll, apcoll
    real(kind=dbPc), allocatable, dimension(:) :: txp, typ, tzp, tdp, tup, tvp, twp, tfk, tfz, tfh, tfg, tft
    !
    nf = mod(last, nnf)
    ns = mod(last, nns)
    nc = mod(last, nnc)
    np = mod(last, nnkl)
    nsf = mod(last, nnsf)
    nfi = mod(last, nnfi)
    nfx = mod(last, nnfx)
    nfii = (last - nfi)/nnfi
    write (ctemp, '(i3)') nfii
    if (nfi == 0) then
        if (myID == 0) then
            open (unit=32, file='./particle_loc/particle_loc'//trim(adjustl(ctemp))//'.plt')
            write (32, "(A82)") 'variables = "XP", "YP", "ZP", "DP", "UP", "VP", "WP", "FK", "FZ", "FH", "FG", "FT"'
            close (32)
            !
            open (unit=33, file='./surface/surface'//trim(adjustl(ctemp))//'.plt')
            write (33, *) 'variables = "X", "Y", "Z", "DP"'
            close (33)
            !
            open (unit=42, file='./field3D/field3D'//trim(adjustl(ctemp))//'.plt')
            write (42, *) 'variables = "X", "Y", "Z", "concentration"'
            close (42)
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
    if (ns == 0) then
        call MPI_ALLREDUCE(ttpcoll, tpcoll, mz, realType, MPI_SUM, comm1D, ierr)
        apcoll = tpcoll/nNodes
        if (myID == 0) then
            open (unit=43, position='append', file='htao.plt')
            write (43, *) 'zone', ' T = "', time, '"'
            do k = 1, mz
                write (43, "(5E15.7)") z(k), thtao(k), thtaop(k), tahff(k), thru(k), apcoll(k)
            end do
            close (43)
        end if
        thtaop = 0.0
        thtao = 0.0
        tahff = 0.0
        thru = 0.0
        ttpcoll = 0.0
    end if
    !
    if (ipar == 1) then
        call MPI_ALLREDUCE(pNum, tnnp, 1, intType, MPI_SUM, comm1D, ierr)
        if (nc == 0) then
            if (myID == 0) then
                open (unit=31, position='append', file='particle_num.plt')
                write (31, "(5E15.7)") time, real(tnnp)
                close (31)
            end if
        end if
    end if
    !
    allocate (txp(tnnp), typ(tnnp), tzp(tnnp), tdp(tnnp), tup(tnnp), tvp(tnnp), twp(tnnp), tfk(tnnp), tfz(tnnp), &
              tfh(tnnp), tfg(tnnp), tft(tnnp))
    if (ipar == 1) then
        if (np == 0) then
            if (last >= pistart) then
                displs(1) = 0
                call MPI_GATHER(pNum, 1, intType, cnt, 1, intType, 0, comm1D, ierr)
                do i = 2, nNodes
                    displs(i) = displs(i - 1) + cnt(i - 1)
                end do
                call MPI_GATHERV(xp, pNum, realType, txp, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(yp, pNum, realType, typ, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(zp, pNum, realType, tzp, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(dp, pNum, realType, tdp, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(up, pNum, realType, tup, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(vp, pNum, realType, tvp, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(wp, pNum, realType, twp, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(fk, pNum, realType, tfk, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(fz, pNum, realType, tfz, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(fh, pNum, realType, tfh, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(fg, pNum, realType, tfg, cnt, displs, realType, 0, comm1D, ierr)
                call MPI_GATHERV(ft, pNum, realType, tft, cnt, displs, realType, 0, comm1D, ierr)
                if (myID == 0) then
                    open (unit=32, position='append', file='./particle_loc/particle_loc'//trim(adjustl(ctemp))//'.plt')
                    write (32, *) 'zone', ' T = "', time, '"'
                    do n = 1, tnnp
                        write (32, "(5E15.7)") txp(n), typ(n), tzp(n), tdp(n), tup(n), tvp(n), &
                            twp(n), tfk(n), tfz(n), tfh(n), tfg(n), tft(n)
                    end do
                    close (32)
                end if
            end if
        end if
    end if
    !
    if (nsf == 0) then
        ! Gather all kx to myID=0
        call gatherx(comm1D, mkxNode, mkx, kx, tpx)
        ! Gather all kz and bedPD to myID=0
        do j = 1, mky
            do i = 2, mkxNode - 1
                ij = (i - 1) + (j - 1)*(mkxNode - 2)
                apz(ij) = kz(i, j)
                apz4(ij) = bedPD(i, j)
            end do
        end do
        numa = (mkxNode - 2)*mky
        call MPI_ALLGATHER(apz, numa, realType, tapz, numa, realType, comm1D, ierr)
        call MPI_ALLGATHER(apz4, numa, realType, tapz4, numa, realType, comm1D, ierr)
        do j = 1, mky
            do i = 2, mkx - 1
                blk = (i - 2)/(mkxNode - 2)
                ij = (i - 1) - blk*(mkxNode - 2) + (j - 1)*(mkxNode - 2) + blk*(mkxNode - 2)*mky
                tpz(i, j) = tapz(ij)
                tpz4(i, j) = tapz4(ij)
            end do
        end do
        if (myID == 0) then
            open (unit=33, position='append', file='./surface/surface'//trim(adjustl(ctemp))//'.plt')
            write (33, *) 'zone', ' T = "', time, '"'
            write (33, *) 'i=', mkx - 2, ' j=', mky - 2, ' datapacking=point'
            do j = 2, mky - 1
                do i = 2, mkx - 1
                    write (33, "(5E15.7)") tpx(i), ky(j), tpz(i, j), tpz4(i, j)
                end do
            end do
            close (33)
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
    call MPI_ALLREDUCE(vpin, tvpin, 3, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(vpout, tvpout, 3, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(norm_vpin, tnorm_vpin, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(norm_vpout, tnorm_vpout, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(vvpin, tvvpin, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(vvpout, tvvpout, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(mpin, tmpin, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(mpout, tmpout, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(npin, tnpin, 1, realType, MPI_SUM, comm1D, ierr)
    call MPI_ALLREDUCE(npout, tnpout, 1, realType, MPI_SUM, comm1D, ierr)
    if (tnpin <= 0.0) then
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
    if (tnpout <= 0.0) then
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
    if (nfx == 0) then
        call MPI_ALLREDUCE(uflx, tuflx, 1, realType, MPI_SUM, comm1D, ierr)
        call MPI_ALLREDUCE(wflx, twflx, 1, realType, MPI_SUM, comm1D, ierr)
        call MPI_ALLREDUCE(uflxz, tuflxz, mz, realType, MPI_SUM, comm1D, ierr)
        call MPI_ALLREDUCE(wflxz, twflxz, mz, realType, MPI_SUM, comm1D, ierr)
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
        if (myID == 0) then
            open (unit=35, position='append', file='average_flux.plt')
            if (twflx /= 0.) then
                write (35, "(5E15.7)") time, tuflx, twflx, tuflx/twflx
            else
                write (35, "(5E15.7)") time, tuflx, twflx, 0.0
            end if
            close (35)
            !
            open (unit=36, position='append', file='flux_vs_height.plt')
            write (36, *) 'zone', ' T = "', time, '"'
            do k = 1, mz
                write (36, "(5E15.7)") z(k), tuflxz(k), twflxz(k)
            end do
            close (36)
            !
            open (unit=39, position='append', file='vin.plt')
            write (39, "(5E15.7)") time, tot_vpin(1), tot_vpin(2), tot_vpin(3), tot_nvpin
            close (39)
            !
            open (unit=46, position='append', file='vout.plt')
            write (46, "(5E15.7)") time, tot_vpout(1), tot_vpout(2), tot_vpout(3), tot_nvpout
            close (46)
            !
            open (unit=44, position='append', file='eminout.plt')
            write (44, "(5E15.7)") time, tot_vvpin, tot_vvpout, tot_mpin, tot_mpout
            close (44)
            !
            open (unit=45, position='append', file='numinout.plt')
            write (45, "(5E15.7)") time, tot_npin, tot_npout
            close (45)
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
    deallocate (txp, typ, tzp, tdp, tup, tvp, twp, tfk, tfz, tfh, tfg, tft)
end subroutine output

function closestPoint2Triangle(pointP, vertex1, vertex2, vertex3)
    implicit none
    ! public
    real(kind=dbPc) :: arebound
    real(kind=dbPc), intent(in) :: alpha, beta
    real(kind=dbPc), intent(in) :: angin, dp2
    ! local
    integer :: n, iii
    real(kind=dbPc) :: x, y
    real(kind=dbPc) :: r1, r2
    real(kind=dbPc) :: pdf
    real(kind=dbPc) :: gama
    real(kind=dbPc) :: xMax, xmin, xmid
    real(kind=dbPc) :: da
    real(kind=dbPc), parameter :: pi = acos(-1.0)
    !
end function closestPoint2Triangle

function ffd(upp, ufp, ddp, nu, rho, rhos)
    implicit none
    ! public
    real(kind=dbPc) :: ffd
    real(kind=dbPc), intent(in) :: upp, ufp
    real(kind=dbPc), intent(in) :: ddp
    real(kind=dbPc), intent(in) :: nu, rho, rhos
    ! local
    real(kind=dbPc) :: cd
    real(kind=dbPc) :: rep, frep
    real(kind=dbPc) :: beta
    real(kind=dbPc) :: mp
    real(kind=dbPc) :: ttp
    real(kind=dbPc), parameter :: pi = acos(-1.0)
    !
    rep = abs(upp - ufp)*ddp/nu
    if (rep == 0.) then
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
        cd = (sqrt(0.5) + sqrt(24.0/rep))**2
        ffd = -pi/8.*cd*rho*ddp**2*abs(upp - ufp)*(upp - ufp)
        !mp = rhos*pi*ddp**3/6.0
        !ttp = rhos*ddp**2/18.0/rho/nu
        !ffd = mp*(ufp-upp)/ttp*frep
    end if
end function ffd

function normal(mmu, sigma)
    implicit none
    ! public
    real(kind=dbPc) :: normal
    real(kind=dbPc), intent(in) :: mmu, sigma
    !local
    integer :: flg
    real(kind=dbPc) :: pi, u1, u2, y1, y2
    save flg
    data flg/0/
    parameter(pi=acos(-1.))
    call random_number(u1)
    call random_number(u2)
    if (flg == 0) then
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
    real(kind=dbPc) :: expdev
    real(kind=dbPc), intent(in) :: lambda
    ! local
    real(kind=dbPc) :: pv
    !
    do while (.true.)
        call random_number(pv)
        if (pv < 1.) then
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
    real(kind=dbPc) :: myerfc
    real(kind=dbPc), intent(in) :: x
    ! function
    real(kind=dbPc) :: gammp, gammq
    ! local
    real(kind=dbPc) :: a
    !
    a = 0.5
    if (x < 0.0) then
        myerfc = 1.0 + gammp(a, x**2)
    else
        myerfc = gammq(a, x**2)
    end if
end function myerfc

function gammq(a, x)
    implicit none
    ! public
    real(kind=dbPc) :: gammq
    real(kind=dbPc), intent(in) :: a
    real(kind=dbPc), intent(in) :: x
    ! local
    integer :: igser, igcf
    real(kind=dbPc) :: gammcf, gamser
    !
    if (x < 0. .or. a <= 0.) return
    if (x < a + 1.0) then
        call gser(gamser, a, x, igser)
        if (igser == 0) then
            gammq = 1.0 - gamser
        else
            gammq = 1.0
        end if
    else
        call gcf(gammcf, a, x, igcf)
        if (igcf == 0) then
            gammq = gammcf
        else
            gammq = 1.0
        end if
    end if
end function gammq

function gammp(a, x)
    implicit none
    ! public
    real(kind=dbPc) :: gammp
    real(kind=dbPc), intent(in) :: a
    real(kind=dbPc), intent(in) :: x
    ! local
    integer :: igser, igcf
    real(kind=dbPc) :: gammcf, gamser
    !
    if (x < 0. .or. a <= 0.) return
    if (x < a + 1.0) then
        call gser(gamser, a, x, igser)
        if (igser == 0) then
            gammp = gamser
        else
            gammp = 1.0
        end if
    else
        call gcf(gammcf, a, x, igcf)
        if (igcf == 0) then
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
    real(kind=dbPc) :: gamser
    real(kind=dbPc), intent(in) :: a
    real(kind=dbPc), intent(in) :: x
    ! function
    real(kind=dbPc) :: gammln
    ! local
    integer, parameter :: itmax = 100
    integer :: n
    real(kind=dbPc), parameter :: eps = 3.0e-7
    real(kind=dbPc) :: gln
    real(kind=dbPc) :: ap
    real(kind=dbPc) :: del
    real(kind=dbPc) :: asum
    !
    gln = gammln(a)
    if (x <= 0.) then
        if (x < 0.) return
        gamser = 0.
        return
    end if
    ap = a
    asum = 1.0/a
    del = asum
    do n = 1, itmax
        ap = ap + 1.0
        del = del*x/ap
        asum = asum + del
        if (abs(del) < abs(asum)*eps) goto 4000
    end do
    print *, 'iter num > itmax in gser'
    igser = 1
    return
4000 continue
    gamser = asum*exp(-x + a*log(x) - gln)
    igser = 0
end subroutine gser

subroutine gcf(gammcf, a, x, igcf)
    implicit none
    ! public
    integer :: igcf
    real(kind=dbPc) :: gammcf
    real(kind=dbPc), intent(in) :: a
    real(kind=dbPc), intent(in) :: x
    ! function
    real(kind=dbPc) :: gammln
    ! local
    integer, parameter :: itmax = 100
    integer :: i
    real(kind=dbPc) :: gln
    real(kind=dbPc), parameter :: eps = 3.0e-7
    real(kind=dbPc), parameter :: fpmin = 1.0e-30
    real(kind=dbPc) :: an, b, c, d
    real(kind=dbPc) :: del
    real(kind=dbPc) :: h
    !
    gln = gammln(a)
    b = x + 1.0 - a
    c = 1.0/fpmin
    d = 1.0/b
    h = d
    do i = 1, itmax
        an = -i*(i - a)
        b = b + 2.0
        d = an*d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0/d
        del = d*c
        h = h*del
        if (abs(del - 1.0) < eps) goto 4001
    end do
    print *, 'iter num > itmax in gcf'
    igcf = 1
    return
4001 continue
    gammcf = exp(-x + a*log(x) - gln)*h
    igcf = 0
end subroutine gcf

function gammln(xx)
    implicit none
    ! public
    real(kind=dbPc) :: gammln
    real(kind=dbPc), intent(in) :: xx
    ! local
    integer :: j
    real(kind=dbPc) :: ser
    real(kind=dbPc) :: tmp
    real(kind=dbPc) :: x, y
    real(kind=dbPc), save :: cof(6)
    real(kind=dbPc), save :: stp
    data cof, stp/76.18009172947146d0, -86.50532032941677d0, &
        24.01409824083091d0, -1.231739572450155d0, &
        .1208650973866179d-2, -.5395239384953d-5, &
        2.5066282746310005d0/
    !
    x = xx
    y = x
    tmp = x + 5.5d0
    tmp = (x + 0.5d0)*log(tmp) - tmp
    ser = 1.000000000190015d0
    do j = 1, 6
        y = y + 1.d0
        ser = ser + cof(j)/y
    end do
    gammln = tmp + log(stp*ser/x)
end function gammln

function arebound(alpha, beta, angin, dp2)
    implicit none
    ! public
    real(kind=dbPc) :: arebound
    real(kind=dbPc), intent(in) :: alpha, beta
    real(kind=dbPc), intent(in) :: angin, dp2
    ! local
    integer :: n, iii
    real(kind=dbPc) :: x, y
    real(kind=dbPc) :: r1, r2
    real(kind=dbPc) :: pdf
    real(kind=dbPc) :: gama
    real(kind=dbPc) :: xMax, xmin, xmid
    real(kind=dbPc) :: da
    real(kind=dbPc), parameter :: pi = acos(-1.0)
    !
    gama = 4.0/9.0*beta**2/(alpha + beta)**2/dp2
    xmin = -angin
    xMax = sqrt(angin/gama)*2.0 - angin
    if (xMax > pi) then
        xMax = pi
    end if
    da = (xMax - xmin)/500.0
    y = 1.0
    pdf = 0.0
    n = 0
    do
        call random_number(r1)
        call random_number(r2)
        n = n + 1
        if (n > 10000) then
            print *, 'arebound, n>10000', dp2
            x = 60.0/180.0*pi
            exit
        end if
        iii = int(r1*(xMax - xmin)/da)
        x = dfloat(iii)*da + 0.5*da + xmin
        y = r2
        pdf = gama*(x + angin)/angin*log(2.0/gama*angin/(x + angin)**2)
        pdf = pdf*da
        if (y <= pdf) exit
    end do
    arebound = x + angin
end function arebound
