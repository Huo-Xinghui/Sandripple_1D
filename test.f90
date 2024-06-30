program gridGeneration
    implicit none
    integer, parameter :: nz = 20
    real, parameter :: zMax = 10.0, zUni = 3.0, zDiffMin = 0.1
    real :: zDiffMax, zGridCtrl, zCurrent, zDiffCurrent
    real, dimension(nz) :: vectorGrid, zDiff
    integer :: i

    ! Calculate the maximum grid spacing for below and above zUni
    zDiffMax = zMax/real(nz - 1)

    ! Define the refinement ratio for below zUni
    zGridCtrl = 1.5 ! This should be chosen such that the grid requirements are met

    ! Initialize grid generation
    zCurrent = 0.0
    vectorGrid(1) = zCurrent
    zDiff(1) = zDiffMin ! Initial grid spacing, could be refined later

    ! Generate vectorGrid array (below zUni part, exponential refinement)
    do i = 2, nz
        if (zCurrent < zUni) then
            ! Exponential refinement towards zUni, but not exceeding zDiffMax
            zDiffCurrent = zDiffMin*zGridCtrl**(i - 2)
            zDiffCurrent = min(zDiffCurrent, zDiffMax)
        else
            ! Uniform spacing above zUni
            zDiffCurrent = zDiffMax
        end if

        zCurrent = zCurrent + zDiffCurrent
        vectorGrid(i) = zCurrent
        zDiff(i) = zDiffCurrent ! Store the current grid spacing in zDiff array

        ! If reached or exceeded zMax, adjust the last grid point and exit
        if (zCurrent >= zMax) then
            vectorGrid(nz) = zMax
            exit
        end if
    end do

    ! If the loop exited before reaching nz, adjust the grid spacing for the last point
    if (i < nz) then
        zDiff(nz) = vectorGrid(nz) - vectorGrid(nz - 1)
    end if

    ! Ensure the first and last grid points are as specified
    vectorGrid(1) = 0.0
    vectorGrid(nz) = zMax

    ! Output vectorGrid and zDiff arrays
    do i = 1, nz
        print *, 'vectorGrid(', i, ') = ', vectorGrid(i), ' zDiff(', i, ') = ', zDiff(i)
    end do
end program gridGeneration

program particleSizeDistribution
    use someModule ! 假设npdf, dpa, dpStddDev和dbPc在此模块中定义
    implicit none
    integer :: nParticles, iParticle
    real(kind=dbPc), dimension(npdf) :: histogram
    real(kind=dbPc) :: randomNum, cumulativeProb, particleSize

    ! 生成正态分布直方图
    call normalDistHistogram(histogram)

    ! 归一化直方图
    histogram = histogram/sum(histogram)

    ! 设定颗粒数量
    nParticles = 1000

    ! 分配粒径给每个颗粒
    do iParticle = 1, nParticles
        call random_number(randomNum) ! 生成0到1之间的随机数
        cumulativeProb = 0.0
        do i = 1, npdf
            cumulativeProb = cumulativeProb + histogram(i)
            if (cumulativeProb >= randomNum) then
                particleSize = binStart + (i - 0.5)*binWidth
                ! 使用particleSize做一些事情，比如打印出来
                print *, 'Particle', iParticle, 'Size:', particleSize
                exit
            end if
        end do
    end do
end program particleSizeDistribution

program particle_linked_list
    implicit none
    integer, parameter :: dbPc = selected_real_kind(15, 307) ! 双精度

    type particle
        real(kind=dbPc) :: position(3)
        real(kind=dbPc) :: velocity(3)
        real(kind=dbPc) :: diameter
    end type particle

    type particleLink
        type(particle) :: data
        type(particleLink), pointer :: next => null()
    end type particleLink

    type(particleLink), pointer :: head => null(), tail => null(), current => null()
    integer :: pNum, i

    ! 假设pNum由用户输入或迭代过程动态确定
    print *, "请输入颗粒数量："
    read *, pNum

    ! 创建链表
    do i = 1, pNum
        allocate (current)
        current%data = particle([real(i, kind=dbPc), real(i + 1, kind=dbPc), real(i + 2, kind=dbPc)], &
                                [0.1*real(i, kind=dbPc), 0.1*real(i + 1, kind=dbPc), 0.1*real(i + 2, kind=dbPc)], &
                                0.01*real(i, kind=dbPc))
        if (.not. associated(head)) then
            head => current
            tail => current
        else
            tail%next => current
            tail => current
        end if
    end do

    ! 遍历链表并打印颗粒信息
    current => head
    while(associated(current))
    print *, "颗粒位置：", current%data%position
    print *, "颗粒速度：", current%data%velocity
    print *, "颗粒直径：", current%data%diameter
    current => current%next
    end while

    ! 清理链表
    current => head
    while(associated(current))
    tail => current
    current => current%next
    deallocate (tail)
    end while

end program particle_linked_list

! 假设你要在链表末尾添加一个新颗粒
allocate (current)
current%data = particle([新位置x, 新位置y, 新位置z], [新速度x, 新速度y, 新速度z], 新直径)
if (.not. associated(head)) then
    head => current
    tail => current
else
    tail%next => current
    tail => current
end if

! 假设你要删除具有特定属性的颗粒，例如位置为(x0, y0, z0)的颗粒
current => head
previous => null()
while(associated(current))
if (all(current%data%position == [x0, y0, z0])) then
    if (associated(previous)) then
        previous%next => current%next
    else
        head => current%next
    end if
    if (associated(current, tail)) then
        tail => previous
    end if
    deallocate (current)
    ! 跳出循环或设置标志以避免重复删除
    exit
else
    previous => current
    current => current%next
end if
end while
