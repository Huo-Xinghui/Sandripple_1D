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
