module fluid
  use setup
  use public_val
  use vectors
  use immersed_boundary
  implicit none
  save
  ! out
  real(dpc), dimension(nxdim, my, mz), target :: u, v, w, p, nut
  real(dpc), dimension(nxdim, my, mz), target :: uf, vf, wf
  real(dpc), dimension(nxdim, my, mz), target :: dpdx, dpdy, dpdz
  real(dpc), dimension(nxdim, my, mz, 3), target :: ctval
  real(dpc), dimension(nxdim, my, mz), target :: utao
  real(dpc), dimension(xpdim, ypnum), target :: utaosf
  ! in
  real(dpc), dimension(nxdim, my, mz) :: fptx, fpty, fptz
  real(dpc), dimension(nxdim, my, mz) :: phirho
  ! private
  real(dpc), private, dimension(nxdim, my, mz) :: ufx, vfx, wfx
  real(dpc), private, dimension(nxdim, my, mz) :: uc, vc, wc
  real(dpc), private, dimension(nxdim, my, mz) :: fmmx, flmx
  real(dpc), private, dimension(nxdim, my, mz) :: cs2
  real(dpc), private, dimension(nxdim, my, mz, 3) :: flox, difx, sorx
  real(dpc), private, dimension(nxdim, my, mz, 3) :: floxx, difxx, sorxx
  real(dpc), private, dimension(nxdim, my, mz, 3) :: rhsx
  real(dpc), private, dimension(nxdim, my, mz) :: ux, vx, wx
  !
  !
  real(dpc) :: ttat, atat
  real(dpc), dimension(nxdim, my, mz) :: rhsux, rhsuy
  real(dpc), dimension(nxdim, my, mz) :: rhsvx, rhsvy
  real(dpc), dimension(nxdim, my, mz) :: rhswx, rhswy
  real(dpc), dimension(3) :: avepg
  real(dpc), dimension(3) :: tutao
  real(dpc), dimension(3) :: aveu
  real(dpc), dimension(mz) :: ttau, atau
  real(dpc), dimension(mz) :: ttauu, atauu
  real(dpc), dimension(mz) :: ttavv, atavv
  real(dpc), dimension(mz) :: ttaww, ataww
  real(dpc), dimension(mz) :: ttauw, atauw
  real(dpc), dimension(mz) :: ttac, atac
contains
  subroutine fluid_init
    ! initialize fluid field
    use fgrid, only : z
    use grid_exchange, only : ubdy
    implicit none
    ! local
    integer :: i, j, k
    real(dpc) :: rr
    !
    call random_seed()
    ! u v w
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    if ((z(k)-zb)>z0) then
      call random_number(rr)
      u(i, j, k) = wind/kapa*log((z(k)-zb)/z0)*(1.0d0+(-2.0d0*rr+1.0d0)*tld)
      v(i, j, k) = cZero + (-2.0d0*rr+1.0d0)*tld
      w(i, j, k) = cZero + (-2.0d0*rr+1.0d0)*tld
      ufx(i, j, k) = wind/kapa*log((z(k)-zb)/z0)
      vfx(i, j, k) = cZero
      wfx(i, j, k) = cZero
    else
      u(i, j, k) = cZero
      v(i, j, k) = cZero
      w(i, j, k) = cZero
      ufx(i, j, k) = cZero
      vfx(i, j, k) = cZero
      wfx(i, j, k) = cZero
    end if
    end do
    end do
    end do
    call ubdy(u, v, w)
    call ubdy(ufx, vfx, wfx)
    ! p
    p = cZero
    dpdx = cZero
    dpdy = cZero
    dpdz = cZero
    ! LES
    nut = cZero
    flmx = cZero
    fmmx = cZero
    ! particle
    phirho = cOne
    ! others
    atac = cZero
    fptx = cZero
    fpty = cZero
    fptz = cZero
    dpaa = dpa
    last = 1
    time = cZero
    ttat = cZero
    ttau = cZero
    ttauu = cZero
    ttavv = cZero
    ttaww = cZero
    ttauw = cZero
    ttac = cZero
    nnp = 0
    uflx = cZero
    wflx = cZero
    uflxz = cZero
    wflxz = cZero
  end subroutine fluid_init

  subroutine dpdn
    ! calculate dpdx dpdy dpdz
    use grid_exchange, only : pbdy
    implicit none
    ! local
    integer :: i, j, k
    !
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    dpdx(i, j, k) = (p(i, j, k)-p(i-1, j, k))/xdif
    dpdy(i, j, k) = (p(i, j, k)-p(i, j-1, k))/ydif
    dpdz(i, j, k) = (p(i, j, k)-p(i, j, k-1))/zdif
    end do
    end do
    end do
    call pbdy(p, dpdx, dpdy, dpdz)
  end subroutine dpdn

  subroutine p_interpolation
    ! fix dpdx dpdy dpdz near IM-boundary
    implicit none
    ! local
    integer :: n
    real(dpc), dimension(:, :, :), pointer :: valp
    type (vector), dimension(:), pointer :: norp
    real(dpc), dimension(:), pointer :: ipdp
    !
    do n = 1, nnif
    valp => dpdx
    norp => norpu
    ipdp => ipdisu
    call interpolation(1, n, norp, ipdp, valp)
    valp => dpdy
    norp => norpv
    ipdp => ipdisv
    call interpolation(2, n, norp, ipdp, valp)
    valp => dpdz
    norp => norpw
    ipdp => ipdisw
    call interpolation(3, n, norp, ipdp, valp)
    end do
  end subroutine p_interpolation

  subroutine centu
    ! interpolate the velocity to the center of control volume
    use grid_exchange, only : ucbdy
    implicit none
    ! local
    integer i, j, k
    !
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    uc(i, j, k) = 0.5d0*(u(i, j, k) + u(i+1, j, k))
    vc(i, j, k) = 0.5d0*(v(i, j, k) + v(i, j+1, k))
    wc(i, j, k) = 0.5d0*(w(i, j, k) + w(i, j, k+1))
    end do
    end do
    end do
    call ucbdy(uc, vc, wc)
    ctval(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, 1) = uc
    ctval(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, 2) = vc
    ctval(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, 3) = wc
  end subroutine centu

  subroutine nut_cal
    ! calculate nut
    ! Charles Meneveau et al., 1996
    use grid_exchange, only : sbdy, lmbdy
    use fgrid, only : x, y, z
    implicit none
    ! local
    integer :: i, j, k, ii, jj, kk
    real(dpc) :: d2
    real(dpc) :: fkm, fkl
    real(dpc) :: tn, eps
    real(dpc) :: mm, lm
    real(dpc), dimension(3) :: pLoc
    real(dpc), dimension(nxdim, my, mz) :: s11, s12, s13, s22, s23, s33
    real(dpc), dimension(nxdim, my, mz) :: ss11, ss12, ss13, ss22, ss23, ss33
    real(dpc), dimension(nxdim, my, mz) :: ts11, ts12, ts13, ts22, ts23, ts33
    real(dpc), dimension(nxdim, my, mz) :: tss11, tss12, tss13, tss22, tss23, tss33
    real(dpc), dimension(nxdim, my, mz) :: smd, tsmd
    real(dpc), dimension(nxdim, my, mz) :: tuc, tvc, twc
    !
    ! sij=Sij, smd=|S|
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    s11(i, j, k) = (u(i+1, j, k)-u(i, j, k))/xdif
    s12(i, j, k) = ((uc(i, j+1, k)-uc(i, j-1, k))/ydif+(vc(i+1, j, k)-vc(i-1, j, k))/xdif)/4.0d0
    s13(i, j, k) = ((uc(i, j, k+1)-uc(i, j, k-1))/zdif+(wc(i+1, j, k)-wc(i-1, j, k))/xdif)/4.0d0
    s22(i, j, k) = (v(i, j+1, k)-v(i, j, k))/ydif
    s23(i, j, k) = ((vc(i, j, k+1)-vc(i, j, k-1))/zdif+(wc(i, j+1, k)-wc(i, j-1, k))/ydif)/4.0d0
    s33(i, j, k) = (w(i, j, k+1)-w(i, j, k))/zdif
    smd(i, j, k) = sqrt(2.0d0*(s11(i, j, k)**2+2.0d0*(s12(i, j, k)**2)+2.0d0*(s13(i, j, k)**2)+ &
      s22(i, j, k)**2+2.0d0*(s23(i, j, k)**2)+s33(i, j, k)**2))
    end do
    end do
    end do
    call sbdy(uc, vc, wc, u, v, w, s11, s12, s13, s22, s23, s33, smd)
    ss11 = smd*s11
    ss12 = smd*s12
    ss13 = smd*s13
    ss22 = smd*s22
    ss23 = smd*s23
    ss33 = smd*s33
    ! tX=hat(X)
    call tstf(s11, ts11)
    call tstf(s12, ts12)
    call tstf(s13, ts13)
    call tstf(s22, ts22)
    call tstf(s23, ts23)
    call tstf(s33, ts33)
    call tstf(ss11, tss11)
    call tstf(ss12, tss12)
    call tstf(ss13, tss13)
    call tstf(ss22, tss22)
    call tstf(ss23, tss23)
    call tstf(ss33, tss33)
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    tsmd(i, j, k) = sqrt(2.0d0*(ts11(i, j, k)**2+2.0d0*(ts12(i, j, k)**2)+2.0d0*(ts13(i, j, k)**2)+ &
      ts22(i, j, k)**2+2.0d0*(ts23(i, j, k)**2)+ts33(i, j, k)**2))
    end do
    end do
    end do
    d2 = (xdif*ydif*zdif)**(2.0d0/3.0d0)
    ! ssij=Mij=2*Delta^2*(hat(|S|Sij)-4|hat(S)|hat(Sij))
    ss11 = 2.0d0*d2*(tss11-4.0d0*tsmd*ts11)
    ss12 = 2.0d0*d2*(tss12-4.0d0*tsmd*ts12)
    ss13 = 2.0d0*d2*(tss13-4.0d0*tsmd*ts13)
    ss22 = 2.0d0*d2*(tss22-4.0d0*tsmd*ts22)
    ss23 = 2.0d0*d2*(tss23-4.0d0*tsmd*ts23)
    ss33 = 2.0d0*d2*(tss33-4.0d0*tsmd*ts33)
    ! s11=uu, s12=uv, s13=uw, s22=vv, s23=vw, s33=ww 
    s11 = uc*uc
    s12 = uc*vc
    s13 = uc*wc
    s22 = vc*vc
    s23 = vc*wc
    s33 = wc*wc
    call tstf(uc, tuc)
    call tstf(vc, tvc)
    call tstf(wc, twc)
    call tstf(s11, ts11)
    call tstf(s12, ts12)
    call tstf(s13, ts13)
    call tstf(s22, ts22)
    call tstf(s23, ts23)
    call tstf(s33, ts33)
    ! tssij=Lij=hat(UiUj)-hat(Ui)hat(Uj)
    tss11 = ts11 - tuc*tuc
    tss12 = ts12 - tuc*tvc
    tss13 = ts13 - tuc*twc
    tss22 = ts22 - tvc*tvc
    tss23 = ts23 - tvc*twc
    tss33 = ts33 - twc*twc
    ! s11=f_MM, s12=f_LM
    if (last==1) then
      s11 = ss11*ss11 + 2.0d0*ss12*ss12 + 2.0d0*ss13*ss13 + ss22*ss22 + 2.0d0*ss23*ss23 + ss33*ss33
      s12 = 0.0256d0*s11
    else
      do k = sz-1, ez+1
      do j = sy-1, ey+1
      do i = sx-1, ex+1
      pLoc(1) = x(i) - uc(i, j, k)*dt
      pLoc(2) = y(j) - vc(i, j, k)*dt
      pLoc(3) = z(k) - wc(i, j, k)*dt
      if (pLoc(1)<x(sx-1)) pLoc(1) = x(sx-1)
      if (pLoc(1)>x(ex+1)) pLoc(1) = x(ex+1)
      if (pLoc(2)<y(sy-1)) pLoc(2) = y(sy-1)
      if (pLoc(2)>y(ey+1)) pLoc(2) = y(ey+1)
      if (pLoc(3)<z(sz-1)) pLoc(3) = z(sz-1)
      if (pLoc(3)>z(ez+1)) pLoc(3) = z(ez+1)
      ii = int((pLoc(1)-x(sx-1))/xdif) + (sx-1)
      jj = int((pLoc(2)-y(sy-1))/ydif) + (sy-1)
      kk = int((pLoc(3)-z(sz-1))/zdif) + (sz-1)
      if (ii>ex) ii = ex
      if (jj>ey) jj = ey
      if (kk>ez) kk = ez
      call val_in_cell(x(ii), x(ii+1), y(jj), y(jj+1), z(kk), z(kk+1), &
        fmmx(ii, jj, kk), fmmx(ii, jj+1, kk), fmmx(ii+1, jj, kk), fmmx(ii+1, jj+1, kk), &
        fmmx(ii, jj, kk+1), fmmx(ii, jj+1, kk+1), fmmx(ii+1, jj, kk+1), fmmx(ii+1, jj+1, kk+1), pLoc, fkm)
      call val_in_cell(x(ii), x(ii+1), y(jj), y(jj+1), z(kk), z(kk+1), &
        flmx(ii, jj, kk), flmx(ii, jj+1, kk), flmx(ii+1, jj, kk), flmx(ii+1, jj+1, kk), &
        flmx(ii, jj, kk+1), flmx(ii, jj+1, kk+1), flmx(ii+1, jj, kk+1), flmx(ii+1, jj+1, kk+1), pLoc, fkl)
      tn = 1.5*sqrt(d2)*((flmx(i, j, k)*fmmx(i, j, k))**(-1.0d0/8.0d0))
      if (abs(tn)>1.0d-8) then
        eps = (dt/tn/(cOne+dt/tn)
      else
        eps = cZero
      end if
      mm = ss11(i, j, k)**2 + 2.0d0*(ss12(i, j, k)**2) + 2.0d0*(ss13(i, j, k)**2) + &
        ss22(i, j, k)**2 + 2.0d0*(ss23(i, j, k)**2) + ss33(i, j, k)**2
      lm = tss11(i, j, k)*ss11(i, j, k) + 2.0d0*tss12(i, j, k)*ss12(i, j, k) + &
        2.0d0*tss13(i, j, k)*ss13(i, j, k) + tss22(i, j, k)*ss22(i, j, k) + &
        2.0d0*tss23(i, j, k)*ss23(i, j, k) + tss33(i, j, k)*ss33(i, j, k)
      s11(i, j, k) = eps*mm + (cOne-eps)*fkm 
      s12(i, j, k) = eps*lm + (cOne-eps)*fkl 
      if (s12(i, j, k)<cZero) s12(i, j, k) = cZero
      end do
      end do
      end do
    end if
    fmmx = s11
    flmx = s12
    call lmbdy(fmmx, flmx)
    do k = sz-1, ez+1
    do j = sy-1, ey+1
    do i = sx-1, ex+1
    if (fmmx(i, j, k)==cZero) then
      cs2(i, j, k) = cZero
      nut(i, j, k) = cZero
    else
      cs2(i, j, k) = flmx(i, j, k)/fmmx(i, j, k)
      nut(i, j, k) = cs2*d2*smd(i, j, k)
    end if
    end do
    end do
    end do
  end subroutine nut_cal

  subroutine tstf(s, ts)
    ! test filter
    implicit none
    include "mpif.h"
    ! public
    real(dpc), intent(in), dimension(nxdim, my, mz) :: s
    real(dpc), dimension(nxdim, my, mz) :: ts
    ! local
    integer i, j, k
    real(dpc), dimension(nxdim, my, mz) :: tts
    !
    do k = sz-1, ez+1
    do j = sy-1, ey+1
    do i = sx, ex
    ts(i, j, k) = 0.25d0*s(i-1, j, k) + 0.5d0*s(i, j, k) + 0.25d0*s(i+1, j, k)
    end do
    end do
    end do
    ts(sx-1, sy-1:ey+1, sz-1:ez+1) = 0.5d0*s(sx-1, sy-1:ey+1, sz-1:ez+1) + 0.5d0*s(sx, sy-1:ey+1, sz-1:ez+1)
    ts(ex+1, sy-1:ey+1, sz-1:ez+1) = 0.5d0*s(ex, sy-1:ey+1, sz-1:ez+1) + 0.5d0*s(ex+1, sy-1:ey+1, sz-1:ez+1)
    do k = sz-1, ez+1
    do j = sy, ey
    do i = sx-1, ex+1
    tts(i, j, k) = 0.25d0*ts(i, j-1, k) + 0.5d0*ts(i, j, k) + 0.25d0*ts(i, j+1, k)
    end do
    end do
    end do
    tts(sx-1:ex+1, sy-1, sz-1:ez+1) = 0.5d0*ts(sx-1:ex+1, sy-1, sz-1:ez+1) + 0.5d0*ts(sx-1:ex+1, sy, sz-1:ez+1)
    tts(sx-1:ex+1, ey+1, sz-1:ez+1) = 0.5d0*ts(sx-1:ex+1, ey, sz-1:ez+1) + 0.5d0*ts(sx-1:ex+1, ey+1, sz-1:ez+1)
    do k = sz, ez
    do j = sy-1, ey+1
    do i = sx-1, ex+1
    ts(i, j, k) = 0.25d0*tts(i, j, k-1) + 0.5d0*tts(i, j, k) + 0.25d0*tts(i, j, k+1)
    end do
    end do
    end do
    ts(sx-1:ex+1, sy-1:ey+1, sz-1) = 0.5d0*s(sx-1:ex+1, sy-1:ey+1, sz-1) + 0.5d0*s(sx-1:ex+1, sy-1:ey+1, sz)
    ts(sx-1:ex+1, sy-1:ey+1, ez+1) = 0.5d0*s(sx-1:ex+1, sy-1:ey+1, ez) + 0.5d0*s(sx-1:ex+1, sy-1:ey+1, ez+1)
  end subroutine tstf

  subroutine val_in_cell(x1, x2, y1, y2, z1, z2, u1, u2, u3, u4, u5, u6, u7, u8, loc_p, val_p)
    ! value in cell at location loc_p
    ! x1=x(i) x2=x(i+1)
    ! y1=y(j) y2=y(j+1)
    ! z1=z(k) z2=z(k+1)
    ! u1=u(i, j, k), u2=u(i, j+1, k), u3=u(i+1, j, k), u4=u(i+1, j+1, k)
    ! u5=u(i, j, k+1), u6=u(i, j+1, k+1), u7=u(i+1, j, k+1), u8=u(i+1, j+1, k+1)
    use array_cal
    implicit none
    ! public
    real(dpc), intent(in) :: x1, x2, y1, y2, z1, z2
    real(dpc), intent(in) :: u1, u2, u3, u4, u5, u6, u7, u8
    real(dpc), intent(in), dimension(3) :: loc_p
    real(dpc) :: val_p
    !local
    real(dpc) :: a, b, c, d, e, f
    real(dpc) :: xu1, xu2, xu3, xu4
    !
    a = abs((loc_p(1)-x1)/(x2-x1))
    if (a>cOne) a = cOne
    b = cOne - a
    xu1 = u1*b + u3*a
    xu2 = u2*b + u4*a
    xu3 = u5*b + u7*a
    xu4 = u6*b + u8*a
    c = abs((loc_p(2)-y1)/(y2-y1))
    if (c>cOne) c = cOne
    d = cOne - c
    yu1 = xu1*d + xu2*c
    yu2 = xu3*d + xu4*c
    e = abs((loc_p(3)-z1)/(z2-z1))
    if (e>cOne) e = cOne
    f = cOne - e
    uu = yu1*f + yu2*e
  end subroutine val_in_cell

  subroutine ux_cal
    ! calculate auxiliary velocity field
    use grid_exchange, only : ubdy
    implicit none
    ! local
    integer :: i, j, k, nnk
    real(dpc), dimension(nxdim, my, mz) :: flo, dif, sor
    !
    level1: do nnk = 1, 3
    call conv(nnk, flo)
    call diff(nnk, dif)
    call sour(nnk, sor)
    if (last==1) then
      floxx(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, nnk) = flo
      difxx(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, nnk) = dif
      sorxx(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, nnk) = sor
    else
      floxx = flox
      difxx = difx
      sorxx = sorx
    end if
    flox(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, nnk) = flo
    difx(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, nnk) = dif
    sorx(sx-1:ex+1, sy-1:ey+1, sz-1:ez+1, nnk) = sor
    !
    level2: do k = sz, ez
    level3: do j = sy, ey
    level4: do i = sx, ex
    if (nnk==1) then
      rhsx(i, j, k, nnk) = -flox(i, j, k, nnk) + difx(i, j, k, nnk) + sorx(i, j, k,nnk) - &
        dpdx(i, j, k)/(rho*phirho(i, j, k))
      ux(i, j, k) = u(i, j, k) + dt*rhsx(i, j, k, nnk)
    else if (nnk==2) then
      rhsx(i, j, k, nnk) = -flox(i, j, k, nnk) + difx(i, j, k, nnk) + sorx(i, j, k,nnk) - &
        dpdy(i, j, k)/(rho*phirho(i, j, k))
      vx(i, j, k) = v(i, j, k) + dt*rhsx(i, j, k, nnk)
    else if (nnk==3) then
      rhsx(i, j, k, nnk) = -flox(i, j, k, nnk) + difx(i, j, k, nnk) + sorx(i, j, k,nnk) - &
        dpdz(i, j, k)/(rho*phirho(i, j, k))
      wx(i, j, k) = w(i, j, k) + dt*rhsx(i, j, k, nnk)
    end if
    end do level4
    end do level3
    end do level2
    end do level1
    call ubdy(ux, vx, wx)
    do k = sz-1, ez+1
    do j = sy-1, ey+1
    do i = sx-1, ex+1
    if (fu(i, j, k)==1) then
      uf(i, j, k) = cZero
    else
      uf(i, j, k) = ux(i, j, k)
    end if
    if (fv(i, j, k)==1) then
      vf(i, j, k) = cZero
    else
      vf(i, j, k) = vx(i, j, k)
    end if
    if (fw(i, j, k)==1) then
      wf(i, j, k) = cZero
    else
      wf(i, j, k) = wx(i, j, k)
    end if
    end do
    end do
    end do
  end subroutine ux_cal

  subroutine conv(tag, floRes)
    !calculate convection terms
    implicit none
    ! public
    integer, intent(in) :: tag
    real(dpc), dimension(nxdim, my, mz) :: floRes
    ! local
    integer i, j, k
    !
    if (tag==1) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      floRes(i, j, k) = (uc(i, j, k)**2-uc(i-1, j, k)**2)/xdif + &
        ((u(i, j+1, k)+u(i, j, k))/2.0d0*(v(i, j+1, k)+v(i-1, j+1, k))/2.0d0- &
        (u(i, j, k)+u(i, j-1, k))/2.0d0*(v(i, j, k)+v(i-1, j, k))/2.0d0)/ydif + &
        ((u(i, j, k+1)+u(i, j, k))/2.0d0*(w(i, j, k+1)+w(i-1, j, k+1))/2.0d0- &
        (u(i, j, k)+u(i, j, k-1))/2.0d0*(w(i, j, k)+w(i-1, j, k))/2.0d0)/zdif
      end do
      end do
      end do
    else if (tag==2) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      floRes(i, j, k) = (vc(i, j, k)**2-vc(i, j-1, k)**2)/ydif + &
        ((v(i+1, j, k)+v(i, j, k))/2.0d0*(u(i+1, j, k)+u(i+1, j-1, k))/2.0d0- &
        (v(i, j, k)+v(i-1, j, k))/2.0d0*(u(i, j, k)+u(i, j-1, k))/2.0d0)/xdif + &
        ((v(i, j, k+1)+v(i, j, k))/2.0d0*(w(i, j, k+1)+w(i, j-1, k+1))/2.0d0- &
        (v(i, j, k)+v(i, j, k-1))/2.0d0*(w(i, j, k)+w(i, j-1, k))/2.0d0)/zdif
      end do
      end do
      end do
    else if (tag==3) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      floRes(i, j, k) = (wc(i, j, k)**2-wc(i, j, k-1)**2)/zdif + &
        ((w(i, j+1, k)+w(i, j, k))/2.0d0*(v(i, j+1, k)+v(i, j+1, k-1))/2.0d0- &
        (w(i, j, k)+w(i, j-1, k))/2.0d0*(v(i, j, k)+v(i, j, k-1))/2.0d0)/ydif + &
        ((w(i+1, j, k)+w(i, j, k))/2.0d0*(u(i+1, j, k)+u(i+1, j, k-1))/2.0d0- &
        (w(i, j, k)+w(i-1, j, k))/2.0d0*(u(i, j, k)+u(i, j, k-1))/2.0d0)/xdif
      end do
      end do
      end do
    else
      print*, 'tag error in conv!'
      stop
    end if
  end subroutine conv

  subroutine diff(tag, difRes)
    ! calculate diffusion terms
    implicit none
    ! public
    integer :: tag
    ! local
    integer i, j, k
    real(dpc) :: nute, nutw, nutn, nuts, nutt, nutb
    !
    if (tag==1) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      nute = nut(i, j, k)
      nutw = nut(i-1, j, k)
      nutn = (nut(i, j, k)+nut(i, j+1, k)+nut(i-1, j+1, k)+nut(i-1, j, k))/4.0d0
      nuts = (nut(i, j, k)+nut(i, j-1, k)+nut(i-1, j-1, k)+nut(i-1, j, k))/4.0d0
      nutt = (nut(i, j, k)+nut(i, j, k+1)+nut(i-1, j, k+1)+nut(i-1, j, k))/4.0d0
      nutb = (nut(i, j, k)+nut(i, j, k-1)+nut(i-1, j, k-1)+nut(i-1, j, k))/4.0d0
      difRes(i, j, k) = ((nu+nute)*(u(i+1, j, k)-u(i, j, k))/xdif- &
        (nu+nutw)*(u(i, j, k)-u(i-1, j, k))/xdif)/xdif + &
        ((nu+nutn)*(u(i, j+1, k)-u(i, j, k))/ydif- &
        (nu+nuts)*(u(i, j, k)-u(i, j-1, k))/ydif)/ydif + &
        ((nu+nutt)*(u(i, j, k+1)-u(i, j, k))/zdif- &
        (nu+nutb)*(u(i, j, k)-u(i, j, k-1))/zdif)/zdif
      end do
      end do
      end do
    else if (tag==2) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      nute = (nut(i, j, k)+nut(i+1, j, k)+nut(i+1, j-1, k)+nut(i ,j-1, k))/4.0d0
      nutw = (nut(i, j, k)+nut(i-1, j, k)+nut(i-1, j-1, k)+nut(i ,j-1, k))/4.0d0
      nutn = nut(i, j, k)
      nuts = nut(i, j-1, k)
      nutt = (nut(i, j, k)+nut(i, j, k+1)+nut(i, j-1, k+1)+nut(i, j-1, k))/4.0d0
      nutb = (nut(i, j, k)+nut(i, j, k-1)+nut(i, j-1, k-1)+nut(i, j-1, k))/4.0d0
      difRes(i, j, k) = ((nu+nute)*(v(i+1, j, k)-v(i, j, k))/xdif- &
        (nu+nutw)*(v(i, j, k)-v(i-1, j, k))/xdif)/xdif + &
        ((nu+nutn)*(v(i, j+1, k)-v(i, j, k))/ydif- &
        (nu+nuts)*(v(i, j, k)-v(i, j-1, k))/ydif)/ydif + &
        ((nu+nutt)*(v(i, j, k+1)-v(i, j, k))/zdif- &
        (nu+nutb)*(v(i, j, k)-v(i, j, k-1))/zdif)/zdif
      end do
      end do
      end do
    else if (tag==3) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      nute = (nut(i, j, k)+nut(i+1, j, k)+nut(i+1, j, k-1)+nut(i ,j, k-1))/4.0d0
      nutw = (nut(i, j, k)+nut(i-1, j, k)+nut(i-1, j, k-1)+nut(i ,j, k-1))/4.0d0
      nutn = (nut(i, j, k)+nut(i, j+1, k)+nut(i, j+1, k-1)+nut(i, j, k-1))/4.0d0
      nuts = (nut(i, j, k)+nut(i, j-1, k)+nut(i, j-1, k-1)+nut(i, j, k-1))/4.0d0
      nutt = nut(i, j, k)
      nutb = nut(i, j, k-1)
      difRes(i, j, k) = ((nu+nute)*(w(i+1, j, k)-w(i, j, k))/xdif- &
        (nu+nutw)*(w(i, j, k)-w(i-1, j, k))/xdif)/xdif + &
        ((nu+nutn)*(w(i, j+1, k)-w(i, j, k))/ydif- &
        (nu+nuts)*(w(i, j, k)-w(i, j-1, k))/ydif)/ydif + &
        ((nu+nutt)*(w(i, j, k+1)-w(i, j, k))/zdif- &
        (nu+nutb)*(w(i, j, k)-w(i, j, k-1))/zdif)/zdif
      end do
      end do
      end do
    else
      print*, 'tag error in diff!'
      stop
    end if
  end subroutine diff

  subroutine sour(tag, sorRes)
    ! calculate source terms
    implicit none
    ! public
    integer :: tag
    ! local
    integer i, j, k
    real(dpc) :: nute, nutw, nutn, nuts, nutt, nutb
    real(dpc) :: fppx, fppy, fppz
    !
    if (ira==0) return
    if (tag==1) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      fppx = (fptx(i, j, k) + fptx(i-1, j, k))/2.0d0
      nute = nut(i, j, k)
      nutw = nut(i-1, j, k)
      nutn = (nut(i, j, k)+nut(i, j+1, k)+nut(i-1, j+1, k)+nut(i-1, j, k))/4.0d0
      nuts = (nut(i, j, k)+nut(i, j-1, k)+nut(i-1, j-1, k)+nut(i-1, j, k))/4.0d0
      nutt = (nut(i, j, k)+nut(i, j, k+1)+nut(i-1, j, k+1)+nut(i-1, j, k))/4.0d0
      nutb = (nut(i, j, k)+nut(i, j, k-1)+nut(i-1, j, k-1)+nut(i-1, j, k))/4.0d0
      sorRes(i, j, k) = (nute*(u(i+1, j, k)-u(i, j, k))/xdif- &
        nutw*(u(i, j, k)-u(i-1, j, k))/xdif)/xdif + &
        (nutn*(v(i, j+1, k)-v(i-1, j+1, k))/xdif- &
        nuts*(v(i, j, k)-v(i-1, j, k))/xdif)/ydif + &
        (nutt*(w(i, j, k+1)-w(i-1, j, k+1))/xdif- &
        nutb*(w(i, j, k)-w(i-1, j, k))/xdif)/zdif + fppx/(xdif*ydif*zdif*rho*phirho(i, j, k))
      end do
      end do
      end do
    else if (tag==2) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      fppy = (fpty(i, j, k) + fpty(i, j-1, k))/2.0d0
      nute = (nut(i, j, k)+nut(i+1, j, k)+nut(i+1, j-1, k)+nut(i ,j-1, k))/4.0d0
      nutw = (nut(i, j, k)+nut(i-1, j, k)+nut(i-1, j-1, k)+nut(i ,j-1, k))/4.0d0
      nutn = nut(i, j, k)
      nuts = nut(i, j-1, k)
      nutt = (nut(i, j, k)+nut(i, j, k+1)+nut(i, j-1, k+1)+nut(i, j-1, k))/4.0d0
      nutb = (nut(i, j, k)+nut(i, j, k-1)+nut(i, j-1, k-1)+nut(i, j-1, k))/4.0d0
      sorRes(i, j, k) = (nute*(u(i+1, j, k)-u(i+1, j-1, k))/ydif- &
        nutw*(u(i, j, k)-u(i, j-1, k))/ydif)/xdif + &
        (nutn*(v(i, j+1, k)-v(i, j, k))/ydif- &
        nuts*(v(i, j, k)-v(i, j-1, k))/ydif)/ydif + &
        (nutt*(w(i, j, k+1)-w(i, j-1, k+1))/ydif- &
        nutb*(w(i, j, k)-w(i, j-1, k))/ydif)/zdif + fppy/(xdif*ydif*zdif*rho*phirho(i, j, k))
      end do
      end do
      end do
    else if (tag==3) then
      do k = sz, ez
      do j = sy, ey
      do i = sx, ex
      fppz = (fptz(i, j, k) + fptz(i, j, k-1))/2.0d0
      nute = (nut(i, j, k)+nut(i+1, j, k)+nut(i+1, j, k-1)+nut(i ,j, k-1))/4.0d0
      nutw = (nut(i, j, k)+nut(i-1, j, k)+nut(i-1, j, k-1)+nut(i ,j, k-1))/4.0d0
      nutn = (nut(i, j, k)+nut(i, j+1, k)+nut(i, j+1, k-1)+nut(i, j, k-1))/4.0d0
      nuts = (nut(i, j, k)+nut(i, j-1, k)+nut(i, j-1, k-1)+nut(i, j, k-1))/4.0d0
      nutt = nut(i, j, k)
      nutb = nut(i, j, k-1)
      sorRes(i, j, k) = (nute*(u(i+1, j, k)-u(i+1, j, k-1))/zdif- &
        nutw*(u(i, j, k)-u(i, j, k-1))/zdif)/xdif + &
        (nutn*(v(i, j+1, k)-v(i, j+1, k-1))/zdif- &
        nuts*(v(i, j, k)-v(i, j, k-1))/zdif)/ydif + &
        (nutt*(w(i, j, k+1)-w(i, j, k))/zdif- &
        nutb*(w(i, j, k)-w(i, j, k-1))/zdif)/zdif + fppz/(xdif*ydif*zdif*rho*phirho(i, j, k))
      end do
      end do
      end do
    else
      print*, 'tag error in sour!'
      stop
    end if
  end subroutine sour

  subroutine u_interpolation
    ! fix nut u v w near IM-boundary
    use grid_exchange, only : ubdy, nutbdy
    implicit none
    ! local
    integer :: n
    real(dpc), dimension(:, :, :), pointer :: valp
    real(dpc), dimension(:, :, :, :), pointer :: ctvalp
    real(dpc), dimension(:, :, :), pointer :: utaop
    real(dpc), dimension(:, :), pointer :: utaosfp
    type (vector), dimension(:), pointer :: norp
    real(dpc), dimension(:), pointer :: ipdp
    !
    utao = cZero
    utaosf = cZero
    utaop => utao
    utaosfp => utaosf
    ctvalp => ctval
    do n = 1, nnif
    valp => nut
    norp => norpp
    ipdp => ipdisw
    call center_interpolation(n, norp, ipdp, ctvalp, valp, utaop, utaosfp)
    valp => uf
    norp => norpu
    ipdp => ipdisu
    call interpolation(4, n, norp, ipdp, valp, utaop)
    valp => vf
    norp => norpv
    ipdp => ipdisv
    call interpolation(5, n, norp, ipdp, valp, utaop)
    valp => wf
    norp => norpw
    ipdp => ipdisw
    call interpolation(6, n, norp, ipdp, valp, utaop)
    end do
    call nutbdy(nut)
    call ubdy(uf, vf, wf)
  end subroutine u_interpolation

  subroutine u_cal
    ! calculate intermediate u v w
    implicit none
    ! local
    integer :: i, j, k, nnk
    integer, dimension(:, :, :), pointer :: flgp
    real(dpc) :: apha, fx
    real(dpc), dimension(nxdim, my, mz) :: flo, dif, sor
    real(dpc), dimension(mz) :: aw, ae, as, an, ab, at, ap, rhs
    real(dpc), dimension(mz) :: omga, gama, beta
    real(dpc), dimension(nxdim, my, mz, 3):: con
    real(dpc), dimension(:, :, :), pointer :: valp
    real(dpc), dimension(:, :, :), pointer :: valfp
    real(dpc), dimension(:, :, :), pointer :: dpdp
    !
    apha = 1.0d0
    level1: do nnk = 1, 3
    if (nnk==1) then
      flgp => ffu
      valp => u
      valfp => uf
      dpdp => dpdx
    else if (nnk==2) then
      flgp => ffv
      valp => v
      valfp => vf
      dpdp => dpdy
    else if (nnk==3) then
      flgp => ffw
      valp => w
      valfp => wf
      dpdp => dpdz
    end if
    level2: do j = sy, ey
    level3: do i = sx, ex
    level4: do k = sz, ez
    if (flgp(i, j, k)/=0) then
      fx = (valfp(i, j, k)-valp(i, j, k))/dt - rhsx(i, j, k, nnk)
      dpdp(i, j, k) = cZero
    else
      fx = cZero
    end if
    con(i, j, k, nnk) = -1.5*flox(i, j, k, nnk) + 0.5d0*floxx(i, j, k, nnk) + &
      1.5*sorx(i, j, k, nnk) - 0.5d0*sorxx(i, j, k, nnk) + 0.5d0*difx(i, j, k, nnk) - &
      dpdp(i, j, k)/(rho*phirho(i, j, k)) + fx
    if (nnk==1) then
      aw(k) = 0.5d0*(nu+nut(i-1, j, k))/xdif**2*dt
      ae(k) = 0.5d0*(nu+nut(i, j, k))/xdif**2*dt
      as(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j-1, k)+nut(i-1, j-1, k)+nut(i-1, j, k)))/ydif**2*dt
      an(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j+1, k)+nut(i-1, j+1, k)+nut(i-1, j, k)))/ydif**2*dt
      ab(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j, k-1)+nut(i-1, j, k-1)+nut(i-1, j, k)))/zdif**2*dt
      at(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j, k+1)+nut(i-1, j, k+1)+nut(i-1, j, k)))/zdif**2*dt
    else if (nnk==2) then
      aw(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i-1, j, k)+nut(i-1, j-1, k)+nut(i, j-1, k)))/xdif**2*dt
      ae(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i+1, j, k)+nut(i+1, j-1, k)+nut(i, j-1, k)))/xdif**2*dt
      as(k) = 0.5d0*(nu+nut(i, j-1, k))/ydif**2*dt
      an(k) = 0.5d0*(nu+nut(i, j, k))/ydif**2*dt
      ab(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j, k-1)+nut(i, j-1, k-1)+nut(i, j-1, k)))/zdif**2*dt
      at(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j, k+1)+nut(i, j-1, k+1)+nut(i, j-1, k)))/zdif**2*dt
    else if (nnk==3) then
      aw(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i-1, j, k)+nut(i-1, j, k-1)+nut(i, j, k-1)))/xdif**2*dt
      ae(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i+1, j, k)+nut(i+1, j, k-1)+nut(i, j, k-1)))/xdif**2*dt
      as(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j-1, k)+nut(i, j-1, k-1)+nut(i, j, k-1)))/ydif**2*dt
      an(k) = 0.5d0*(nu+0.25d0*(nut(i, j, k)+nut(i, j+1, k)+nut(i, j+1, k-1)+nut(i, j, k-1)))/ydif**2*dt
      ab(k) = 0.5d0*(nu+nut(i, j, k-1))/zdif**2*dt
      at(k) = 0.5d0*(nu+nut(i, j, k))/zdif**2*dt
    end if
    ap(k) = 1.0d0 + ab(k) + at(k) + aw(k) + ae(k) + as(k) + an(k)
    rhs(k) = (con(i, j, k, nnk)*dt + valp(i, j, k) + &
      aw(k)*valp(i-1, j, k)+ae(k)*valp(i+1, j, k)+as(k)*valp(i, j-1, k)+an(k)*valp(i, j+1, k))*apha + &
      (1.0d0-apha)*valp(i, j, k)*ap(k)
    end do level4
    if (nnk==1) then
      ab(ez+1) = 0.5d0*(nu+0.25d0*(nut(i, j, ez+1)+nut(i, j, ez)+nut(i-1, j, ez)+nut(i-1, j, ez+1)))/zdif**2*dt
      at(sz-1) = 0.5d0*(nu+0.25d0*(nut(i, j, sz-1)+nut(i, j, sz)+nut(i-1, j, sz)+nut(i-1, j, sz-1)))/zdif**2*dt
    else if (nnk==2) then
      ab(ez+1) = 0.5d0*(nu+0.25d0*(nut(i, j, ez+1)+nut(i, j, ez)+nut(i, j-1, ez)+nut(i, j-1, ez+1)))/zdif**2*dt
      at(sz-1) = 0.5d0*(nu+0.25d0*(nut(i, j, sz-1)+nut(i, j, sz)+nut(i, j-1, sz)+nut(i, j-1, sz-1)))/zdif**2*dt
    else if (nnk==3) then
      ab(ez+1) = 0.5d0*(nu+nut(i, j, ez))/zdif**2*dt
      at(sz-1) = 0.5d0*(nu+nut(i, j, sz-1))/zdif**2*dt
    end if
    ap(sz-1) = 1.0d0 + ab(sz-1) + at(sz-1) + aw(sz-1) + ae(sz-1) + as(sz-1) + an(sz-1)
    ap(ez+1) = 1.0d0 + ab(ez+1) + at(ez+1) + aw(ez+1) + ae(ez+1) + as(ez+1) + an(ez+1)
    rhs(sz-1) = (con(i, j, sz-1, nnk)*dt + valp(i, j, sz-1) + &
      aw(sz-1)*valp(i-1, j, sz-1)+ae(sz-1)*valp(i+1, j, sz-1)+ &
      as(sz-1)*valp(i, j-1, sz-1)+an(sz-1)*valp(i, j+1, sz-1))*apha + &
      (1.0d0-apha)*valp(i, j, sz-1)*ap(sz-1)
    rhs(ez+1) = (con(i, j, ez+1, nnk)*dt + valp(i, j, ez+1) + &
      aw(ez+1)*valp(i-1, j, ez+1)+ae(ez+1)*valp(i+1, j, ez+1)+ &
      as(ez+1)*valp(i, j-1, ez+1)+an(ez+1)*valp(i, j+1, ez+1))*apha + &
      (1.0d0-apha)*valp(i, j, ez+1)*ap(ez+1)
    omga(sz-1) = -at(sz-1)
    gama(sz-1) = rhs(sz-1)
    level4: do k = sz, ez+1
    beta(k) = -ab(k)/omga(k-1)
    omga(k) = ap(k) + at(k-1)*beta(k)
    gama(k) = rhs(k) - beta(k)*gama(k-1)
    end do level4
    valp(i, j, ez+1) = gama(ez+1)/omga(ez+1) !Line G-S & Line SOR
    level4: do k = ez, sz-1, -1
    valp(i, j, k) = (gama(k)+at(k)*valp(i, j, k+1))/omga(k) !Line G-S & Line SOR
    if (isnan(valp(i, j, k))) then
      print*, 'valp = nan when nnk=', nnk, 'omga=', omga(k), 'gama=', gama(k)
      stop
    end if
    end do level4
    end do level3
    end do level2
    end do level1
  end subroutine u_cal

  subroutine ux_correction
    use grid_exchange, only : ubdy
    implicit none
    ! local
    integer i, j, k
    !
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    u(i, j, k) = u(i, j, k) + dpdx(i, j, k)/(rho*phirho(i, j, k))*dt
    v(i, j, k) = v(i, j, k) + dpdy(i, j ,k)/(rho*phirho(i, j, k))*dt
    w(i, j, k) = w(i, j, k) + dpdz(i, j, k)/(rho*phirho(i, j, k))*dt
    if (isnan(u(i, j, k)) .or. u(i, j, k)>1.0d8) then
      print*, dpdx(i, j, k), u(i, j, k), 'uxcor_u'
      stop
    end if
    if (isnan(v(i, j, k)) .or. v(i, j, k)>1.0d8) then
      print*, dpdy(i, j, k), v(i, j, k), 'uxcor_v'
      stop
    end if
    if (isnan(w(i, j, k)) .or. w(i, j, k)>1.0d8) then
      print*, dpdz(i, j, k), w(i, j, k), 'uxcor_w'
      stop
    end if
    end do
    end do
    end do
    call ubdy(u, v, w)
  end subroutine ux_correction

  subroutine u_correction
    use grid_exchange, only : ubdy
    implicit none
    ! local
    integer i, j, k
    !
    do k = sz, ez
    do j = sy, ey
    do i = sx, ex
    u(i, j, k) = u(i, j, k) - (p(i, j, k)-p(i-1, j, k))/xdif/(rho*phirho(i, j, k))*dt
    v(i, j, k) = v(i, j, k) - (p(i, j, k)-p(i, j-1, k))/ydif/(rho*phirho(i, j, k))*dt
    w(i, j, k) = w(i, j, k) - (p(i, j, k)-p(i, j, k-1))/zdif/(rho*phirho(i, j, k))*dt
    if (isnan(u(i, j, k)) .or. u(i, j, k)>1.0d8) then
      print*, p(i, j, k)-p(i-1, j, k), 'ucor_u'
      stop
    end if
    if (isnan(v(i, j, k)) .or. v(i, j, k)>1.0d8) then
      print*, p(i, j, k)-p(i, j-1, k), 'ucor_v'
      stop
    end if
    if (isnan(w(i, j, k)) .or. w(i, j, k)>1.0d8) then
      print*, p(i, j, k)-p(i, j, k-1), 'ucor_w'
      stop
    end if
    end do
    end do
    end do
    call ubdy(u, v, w)
  end subroutine u_correction
end module fluid
