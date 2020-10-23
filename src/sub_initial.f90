!-----------------------------------------------------------------------
subroutine sub_initial
!-----------------------------------------------------------------------
  use mod_variables
  implicit none
!
  integer :: ic,jc,icen,jcen,ib,ierr,jmax_tmp,ios,na,nb
  real(8) :: rblank
  real(4) :: xout(imax),yout(jmax)
  integer :: nblock = 1
  real(4) :: re4,step4,time4
  type(variances) vars
!-----------------------------------------------------------------------
!
  ostep     = 0
  tstep     = 0.0
!
!-Variables
  allocate(ustg(imx1,jmax),stat=ierr)
  allocate(vstg(imax,jmx1),stat=ierr)
  allocate(udlt(imx1,jmax),stat=ierr)
  allocate(vdlt(imax,jmx1),stat=ierr)
  allocate(pcnt(imax,jmax),stat=ierr)
  allocate(divg(imax,jmax),stat=ierr)
  allocate(iblk(imax,jmax),stat=ierr)
!
  allocate(uinf(jmax)     ,stat=ierr)
!
  allocate(xcen(imx1)     ,stat=ierr)
  allocate(ycen(jmx1)     ,stat=ierr)
!
  allocate(q4(imax,jmax,4),stat=ierr)
!
! Specify the number of points per one side of square cylinder
  jobj = int(0.100*jmax)  ! 8 points per one side when jmax = 80
! jobj = int(0.064*jmax)  ! jmax = 320: the original setting 2018/12/26
  iobj = jobj
  dx   = 1.0/dble(jobj)
!
  do ic=1,imx1
    xcen(ic) = dx* dble(ic-1-0.25*imax)
  enddo
  do ic=1,imax
    xout(ic) = dx*(dble(ic-0.25*imax)-0.5d0)
  enddo
  do jc=1,jmx1
    ycen(jc) = dx* dble(jc-1-0.50*jmax)
  enddo
  do jc=1,jmax
    yout(jc) = dx*(dble(jc-0.50*jmax)-0.5d0)
  enddo
!
!-Definition of object (a rectangular cylinder)
  icen = int(0.25d0*dble(imax))
  jcen = int(0.50d0*dble(jmax))
!
  do jc=1,jmax
  do ic=1,imax
    if(abs(xout(ic))<=0.5d0 .and. abs(yout(jc))<=0.5d0) then
      iblk(ic,jc) = 0
    else
      iblk(ic,jc) = 1
    endif
  enddo
  enddo
!
  do jc=2,jmax-1
  do ic=2,imax-1
    if(iblk(ic,jc)==0) then
      if(iblk(ic-1,jc)==1 .and. iblk(ic,jc-1)==1) then
        corner(1,1) = ic
        corner(2,1) = jc
      endif
      if(iblk(ic+1,jc)==1 .and. iblk(ic,jc-1)==1) then
        corner(1,2) = ic
        corner(2,2) = jc
      endif
      if(iblk(ic-1,jc)==1 .and. iblk(ic,jc+1)==1) then
        corner(1,3) = ic
        corner(2,3) = jc
      endif
      if(iblk(ic+1,jc)==1 .and. iblk(ic,jc+1)==1) then
        corner(1,4) = ic
        corner(2,4) = jc
      endif
    endif
  enddo
  enddo
!
  do jc=1,jmax
  do ic=1,imax
    rblank = real(iblk(ic,jc))
    ustg(ic  ,jc) = u_inf*rblank
    vstg(ic,jc  ) = v_inf*rblank
    ustg(ic+1,jc) = u_inf*rblank
    vstg(ic,jc+1) = v_inf*rblank
!
    udlt(ic  ,jc) = 0.0
    vdlt(ic,jc  ) = 0.0
    udlt(ic+1,jc) = 0.0
    vdlt(ic,jc+1) = 0.0
  enddo
  enddo
!
! Inflow velocity
  uinf(:) = u_inf
!
  do jc=1,jmax
  do ic=1,imax
    pcnt(ic,jc) = 0.0
    divg(ic,jc) = 0.0
  enddo
  enddo
!
!!cfl = 0.2
!!call sub_timestep(dt)
  call sub_cflnum(cfl,dt)
  call sub_bc_outer(ustg,vstg,pcnt)
  call sub_bc_wall(ustg,vstg,pcnt)
!
!
! Output mesh
  if(ioutput==1) then
    open(30,file='./output/mesh.g',form='unformatted',status='replace',action='write')
    nblock = 1
    write(30) nblock
    write(30) (imax,jmax,ib=1,nblock)
    do ib=1,nblock
      write(30) ((xout(ic)   ,ic=1,imax),jc=1,jmax), &
                ((yout(jc)   ,ic=1,imax),jc=1,jmax), &
                ((iblk(ic,jc),ic=1,imax),jc=1,jmax)
    enddo
    close(30)
    open(30,file='./statis/mesh.g',form='unformatted',status='replace',action='write')
    nblock = 1
    write(30) nblock
    write(30) (imax,jmax,ib=1,nblock)
    do ib=1,nblock
      write(30) ((xout(ic)   ,ic=1,imax),jc=1,jmax), &
                ((yout(jc)   ,ic=1,imax),jc=1,jmax), &
                ((iblk(ic,jc),ic=1,imax),jc=1,jmax)
    enddo
    close(30)
  endif
!
  return
end subroutine sub_initial
!-----------------------------------------------------------------------
!
!
!
! ----------------------------------------------------------------------
subroutine sub_restart
!-----------------------------------------------------------------------
  use mod_variables
  implicit none

  integer :: ic,jc,m,imx,jmx,ib,nblock
  real(4) :: re4,step4,time4
!-----------------------------------------------------------------------
!
!
!-Read an instantaneous field
  open(15,file='./restart.q',form='unformatted',action='read',status='old')
  nblock = 1
  read(15) nblock
  read(15) (imx,jmx,ib=1,nblock)

  ! consistency check >>>>>>
  if(imx/=imax .or. jmx/=jmax) stop 'stop imx/=imax or jmx/=jmax'

  do ib=1,nblock
    read(15) re4,step4,re4,time4
    read(15) (((q4(ic,jc,m),ic=1,imax),jc=1,jmax),m=1,4)
  enddo
  close(15)
!
  ostep = int(step4)
  tstep = real(time4)
!
  do jc=1,jmax
    do ic=2,imax
      ustg(ic,jc) = real(q4(ic,jc,2)+q4(ic-1,jc,2))*0.5d0
    enddo
    ustg(1   ,jc) = real(q4(1   ,jc,2))
    ustg(imx1,jc) = real(q4(imax,jc,2))
  enddo
  do ic=1,imax
    do jc=2,jmax
      vstg(ic,jc) = real(q4(ic,jc,3)+q4(ic,jc-1,3))*0.5d0
    enddo
    vstg(ic,1   ) = real(q4(ic,1   ,3))
    vstg(ic,jmx1) = real(q4(ic,jmax,3))
  enddo
  do jc=1,jmax
  do ic=1,imax
    pcnt(ic,jc) = real(q4(ic,jc,4))
  enddo
  enddo
!
! Inlet boundary condition
  do jc=1,jmax
    uinf(jc) = ustg(1,jc)
  enddo
!
  return
end subroutine sub_restart
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
subroutine sub_addvtx(xc,yc,ga)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,xcen,ycen,ustg,vstg,iblk
  implicit none
!
  integer :: ic,jc,icen,jcen,ib,nblock,ierr,jmax_tmp,ios
  real(8),parameter :: rc = 0.5d0
!!real(8),parameter :: ga = 1.0d0 
  real(8) :: rblank,gp,xc,yc,rr,yl,ga,pi
!-----------------------------------------------------------------------
!
!
  pi = 4.d0*datan(1.d0)
  gp = ga/(2.d0*pi)
  yl = ycen(jmax)-ycen(1)
!
  do jc=1,jmax
  do ic=1,imax
    rblank = real(iblk(ic,jc))
!   Main vortex
    rr     = (xcen(ic)-xc)**2+(ycen(jc)-yc)**2+rc**2
    ustg(ic,jc) = ustg(ic,jc)-(ga/rr)*(ycen(jc)-yc)*rblank
    vstg(ic,jc) = vstg(ic,jc)+(ga/rr)*(xcen(ic)-xc)*rblank
!   Lower mirror vortex
    rr     = (xcen(ic)-xc)**2+(ycen(jc)-yc-yl)**2+rc**2
    ustg(ic,jc) = ustg(ic,jc)+(ga/rr)*(ycen(jc)-yc)*rblank
    vstg(ic,jc) = vstg(ic,jc)-(ga/rr)*(xcen(ic)-xc)*rblank
!   Upper mirror vortex
    rr     = (xcen(ic)-xc)**2+(ycen(jc)-yc+yl)**2+rc**2
    ustg(ic,jc) = ustg(ic,jc)+(ga/rr)*(ycen(jc)-yc)*rblank
    vstg(ic,jc) = vstg(ic,jc)-(ga/rr)*(xcen(ic)-xc)*rblank
  enddo
  enddo
!
  return
end subroutine sub_addvtx
!-----------------------------------------------------------------------
