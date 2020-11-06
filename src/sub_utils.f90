!-----------------------------------------------------------------------
! Generate initial ensemble
!-----------------------------------------------------------------------
subroutine sub_ensemble_gen(vars)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,dt,variances,n_var,n_ens, &
                           ef,ave,ustg,vstg,pcnt,uini,vini
  use m_random3
  use m_ranmean3
  implicit none
!
  type(variances) :: vars
  integer(4):: seed 
  integer :: i,j,k,ic,jc,iens,istep,ind2,ierr
  integer :: oneperiod,steplag
  real(8) :: xc,yc,ga,w1i
  real(8),allocatable :: wn(:,:)
  real(8),allocatable :: XX(:,:),SS(:,:),w(:),rv1(:),pmode(:,:)
!-----------------------------------------------------------------------
!
!
! allocate(wn(n_var,n_ens))
  allocate(wn(n_var,1))
  if(vars%ivr==1) then
   !call random3(wn,n_var,n_ens)            ! N(0,var%ini) white noise
    do i=1,n_ens
      call random3(wn,n_var,1)              ! N(0,var%ini) white noise
      do j=1,n_var
       !ef(j,i) = ef(j,i)+dsqrt(vars%ini)*wn(j,i)
        ef(j,i) = ef(j,i)+dsqrt(vars%ini)*wn(j,1)
      enddo
    enddo
  elseif(vars%ivr==2) then
    oneperiod = int(7.d0/dt)
    steplag   = int(oneperiod/(n_ens))
    iens      = 0
    do istep = 1,oneperiod
!-----RHS
      call sub_bc_outer(ustg,vstg,pcnt)
      call sub_bc_wall(ustg,vstg,pcnt)
      call sub_rhs
!-----Pressure
      call sub_HSMAC(ustg,vstg,pcnt)
      if(mod(istep,steplag)==0 .or. istep==1) then
        iens = iens+1
        if(iens>n_ens+1) exit
        ind2 = 0
        do jc=1,jmax
        do ic=1,imx1
          ind2 = ind2+1
          ef(ind2,iens) = ustg(ic,jc)
        enddo
        enddo
        do jc=1,jmx1
        do ic=1,imax
          ind2 = ind2+1
          ef(ind2,iens) = vstg(ic,jc)
        enddo
        enddo
      endif
    enddo
  elseif(vars%ivr==3) then
    deallocate(wn)
    allocate(wn(n_ens,3))
    call random_number(wn(:,1))
    call random_number(wn(:,2))
    call random_number(wn(:,3))
    do iens=1,n_ens
      ustg(1:imx1,1:jmax) = uini(1:imx1,1:jmax)
      vstg(1:imax,1:jmx1) = vini(1:imax,1:jmx1)
      xc = 4.d0*wn(iens,1)-6.0d0   ![-6.0,-2.0]
      yc = 1.d0*wn(iens,2)-0.5d0   ![-0.5, 0.5]
      write(*,*) "Initial vortex position (x,y):",xc,yc
!!    ga = 0.6d0-0.2*(wn(iens,3))
      ga = 1.d0-0.2*(wn(iens,3))   ![ 0.8, 1.0]
      call sub_addvtx(xc,yc,ga)
      ind2 = 0
      do jc=1,jmax
      do ic=1,imx1
        ind2 = ind2+1
        ef(ind2,iens) = ustg(ic,jc)
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        ind2 = ind2+1
        ef(ind2,iens) = vstg(ic,jc)
      enddo
      enddo
    enddo
  elseif(vars%ivr==4) then
    oneperiod = int(7.d0/dt)
    steplag   = int(oneperiod/(n_ens))
    iens      = 0
    do istep = 1,oneperiod
!-----RHS
      call sub_bc_outer(ustg,vstg,pcnt)
      call sub_bc_wall(ustg,vstg,pcnt)
      call sub_rhs
!-----Pressure
      call sub_HSMAC(ustg,vstg,pcnt)
!
!     if(mod(istep,steplag)==0 .or. istep==1) then
      if(mod(istep,steplag)==0) then
        iens = iens+1
        if(iens>=n_ens+1) exit
        ind2 = 0
        do jc=1,jmax
        do ic=1,imx1
          ind2 = ind2+1
          ef(ind2,iens) = ustg(ic,jc)
        enddo
        enddo
        do jc=1,jmx1
        do ic=1,imax
          ind2 = ind2+1
          ef(ind2,iens) = vstg(ic,jc)
        enddo
        enddo
      endif
    enddo
    call ranmean3(ef,ave,n_var,n_ens)
    do j=1,n_ens
      do i=1,n_var
        ef(i,j) = ef(i,j)-ave(i)
      enddo
    enddo
    allocate(XX(n_ens,n_ens),SS(n_ens,n_ens),w(n_ens))
    allocate(rv1(n_ens),pmode(n_var,n_ens))
   !Compose a matrix (X^T)(X) for offline POD
    do k=1,n_ens
      do j=1,n_ens
        XX(j,k) = 0.d0
        do i=1,n_var
          XX(j,k) = XX(j,k)+ef(i,j)*ef(i,k)
        enddo
      enddo
    enddo
    w(:)=0.d0 ; SS(:,:) = 0.d0
    call svd(n_ens,n_ens,n_ens,XX,w,.true.,XX,.true.,SS,ierr,rv1)
    do k=1,n_ens
      w1i = 1.d0/dsqrt(w(k))
      do i=1,n_var
        pmode(i,k) = 0.d0
        do j=1,n_ens
          pmode(i,k) = pmode(i,k)+ef(i,j)*XX(j,k)
        enddo
        pmode(i,k) = pmode(i,k)*w1i
      enddo
     !write(*,*) w(k),minval(pmode(:,k)),maxval(pmode(:,k))
    enddo
    do k=1,n_ens
      do i=1,n_var
        ef(i,k) = ef(i,n_ens+1)+dsqrt(vars%ini)*pmode(i,k)
      enddo
    enddo
    deallocate(XX,SS,w,rv1,pmode)
!
  endif
  deallocate(wn)
!
  return
end subroutine sub_ensemble_gen
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Calculate summation of energy in entire flow field
!-----------------------------------------------------------------------
subroutine sub_energy(sum)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,dx,ustg,vstg
  implicit none
!
  integer :: ic,jc
  real(8) :: vol
  real(8) :: sum 
!-----------------------------------------------------------------------
!
  sum = 0.0
  vol = dx**2
!
  do jc = 1,jmax
  do ic = 1,imax
    sum = sum+vol*0.5*(ustg(ic,jc)**2+vstg(ic,jc)**2)
  enddo
  enddo
!
  return
end subroutine sub_energy
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Define time step size based on CFL number
!-----------------------------------------------------------------------
subroutine sub_timestep(dt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,ustg,vstg,dx,cfl
  implicit none
!
  integer :: ic,jc
  real(8) :: dt
!-----------------------------------------------------------------------
!
  dt  = 1.0
!
  do jc = 1,jmax
  do ic = 1,imax
    dt = min(dt,cfl*dx/(abs(ustg(ic,jc))+abs(vstg(ic,jc))))
  enddo
  enddo
!
  return
end subroutine sub_timestep
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Calculate maximum value of local CFL number
!-----------------------------------------------------------------------
subroutine sub_cflnum(cfl,dt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,ustg,vstg,dx
  implicit none
!
  integer :: ic,jc
  real(8) :: cfl,dt,dxi
!-----------------------------------------------------------------------
!
  cfl = 0.0
  dxi = 1.0/dx
!
  do jc = 1,jmax
  do ic = 1,imax
    cfl = max(cfl,(abs(ustg(ic,jc))+abs(vstg(ic,jc)))*dt*dxi)
  enddo
  enddo
!
  return
end subroutine sub_cflnum
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Probing u-velocity
!-----------------------------------------------------------------------
subroutine sub_probe(istep,iter,char)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,xcen,ycen,ustg,vstg,iprb,jprb,prob,nstep,niter,ihx
  implicit none
!
  real(8),parameter    :: u_inf = 1.0

  integer :: ic,jc,istep,ierr,iter,ip
  real(8) :: xx,yy,xt,yt,dist,dmin
  character(len=3) :: char
!-----------------------------------------------------------------------
!
!
  if(istep==1 .and. iter==1) then
!
   if(allocated(prob)) deallocate(prob)
   allocate(prob(4,nstep),stat=ierr)
   prob(:,:) = 0.0
!
!
!   Before: (3.0,0.0), (3.0,1.0), (3.0,1.0)
!   Position 1
    xt   =-2.5d0 
    yt   =-0.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(1) = ic
          jprb(1) = jc
        endif
!     endif
    enddo
    enddo
!   Position 2
    xt   = 0.0d0 
    yt   =-2.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(2) = ic
          jprb(2) = jc
        endif
!     endif
    enddo
    enddo
!   Position 3
    xt   = 2.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(3) = ic
          jprb(3) = jc
        endif
!     endif
    enddo
    enddo
!   Position 4
    xt   = 7.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(4) = ic
          jprb(4) = jc
        endif
!     endif
    enddo
    enddo
    jprb(4) = jprb(3)
!
    if(char=="ref") then
      open(60,file='probe_ref.dat',form='formatted',status='replace',action='write')
    elseif(char=="est") then
      open(60,file='probe_est.dat',form='formatted',status='replace',action='write')
    endif
    write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
    write(60,'(" iter   step      ref1        ref2        ref3        ref4")')
    close(60)
!
  endif
!
! prob(1,istep) = dsqrt(ustg(iprb(1),jprb(1))**2+vstg(iprb(1),jprb(1))**2)
! prob(2,istep) = dsqrt(ustg(iprb(2),jprb(2))**2+vstg(iprb(2),jprb(2))**2)
! prob(3,istep) = dsqrt(ustg(iprb(3),jprb(3))**2+vstg(iprb(3),jprb(3))**2)
  prob(1,istep) = ustg(iprb(1),jprb(1))
  prob(2,istep) = ustg(iprb(2),jprb(2))
  prob(3,istep) = ustg(iprb(3),jprb(3))
  prob(4,istep) = ustg(iprb(4),jprb(4))
!
  if(istep==nstep) then
    if(char=="ref") then
      open(60,file='probe_ref.dat',form='formatted',position='append')
    elseif(char=="est") then
      open(60,file='probe_est.dat',form='formatted',position='append')
    endif
    do ic=1,nstep
      write(60,'(2i6,4f12.7)') iter,ic+nstep*(iter-1),prob(1,ic),prob(2,ic),prob(3,ic),prob(4,ic)
    enddo
    close(60)
  endif
!
  return
end subroutine sub_probe
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Probing u-velocity
!-----------------------------------------------------------------------
subroutine sub_probe_ens(iter)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,xcen,ycen,ustg,vstg,nstep,niter, &
                           iprb,jprb,prob,                                      &
                           n_ens,n_mes,rm,ihx,uens,vens,uref,vref
  implicit none
!
  integer :: ic,jc,istep,ierr,iter,ip,iens
  real(8) :: ave(4,nstep),var(4,nstep),std(4,nstep),dist,dmin,xx,yy,xt,yt
!-----------------------------------------------------------------------
!
!
  if(iter==1) then
!
!   Position 1
    xt   =-2.5d0 
    yt   =-0.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(1) = ic
          jprb(1) = jc
        endif
!     endif
    enddo
    enddo
!   Position 2
    xt   = 0.0d0 
    yt   =-2.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(2) = ic
          jprb(2) = jc
        endif
!     endif
    enddo
    enddo
!   Position 3
    xt   = 2.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(3) = ic
          jprb(3) = jc
        endif
!     endif
    enddo
    enddo
!   Position 4
    xt   = 7.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(4) = ic
          jprb(4) = jc
        endif
!     endif
    enddo
    enddo
    jprb(4) = jprb(3)
!
    open(60,file='probe_ens.dat',form='formatted',status='replace',action='write')
    write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
    write(60,'(" iter   step      ave1        std1        ave2        std2        ave3        std3        ave4        std4")')
    close(60)
    open(60,file='probe_mes.dat',form='formatted',status='replace',action='write')
    write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
    write(60,'(" iter   step      mes1        mes2        mes3        mes4")')
    close(60)
!
  endif
!
  do istep=1,nstep
    ave(1:4,istep) = 0.d0
    var(1:4,istep) = 0.d0
    std(1:4,istep) = 0.d0
    do iens=1,n_ens
      ave(1,istep) = ave(1,istep)+uens(iprb(1),jprb(1),istep,iens)
      ave(2,istep) = ave(2,istep)+uens(iprb(2),jprb(2),istep,iens)
      ave(3,istep) = ave(3,istep)+uens(iprb(3),jprb(3),istep,iens)
      ave(4,istep) = ave(4,istep)+uens(iprb(4),jprb(4),istep,iens)
    enddo
    ave(1:4,istep) = ave(1:4,istep)/dble(n_ens)
    do iens=1,n_ens
      var(1,istep) = var(1,istep)+(uens(iprb(1),jprb(1),istep,iens)-ave(1,istep))**2
      var(2,istep) = var(2,istep)+(uens(iprb(2),jprb(2),istep,iens)-ave(2,istep))**2
      var(3,istep) = var(3,istep)+(uens(iprb(3),jprb(3),istep,iens)-ave(3,istep))**2
      var(4,istep) = var(4,istep)+(uens(iprb(4),jprb(4),istep,iens)-ave(4,istep))**2
    enddo
    var(1:4,istep) = var(1:4,istep)/dble(n_ens)
    std(1:4,istep) = dsqrt(var(1:4,istep))
  enddo
!
  open(60,file='probe_ens.dat',form='formatted',position='append')
  do istep=1,nstep
    write(60,'(2i6,8f12.7)') iter,istep+nstep*(iter-1),(ave(ip,istep),std(ip,istep),ip=1,4)
  enddo
  close(60)
  open(60,file='probe_mes.dat',form='formatted',position='append')
  write(60,'(2i6,4f12.7)') iter,nstep*iter,(uref(iprb(ip),jprb(ip),nstep*iter),ip=1,4)
  close(60)
!
  return
end subroutine sub_probe_ens
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
! Probing u-velocity
!-----------------------------------------------------------------------
subroutine sub_probe_ens_all(iter)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,xcen,ycen,ustg,vstg,nstep,niter, &
                           iprb,jprb,prob,                                      &
                           n_ens,n_mes,rm,ihx,uens,vens,uref,vref
  implicit none
!
  integer :: ic,jc,istep,ierr,iter,ip,iens
  real(8) :: ave(4,nstep),var(4,nstep),std(4,nstep),dist,dmin,xx,yy,xt,yt
!-----------------------------------------------------------------------
!
!
  if(iter==1) then
!
!   Position 1
    xt   =-2.5d0 
    yt   =-0.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(1) = ic
          jprb(1) = jc
        endif
!     endif
    enddo
    enddo
!   Position 2
    xt   = 0.0d0 
    yt   =-2.5d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(2) = ic
          jprb(2) = jc
        endif
!     endif
    enddo
    enddo
!   Position 3
    xt   = 2.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(3) = ic
          jprb(3) = jc
        endif
!     endif
    enddo
    enddo
!   Position 4
    xt   =10.5d0 
    yt   = 0.0d0
    dmin = 9.d0
    do jc=1,jmax
    do ic=1,imx1
!     if(ihx(ic,jc)==1) then      
        xx   = xcen(ic)
        yy   = 0.5d0*(ycen(jc)+ycen(jc+1))
        dist = (xx-xt)**2+(yy-yt)**2
        if(dist<=dmin) then
          dmin    = dist
          iprb(4) = ic
          jprb(4) = jc
        endif
!     endif
    enddo
    enddo
    jprb(4) = jprb(3)
!
    open(60,file='probe_ens.dat',form='formatted',status='replace',action='write')
    write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
    write(60,'(" iter   step      ave1        std1        ave2        std2        ave3        std3        ave4        std4")')
    close(60)
    open(60,file='probe_mes.dat',form='formatted',status='replace',action='write')
    write(60,'(" pos(x,y)",8f12.7)') (xcen(iprb(ip)),ycen(jprb(ip)),ip=1,4)
    write(60,'(" iter   step      mes1        mes2        mes3        mes4")')
    close(60)
!
  endif
!
  do istep=1,nstep
    ave(1:4,istep) = 0.d0
    var(1:4,istep) = 0.d0
    std(1:4,istep) = 0.d0
    do iens=1,n_ens
      ave(1,istep) = ave(1,istep)+uens(iprb(1),jprb(1),istep+nstep*(iter-1),iens)
      ave(2,istep) = ave(2,istep)+uens(iprb(2),jprb(2),istep+nstep*(iter-1),iens)
      ave(3,istep) = ave(3,istep)+uens(iprb(3),jprb(3),istep+nstep*(iter-1),iens)
      ave(4,istep) = ave(4,istep)+uens(iprb(4),jprb(4),istep+nstep*(iter-1),iens)
    enddo
    ave(1:4,istep) = ave(1:4,istep)/dble(n_ens)
    do iens=1,n_ens
      var(1,istep) = var(1,istep)+(uens(iprb(1),jprb(1),istep+nstep*(iter-1),iens)-ave(1,istep))**2
      var(2,istep) = var(2,istep)+(uens(iprb(2),jprb(2),istep+nstep*(iter-1),iens)-ave(2,istep))**2
      var(3,istep) = var(3,istep)+(uens(iprb(3),jprb(3),istep+nstep*(iter-1),iens)-ave(3,istep))**2
      var(4,istep) = var(4,istep)+(uens(iprb(4),jprb(4),istep+nstep*(iter-1),iens)-ave(4,istep))**2
    enddo
    var(1:4,istep) = var(1:4,istep)/dble(n_ens)
    std(1:4,istep) = dsqrt(var(1:4,istep))
  enddo
!
  open(60,file='probe_ens.dat',form='formatted',position='append')
  do istep=1,nstep
    write(60,'(2i6,8f12.7)') iter,istep+nstep*(iter-1),(ave(ip,istep),std(ip,istep),ip=1,4)
  enddo
  close(60)
  open(60,file='probe_mes.dat',form='formatted',position='append')
  write(60,'(2i6,4f12.7)') iter,nstep*iter,(uref(iprb(ip),jprb(ip),nstep*iter),ip=1,4)
  close(60)
!
  return
end subroutine sub_probe_ens_all
!-----------------------------------------------------------------------
