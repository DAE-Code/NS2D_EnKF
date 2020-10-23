!-----------------------------------------------------------------------
!
!   Ensemble Kalman filter (EnKF) coupled with 2D Navier-Stokes code
!   prepared for the book "Data Assimilation Fluid Science"
!
!   Further information is avairable on: https://github.com/DAE-Code
!
!-----------------------------------------------------------------------
program NS2D_EnKF 
  use mod_variables
  use m_random3
  use m_ranmean3
  use m_ranvar3
  implicit none
!
  type(variances) vars
  integer :: i,j,k,m,n,lc,kc,iens,ii,icheck
  real(8) :: fJ_ens,fJ_cen,fJ_all,fJalll,mvi
  real(8) :: sumenergy
  integer :: iter,ios,reclen
  integer :: ic,jc,ib,indx,ind2
  integer :: nblock = 1
  real(4) :: re4,step4,time4
  real(4) :: time1,tnow
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reading initial parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  open(28,file='EnKF.inp',form='formatted')
  read(28,*) n_prb               ! Karman vortex or vortex advection
  read(28,*) n_kfs               ! 1:EnKF(mes), 2:EnKS(mes), 3:EnKF(ens), 4:EnKS(ens)
  read(28,*) n_flt               ! Number of filtering
  read(28,*) n_ens               ! Ensemble size
  read(28,*) vars%ivr            ! Ensemble generation methods
  read(28,*) vars%ini            ! Initial error variance
  read(28,*) vars%ref            ! Error variance of pseudo measurement
  read(28,*) vars%mes            ! Measurement error variance
  read(28,*) vars%col            ! Covariance localization distance
  read(28,*) vars%coi            ! Covariance inflation factor
  read(28,*) iskip,jskip         ! Every iskip & jskip for measuremnt 
  read(28,*) mlft,mrht           ! Left and right of measurement area
  read(28,*) mbot,mtop           ! Bottom and top of measurement area
  read(28,*)
  read(28,*) nstep               ! Number of time-steps
  read(28,*) Re                  ! Reynolds number
  read(28,*) jmax                ! Number of mesh: jmax (imax x jmax, imax = 3*jmax)
  read(28,*) dt                  ! Time step
  read(28,*) irestart            ! 0:initial, 1:restart
  read(28,*) ioutput             ! 0:No output to screen,1:otherwise
  read(28,*) iskip_plot          ! Output plot3d interval
  read(28,*) iskip_hist          ! Output history interval
  close(28)
  imax = 3*jmax
  imx1 = imax+1
  jmx1 = jmax+1
!
  call system('mkdir output')
  call system('mkdir statis')
  call cpu_time(time1)
!
!
!-----------------------------------------------------------------------
!-Initial setup
!-----------------------------------------------------------------------
  call sub_initial
!
  if(irestart==1) then
    if(ioutput==1) write(*,*) 'Restart'
    call sub_restart
    tstep_ini = tstep
    ostep_ini = ostep
  else
    if(ioutput==1) write(*,*) 'Impulse start'
    tstep = 0.d0
    ostep = 0
  endif
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Define measurement operator
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call sub_measure(vars)
!
! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Allocate arrays for EnKF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  n_var = imx1*jmax+imax*jmx1    ! velocity components (u,v)
  niter = n_flt
  allocate(ef(n_var,n_ens+1))    ! ensemble of state (n_ens+1: central prediction)
  allocate(pm(n_mes,n_ens+1))    ! ensemble of state (n_ens+1: central prediction)
  allocate(rm(n_mes))
  allocate(cf(n_var))
  allocate(rf(n_var))
  allocate(cl(n_var))
  allocate(x1d(n_var))
  allocate(y1d(n_var))
  allocate(ave(n_var),var(n_var),dif(n_var)) 
  allocate(X5(n_ens,n_ens))
!
  allocate(uini(1:imx1,1:jmax), &
           vini(1:imax,1:jmx1), &
           pini(1:imax,1:jmax))
  uini(1:imx1,1:jmax) = ustg(1:imx1,1:jmax)
  vini(1:imax,1:jmx1) = vstg(1:imax,1:jmx1)
  pini(1:imax,1:jmax) = pcnt(1:imax,1:jmax)
  allocate(uref(1:imx1,1:jmax,nstep*niter), &
           vref(1:imax,1:jmx1,nstep*niter), &
           pref(1:imax,1:jmax,nstep*niter))
  uref(1:imx1,1:jmax,:) = 0.d0
  vref(1:imax,1:jmx1,:) = 0.d0
  pref(1:imax,1:jmax,:) = 0.d0
  allocate(uens(1:imx1,1:jmax,nstep,n_ens+1), &
           vens(1:imax,1:jmx1,nstep,n_ens+1), &
           pens(1:imax,1:jmax,nstep,n_ens+1))
  uens(1:imx1,1:jmax,:,:) = 0.d0
  vens(1:imax,1:jmx1,:,:) = 0.d0
  pens(1:imax,1:jmax,:,:) = 0.d0
  allocate(uwns(1:imx1,1:jmx1,niter), &
           vwns(1:imx1,1:jmx1,niter))
  do iter=1,niter
    call random3(uwns(:,:,iter),imx1,jmx1)
    call random3(vwns(:,:,iter),imx1,jmx1)
    uwns(:,:,iter) = dsqrt(vars%ref)*uwns(:,:,iter)
    vwns(:,:,iter) = dsqrt(vars%ref)*vwns(:,:,iter)
  enddo
!-----------------------------------------------------------------------
!-Initial setup end
!-----------------------------------------------------------------------
!
!
  write(*,*) 
  write(*,*) '----------------------------------------------------------'
  write(*,*) 'EnKF(1,3) or EnKS(2,4)    = ',n_kfs
  write(*,*) 'Number of filtering       = ',n_flt
  write(*,*) 'Number of ens members     = ',n_ens
  write(*,*) 'WN(1),TL(2),VTX(3),POD(4) = ',vars%ivr
  write(*,*) 'Initial variance          = ',vars%ini
  write(*,*) 'Variance of pseudo obs    = ',vars%ref
  write(*,*) 'Observation variance      = ',vars%mes
  write(*,*) 'Number of time-step       = ',nstep
  write(*,*) 'Mesh points               = ',imax,jmax
  write(*,*) 'Time-step size            = ',dt
  write(*,*) 'Grid spacing              = ',dx
  write(*,*) 'Reynolds number           = ',Re
  write(*,*) '----------------------------------------------------------'
  write(*,*)
!
  open(87,file='./output/history.txt',form='formatted')
  write(87,'(a8,a5,10a23)') 'istep','itrp','divmax','sumenergy','dt'
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Forward run to generate pseudo measurement data for whole time period
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(n_prb==2) call sub_addvtx(-4.d0,0.d0,1.d0)
  do iter=1,niter
    do istep = 1,nstep
!
      tstep = tstep+dt ! current time
      ostep = ostep+1  ! current time-step
!
      uref(1:imx1,1:jmax,istep+nstep*(iter-1)) = ustg(1:imx1,1:jmax)
      vref(1:imax,1:jmx1,istep+nstep*(iter-1)) = vstg(1:imax,1:jmx1)
      pref(1:imax,1:jmax,istep+nstep*(iter-1)) = pcnt(1:imax,1:jmax)
!
!-----RHS
      call sub_bc_outer(ustg,vstg,pcnt)
      call sub_bc_wall(ustg,vstg,pcnt)
      call sub_rhs
!
!-----Pressure
      call sub_HSMAC(ustg,vstg,pcnt)
!
      call sub_probe(istep,iter,"ref")
!
!-----Standard output to screen
      if(mod(istep,iskip_hist)==0) then
        call sub_energy(sumenergy)
        call sub_cflnum(cfl,dt)
        if(ioutput==1) write(87,'(i8,i5,10e23.15)') istep,itrp,divmax,sumenergy,dt
!
        call cpu_time(tnow)
        if(ioutput==1) then
          write(*,'("----------------------------------------------------------")') 
          write(*,'(" Time step     =",i11,  ",   Time          =",f11.5)') &
                                       ostep,                      tstep
          write(*,'(" Sum energy    =",e11.4,",   Time step     =",f11.5)') &
                                       sumenergy,                  dt
          write(*,'(" Divmax        =",e11.4,",   CFL number    =",e11.4)') &
                                       divmax,                     cfl               
          write(*,'(" Pressure iter =",i11,  ",   Elapsed time  =",e11.4)') &
                                      itrp,       tnow-time1
        endif
      endif
!
!-----Plot3D & Tecplot output
      if(mod(istep+nstep*(iter-1),iskip_plot)==0) then
        call sub_p3dwrite(ostep,tstep,'ref')
      endif
!
    enddo
  enddo
!
!!call sub_restart
  ustg(1:imx1,1:jmax) = uini(1:imx1,1:jmax)
  vstg(1:imax,1:jmx1) = vini(1:imax,1:jmx1)
  pcnt(1:imax,1:jmax) = pini(1:imax,1:jmax)
  tstep = tstep_ini
  ostep = ostep_ini
  if(n_prb==1) then
    do istep = 1,int(7.d0/dt)/2   ! Shift half period of vortex shedding
!
      tstep = tstep+dt ! current time
      ostep = ostep+1  ! current time-step
!
!-----RHS
      call sub_bc_outer(ustg,vstg,pcnt)
      call sub_bc_wall(ustg,vstg,pcnt)
      call sub_rhs
!
!-----Pressure
      call sub_HSMAC(ustg,vstg,pcnt)
!
    enddo
  endif
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Reading initial case (used as a center of an ensemble)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do iens=1,n_ens+1   ! n_ens+1 for central prediction (not used)
    indx = 0
    do jc=1,jmax
    do ic=1,imx1
      indx = indx+1
      ef(indx,iens) = ustg(ic,jc)
      x1d(indx)     = xcen(ic)
      y1d(indx)     = ycen(jc)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      indx = indx+1
      ef(indx,iens) = vstg(ic,jc)
      x1d(indx)     = xcen(ic)
      y1d(indx)     = ycen(jc)
    enddo
    enddo
  enddo
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> Generate inital ensemble
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call sub_ensemble_gen(vars)
!
  do iens=1,n_ens+1
    ind2 = 0
    do jc=1,jmax
    do ic=1,imx1
      ind2 = ind2+1
      uens(ic,jc,nstep,iens) = ef(ind2,iens)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      ind2 = ind2+1
      vens(ic,jc,nstep,iens) = ef(ind2,iens)
    enddo
    enddo
    do jc=1,jmax
    do ic=1,imax
      pens(ic,jc,nstep,iens) = pcnt(ic,jc)
    enddo
    enddo
  enddo
  ustg(1:imx1,1:jmax) =  uens(1:imx1,1:jmax,nstep,n_ens+1)
  vstg(1:imax,1:jmx1) =  vens(1:imax,1:jmx1,nstep,n_ens+1)
  pcnt(1:imax,1:jmax) =  pens(1:imax,1:jmax,nstep,n_ens+1)
  istep = 0
  call sub_p3dwrite(ostep_ini,tstep_ini,'res')
!
!-Open files
  open(10,file='RMSE.dat')
  write(10,'("  Step       RMSE")') 
!
!-Open scratch file to use it as temporary memoery
! EnKS only for u-velocity component to save memory 
  if(n_kfs==2 .or. n_kfs==4) then
    inquire(iolength=reclen) iter,uens     !!,vens,pens
    open(50,form='unformatted',access='direct',recl=reclen,status='scratch')
  endif
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!---> Filtering loop
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  fJalll=0.d0
  do iter=1,niter
!
!---Forward run of each ensemble member
    do iens=1,n_ens+1
!
      ustg(1:imx1,1:jmax) = uens(1:imx1,1:jmax,nstep,iens)
      vstg(1:imax,1:jmx1) = vens(1:imax,1:jmx1,nstep,iens)
      pcnt(1:imax,1:jmax) = pens(1:imax,1:jmax,nstep,iens)
!
      tstep = tstep_ini+nstep*(iter-1)*dt ! current time
      ostep = ostep_ini+nstep*(iter-1)
      do istep = 1,nstep
!
        tstep = tstep+dt ! current time
        ostep = ostep+1  ! current time-step
!
!-------RHS
        call sub_bc_outer(ustg,vstg,pcnt)
        call sub_bc_wall(ustg,vstg,pcnt)
        call sub_rhs
!
!-------Pressure
        call sub_HSMAC(ustg,vstg,pcnt)
!
!       if(iens==n_ens+1) call sub_probe(istep,iter,"est")
!
        uens(1:imx1,1:jmax,istep,iens) = ustg(1:imx1,1:jmax)
        vens(1:imax,1:jmx1,istep,iens) = vstg(1:imax,1:jmx1)
        pens(1:imax,1:jmax,istep,iens) = pcnt(1:imax,1:jmax)
!
!-------Plot3D & Tecplot output
        if(mod(nstep*(iter-1)+istep,iskip_plot)==0 .and. iens==n_ens+1) then
          call sub_p3dwrite(ostep,tstep,'res')
        endif
!
      enddo
!
!-----Ensemble data 
      indx = 0
      ind2 = 0
      do jc=1,jmax
      do ic=1,imx1
        if(ihx(ic,jc)==1) then
          indx = indx+1
          pm(indx,iens) = ustg(ic,jc)
        endif
        ind2 = ind2+1
        ef(ind2,iens) = ustg(ic,jc)
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        if(ihy(ic,jc)==1) then
          indx = indx+1
          pm(indx,iens) = vstg(ic,jc)
        endif
        ind2 = ind2+1
        ef(ind2,iens) = vstg(ic,jc)
      enddo
      enddo
!
    enddo
!
!   Write ensembe
!   EnKS only for u-velocity component to save memory 
    if(n_kfs==2 .or. n_kfs==4) then
      write(50,rec=iter) iter,uens    !!,vens,pens
    endif
!
!---Pesudo measurement data 
    indx = 0
    ind2 = 0
    do jc=1,jmax
    do ic=1,imx1
      if(ihx(ic,jc)==1) then
        indx = indx+1
        rm(indx) = uref(ic,jc,nstep*iter)+uwns(ic,jc,iter)
      endif
      ind2 = ind2+1
      rf(ind2) = uref(ic,jc,nstep*iter)
    enddo
    enddo
    do jc=1,jmx1
    do ic=1,imax
      if(ihy(ic,jc)==1) then
        indx = indx+1
        rm(indx) = vref(ic,jc,nstep*iter)+vwns(ic,jc,iter)
      endif
      ind2 = ind2+1
      rf(ind2) = vref(ic,jc,nstep*iter)
    enddo
    enddo
!
!
!---Analysis before EnKF
    call ranmean3(ef,ave,n_var,n_ens)
    call ranvar3(ef,ave,var,n_var,n_ens)
    indx = 0
    do jc=1,jmax
    do ic=1,imax
      indx = indx+1
      q4(ic,jc,1) = real(ave(indx))          ! u-velocity only
      q4(ic,jc,2) = real(var(indx))          ! real(pm(indx,n_ens+1)-rm(indx))
      q4(ic,jc,3) = real( ef(indx,n_ens+1))
      q4(ic,jc,4) = real( rf(indx))
    enddo
    indx = indx+1  ! for i=imx1
    enddo
    write(cicount,'(i6.6)') nstep*iter
    open(15,file='./statis/stat_'//cicount//'.q',form='unformatted')
    write(15) nblock
    write(15) (imax,jmax,ib=1,nblock)
    do ib=1,nblock
      write(15) re4,step4,re4,time4
      write(15) (((q4(ic,jc,m),ic=1,imax),jc=1,jmax),m=1,4)
    enddo
    close(15)
!
!
!---Evaluate cost function & RMSE--------------------------------------
    call ranmean3(pm,ave,n_mes,n_ens)
    call ranvar3(pm,ave,var,n_mes,n_ens)
    fJ_ens = 0.d0
    fJ_cen = 0.d0
    do i=1,n_mes
      mvi    = 1.d0/vars%mes
      fJ_ens = fJ_ens+0.5d0*(ave(i)       -rm(i))**2*mvi
      fJ_cen = fJ_cen+0.5d0*(pm(i,n_ens+1)-rm(i))**2*mvi
    enddo
    write(* ,'(" Step:",i5,",  Ens_cost:",e14.7,",  Cen_cost:",e14.7)') nstep*iter,fJ_ens,fJ_cen
!
!   Ensemble average for calculating RMSE
    do istep = 1,nstep
      do jc=1,jmax
      do ic=1,imx1
        uens(ic,jc,istep,n_ens+1) = sum(uens(ic,jc,istep,1:n_ens))/dble(n_ens)
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        vens(ic,jc,istep,n_ens+1) = sum(vens(ic,jc,istep,1:n_ens))/dble(n_ens)
      enddo
      enddo
    enddo
!
    do istep = 1,nstep
      fJ_all = 0.d0
      do jc=1,jmax
      do ic=1,imx1
        fJ_all = fJ_all+(uens(ic,jc,istep,n_ens+1)-uref(ic,jc,istep+nstep*(iter-1)))**2
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        fJ_all = fJ_all+(vens(ic,jc,istep,n_ens+1)-vref(ic,jc,istep+nstep*(iter-1)))**2
      enddo
      enddo
      write(10,'(i6,e15.7)') istep+nstep*(iter-1),dsqrt(fJ_all/dble(imx1*jmax+imax*jmx1))
      fJalll = fJalll+fJ_all
    enddo
!
!
!---Kalman filter analysis---------------------------------------------
    if(n_kfs<=2) then
      call sub_kfilter_mes(vars,n_mes,n_ens,n_var,ef,ef(:,n_ens+1),pm,pm(:,n_ens+1),rm,ih,x1d,y1d,X5)
    else
      call sub_kfilter_ens(vars,n_mes,n_ens,n_var,ef,ef(:,n_ens+1),pm,pm(:,n_ens+1),rm,ih,x1d,y1d,X5)
    endif
!
!
!   EnKS only for u-velocity component to save memory 
    if(n_kfs==2 .or. n_kfs==4) then
      do ii=1,iter
        read(50,rec=ii) icheck,uens   !!,vens,pens
        do i=1,nstep
          do iens=1,n_ens
            ustg(1:imx1,1:jmax) = 0.d0
     !!     vstg(1:imax,1:jmx1) = 0.d0
            do jc=1,jmax
            do ic=1,imx1
              do j=1,n_ens
                ustg(ic,jc) = ustg(ic,jc)+uens(ic,jc,i,j)*X5(j,iens)
              enddo
            enddo
            enddo
     !!     do jc=1,jmx1
     !!     do ic=1,imax
     !!       do j=1,n_ens
     !!         vstg(ic,jc) = vstg(ic,jc)+vens(ic,jc,i,j)*X5(j,iens)
     !!       enddo
     !!     enddo
     !!     enddo
            uens(1:imx1,1:jmax,i,iens) = ustg(1:imx1,1:jmax)
     !!     vens(1:imax,1:jmx1,i,iens) = vstg(1:imax,1:jmx1)
          enddo
        enddo
        write(50,rec=ii) icheck,uens  !!,vens,pens
      enddo
    else
      call sub_probe_ens(iter)
    endif
!
!
!---Update ensemble
    do iens=1,n_ens+1
      ind2 = 0
      do jc=1,jmax
      do ic=1,imx1
        ind2 = ind2+1
        uens(ic,jc,nstep,iens) = ef(ind2,iens)
      enddo
      enddo
      do jc=1,jmx1
      do ic=1,imax
        ind2 = ind2+1
        vens(ic,jc,nstep,iens) = ef(ind2,iens)
      enddo
      enddo
    enddo
!
  enddo
999 continue
!
! 
! EnKS only for u-velocity component to save memory 
  if(n_kfs==2 .or. n_kfs==4) then
    do iter=1,niter
      read(50,rec=iter) icheck,uens   !!,vens,pens
      call sub_probe_ens(iter)
    enddo
    close(50)
  endif
!
  close(10)
  fJalll = dsqrt(fJalll/dble(nstep*niter*(imx1*jmax+imax*jmx1)))
  open(10,file='RMSE_all.dat')
  write(10,'("RMSE:",e15.7)') fJalll
  write(* ,'("RMSE:",e15.7)') fJalll
  close(10)
!
  deallocate(ave,var,dif) 
  deallocate(ef,cf,rf,pm,rm)
!
  call cpu_time(tnow)
  write(*,'("Program end",f15.3," [sec]")') tnow-time1
!
end program NS2D_EnKF
!-----------------------------------------------------------------------
