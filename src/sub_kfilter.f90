!------------------------------------------------------------------------------
!
!   Ensemble Kalman Filter (EnKF) based on:
!
!   G. Evensen, "The Ensemble Kalman Filter: Theoretical Formulation and
!   Practical Implementation," Ocean Dynamics, Vol.53, pp.343-367, 2003.
!
!   <<< n_mes x n_mes matrix inversion >>>
!
!------------------------------------------------------------------------------
subroutine sub_kfilter_mes(vars,n_mes,n_ens,n_var,ef,cf,pm,cm,rm,ih,x1d,y1d,X5)
  use mod_variables,only : variances
  use m_ranmean3
  use m_random3
  implicit none
!
  type(variances),intent(in)    :: vars                   ! System and mesurement variances
!
  real(8)        ,intent(inout) :: ef(n_var,n_ens)        ! State vector matrix 
                                                          ! -- state vectors for all samples
  real(8)        ,intent(inout) :: cf(n_var)              ! State vector matrix 
  real(8)        ,intent(inout) :: pm(n_mes,n_ens)        ! Pseudo observations
  real(8)        ,intent(inout) :: cm(n_mes)              ! Central case
  real(8)        ,intent(in)    :: rm(n_mes)              ! Real observations
  real(8)        ,intent(in)    :: x1d(n_var)             ! Central state vector
  real(8)        ,intent(in)    :: y1d(n_var)             ! Central state vector
  integer        ,intent(in)    :: ih(n_mes)              ! Real observations
!
  integer                       :: i,j,m,n_mes,n_ens,n_var,ierr
  real(8)                       :: alpha,disn,localop
!
  real(8)                       :: efl(n_var,n_ens)       ! Localized state vector matrix 
  real(8)                       :: HVHR(n_mes,n_mes)      ! (HVH+R)
  real(8)                       :: HVHRi(n_mes,n_mes)     ! (HVH+R)^(-1)
  real(8)                       :: VH(n_var,n_mes)        ! (VH)
  real(8)                       :: w(n_var)               ! Average state
  real(8)                       :: v(n_mes)               ! Average measurement
  real(8)                       :: R(n_mes,n_mes)         ! Measurement error matrix (R)
  real(8)                       :: Le(n_var)              ! Localization function
!
  real(8)                       :: pml(n_mes,n_ens)       ! Work variables
  real(8)                       :: d(n_mes)
  real(8)                       :: WW(n_mes)
  real(8)                       :: UU(n_mes,n_mes)
  real(8)                       :: VV(n_mes,n_mes)
  real(8)                       :: rv1(n_mes)
!
  real(8)                       :: X5(n_ens,n_ens)
!
!
!-Measurement error covariance matrix (R) and measuremnt operator matrix (H)
  R(:,:) = 0.d0
  do m=1,n_mes
    R(m,m) = vars%mes      ! n_mes x n_mes diagonal obs error covariance matrix
  enddo
!
!-Average of sample states (w)
  call ranmean3(ef,w,n_var,n_ens)
  call ranmean3(pm,v,n_mes,n_ens)
!
!-Deviation from the average (ef-w),
!-where covariance matrix (V) is represented by V = (1/(n_ens-1))(ef-w)(ef-w)^t
  do i=1,n_ens
    efl(:,i) = ef(:,i)-w(:)
  enddo
  do i=1,n_ens
    pml(:,i) = pm(:,i)-v(:)
  enddo
!
!-Covariance inflation
  do i=1,n_ens
    efl(:,i) = efl(:,i)*vars%coi 
  enddo
  do i=1,n_ens
    pml(:,i) = pml(:,i)*vars%coi
  enddo
!
!-Covariance localization
  if(vars%col>0.d0) then
    Le(:) = 0.d0
    do i=1,n_mes
      do j=1,n_var
        disn = 2.d0*dsqrt((x1d(j)-x1d(ih(i)))**2+(y1d(j)-y1d(ih(i)))**2)/dble(vars%col)
        if(disn<=1.d0) then
          localop = 1.d0-0.25d0*disn**5+0.5*disn**4+(5.d0/8.d0)*disn**3 &
                   -(5.d0/3.d0)*disn**2
        elseif(1.d0<disn .and. disn<=2.d0) then 
          localop = (1.d0/12.d0)*disn**5-0.5*disn**4+(5.d0/8.d0)*disn**3 &
                   +(5.d0/3.d0)*disn**2-5.d0*disn+4.d0-(2.d0/3.d0)/disn
        else
          localop = 0.d0
        endif
        Le(j) = Le(j)+localop
      enddo
    enddo
    do j=1,n_var
      Le(j) = dmin1(Le(j),1.d0)
    enddo
    do j=1,n_ens
      do i=1,n_var
        efl(i,j) = efl(i,j)*Le(i)
      enddo
    enddo
  endif
!
!-Calculation of (covariance matrix) x (mesurement operator)^T --> (VH)
  VH  = matmul(efl,transpose(pml))
  VH = (1.d0/dble(n_ens-1))*VH
!
!-Calculation of a matrix (HVH^t+R) to be inversed
  HVHR = (1.d0/dble(n_ens-1))*matmul(pml,transpose(pml))+R
!
!-Matrix inversion (HVH^t+R)^(-1)
  call svd(n_mes,n_mes,n_mes,HVHR,WW,.true.,UU,.true.,VV,ierr,rv1)
  do j=1,n_mes  
    do i=1,n_mes  
      VV(i,j) = VV(i,j)/WW(j)
    enddo
  enddo
  HVHRi = matmul(VV,transpose(UU))
!
!-Filtering of central state vector
  d = rm-cm
! (HVH^t+R)^(-1) (y-Hx)
  d = matmul(HVHRi,d) 
! x_new = x + VH^t (HVH^t+R)^(-1) (y-Hx)
  cf = cf+matmul(VH,d)
!
!-Filtering of each ensemble member
  pml = pml/dble(n_ens-1) 
  do j=1,n_ens
!
    call random3(d,n_mes,1)
!   (y-Hx+w)
    do i=1,n_mes
      d(i) = rm(i)+dsqrt(vars%mes)*d(i)-pm(i,j)
    enddo
!   (HVH^t+R)^(-1) (y-Hx+w)
    d = matmul(HVHRi,d) 
!   x_new = x + VH^t (HVH^t+R)^(-1) (y-Hx+w)
    ef(:,j) = ef(:,j)+matmul(VH,d)
!
!   X5(n_ens,n_ens) matrix
    X5(:,j) = matmul(transpose(pml),d)
!
  enddo
!
  do j=1,n_ens
    X5(j,j) = X5(j,j)+1.d0
  enddo
!
  return
end subroutine sub_kfilter_mes
!------------------------------------------------------------------------------
!
!
!
!------------------------------------------------------------------------------
!
!   Ensemble Kalman Filter (EnKF) based on:
!
!   G. Evensen, "The Ensemble Kalman Filter: Theoretical Formulation and
!   Practical Implementation," Ocean Dynamics, Vol.53, pp.343-367, 2003.
!
!   <<< n_ens x n_ens matrix inversion >>>
!
!------------------------------------------------------------------------------
subroutine sub_kfilter_ens(vars,n_mes,n_ens,n_var,ef,cf,pm,cm,rm,ih,x1d,y1d,X5)
  use mod_variables,only : variances
  use m_ranmean3
  use m_random3
  implicit none
!
  type(variances),intent(in)    :: vars                   ! System and mesurement variances
  real(8),        intent(inout) :: ef(n_var,n_ens)        ! State vector matrix 
                                                          ! -- state vectors for all samples
  real(8),        intent(inout) :: cf(n_var)              ! Central state vector
  real(8),        intent(inout) :: pm(n_mes,n_ens)        ! Pseudo observations
  real(8),        intent(inout) :: cm(n_mes)              ! Central pseudo observations
  real(8),        intent(in)    :: x1d(n_var)             ! Central state vector
  real(8),        intent(in)    :: y1d(n_var)             ! Central state vector
  real(8),        intent(in)    :: rm(n_mes)              ! Real observations
  integer,        intent(in)    :: ih(n_mes)              ! Real observations
!
  integer                       :: i,j,ierr,m,n_mes,n_ens,n_var
  real(8)                       :: alpha,disn,localop
!
  real(8)                       :: efl(n_var,n_ens)       ! Localized state vector matrix 
  real(8)                       :: w(n_var)               ! Average state
  real(8)                       :: v(n_mes)               ! Average measurement
  real(8)                       :: R(n_mes,n_mes)         ! Measurement error matrix (R)
  real(8)                       :: Le(n_var)              ! Local
  real(8)                       :: Lm(n_mes,n_mes)        ! Local
  real(8)                       :: af(n_var)              ! Average state vector
!
  real(8)                       :: pml(n_mes,n_ens)       ! Work variables
  real(8)                       :: d (n_mes)
  real(8)                       :: dd(n_ens)
!
!-For (n_ens)x(n_ens) matrix inversion
  real(8)                       :: IHEREH(n_ens,n_ens)    ! (I+(HE)^T R^(-1) HE)
  real(8)                       :: IHEREHi(n_ens,n_ens)   ! (I+(HE)^T R^(-1) HE)^(-1)
  real(8)                       :: Ri(n_mes,n_mes)        ! Inverse of measurement error matrix (R)
  real(8)                       :: II(n_ens,n_ens)        ! n_ens x n_ens identity matrix 
  real(8)                       :: WW(n_ens)
  real(8)                       :: UU(n_ens,n_ens)
  real(8)                       :: VV(n_ens,n_ens)
  real(8)                       :: rv1(n_ens)
!
  real(8)                       :: X5(n_ens,n_ens)
!
!
!-Measurement error covariance matrix (R) and measuremnt operator matrix (H)
  R (:,:) = 0.d0
  Ri(:,:) = 0.d0
  do m=1,n_mes
    R (m,m) = vars%mes     ! n_mes x n_mes diagonal obs error covariance matrix
    Ri(m,m) = 1.d0/R(m,m)  ! Inverse of n_mes x n_mes diagonal obs error covariance matrix
  enddo
!
  II(:,:) = 0.d0
  do m=1,n_ens
    II(m,m) = 1.d0         ! n_ens x n_ens identity matrix
  enddo
!
!-Average of sample states (w)
  call ranmean3(ef,w,n_var,n_ens)
  call ranmean3(pm,v,n_mes,n_ens)
!
!-Deviation from the average efl = (ef-w)/sqrt(n_ens-1),
!-where covariance matrix (V) is represented by V = (efl)(efl)^t
  do i=1,n_ens
    efl(:,i) = ef(:,i)-w(:)
  enddo
  do i=1,n_ens
    pml(:,i) = pm(:,i)-v(:)
  enddo
!
!-Covariance inflation
  do i=1,n_ens
    efl(:,i) = efl(:,i)*vars%coi
  enddo
  do i=1,n_ens
    pml(:,i) = pml(:,i)*vars%coi
  enddo
!
!-Covariance localization
  if(vars%col>0.d0) then
    Le(:) = 0.d0
    do i=1,n_mes
      do j=1,n_var
        disn = 2.d0*dsqrt((x1d(j)-x1d(ih(i)))**2+(y1d(j)-y1d(ih(i)))**2)/dble(vars%col)
        if(disn<=1.d0) then
          localop = 1.d0-0.25d0*disn**5+0.5*disn**4+(5.d0/8.d0)*disn**3 &
                   -(5.d0/3.d0)*disn**2
        elseif(1.d0<disn .and. disn<=2.d0) then 
          localop = (1.d0/12.d0)*disn**5-0.5*disn**4+(5.d0/8.d0)*disn**3 &
                   +(5.d0/3.d0)*disn**2-5.d0*disn+4.d0-(2.d0/3.d0)/disn
        else
          localop = 0.d0
        endif
        Le(j) = Le(j)+localop
      enddo
    enddo
    do j=1,n_var
      Le(j) = dmin1(Le(j),1.d0)
    enddo
    do j=1,n_ens
      do i=1,n_var
        efl(i,j) = efl(i,j)*Le(i)
      enddo
    enddo
  endif
!
!-Calculation of a matrix (I+(HE)^T R^(-1) HE), which needs to be inversed
  IHEREH=II+matmul(matmul(transpose(pml),Ri),pml)/dble(n_ens-1)
!
!-Matrix inversion (HVH^t+R)^(-1)            : n_mes x n_mes (sub_kfilter_mes)
!-Matrix inversion (I+(HE)^t R^(-1) HE)^(-1) : n_ens x n_ens (here)
  call svd(n_ens,n_ens,n_ens,IHEREH,WW,.true.,UU,.true.,VV,ierr,rv1)
  do j=1,n_ens  
    do i=1,n_ens   
      VV(i,j) = VV(i,j)/WW(j)
    enddo
  enddo
  IHEREHi = matmul(VV,transpose(UU))
!
!-Filtering of central state vector
! (y-Hx)
  d = rm-cm 
! (I+(HE)^T R^(-1) HE)^(-1) (HE)^(T) R^(-1) (y-Hx)
  dd = matmul(matmul(IHEREHi,matmul(transpose(pml),Ri)),d)
! x_new = x + E (I+(HE)^T R^(-1) HE)^(-1) (HE)^(T) R^(-1) (y-Hx)
  cf = cf+matmul(efl,dd)/dble(n_ens-1)
!
!-Filtering of each sample
  do j=1,n_ens
!
    call random3(d,n_mes,1)
!   (y-Hx+w)
    d       = rm+dsqrt(vars%mes)*d-pm(:,j)
!   (I+(HE)^T R^(-1) HE)^(-1) (HE)^(T) R^(-1) (y-Hx)
    dd      = matmul(matmul(IHEREHi,matmul(transpose(pml),Ri)),d)
!   x_new = x + E (I+(HE)^T R^(-1) HE)^(-1) (HE)^(T) R^(-1) (y-Hx)
    ef(:,j) = ef(:,j)+matmul(efl,dd)/dble(n_ens-1)
!
!   X5(n_ens,n_ens) matrix
    X5(:,j) = dd(:)/dble(n_ens-1) 
!
  enddo
!
  do j=1,n_ens
    X5(j,j) = X5(j,j)+1.d0
  enddo
!
  return
end subroutine sub_kfilter_ens
!------------------------------------------------------------------------------


