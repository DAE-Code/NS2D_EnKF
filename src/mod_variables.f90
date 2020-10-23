module mod_variables
!
  integer             :: n_prb              ! problem (1:Karman vortex, 2:vortex advection)
  integer             :: n_kfs              ! EnKF(1) or EnKS(2)
  integer             :: n_ens              ! size of ensemble 
  integer             :: n_var              ! size of state vector
  integer             :: n_mes              ! size of measurement vector
  integer             :: n_flt              ! number of filtering
!
  type variances
    integer           :: ivr
    real(8)           :: ini
    real(8)           :: ref                ! Error variance of pseudo measurement
    real(8)           :: mes                ! Measurement error variance
    real(8)           :: col
    real(8)           :: coi
  end type variances
!
  integer             :: iskip              ! Every iskip-th mesh point as measuremnt 
  integer             :: jskip              ! Every jskip-th mesh point as measuremnt 
  real(8)             :: mbot               ! Bottom of measurement area
  real(8)             :: mtop               ! Top of measurement area
  real(8)             :: mlft               ! Left corner of measurement area
  real(8)             :: mrht               ! Right corner of measurement area
!
  real(8),allocatable :: ef(:,:)            ! Ensemble forecast (prediction)
  real(8),allocatable :: pm(:,:)            ! Central forecast (prediction)
  real(8),allocatable :: rm(:)              ! Central forecast (prediction)
  real(8),allocatable :: cf(:)              ! Central forecast (prediction)
  real(8),allocatable :: rf(:)              ! Reference (prediction)
  real(8),allocatable :: cl(:)              ! Covariance localization
  integer,allocatable :: ih (:)             ! Index for measurement points (1D to 1D)
  integer,allocatable :: ihx(:,:)           ! Index for measurement points (2D mesh)
  integer,allocatable :: ihy(:,:)           ! Index for measurement points (2D mesh)
  real(8),allocatable :: x1d(:)             ! Covariance localization
  real(8),allocatable :: y1d(:)             ! Covariance localization
  real(8),allocatable :: X5(:,:)            ! X5 matrix of Evensen(2003)
!
  real(8),allocatable :: uini(:,:)
  real(8),allocatable :: vini(:,:)
  real(8),allocatable :: pini(:,:)
  real(8),allocatable :: uref(:,:,:)
  real(8),allocatable :: vref(:,:,:)
  real(8),allocatable :: pref(:,:,:)
  real(8),allocatable :: uwns(:,:,:)
  real(8),allocatable :: vwns(:,:,:)
  real(4),allocatable :: uens(:,:,:,:)
  real(4),allocatable :: vens(:,:,:,:)
  real(4),allocatable :: pens(:,:,:,:)
!
  real(8),allocatable :: ave(:)             ! Average of predicted state
  real(8),allocatable :: var(:)             ! Variance of predicted state
  real(8),allocatable :: dif(:)             ! Difference of predicted state from 'cf'
!
! Variables for NS2D
  real(8),parameter   :: u_inf = 1.d0
  real(8),parameter   :: v_inf = 0.d0
  real(8),parameter   :: w_inf = 0.d0
  real(8),parameter   :: p_inf = 0.d0
!
  real(8)             :: dt
  real(8)             :: dx
  real(8)             :: re
! 
  integer             :: itrp               ! # of iteration
  real(8)             :: divmax
!
  real(8)             :: cfl
  integer             :: nstep              ! # of sub time-steps
  integer             :: niter              ! # of Kalman filtering
  real(8)             :: tstep              ! current time
  real(8)             :: tstep_ini          ! initial time
  integer             :: istep              ! time step
  integer             :: ostep              ! current step
  integer             :: ostep_ini          ! initial step
  integer             :: irestart           ! restart flag
  integer             :: ioutput            ! output flag
  integer             :: iskip_plot
  integer             :: iskip_hist
!
  integer             :: iobj               ! object size
  integer             :: jobj               !
  integer             :: imax               ! # of mesh points
  integer             :: jmax               !
  integer             :: imx1               !
  integer             :: jmx1               !
  integer             :: corner(2,4)
  integer             :: iprb(4)
  integer             :: jprb(4)
!
  character(len=6)    :: cicount = '000000'
!
  real(8),allocatable :: xcen(:)
  real(8),allocatable :: ycen(:)
!
  real(8),allocatable :: ustg(:,:)
  real(8),allocatable :: vstg(:,:)
  real(8),allocatable :: udlt(:,:)
  real(8),allocatable :: vdlt(:,:)
  real(8),allocatable :: pcnt(:,:)
  real(8),allocatable :: divg(:,:)
!
  integer,allocatable :: iblk(:,:)
  real(8),allocatable :: uinf(:)
  real(8),allocatable :: prob(:,:)
  real(4),allocatable :: q4(:,:,:)
!
end module mod_variables
