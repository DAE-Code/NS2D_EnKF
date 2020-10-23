subroutine sub_measure(vars)
  use mod_variables,only : imax,jmax,imx1,jmx1,n_prb,n_mes,variances,  &
                           xcen,ycen,ih,ihx,ihy,iskip,jskip,mbot,mtop, &
                           mlft,mrht
  implicit none
!
  type(variances) :: vars
  integer :: ic,jc,indx
  integer :: iht(imax,jmax)
!
!
  iht(:,:) = 0
!
  do jc=1,jmax,iskip
  do ic=1,imax,jskip
    if(mlft<xcen(ic) .and. xcen(ic)<mrht .and. mbot<ycen(jc) .and. ycen(jc)<mtop) then
      iht(ic,jc) = 1
    else
      iht(ic,jc) = 0
    endif
  enddo
  enddo
!
  allocate(ihx(imx1,jmax),ihy(imax,jmx1))
  ihx(:,:) = 0
  ihy(:,:) = 0
!
  do jc=1,jmax
  do ic=1,imax
    if(iht(ic,jc)==1) then
      ihx(ic  ,jc) = 1
      ihx(ic+1,jc) = 1
      ihy(ic,jc  ) = 1
      ihy(ic,jc+1) = 1
    endif
  enddo
  enddo
  n_mes = 0
  do jc=1,jmax
  do ic=1,imx1
    if(ihx(ic,jc)==1) then
      n_mes = n_mes+1
    endif
  enddo
  enddo
  do jc=1,jmx1
  do ic=1,imax
    if(ihy(ic,jc)==1) then
      n_mes = n_mes+1
    endif
  enddo
  enddo
!
!
  allocate(ih(n_mes))
  ih(:) = 0
  n_mes = 0
  indx  = 0
  do jc=1,jmax
  do ic=1,imx1
    indx = indx+1
    if(ihx(ic,jc)==1) then
      n_mes = n_mes+1
      ih(n_mes) = indx
    endif
  enddo
  enddo
  do jc=1,jmx1
  do ic=1,imax
    indx = indx+1
    if(ihy(ic,jc)==1) then
      n_mes = n_mes+1
      ih(n_mes) = indx
    endif
  enddo
  enddo
  write(*,*)
  write(*,*) "Number of measurement points:",n_mes
!
  return
end subroutine sub_measure
