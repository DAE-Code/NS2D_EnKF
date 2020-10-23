!-----------------------------------------------------------------------
subroutine sub_bc_outer(ustg,vstg,pcnt)
!-----------------------------------------------------------------------
  use mod_variables,only : imax,jmax,imx1,jmx1,     &
                           u_inf,v_inf,p_inf,uinf,  &
                           dt,dx   
  implicit none
!
  integer :: ic,jc,kc,is,js,ks
  real(8) :: uave,dtdxi
  real(8) :: ustg(imx1,jmax)
  real(8) :: vstg(imax,jmx1)
  real(8) :: pcnt(imax,jmax)
!-----------------------------------------------------------------------
!
!
! Boundaries in y-direction: Slip boundary
  do ic=3,imax-1
    ustg(ic,1:2        ) = ustg(ic,3     )
    ustg(ic,jmax-1:jmax) = ustg(ic,jmax-2)
  enddo
  do ic=3,imax-2
    vstg(ic,1          ) =-vstg(ic,4     )
    vstg(ic,2          ) =-vstg(ic,3     )
    vstg(ic,jmx1       ) =-vstg(ic,jmx1-3)
    vstg(ic,jmx1-1     ) =-vstg(ic,jmx1-2)
  enddo
  do ic=2,imax-1
    pcnt(ic,1          ) = pcnt(ic,2     )
    pcnt(ic,jmax       ) = pcnt(ic,jmax-1)
  enddo
!
!
! Boundaries in x-direction
  do jc=1,jmax
    ustg(1:2        ,jc) = uinf(jc)
    ustg(imx1-1:imx1,jc) = ustg(imx1-2,jc)
    pcnt(1          ,jc) = pcnt(2     ,jc)
    pcnt(imax       ,jc) = pcnt(imax-1,jc)
  enddo
  do jc=1,jmx1
    vstg(1:2        ,jc) = 0.d0
    vstg(imax-1:imax,jc) = vstg(imax-2,jc)
  enddo
!
  return
end subroutine sub_bc_outer
!-----------------------------------------------------------------------
