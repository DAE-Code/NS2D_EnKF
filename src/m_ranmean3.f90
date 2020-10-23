module m_ranmean3
contains
subroutine ranmean3(ens,ave,n_var,n_ens)
! Return ensemble mean vector
!
  implicit none 
  integer             :: j
  integer,intent(in)  :: n_var,n_ens
  real(8),intent(in)  :: ens(n_var,n_ens)
  real(8),intent(out) :: ave(n_var)
!
  ave(:)=ens(:,1)
  do j=2,n_ens
    ave(:)=ave(:)+ens(:,j)
  enddo
  ave(:)=ave(:)/dble(n_ens)
!
end subroutine ranmean3
end module
