module m_ranvar3
contains
subroutine ranvar3(ens,ave,var,n_var,n_ens)
! Return ensemble variance vector
!
  implicit none
  integer             :: j
  integer,intent(in)  :: n_var,n_ens
  real(8),intent(in)  :: ens(n_var,n_ens)
  real(8),intent(in)  :: ave(n_var)
  real(8),intent(out) :: var(n_var)
!
  var(:)=0.d0
  do j=1,n_ens
    var(:)=var(:)+(ens(:,j)-ave(:))*(ens(:,j)-ave(:))
  enddo
  var(:)=var(:)/dble(n_ens-1)
!
end subroutine ranvar3
end module
