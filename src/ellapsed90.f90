! Computation ellapsed time in days
! TD - time difference in days


subroutine ellapsed90(TD)
  implicit none
  double precision, intent(inout) :: TD
  double precision :: HPC
  double precision, parameter :: PC=366.0
  HPC = PC/2
  ! Compute ellapsed time difference
  TD = sqrt(TD * TD)
  ! Correct ellapsed time with year cycle
  if (TD > HPC) then
     TD = PC - TD
  end if
  TD = abs(TD)
end subroutine ellapsed90
