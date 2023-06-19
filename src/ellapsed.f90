! Computation ellapsed time in days
! TD - time difference in days

double precision function ellapsed(X)
  implicit none
  double precision, intent(in) :: X
  double precision, parameter :: PC=366.0, HPC = PC/2
  ! Compute ellapsed time difference
  ellapsed = sqrt(X * X)
  ! Correct ellapsed time with year cycle
  if (ellapsed > HPC) then
     ellapsed = PC - ellapsed
  end if
  ellapsed = abs(ellapsed)
end function ellapsed
