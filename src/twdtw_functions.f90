! Compute TWDTW distance using logistic weight
! XM - matrix with the time series (N,D)
! YM - matrix with the temporal profile (M,D)
! N  - Number of rows in CM, DM, and VM - time series
! M  - Number of columns CM, DM, and VM - temporal profile
! D  - Number of spectral dimensions including time in XM and YM
! I  - Single point in the time series to calculate the local distance
! J  - Single point in the temporal profile to calculate the local distance
! A  - Time-Weight parameter alpha
! B  - Time-Weight parameter beta


double precision function distance(YM, XM, N, M, D, I, J, TW, TD)
  implicit none
  integer, intent(in) :: N, M, D, I, J
  double precision, intent(in) :: XM(M,D), YM(N,D), TW(2), TD
  double precision :: BD, CD
  integer :: K
  double precision :: dist

  CD = 0.0
  do K = 2, D
     BD = YM(I,K) - XM(J,K)
     CD = CD + (BD * BD)
  end do

  dist = sqrt(CD)
  dist = dist + 1.0 / (1.0 + exp(TW(1) * (TD - TW(2))))

  distance = dist
end function distance


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
