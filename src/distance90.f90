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


real function distancen(YM, XM, N, M, D, I, J, TW, TD)
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, intent(in) :: N, M, D, I, J
  double precision, dimension(:,:), intent(in) :: XM, YM
  double precision, dimension(2), intent(in) :: TW
  double precision, intent(in) :: TD
  double precision :: BD, CD
  integer :: K
  distancen = ieee_value(0.0, ieee_quiet_nan)
  CD = 0.0
  do K = 2, D
     BD = YM(I,K) - XM(J,K)
     CD = CD + (BD * BD)
  end do
  distancen = sqrt(CD) + 1.0 / (1.0 + exp(TW(1) * (TD - TW(2))))
end function distancen
