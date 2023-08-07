! This is a Fortran function that applies a logistic transformation on an input time difference (TD) and distance (DIST)
! using a pair of weights (TW1, TW2). It is intended to be called from C, and thus uses the iso_c_binding module
! to ensure compatibility.

! Args:
! DIST: A double precision value representing the distance parameter in the logistic transformation.
! TD: A double precision value representing the time difference parameter in the logistic transformation.
! TW1: A double precision value representing the first weight parameter in the logistic transformation.
! TW2: A double precision value representing the second weight parameter in the logistic transformation.

! Returns:
! A double precision value that is the result of applying the logistic transformation on DIST and TD
! using the weights TW1 and TW2.
function logistic_tw(DIST, TD, TW1, TW2) bind(C, name="logistic_tw") result(res)
  use iso_c_binding
  implicit none
  real(c_double), intent(in) :: DIST, TD, TW1, TW2
  real(c_double) :: res
  res = DIST + 1.0d0 / (1.0d0 + exp(-TW1 * (TD - TW2)))
end function logistic_tw






! Computation of TWDTW
! XM - matrix with the time series (N,D)
! YM - matrix with the temporal profile (M,D)
! CM - Output cumulative cost matrix (N+1,M)
! DM - Direction matrix (N+1,M)
! VM - Starting points matrix (N+1,M)
! N  - Number of observations in the time series YM
! M  - Number of observations in the time series XM
! D  - Number of spectral dimensions including time in XM and YM
! TW - Time-Weight parameters alpha and beta
! LB - Maximum elapsed time to constrain TWDTW calculation
! JB - Output array of starting points
! CL - The length of the time cycle
! callback_func - A time-weight fucntion
subroutine twdtwf90(XM, YM, CM, DM, VM, N, M, D, TW, LB, JB, CL, callback_func) bind(C, name="twdtwf90")
  use, intrinsic :: ieee_arithmetic
  use iso_c_binding
  implicit none
  interface
    function callback_func(x, y, z, w) bind(C)
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), intent(in) :: x, y, z, w
      real(c_double) :: callback_func
    end function callback_func
  end interface

  ! I/O Variables
  integer(c_int), intent(in) :: N, M, D
  integer(c_int), intent(out) :: DM(N+1,M), VM(N+1,M), JB(N)
  real(c_double), intent(in) :: XM(M,D), YM(N,D), TW(2), LB, CL
  real(c_double), intent(out) :: CM(N+1,M)
  ! Internals
  real(c_double) :: CP, ST, TD, DIST
  integer(c_int) :: I, J, K, ZERO=0, ONE=1, JM, VM_value
  real(c_double) :: NAN, INF
  NAN = ieee_value(0.0, ieee_quiet_nan)
  INF = ieee_value(0.0, ieee_positive_inf)

  VM(1,1) = 1

  ! Initialize the first row and col of the matrices
  do I = 2, N+1
    TD = abs(YM(I-1,1) - XM(1,1))
    TD = min(TD, CL - TD)
    DIST = 0.0
    do K = 2, D
      DIST = DIST + (YM(I-1,K) - XM(1,K))**2
    end do
    CM(I,1) = CM(I-1,1) + callback_func(sqrt(DIST), TD, TW(1), TW(2))
    DM(I,1) = 3
    VM(I,1) = 1
  end do

  ! Compute cumulative cost matrix
  J = 2
  do while ( J .le. M )
    VM(1,J) = J
    I = 2
    do while ( I .le. N+1 )
      ! Calculate local distance
      ! the call takes I-1 because local matrix has an additional row at the beginning
      TD = abs(YM(I-1,1) - XM(J,1))
      TD = min(TD, CL - TD)
      if (TD.gt.LB) then
        CM(I,J) = INF
        DM(I,J) = -ONE
        VM(I,J) = ZERO
      else
        DIST = 0.0
        do K = 2, D
          DIST = DIST + (YM(I-1,K) - XM(J,K))**2
        end do
        CP = callback_func(sqrt(DIST), TD, TW(1), TW(2))
        CM(I, J) = CP + CM(I - 1, J - 1)
        DM(I,J) = ONE
        VM(I,J) = VM(I - 1, J - 1)

        ST = CP + CM(I    , J - 1)
        if (ST < CM(I, J)) then
              DM(I,J) = 2
              CM(I, J) = ST
              VM(I,J) = VM(I, J - 1)
        endif

        ST = CP + CM(I - 1, J    )
        if (ST < CM(I, J)) then
              DM(I,J) = 3
              CM(I, J) = ST
              VM(I,J) = VM(I - 1, J)
        endif

      end if
      I = I + 1
    end do
    J = J + 1
  end do
  J = 1
  K = ZERO
  do while ( J .le. M )
    VM_value = VM(N+1,J)
    if (VM_value.ne.ZERO) then
      if (K.eq.ZERO) then
        K = 1
        JB(K) = J
        JM = VM_value
      elseif (VM_value.ne.JM) then
        K = K + 1
        JB(K) = J
        JM = VM_value
      elseif (CM(N+1,J).lt.CM(N+1,JB(K))) then
        JB(K) = J
      end if
    end if
    J = J + 1
  end do
end subroutine twdtwf90


! Computation of TWDTW
! XM - matrix with the time series (N,D)
! YM - matrix with the temporal profile (M,D)
! CM - Output cumulative cost matrix (N+1,M)
! DM - Direction matrix (N+1,M)
! VM - Starting points matrix (N+1,M)
! N  - Number of observations in the time series YM
! M  - Number of observations in the time series XM
! D  - Number of spectral dimensions including time in XM and YM
! TW - Time-Weight parameters alpha and beta
! LB - Maximum elapsed time to constrain TWDTW calculation
! JB - Output array of starting points
! CL - The length of the time cycle
subroutine twdtwf90gt(XM, YM, CM, DM, VM, N, M, D, TW, LB, JB, CL)
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, intent(in) :: N, M, D
  integer :: DM(N+1,M), VM(N+1,M), JB(N)
  double precision, intent(in) :: XM(M,D), YM(N,D), TW(2), LB, CL
  double precision :: CM(N+1,M), CP, ST, TD, DIST
  integer :: I, J, K, ZERO, ONE, JM, VM_value
  parameter(ZERO=0, ONE=1)
  double precision :: NAN, INF
  NAN = ieee_value(0.0, ieee_quiet_nan)
  INF = ieee_value(0.0, ieee_positive_inf)
  VM(1,1) = 1

  ! Initialize the first row and column of the matrices
  do 21 I = 2, N+1
     TD = abs(YM(I-1,1) - XM(1,1))
     TD = min(TD, CL - TD)
     DIST = 0.0
     do K = 2, D
       DIST = DIST + (YM(I-1,K) - XM(1,K))**2
     end do
     CM(I,1) = CM(I-1,1) + sqrt(DIST) + 1.0 / (1.0 + exp(-TW(1) * (TD - TW(2))))
     DM(I,1) = 3
     VM(I,1) = 1
21 continue

  ! Compute cumulative cost matrix
  J = 2
  do 32 while (J .le. M)
     VM(1,J) = J
     I = 2
     do 22 while (I .le. N+1)

        TD = abs(YM(I-1,1) - XM(J,1))
        TD = min(TD, CL - TD)

        if (TD.gt.LB) then
           CM(I,J) =  INF
           DM(I,J) = -ONE
           VM(I,J) = ZERO
           goto 44
        endif

        DIST = 0.0
        do K = 2, D
          DIST = DIST + (YM(I-1,K) - XM(J,K))**2
        end do
        CP = sqrt(DIST) + 1.0 / (1.0 + exp(-TW(1) * (TD - TW(2))))
        CM(I, J) = CP + CM(I - 1, J - 1)
        DM(I,J) = ONE
        VM(I,J) = VM(I - 1, J - 1)

        ST = CP + CM(I    , J - 1)
        if (ST < CM(I, J)) then
              DM(I,J) = 2
              CM(I, J) = ST
              VM(I,J) = VM(I, J - 1)
        endif

        ST = CP + CM(I - 1, J    )
        if (ST < CM(I, J)) then
              DM(I,J) = 3
              CM(I, J) = ST
              VM(I,J) = VM(I - 1, J)
        endif

44    continue
        I = I + 1
22    continue
     J = J + 1
32 continue
99 continue

  J = 1
  K = ZERO
  do 69 while (J .le. M)
     VM_value = VM(N+1,J)
     if (VM_value /= ZERO) then
        if (K == ZERO) then
           K = 1
           JB(K) = J
           JM = VM_value
           goto 68
        endif
        if (VM_value /= JM) then
           K = K + 1
           JB(K) = J
           JM = VM_value
           goto 68
        endif

        if (CM(N+1,J) < CM(N+1,JB(K))) then
           JB(K) = J
           goto 68
        endif
     endif
68    continue
     J = J + 1
69 continue

end subroutine twdtwf90gt
