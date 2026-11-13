! Time weights for the batch computation of TWDTW, see twdtwknn
!
! Computes the logistic time weight between each observation of the temporal
! patterns and each date of the time series. As the time series share the same
! dates, the weights are computed only once and reused for all time series.
! Weights for an elapsed time greater than LB are set to infinity, use an
! infinite LB to compute the weights without this constraint.
!
! TX - Dates of the time series in units of the time cycle (NT)
! NT - Number of dates of the time series
! TY - Dates of the observations of the temporal patterns in units of the time cycle (NY)
! NY - Number of observations of all temporal patterns
! TW - Time-Weight parameters alpha and beta
! LB - Maximum elapsed time to constrain TWDTW calculation
! CL - The length of the time cycle
! W  - Output time weights (NT,NY)
subroutine twdtwknn_weights(TX, NT, TY, NY, TW, LB, CL, W) bind(C, name="twdtwknn_weights")
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_c_binding
  implicit none

  ! I/O Variables
  integer(c_int), intent(in) :: NT, NY
  real(c_double), intent(in) :: TX(NT), TY(NY), TW(2), LB, CL
  real(c_double), intent(out) :: W(NT,NY)
  ! Internals
  real(c_double) :: TD, INF
  integer(c_int) :: I, J
  INF = ieee_value(1.0d0, ieee_positive_inf)

  do I = 1, NY
    do J = 1, NT
      TD = abs(TY(I) - TX(J))
      TD = min(TD, CL - TD)
      if (TD.gt.LB) then
        W(J,I) = INF
      else
        W(J,I) = 1.0d0 / (1.0d0 + exp(-TW(1) * (TD - TW(2))))
      end if
    end do
  end do
end subroutine twdtwknn_weights


! Batch computation of the TWDTW k-nearest neighbours
!
! Computes the TWDTW distance between each time series in XM and each temporal
! pattern in YM, and keeps for each time series the K patterns with the lowest
! distance. The distance is the same as the lowest distance computed by twdtwf90.
!
! All time series have the same dates. Observations of a time series with a
! missing value (NaN) in any dimension are skipped. The temporal patterns are
! stacked by rows in YM and can have different number of observations.
!
! Only the cumulative cost of one row of the TWDTW matrix is stored, i.e. the
! matrices CM, DM, and VM of twdtwf90 are not available. The lowest cumulative
! cost of a row never decreases from one row to the next. Therefore, the
! computation for a pattern stops as soon as the lowest cost of a row is not
! lower than the distance to the current K-th nearest pattern.
!
! XM   - Array with the time series (NS,NT,ND), missing values are NaN
! NS   - Number of time series
! NT   - Number of dates of the time series
! ND   - Number of spectral dimensions in XM and YM, time not included
! P1   - First time series in XM to be processed
! P2   - Last time series in XM to be processed
! YM   - Matrix with the temporal patterns stacked by rows (NY,ND), without missing values
! YO   - Row offset of each temporal pattern in YM (NP+1), pattern I is in rows YO(I)+1 to YO(I+1)
! NY   - Number of observations of all temporal patterns
! NP   - Number of temporal patterns
! WU   - Time weights without maximum elapsed time constraint (NT,NY), see twdtwknn_weights
! WL   - Time weights constrained by the maximum elapsed time (NT,NY), see twdtwknn_weights
! K    - Number of nearest neighbours
! KI   - Output index of the K nearest patterns (NS,K), not modified for time series without observations
! KD   - Output distance to the K nearest patterns (NS,K), not modified for time series without observations
! IERR - Output error flag, zero if the work arrays were allocated
subroutine twdtwknn(XM, NS, NT, ND, P1, P2, YM, YO, NY, NP, WU, WL, K, KI, KD, IERR) bind(C, name="twdtwknn")
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_c_binding
  implicit none

  ! I/O Variables
  integer(c_int), intent(in) :: NS, NT, ND, P1, P2, NY, NP, K
  integer(c_int), intent(in) :: YO(NP+1)
  real(c_double), intent(in) :: XM(NS,NT,ND), YM(NY,ND), WU(NT,NY), WL(NT,NY)
  integer(c_int), intent(inout) :: KI(NS,K)
  real(c_double), intent(inout) :: KD(NS,K)
  integer(c_int), intent(out) :: IERR
  ! Internals
  ! XB - Observations of a time series without missing values (M,ND)
  ! JX - Index of the dates of the observations in XB (M)
  ! DT - Local cost between an observation of a pattern and the time series (M)
  ! CR - Cumulative cost of a row of the TWDTW matrix (M)
  ! BD - Distance to the nearest patterns of a time series in ascending order (K)
  ! BI - Index of the nearest patterns of a time series (K)
  ! CP, CMIN   - Cost of a cell and lowest cost of a row
  ! DIAG, LEFT, UP - Cumulative cost of the neighbouring cells
  real(c_double), allocatable :: XB(:,:), DT(:), CR(:), BD(:)
  integer(c_int), allocatable :: JX(:), BI(:)
  real(c_double) :: XV, YV, CP, CMIN, DIAG, LEFT, UP
  integer(c_int) :: P, I, J, L, M, IP, NB, POS
  logical :: COMPLETE, ABANDON

  allocate(XB(NT,ND), DT(NT), CR(NT), BD(K), JX(NT), BI(K), stat=IERR)
  if (IERR.ne.0) return

  do P = P1, P2

    ! Gather the observations of the time series without missing values
    M = 0
    do J = 1, NT
      COMPLETE = .true.
      do L = 1, ND
        XV = XM(P,J,L)
        if (ieee_is_nan(XV)) COMPLETE = .false.
        XB(M+1,L) = XV
      end do
      if (COMPLETE) then
        M = M + 1
        JX(M) = J
      end if
    end do

    ! Time series without observations keep the initial values of KI and KD
    if (M.eq.0) cycle

    NB = 0
    do IP = 1, NP

      ! Patterns without observations cannot be a neighbour
      if (YO(IP+1).le.YO(IP)) cycle

      ! Row zero of the TWDTW matrix, a match can start at any observation of the time series
      do J = 1, M
        CR(J) = 0.0d0
      end do
      CMIN = 0.0d0

      ABANDON = .false.
      do I = YO(IP) + 1, YO(IP+1)

        ! Calculate local distance between observation I of the pattern and the time series
        do J = 1, M
          DT(J) = 0.0d0
        end do
        do L = 1, ND
          YV = YM(I,L)
          do J = 1, M
            DT(J) = DT(J) + (YV - XB(J,L))**2
          end do
        end do
        do J = 1, M
          DT(J) = sqrt(DT(J))
        end do

        ! Compute cumulative cost of the row, CR has the cost of the previous row
        ! The first column has only the cost from the previous row and is not
        ! constrained by the maximum elapsed time, consistent with twdtwf90
        DIAG = CR(1)
        LEFT = DT(1) + WU(JX(1),I) + DIAG
        CR(1) = LEFT
        CMIN = LEFT
        do J = 2, M
          UP = CR(J)
          CP = DT(J) + WL(JX(J),I)
          LEFT = CP + min(DIAG, LEFT, UP)
          DIAG = UP
          CR(J) = LEFT
          CMIN = min(CMIN, LEFT)
        end do

        ! Stop if the pattern cannot be among the K nearest
        if (NB.eq.K) then
          if (CMIN.ge.BD(K)) then
            ABANDON = .true.
            exit
          end if
        end if

      end do
      if (ABANDON) cycle

      ! CMIN is the lowest cost of the last row, i.e. the TWDTW distance
      ! Insert the pattern in the list of nearest patterns sorted by distance,
      ! patterns with the same distance are kept in the order of their index
      if (NB.lt.K) NB = NB + 1
      POS = NB
      do while (POS.gt.1)
        if (BD(POS-1).le.CMIN) exit
        BD(POS) = BD(POS-1)
        BI(POS) = BI(POS-1)
        POS = POS - 1
      end do
      BD(POS) = CMIN
      BI(POS) = IP

    end do

    do J = 1, NB
      KI(P,J) = BI(J)
      KD(P,J) = BD(J)
    end do

  end do

  deallocate(XB, DT, CR, BD, JX, BI)
end subroutine twdtwknn


! Local distance for the batch computation of TWDTW with a time-weight function
!
! Computes the Euclidean distance between each observation of a block of time
! series and each observation of one temporal pattern. This is the distance
! that a time-weight function receives, see twdtwknn_cost. Observations of a
! time series with a missing value (NaN) in any dimension are flagged in OK
! and have distance zero.
!
! XM - Array with the time series (NS,NT,ND), missing values are NaN
! NS - Number of time series
! NT - Number of dates of the time series
! ND - Number of spectral dimensions in XM and YM, time not included
! P1 - First time series in XM to be processed
! P2 - Last time series in XM to be processed
! YM - Matrix with the temporal patterns stacked by rows (NY,ND), without missing values
! NY - Number of observations of all temporal patterns
! Y1 - First row of the temporal pattern in YM
! Y2 - Last row of the temporal pattern in YM
! DM - Output local distance (P2-P1+1,NT,Y2-Y1+1)
! OK - Output flag of the observations of the time series (P2-P1+1,NT), one if complete and zero otherwise
subroutine twdtwknn_dist(XM, NS, NT, ND, P1, P2, YM, NY, Y1, Y2, DM, OK) bind(C, name="twdtwknn_dist")
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_c_binding
  implicit none

  ! I/O Variables
  integer(c_int), intent(in) :: NS, NT, ND, P1, P2, NY, Y1, Y2
  real(c_double), intent(in) :: XM(NS,NT,ND), YM(NY,ND)
  real(c_double), intent(out) :: DM(P2-P1+1,NT,Y2-Y1+1)
  integer(c_int), intent(out) :: OK(P2-P1+1,NT)
  ! Internals
  real(c_double) :: DIST
  integer(c_int) :: P, I, J, L

  ! Flag the observations of the time series without missing values
  do J = 1, NT
    do P = P1, P2
      OK(P-P1+1,J) = 1
      do L = 1, ND
        if (ieee_is_nan(XM(P,J,L))) OK(P-P1+1,J) = 0
      end do
    end do
  end do

  ! Calculate local distance as in twdtwf90
  do I = Y1, Y2
    do J = 1, NT
      do P = P1, P2
        if (OK(P-P1+1,J).eq.1) then
          DIST = 0.0d0
          do L = 1, ND
            DIST = DIST + (YM(I,L) - XM(P,J,L))**2
          end do
          DM(P-P1+1,J,I-Y1+1) = sqrt(DIST)
        else
          DM(P-P1+1,J,I-Y1+1) = 0.0d0
        end if
      end do
    end do
  end do
end subroutine twdtwknn_dist


! Batch computation of TWDTW for one temporal pattern given the local cost
!
! Computes the TWDTW distance between each time series of a block and one
! temporal pattern, and updates the K patterns with the lowest distance of
! each time series. CM has the local cost, i.e. the value of a time-weight
! function for the distance and the elapsed time between each observation of
! the time series and each observation of the pattern. Call this subroutine
! once for each pattern, in the order of the patterns.
!
! The cumulative cost is computed as in twdtwknn. However, the computation
! never stops before the last observation of the pattern, because the cost
! given by a time-weight function can be negative.
!
! CM   - Local cost (NB,NT,NI)
! NB   - Number of time series in the block
! NT   - Number of dates of the time series
! NI   - Number of observations of the temporal pattern
! OK   - Flag of the observations of the time series (NB,NT), one if complete, see twdtwknn_dist
! EL   - Elapsed time between each date of the time series and each observation of the pattern (NT,NI)
! LB   - Maximum elapsed time to constrain TWDTW calculation
! IP   - Index of the temporal pattern
! K    - Number of nearest neighbours
! NS   - Number of time series, i.e. the number of rows of KN, KI, and KD
! P1   - Row of the first time series of the block in KN, KI, and KD
! KN   - Number of nearest patterns found for each time series (NS), updated
! KI   - Index of the K nearest patterns (NS,K), updated
! KD   - Distance to the K nearest patterns in ascending order (NS,K), updated
! IERR - Output error flag, zero if the work arrays were allocated
subroutine twdtwknn_cost(CM, NB, NT, NI, OK, EL, LB, IP, K, NS, P1, KN, KI, KD, IERR) bind(C, name="twdtwknn_cost")
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_c_binding
  implicit none

  ! I/O Variables
  integer(c_int), intent(in) :: NB, NT, NI, IP, K, NS, P1
  integer(c_int), intent(in) :: OK(NB,NT)
  real(c_double), intent(in) :: CM(NB,NT,NI), EL(NT,NI), LB
  integer(c_int), intent(inout) :: KN(NS), KI(NS,K)
  real(c_double), intent(inout) :: KD(NS,K)
  integer(c_int), intent(out) :: IERR
  ! Internals
  ! JX - Index of the dates of the complete observations of a time series (M)
  ! CR - Cumulative cost of a row of the TWDTW matrix (M)
  ! CMIN - Lowest cost of a row
  ! DIAG, LEFT, UP - Cumulative cost of the neighbouring cells
  real(c_double), allocatable :: CR(:)
  integer(c_int), allocatable :: JX(:)
  real(c_double) :: CMIN, DIAG, LEFT, UP, INF
  integer(c_int) :: P, R, I, J, M, NN, POS

  allocate(CR(NT), JX(NT), stat=IERR)
  if (IERR.ne.0) return

  INF = ieee_value(1.0d0, ieee_positive_inf)

  do P = 1, NB

    ! Gather the complete observations of the time series
    M = 0
    do J = 1, NT
      if (OK(P,J).eq.1) then
        M = M + 1
        JX(M) = J
      end if
    end do

    ! Time series without observations have no neighbours
    if (M.eq.0) cycle

    ! Row zero of the TWDTW matrix, a match can start at any observation of the time series
    do J = 1, M
      CR(J) = 0.0d0
    end do
    CMIN = 0.0d0

    do I = 1, NI

      ! Compute cumulative cost of the row, CR has the cost of the previous row
      ! The first column has only the cost from the previous row and is not
      ! constrained by the maximum elapsed time, consistent with twdtwf90
      DIAG = CR(1)
      LEFT = CM(P,JX(1),I) + DIAG
      CR(1) = LEFT
      CMIN = LEFT
      do J = 2, M
        UP = CR(J)
        if (EL(JX(J),I).gt.LB) then
          LEFT = INF
        else
          LEFT = CM(P,JX(J),I) + min(DIAG, LEFT, UP)
        end if
        DIAG = UP
        CR(J) = LEFT
        CMIN = min(CMIN, LEFT)
      end do

    end do

    ! CMIN is the lowest cost of the last row, i.e. the TWDTW distance
    ! A distance that is not a number cannot be compared
    if (ieee_is_nan(CMIN)) cycle

    ! Insert the pattern in the list of nearest patterns sorted by distance,
    ! patterns with the same distance are kept in the order of their index
    R = P1 + P - 1
    NN = KN(R)
    if (NN.eq.K) then
      if (CMIN.ge.KD(R,K)) cycle
    else
      NN = NN + 1
    end if
    POS = NN
    do while (POS.gt.1)
      if (KD(R,POS-1).le.CMIN) exit
      KD(R,POS) = KD(R,POS-1)
      KI(R,POS) = KI(R,POS-1)
      POS = POS - 1
    end do
    KD(R,POS) = CMIN
    KI(R,POS) = IP
    KN(R) = NN

  end do

  deallocate(CR, JX)
end subroutine twdtwknn_cost
