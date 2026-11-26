! Time weights for the batch computation of TWDTW, see twdtwnearest
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
subroutine twdtwnearest_weights(TX, NT, TY, NY, TW, LB, CL, W) bind(C, name="twdtwnearest_weights")
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
end subroutine twdtwnearest_weights


! Near dates of the observations of the temporal patterns, see twdtwnearest
!
! For each observation of the temporal patterns, lists the dates of the time
! series with a low time weight, i.e. the dates near the observation in time,
! and computes the lowest time weight of the other dates. A date is near if
! its weight is not greater than WMIN + TAU * (WMAX - WMIN), where WMIN and
! WMAX are the lowest and the highest finite weight of the observation.
!
! The local cost of an observation of a pattern and a date of the time series
! is the distance of their values, which is never negative, plus the time
! weight. Therefore, the lowest cost of the near dates, or the lowest weight
! of the other dates, is a lower bound of the cost of the observation in any
! alignment. This holds for any choice of the near dates.
!
! Call this subroutine with NM equal to zero to get the number of near dates
! in NO(NY+1), and then with NEAR of that size.
!
! WL   - Time weights constrained by the maximum elapsed time (NT,NY), see twdtwnearest_weights
! NT   - Number of dates of the time series
! NY   - Number of observations of all temporal patterns
! NM   - Size of NEAR, near dates beyond this size are counted and not stored
! NO   - Output offset of each observation in NEAR (NY+1), the near dates of observation I are NEAR(NO(I)+1) to NEAR(NO(I+1))
! NEAR - Output index of the near dates of the observations (NM)
! FW   - Output lowest time weight of the dates of each observation that are not near (NY), infinity if there is none
subroutine twdtwnearest_near(WL, NT, NY, NM, NO, NEAR, FW) bind(C, name="twdtwnearest_near")
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_c_binding
  implicit none

  ! I/O Variables
  integer(c_int), intent(in) :: NT, NY, NM
  real(c_double), intent(in) :: WL(NT,NY)
  integer(c_int), intent(out) :: NO(NY+1)
  integer(c_int), intent(inout) :: NEAR(*)
  real(c_double), intent(out) :: FW(NY)
  ! Internals
  real(c_double), parameter :: TAU = 0.15d0
  real(c_double) :: WMIN, WMAX, THR, INF
  integer(c_int) :: I, J, N
  INF = ieee_value(1.0d0, ieee_positive_inf)

  N = 0
  NO(1) = 0
  do I = 1, NY
    WMIN = INF
    WMAX = -INF
    do J = 1, NT
      if (WL(J,I).lt.INF) then
        WMIN = min(WMIN, WL(J,I))
        WMAX = max(WMAX, WL(J,I))
      end if
    end do
    FW(I) = INF
    if (WMIN.lt.INF) then
      THR = WMIN + TAU * (WMAX - WMIN)
      do J = 1, NT
        if (WL(J,I).le.THR) then
          N = N + 1
          if (N.le.NM) NEAR(N) = J
        else
          FW(I) = min(FW(I), WL(J,I))
        end if
      end do
    end if
    NO(I+1) = N
  end do
end subroutine twdtwnearest_near


! Batch computation of the TWDTW k-nearest neighbours
!
! Computes for each time series in XM the K temporal patterns in YM with the
! lowest TWDTW distance. The distance is the same as the lowest distance
! computed by twdtwf90. Patterns with the same distance are in the order of
! their index.
!
! All time series have the same dates. Observations of a time series with a
! missing value (NaN) in any dimension are skipped. The temporal patterns are
! stacked by rows in YM and can have different number of observations.
!
! Only the cumulative cost of one row of the TWDTW matrix is stored, i.e. the
! matrices CM, DM, and VM of twdtwf90 are not available. The lowest cumulative
! cost of a row never decreases from one row to the next. Therefore, the
! computation for a pattern stops as soon as the lowest cost of a row shows
! that the pattern cannot be among the K nearest.
!
! With USELB equal to one, most patterns are skipped without computing their
! TWDTW matrix. The sum of the lower bounds of the cost of the observations of
! a pattern, see twdtwnearest_near, is a lower bound of its TWDTW distance. A
! pattern with a bound greater than the distance to the current K-th nearest
! pattern cannot be among the K nearest. The patterns are processed in three steps:
!   1. The nearest patterns of the previous time series, which are often also
!      near the current time series, e.g. for neighbouring pixels of an image.
!   2. The bound of each other pattern, added row by row until it is greater
!      than the distance to the K-th nearest pattern. A pattern with a bound
!      lower than the bound of all patterns before it is processed at once.
!   3. The remaining patterns in ascending order of their bound, until the
!      bound is greater than the distance to the K-th nearest pattern.
! The result does not depend on the order in which the patterns are processed,
! and it is the same with USELB equal to zero.
!
! XM    - Array with the time series (NS,NT,ND), missing values are NaN
! NS    - Number of time series
! NT    - Number of dates of the time series
! ND    - Number of spectral dimensions in XM and YM, time not included
! P1    - First time series in XM to be processed
! P2    - Last time series in XM to be processed
! YM    - Matrix with the temporal patterns stacked by rows (NY,ND), without missing values
! YO    - Row offset of each temporal pattern in YM (NP+1), pattern I is in rows YO(I)+1 to YO(I+1)
! NY    - Number of observations of all temporal patterns
! NP    - Number of temporal patterns
! WL    - Time weights constrained by the maximum elapsed time (NT,NY), see twdtwnearest_weights
! USELB - One to skip patterns using the lower bound of their distance, zero otherwise
! NO    - Offset of each observation of the patterns in NEAR (NY+1), see twdtwnearest_near, not used if USELB is zero
! NEAR  - Index of the near dates of the observations of the patterns, not used if USELB is zero
! FW    - Lowest time weight of the dates that are not near (NY), not used if USELB is zero
! K     - Number of nearest neighbours
! KI    - Output index of the K nearest patterns (NS,K), not modified for time series without observations.
!         The neighbours of time series P1-1 are read, if they were already computed
! KD    - Output distance to the K nearest patterns (NS,K), not modified for time series without observations
! IERR  - Output error flag, zero if the work arrays were allocated
subroutine twdtwnearest(XM, NS, NT, ND, P1, P2, YM, YO, NY, NP, WL, USELB, NO, NEAR, FW, K, KI, KD, IERR) &
    bind(C, name="twdtwnearest")
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_c_binding
  implicit none

  ! I/O Variables
  integer(c_int), intent(in) :: NS, NT, ND, P1, P2, NY, NP, USELB, K
  integer(c_int), intent(in) :: YO(NP+1), NO(*), NEAR(*)
  real(c_double), intent(in) :: XM(NS,NT,ND), YM(NY,ND), WL(NT,NY), FW(*)
  integer(c_int), intent(inout) :: KI(NS,K)
  real(c_double), intent(inout) :: KD(NS,K)
  integer(c_int), intent(out) :: IERR
  ! Internals
  ! XB - Observations of a time series without missing values (M,ND)
  ! JX - Index of the dates of the observations in XB (M)
  ! PX - Position in XB of each date of the time series, zero if the observation is missing (NT)
  ! DT - Local cost between an observation of a pattern and the time series (M)
  ! CR - Cumulative cost of a row of the TWDTW matrix (M)
  ! BD - Distance to the nearest patterns of a time series in ascending order (K)
  ! BI - Index of the nearest patterns of a time series (K)
  ! FIRST - Nearest patterns of the previous time series (NF)
  ! SEEN  - Last time series for which a pattern was processed in step 1 (NP)
  ! RB - Lower bound of the cost of each observation of the patterns, replaced by
  !      the bound of the cost of the observations after it before the pattern is processed (NY)
  ! PB - Lower bound of the distance of the patterns that remain for step 3 (NP)
  ! SV - Patterns that remain for step 3, a heap with the lowest bound first (NSV)
  ! LBP, LOW   - Lower bound of the distance of a pattern and the lowest bound of a time series
  ! CP, CMIN   - Cost of a cell and lowest cost of a row
  ! DIAG, LEFT, UP - Cumulative cost of the neighbouring cells
  ! SAFE - Margin for the rounding of the bounds, which are sums in a different order than the distance
  real(c_double), allocatable :: XB(:,:), DT(:), CR(:), BD(:), RB(:), PB(:)
  integer(c_int), allocatable :: JX(:), PX(:), BI(:), FIRST(:), SEEN(:), SV(:)
  real(c_double) :: XV, YV, CP, CMIN, DIAG, LEFT, UP, INF, RMIN, S, LBP, LOW, REM
  real(c_double), parameter :: SAFE = 1.0d0 - 1.0d-9
  integer(c_int) :: P, I, J, L, M, IP, NB, POS, Q, JP, JD, NF, STEP, STAGE, NSV, NR
  logical :: COMPLETE, ABANDON, BOUND

  ! The work arrays of the bounds are only needed with USELB
  NR = 1
  if (USELB.eq.1) NR = NY
  allocate(XB(NT,ND), DT(NT), CR(NT), BD(K), RB(NR), PB(NP), JX(NT), PX(NT), BI(K), FIRST(K), SEEN(NP), SV(NP), stat=IERR)
  if (IERR.ne.0) return

  INF = ieee_value(1.0d0, ieee_positive_inf)
  SEEN = 0

  ! Nearest patterns of the time series before the first one, if available
  NF = 0
  if (P1.gt.1) then
    do J = 1, K
      if (KI(P1-1,J).lt.1 .or. KI(P1-1,J).gt.NP) exit
      FIRST(J) = KI(P1-1,J)
      NF = J
    end do
  end if

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
        PX(J) = M
      else
        PX(J) = 0
      end if
    end do

    ! Time series without observations keep the initial values of KI and KD
    if (M.eq.0) cycle

    NB = 0
    NSV = 0
    LOW = INF
    STAGE = 1
    STEP = 0
    patterns: do

      ! Select the next pattern IP to be processed
      BOUND = .false.
      if (STAGE.eq.1) then

        ! Step 1: the nearest patterns of the previous time series
        STEP = STEP + 1
        if (STEP.gt.NF) then
          STAGE = 2
          STEP = 0
          cycle patterns
        end if
        IP = FIRST(STEP)
        SEEN(IP) = P

      else if (STAGE.eq.2) then

        ! Step 2: the other patterns in the order of their index
        STEP = STEP + 1
        if (STEP.gt.NP) then
          STAGE = 3
          do I = NSV / 2, 1, -1
            call sift_down(I)
          end do
          cycle patterns
        end if
        IP = STEP
        if (SEEN(IP).eq.P) cycle patterns

        ! Patterns without observations cannot be a neighbour
        if (YO(IP+1).le.YO(IP)) cycle patterns

        ! Lower bound of the distance, as soon as K patterns with a finite distance were found
        if (USELB.eq.1 .and. NB.eq.K) then
          if (BD(K).lt.INF) then
            LBP = 0.0d0
            do I = YO(IP) + 1, YO(IP+1)
              RMIN = FW(I)
              do Q = NO(I) + 1, NO(I+1)
                JD = NEAR(Q)
                JP = PX(JD)
                if (JP.eq.0) cycle
                S = 0.0d0
                do L = 1, ND
                  S = S + (YM(I,L) - XB(JP,L))**2
                end do
                RMIN = min(RMIN, sqrt(S) + WL(JD,I))
              end do
              RB(I) = RMIN
              LBP = LBP + RMIN
              ! The pattern cannot be among the K nearest
              if (LBP * SAFE .gt. BD(K)) cycle patterns
            end do
            if (LBP.ge.LOW) then
              ! Keep the pattern for step 3
              NSV = NSV + 1
              SV(NSV) = IP
              PB(IP) = LBP
              cycle patterns
            end if
            ! The most promising pattern so far is processed at once,
            ! which can lower the distance to the K-th nearest pattern
            LOW = LBP
            BOUND = .true.
          end if
        end if

      else

        ! Step 3: the remaining patterns in ascending order of their bound
        if (NSV.eq.0) exit patterns
        IP = SV(1)
        if (PB(IP) * SAFE .gt. BD(K)) exit patterns
        SV(1) = SV(NSV)
        NSV = NSV - 1
        call sift_down(1)
        BOUND = .true.

      end if

      ! Replace the bound of each observation of the pattern by the bound
      ! of the cost of the observations after it
      if (BOUND) then
        REM = 0.0d0
        do I = YO(IP+1), YO(IP) + 1, -1
          S = RB(I)
          RB(I) = REM
          REM = REM + S
        end do
      end if

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
        ! The first column has only the cost from the previous row. It is
        ! constrained by the maximum elapsed time like the other columns,
        ! consistent with twdtwf90
        DIAG = CR(1)
        LEFT = DT(1) + WL(JX(1),I) + DIAG
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

        ! Stop if the pattern cannot be among the K nearest: its distance is
        ! not lower than the lowest cost of the row, and not lower than this
        ! cost plus the bound of the cost of the remaining observations
        if (NB.eq.K) then
          if (CMIN.gt.BD(K)) then
            ABANDON = .true.
          else if (CMIN.ge.BD(K) .and. IP.gt.BI(K)) then
            ABANDON = .true.
          else if (BOUND) then
            if ((CMIN + RB(I)) * SAFE .gt. BD(K)) ABANDON = .true.
          end if
          if (ABANDON) exit
        end if

      end do
      if (ABANDON) cycle patterns

      ! CMIN is the lowest cost of the last row, i.e. the TWDTW distance
      ! Insert the pattern in the list of nearest patterns sorted by distance,
      ! patterns with the same distance are sorted by their index
      if (NB.lt.K) NB = NB + 1
      POS = NB
      do while (POS.gt.1)
        if (BD(POS-1).lt.CMIN) exit
        if (BD(POS-1).le.CMIN .and. BI(POS-1).lt.IP) exit
        BD(POS) = BD(POS-1)
        BI(POS) = BI(POS-1)
        POS = POS - 1
      end do
      BD(POS) = CMIN
      BI(POS) = IP

    end do patterns

    do J = 1, NB
      KI(P,J) = BI(J)
      KD(P,J) = BD(J)
      FIRST(J) = BI(J)
    end do
    NF = NB

  end do

  deallocate(XB, DT, CR, BD, RB, PB, JX, PX, BI, FIRST, SEEN, SV)

contains

  ! Moves the pattern at position START of the heap SV down to its place,
  ! such that each pattern has a bound not greater than the patterns below it
  subroutine sift_down(START)
    integer(c_int), intent(in) :: START
    integer(c_int) :: A, B, T
    A = START
    do
      B = 2 * A
      if (B.gt.NSV) exit
      if (B.lt.NSV) then
        if (PB(SV(B+1)).lt.PB(SV(B))) B = B + 1
      end if
      if (PB(SV(A)).le.PB(SV(B))) exit
      T = SV(A)
      SV(A) = SV(B)
      SV(B) = T
      A = B
    end do
  end subroutine sift_down

end subroutine twdtwnearest


! Local distance for the batch computation of TWDTW with a time-weight function
!
! Computes the Euclidean distance between each observation of a block of time
! series and each observation of one temporal pattern. This is the distance
! that a time-weight function receives, see twdtwnearest_cost. Observations of a
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
subroutine twdtwnearest_dist(XM, NS, NT, ND, P1, P2, YM, NY, Y1, Y2, DM, OK) bind(C, name="twdtwnearest_dist")
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
end subroutine twdtwnearest_dist


! Batch computation of TWDTW for one temporal pattern given the local cost
!
! Computes the TWDTW distance between each time series of a block and one
! temporal pattern, and updates the K patterns with the lowest distance of
! each time series. CM has the local cost, i.e. the value of a time-weight
! function for the distance and the elapsed time between each observation of
! the time series and each observation of the pattern. Call this subroutine
! once for each pattern, in the order of the patterns.
!
! The cumulative cost is computed as in twdtwnearest. However, the computation
! never stops before the last observation of the pattern, because the cost
! given by a time-weight function can be negative.
!
! CM   - Local cost (NB,NT,NI)
! NB   - Number of time series in the block
! NT   - Number of dates of the time series
! NI   - Number of observations of the temporal pattern
! OK   - Flag of the observations of the time series (NB,NT), one if complete, see twdtwnearest_dist
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
subroutine twdtwnearest_cost(CM, NB, NT, NI, OK, EL, LB, IP, K, NS, P1, KN, KI, KD, IERR) bind(C, name="twdtwnearest_cost")
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
      ! The first column has only the cost from the previous row. It is
      ! constrained by the maximum elapsed time like the other columns,
      ! consistent with twdtwf90
      DIAG = CR(1)
      if (EL(JX(1),I).gt.LB) then
        LEFT = INF
      else
        LEFT = CM(P,JX(1),I) + DIAG
      end if
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
end subroutine twdtwnearest_cost
