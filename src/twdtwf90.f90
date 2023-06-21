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


! Computation of TWDTW cost matrix

! XM - matrix with the time series (N,D)
! YM - matrix with the temporal profile (M,D)
! CM - Output cumulative cost matrix
! DM - Direction matrix
! VM - Starting points matrix
! SM - Matrix of step patterns
! N  - Number of rows in CM, DM, and VM - time series
! M  - Number of columns CM, DM, and VM - temporal profile
! D  - Number of spectral dimensions including time in XM and YM
! NS - Number of rows in SM
! TW - Time-Weight parameters alpha and beta
! LB - Constrain TWDTW calculation to band given by TW(2)
subroutine twdtwf90(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB)
  ! I/O Variables
  use, intrinsic :: ieee_arithmetic
  implicit none
  double precision :: ellapsed, distance
  integer, intent(in) :: N, M, D, NS
  integer, intent(in) :: SM(NS,4)
  integer, intent(out) :: DM(N+1,M), VM(N+1,M), JB(N)
  double precision, intent(in) :: XM(M,D), YM(N,D), TW(2)
  double precision, intent(out) :: CM(N+1,M)
  logical, intent(in) :: LB
  ! Internals
  double precision :: W, CP(NS), VMIN, A, B, TD
  integer :: I, J, IL(NS), JL(NS), K, PK, KMIN, ZERO=0, ONE=1, JM, JLMIN, ILMIN
  double precision :: NAN, INF
  NAN = ieee_value(0.0, ieee_quiet_nan)
  INF = ieee_value(0.0, ieee_positive_inf)

  VM(1,1) = 1

  ! Initialize the first row and col of the matrices
  do I = 2, N+1
    TD = ellapsed(YM(I-1,1) - XM(1,1))
    CM(I,1) = CM(I-1,1) + distance(YM, XM, N, M, D, I-1, 1, TW, TD)
    DM(I,1) = 3
    VM(I,1) = 1
  end do

  do J = 2, M
    TD = ellapsed(YM(2,1) - XM(J,1))
    CM(2,J) = CM(2,J-1) + distance(YM, XM, N, M, D, 1, J, TW, TD)
    DM(1,J) = 2
    VM(1,J) = J
  end do

  ! Compute cumulative cost matrix
  J = 2
  do while ( J .le. M )
    I = 2
    do while ( I .le. N+1 )
      ! Calculate local distance
      ! the call takes I-1 because local matrix has an additional row at the beginning
      TD = ellapsed(YM(I-1,1) - XM(J,1))
      if (LB.and.(TD.gt.TW(2))) then
        CM(I,J) = INF
        DM(I,J) = -ONE
        VM(I,J) = ZERO
      else
        CM(I,J) = distance(YM, XM, N, M, D, I-1, J, TW, TD)
      end if
      ! Initialize list of step cost
      CP(:) = NAN
      do K = 1, NS
        PK = SM(K,1)
        IL(K) = I - SM(K,2)
        JL(K) = J - SM(K,3)
        if ((IL(K).gt.ZERO).and.(JL(K).gt.ZERO)) then
          W = SM(K,4)
          if (W.eq.-ONE) then
            CP(PK) = CM(IL(K),JL(K))
          else
            CP(PK) = CP(PK) + CM(IL(K),JL(K))*W
          end if
        end if
      end do
      KMIN = -ONE
      VMIN = INF
      do K = 1, NS
        PK = SM(K,1)
        if (CP(PK).eq.CP(PK).and.CP(PK).lt.VMIN) then
          KMIN = PK
          VMIN = CP(PK)
          ILMIN = IL(K)
          JLMIN = JL(K)
        end if
      end do
      if (KMIN.gt.-ONE) then
        CM(I,J) = VMIN
        DM(I,J) = KMIN
        VM(I,J) = VM(ILMIN, JLMIN)
      end if
      I = I + 1
    end do
    J = J + 1
  end do
  J = 1
  K = ZERO
  do while ( J .le. M )
    if (VM(N+1,J).ne.ZERO) then
      if (K.eq.ZERO) then
        K = 1
        JB(K) = J
        JM = VM(N+1,J)
      elseif (VM(N+1,J).ne.JM) then
        K = K + 1
        JB(K) = J
        JM = VM(N+1,J)
      elseif (CM(N+1,J).lt.CM(N+1,JB(K))) then
        JB(K) = J
      end if
    end if
    J = J + 1
  end do
end subroutine twdtwf90


! Computation of TWDTW cost matrix
! XM - matrix with the time series (N,D)
! YM - matrix with the temporal profile (M,D)
! CM - Output cumulative cost matrix
! DM - Direction matrix
! VM - Starting points matrix
! SM - Matrix of step patterns
! N  - Number of rows in CM, DM, and VM - time series
! M  - Number of columns CM, DM, and VM - temporal profile
! D  - Number of spectral dimensions including time in XM and YM
! NS - Number of rows in SM
! TW - Time-Weight parameters alpha and beta
! LB - Constrain TWDTW calculation to band given by TW(2)
subroutine twdtwf90gt(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB)
  use, intrinsic :: ieee_arithmetic
  implicit none
  double precision :: ellapsed, distance
  integer, intent(in) :: N, M, D, NS
  integer :: SM(NS,4), DM(N+1,M), VM(N+1,M), JB(N)
  double precision, intent(in) :: XM(M,D), YM(N,D), TW(2)
  logical, intent(in) :: LB
  double precision :: CM(N+1,M), W, CP(NS), VMIN, A, B, TD
  integer :: I, J, IL(NS), JL(NS), K, PK, KMIN, ZERO, ONE, JM, ILMIN, JLMIN, IML
  parameter(ZERO=0, ONE=1)
  double precision :: NAN, INF
  NAN = ieee_value(0.0, ieee_quiet_nan)
  INF = ieee_value(0.0, ieee_positive_inf)
  IML = 1
  VM(1,1) = 1

  ! Initialize the first row and column of the matrices
  do 21 I = 2, N+1
     TD = ellapsed(YM(I-1,1) - XM(1,1))
     CM(I,1) = CM(I-1,1) + distance(YM, XM, N, M, D, I-1, 1, TW, TD)
     DM(I,1) = 3
     VM(I,1) = 1
21 continue

  do 31 J = 2, M
     TD = ellapsed(YM(2,1) - XM(J,1))
     CM(2,J) = CM(2,J-1) + distance(YM, XM, N, M, D, 1, J, TW, TD)
     DM(1,J) = 2
     VM(1,J) = J
31 continue

  ! Compute cumulative cost matrix
  J = 2
  do 32 while (J .le. M)
     I = 2
     do 22 while (I .le. N+1)
        TD = ellapsed(YM(I-1,1) - XM(J,1))
        if (LB .and. (TD > TW(2))) then
           CM(I,J) = INF
           DM(I,J) = -ONE
           VM(I,J) = ZERO
           goto 44
        else
           CM(I,J) = distance(YM, XM, N, M, D, I-1, J, TW, TD)
        endif

        ! Initialize list of step cost
        do 10 K = 1, NS
           CP(K) = NAN
10    continue

        do 11 K = 1, NS
           PK = SM(K,1)
           IL(K) = I - SM(K,2)
           JL(K) = J - SM(K,3)
           if ((IL(K) > ZERO) .and. (JL(K) > ZERO)) then
              W = SM(K,4)
              if (W .eq. -ONE) then
                 CP(PK) = CM(IL(K),JL(K))
              else
                 CP(PK) = CP(PK) + CM(IL(K),JL(K)) * W
              endif
           endif
11    continue

        KMIN = -ONE
        VMIN = INF
        ILMIN = -ONE
        JLMIN = -ONE
        do 12 K = 1, NS
           PK = SM(K,1)
           if ((CP(PK) == CP(PK)) .and. (CP(PK) < VMIN)) then
              KMIN = PK
              VMIN = CP(PK)
              ILMIN = IL(K)
              JLMIN = JL(K)
           endif
12    continue

        if (KMIN > -ONE) then
           CM(I,J) = VMIN
           DM(I,J) = KMIN
           VM(I,J) = VM(ILMIN, JLMIN)
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
     if (VM(N+1,J) /= ZERO) then
        if (K == ZERO) then
           K = 1
           JB(K) = J
           JM = VM(N+1,J)
           goto 68
        endif
        if (VM(N+1,J) /= JM) then
           K = K + 1
           JB(K) = J
           JM = VM(N+1,J)
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
