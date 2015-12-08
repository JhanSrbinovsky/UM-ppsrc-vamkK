! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE LSP_FOCWWIL---------------------------------------------

!    Purpose: Calculate from temperature the Fraction Of Cloud Water
!             Which Is Liquid.
!       NOTE: Operates within range 0 to -9 deg.C based upon MRF
!             observational analysis. Not robust to changes in TM or T0C


!    Programming standard: Unified Model Documentation Paper No 4,
!                          Version 1, dated 12/9/89.

!    Logical component covered: Part of P26.

!    System task:

!    Documentation: Unified Model Documentation Paper No 26: Eq 26.50.

!    Called by components P26, P23.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!    Arguments:---------------------------------------------------------

      SUBROUTINE lsp_focwwil(                                           &
       t,points,rocwwil                                                 &
      )
      USE water_constants_mod, ONLY: tm
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
      INTEGER                                                           &
                       ! Input integer scalar :-
       points          ! IN Number of points to be processed.
      REAL                                                              &
                       ! Input real arrays :-
       t(points)       ! IN Temperature at this level (K).
      REAL                                                              &
                       ! Updated real arrays :-
       rocwwil(points) ! OUT Ratio Of Cloud Water Which Is Liquid.
!     External subprogram called :-
!     EXTERNAL None.

!-----------------------------------------------------------------------
!  Common, then local, physical constants.
!-----------------------------------------------------------------------
      REAL                                                              &
       tstart                                                           &
                       ! Temperature at which ROCWWIL reaches 1.
      ,trange          ! Temperature range over which 0 < ROCWWIL < 1.
      PARAMETER(tstart=tm,                                              &
                trange=9.0)
!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
                       ! Real workspace. At end of DO loop, contains :-
       tfoc            ! T(I) within DO loop. Allows routines to call
!                        LSP_FOCWWIL(WORK1, POINTS, WORK1) to save space
!  (b) Others.
      INTEGER i       ! Loop counter (horizontal field index).

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('LSP_FOCWWIL',zhook_in,zhook_handle)
      DO  i = 1, points

        tfoc = t(i)
!-----------------------------------------------------------------------
!  0. Calculate fraction of cloud water which is liquid (FL),
!     according to equation P26.50.
!-----------------------------------------------------------------------
        IF (tfoc  <=  (tstart - trange)) THEN
!       Low temperatures, cloud water all frozen------------------------
          rocwwil(i) = 0.0

        ELSE IF (tfoc  <   tstart) THEN
!       Intermediate temperatures---------------------------------------
          rocwwil(i) = (tfoc - tstart + trange) / trange

        ELSE
!       High temperatures, cloud water all liquid-----------------------
          rocwwil(i) = 1.0

        END IF

      END DO ! Loop over points

      IF (lhook) CALL dr_hook('LSP_FOCWWIL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE lsp_focwwil
