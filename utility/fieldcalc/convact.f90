! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine used in Fieldcalc CB actions

!=======================================================================

SUBROUTINE ConvAct( NumLevs ,       &  ! in
                    Factor,         &  ! in
                    CPNRT_Fld,      &  ! in
                    BlkCld_Flds,    &  ! in
                    ConCld_Flds,    &  ! in
                    Ptheta_Flds,    &  ! in
                    Ttheta_Flds,    &  ! in
                    P_CBB_Fld,      &  ! out
                    P_CBT_Fld,      &  ! out
                    I_CBB_Fld,      &  ! out
                    I_CBT_Fld,      &  ! out
                    P_ECBB_Fld,     &  ! out
                    P_ECBT_Fld,     &  ! out
                    I_ECBB_Fld,     &  ! out
                    I_ECBT_Fld,     &  ! out
                    CBHorE_Fld,     &  ! out
                    ErrorStatus)       ! inout
!
! Description:
!   Calculates CB Action fields
!
! Method:
!  Subroutine ConvAct takes fields of convective cloud combined with the
!  field of convective precipitation rate to produce a mask of Cumulonimbus
!  activity. This Cb mask is then compared to convective cloud top and base
!  pressures to give Cb cloud top and base pressures. These fields can later
!  be converted to heights (in Kft) using the ICAO_HT action.
!
!
!  Contents:
!    P_CBB_Fld= Will contain Cb base pressure
!    P_CBT_Fld= Will contain Cb top pressure
!    I_CBB_Fld= In this subroutine, used to contain CB mask
!                (after this subroutine has finished this will
!                later be used to contain the Cb base height)
!    I_CBT_Fld= Used as a dummy/spare field space for adjusting fields
!                (after this subroutine has finished this will
!                later be used to contain Cb top height)
!    P_ECBB_Fld= Will contain embedded Cb base pressure
!    P_ECBT_Fld= Will contain embedded Cb top pressure
!    I_ECBB_Fld= In this subroutine, used to contain embedded CB mask
!                (after this subroutine has finished this
!                later be used to contain embedded Cb base height)
!    I_ECBT_Fld= Used as a spare field space for adjusting fields
!                (after this subroutine has finished this will
!                later be used to contain embedded Cb top height)
!    CBHorE_Fld= Will contain Cb horizontal extent (as index)
!
!  Input :
!    CPNRT_Fld : Convective Precipitation Rate  (5/205)
!
!    BlkCld_Flds  : Bulk Cloud Fraction     (0/266)
!    ConCld_Flds  : Convective Cloud Amount (5/212)
!    Ptheta_Flds  : Theta level Pressure    (0/408)
!    Ttheta_Flds  : Theta level Temperature (16/004)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!

USE IO_mod, ONLY:        &
    PP_Header_type,      &
    PP_Field_type
USE Err_Mod, ONLY:       &
    StatusOK
USE FldCodes_mod, ONLY:  &
ST_P_CBB,   MO8_P_CBB,  PP_P_CBB,  VC_P_CBB, &
ST_P_CBT,   MO8_P_CBT,  PP_P_CBT,  VC_P_CBT, &
ST_I_CBB,   MO8_I_CBB,  PP_I_CBB,  VC_I_CBB, &
ST_I_CBT,   MO8_I_CBT,  PP_I_CBT,  VC_I_CBT, &
ST_P_ECBB,  MO8_P_ECBB, PP_P_ECBB, VC_P_ECBB,&
ST_P_ECBT,  MO8_P_ECBT, PP_P_ECBT, VC_P_ECBT,&
ST_I_ECBB,  MO8_I_ECBB, PP_I_ECBB, VC_I_ECBB,&
ST_I_ECBT,  MO8_I_ECBT, PP_I_ECBT, VC_I_ECBT,&
ST_CBHorE,  MO8_CBHorE, PP_CBHorE, VC_CBHorE,&
ST_CPNRT,   MO8_CPNRT,  PP_CPNRT,  VC_CPNRT, &
LV_Special

IMPLICIT None


INTEGER, INTENT(IN)                :: Numlevs         !No. of levels
REAL, INTENT(IN)                   :: Factor          ! Options in
INTEGER, INTENT(INOUT)             :: Errorstatus
TYPE(PP_Field_type), INTENT(IN)    :: CPNRT_Fld            ! in
TYPE(PP_Field_type), INTENT(IN)    :: BlkCld_Flds(NumLevs) ! in
TYPE(PP_Field_type), INTENT(IN)    :: ConCld_Flds(NumLevs) ! in
TYPE(PP_Field_type), INTENT(IN)    :: Ptheta_Flds(NumLevs) ! in
TYPE(PP_Field_type), INTENT(IN)    :: Ttheta_Flds(NumLevs) ! in
TYPE(PP_Field_type), INTENT(OUT)   :: P_CBB_Fld            ! out
TYPE(PP_Field_type), INTENT(OUT)   :: P_CBT_Fld            ! out
TYPE(PP_Field_type), INTENT(OUT)   :: I_CBB_Fld            ! out
TYPE(PP_Field_type), INTENT(OUT)   :: I_CBT_Fld            ! out
TYPE(PP_Field_type), INTENT(OUT)   :: P_ECBB_Fld           ! out
TYPE(PP_Field_type), INTENT(OUT)   :: P_ECBT_Fld           ! out
TYPE(PP_Field_type), INTENT(OUT)   :: I_ECBB_Fld           ! out
TYPE(PP_Field_type), INTENT(OUT)   :: I_ECBT_Fld           ! out
TYPE(PP_Field_type), INTENT(OUT)   :: CBHorE_Fld            ! out


! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ConvAct"
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------


! Local Variables:
INTEGER :: i, j, k         ! Loop counters
INTEGER :: ii, jj          ! Loop counters in optimisation loops
INTEGER :: z               ! Loop counter for field index
REAL    :: pcnt            ! Holds value when determining Cb horiz extent

! Parameters used in connection with the subroutine Fill_Pressure
INTEGER :: In_field             ! index for input field to Fill_Pressure
INTEGER :: In_mask              ! index for input mask field to Fill_Pressure
INTEGER :: Out_field            ! index for input field to Fill_Pressure
INTEGER :: Nofill               ! flag, set to 1 to call Fill_Pressure to
                                !   fill in values in area of mask

! Variables minnm, minnmb, diffn, diffnb below used in to check cloud tops
! are above cloud bases
!
REAL :: diffn                   ! difference of CB base pressure
                                !  and CB top pressure
REAL :: diffnb                  ! difference of embedded CB base pressure
                                !  and embedded CB top pressure
REAL :: minnm                   ! indicates CB base is above CB top if negative
REAL :: minnmb                  ! indicates embedded CB base is above embedded
                                !  CB top if negative

!  Variables p & q below are used to contain field values in a 3x3 box
!  around the point of interest (which is p5/q5)
!
!  p1,p2,p3,
!  p4,p5,p6,
!  p7,p8,p9,
!
!  q1,q2,q3,
!  q4,q5,q6,
!  q7,q8,q9
!
REAL :: p1, p2, p3, p4, p5, p6, p7, p8, p9
REAL :: q1, q2, q3, q4, q5, q6, q7, q8, q9


! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Deallocate previous data.
IF (ASSOCIATED(P_CBB_Fld % RData)) THEN
  DEALLOCATE(P_CBB_Fld % RData)
END IF
IF (ASSOCIATED(P_CBT_Fld % RData)) THEN
  DEALLOCATE(P_CBT_Fld % RData)
END IF
IF (ASSOCIATED(I_CBB_Fld % RData)) THEN
  DEALLOCATE(I_CBB_Fld % RData)
END IF
IF (ASSOCIATED(I_CBT_Fld % RData)) THEN
  DEALLOCATE(I_CBT_Fld % RData)
END IF
IF (ASSOCIATED(P_ECBB_Fld % RData)) THEN
  DEALLOCATE(P_ECBB_Fld % RData)
END IF
IF (ASSOCIATED(P_ECBT_Fld % RData)) THEN
  DEALLOCATE(P_ECBT_Fld % RData)
END IF
IF (ASSOCIATED(I_ECBB_Fld % RData)) THEN
  DEALLOCATE(I_ECBB_Fld % RData)
END IF
IF (ASSOCIATED(I_ECBT_Fld % RData)) THEN
  DEALLOCATE(I_ECBT_Fld % RData)
END IF
IF (ASSOCIATED(CBHorE_Fld % RData)) THEN
  DEALLOCATE(CBHorE_Fld % RData)
END IF

!---- Allocate output fields and set fieldcodes etc. in headers ----
P_CBB_Fld % LookupPos      = Ptheta_Flds(1) % LookupPos
P_CBB_Fld % Hdr            = Ptheta_Flds(1) % Hdr
P_CBB_Fld % Hdr % STCode   = ST_P_CBB
P_CBB_Fld % Hdr % MO8Type  = MO8_P_CBB
P_CBB_Fld % Hdr % MO8Level = LV_Special
P_CBB_Fld % Hdr % PPcode   = PP_P_CBB
P_CBB_Fld % Hdr % LBVC     = VC_P_CBB
P_CBB_Fld % Hdr % BMDI     = RMDI
! Reset REAL part of header to correct single level info.
P_CBB_Fld % Hdr % BULev    = 0.0
P_CBB_Fld % Hdr % BHULev   = 0.0
P_CBB_Fld % Hdr % RLevel   = -1.0
P_CBB_Fld % Hdr % BHLev    = 0.0
P_CBB_Fld % Hdr % BHRLev   = 0.0
ALLOCATE( P_CBB_Fld % RData(P_CBB_Fld % Hdr % NumCols,     &
                            P_CBB_Fld % Hdr % NumRows) )
P_CBB_Fld % RData(:,:) = RMDI
! --
P_CBT_Fld % LookupPos      = P_CBB_Fld % LookupPos
P_CBT_Fld % Hdr            = P_CBB_Fld % Hdr
P_CBT_Fld % Hdr % STCode   = ST_P_CBT
P_CBT_Fld % Hdr % MO8Type  = MO8_P_CBT
P_CBT_Fld % Hdr % MO8Level = LV_Special
P_CBT_Fld % Hdr % PPcode   = PP_P_CBT
P_CBT_Fld % Hdr % LBVC     = VC_P_CBT
P_CBT_Fld % Hdr % BMDI     = RMDI
ALLOCATE( P_CBT_Fld % RData(P_CBT_Fld % Hdr % NumCols,     &
                            P_CBT_Fld % Hdr % NumRows) )
P_CBT_Fld % RData(:,:) = RMDI
! --
P_CBB_Fld % LookupPos      = P_CBB_Fld % LookupPos
I_CBB_Fld % Hdr            = P_CBB_Fld % Hdr
I_CBB_Fld % Hdr % STCode   = ST_I_CBB
I_CBB_Fld % Hdr % MO8Type  = MO8_I_CBB
I_CBB_Fld % Hdr % MO8Level = LV_Special
I_CBB_Fld % Hdr % PPcode   = PP_I_CBB
I_CBB_Fld % Hdr % LBVC     = VC_I_CBB
I_CBB_Fld % Hdr % BMDI     = RMDI
ALLOCATE( I_CBB_Fld % RData(I_CBB_Fld % Hdr % NumCols,     &
                            I_CBB_Fld % Hdr % NumRows) )
I_CBB_Fld % RData(:,:) = RMDI
! --
I_CBT_Fld % LookupPos      = P_CBB_Fld % LookupPos
I_CBT_Fld % Hdr            = P_CBB_Fld % Hdr
I_CBT_Fld % Hdr % STCode   = ST_I_CBT
I_CBT_Fld % Hdr % MO8Type  = MO8_I_CBT
I_CBT_Fld % Hdr % MO8Level = LV_Special
I_CBT_Fld % Hdr % PPcode   = PP_I_CBT
I_CBT_Fld % Hdr % LBVC     = VC_I_CBT
I_CBT_Fld % Hdr % BMDI     = RMDI
ALLOCATE( I_CBT_Fld % RData(I_CBT_Fld % Hdr % NumCols,     &
                            I_CBT_Fld % Hdr % NumRows) )
I_CBT_Fld % RData(:,:) = RMDI
!--
P_ECBB_Fld % LookupPos      = P_CBB_Fld % LookupPos
P_ECBB_Fld % Hdr            = P_CBB_Fld % Hdr
P_ECBB_Fld % Hdr % STCode   = ST_P_ECBB
P_ECBB_Fld % Hdr % MO8Type  = MO8_P_ECBB
P_ECBB_Fld % Hdr % MO8Level = LV_Special
P_ECBB_Fld % Hdr % PPcode   = PP_P_ECBB
P_ECBB_Fld % Hdr % LBVC     = VC_P_ECBB
P_ECBB_Fld % Hdr % BMDI     = RMDI
ALLOCATE( P_ECBB_Fld % RData(P_ECBB_Fld % Hdr % NumCols,     &
                             P_ECBB_Fld % Hdr % NumRows) )
P_ECBB_Fld % RData(:,:) = RMDI
! --
P_ECBT_Fld % LookupPos      = P_CBB_Fld % LookupPos
P_ECBT_Fld % Hdr            = P_CBB_Fld % Hdr
P_ECBT_Fld % Hdr % STCode   = ST_P_ECBT
P_ECBT_Fld % Hdr % MO8Type  = MO8_P_ECBT
P_ECBT_Fld % Hdr % MO8Level = LV_Special
P_ECBT_Fld % Hdr % PPcode   = PP_P_ECBT
P_ECBT_Fld % Hdr % LBVC     = VC_P_ECBT
P_ECBT_Fld % Hdr % BMDI     = RMDI
ALLOCATE( P_ECBT_Fld % RData(P_ECBT_Fld % Hdr % NumCols,     &
                             P_ECBT_Fld % Hdr % NumRows) )
P_ECBT_Fld % RData(:,:) = RMDI
! --
I_ECBB_Fld % LookupPos      = P_CBB_Fld % LookupPos
I_ECBB_Fld % Hdr            = P_CBB_Fld % Hdr
I_ECBB_Fld % Hdr % STCode   = ST_I_ECBB
I_ECBB_Fld % Hdr % MO8Type  = MO8_I_ECBB
I_ECBB_Fld % Hdr % MO8Level = LV_Special
I_ECBB_Fld % Hdr % PPcode   = PP_I_ECBB
I_ECBB_Fld % Hdr % LBVC     = VC_I_ECBB
I_ECBB_Fld % Hdr % BMDI     = RMDI
ALLOCATE( I_ECBB_Fld % RData(I_ECBB_Fld % Hdr % NumCols,     &
                             I_ECBB_Fld % Hdr % NumRows) )
I_ECBB_Fld % RData(:,:) = RMDI
! --
I_ECBT_Fld % LookupPos      = P_CBB_Fld % LookupPos
I_ECBT_Fld % Hdr            = P_CBB_Fld % Hdr
I_ECBT_Fld % Hdr % STCode   = ST_I_ECBT
I_ECBT_Fld % Hdr % MO8Type  = MO8_I_ECBT
I_ECBT_Fld % Hdr % MO8Level = LV_Special
I_ECBT_Fld % Hdr % PPcode   = PP_I_ECBT
I_ECBT_Fld % Hdr % LBVC     = VC_I_ECBT
I_ECBT_Fld % Hdr % BMDI     = RMDI
ALLOCATE( I_ECBT_Fld % RData(I_ECBT_Fld % Hdr % NumCols,     &
                             I_ECBT_Fld % Hdr % NumRows) )
I_ECBT_Fld % RData(:,:) = RMDI
! --
CBHorE_Fld % LookupPos      = P_CBB_Fld % LookupPos
CBHorE_Fld % Hdr            = P_CBB_Fld % Hdr
CBHorE_Fld % Hdr % STCode   = ST_CBHorE
CBHorE_Fld % Hdr % MO8Type  = MO8_CBHorE
CBHorE_Fld % Hdr % MO8Level = LV_Special
CBHorE_Fld % Hdr % PPcode   = PP_CBHorE
CBHorE_Fld % Hdr % LBVC     = VC_CBHorE
CBHorE_Fld % Hdr % BMDI     = RMDI
ALLOCATE( CBHorE_Fld % RData(CBHorE_Fld % Hdr % NumCols,     &
                             CBHorE_Fld % Hdr % NumRows) )
CBHorE_Fld % RData(:,:) = RMDI

!---- Generate initial Cb and layer cloud masks as follows:
!
! If convective PPTN rate greater than .00025 then Cb. Put into I_CBB_Fld.
! If layer cloud > .95 on any model level then layer cloud. Put into I_ECBB_Fld.
!   This mask will be used to generate embedded Cb mask later.
!
DO j = 1, P_CBB_Fld % Hdr % NumRows
  DO i = 1, P_CBB_Fld % Hdr % NumCols
    IF ( CPNRT_Fld % RData(i,j)   > 0.00025    )  THEN
       I_CBB_Fld % RData(i,j) = 1.
    ELSE
       I_CBB_Fld % RData(i,j) = RMDI
    END IF
  END DO
END DO
DO k = 1, NumLevs
  DO j = 1, I_ECBB_Fld % Hdr % NumRows
    DO i = 1, I_ECBB_Fld % Hdr % NumCols
      IF  ( BlkCld_Flds(k) % RData(i,j)   > 0.95  )  THEN
        I_ECBB_Fld % RData(i,j) = 1.
      END IF
    END DO
  END DO
END DO


!---  Adjust mask fields using subroutine FILL_N_DSPEC to fill in small gaps
!     and smooth edges. Subroutine iterates 30 times.
!
! Note: For Cb mask, need to copy mask field into I_CBT_Fld first
!       as this is used in the iteration process in FILL_N_DSPEC
!       but don't need to do this for the layer cloud mask.

! Cb mask
I_CBT_Fld % RData(:,:) = I_CBB_Fld % RData(:,:)

IF (Factor /= -1.0) THEN
! DEPENDS ON: fill_n_dspec
  Call FILL_N_DSPEC(30,I_CBB_Fld,I_CBT_Fld,ErrorStatus)
END IF

I_CBT_Fld % RData(:,:) = I_CBB_Fld % RData(:,:)


! Embedded Cb mask
IF (Factor /= -1.0) THEN
! DEPENDS ON: fill_n_dspec
  Call FILL_N_DSPEC(30,I_ECBB_Fld,I_ECBT_Fld,ErrorStatus)
END IF
I_ECBT_Fld % RData(:,:) = I_ECBB_Fld % RData(:,:)
I_ECBB_Fld % RData(:,:) = RMDI



!---- Identify embedded Cbs and generate embedded Cb mask
!
! There are no output fields of layer cloud only with which to identify
! embedded Cbs. In order to classify Cbs as embedded the following
! criteria are used. A cell which has a Cb will be considered as
! embedded if at least one of the 8 adjacent cells is not a Cb AND is
! not clear of any cloud.
!
! I_CBB_Fld contains Cb mask
! I_ECBT_Fld contains layer cloud mask
! Output I_ECBB_Fld will contain the embedded Cb mask

DO j = 2, P_CBB_Fld % Hdr % NumRows  -1
  DO i = 2, P_CBB_Fld % Hdr % NumCols -1
    p1 = I_CBB_Fld % RData(i-1,j+1)
    p2 = I_CBB_Fld % RData(i,j+1)
    p3 = I_CBB_Fld % RData(i+1,j+1)
    p4 = I_CBB_Fld % RData(i-1,j)
    p5 = I_CBB_Fld % RData(i,j)
    p6 = I_CBB_Fld % RData(i+1,j)
    p7 = I_CBB_Fld % RData(i-1,j-1)
    p8 = I_CBB_Fld % RData(i,j-1)
    p9 = I_CBB_Fld % RData(i+1,j-1)
    q1 = I_ECBT_Fld % RData(i-1,j+1)
    q2 = I_ECBT_Fld % RData(i,j+1)
    q3 = I_ECBT_Fld % RData(i+1,j+1)
    q4 = I_ECBT_Fld % RData(i-1,j)
    q5 = I_ECBT_Fld % RData(i,j)
    q6 = I_ECBT_Fld % RData(i+1,j)
    q7 = I_ECBT_Fld % RData(i-1,j-1)
    q8 = I_ECBT_Fld % RData(i,j-1)
    q9 = I_ECBT_Fld % RData(i+1,j-1)
    IF ((((q1 > 0.) .AND. (p1 == RMDI))  .OR. &
         ((q2 > 0.) .AND. (p2 == RMDI))  .OR. &
         ((q3 > 0.) .AND. (p3 == RMDI))  .OR. &
         ((q4 > 0.) .AND. (p4 == RMDI))  .OR. &
         ((q6 > 0.) .AND. (p6 == RMDI))  .OR. &
         ((q7 > 0.) .AND. (p7 == RMDI))  .OR. &
         ((q8 > 0.) .AND. (p8 == RMDI))  .OR. &
         ((q9 > 0.) .AND. (p9 == RMDI)) ) .AND. (p5 == 1.) ) THEN
      I_ECBB_Fld % RData(i,j) = 1.
    END IF
  END DO
END DO



!---- Recheck embedded Cb mask ----
!
! If cell is marked as having an embedded Cb, check the surrounding
! 8 cells in turn - if they contain Cbs, mark them as embedded Cbs too.
!
!Comment on optimisation below:
!kk Optimisation exchanged I and J loop to avoid bank conflicts
!kk In the I loop, there are dependencies in I_ECBB_Fld%RData which
!kk prevent vectorization. Running this loop with stride 3 allows
!kk vectorization, but slightly changes the algorithm .

DO k=1,50
  DO jj=0,2
    DO j = jj+2, P_CBB_Fld % Hdr % NumRows  -1, 3
      DO ii=0,2
        DO i = ii+2, P_CBB_Fld % Hdr % NumCols -1, 3
          p1=I_CBB_Fld % RData(i-1,j+1)
          p2=I_CBB_Fld % RData(i  ,j+1)
          p3=I_CBB_Fld % RData(i+1,j+1)
          p4=I_CBB_Fld % RData(i-1,j  )
          p6=I_CBB_Fld % RData(i+1,j  )
          p7=I_CBB_Fld % RData(i-1,j-1)
          p8=I_CBB_Fld % RData(i  ,j-1)
          p9=I_CBB_Fld % RData(i+1,j-1)

          q5=I_ECBB_Fld % RData(i  ,j  )
          IF ( q5 ==1. ) THEN
            IF ( p1 == 1. ) THEN
              I_ECBB_Fld % RData(i-1,j+1) = 1.
            END IF
            IF ( p2 == 1. ) THEN
              I_ECBB_Fld % RData(i,j+1) = 1.
            END IF
            IF ( p3 == 1. ) THEN
              I_ECBB_Fld % RData(i+1,j+1) = 1.
            END IF
            IF ( p4 == 1. ) THEN
              I_ECBB_Fld % RData(i-1,j) = 1.
            END IF
            IF ( p6 == 1. ) THEN
              I_ECBB_Fld % RData(i+1,j) = 1.
            END IF
            IF ( p7 == 1. ) THEN
              I_ECBB_Fld % RData(i-1,j-1) = 1.
            END IF
            IF ( p8 == 1. ) THEN
              I_ECBB_Fld % RData(i,j-1) = 1.
            END IF
            IF ( p9 == 1. ) THEN
              I_ECBB_Fld % RData(i+1,j-1) = 1.
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
END DO




!---- Identify pressure at base and top of isolated and embedded Cbs ----
! Use:
! I_CBB_Fld  containing Cb mask
! I_ECBB_Fld  containing Cb mask
! ConCld_Flds(k)  containing convective cloud fraction at level k
! Ptheta_Flds(k) containing pressure at level k
!
! Output:
! P_CBB_Fld will contain Cb base pressure
! P_CBT_Fld will contain Cb top pressure
! P_ECBB_Fld will contain embedded Cb base pressure
! P_ECBT_Fld will contain embedded Cb top pressure
!
! Do this by looping over model levels bottom to top.
!
! For cloud bases:
!  If Cb mask indicates Cb in this cell, and convective cloud amount at this
!  level > 0, and no base has yet been identified, set the Cb base to this
!  level. Repeat for embedded Cb using embedded Cb mask.
!
! For cloud tops:
!  If Cb mask indicates Cb in this cell, and convective cloud amount at this
!  level > 0, and cloud base has been set and is lower than or at this level
!  then set the cloud top to be this level. This is overwritten each time the
!  model level increases and these conditions are met, until the conditions
!  aren't met anymore. Repeat for embedded Cb using embedded Cb mask.
!

DO k= 1, numlevs
  DO j = 1, P_CBB_Fld % Hdr % NumRows
    DO i = 1, P_CBB_Fld % Hdr % NumCols
!  find cloud Bases
      IF ( ( ConCld_Flds(k) % RData(i,j)  /= 0. )  .AND.  &
           ( I_CBB_Fld % RData(i,j)     /= RMDI )  .AND.  &
           ( P_CBB_Fld % RData(i,j)     == RMDI )  )  THEN
        P_CBB_Fld % RData(i,j) =  Ptheta_Flds(k) % RData(i,j)
      END IF

      IF ( ( ConCld_Flds(k) % RData(i,j)  /= 0. )  .AND.  &
           ( I_ECBB_Fld % RData(i,j)     /= RMDI )  .AND.  &
           ( P_ECBB_Fld % RData(i,j)     == RMDI )  )  THEN
        P_ECBB_Fld % RData(i,j) = Ptheta_Flds(k) % RData(i,j)

      END IF

! Find cloud tops
      IF ( ( ConCld_Flds(k) % RData(i,j)  /= 0. ) .AND.  &
           ( I_CBB_Fld % RData(i,j)     /= RMDI ) .AND.  &
           ( P_CBB_Fld % RData(i,j) >= Ptheta_Flds(k) % RData(i,j))) THEN
        P_CBT_Fld % RData(i,j) = Ptheta_Flds(k) % RData(i,j)
      END IF

      IF ( ( ConCld_Flds(k) % RData(i,j)  /= 0. ) .AND.  &
           ( I_ECBB_Fld % RData(i,j)     /= RMDI ) .AND.  &
           ( P_ECBB_Fld % RData(i,j) >= Ptheta_Flds(k) % RData(i,j))) THEN
        P_ECBT_Fld % RData(i,j) = Ptheta_Flds(k) % RData(i,j)
      END IF
    END DO
  END DO
END DO



!---- Test that Cb / embedded Cb tops are above bases ----
!
! Do this by calculating the difference in pressure between the cloud base
!  and top. If this is negative, the cloud base is above cloud top. In this
!  case set the flag minnm (for isolated Cb) or minnmb (for embedded Cb)
!  to be equal to the difference. This is then used to reset all values in the
!  Cb base pressure or embedded Cb base pressure to the min value. (I.e. if
!  any one point is wrong, all values in the field are set to the min value.)

minnm=1.
minnmb=1.

DO j = 1, P_CBB_Fld % Hdr % NumRows
   DO i = 1, P_CBB_Fld % Hdr % NumCols
     diffn = P_CBB_Fld % RData(i,j) - P_CBT_Fld % RData(i,j)
     IF (diffn < 0.) THEN
       minnm = diffn
     END IF
     diffnb = P_ECBB_Fld % RData(i,j) - P_ECBT_Fld % RData(i,j)
     IF (diffnb < 0.) THEN
       minnmb = diffnb
     END IF
  END DO
END DO

DO j = 1, P_CBB_Fld % Hdr % NumRows
   DO i = 1, P_CBB_Fld % Hdr % NumCols
     IF (minnm < 0.) THEN
       P_CBB_Fld % RData(i,j) = minnm
     END IF
     IF (minnmb < 0.) THEN
       P_ECBB_Fld % RData(i,j) = minnmb
     END IF
  END DO
END DO



!---- Calculate Cb horizontal extent ----
! Input:
!   I_CBB_Fld    contains Cb mask
!   ConCld_Flds(k) contains convective cloud fraction at level k
!
! Output:
!   CBHorE_Fld    containing Cb horizontal extent, given as value of convective
!                 cloud fraction.
!
!  For each gridpoint:
!  1) Set CBHorE_Fld=largest convective cloud fraction in the column
!

DO k = 1, NumLevs
  DO j = 1, P_CBB_Fld % Hdr % NumRows
    DO i = 1, P_CBB_Fld % Hdr % NumCols
      IF ( (ConCld_Flds(k) % RData(i,j) > CBHorE_Fld % RData(i,j)) .AND. &
           (I_CBB_Fld % RData(i,j) > 0. )  ) THEN
        CBHorE_Fld % RData(i,j) = ConCld_Flds(k) % RData(i,j)
      END IF
    END DO
  END DO
END DO

DO j = 1, P_CBB_Fld % Hdr % NumRows
  DO i = 1, P_CBB_Fld % Hdr % NumCols
    pcnt = RMDI
    IF ( CBHorE_Fld % RData(i,j) > 0.0 ) THEN 
      pcnt = CBHorE_Fld % RData(i,j)
    END IF
    CBHorE_Fld % RData(i,j) = pcnt
  END DO
END DO

!---- Fill in field values ----
!
! If Nofill=1 use Fill_Pressure to check/fill in missing output field
! values where Cb / embedded Cb mask indicates presence of Cb.
!
! Inputs/output for subroutine Fill_Pressure:
! In_field  - index of input field to be filled in. This will be adjusted
!               within subroutine. On output will hold the adjusted field.
! In_mask   - index of mask field used for comparison (will be either 3 or 7)
! Out_field - index of dummy field which will be overwritten and used within
!               subroutine to compare field values.
!
! On output, the adjusted field will be in both In_field and
! Out_field, and either could be used. The code below uses
! In_field for the adjusted field and Out_field as a
! dummy field which is overwritten.
!

Nofill=1
IF ( Nofill == 1 ) THEN


!   Fill in Cloud Base Pressure values within area of Cb mask.
! P_CBB_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_CBB_Fld : Mask field
! I_CBT_Fld : Dummy field, will be overwritten in the subroutine
!             and replaced by adjusted field

! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(P_CBB_Fld,I_CBB_Fld,I_CBT_Fld,ErrorStatus)


!   Fill in Cloud Top Pressure values within area of Cb mask.
! P_CBT_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_CBB_Fld : Mask field
! I_CBT_Fld : Dummy field, will be overwritten in the subroutine
!             and replaced by adjusted field

! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(P_CBT_Fld,I_CBB_Fld,I_CBT_Fld,ErrorStatus)


!   Fill in Cloud Base Pressure values within area of embedded Cb mask.
! P_ECBB_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_ECBB_Fld : Mask field
! I_ECBT_Fld : Dummy field, will be overwritten in the subroutine
!              and replaced by adjusted field

! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(P_ECBB_Fld,I_ECBB_Fld,I_ECBT_Fld,ErrorStatus)


!   Fill in Cloud Top Pressure values within area of embedded Cb mask.

! P_ECBT_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_ECBB_Fld : Mask field
! I_ECBT_Fld : Dummy field, will be overwritten in the subroutine
!              and replaced by adjusted field

! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(P_ECBT_Fld,I_ECBB_Fld,I_ECBT_Fld,ErrorStatus)


!   Fill in Cb horizontal extent values within area of Cb mask.
! CBHorE_Fld : Field to be adjusted (will be replaced by adjusted field)
! I_CBB_Fld  : Mask field
! I_ECBT_Fld : Dummy field, will be overwritten in the subroutine
!              and replaced by adjusted field

! DEPENDS ON: fill_pressure
  CALL Fill_Pressure(CBHorE_Fld,I_CBB_Fld,I_ECBT_Fld,ErrorStatus)


END IF


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE ConvAct
