! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate contribution to visibility from mineral dust.
! Uses this to calculate vis in dust and total vis (inc. dust and ppn)

SUBROUTINE DUST_VIS (vis_incppn,    & !in
                     tscreen,       & !in
                     pstar,         & !in
                     dust1,         & !in
                     dust2,         & !in
                     dust3,         & !in
                     dust4,         & !in
                     dust5,         & !in
                     dust6,         & !in   
                     vis_tot,       & !inout                  
                     vis_dust,      & !inout
                     ErrorStatus )    !inout

! Description: Calculates visibility contribution from dust. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!

USE atmos_constants_mod, ONLY: r
    
USE visbty_constants_mod, ONLY:     &
  LnLiminalContrast, RecipVisAir

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:             &
  ST_VisTot, MO8_VisTot, PP_VisTot, &
  ST_VisDust,MO8_VisDust,PP_VisDust

IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: vis_incppn  ! 1.5m Vis incl. ppn
TYPE(PP_Field_type), INTENT(IN) :: tscreen     ! 1.5m Temperature
TYPE(PP_Field_type), INTENT(IN) :: pstar       ! Surface pressure
TYPE(PP_Field_type), INTENT(IN) :: dust1       ! Level 1 div 1 dust MMR
TYPE(PP_Field_type), INTENT(IN) :: dust2       ! Level 1 div 2 dust MMR
TYPE(PP_Field_type), INTENT(IN) :: dust3       ! Level 1 div 3 dust MMR
TYPE(PP_Field_type), INTENT(IN) :: dust4       ! Level 1 div 4 dust MMR
TYPE(PP_Field_type), INTENT(IN) :: dust5       ! Level 1 div 5 dust MMR
TYPE(PP_Field_type), INTENT(IN) :: dust6       ! Level 1 div 6 dust MMR
TYPE(PP_Field_type), INTENT(INOUT) :: vis_tot  ! 1.5m Vis incl. ppn and dust (m)
TYPE(PP_Field_type), INTENT(INOUT) :: vis_dust ! 1.5m Vis in dust only (m)

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DUST_VIS"

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
INTEGER :: i,j           ! Loop counters

! extinction in clean air
REAL :: beta_air
! specific extinction coefficients at 550nm for each bin
REAL, PARAMETER, DIMENSION(6)  :: k_ext=(/ 652.797, 3626.70, 979.290,  &
                                           260.675, 77.6338, 23.9075 /)
! air density
REAL :: rho
! extinction due to dust
REAL :: beta_dust (vis_incppn%Hdr%NumCols, vis_incppn%Hdr%NumRows)
! running total of the extinction
REAL :: beta_tot (vis_incppn%Hdr%NumCols, vis_incppn%Hdr%NumRows)

! End of header -----------------------------------------------

CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( vis_tot % RData ) ) THEN
  DEALLOCATE( vis_tot % RData )
END IF
IF ( ASSOCIATED( vis_dust % RData ) ) THEN
  DEALLOCATE( vis_dust % RData )
END IF

vis_tot % Hdr = vis_incppn % Hdr
vis_tot % Hdr % BMDI     = RMDI
vis_tot % Hdr % PPCode   =  PP_VisTot
vis_tot % Hdr % MO8Type  = MO8_VisTot
vis_tot % Hdr % STCode   =  ST_VisTot

vis_dust % Hdr = vis_tot % Hdr
vis_dust % Hdr % PPCode   =  PP_VisDust
vis_dust % Hdr % MO8Type  = MO8_VisDust
vis_dust % Hdr % STCode   =  ST_VisDust

ALLOCATE( vis_tot % RData(vis_tot % Hdr % NumCols, &
                        vis_tot % Hdr % NumRows) )
ALLOCATE( vis_dust % RData(vis_dust % Hdr % NumCols, &
                        vis_dust % Hdr % NumRows) )

! Calculate the extinction in clean air from RecipVisAir
beta_air = -LnLiminalContrast * RecipVisAir

! Calculate the extinction due to dust

beta_dust(:,:)=0.0
DO j = 1, vis_tot % Hdr % NumRows
  DO i = 1, vis_tot % Hdr % NumCols

    rho = pstar%RData(i,j) / (tscreen%RData(i,j) * r)

    beta_dust(i,j)= rho * (             dust1%RData(i,j)*k_ext(1) + &
            dust2%RData(i,j)*k_ext(2) + dust3%RData(i,j)*k_ext(3) + &
            dust4%RData(i,j)*k_ext(4) + dust5%RData(i,j)*k_ext(5) + &
                                        dust6%RData(i,j)*k_ext(6))
  END DO
END DO

! Calculate fields

DO j = 1, vis_tot % Hdr % NumRows
  DO i = 1, vis_tot % Hdr % NumCols

!   Invert visibility to calculate total extinction
    beta_tot(i,j) = -LnLiminalContrast / vis_incppn%RData(i,j)

!   Add the extinction from dust to the total
    beta_tot(i,j) = beta_tot(i,j) + beta_dust(i,j)

!   Invert back to get visibilities 
    vis_tot%RData(i,j)  = -LnLiminalContrast / beta_tot(i,j)
!   Include small contribution from air to limit vis model's max value
    vis_dust%RData(i,j) = -LnLiminalContrast / (beta_dust(i,j) + beta_air)
  END DO
END DO

9999 CONTINUE

CALL Timer( RoutineName, 4 )

END SUBROUTINE DUST_VIS
