! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+Decode the STASH level code

Module Rcf_Level_Code_Mod

! Description:
!   Sets ILOUT to an appropriate level size according to the value of
!   ILIN.
!
!******************************************************************
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
! Any Changes to this routine must be accompanied with equivalent
! changes to levcod.F90 and 
!******************************************************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.

Contains

SUBROUTINE Rcf_Level_Code( ILIN, ILOUT, Grid )

Use Rcf_Model_Mod, Only :  &
    STLevGWDrag,           &
    BotVDiffLev,           &
    TopVDiffLev


Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

! The following are retired from the STASHmaster.  Using pseudo levels instead.
USE rad_input_mod, ONLY:    &
    H_SWBands,              &
    H_LWBands

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_HeadAddress_Mod, Only : &
    FH_GridStagger_C,           &
    FH_GridStagger_Endgame

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)    :: ILIN           ! Model level code
Integer, Intent(Out)   :: ILOUT          ! The actual level

Type (Grid_Type), Intent(In) :: Grid     ! The grid that decoded
                                         ! values should correspond to


! Local variables
Integer                      :: ErrorStatus
Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName = 'Rcf_Level_Code'


Select Case ( ILIN )
  Case ( 1 )                            ! First atmos level
    ILOUT=1

  Case ( 2 )                            ! Top atmos level
    ILOUT= Grid % MODEL_LEVELS

  Case ( 3 )                            ! Top wet level
    ILOUT= Grid % WET_LEVELS

  Case ( 4 )
    ILOUT= Grid % MODEL_LEVELS - 1

  Case ( 5 )                            ! First boundary layer level
    ILOUT=MIN(1,Grid % BL_LEVELS)

  Case ( 6 )                            ! Last boundary layer level
    ILOUT=Grid % BL_LEVELS

  Case ( 7 )
    ILOUT= Grid % BL_LEVELS+1

  Case ( 8 )                            ! First soil level
    ILOUT=MIN(1, Grid % ST_LEVELS)

  Case ( 9 )                            ! Last soil level
    ILOUT= Grid % ST_LEVELS

  Case ( 10 )                           ! First tracer level
    ILOUT= Grid % MODEL_LEVELS - Grid % TR_LEVELS+1

  Case ( 11 )                           ! Last tracer level
    ILOUT= Grid % MODEL_LEVELS

  Case ( 12 )
    ILOUT= Grid % MODEL_LEVELS+1

  Case ( 13 )                           ! First gravity wave drag level
    ILOUT=StLevGWdrag

  Case ( 14 )                           ! Last gravity wave drag level
    ILOUT= Grid % MODEL_LEVELS

  Case ( 15 )
    ILOUT=BotVDiffLev

  Case ( 16 )
    ILOUT=TopVDiffLev-1

  Case ( 17 )
    ILOUT=TopVDiffLev

  Case ( 18 )
    ILOUT= Grid % BL_LEVELS-1

  Case ( 19 )
    ILOUT= Grid % MODEL_LEVELS+1

  Case ( 20 )
    ILOUT=MIN(2, Grid % ST_LEVELS)

  ! Ocean removed at vn7.0 so this is redundant
  Case ( 21 )
    ILOUT=1
  
  ! Ocean removed at vn7.0 so this is redundant
  Case ( 22 )
    ILOUT=0

  Case ( 23 )
    ILOUT= Grid % OZONE_LEVELS

  Case ( 24 )
    ILOUT= Grid % MODEL_LEVELS*H_SWBANDS

  Case ( 25 )
    ILOUT=( Grid % MODEL_LEVELS+1)*H_SWBANDS

  Case ( 26 )
    ILOUT= Grid % WET_LEVELS*H_SWBANDS

  Case ( 27 )
    ILOUT= Grid % MODEL_LEVELS*H_LWBANDS

  Case ( 28 )
    ILOUT=( Grid % MODEL_LEVELS+1)*H_LWBANDS

  Case ( 29 )
    ILOUT= Grid % WET_LEVELS*H_LWBANDS

  Case ( 30 )
    ILOUT=2

  Case ( 32 )
    ILOUT=H_SWBANDS

  Case ( 33 )
    ILOUT=H_LWBANDS

  Case ( 34 )
    ILOUT= Grid % SM_LEVELS

  Case ( 35 )
    ILOUT= Grid % CLOUD_LEVELS

  Case ( 36 )                       ! Wave model first level (direction)
    ILOUT=1

  Case ( 37 )                       ! Wave model last level (direction)
                                    ! No wave model in rcf
!    ILOUT=NANG
     ILOUT=0

  Case ( 38 )                       ! Surace theta level
    ILOUT=0

  Case ( 39 )
    ! Number of ISCCP simulator levels
    ILOUT=7

  Case ( 40 )
! Fields which are different between ENDGAME and New Dynamics regarding the
! whether it starts at theta level 0 in ENDGAME but 1 in New Dynamics.  Use
! grid stagger value to identify endgame grid - not perfect.
    IF (grid % grid_stagger == FH_GridStagger_Endgame) THEN
      ILOUT=0
    ELSE
      ILOUT=1
    END IF

  Case Default
    WRITE(Cmessage, '(A, I5)') 'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND: ',ILIN
    ErrorStatus=1
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

RETURN
END Subroutine Rcf_Level_Code
End Module Rcf_Level_Code_Mod

