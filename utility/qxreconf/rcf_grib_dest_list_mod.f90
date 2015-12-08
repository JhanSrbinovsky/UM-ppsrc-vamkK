! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

MODULE Rcf_Grib_Dest_List_Mod

! Description: Routine to DeAllocate all entries in a linked list
!              and nullify the pointers which locate the 'ends'.
!
! Method:
!         For lists with more than 1 member-
!           Step through the list deallocating the previous entry
!           Deallocate the end
!           Nullify the pointers to the end
!         For lists with only 1 member-
!           Deallocate the member.
!           Nullify the pointers.

CONTAINS

  SUBROUTINE Rcf_Grib_Dest_List(List_Header)

    USE Rcf_GRIB_Block_Params_Mod, ONLY :                                      &
        List_Marker,                                                           &
        Grib_Record

    USE EReport_Mod, ONLY :                                                    &
        EReport

    USE UM_ParVars, ONLY :                                                     &
        mype

    USE PrintStatus_mod, ONLY :                                                &
        PrintStatus,                                                           &
        PrStatus_Diag                   ! =4 Extra Diagnostic output

    USE lookup_addresses

    IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
    TYPE (List_Marker)               :: List_Header

! Comdecks
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
! contains LBLREC (amongst others)

! Local constants
    CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_GRIB_Dest_List'

! Local variables

    TYPE (Grib_Record), POINTER      :: Current

    CHARACTER (LEN=80)               :: Cmessage(2)   ! used for EReport
    INTEGER                          :: ErrorStatus   ! used for EReport

    INTEGER                          :: cnter

!=======================================================================
!  Routine Code Start :
!=======================================================================

! Check list has _some_ members.
    IF (ASSOCIATED (List_Header % Begin )) THEN

      cnter = 0

  ! Lists with more than 1 member
      IF ( List_Header % LstCount  > 1 ) THEN
        Current => List_Header % Begin

        DO While (ASSOCIATED (Current % Next))
      ! Deallocate array attached to pointer
          IF (ASSOCIATED(Current % VertCoords)) THEN
            DEALLOCATE(Current % VertCoords)
          END IF
          Current => Current % Next
          DEALLOCATE (Current % Prev)
          cnter = cnter + 1
        END DO
      END IF

  ! do for all lists
      NULLIFY (Current)
      DEALLOCATE (List_Header % End)

      IF ( cnter < List_Header % LstCount -1 ) THEN
        WRITE(Cmessage(1),'(A,I2)') 'Destroyed less list entries than I&
        & thought existed'
        ErrorStatus = 10
        CALL EReport( RoutineName, ErrorStatus, Cmessage(1) )
      ELSE
        IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
          WRITE (6,*) "Destroyed a list containing ",cnter +1, " members."
        END IF
      END IF

    END IF ! list had members

    NULLIFY (List_Header % Begin)
    NULLIFY (List_Header % End)
    List_Header % LstCount = 0

    RETURN

  END SUBROUTINE Rcf_Grib_Dest_List
END MODULE Rcf_Grib_Dest_List_Mod
