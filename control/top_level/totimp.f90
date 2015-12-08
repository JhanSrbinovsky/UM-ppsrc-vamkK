! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Convert multiple of specified time period to no. of timesteps
! Function Interface:

      INTEGER FUNCTION TOTIMP(PERIOD,UNIT,MDL)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim, dumpfreqim
      USE submodel_mod, ONLY: n_internal_model_max
      IMPLICIT NONE

! Description: Get no. of timesteps from time period information taken
!  from STASH profiles to convert to no. of timesteps required to
!  generate STASH list contents for controlling diagnostic output times.
!
! Method: Simple conversion of specified time periods. Illegal
!  combinations are returned as -999 to be trapped by calling routine.
!
! Note: Negative periods are allowed since the UMUI allows time profiles 
! with a start and end date indicating when the diagnostic is to be 
! written out: Model runs stopped by the operator and restarted may result 
! in a start date for the diagnostic which is earlier than the model basis 
! time.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!

! Function arguments:
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN)         :: PERIOD ! multiples of time period
      CHARACTER(LEN=2), INTENT(IN):: UNIT   ! descriptor for time period
      INTEGER,  INTENT(IN)        :: MDL    ! internal model

! Local scalars:
      INTEGER, PARAMETER::      TOTIMP_ERROR=-999 ! error marker
      REAL FAC

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header ----------------------------------------------------

        IF (lhook) CALL dr_hook('TOTIMP',zhook_in,zhook_handle)

        FAC= REAL(steps_per_periodim(MDL))/REAL(secs_per_periodim(MDL))

        IF      (UNIT == 'T ') THEN       ! timesteps
          TOTIMP=PERIOD
        ELSE IF (UNIT == 'H ') THEN       ! hours
! Convert (negative or positive) hours into timesteps using nearest integer. 
          TOTIMP=NINT(PERIOD*FAC*3600.0)

        ELSE IF (UNIT == 'DA') THEN       ! days
! Convert (negative or positive) days into timesteps using nearest integer. 
          TOTIMP=NINT(PERIOD*FAC*3600.0*24.0)

        ELSE IF (UNIT == 'DU') THEN       ! dump periods
          IF (dumpfreqim(MDL) == 0) THEN
            WRITE(6,*)'TOTIMP:IRREGULAR DUMPS FOR DUMP FREQUENCY'
            TOTIMP= TOTIMP_ERROR
          ELSE
            TOTIMP=PERIOD*dumpfreqim(MDL)
          END IF
        ELSE                              ! illegal unit
          TOTIMP= TOTIMP_ERROR
          WRITE(6,*)'TOTIMP: UNEXPECTED TIME UNIT=',UNIT
        END IF
! Note: the special case of period  ==  -1 (indefinite) is trapped by
! the calling routine, otherwise it would be necessary to include lines
! ELSE IF(PERIOD  ==  -1) THEN TOTIMP= -1 here.


      IF (lhook) CALL dr_hook('TOTIMP',zhook_out,zhook_handle)
      RETURN

      END FUNCTION TOTIMP
