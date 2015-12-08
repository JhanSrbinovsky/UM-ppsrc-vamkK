! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine INTF_NEW_FILES -----------------------------------------
!
! Purpose: To test whether a new boundary data file needs to be
!           opened
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
      subroutine intf_new_files(first_unit, last_unit, max_n_intf, im,  &
     &    TYPE_LETTER_1, FT_OUTPUT, INTF_FREQ_HR, INTF_FREQ_MN,         &
     &    INTF_FREQ_SC, FT_STEPS, STEP, FT_FIRSTSTEP, INTERFACE_STEPS,  &
     &    LNEWBND )

! Purpose: determines new output interface files for a submodel.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      implicit none
      integer first_unit ! IN first unit to test
      integer last_unit  ! IN last unit to test
      integer max_n_intf ! IN number of interface files
      integer im         ! IN  sub-model identifier
      CHARACTER(LEN=1) TYPE_LETTER_1(20:last_unit) ! IN
      CHARACTER(LEN=1) FT_OUTPUT(20:last_unit)     ! IN
      integer INTF_FREQ_HR(max_n_intf)     ! IN
      integer INTF_FREQ_MN(max_n_intf)     ! IN
      integer INTF_FREQ_SC(max_n_intf)     ! IN
      integer FT_STEPS(20:last_unit)          ! IN
      integer STEP                         ! IN model step no.
      integer FT_FIRSTSTEP(20:last_unit)      ! IN
      integer INTERFACE_STEPS(max_n_intf)  ! IN
      logical LNEWBND(max_n_intf)          ! OUT
!-----------------------------------------------------------------------
!L Declaration of local variables
      integer iunit
!     logical ll_intf_type   ! OUT T => file is an output interface file
      integer jintf          ! boundary file area number
      INTEGER INTF_FREQ_SECS  ! Interface frequency in seconds

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!-----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('INTF_NEW_FILES',zhook_in,zhook_handle)
      do iunit = first_unit, last_unit

       if (type_letter_1(iunit) == 'b') then  !  Boundary file

         IF (FT_OUTPUT(IUNIT) == 'Y') THEN ! Intf. data output?
! DEPENDS ON: intf_area
          call intf_area ( im, iunit, JINTF)

      INTF_FREQ_SECS=3600*INTF_FREQ_HR(JINTF) +                         &
     &  60*INTF_FREQ_MN(JINTF)+INTF_FREQ_SC(JINTF)
      IF (INTF_FREQ_SECS >  0) THEN

            IF (STEP == 0) THEN

              LNEWBND(JINTF) = .TRUE. ! New intf data file required at
                                      ! first entry to ININTF1
            ELSE

              IF (FT_STEPS(IUNIT) == 0) LNEWBND(JINTF) = .FALSE. !False
                                         ! if incomplete single file

              IF (FT_STEPS(IUNIT) >  0) LNEWBND(JINTF) = .NOT.(         &
!               step = first timestep to get boundary data
     &          (STEP-FT_FIRSTSTEP(IUNIT) == 0 .OR.                     &
!               step = timestep to start new file
     &          MOD(STEP-FT_FIRSTSTEP(IUNIT),FT_STEPS(IUNIT)) /= 0)     &
     &          .AND.                                                   &
     &     STEP >  FT_FIRSTSTEP(IUNIT)-INTERFACE_STEPS(JINTF))
!                 ! False if incomplete file in sequence
            END IF  ! STEP

           ENDIF  !  INTF_FREQ_HR
!
          ELSE    !  FT_OUTPUT(IUNIT)
!
!  Possible place for setting switches for reading in interface data
!
          ENDIF   !  FT_OUTPUT
        ENDIF  !  TYPE_LETTER_1

      ENDDO   ! IUNIT

      IF (lhook) CALL dr_hook('INTF_NEW_FILES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE intf_new_files
