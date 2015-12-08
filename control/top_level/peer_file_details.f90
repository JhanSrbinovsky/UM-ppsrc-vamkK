! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routines to initialise and write to PEER output files which can then
! be analysed with the PEER utility.
!
! Subroutine interface:
      SUBROUTINE PEER_FILE_DETAILS(                                     &
     &  PUNIT,PEER_FILENAME,SUNIT,SUMM_FILENAME                         &
     &  )

! Description:
!   Returns unit number and filename for sending output
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Subroutine Arguments:


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE nlstcall_mod, ONLY : expt_id, &
                               job_id

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

      INTEGER                                                           &
     &  PUNIT                                                           &
                                ! OUT: Unit number for peer file
     &  ,SUNIT                  ! OUT: Unit number for summary file

      CHARACTER(LEN=30)                                                      &
     &  PEER_FILENAME                                                   &
                                ! OUT: Peer filename
     &  ,SUMM_FILENAME          ! OUT: Summary filename

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Create filename for output from this pe
      IF (lhook) CALL dr_hook('PEER_FILE_DETAILS',zhook_in,zhook_handle)
      WRITE(PEER_FILENAME,10)mype
      PEER_FILENAME(1:4)=expt_id
      PEER_FILENAME(5:5)=job_id
      WRITE(SUMM_FILENAME,20)mype
      SUMM_FILENAME(1:4)=expt_id
      SUMM_FILENAME(5:5)=job_id
 10   FORMAT('......peer.pe',I3.3)
 20   FORMAT('......summ.pe',I3.3)

      PUNIT=58
      SUNIT=59

      IF (lhook) CALL dr_hook('PEER_FILE_DETAILS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PEER_FILE_DETAILS
