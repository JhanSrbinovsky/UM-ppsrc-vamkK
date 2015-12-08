! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Provides convenience routines for async stash operations

MODULE stwork_aux

  USE IOS_Stash
  USE IOS_Stash_Common
  USE IOS_Constants
  USE IOS_Model_Geometry, ONLY :                                               &
      maxFieldDomain
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  CHARACTER (LEN=132) :: stwork_aux_Message

  ! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS

  SUBROUTINE completedBuffer(s)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: s
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:COMPLETEDBUFFER',                                         &
        zhook_in,zhook_handle)

    ! Release data buffers to the IO server
    CALL ios_stash_dispatch( s )
    ! Send the control message to the IO Server.
    IF (model_rank == 0)                                                       &
        CALL ios_stash_write_pp_data( s )

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:COMPLETEDBUFFER',                                         &
        zhook_out,zhook_handle)

  END SUBROUTINE completedBuffer


  SUBROUTINE flushPendingStash()
    IMPLICIT NONE
    INTEGER         :: i
    INTEGER         :: state
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:FLUSHPENDINGSTASH',                                       &
        zhook_in,zhook_handle)

    IF (isUsingAsyncStash()) THEN
      DO i=1,IOS_AsyncNumSlots
        state=IOS_getSlotState(i)

        IF (slot(i)%bufferType==IOS_Action_StashWritePPData) THEN
          SELECT CASE (state)

          CASE (ios_queue_slot_partfilled)
            CALL completedBuffer(i)
          CASE (ios_queue_slot_initialized)
            IF (IOS_Verbosity >= IOS_PrStatus_Oper )                           &
                WRITE(6,'(A,I3,A)')                                            &
                'stwork_aux: flush ',i,' initialised (resetting)'
            CALL resetBuffer(i)
          CASE (ios_queue_slot_dispatched)

          CASE (ios_queue_slot_unused)

          CASE (ios_queue_slot_ready)

          CASE DEFAULT
            WRITE(stwork_aux_message,'(A,I8)')                                 &
                'Slot :',i,' is in unknown state ',state
            CALL IOS_Ereport( 'stwork_aux:flushPending',                       &
                -10, stwork_aux_Message )
          END SELECT
        END IF
      END DO
    END IF

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:FLUSHPENDINGSTASH',                                       &
        zhook_out,zhook_handle)

  END SUBROUTINE flushPendingStash


  SUBROUTINE flushPendingStashForUnit(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER             :: i
    INTEGER             :: state
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:FLUSHPENDINGSTASHFORUNIT',                                &
        zhook_in,zhook_handle)

    IF (isUsingAsyncStash()) THEN
      DO i=1,IOS_AsyncNumSlots
        IF (unit==IOS_getSlotUnit(i)) THEN
          state=IOS_getSlotState(i)
          SELECT CASE (state)

          CASE (ios_queue_slot_partfilled)
            CALL completedBuffer(i)
          CASE (ios_queue_slot_initialized)
            IF (IOS_Verbosity >= IOS_PrStatus_Oper )                           &
                WRITE(6,'(A,I3,A)')                                            &
                'stwork_aux: flush ',i,                                        &
                ' initialised (resetting)'
            CALL resetBuffer(i)
          CASE (ios_queue_slot_dispatched)

          CASE (ios_queue_slot_unused)

          CASE (ios_queue_slot_ready)

          CASE DEFAULT
            WRITE(stwork_aux_message,'(A,I8)')                                 &
                'Slot :',i,' is in unknown state ',state
            CALL IOS_Ereport( 'stwork_aux:flushPending',                       &
                -10, stwork_aux_Message )
          END SELECT
        END IF
      END DO
    END IF

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:FLUSHPENDINGSTASHFORUNIT',                                &
        zhook_out,zhook_handle)

  END SUBROUTINE flushPendingStashForUnit


  FUNCTION getSlotForNextLevel(unit) RESULT (s)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER             :: s
    INTEGER             :: current_levels_packed
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:GETSLOTFORNEXTLEVEL',                                     &
        zhook_in,zhook_handle)

    s=-1

    IF (isUsingAsyncStash())THEN

      s=IOS_getCurrentSlot(unit)
      IF (s<0) THEN

        CALL ios_stash_next_buffer( unit, s )
        IF (s==-1)THEN
          IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
            WRITE(6,'(A)')                                                     &
                'IOS: Warning: stwork_aux: No new slot returned'
            WRITE(6,'(A,A)')'IOS: Warning: stwork_aux: ',                      &
                '  flushing everything - expensive '
          END IF

          ! next buffer failed, flush stuff and retry...
          CALL flushPendingStash()
          CALL ios_stash_next_buffer( unit, s )

          IF (s==-1)THEN
            WRITE(stwork_aux_message,'(A)')                                    &
                'Failure tring to locate a buffer (code fault)'
            CALL IOS_Ereport( 'stwork_aux:slotfornextlevel',                   &
                10, stwork_aux_Message )
          END IF

        END IF

        IF (IOS_Verbosity>=IOS_prstatus_diag)                                  &
            WRITE(6,'(A,A,I3,A,I3)')'IOS: Info: stwork_aux: ',                 &
            'No current slot for unit, ',unit,                                 &
            ' new slot is ',s

      ELSE
        current_levels_packed=IOS_getLevsInPack(s)
        IF (current_levels_packed >= IOS_AsyncMaxFieldsInPack )THEN
          ! This buffer is full need a new buffer

          ! Dispatch the old one
          CALL completedBuffer(s)

          !Get a new one
          CALL ios_stash_next_buffer( unit, s )

          IF (s==-1)THEN
            WRITE(stwork_aux_message,'(A)')                                    &
                'Failure trying to locate a buffer (code fault)'
            CALL IOS_Ereport( 'stwork_aux:slotfornextlevel',                   &
                10, stwork_aux_Message )
          END IF

          IF (IOS_Verbosity>=IOS_prstatus_diag)                                &
              WRITE(6,'(A,A,I3,A,I3)')'IOS: Info: stwork_aux: ',               &
              'Buffer for unit, ',unit,                                        &
              ' was full, new slot is ',s

        END IF
      END IF
    END IF

    IF (lhook) CALL dr_hook                                                    &
        ('STWORK_AUX:GETSLOTFORNEXTLEVEL',                                     &
        zhook_out,zhook_handle)

  END FUNCTION getSlotForNextLevel

END MODULE stwork_aux
