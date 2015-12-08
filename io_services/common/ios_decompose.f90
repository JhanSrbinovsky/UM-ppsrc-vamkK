! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Description:
!
! Auxiliary functions to assist in decomposing the IOS server.

MODULE IOS_decompose

  IMPLICIT NONE

CONTAINS

  SUBROUTINE distribute_range( &
      global_start,global_end, &
      local_start ,local_end,  &
      local_bin, bins)

! Method:
! Deterministically distribute a number range 
! global_start...global_end over a number of bins
! and return the range owned by a particular local_bin.
!
! Typically bins may be a number of processors, local_bin will be
! my rank in that set, such that local_start - local_end 
! is the processor local domain.

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: global_start
    INTEGER, INTENT(IN)  :: global_end
    INTEGER, INTENT(IN)  :: bins
    INTEGER, INTENT(IN)  :: local_bin
    INTEGER, INTENT(OUT) :: local_start
    INTEGER, INTENT(OUT) :: local_end
    
    INTEGER              :: base
    INTEGER              :: items
    INTEGER              :: remainder
    INTEGER              :: local_items
    INTEGER              :: num_below

    items=global_end-global_start+1

    base=items/bins
    remainder=items-base*bins

    local_items=base

    IF ( local_bin < remainder ) local_items=local_items+1

    num_below=MIN(remainder,local_bin)

    local_start = local_bin*base+1+num_below
    local_end   = local_start+local_items-1

  END SUBROUTINE distribute_range


  FUNCTION locate_in_range(    &
      global_start,global_end, &
      item, bins)              &
      RESULT (local_bin)

! Method: 
! An inefficient function to determine which bin in a 
! distributed range (as specified by distribute_range)
! owns a particular number in the global range.

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: global_start
    INTEGER, INTENT(IN)  :: global_end
    INTEGER, INTENT(IN)  :: bins
    INTEGER, INTENT(IN)  :: item

    INTEGER              :: local_bin    
    INTEGER              :: local_start
    INTEGER              :: local_end
    INTEGER              :: i

    DO i=0,bins-1

       CALL distribute_range(       &
           global_start,global_end, &
           local_start,local_end,   &
           i, bins)

       IF (item >= local_start .AND. item <= local_end) THEN
         local_bin=i
         RETURN
       END IF

    END DO

  END FUNCTION locate_in_range

END MODULE IOS_decompose
