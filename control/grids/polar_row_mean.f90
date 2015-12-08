! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_Row_Mean
! Purpose:
!         Resets polar values to mean over row to remove any
!         deviations due to rounding error.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
MODULE polar_row_mean_mod

IMPLICIT NONE

CONTAINS



Subroutine polar_row_mean(                                             & 
                           field,                                      &
                           ini_start,ini_end,                          &
                           inj_start,inj_end,                          &
                           ink_start,ink_end,                          &
                           global_row_length,                          &
                           n_proc, n_procy, proc_row_group,            &
                           at_extremity)

USE global_2d_sums_mod, ONLY : global_2d_sums
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
IMPLICIT NONE



INTEGER, INTENT(IN)  ::  &
  ini_start,ini_end      & ! in row_length array bounds
, inj_start,inj_end      & ! in rows array bounds
, ink_start,ink_end      & ! in levels array bounds
, n_proc                 & ! Number of procs
, n_procy                & ! Number of procs N/S
, global_row_length      & ! number of points per global row
, proc_row_group           ! Group id for processors on the same row

LOGICAL  ::              &
  at_extremity(4)          ! Indicates if this processor is at north,
                           ! south, east or west of the processor grid


REAL, INTENT(INOUT) ::                            &
field (ini_start:ini_end, inj_start:inj_end,      &
                          ink_start:ink_end )

! Local Variables.

INTEGER  ::              &
 i,k                     & ! Loop indices
,info, levels, rlength


REAL     ::                                       &
polar_row (ini_start:ini_end,ink_start:ink_end),  &
sum       (ink_start:ink_end)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('POLAR_ROW_MEAN',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1.   Reset polar values to mean over the row.
! ----------------------------------------------------------------------

rlength= ini_end-ini_start+1
levels=ink_end-ink_start+1

! variables on all model levels.

  IF (n_procy  >   1) THEN 

     IF(at_extremity(PNorth)) THEN
          DO k = ink_start,ink_end
            DO i = ini_start,ini_end
              polar_row(i,k) = field(i,inj_end,k)
            END DO
          END DO

     ELSE IF(at_extremity(PSouth)) THEN
          DO k = ink_start,ink_end
            DO i = ini_start,ini_end
              polar_row(i,k) = field(i,inj_start,k)
            END DO
          END DO

     ELSE
          DO k = ink_start,ink_end
            DO i = ini_start,ini_end
              polar_row(i,k) = 0.
            END DO
          END DO
          
     END IF

     
     CALL global_2d_sums(polar_row, rlength, 1, 0, 0, levels, &
                         sum, proc_row_group)

     DO i = ink_start,ink_end 
        sum(i) = sum(i) / global_row_length
     END DO

     IF(at_extremity(PNorth)) THEN
          DO k = ink_start,ink_end 
            DO i = ini_start,ini_end
              field(i,inj_end,k) = sum(k)
            END DO
          END DO
         

     ELSE IF(at_extremity(PSouth)) THEN

          DO k = ink_start,ink_end 
            DO i = ini_start,ini_end
              field(i,inj_start,k) = sum(k)
            END DO
          END DO

     END IF

  ELSE


! One processor N/S so must have both North and South
! Do north pole

     DO k = ink_start,ink_end 
       DO i = ini_start,ini_end
         polar_row(i,k) = field(i,inj_end,k)
       END DO
     END DO
        
     CALL global_2d_sums(polar_row, rlength, 1, 0, 0, levels, &
                         sum, proc_row_group)

     DO i = ink_start,ink_end 
       sum(i) = sum(i) / global_row_length
     END DO

     DO k = ink_start,ink_end 
       DO i = ini_start,ini_end
         field(i,inj_end,k) = sum(k)
       END DO
     END DO



! Do south pole

     DO k = ink_start,ink_end
       DO i = ini_start,ini_end
         polar_row(i,k) = field(i,inj_start,k)
       END DO
     END DO
      

     CALL global_2d_sums(polar_row, rlength, 1, 0, 0, levels, &
                         sum, proc_row_group)

     DO i = ink_start,ink_end
       sum(i) = sum(i) / global_row_length
     END DO

     DO k = ink_start,ink_end
       DO i = ini_start,ini_end
         field(i,inj_start,k) = sum(k)
       END DO
     END DO
        
  END IF

! End of routine
IF (lhook) CALL dr_hook('POLAR_ROW_MEAN',zhook_out,zhook_handle)
RETURN
END SUBROUTINE Polar_Row_Mean
END MODULE polar_row_mean_mod
