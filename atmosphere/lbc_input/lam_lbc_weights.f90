
! *****************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*********************************

MODULE lam_lbc_weights_mod
IMPLICIT NONE
REAL, SAVE, ALLOCATABLE :: ns_weights(:,:), ew_weights(:,:)

CONTAINS

   SUBROUTINE init_lbc_wghts(row_length,rows,rimwidth,halo_i, halo_j,          &
                             at_extremity, rimweights                          &
                            )
   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim
   USE UM_ParParams
   IMPLICIT NONE
!
! Description: Precalculate RIM weights for LBC data
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

   INTEGER, INTENT(IN) :: row_length, rows, rimwidth
   INTEGER, INTENT(IN) :: halo_i, halo_j
   REAL,    INTENT(IN) :: rimweights(rimwidth)
   LOGICAL, INTENT(IN) :: at_extremity(4)

! Local variables

   INTEGER :: i, j, rim
   INTEGER :: row_end_pt, row_start_pt, last_row, first_row

   INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
   INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
   REAL(KIND=jprb)               :: zhook_handle

   IF (lhook) CALL dr_hook('EG_BALANCE_LBC_VALUES',zhook_in,zhook_handle)

   IF (at_extremity(PNorth) .OR. at_extremity(PSouth)) THEN

      ALLOCATE( ns_weights(1-halo_i:row_length+halo_i,1-halo_j:rimwidth) )

! Initialise the ns_weights array

! The North/South edge

      DO j = 1-halo_j, 0
         DO i = 1-halo_i, row_length+halo_i
            ns_weights(i,j) = 1.0
         END DO
      END DO

! The Western corner

      IF (at_extremity(PWest)) THEN
         DO j = 1, rimwidth
            DO i = 1-halo_i, 0
               ns_weights(i,j) = 1.0
            END DO
         END DO
      END IF

! The Eastern corner

      IF (at_extremity(PEast)) THEN
         DO j = 1, rimwidth
            DO i = row_length+1, row_length+halo_i
               ns_weights(i,j) = 1.0
            END DO
         END DO
      END IF

! The Western corner (South : m, North : c)

      IF (at_extremity(PWest)) THEN
         DO j = 1, rimwidth
            DO i = 1, rimwidth
               rim = MIN(i,j)

               IF (rim  <=  rimwidth) THEN
                 ns_weights(i,j) = rimweights(rim)
               ELSE
                 ns_weights(i,j) = 0.0
               END IF

            END DO
         END DO
      END IF

! The Eastern corner (South : o, North : e)

      IF (at_extremity(PEast)) THEN
         DO j = 1, rimwidth
            DO i = row_length-rimwidth+1, row_length
               rim   = MIN(row_length+1-i,j)

               IF (rim  <=  rimwidth) THEN
                 ns_weights(i,j) = rimweights(rim)
               ELSE
                 ns_weights(i,j) = 0.0
               END IF

            END DO
         END DO
      END IF

! The bit between the corners (South : n, North : d)

      IF (at_extremity(PEast)) THEN
         row_end_pt = row_length-rimwidth
      ELSE
         row_end_pt = row_length+halo_i
      END IF

      IF (at_extremity(Pwest)) THEN
         row_start_pt = rimwidth+1
      ELSE
         row_start_pt = 1-halo_i
      END IF

      DO j = 1, rimwidth
         DO i = row_start_pt, row_end_pt
            IF (j  <=  rimwidth) THEN
              ns_weights(i,j) = rimweights(j)
            ELSE
              ns_weights(i,j) = 0.0
            END IF
         END DO
      END DO

   END IF

   IF (at_extremity(PWest) .OR. at_extremity(PEast)) THEN

      ALLOCATE( ew_weights(1-halo_i:rimwidth,1-halo_j:rows+halo_j) )

! We will need the ew_weights array

      IF (at_extremity(PSouth)) THEN
         first_row = rimwidth+1
      ELSE
         first_row = 1-halo_j
      END IF

      IF (at_extremity(PNorth)) THEN
         last_row = rows-rimwidth
      ELSE
         last_row = rows+halo_j
      END IF

! First the halo area (West : g , East k)

      DO j = first_row, last_row
         DO i = 1-halo_i, 0
            ew_weights(i,j) = 1.0
         END DO
      END DO

! And now the boundary area (West : h, East : j)

      DO j=first_row, last_row
         DO i = 1, rimwidth
            IF (i <= rimwidth) THEN
               ew_weights(i,j) = rimweights(i)
            ELSE
               ew_weights(i,j) = 0.0
            END IF
         END DO
      END DO

   END IF

   IF (lhook) CALL dr_hook('EG_BALANCE_LBC_VALUES',zhook_out,zhook_handle)

   END SUBROUTINE init_lbc_wghts
END MODULE lam_lbc_weights_mod
