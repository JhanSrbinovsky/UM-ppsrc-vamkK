! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Initialises accumulated carbon fluxes to zero if new calling period
!
! Subroutine Interface:
      SUBROUTINE INIT_ACC(LAND_PTS,                                     &
     &       NPP_PFT_ACC,G_LEAF_PHEN_PFT_ACC,                           &
     &       RESP_W_PFT_ACC,RESP_S_ACC,ICODE,CMESSAGE)


  USE nstypes
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!   Resets accumulation prognostics to zero if a new TRIFFID calling
!   period is starting.  This routine is needed when starting an NRUN
!   from an initial dump created in either of the following situations:
!
!   i)  Initial dump created from a non-TRIFFID run
!
!   ii) Initial dump created in a TRIFFID run mid-way through a TRIFFID
!       calling period.  The NRUN may re-start at the same point within
!       this calling period and continue with the accumulation already
!       part-completed in this dump; in this case this routine will not
!       be used.  Alternatively, the NRUN may start a new calling
!       period, in which case the accumulation must begin; this routine
!       allows this by re-setting the relevant prognostics to zero.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Vegetation
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Arguments


      INTEGER                                                           &
     & LAND_PTS                            ! IN number of land points


      REAL                                                              &
     & NPP_PFT_ACC(LAND_PTS,NPFT)                                       &
                                         !INOUT Accumulated NPP on
!                                        !      Plant Functional Types
     &,G_LEAF_PHEN_PFT_ACC(LAND_PTS,NPFT)                               &
                                         !INOUT Accum. phenological
!                                        !      leaf turnover rate PFTs
     &,RESP_W_PFT_ACC(LAND_PTS,NPFT)                                    &
                                         !INOUT Accumulated wood
!                                        !      respiration on PFTs
     &,RESP_S_ACC(LAND_PTS,4)         !INOUT


      INTEGER                                                           &
     & L                                                                &
                               ! Loop counter for land points
     &,N                       ! Loop counter for plant functional types

      INTEGER ICODE            ! Work - Internal return code
      CHARACTER(LEN=80) CMESSAGE    ! Work - Internal error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('INIT_ACC',zhook_in,zhook_handle)
      WRITE (6,*)                                                       &
     & 'INIT_ACC: setting accumulation prognostics to zero'

      DO L=1,LAND_PTS
        DO N=1,NPFT
          NPP_PFT_ACC(L,N) = 0.0
          G_LEAF_PHEN_PFT_ACC(L,N) = 0.0
          RESP_W_PFT_ACC(L,N) = 0.0
        ENDDO
      ENDDO
      DO N=1,4
        DO L=1,LAND_PTS
          RESP_S_ACC(L,N) = 0.0
        ENDDO
      ENDDO

      IF (lhook) CALL dr_hook('INIT_ACC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INIT_ACC
