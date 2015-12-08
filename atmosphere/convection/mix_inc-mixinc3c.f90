! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  SUBROUTINE MIX_INC------------------------------------------------
!LL
!LL  PURPOSE : TO WELL-MIX CONVECTIVE INCREMENTS IN THE BOUNDARY LAYER
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
      SUBROUTINE MIX_INC (NP_FIELD,NPNTS,NCPTS,NLEV,NBL,NWML,           &
                          DTHBYDT,DQBYDT,DUBYDT,DVBYDT,L_TRACER,NTRA,   &
                          DTRABYDT,P_LAYER_BOUNDARIES,P_LAYER_CENTRES,  &
                          INDEX4)

      Use cv_run_mod, Only:                                             &
          l_mom

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
!
!----------------------------------------------------------------------
! VECTOR LENGTH AND LOOP COUNTERS
!----------------------------------------------------------------------
!

      INTEGER NP_FIELD            ! LENGTH OF DATA (ALSO USED TO
                                  ! SPECIFY STARTING POINT OF
                                  ! DATA PASSED IN)
!
      INTEGER NPNTS               ! IN Full VECTOR LENGTH
!
      INTEGER NCPTS               ! IN Vector length
!
      INTEGER NLEV                ! IN Number of model levels
!
      INTEGER NBL                 ! IN Number of model layers
!                                 !    potentially in the boundary layer
!
      INTEGER NTRA                ! IN Number of tracer fields
!
      INTEGER I,K,ITRA            ! Loop counters
!
      LOGICAL L_TRACER            ! IN Logical switch for inclusion
!                                 !    of convective mixing of tracers
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      INTEGER INDEX4(NPNTS)
      INTEGER NWML(NPNTS)         ! IN Number of model layers for which
!                                 !    increments are to be well-mixed.
!
      REAL P_LAYER_BOUNDARIES(NP_FIELD,0:NLEV)
      REAL P_LAYER_CENTRES(NP_FIELD,0:NLEV)
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!----------------------------------------------------------------------
!
      REAL DTHBYDT(NP_FIELD,NLEV) ! INOUT
!                                 ! IN  Increment to potential
!                                 !     temperature due to convection
!                                 ! OUT Modified increment
!
      REAL DQBYDT(NP_FIELD,NLEV)  ! INOUT
!                                 ! IN  Increment to mixing ratio
!                                 !     due to convection
!                                 ! OUT Modified increment
!
      REAL DUBYDT(NP_FIELD,NBL)   ! INOUT
!                                 ! IN  Increment to x-component of
!                                 !     wind due to convection
!                                 ! OUT Modified increment
!
      REAL DVBYDT(NP_FIELD,NBL)   ! INOUT
!                                 ! IN  Increment to y-component of
!                                 !     wind due to convection
!                                 ! OUT Modified increment
!
      REAL DTRABYDT(NPNTS,NLEV,                                         &
                                ! INOUT
     &              NTRA)         ! IN  Increments to tracers
!                                 !     due to convection
!                                 ! OUT Modified increment
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE LOCALLY DEFINED
!----------------------------------------------------------------------
!
      REAL THSUM(NCPTS)           ! SUMMATION OF INCREMENTS TO
                                  ! POTENTIAL TEMPERATURE DUE TO
                                  ! CONVECTION IN THE VERTICAL,
                                  ! WEIGHTED ACCORDING TO MASS.
!
      REAL QSUM(NCPTS)            ! SUMMATION OF INCREMENTS TO
                                  ! MIXING RATIO DUE TO CONVECTION
                                  ! IN THE VERTICAL, WEIGHTED
                                  ! ACCORDING TO MASS.
!
      REAL USUM(NCPTS)            ! SUMMATION OF INCREMENTS TO
                                  ! X_COMPONENT OF WIND DUE TO
                                  ! CONVECTION IN THE VERTICAL,
                                  ! WEIGHTED ACCORDING TO MASS.
!
      REAL VSUM(NCPTS)            ! SUMMATION OF INCREMENTS TO
                                  ! Y_COMPONENT OF WIND DUE TO
                                  ! CONVECTION IN THE VERTICAL,
                                  ! WEIGHTED ACCORDING TO MASS.
!
      REAL TRASUM(NCPTS)          ! SUMMATION OF INCREMENTS TO
                                  ! TRACERS DUE TO
                                  ! CONVECTION IN THE VERTICAL,
                                  ! WEIGHTED ACCORDING TO MASS.
!
      REAL DELPSUM(NCPTS)         ! Summation of model layer thicknesses
                                  ! with height. (Pa)
      REAL DELPSUM_UV(NCPTS)      ! Summation of model layer thicknesses
                                  ! with height. (Pa)
!
      REAL DELPK(NCPTS,NBL)       ! DIFFERENCE IN PRESSURE ACROSS A
                                  ! LAYER (Pa)
      REAL DELPK_UV(NCPTS,NBL)    ! DIFFERENCE IN PRESSURE ACROSS A
                                  ! LAYER (Pa)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!
!*----------------------------------------------------------------------
!
!L----------------------------------------------------------------------
!L  Sum up mass weighted rates of change of variables
!L  and layer thicknesses for model layers 1 to NWML.
!L----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('MIX_INC',zhook_in,zhook_handle)
      DO I = 1,NCPTS
        DELPK(I,1) = P_LAYER_BOUNDARIES(INDEX4(I),0) -                  &
     &               P_LAYER_BOUNDARIES(INDEX4(I),1)
        THSUM(I) = DTHBYDT(INDEX4(I),1) * DELPK(I,1)
        QSUM(I) = DQBYDT(INDEX4(I),1) * DELPK(I,1)
        DELPSUM(I) = DELPK(I,1)
      ENDDO
!
      DO K = 2,NBL
        DO I = 1,NCPTS
          IF (K  <=  NWML(INDEX4(I))) THEN
            DELPK(I,K) = P_LAYER_BOUNDARIES(INDEX4(I),k-1) -            &
     &                   P_LAYER_BOUNDARIES(INDEX4(I),k)
            THSUM(I) =                                                  &
     &        THSUM(I) + DTHBYDT(INDEX4(I),K) * DELPK(I,K)
            QSUM(I) =                                                   &
     &         QSUM(I) + DQBYDT(INDEX4(I),K) * DELPK(I,K)
            DELPSUM(I) = DELPSUM(I) + DELPK(I,K)
          ENDIF
        ENDDO
      ENDDO
!
!L----------------------------------------------------------------------
!L Reset potential temperature and humidity increments in layers 1 to
!L NWML to mean values over these layers.
!L----------------------------------------------------------------------
!
      DO K = 1,NBL
        DO I = 1,NCPTS
          IF (K  <=  NWML(INDEX4(I))) THEN
            DTHBYDT(INDEX4(I),K) = THSUM(I) / DELPSUM(I)
            DQBYDT(INDEX4(I),K) = QSUM(I) / DELPSUM(I)
          ENDIF
        ENDDO
      ENDDO
!
      IF (L_MOM) THEN
!
!L----------------------------------------------------------------------
!L  Sum up mass weighted rates of change of variables
!L  for model layers 1 to NWML.
!L----------------------------------------------------------------------
!
        DO I = 1,NCPTS
          DELPK_UV(I,1) = P_LAYER_CENTRES(INDEX4(I),0) -                &
     &                    P_LAYER_CENTRES(INDEX4(I),1)
          DELPSUM_UV(I) = DELPK_UV(I,1)
          USUM(I) = DUBYDT(INDEX4(I),1) * DELPK_UV(I,1)
          VSUM(I) = DVBYDT(INDEX4(I),1) * DELPK_UV(I,1)
        ENDDO
!
        DO K = 2,NBL
          DO I = 1,NCPTS
            IF (K  <=  NWML(INDEX4(I))) THEN
              DELPK_UV(I,K) = P_LAYER_CENTRES(INDEX4(I),K-1) -          &
     &                        P_LAYER_CENTRES(INDEX4(I),K)
              DELPSUM_UV(I) = DELPSUM_UV(I) + DELPK_UV(I,K)
              USUM(I) = USUM(I) + DUBYDT(INDEX4(I),K) * DELPK_UV(I,K)
              VSUM(I) = VSUM(I) + DVBYDT(INDEX4(I),K) * DELPK_UV(I,K)
            ENDIF
          ENDDO
        ENDDO
!
!L----------------------------------------------------------------------
!L Reset wind component increments in layers 1 to NWML to mean values
!L over these layers.
!L----------------------------------------------------------------------
!
        DO K = 1,NBL
          DO I = 1,NCPTS
            IF (K  <=  NWML(INDEX4(I))) THEN
              DUBYDT(INDEX4(I),K) = USUM(I) / DELPSUM_UV(I)
              DVBYDT(INDEX4(I),K) = VSUM(I) / DELPSUM_UV(I)
            ENDIF
          ENDDO
        ENDDO
!
      ENDIF  ! L_MOM
!
      IF (L_TRACER) THEN
!
!L----------------------------------------------------------------------
!L  Sum up mass weighted rates of change of each tracer in turn
!L  for model layers 1 to NWML.
!L----------------------------------------------------------------------
!
        DO ITRA = 1,NTRA
          DO I = 1,NCPTS
            TRASUM(I) = DTRABYDT(INDEX4(I),1,ITRA) * DELPK(I,1)
          ENDDO
!
          DO K = 2,NBL
            DO I = 1,NCPTS
              IF (K  <=  NWML(INDEX4(I))) THEN
                TRASUM(I) = TRASUM(I) +                                 &
     &                        DTRABYDT(INDEX4(I),K,ITRA) * DELPK(I,K)
              ENDIF
            ENDDO
          ENDDO
!
!L----------------------------------------------------------------------
!L Reset tracer increments in layers 1 to NWML to mean values
!L over these layers.
!L----------------------------------------------------------------------
!
          DO K = 1,NBL
            DO I = 1,NCPTS
              IF (K  <=  NWML(INDEX4(I))) THEN
                DTRABYDT(INDEX4(I),K,ITRA) = TRASUM(I) / DELPSUM(I)
              ENDIF
            ENDDO
          ENDDO
        ENDDO  ! Loop over tracers
      ENDIF  ! L_TRACER
!
      IF (lhook) CALL dr_hook('MIX_INC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE MIX_INC
