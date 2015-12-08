! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine GRAVSETT ----------------------------------------------
!
! Purpose: To perform gravitational settlement of tracer particles
!          down to the lowest layer of the model.
!          This version allows tracers to fall through 1 or 2 layers.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!
! Code description:
!  Language: FORTRAN90
!  Programming standard: UMDP 3 
!
! System components covered:
!
! System task:
!
!Documentation: Not yet available
!
!-----------------------------------------------------------------------
SUBROUTINE gravsett(                                                    &
          diam,rhop,p_layer_centres,p_layer_boundaries,t,tracfld,drydep)

  USE atm_fields_bounds_mod, ONLY: pdims, tdims
  USE atmos_constants_mod, ONLY: r
  USE timestep_mod, ONLY: timestep
  USE earth_constants_mod, ONLY: g
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

  REAL, INTENT(IN) :: diam                !IN tracer particle diameter
  REAL, INTENT(IN) :: rhop                !IN tracer particle density
  REAL, INTENT(IN) :: p_layer_centres(pdims%i_start:pdims%i_end,        &
                          pdims%j_start:pdims%j_end,0:pdims%k_end)!IN
  REAL, INTENT(IN) :: p_layer_boundaries(pdims%i_start:pdims%i_end,     &
                             pdims%j_start:pdims%j_end,0:pdims%k_end)!IN
  REAL, INTENT(IN) :: t(tdims%i_start:tdims%i_end,                      &
                 tdims%j_start:tdims%j_end,1:tdims%k_end)!IN temperature
 
  REAL, INTENT(INOUT) :: tracfld(pdims%i_start:pdims%i_end,             &
                         pdims%j_start:pdims%j_end,                     &
                         pdims%k_start:pdims%k_end) !IN/OUT tracer field

  REAL, INTENT(INOUT) :: drydep(pdims%i_start:pdims%i_end,              &
                                pdims%j_start:pdims%j_end)
                                   !IN/OUT dep flux from
                                   !       layer2(kg m-2 s-1)

! Local variables

  INTEGER :: k                  !LOC loop counter for levels
  INTEGER :: j                  !LOC loop counter for points
  INTEGER :: i                  !LOC loop counter for points

  REAL :: vrhoctimestep(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !  v*rho*tracer*deltat @lev
  REAL :: rhok2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !  rho(lev+2)
  REAL :: rhok1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !  rho(lev+1)
  REAL :: rhok(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! rho(lev)
  REAL :: dzk(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! thickness of layer lev
  REAL :: dzk1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! thickness of layer lev+1
  REAL :: dzk2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! thickness of layer lev+2
  REAL :: vstokes(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  pdims%k_start:pdims%k_end)!deposition velocity
!                                         (vstokes corrected)
  REAL :: massout2k2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !flux falling 2 levs from lev k+2
  REAL :: massout1k2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !flux falling 1 levs from lev k+2
  REAL :: massout2k1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !flux falling 2 levs from lev k+1
  REAL :: massout1k1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !flux falling 1 levs from lev k+1
  REAL :: massout2k(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !flux falling 2 levs from lev k
  REAL :: massout1k(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  !flux falling 1 levs from lev k
  REAL :: dummy1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                 pdims%k_start:pdims%k_end) !
  REAL :: dummy2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,   &
                 pdims%k_start:pdims%k_end) !

  INTEGER :: jj, j_block ! omp blocking variables

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('GRAVSETT',zhook_in,zhook_handle)

  j_block = 4

! Calculate settlement velocity

!      CALL VGRAV(PFIELD,NLEVS,DIAM,RHOP,PSTAR,AK,BK,T,V,DUMMY1,DUMMY2,
!     &           FIRST_POINT,LAST_POINT)

! DEPENDS ON: vgrav
  CALL vgrav(                                                           &
   pdims%k_end,diam,rhop,                                               &
   p_layer_centres(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,pdims%k_start:pdims%k_end),&
   t,vstokes,dummy1,dummy2)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,jj)

! Calculate new tracer mixing ratios

! Initialise deposition flux to zero
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        drydep(i,j)=0.
      END DO
    END DO
!$OMP END DO

! Level 1 (K at start of loop)

!  Do not allow tracer to settle directly from level 1
!  but treat this in parallel with tracer mixing in bl_trmix_dd.

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rhok(i,j)=p_layer_centres(i,j,1)/(r*t(i,j,1))
        dzk(i,j)=(p_layer_boundaries(i,j,0)-p_layer_boundaries(i,j,1))  &
                /(rhok(i,j)*g)
        massout2k(i,j)=0.
        massout1k(i,j)=0.
      END DO               !I
    END DO                !J
!$OMP END DO


! Level 2 (K+1 at start of loop)
!   NB  deposit tracer direct to ground from lev 2 if V high enough

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      rhok1(i,j)=p_layer_centres(i,j,2)/(r*t(i,j,2))
      dzk1(i,j)=(p_layer_boundaries(i,j,1)-                             &
       p_layer_boundaries(i,j,2))/(rhok1(i,j)*g)

!   check for deposition :
      IF (vstokes(i,j,2)*timestep  >   dzk(i,j)) THEN
!       some tracer deposited onto ground

        IF (vstokes(i,j,2)*timestep  >   dzk1(i,j)+dzk(i,j)) THEN
!         all deposited to ground
          massout2k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*dzk1(i,j)
          massout1k1(i,j)=0.
        ELSE IF ( vstokes(i,j,2)*timestep  >   dzk1(i,j) ) THEN
!         some deposited to ground, some to layer 1

          massout2k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*                    &
           (vstokes(i,j,2)*timestep-dzk(i,j))
          massout1k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*                    &
           dzk1(i,j)-massout2k1(i,j)
        ELSE
!         some deposited to ground, some to layer1, some left in layer2
          massout2k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*                    &
           (vstokes(i,j,2)*timestep-dzk(i,j))
          massout1k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*dzk(i,j)
        END IF

        drydep(i,j)=drydep(i,j)                                         &
                   +massout2k1(i,j)/(timestep)

      ELSE
!         only falls into layer 1
        massout2k1(i,j)=0.
        IF ( vstokes(i,j,2)*timestep  >   dzk1(i,j)) THEN
!           all falls into layer 1
          massout1k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*dzk1(i,j)
        ELSE
!           some to layer 1 , some left in layer2
          massout1k1(i,j)=rhok1(i,j)*tracfld(i,j,2)*vstokes(i,j,2)      &
           *timestep
        END IF

      END IF

    END DO                !END I LOOP
  END DO                 !END J LOOP
!$OMP END DO



! Main loop through levels, from bottom up
!$OMP DO SCHEDULE(STATIC)
  DO jj=pdims%j_start, pdims%j_end, j_block
    DO k = pdims%k_start, pdims%k_end-2

      DO j = jj, MIN(jj+j_block-1,pdims%j_end)
        DO i = pdims%i_start, pdims%i_end

          rhok2(i,j)=p_layer_centres(i,j,k+2)/(r*t(i,j,k+2))
          dzk2(i,j)=(p_layer_boundaries(i,j,k+1)-                       &
               p_layer_boundaries(i,j,k+2))/(rhok2(i,j)*g)
          
!       Calculate mass of tracer falling between levels

!        limit fall to 2 levs
          IF (vstokes(i,j,k+2)*timestep >  (dzk1(i,j)+dzk(i,j)))        &
               vstokes(i,j,k+2)=(dzk1(i,j)+dzk(i,j))/timestep
          
!          check how far tracer falls:
          IF ( vstokes(i,j,k+2)*timestep  >   dzk1(i,j) ) THEN
!          it falls through more than 1 layer
            IF ( vstokes(i,j,k+2)*timestep  >                           &
                 (dzk2(i,j)+dzk1(i,j)) ) THEN
!             all into layer k
              massout2k2(i,j)=rhok2(i,j)*tracfld(i,j,k+2)*dzk2(i,j)
              massout1k2(i,j)=0.
            ELSE IF ( vstokes(i,j,k+2)*timestep  >   dzk2(i,j) ) THEN
!            some into k+1, some into k
              massout2k2(i,j)=                                          &
                   rhok2(i,j)*tracfld(i,j,k+2)*(vstokes(i,j,k+2)*       &
                   timestep-dzk1(i,j))
              massout1k2(i,j)=                                          &
                   rhok2(i,j)*tracfld(i,j,k+2)*dzk2(i,j)-massout2k2(i,j)
            ELSE
!            some left in k+2, some into k+1, some into k
              massout2k2(i,j)=                                          &
                   rhok2(i,j)*tracfld(i,j,k+2)*(vstokes(i,j,k+2)*       &
                   timestep-dzk1(i,j))
              massout1k2(i,j)=rhok2(i,j)*tracfld(i,j,k+2)*dzk1(i,j)
            END IF
            
          ELSE
!          falls no more than 1 layer
            massout2k2(i,j)=0.
            IF (vstokes(i,j,k+2)*timestep  >   dzk2(i,j)) THEN
!            all falls into layer k+1
              massout1k2(i,j)=rhok2(i,j)*tracfld(i,j,k+2)*dzk2(i,j)
            ELSE
!            some falls into k+1, some left in k+2
              massout1k2(i,j)=rhok2(i,j)*tracfld(i,j,k+2)*              &
                   vstokes(i,j,k+2)*timestep
            END IF
            
          END IF

! Update tracer field

          tracfld(i,j,k)=tracfld(i,j,k)+(massout2k2(i,j)+               &
               massout1k1(i,j)-massout2k(i,j)-massout1k(i,j))/          &
               (rhok(i,j)*dzk(i,j))
          
! Put k+2 vals in k+1's & k+1's in k's
          massout1k(i,j)=massout1k1(i,j)
          massout1k1(i,j)=massout1k2(i,j)
          massout2k(i,j)=massout2k1(i,j)
          massout2k1(i,j)=massout2k2(i,j)
          dzk(i,j)=dzk1(i,j)
          dzk1(i,j)=dzk2(i,j)
          rhok(i,j)=rhok1(i,j)
          rhok1(i,j)=rhok2(i,j)
          
        END DO           !END I LOOP
      END DO            !END J LOOP

    END DO              !END K LOOP
  END DO                !END JJ LOOP
!$OMP END DO

! Top 2 levels

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      tracfld(i,j,pdims%k_end-1)=tracfld(i,j,pdims%k_end-1)+            &
       (massout1k1(i,j)-massout2k(i,j)-massout1k(i,j))/                 &
       (rhok(i,j)*dzk(i,j))
      tracfld(i,j,pdims%k_end)=tracfld(i,j,pdims%k_end)-                &
       (massout2k1(i,j)+massout1k1(i,j))/                               &
                  (rhok1(i,j)*dzk1(i,j))

    END DO       !I
  END DO        !J
!$OMP END DO

!$OMP END PARALLEL 

  IF (lhook) CALL dr_hook('GRAVSETT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE gravsett
