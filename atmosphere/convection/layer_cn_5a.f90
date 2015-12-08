! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates layer dependent constants for model layer k

      SUBROUTINE LAYER_CN(k,npnts,nlev                                 &
      ,                   mdet_on, ent_on                              &
      ,                   ntml,ntpar                                   &
      ,                   l_shallow,l_congestus,l_deep                 &
      ,                   bconv,bwk,bwkp1                              &
      ,                   exner_layer_boundaries                       &
      ,                   exner_layer_centres                          &
      ,                   p_layer_boundaries,p_layer_centres           &
      ,                   recip_pstar,entrain_coef,rhum                &
      ,                   zk, zkp12, zkp1                              &
      ,                   thek, qek,qsek, thekp1,qekp1,qsekp1          &
      ,                   thpk,qpk ,qs_cb, ekm14                       &
      ,                   pkp1,delpkp1,exkp1                           &
      ,                   pk,delpk,delpkp12,exk,delexkp1               &
      ,                   delp_uv_k, delp_uv_kp1                       &
      ,                   ekp14,ekp34,amdetk                           &
                         )

      USE atmos_constants_mod, ONLY: cp

      Use cv_run_mod, Only:                                            &
          ent_fac_dp, ent_fac_md, amdet_fac, ent_opt_dp, ent_opt_md,   &
          ent_dp_power, ent_md_power

      USE cv_param_mod, ONLY: ae2

      Use cv_dependent_switch_mod, Only:                               &
          l_var_entrain, l_new_det, l_const_ent, sh_grey
! , l_rh_dep, l_rh_dep2  May want to add to list for future use

      USE earth_constants_mod, ONLY: g
      USE entcoef_mod, ONLY: entcoef

      USE water_constants_mod, ONLY: lc, lf
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
 
!----------------------------------------------------------------------
! Description:
!   Calculates layer dependent constants for layer K i.e.
!            pressure
!            layer thickness
!            entrainment coefficients
!            detrainment coefficients
!
! Method:
!   See Unified Model documentation paper 27.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.0 programming standards.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

! Vector lengths and loop counters

      Integer, intent(in) ::   &
       k                       & ! Present model layer 
      ,npnts                   & ! Vector length
      ,nlev                      ! Number of model levels

! Switches

      Integer, intent(in) ::   &
       mdet_on                 & ! flag for adaptive mixing detrainment
                                 !  on = 1, off = 0
      ,ent_on                    ! flag for adaptive entrainment
                                 !  on = 1, off = 0

      Integer, intent(in) ::   &
       ntml(npnts)             & ! Number of levels in the surface-based
                                 ! turbulently mixed layer
      ,ntpar(npnts)             ! Top of initial parcel ascent
      
! Logical Switches

      Logical, intent(in) ::   &
       l_shallow               & ! indicator all points are shallow
      ,l_congestus             & ! indicator all points are congestus
      ,l_deep                    ! indicator all points are deep

      Logical, intent(in) :: &
       bconv(npnts)       & ! mask for points at which convection takes place
      ,bwk(npnts)         & ! mask for points where condensate is liquid on k
      ,bwkp1(npnts)         ! mask for points where condensate is liquid on k+1

! Field on model levels

      Real,  intent(in) ::                   &
       exner_layer_boundaries(npnts,0:nlev)  &
                                       ! Exner function at layer boundary
                                       ! starting at level k-1/2
      ,exner_layer_centres(npnts,0:nlev)     &
                                       ! Exner function at layer  centre
      ,p_layer_centres(npnts,0:nlev)         &
                                       ! Pressure at layer  centre (Pa)
      ,p_layer_boundaries(npnts,0:nlev)      
                                       ! Pressure at layer  boundary (Pa)
! Field on a level
      Real,  intent(in) :: &
       recip_pstar(npnts)  & ! Reciprocal of pstar array (1/Pa)
      ,entrain_coef(npnts) & ! entrainment coefficients
      ,rhum(npnts)         & ! Relative humidity at level K
      ,zk(npnts)           & ! height on k
      ,zkp12(npnts)        & ! height on k+1/2
      ,zkp1(npnts)         & ! height on k+1
      ,thek(npnts)         & ! theta for environment on k
      ,qek(npnts)          & ! q for environment on k
      ,qsek(npnts)         & ! q sat for environment on k
      ,thekp1(npnts)       & ! theta for environment on k+1
      ,qekp1(npnts)        & ! q for environment on k+1
      ,qsekp1(npnts)       & ! q sat for environment on k+1
      ,thpk(npnts)         & ! parcel theta on k
      ,qpk(npnts)          & ! parcel q on k
      ,qs_cb(npnts)        & ! qsat at cloud base 
      ,ekm14(npnts)          ! ek14 from previous pass i.e. ek for k-1+1/4 

! Information for layer k+1

      Real,  intent(inout) ::  &
       pkp1(npnts)             & ! pressure at layer K+1 (Pa) 
      ,delpkp1(npnts)          & ! thickness of layer K+1 (Pa) 
      ,exkp1(npnts)              ! Exner function at level K+1

! Information for layer k

      Real,  intent(out) ::    &
       pk(npnts)               & ! pressure at layer K (Pa) 
      ,delpk(npnts)            & ! thickness of layer K (Pa) 
      ,delpkp12(npnts)         & ! thickness between layer K & K+1 (Pa) 
      ,exk(npnts)              & ! Exner function at level K
      ,delexkp1(npnts)         & ! Difference in Exner function between
                                 ! K+3/2 and K+1/2
      ,delp_uv_k(npnts)        & ! thickness of uv layer K (Pa)
      ,delp_uv_kp1(npnts)      & ! thickness of uv layer K+1 (Pa)
      ,ekp14(npnts)            & ! entrainment coefficient at level k+1/4
                                 ! multiplied by appropriate layer thickness
      ,ekp34(npnts)            & ! entrainment coefficient at level k+3/4
                                 ! multiplied by appropriate layer thickness
      ,amdetk(npnts)             ! mixing detrainment coefficient at level k
                                 ! multiplied by appropriate layer thickness

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------

      Integer :: I             ! Counter for do loop
!
      Real ::                &
       Aekp14,Aekp34         & ! Used in calculation of entrainment rate
      ,aekp14_md,aekp34_md   & ! Used in calculation of entrainment rate (mid)
      ,RENT                  & ! parameter controlling fractional entrainment
      ,EL                    & ! latent heat of condensation(/freezing)
      ,gzk                   & ! g*z(k)
      ,cpexnerk              & ! cp*exner(k)
      ,gzkp1                 & ! g*z(k+1)
      ,cpexnerkp1_th         & ! cp*exner(k+1)*theta(k)
      ,hsat_ek(npnts)        &
      ,hsat_ekp1(npnts)      & ! saturated moist static energies
      ,hsat_ekp12(npnts)     &
      ,h_pk(npnts)           & ! moist static energy
      ,h_pkp12(npnts)        &
      ,h_ek(npnts)           &
      ,h_ekp1(npnts)         &
      ,h_ekp12(npnts)        &
      ,delta_hsat_env(npnts) & ! hsat_ekp1 - hsat_ek
      ,dhpbydp(npnts)        &
      ,dhebydp(npnts)        &
      ,epsilonp14(npnts)     & ! fractional entrainment coef calc from hsat
      ,epsilonp34(npnts)     & ! fractional entrainment coef calc from hsat
      ,ekp14_ad(npnts)       & ! EPS * layer thickness for 1/4 layer
      ,ekp34_ad(npnts)       & ! EPS * layer thickness for 3/4 layer
      ,delp_cld(npnts)       & ! thickness of cloud layer (Pa)
      ,delpkp14(npnts)       & ! thickness of layer for k+1/4 (Pa)
      ,delpkp34(npnts)       & ! thickness of layer for k+3/4 (Pa)
      ,amdet_func(npnts)       ! detrainment function of RH
      
      Real ::                &
       ratio14               &
      ,ratio34               &
      ,factor                & 
      ,factor1               & 
      ,factor2 


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('LAYER_CN',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Set constant ae used in calculation of entrainment and detrainment
! rates depending upon level.
!----------------------------------------------------------------------
!
      aekp14 = ent_fac_dp*ae2
      aekp34 = ent_fac_dp*ae2
      aekp14_md = ent_fac_md*ae2
      aekp34_md = ent_fac_md*ae2
 

!      write(6,*) ' layer_cn', l_deep,l_shallow
!      write(6,*) ' entrain_coef '
!      write(6,*) (entrain_coef(i),i=1,npnts)

!---------------------------------------------------------------------
! Initialise various arrays holding layer thicknesses etc.
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Calculate pk and delpk - if K = 1 (lowest model layer) then
! values for previous pass through routine at (k-1)+1 are taken.
! Calculate Exner functions at mid-layers K and K+1, and difference
! of exner function across layer K
!---------------------------------------------------------------------

      IF (K == 2) THEN
        DO i=1,npnts
          pk(i)    = p_layer_centres(i,k)
          delpk(i) = p_layer_boundaries(i,k-1) - p_layer_boundaries(i,k)
          exk(i)   = exner_layer_centres(i,k)
        END DO 
      ELSE
        DO i=1,npnts
          pk(i)    = pkp1(i)
          delpk(i) = delpkp1(i)
          exk(i)   = exkp1(i)
        END DO 
      END IF

!---------------------------------------------------------------------
! Calculate pkp1, delpkp1 and delpk+1/2
!---------------------------------------------------------------------
      DO i=1,npnts
        pkp1(i)     = p_layer_centres(i,k+1)
        delpkp1(i)  = p_layer_boundaries(i,k) -                      &
                                       p_layer_boundaries(i,k+1)
        delpkp12(i) = pk(i) - pkp1(i)
        exkp1(i)    = exner_layer_centres(i,k+1)
        delexkp1(i) = exner_layer_boundaries(i,k)-                   &
                                     exner_layer_boundaries(i,k+1)

!---------------------------------------------------------------------
! Calculate delp_uv_k and delp_uv_kp1
! NB: Ensure positive by taking - of difference
!---------------------------------------------------------------------
!
        delp_uv_k(i)   = p_layer_centres(i,k-1) - p_layer_centres(i,k)
        delp_uv_kp1(i) = p_layer_centres(i,k) - p_layer_centres(i,k+1)

! Thickness of layer for k+1/4
        delpkp14(i)    = pk(i) - p_layer_boundaries(i,k)

! Thickness of layer for k+3/4
        delpkp34(i)    = p_layer_boundaries(i,k) - pkp1(i)
      END DO


! For use later in this routine only 
! Cloud thickness  - only used for shallow convection

      IF (l_shallow) THEN
        DO i=1,npnts
          delp_cld(i)    = p_layer_boundaries(i,ntml(i))               &
                                   -p_layer_boundaries(i,ntpar(i))
        END DO
      END IF
!----------------------------------------------------------------------
! Adaptive entrainment option 
! UM documentation paper 27 Section 7.3
!----------------------------------------------------------------------

      IF (ent_on  ==  1 ) THEN

        rent=0.5

! Initialise everything in sight!,  except ekm14, obviously

        DO i=1, npnts
          hsat_ek(i)=0.0
          hsat_ekp1(i)=0.0
          hsat_ekp12(i)=0.0
          h_pk(i)=0.0
          h_pkp12(i)=0.0
          h_ek(i)=0.0
          h_ekp1(i)=0.0
          h_ekp12(i)=0.0
          delta_hsat_env(i)=0.0
          dhpbydp(i)=0.0
          dhebydp(i)=0.0
          epsilonp14(i)=0.0
          epsilonP34(i)=0.0
          ekp14_ad(i)=0.0
          ekp34_ad(i)=0.0
        END DO
        DO i=1,npnts

!Calculate hsat in case where adaptive entrainment is switched on

          IF(bwk(i)) THEN
            EL = LC
          ELSE
            EL = LC+LF
          END IF
          gzk       = g*zk(i)
          cpexnerk  = cp * exner_layer_centres(i,k)
          hsat_ek(i)= cpexnerk * thek(i) + EL * qsek(i) + gzk
          h_ek(i)   = cpexnerk * thek(i) + EL * qek(i)  + gzk
          h_pk(i)   = cpexnerk * thpk(i) + EL * qpk(i)  + gzk

          IF(bwkp1(i)) THEN
            EL = LC
          ELSE
            EL = LC+LF
          END IF
          gzkp1         = g*zkp1(i)
          cpexnerkp1_th = cp * exner_layer_centres(i,k+1)* thekp1(i)
          hsat_ekp1(i)= cpexnerkp1_th  + EL * qsekp1(i) + gzkp1
          h_ekp1(i)   = cpexnerkp1_th  + EL * qekp1(i)  + gzkp1

          delta_hsat_env(i)=(hsat_ekp1(i)-hsat_ek(i))


!version below has stability criterion, but doesn't really help
!          delta_hsat_env(i)=MIN((hsat_ekp1(i)-hsat_ek(i)), 0.0)
!Not sure precisely what to do with this- if delta_hsat_env
!before correction was positive, convection shouldn't be taking
!place at all, never mind the value of the entrainment coefficient.


          hsat_ekp12(i) = hsat_ek(i)+ delpkp14(i)/delpkp12(i) *         &
                                              delta_hsat_env(i)

       !NB already inside loop over npnts

          dhpbydp(i) = (1.0 - RENT) * delta_hsat_env(i)/                &
                                              (-1.0*delpkp12(i))
          h_pkp12(i) = h_pk(i)  -delpkp14(i)* dhpbydp(i)

          dhebydp(i) = (h_ekp1(i) - h_ek(i))/(-1.0*delpkp12(i)) 

          h_ekp12(i) = h_ek(i)  - delpkp14(i)* dhebydp(i)

        END DO    
      END IF      ! test on adaptive entrainment

! ---------------------------------------------------------------------
! Calculate entrainment coefficients multiplied by approppriate
! layer thickness. Same calculation required whether adaptive 
! entrainment or not.  Stop exponential decrease of entrainment rate 
! with height above NTPAR to ensure massflux continues to decrease 
! significantly with height (characteristic of shallow)
! (UM Doc 27 section 2C, equation 14)
! ---------------------------------------------------------------------
!
      IF (l_shallow) THEN

        IF (l_var_entrain) THEN      ! new variable entrainment
          If  (l_const_ent) THEN         ! no height dependence
            DO i=1,npnts
              IF (entrain_coef(i)  > 0.0) THEN

                ekp14(I) = entrain_coef(i) * (zkp12(i) - zk(i))
                ekp34(I) = entrain_coef(i) * (zkp1(i)-zkp12(i))

              ELSE         ! use equilibrium values
                ekp14(i) = delpkp14(i) * 0.03                          &
                        * EXP( -1.0*MIN( 1.0,                          &
                 ( p_layer_boundaries(I,ntml(i))-pk(i) )/ delp_cld(i) )&
                             )             / delp_cld(i) 
                ekp34(i) = delpkp34(i) * 0.03                          &
                        * EXP( -1.0*MIN( 1.0,                          &
                 ( p_layer_boundaries(I,ntml(i))                       &
                           - p_layer_boundaries(i,k) )/ delp_cld(i) )  &
                             )             / delp_cld(i)

              END IF
              IF (K  ==  ntml(i)) ekp14(i) = 0.0
            END DO 
          ELSE           ! rates vary with height
            DO i=1,npnts
              IF (entrain_coef(i)  > 0.0) THEN

! Possible future use commented outto reduce CPU cost
!                IF (l_rh_dep) THEN 
!                  IF (rhum(i) > 1.0) THEN 
!                    factor = 1.0
!                  ELSE IF (rhum(i) <0.05 ) THEN
!                    factor = 20.    ! is this reasonable ?    
!                  ELSE
!                    factor = 1.0/rhum(i)
!                  END IF 
!                ELSE
                  factor = 1.0
!                END IF 
                ekp34(I) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))  &
                                         *factor / (zkp12(i)+zkp1(i))

                ekp14(I) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))   &
                                         *factor /(zkp12(i) + zk(i))
              ELSE         ! use equilibrium values
                ekp14(i) = delpkp14(i) * 0.03                          &
                        * EXP( -1.0*MIN( 1.0,                          &
                 ( p_layer_boundaries(I,ntml(i))-pk(i) )/ delp_cld(i) )&
                             )              / delp_cld(i) 

                ekp34(i) = delpkp34(i) * 0.03                          &
                        * EXP( -1.0*MIN( 1.0,                          &
                 ( p_layer_boundaries(I,ntml(i))                       &
                           - p_layer_boundaries(i,k) )/ delp_cld(i) )  &
                             )              / delp_cld(i)

              END IF    ! test on entrain_coef
              IF (K  ==  ntml(i)) ekp14(i) = 0.0
            END DO 
          END IF   !  test on l_const_ent

        ELSE IF (sh_grey == 1) THEN
          ! increase cloud base entrainment for grey zone param 
          ! - gives better match to LES
          DO i=1,npnts

            IF (K  >   ntml(i) ) THEN

              ekp14(i) = delpkp14(i) * 0.06                            &
                        * EXP( -1.5*MIN( 1.0,                          &
                ( p_layer_boundaries(i,ntml(i))-pk(i) )/ delp_cld(i) ) &
                             )             / delp_cld(i) 
            ELSE

              ekp14(i) = entcoef * ae2 * pk(i) * delpkp14(i) *         &
                        recip_pstar(i) * recip_pstar(i)
  
            END IF    ! level
            IF (K  ==  ntml(i)) ekp14(i) = 0.0
            IF ( K  >=  ntml(i) ) THEN

              ekp34(i) = delpkp34(i) * 0.06                            &
                         * EXP( -1.5*MIN( 1.0,                         &
                ( p_layer_boundaries(i,ntml(i))                        &
                            - p_layer_boundaries(i,k) )/ delp_cld(i) ) &
                              )             / delp_cld(i)

            ELSE ! Deep or mid level (or levels not used)        

              ekp34(i) = entcoef * aekp34 * (p_layer_boundaries(i,k)) *&
                        delpkp34(i) *                                  &
                         recip_pstar(i) * recip_pstar(i)

            END IF   ! type of convection and level
          END DO 
        ELSE             ! original code
          DO i=1,npnts

            IF (K  >   ntml(i) ) THEN

              ekp14(i) = delpkp14(i) * 0.03                            &
                        * EXP( -1.0*MIN( 1.0,                          &
                ( p_layer_boundaries(i,ntml(i))-pk(i) )/ delp_cld(i) ) &
                             )             / delp_cld(i) 
            ELSE

              ekp14(i) = entcoef * ae2 * pk(i) * delpkp14(i) *         &
                        recip_pstar(i) * recip_pstar(i)
  
            END IF    ! level
            IF (K  ==  ntml(i)) ekp14(i) = 0.0
            IF ( K  >=  ntml(i) ) THEN

              ekp34(i) = delpkp34(i) * 0.03                            &
                         * EXP( -1.0*MIN( 1.0,                         &
                ( p_layer_boundaries(i,ntml(i))                        &
                            - p_layer_boundaries(i,k) )/ delp_cld(i) ) &
                              )             / delp_cld(i)

            ELSE ! Deep or mid level (or levels not used)        

              ekp34(i) = entcoef * aekp34 * (p_layer_boundaries(i,k)) *&
                        delpkp34(i) *                                  &
                         recip_pstar(i) * recip_pstar(i)

            END IF   ! type of convection and level
          END DO 
        END IF      ! l_var_entrain and sh_grey eq 1
 
      !-----------------------------------------------------------------------
      ELSE IF (L_congestus) THEN
      !-----------------------------------------------------------------------

         IF (l_var_entrain) THEN      ! new variable entrainment

          IF (l_const_ent) THEN  

            DO i=1,npnts
              IF (entrain_coef(i) > 0.0 ) THEN

                ekp34(I) = entrain_coef(i) *(zkp1(i)-zkp12(i))
                ekp14(I) = entrain_coef(i) *(zkp12(i) - zk(i))
              ELSE       ! use equilibrium values
! 0.5/z rates
                ekp34(I) =0.5* 2.0 * (zkp1(i)-zkp12(i))/                     &
                                                (zkp12(i)+zkp1(i))
                ekp14(I) =0.5* 2.0 * (zkp12(i) - zk(i))/                     &
                                                (zkp12(i) + zk(i))
              END IF 
              IF (K  ==  ntml(i)) ekp14(i) = 0.0
            END DO 
          ELSE   ! rates vary with height
            DO i=1,npnts
              IF (entrain_coef(i) > 0.0 ) THEN
                ekp34(I) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))/        &
                                                 (zkp12(i)+zkp1(i)) 

                ekp14(I) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))/        &
                                                (zkp12(i) + zk(i))

              ELSE       ! use equilibrium values
! 0.5/z rates
                ekp34(I) = 0.5*2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
                ekp14(I) = 0.5*2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
              END IF 
              IF (K  ==  ntml(i)) ekp14(i) = 0.0
            END DO 
          END IF       ! l_const_ent    

        ELSE           ! original rates

          DO i=1,npnts
          
            IF (K  >   ntml(i) ) THEN
! 0.5/z rates
                ekp34(I) = 0.5*2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
                ekp14(I) = 0.5*2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
            ELSE   ! should not be used.

                ekp14(i) = 0.0
                ekp34(I) = 0.0

            END IF   ! type of convection and level
          END DO

        END IF    ! test on new variable entrainment

      !-----------------------------------------------------------------------
      ELSE IF (l_deep) THEN ! Deep 
      !-----------------------------------------------------------------------
  

        SELECT CASE (ent_opt_dp) 
        CASE(0)     ! orignal Ap/(p*)^2  style entrainment 
        IF (l_var_entrain) THEN      ! new variable entrainment

          IF (l_const_ent) THEN

            DO i=1,npnts

              IF (entrain_coef(i) > 0.0) THEN

                ekp14(I) = entrain_coef(i) * (zkp12(i) - zk(i))
                ekp34(I) = entrain_coef(i) * (zkp1(i)-zkp12(i))

              ELSE       ! use equilibrium values   
! Possible future use
!                IF (l_rh_dep) THEN 
!                  IF (rhum(i) > 1.0) THEN 
!                   factor = 1.0
!                  ELSE IF (rhum(i) <0.05 ) THEN
!                   factor = 20.    ! is this reasonable ?    
!                  ELSE
!                   factor = 1.0/rhum(i)
!                  END IF 
!                ELSE IF (l_rh_dep2) THEN    ! alternative
!                  IF (rhum(i) > 1.0) THEN 
!                   factor = 0.75
!                  ELSE
!                  factor = 1.75-rhum(i)
!                  END IF 
!                ELSE 
                  factor = 1.0
!                END IF 

                ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *     &
                        recip_pstar(i) * recip_pstar(i)*factor
                ekp34(i) = entcoef * aekp34 * p_layer_boundaries(i,k) * &
                         delpkp34(i) *                                 &
                         recip_pstar(i) * recip_pstar(i)*factor
              END IF 
            END DO 

          ELSE          ! rates vary with height

            DO i=1,npnts

! Possible future use
!                IF (l_rh_dep) THEN 
!                  IF (rhum(i) > 1.0) THEN 
!                    factor = 1.0
!                  ELSE IF (rhum(i) <0.05 ) THEN
!                    factor = 20.    ! is this reasonable ?    
!                  ELSE
!                    factor = 1.0/rhum(i)
!                  END IF 
!                ELSE IF (l_rh_dep2) THEN    ! alternative
!                  IF (rhum(i) > 1.0) THEN 
!                   factor = 0.75
!                  ELSE
!                   factor = 1.75-rhum(i)
!                  END IF 
!                ELSE                  
                  factor = 1.0
!                END IF 

              IF (entrain_coef(i) > 0.0) THEN
                ekp14(I) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))    &
                                   *factor   / (zkp12(i) + zk(i))
                ekp34(I) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))    &
                                   *factor   / (zkp12(i)+zkp1(i))
              ELSE       ! use equilibrium values   
 
                ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *     &
                          recip_pstar(i) * recip_pstar(i)*factor
                ekp34(i) = entcoef * aekp34 * p_layer_boundaries(i,k) * &
                          delpkp34(i) *                                   &
                           recip_pstar(i) * recip_pstar(i)*factor
              END IF    ! test on entrain_coef
            END DO 
          END IF       ! l_const_ent

        ELSE          ! original values
          DO i=1,npnts 
              ekp14(i) = entcoef * aekp14 * pk(i) * delpkp14(i) *        &
                          recip_pstar(i) * recip_pstar(i)
              ekp34(i) = entcoef * aekp34 * (p_layer_boundaries(i,k)) *  &
                          delpkp34(i) *                                  &
                         recip_pstar(i) * recip_pstar(i)
          END DO
        END IF  

        CASE(1)     !  n/z style entrainment  where n = ent_fac

          IF (l_var_entrain) THEN      ! new variable entrainment coming from 
                                       ! diagnosis

            IF (l_const_ent) THEN      ! Constant entrainment rates

              DO i=1,npnts

                IF (entrain_coef(i) > 0.0) THEN

                  ekp14(I) = entrain_coef(i) * (zkp12(i) - zk(i))
                  ekp34(I) = entrain_coef(i) * (zkp1(i)-zkp12(i))

                ELSE       ! use equilibrium values   

                  ekp14(i) = ent_fac_dp*2.0* (zkp12(i) - zk(i))              &
                                            /(zkp12(i) + zk(i)) 
                  ekp34(i) = ent_fac_dp*2.0* (zkp1(i)-zkp12(i))              &
                                            /(zkp12(i)+zkp1(i))

                END IF 
              END DO 

            ELSE          ! rates vary with height

              DO i=1,npnts 

! Possible future use
!                IF (l_rh_dep) THEN 
!                  IF (rhum(i) > 1.0) THEN 
!                    factor = 1.0
!                  ELSE IF (rhum(i) <0.05 ) THEN
!                    factor = 20.    ! is this reasonable ?    
!                  ELSE
!                    factor = 1.0/rhum(i)
!                  END IF 
!                ELSE IF (l_rh_dep2) THEN    ! alternative
!                  IF (rhum(i) > 1.0) THEN 
!                   factor = 0.75
!                  ELSE
!                   factor = 1.75-rhum(i)
!                  END IF 
!                ELSE
                  factor = 1.0
!                END IF 

                IF (entrain_coef(i) > 0.0) THEN
                  factor = entrain_coef(i) *factor
                ELSE       ! use equilibrium values   
                  factor = ent_fac_dp *factor
                END IF 

                ekp14(I) = factor *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i))
                ekp34(I) = factor *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))

              END DO 
            END IF       ! l_const_ent

          ELSE          ! Simple n/Z  entrainment rates
               
            DO i=1,npnts 
              ekp14(i) = ent_fac_dp *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i)) 
              ekp34(i) = ent_fac_dp *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
            END DO

          END IF  

        CASE(2)     ! Higher entrainment near surface (various versions have
                    ! been used during GA3.0 plus testing) current code has
                    ! version known as New 3. Other versions left commented out
                    ! for reference at present.  

          IF (l_var_entrain) THEN    ! new variable entrainment

            DO i=1,npnts
              IF (entrain_coef(i) > 0.0) THEN  ! n/z
                ekp14(I) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))    &
                                    / (zkp12(i) + zk(i))
                ekp34(I) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))    &
                                    / (zkp12(i)+zkp1(i))

              ELSE       ! use  factor * original Ap/(p*)^2 

                IF (pk(i) > 50000.0) THEN
! New 1           factor1=1.0 +3.0*(1.-(100000.0 - pk(i))/50000.0)
! new 2           factor1=1.0 +2.0*(1.-(100000.0 - pk(i))/50000.0)
! New 3 
                  factor1=1.0 +1.25*(1.-(100000.0 - pk(i))/50000.0)
                ELSE 
                  factor1=1.0
                END IF 
                IF (p_layer_boundaries(i,k) > 50000.0) THEN
! new 1           factor2=1.0+                                                &
!                        3.0*(1.0-(100000.0-p_layer_boundaries(i,k))/50000.0)
! new 2           factor2=1.0+                                                &
!                        2.0*(1.0-(100000.0- p_layer_boundaries(i,k))/50000.0)
! New 3   
                  factor2=1.0+                                                &
                          1.25*(1.0-(100000.0- p_layer_boundaries(i,k))/50000.0)
                ELSE 
                  factor2=1.0
                END IF 
                ekp14(i) = factor1*entcoef*aekp14 * pk(i) *                   &
                          delpkp14(i) * recip_pstar(i) * recip_pstar(i)
                ekp34(i) = factor2*entcoef*aekp34 *p_layer_boundaries(i,k)*   &
                          delpkp34(i) * recip_pstar(i) * recip_pstar(i)

              END IF    ! test on entrain_coef
            END DO 

          ELSE          ! factor * original Ap/(p*)^2 
            DO i=1,npnts 
              IF (pk(i) > 50000.0) THEN
! new 1          factor1=1.0 +3.0*(1.-(100000. - pk(i))/50000.)
! new 2          factor1=1.0 +2.0*(1.-(100000.0 - pk(i))/50000.0)
! New 3 
                factor1=1.0 +1.25*(1.-(100000.0 - pk(i))/50000.0)
              ELSE 
                factor1=1.0
              END IF 
              IF (p_layer_boundaries(i,k) > 50000.0) THEN
! New 1         factor2=1.0+3.0*(1.0-(100000.0-p_layer_boundaries(i,k))/50000.0)
! new 2         factor2=1.0+2.0*(1.0-(100000.0-p_layer_boundaries(i,k))/50000.0)
! New 3
                factor2=1.0+                                                &
                          1.25*(1.0-(100000.0- p_layer_boundaries(i,k))/50000.0)
              ELSE 
                factor2=1.0
              END IF 
              ekp14(i) = factor1*entcoef * aekp14 * pk(i) * delpkp14(i) *     &
                          recip_pstar(i) * recip_pstar(i)
              ekp34(i) = factor2*entcoef * aekp34 * p_layer_boundaries(i,k)*  &
                          delpkp34(i) * recip_pstar(i) * recip_pstar(i)
            END DO
          END IF  

        CASE(3)     ! factor * (A/p*)*((p/p*)^m)  style entrainment over sea
                    ! Diurnal over land if active 

          IF (l_var_entrain) THEN    ! new variable entrainment

            DO i=1,npnts
              IF (entrain_coef(i) > 0.0) THEN  ! n/z
                ekp14(I) = entrain_coef(i) *2.0 * (zkp12(i) - zk(i))    &
                                    / (zkp12(i) + zk(i))
                ekp34(I) = entrain_coef(i) *2.0 * (zkp1(i)-zkp12(i))    &
                                    / (zkp12(i)+zkp1(i))

              ELSE       !

                ekp14(i) = entcoef*aekp14*delpkp14(i)*recip_pstar(i) &
                         * ((pk(i)* recip_pstar(i))**ent_dp_power)
                ekp34(i) = entcoef*aekp34*delpkp34(i)*recip_pstar(i) &
                     * ((p_layer_boundaries(i,k)*recip_pstar(i))**ent_dp_power)

              END IF    ! test on entrain_coef
            END DO 

          ELSE       
            DO i=1,npnts 
              ekp14(i) = entcoef*aekp14*delpkp14(i)*recip_pstar(i) &
                         * ((pk(i)* recip_pstar(i))**ent_dp_power)
              ekp34(i) = entcoef*aekp34*delpkp34(i)*recip_pstar(i) &
                     * ((p_layer_boundaries(i,k)*recip_pstar(i))**ent_dp_power)
            END DO
          END IF  

        END SELECT    ! test on ent_opt_dp

      !-----------------------------------------------------------------------
      ELSE ! Mid level 
      !-----------------------------------------------------------------------

        SELECT CASE (ent_opt_md) 
        CASE(0)     ! orignal Ap/(p*)^2  style entrainment 
                
          DO i=1,npnts 
            ekp14(i) = entcoef * aekp14_md * pk(i) * delpkp14(i) *        &
                        recip_pstar(i) * recip_pstar(i)
            ekp34(i) = entcoef * aekp34_md * (p_layer_boundaries(i,k)) *  &
                        delpkp34(i) *                                     &
                        recip_pstar(i) * recip_pstar(i)
          END DO

        CASE(1)     !  n/z style entrainment  where n = ent_fac

          DO i=1,npnts 
            ekp34(i) = ent_fac_md *2.0 * (zkp1(i)-zkp12(i))/(zkp12(i)+zkp1(i))
            ekp14(i) = ent_fac_md *2.0 * (zkp12(i) - zk(i))/(zkp12(i) + zk(i)) 
          END DO

        CASE(2)     ! New 3 profile - higher near surface (option not used 
                    ! during GA3.0 plus testing).
                
          DO i=1,npnts 
            IF (pk(i) > 50000.0) THEN
              factor1=1.0 +1.25*(1.-(100000.0 - pk(i))/50000.0)
            ELSE 
              factor1=1.0
            END IF 
            IF (p_layer_boundaries(i,k) > 50000.0) THEN
              factor2=1.0+1.25*(1.0-(100000.0-p_layer_boundaries(i,k))/50000.0)
            ELSE 
              factor2=1.0
            END IF 
            ekp14(i) = factor1*entcoef * aekp14_md * pk(i) *                   &
                          delpkp14(i) * recip_pstar(i) * recip_pstar(i)
            ekp34(i) = factor2*entcoef * aekp34_md * p_layer_boundaries(i,k) * &
                          delpkp34(i) * recip_pstar(i) * recip_pstar(i)
          END DO

        CASE(3)     ! factor * (A/p*)*((p/p*)^m)  style entrainment over sea

          DO i=1,npnts 
            ekp14(i) = entcoef*aekp14_md*delpkp14(i)*recip_pstar(i) &
                       * ((pk(i)* recip_pstar(i))**ent_md_power)
            ekp34(i) = entcoef*aekp34_md*delpkp34(i)*recip_pstar(i) &
                     * ((p_layer_boundaries(i,k)*recip_pstar(i))**ent_md_power)
          END DO


        END SELECT    ! test on ent_opt

      END IF   ! type of convection 


! ----------------------------------------------------------------------
! Second part of adaptive entrainment calculation
! ----------------------------------------------------------------------

      IF (ent_on  ==  1 ) THEN  

        !-------------------
        !Calculate epsilon14 and epsilon 34
        !-------------------
        DO i=1,npnts 
          If(bconv(i)) THEN
            If((h_pk(i) - h_ek(i))  >   0.0) THEN
              epsilonp14(i)= dhpbydp(i)/MAX((h_pk(i) - h_ek(i)),500.0)

          !500 chosen with view to not letting buoyancy drop too much
          !approx 0.2 cp + 0.2/1000 L
            ELSE
               epsilonp14(i)= 0.0
            END IF
            If((h_pkp12(i) - h_ekp12(i))  >   0.0) THEN
               epsilonp34(i)= dhpbydp(i)/MAX((h_pkp12(i) - h_ekp12(i)),500.0)
            ELSE
               epsilonp34(i)= 0.0
            END IF

            If(((h_pk(i) - h_ek(i))  >   500.0)                        &
                    .and.(dhpbydp(i)  >   1.0e-9)) THEN
              ekp14_ad(i) = MIN((epsilonp14(i) *                       &
                                (pk(i) - p_layer_boundaries(i,k))),1.0)
            ELSE
              ekp14_ad(i) = 0.0
            END IF

            If(((h_pkp12(i) - h_ekp12(i))  >  1.0e-10)                 &
                        .and. (dhpbydp(i)  >  1.0e-10)) THEN
              ekp34_ad(i) = MIN((epsilonp34(i) *                       &
                                (p_layer_boundaries(i,k)  - pkp1(i))),1.0)

!Check on sign of difference
              If((((ekp14_ad(i) - ekm14(i))    <   0.0) .and.          &
                  ((ekp34_ad(i) - ekp14_ad(i)) >   0.0))               &
               .OR.                                                    &
                 (((ekp14_ad(i) - ekm14(i))    >   0.0) .and.          &
                  ((ekp34_ad(i) - ekp14_ad(i)) <   0.0))) THEN

                ekp34_ad(i) = ekp14_ad(i)

              END IF

            ELSE   
              ekp34_ad(i) = 0.0
            END IF     ! test  h gradients

!       reset values of entrainment coefficients to adaptive versions
!        ekp14(i)=MAX(ekp14_ad(i),0.0)
!        ekp34(i)=MAX(ekp34_ad(i),0.0)
!Try an average of the adaptive and original versions
!but using GR as minimum value below which  entrainment doesn't fall
            If(ekp14_ad(i)  >   0.0) THEN
              ekp14(i)=0.5*(ekp14_ad(i)+ekp14(i))
            END IF
            If(ekp34(i)  >   0.0) THEN
              ekp34(i)=0.5*(ekp34_ad(i)+ekp34(i))
            END IF
          END IF        ! bconv test

        END DO  ! npnts

!----------------------------------------------------------------------

      END IF      ! end of adaptive entrainment calculation

! ---------------------------------------------------------------------
!  Calculate mixing detrainment coefficient multiplied by appropriate
!  layer thickness.
!
!  UM Documentation paper 27, section (2C) equation(15)
! ---------------------------------------------------------------------
 
      IF (l_shallow) THEN 
        IF (l_new_det) THEN    ! alter relationship for detrainment
          DO i=1,npnts
           If(K == 1)Then
             amdetk(i) = 0.0
           ELSE IF (K  >=  ntml(i) ) THEN

             IF (entrain_coef(i) > 0.0) THEN    ! alter detrainment
               ! Trying 2* entrainment rates for shallow

                amdetk(i) = (1.0 + 1.0)*(ekp14(i) + ekp34(i))
               !    But set to 1 if RH is 100%    
                IF (rhum(i)  >   1.0) THEN
                  amdetk(i) = 1.0*(ekp14(i) + ekp34(i))
                END IF
             ELSE          ! use original values
              IF (rhum(i)  <=  0.85) THEN
                amdetk(i) = (1.0 + 0.3)*(ekp14(i) + ekp34(i))
              ELSE IF (rhum(i)  >   1.0) THEN
                amdetk(i) = 1.0*(ekp14(i) + ekp34(i))
              ELSE
                amdetk(i) = (1.0 + (0.3/0.15)*(1.0-rhum(i)))            &
                                         *(ekp14(i) + ekp34(i))
              END IF
             END IF
           END IF
          END DO !npnts

        ELSE IF (sh_grey == 1) THEN
          ! increase detrainment and remove RH dependence 
          ! - gives better match to LES of shallow cumulus
          DO i=1,npnts
            IF(K == 1)THEN
              amdetk(i) = 0.0
            ELSE IF (K  >=  ntml(i) ) THEN
              amdetk(i) = 1.5*(ekp14(i) + ekp34(i))
            ELSE
              amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
            END IF
          END DO !npnts

        ELSE  ! original code unaltered

          DO i=1,npnts
            IF(K == 1)THEN
              amdetk(i) = 0.0
            ELSE IF (K  >=  ntml(i) ) THEN
              IF (rhum(i)  <=  0.85) THEN
                amdetk(i) = (1.0 + 0.3)*(ekp14(i) + ekp34(i))
              ELSE IF (rhum(i)  >   1.0) THEN
                amdetk(i) = 1.0*(ekp14(i) + ekp34(i))
              ELSE
                amdetk(i) = (1.0 + (0.3/0.15)*(1.0-rhum(i)))            &
                                         *(ekp14(i) + ekp34(i))
              END IF
            ELSE
             amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
            END IF
          END DO !npnts

        END IF  ! test on l_new det and sh_grey

      ELSE IF (l_congestus) THEN

         !------------------------------------------------
         ! adaptive mixing detrainment used in all cases
         !------------------------------------------------
          DO i=1,npnts
            IF (K == 1) THEN
              amdetk(i) = 0.0
            ELSE IF (K  >=  ntml(i) .and. rhum(i) <= 1.0 ) THEN
              amdetk(i) = (1.0 -rhum(i))*(ekp14(i) + ekp34(i))
            ELSE
              amdetk(i) = 0.0
            END IF
          END DO !npnts



      ELSE IF (l_deep) THEN

!-----------------------------------------------------------------------
! Adaptive mixing detrainment allowed as an option.
! See UM documentation paper section 7.2 equation 7.3
!-----------------------------------------------------------------------

        IF (K == 1) THEN
          DO i=1,npnts
            amdetk(i) = 0.0
          END DO !npnts
        ELSE
          IF (mdet_on  ==  1) THEN

            IF (l_new_det) THEN   ! altered code
              DO i=1,npnts

                amdet_func(i)= amdet_fac*(1.0-rhum(i))  ! original value

                IF (K  >=  ntml(i) ) THEN
                  IF (rhum(i)  <=  1.0) THEN
                    amdetk(i) = amdet_func(i)*(ekp14(i) + ekp34(i))
                  ELSE
                    amdetk(i) = 0.0
                  END IF
                ELSE          ! original mixing detrainment option
                  amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
                END IF
              END DO !npnts

            ELSE          ! orginal code

              DO i=1,npnts
               IF (K  >=  ntml(i) ) THEN
                 IF (rhum(i)  <=  1.0) THEN
                   amdetk(i) = amdet_fac*(ekp14(i) + ekp34(i))*(1-rhum(i))
                 ELSE
                   amdetk(i) = 0.0
                 END IF
               ELSE          ! original mixing detrainment option
                 amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
               END IF
              END DO !npnts
            END IF 
          ELSE     ! not mdet
            IF (l_new_det) THEN 
           ! Alter relationship between detrainment and entrainment
           ! try D=0.77 E
              DO i=1,npnts
                IF (entrain_coef(i) > 0.0) THEN    ! alter detrainment
                  amdetk(i) = (ekp14(i) + ekp34(i)) * 0.77
                ELSE   
                  amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
                END IF
              END DO !npnts
            ELSE         ! Original values
              DO i=1,npnts
                amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
              END DO !npnts
            END IF  
          END IF     ! test on mdet_on
        END IF   ! test on k=1

      ELSE      ! mid_level scheme  (no adaptive mixing detrainment)

        IF (K == 1) THEN
          DO i=1,npnts
            amdetk(i) = 0.0
          END DO !npnts
        ELSE
          DO i=1,npnts
            amdetk(i) = (ekp14(i) + ekp34(i)) * (1.0-1.0/aekp34)
          END DO !npnts
        END IF

      END IF     ! type of convection scheme

!      write(6,*) ' ekp14 '
!      write(6,*) (ekp14(i),i=1,npnts)
!      write(6,*) ' ekp34 '
!      write(6,*) (ekp34(i),i=1,npnts)
!

      IF (lhook) CALL dr_hook('LAYER_CN',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LAYER_CN
