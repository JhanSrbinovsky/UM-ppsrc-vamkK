! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Subroutine to generate sub-columns for McICA

! Method:
!       Allocates arrays to store sub-columns, rearranges cloud fields
!       so that they are in suitable order for generator and calls
!       generator which fills these arrays. Copies cloudy sub-columns 
!       if more are required

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
SUBROUTINE open_cloud_gen(                                              &
!                       Parallel variables
  global_row_length, global_rows                                        &
, me, n_proc, at_extremity                                              &
!                       Model Dimensions
, row_length, rows                                                      &
, n_layer, n_wet_layer                                                  &
, n_profile, p_temp                                                     &
!                       Properties of clouds
, w_cloud1, dp_corr_strat, cct, cloud_levels                            &
!                       Model switches
, l_rad_step_diag, l_rad_step_prog, model_domain                        &
!                       in time stepping information.
, val_year, val_day_number, val_hour, val_minute                        &
, val_second                                                            &
!                       error information
, ierr)

  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE rad_pcf
  USE rand_no_mcica
  USE mcica_mod
  USE sw_control_struct
  USE lw_control_struct
  USE spec_sw_lw
  USE rad_input_mod, ONLY: rad_mcica_sampling, rad_mcica_sigma
  USE cld_generator_mod
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParParams
  USE UM_ParVars, ONLY : g_datastart
  USE level_heights_mod, ONLY: r_theta_levels
  IMPLICIT NONE




!                       Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
    ierr
!           Error flag
!
!                       Parallel Variables
  INTEGER, INTENT(IN) ::                                                &
    global_row_length                                                   &
!           Number of global points on a row
  , global_rows                                                         &
!           Number of global rows
  , n_proc                                                              &
!           Total number of processors
  , me
!           My processor number
!
  LOGICAL, INTENT(IN) ::                                                &
    at_extremity(4)  
!           Indicates if this processor is at north, south
!           east or west of the processor grid
!
!
!                       Dimensions of arrays
  INTEGER, INTENT(IN) ::                                                &
    row_length                                                          &
!           Number of points on a row
  , rows                                                                &
!           Number of rows
  , n_profile                                                           &
!           Size allocated for atmospheric profiles
  , n_layer                                                             &
!           Size allocated for atmospheric layers
  , n_wet_layer
!           number of potentially cloudy layers
!
!                                                        
!                       Information about model time
  INTEGER, INTENT(IN) ::                                                &
    val_year                                                            &
!           time information for current timestep
  , val_day_number                                                      &
!           time information for current timestep
  , val_hour                                                            &
!           time information for current timestep
  , val_minute                                                          &
!           time information for current timestep
  , val_second
!           time information for current timestep

  REAL, INTENT(IN) ::                                                   &
    p_temp(n_profile, 0:n_layer)                                        &
!           Pressure
  , w_cloud1(n_profile, n_wet_layer)
!           Amount of cloud
!
!
!                       Properties of clouds
  REAL, INTENT(IN) ::                                                   &
    dp_corr_strat 
!           Decorrelation pressure scale for large scale cloud!

  INTEGER, INTENT(IN) ::                                                &
    cct(n_profile)
!           Level of top of convective cloud

  INTEGER, INTENT(OUT) ::                                               &
    cloud_levels
!           Number of global cloud levels

!                       Model Switches
  LOGICAL ::                                                            &
    l_rad_step_diag                                                     &
!           true if fast radiation timestep    (3C)
  , l_rad_step_prog
!           true if slow radiation timestep    (3C)
!
  INTEGER ::                                                            &
    model_domain
!           whether the model is global or limited area.             
!
!
!                       Local variables.
  INTEGER ::                                                            &
    i                                                                   &
!           Loop variable
  , j                                                                   &
!           Loop variable
  , k                                                                   &
!           Loop variable
  , l                                                                   &
!           Loop variable
  , info                                                                &
!           Loop variable
  , list(n_profile)                                                     &
!           Global location of each point on the processor.
  , random_dummy_init                                                   &
!           Seed for generating seeds for random numbers used in the
!           generator          
  , random_dummy(n_profile)                                             &
!           Seed for generating seeds for random numbers used in the
!           generator
  , first_row                                                           &
!           index of first row
  , last_row                                                            &
!           index of last row
  , cloud_top
!           Top global cloudy layer

  REAL ::                                                               &
    p(n_profile, n_layer)                                               &
!           Pressure
  , eps                                                                 &
!           small number to prevent division by zero.
  , cloud_scale                                                         &
!           cloud fraction times gridbox size
  , thickness_part(n_profile, n_layer)                                  &
!           part of FSD param related to layer thickness
  , dp_corr_cloud(n_profile,n_layer)                                    &
!           Cloud fraction decorrelation length
  , dp_corr_cond(n_profile,n_layer)                                     &
!           Cloud condensate decorrelation length
  , sigma_qcw(n_profile,n_layer)                                        &
!           Normalized cloud condensate std. dev
  , w_cloud(n_profile, 0:n_layer)
!           Amount of cloud

  REAL,PARAMETER ::                                                     &
    one_third = 1.0/3.0
!           1.0/3.0

  LOGICAL ::                                                            &
    l_layer_clear
!           Flag for layer free of clouds

  REAL, ALLOCATABLE :: temp_rand(:,:,:)
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('OPEN_CLOUD_GEN',zhook_in,zhook_handle)

  IF ((l_rad_step_prog) .OR. (l_rad_step_diag)) THEN

    IF ((sw_control(1)%i_cloud == ip_cloud_mcica) .OR.                  &
        (lw_control(1)%i_cloud == ip_cloud_mcica) .OR.                  &
        (sw_control(1)%i_fsd == ip_fsd_param)     .OR.                  &
        (lw_control(1)%i_fsd == ip_fsd_param)) THEN

      w_cloud=0.0

      DO j=1,n_wet_layer
        DO i=1,n_profile
          IF (w_cloud1(i,j) .gt. cut) THEN
            w_cloud(i,n_layer+1-j)=w_cloud1(i,j)
          END IF
        END DO
      END DO

    END IF

    IF ((sw_control(1)%i_fsd == ip_fsd_param) .OR.                      &
        (lw_control(1)%i_fsd == ip_fsd_param)) THEN

      DO k=1,n_layer
        DO j=1,rows
          DO i=1,row_length 
            l=(j-1)*row_length+i
            thickness_part(l,n_layer+1-k)=(0.001*(r_theta_levels(i,j,k) &
              -r_theta_levels(i,j,k-1)))**0.11
          END DO
        END DO
      END DO

      DO k=1,n_layer
        DO j=1,rows
          DO i=1,row_length
            l=(j-1)*row_length+i
            IF (w_cloud(l,k) == 1.0) THEN
              sigma_qcw(l,k)=0.21*(gridbox_size(i,j)**one_third)          &
                *((((0.016*gridbox_size(i,j))**1.10)+1)**(-0.26))         &
                *thickness_part(l,k)
            ELSE
              cloud_scale=w_cloud(l,k)*gridbox_size(i,j)
              sigma_qcw(l,k)=(0.41-(0.07*w_cloud(l,k)))                 &
                *(cloud_scale**one_third)                               &
                *((((0.016*cloud_scale)**1.10)+1)**(-0.26))             &
                *thickness_part(l,k)
            END IF
          END DO
        END DO
      END DO

      ALLOCATE(cloud_inhom_param_full(n_profile,1:n_layer))

      DO j=1, n_layer
        DO i=1, n_profile
! Calculate scaling parameter, =EXP(MEAN(LOG(water content)))/MEAN(WC)
! from standard deviation of water content, assuming a log-normal
! distribution to simplify the maths. 
          cloud_inhom_param_full(i,j)=(1.0/SQRT((sigma_qcw(i,j)**2)+1))
        END DO
      END DO
 
    ELSE
      DO j=1, n_layer
        DO i=1, n_profile
          sigma_qcw(i,j)=rad_mcica_sigma
        END DO
      END DO
    END IF

      
    IF ((sw_control(1)%i_cloud == ip_cloud_mcica) .OR.                  &
        (lw_control(1)%i_cloud == ip_cloud_mcica)) THEN

      cloud_levels=1
      DO j=1,n_wet_layer
        DO i=1,n_profile
          IF (w_cloud1(i,j) .gt. cut) THEN
            cloud_levels=j
          END IF
        END DO
      END DO


! To obtain reproducible results independent of the
! decomposition of the domain used on an MPP machine a global
! value for the topmost cloudy layer is used.
      CALL gc_imax(1, n_proc, info, cloud_levels)


      cloud_top=n_layer+1-cloud_levels

      IF  (l_rad_step_prog) THEN

        eps=EPSILON(eps)

! Set the SW and LW values of subcol_k (the number of sub-columns 
! each k-term "sees") and subcol_reorder (a reordering of the 
! sub-columns so that each sub-column is equivalently important in 
! the SW and LW).
!
        SELECT CASE (rad_mcica_sampling)

        CASE (ip_mcica_full_sampling)
          subcol_need=tot_subcol_gen
          sw_subcol_k=tot_subcol_gen
          lw_subcol_k=tot_subcol_gen

          ALLOCATE(lw_subcol_reorder(subcol_need))
          DO i=1,subcol_need
            lw_subcol_reorder(i)=i
          END DO

        CASE (ip_mcica_single_sampling)
          subcol_need=subcol_need_single
          sw_subcol_k=1
          lw_subcol_k=1

          ALLOCATE(lw_subcol_reorder(subcol_need))
          DO i=1,subcol_need
            lw_subcol_reorder(i) =                                      &
              MOD(lw_subcol_reorder_single(i),tot_subcol_gen)
          END DO

        CASE (ip_mcica_optimal_sampling)
          subcol_need=subcol_need_optimal
!          sw_subcol_k and lw_subcol_k have been read from data file

          ALLOCATE(lw_subcol_reorder(subcol_need))
          DO i=1,subcol_need
            lw_subcol_reorder(i) =                                      &
              MOD(lw_subcol_reorder_optimal(i),tot_subcol_gen)
          END DO
       
        END SELECT

        ALLOCATE(sw_subcol_reorder(subcol_need))
        DO i=1,subcol_need
          sw_subcol_reorder(i)=i
        END DO

        DO j=1,n_layer
          DO i=1,n_profile
            p(i,n_layer+1-j)=p_temp(i,j)
          END DO
        END DO

        DO j=1, n_layer
          DO i=1, n_profile
            dp_corr_cloud(i,j)=dp_corr_strat
            dp_corr_cond(i,j)=dp_corr_cloud(i,j)*0.5
          END DO
        END DO
     
        ALLOCATE(clw_sub_full(n_profile,cloud_top:n_layer,tot_subcol_gen))
!        ALLOCATE(cic_sub_full(n_profile,cloud_top:n_layer,tot_subcol_gen))
        ALLOCATE(temp_rand(n_profile,1,tot_subcol_gen))
        ALLOCATE(frac_cloudy_full(n_profile))
        ALLOCATE(ncldy(n_profile))

        k = 0
        first_row=0
        last_row=rows-1
        IF (.NOT. l_vatpoles) THEN
          IF (model_domain.eq.mt_global) THEN
            IF (at_extremity(psouth)) THEN
              first_row=1
              DO i=0, row_length-1
                k = k + 1
                list(k)=1
              END DO
            END IF
            IF (at_extremity(pnorth)) THEN
              last_row = rows-2
            END IF
          END IF
        END IF ! vatpoles
        DO j=first_row, last_row
          DO i=0, row_length-1
            k = k + 1
            list(k) = (g_datastart(2,me)+j-1)*global_row_length         &
                       + g_datastart(1,me)+i
          END DO
        END DO
        IF (.NOT. l_vatpoles) THEN
          IF ((model_domain.eq.mt_global).AND.(at_extremity(pnorth))) THEN
            DO i=0, row_length-1
              k = k + 1
              list(k)=(global_rows-1)*global_row_length+1
            END DO
          END IF
        END IF ! vatpoles

        random_dummy_init=(ABS(val_year-2000)*366*24*60*60)             &
          +(val_day_number*24*60*60)+(val_hour*60*60)                   &
          +(val_minute*60)+val_second
     
        DO i=1,n_profile
          random_dummy(i)=random_dummy_init+list(i)
        END DO
    
! The initial values of random dummy are successive integers, which 
! result in random numbers that are close to each other.This first 
! call to mcica_rand_no is to ensure that the seed for each 
! profile is itself random.
  
        CALL mcica_rand_no(random_dummy, temp_rand,n_profile            &
          ,tot_subcol_gen)

        ALLOCATE(rand_seed_x(n_profile,tot_subcol_gen))

        DO i=1,tot_subcol_gen
          CALL mcica_rand_no(random_dummy, temp_rand,n_profile,1)
          DO j=1, n_profile
            rand_seed_x(j,i)=random_dummy(j)
          END DO
        END DO
      
        ALLOCATE(rand_seed_y(n_profile,tot_subcol_gen))

        DO i=1,tot_subcol_gen
          CALL mcica_rand_no(random_dummy, temp_rand,n_profile,1)
          DO j=1, n_profile
            rand_seed_y(j,i)=random_dummy(j)
          END DO
        END DO

! Zero out fields
        DO i = 1,n_profile
          ncldy(i)  = 0
        END DO ! i

        DO k = 1, tot_subcol_gen
          DO j = cloud_top, n_layer
            DO i = 1, n_profile
!              cic_sub_full(i,j,k) = 0.0
              clw_sub_full(i,j,k) = 0.0
            END DO
          END DO
        END DO

! Set the overlap used in the cloud generator
        ioverlap=sw_control(1)%i_overlap

        CALL cld_generator(n_layer, cloud_top, n_profile               &
        , tot_subcol_gen, dp_corr_cloud, dp_corr_cond, sigma_qcw       &
        , w_cloud, p, 1, n_profile)

        IF (ALLOCATED(rand_seed_y)) DEALLOCATE(rand_seed_y)     
        IF (ALLOCATED(rand_seed_x)) DEALLOCATE(rand_seed_x)     
        IF (ALLOCATED(temp_rand)) DEALLOCATE(temp_rand)     

        IF (rad_mcica_sampling /= ip_mcica_full_sampling) THEN
          DO i=1,n_profile
            frac_cloudy_full(i)=REAL(ncldy(i))/REAL(tot_subcol_gen)
          END DO
        ELSE IF (rad_mcica_sampling == ip_mcica_full_sampling) THEN
! In this case we treat the clear sub-columns as cloudy sub-columns
! so the clear-sky fraction is implicit in the summing of the 
! "cloudy" sub-columns
          DO i=1,n_profile
            frac_cloudy_full(i)=1.0
          END DO
        END IF

        IF (rad_mcica_sampling /= ip_mcica_full_sampling) THEN
! For the case where are less cloudy subcolumns than required, 
! copy cloudy values to ensure enough cloudy subcolumns
          DO i=1, n_profile
            IF (ncldy(i) < subcol_need .AND. ncldy(i) > 0) THEN
              DO j=ncldy(i)+1,subcol_need
                DO k=cloud_top,n_layer
                  l=j-ncldy(i)
                  clw_sub_full(i,k,j)=clw_sub_full(i,k,l)
!                  cic_sub_full(i,k,j)=cic_sub_full(i,k,l)
                END DO
              END DO
            END IF
          END DO
        END IF
      END IF! l_rad_step_prog

    ELSE! If not using McICA in either SW or LW 

      l_layer_clear=.TRUE.
      DO j=n_wet_layer,1,-1
        DO i=1,n_profile
          l_layer_clear = l_layer_clear .AND.                             &
            (w_cloud1(i, j) <= 0.0e+00) .AND. (cct(i) < j-1)
        END DO
        IF (.NOT.l_layer_clear) THEN
          cloud_levels = j
          EXIT
        END IF
      END DO


!     To obtain reproducible results independent of the
!     decomposition of the domain used on an MPP machine a global
!     value for the topmost cloudy layer is used.
      CALL gc_imax(1, n_proc, info, cloud_levels)


    END IF! not mcica

  END IF! l_rad_step_prog .OR. l_rad_step_diag

  IF (lhook) CALL dr_hook('OPEN_CLOUD_GEN',zhook_out,zhook_handle)

END SUBROUTINE open_cloud_gen
