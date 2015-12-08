! *****************************COPYRIGHT******************************* 
! 
! (c) [University of Cambridge] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! 
! *****************************COPYRIGHT******************************* 
!                                                                       
! Purpose: Subroutine to overwrite values at top of model            
!          using interpolated 5-day fields from the 2-d model.          
!          Based on STRATF.F from Cambridge TOMCAT model and            
!          modified by Olaf Morgenstern to allow for flexible           
!          positioning of the NOy species.                              
!                                                                       
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTRY_CTL.                              
!                                                                       
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v8.3 programming standards.
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
      SUBROUTINE UKCA_STRATF(i_day_number, row_length, rows,            &
                             model_levels, theta_field_size,            &
                             first_row, global_row_length,              &
                             global_rows, ntracer, sinlat,              &
                             pres, um_ozone3d, p_above, tracer)         

      USE ASAD_MOD,          ONLY: advt
      USE UKCA_TROPOPAUSE,   ONLY: tropopause_level
      USE UKCA_CSPECIES,     ONLY: n_o3, n_o3s, n_hono2, n_ch4, n_no,   &
                                   n_no2, n_no3, n_n2o5, n_ho2no2
      USE UKCA_CONSTANTS,    ONLY: c_o3, c_n, c_ch4, c_hno3, c_no,      &
                                   c_no3, c_no2, c_n2o5, c_hono2,       &
                                   c_ho2no2
      USE ukca_option_mod,   ONLY: L_ukca_use_2dtop, strat2d_dir
      USE parkind1,          ONLY: jprb, jpim
      USE yomhook,           ONLY: lhook, dr_hook
      USE UM_ParVars
      USE Control_Max_Sizes
      USE ereport_mod,       ONLY: ereport
      USE PrintStatus_mod
      USE param2d_mod, ONLY: nolat,nolev,nlphot,ntphot,jpphio,jpphin
      IMPLICIT NONE  

      INTEGER, INTENT(IN) :: row_length         ! No of longitudes            
      INTEGER, INTENT(IN) :: rows               ! No of latitudes             
      INTEGER, INTENT(IN) :: model_levels       ! No of levels                
      INTEGER, INTENT(IN) :: theta_field_size   ! No of spatial points        
      INTEGER, INTENT(IN) :: first_row          ! First global latitude       
      INTEGER, INTENT(IN) :: global_row_length  ! No of global longitudes     
      INTEGER, INTENT(IN) :: global_rows        ! No of global latitudes      
      INTEGER, INTENT(IN) :: ntracer            ! No of chemical tracers
                                                                              
      REAL, INTENT(IN) :: sinlat(theta_field_size)           ! Sine(latitude)
      REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! Model pressures
      REAL, INTENT(IN) :: um_ozone3d(row_length,rows,model_levels) ! UM O3
      REAL, INTENT(IN) :: p_above               ! above which tbc are applied            

! Tracer concentrations in mass mixing ratio
      REAL, INTENT(INOUT) :: tracer(row_length,rows,model_levels,ntracer)
                                                                              
!     Local variables                                                       

      INTEGER :: i_day_number   ! Day number                                 
      INTEGER :: info           ! Tag used in communication
      INTEGER :: i,ij,j,k,l     ! Loop variables
      INTEGER :: m              ! Loop variable
     
      LOGICAL :: MASK(row_length,rows,model_levels) ! mask to identify stratosphere
      
      REAL, PARAMETER :: o3_hno3_ratio = 1.0/1000.0 ! kg[N]/kg[O3] from 
                                                    ! Murphy and Fahey 1994  
 
      REAL :: noy(nolat,nolev)                      ! 2D field interpolated in time
      REAL :: o3(nolat,nolev)                       ! 2D field interpolated in time
      REAL :: ch4(nolat,nolev)                      ! 2D field interpolated in time
      REAL :: o33d(row_length,rows,model_levels)    ! 2D field interpolated onto 3D
      REAL :: noy3d(row_length,rows,model_levels)   ! 2D field interpolated onto 3D
      REAL :: ch43d(row_length,rows,model_levels)   ! 2D field interpolated onto 3D
      REAL :: hno33d(row_length,rows,model_levels)  ! 3D field from fixed o3:hno3 ratio      
      REAL :: noyzm(rows,model_levels)              ! NOy zonal mean            
      REAL :: hno3zm(rows,model_levels)             ! HNO3 zonal mean                  
      REAL :: hno4zm(rows,model_levels)             ! HNO4 zonal mean                  
      REAL :: nozm(rows,model_levels)               ! NOz zonal mean            
      REAL :: no2zm(rows,model_levels)              ! NO2 zonal mean                   
      REAL :: no3zm(rows,model_levels)              ! NO3 zonal mean                   
      REAL :: n2o5zm(rows,model_levels)             ! N2O5 zonal mean           

! Parameter to use UM ancil ozone instead of 2-D O3
      LOGICAL, PARAMETER :: L_use_UMO3 = .TRUE.

! Parameter to use fixed o3 to hno3 ratio rather than 2d NOy
      LOGICAL, PARAMETER :: L_use_O3HNO3ratio = .TRUE.

! Parameter to overwrite stratosphere (fixed no of levels above tropopause)
      LOGICAL, PARAMETER :: L_all_strat   = .TRUE.
      INTEGER, PARAMETER :: no_above_trop1 = 3        ! Suitable for L38/L60
      INTEGER, PARAMETER :: no_above_trop2 = 10       ! Suitable for L63/L85
      INTEGER :: no_above_trop

! Parameter to overwrite CH4 with 2D boundary conditions
      LOGICAL, PARAMETER :: L_overwrite_CH4 = .FALSE.

      INTEGER           :: errcode                    ! Error code: ereport
      CHARACTER(LEN=72) :: cmessage                   ! Error message
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_STRATF',zhook_in,zhook_handle)

! Set number of levels above tropopause depending on vertical resolution
      IF (model_levels == 38 .OR. model_levels == 60) THEN
        no_above_trop = no_above_trop1
      ELSE IF (model_levels == 63 .OR. model_levels == 70 .OR.          &
               model_levels == 85) THEN
        no_above_trop = no_above_trop2
      ELSE
        errcode = 1
        cmessage = 'Levels above tropopause not set at this resolution'
        CALL EREPORT('UKCA_STRATF',errcode,cmessage)
      END IF

      IF (printstatus == Prstatus_Diag) THEN
        WRITE(6,*) 'UKCA_STRATF: Logicals in use:'
        WRITE(6,*) 'L_overwrite_CH4: ',L_overwrite_CH4
        WRITE(6,*) 'L_all_strat: ',L_all_strat
        WRITE(6,*) 'L_use_O3HNO3ratio: ',L_use_O3HNO3ratio
        WRITE(6,*) 'L_use_UMO3: ',L_use_UMO3
        WRITE(6,*) 'L_ukca_use_2dtop: ',L_ukca_use_2dtop
        WRITE(6,*) 'no_above_trop: ',no_above_trop
        WRITE(6,*) ' '
      END IF

      IF (L_ukca_use_2dtop) THEN
! Using data from Cambridge 2D model
! Get 2D data for current time then  interpolate in space

! DEPENDS ON: ukca_2d_bc_read_interp
        CALL UKCA_2D_BC_READ_INTERP(mype,nproc,i_day_number,strat2d_dir, &
                                    noy,o3,ch4)                                                                        
      END IF    ! End of IF(L_ukca_use_2dtop) statement

! Interpolate onto 3-D lats and levs and convert to mass mixing ratio (O3,CH4)
!  or calculate as VMR.                              
      IF (L_use_UMO3) THEN
        o33d(:,:,:) = um_ozone3d(:,:,:)
      ELSE
! DEPENDS ON: ukca_interp
        CALL UKCA_INTERP(o3, pres, row_length, rows, model_levels,      &
                         theta_field_size, sinlat, o33d)
        o33d(:,:,:) = o33d(:,:,:)*c_o3
      END IF

      IF (L_use_O3HNO3ratio) THEN
        hno33d(:,:,:) = o33d(:,:,:)*o3_hno3_ratio*c_hno3/c_n
      ELSE
! DEPENDS ON: ukca_interp
        CALL UKCA_INTERP(noy, pres, row_length, rows, model_levels,     &
                         theta_field_size, sinlat, noy3d)
      END IF

      IF (.NOT. L_overwrite_CH4 ) THEN
! DEPENDS ON: ukca_interp
        CALL UKCA_INTERP(ch4, pres, row_length, rows, model_levels,     &
                         theta_field_size, sinlat, ch43d)
        ch43d(:,:,:) = ch43d(:,:,:)*c_ch4
      ENDIF

      IF (.NOT. L_use_O3HNO3ratio) THEN
! Calculate zonal means for NOy species                            
! DEPENDS ON: ukca_calc_noy_zmeans 
         CALL UKCA_CALC_NOY_ZMEANS(row_length, rows, model_levels,      &
                            global_row_length, global_rows, first_row,  &
                            tracer, noyzm, nozm, no2zm, no3zm,          &
                            n2o5zm, hno3zm, hno4zm)
      END IF

! Overwrite values in tracer array

!       New N fields = zonal 3D field * (2d NOy/3d NOy) - at present only
!       scale to total NOy - the partitioning between NOy species is 
!       not really correct because the species are treated differently
!       in the two models. This can be improved. --- Comments from TOMCAT

      IF (L_use_O3HNO3ratio .AND. (.NOT. L_all_strat)) THEN

! Overwrite o3, hno3, and ch4 vmr at all gridboxes above p_above
!  O3   - either UM ancillary or Cambridge 2d
!  HNO3 - using fixed o3:hno3 ratio
!  CH4  - from Cambridge 2D model

        MASK(:,:,:) = .FALSE.
        WHERE (pres <= p_above) MASK = .TRUE.

        WHERE(MASK(:,:,:))
          tracer(:,:,:,n_o3)    = o33d(:,:,:)
          tracer(:,:,:,n_hono2) = hno33d(:,:,:)
        ENDWHERE
        IF (n_o3s > 0) THEN
          WHERE(MASK(:,:,:))
            tracer(:,:,:,n_o3s) = o33d(:,:,:)
          ENDWHERE
        END IF

        IF (L_overwrite_ch4) THEN
          WHERE(MASK(:,:,:))
            tracer(:,:,:,n_ch4) = ch43d(:,:,:)
          ENDWHERE
        ENDIF

      ELSE IF (L_use_O3HNO3ratio .AND. L_all_strat) THEN

! Overwrite o3, hno3, and ch4 at all gridboxes a fixed
!  number of model levels above the tropopause
!  O3   - either UM ancillary or Cambridge 2d
!  HNO3 - using fixed o3:hno3 ratio
!  CH4  - from Cambridge 2D model

        MASK(:,:,:) = .FALSE.
        DO l = model_levels,1,-1
          MASK(:,:,l) = tropopause_level(:,:)+no_above_trop <= l
        ENDDO

        WHERE (MASK(:,:,:))
          tracer(:,:,:,n_o3)    = o33d(:,:,:)
          tracer(:,:,:,n_hono2) = hno33d(:,:,:)
        ENDWHERE
        IF (n_o3s > 0) THEN
          WHERE (MASK(:,:,:))
            tracer(:,:,:,n_o3s) = o33d(:,:,:)
          ENDWHERE
        END IF

        IF (L_overwrite_ch4) THEN
          WHERE(MASK(:,:,:))
            tracer(:,:,:,n_ch4) = ch43d(:,:,:)
          ENDWHERE
        ENDIF
      ELSE

! Overwrite O3, NOy and CH4 at all gridboxes above p_above
!  O3  - use either UM ancillary or Cambridge 2D
!  NOy - use Cambridge 2D
!  CH4 - use Cambridge 2D (if L_overwrite_ch4)

        MASK(:,:,:) = .FALSE.
        WHERE (pres <= p_above) MASK = .TRUE.

        DO k = 1,ntracer
          SELECT CASE (advt(k))
            CASE ('O3        ')
              WHERE(MASK) tracer(:,:,:,n_o3) = o33d(:,:,:)
            CASE ('O3S       ')
              WHERE(MASK) tracer(:,:,:,n_o3s) = o33d(:,:,:)
            CASE ('CH4       ')
              IF (L_overwrite_ch4) THEN
                WHERE(MASK) tracer(:,:,:,n_ch4) = ch43d(:,:,:)
              END IF
            CASE ('NO        ')
              WHERE(MASK) tracer(:,:,:,n_no) = noy3d(:,:,:)*            &
                  SPREAD(nozm/noyzm,DIM=1,NCOPIES=row_length)*c_no
            CASE ('NO3       ')
              WHERE(MASK) tracer(:,:,:,n_no3) = noy3d(:,:,:)*           &
                  SPREAD(no3zm/noyzm,DIM=1,NCOPIES=row_length)*c_no3
            CASE ('NO2       ')
              WHERE(MASK) tracer(:,:,:,n_no2) = noy3d(:,:,:)*           &
                  SPREAD(no2zm/noyzm,DIM=1,NCOPIES=row_length)*c_no2
            CASE ('N2O5      ')
              WHERE(MASK) tracer(:,:,:,n_n2o5) = noy3d(:,:,:)*          &
                  SPREAD(n2o5zm/noyzm,DIM=1,NCOPIES=row_length)*c_n2o5
            CASE ('HO2NO2    ')
              WHERE(MASK) tracer(:,:,:,n_ho2no2) = noy3d(:,:,:)*        &
                  SPREAD(hno4zm/noyzm,DIM=1,NCOPIES=row_length)*c_ho2no2
            CASE ('HONO2     ')
              WHERE(MASK) tracer(:,:,:,n_hono2) = noy3d(:,:,:)*         &
                  SPREAD(hno3zm/noyzm,DIM=1,NCOPIES=row_length)*c_hono2
          END SELECT
        ENDDO

      END IF    ! End of IF(L_use_O3HNO3ratio and L_all_strat) statement

      IF (lhook) CALL dr_hook('UKCA_STRATF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_STRATF                                      
