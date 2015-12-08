! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!+ data module for switches/options concerned with the cloud scheme.
  ! Description:
  !   Module containing runtime options/data used by the cloud scheme
  !
  ! Method:
  !   Switches and associated data values used by the cloud scheme
  !   are defined here and assigned default values. These may be overridden
  !   by namelist input.
  !   
  !   A description of what each switch or number refers to is provided
  !   with the namelist
  !
  !   Any routine wishing to use these options may do so with the 'USE'
  !   statement.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Large Scale Cloud
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !

MODULE cloud_inputs_mod

  USE missing_data_mod, ONLY: RMDI, IMDI
  USE atmos_max_sizes,  ONLY: model_levels_max

  IMPLICIT NONE

!===========================================================================
! INTEGER options set from RUN_CLOUD namelist
!===========================================================================

 INTEGER :: pc2_falliceshear_method = IMDI
                                        ! Select method for including
                                        ! effects of shear on falling 
                                        ! ice effecting CFF. 
 INTEGER :: cloud_fraction_method   = IMDI
                                        ! Selects total cloud fraction
                                        ! calculation method
 INTEGER :: i_fixbug_pc2_checks =  IMDI ! Options for changing CFL and CF
                                        ! when creating extra
                                        ! QCL when qv>qsat.

 INTEGER :: i_pc2_conv_coupling = IMDI  ! Integer option to determine how the 
                                        ! PC2 cloud scheme and the convection 
                                        ! scheme are coupled and how cloud gets
                                        ! generated by convection. 

 INTEGER :: i_pc2_erosion_method = IMDI ! Select method for calculating
                                        ! PC2 cloud erosion.

 INTEGER :: forced_cu = IMDI            ! Options for representing forced 
                                        ! cumulus: currently 0=off, 1=on

!===========================================================================
! LOGICAL options set from RUN_CLOUD namelist
!===========================================================================

 LOGICAL :: L_eacf = .false.      ! Use empirically adjusted
                                  ! cloud fraction
 LOGICAL :: l_ensure_min_in_cloud_qcf = .FALSE.
                                  ! Reduce CFF when qcf is low to 
                                  ! ensure a minimum qcf/CFF. 
 LOGICAL :: l_fixbug_pc2_qcl_incr = .FALSE. 
                                  ! Prevent QCL incr removing more
                                  ! QCL than there is, down at the 
                                  ! individual subroutine level.
 LOGICAL :: l_fixbug_pc2_mixph = .FALSE.   
                                  ! Collection of bug-fixes related
                                  ! to how PC2 treats mixed phase
                                  ! cloud fraction
 LOGICAL :: l_micro_eros = .FALSE.
                                  ! If false, erosion is done as part of
                                  ! the convection scheme
                                  ! If true, erosion is done as part of 
                                  ! the microphysics scheme

 LOGICAL :: l_cld_area  = .FALSE. ! Controls cloud area parametrization

 LOGICAL :: l_acf_cusack = .FALSE.! Cusack method of cloud area parametrization
                                  ! based on temperature and moisture profiles

 LOGICAL :: l_pc2 = .FALSE.       ! Controls PC2 cloud scheme

 LOGICAL :: l_rhcpt = .FALSE.     ! Controls the use of new RHcrit
                                  ! parametrization option in Sec 9 vn 2A.

 LOGICAL :: l_filter_cloud = .FALSE. 
                                  ! Filters area and combined cloud amounts 
                                  ! to remove sub-visual cirrus prior to  
                                  ! calculating low/med/high and  
                                  ! total cloud amounts. 

!=======================================================================
! REAL values set from RUN_CLOUD namelist
!=======================================================================

 REAL    :: dbsdtbs_turb_0 = RMDI      ! PC2 erosion rate / s-1
 REAL    :: starticeTKelvin = RMDI   ! Temperature at which detrained
                                       ! condensate is assumed to start
                                       ! being in the ice phase (Kelvin)
 REAL    :: alliceTdegC = RMDI        ! Temperature at which all 
                                       ! detrained condensate is 
                                       ! assumed to be ice (deg C)
 REAL    :: cff_spread_rate = RMDI   ! Rate at which CFF spreads out.
 REAL    :: rhcrit(model_levels_max) = RMDI
                                       ! Critical relative humidity

 REAL    :: ice_width = RMDI           ! Specifies the ice content 
                                       ! (in terms of a fraction of qsat_liq)
                                       ! that corresponds to a factor of 
                                       ! two reduction in the width of
                                       ! the vapour distribution in the 
                                       ! liquid-free part of the gridbox.
 REAL    :: tau_thresh = RMDI          ! Value of optical depth below 
                                       ! which ice cloud is set to zero. 

!=======================================================================
! Constants, removed from RUN_CLOUD namelist
!=======================================================================

 INTEGER :: ice_fraction_method = 1     
                                        ! Selects ice cloud fraction
                                        ! calculation method

 REAL    :: overlap_ice_liquid = -1.0    
                                       ! Generic overlap parameter
                                       ! between ice and liquid phases
 REAL    :: ctt_weight = 0.333           
                                       ! Cloud top temperature weight
 REAL    :: t_weight = 0.333   
                                       ! Local temperature weight
 REAL    :: qsat_fixed = 0.1E-3        ! Prescribed qsat value
 REAL    :: sub_cld   = 0.225          ! Scaling factor

!=======================================================================
! Control options, removed from namelist
!=======================================================================

 LOGICAL :: l_pc2_reset = .FALSE. ! Run PC2 scheme diagnostically

 LOGICAL :: l_pc2_lbc = .FALSE.   ! LBCs contain cloud fractions

 LOGICAL :: l_acf_brooks = .FALSE.! Diagnostic Brooks method for cloud 
                                  ! parameterization for non-convective cloud 
!----------------------------------------------------------------------
!         
!----------------------------------------------------------------------

! Define the RUN_CLOUD namelist

 NAMELIST/RUN_Cloud/ rhcrit, L_eacf, cloud_fraction_method, forced_cu,  &
        dbsdtbs_turb_0, pc2_falliceshear_method,                        &
        l_ensure_min_in_cloud_qcf,                                      &
        l_fixbug_pc2_qcl_incr, i_fixbug_pc2_checks, l_fixbug_pc2_mixph, &
        i_pc2_conv_coupling, i_pc2_erosion_method, l_micro_eros,        &
        starticeTKelvin, alliceTdegC, cff_spread_rate, ice_width,       &
        l_cld_area, l_acf_cusack, l_pc2, l_rhcpt,                       &
        l_filter_cloud, tau_thresh

CONTAINS

SUBROUTINE check_run_cloud()

! Description:
!   Subroutine to apply logic checks based on the options selected in 
!   the run_cloud namelist.

! Dr Hook Modules
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod,    ONLY: ereport

IMPLICIT NONE

INTEGER                       :: icode         ! used for ereport
CHARACTER (LEN=120)           :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'check_run_cloud'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CHECK_RUN_CLOUD',zhook_in,zhook_handle)

IF (l_cld_area .AND.                                                     &
     ((.NOT. l_acf_cusack .AND. .NOT. l_acf_brooks) .OR.                 &
      ( l_acf_cusack .AND.  l_acf_brooks ))) THEN
  WRITE (cmessage,'(A89)') 'l_cld_area is set to true: either l_acf_brooks or '&
                            // 'l_acf_cusack must be true, but not both'
  icode = 98
  CALL ereport(RoutineName, icode, cmessage)
END IF

IF (lhook) CALL dr_hook('CHECK_RUN_CLOUD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_cloud

END MODULE cloud_inputs_mod