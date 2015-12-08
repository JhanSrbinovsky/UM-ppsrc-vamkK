! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!+ data module for switches/options concerned with the cloud scheme.
  
  ! Module holding constants for the PC2 cloud scheme
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

MODULE pc2_constants_mod

  IMPLICIT NONE

  !=================================================================
  ! PC2 Cloud Scheme Terms
  !=================================================================

  ! Constants migrated from c_lspmic.h
  ! --------------------------------------------------------------------
  REAL, PARAMETER :: WIND_SHEAR_FACTOR = 1.5E-4
  ! Parameter that governs the rate of spread of ice cloud fraction
  ! due to windshear

  REAL, PARAMETER :: cloud_rounding_tol = 1.0E-12
  ! Tolerance to use when comparing cloud amounts, to allow for
  ! possible effects of rounding error.

  ! Constants migrated from ni_imp_ctl.F90 
  ! --------------------------------------------------------------------
  REAL, PARAMETER :: ls_bl0 = 1.0e-4 
  ! Specified in-plume value of ice content when there is not enough 
  ! information to calculate this parameter another way (kg kg-1)

  ! Erosion options
  ! --------------------------------------------------------------------
  INTEGER, PARAMETER :: pc2eros_exp_rh = 1  
  ! Original method described in Wilson et al. (2008). Uses the 
  ! notion of the relative rate of narrowing of the moisture PDF. 
  ! This rate is related pragmatically to RH, using an ad hoc 
  ! exponential. Change to qcl and CFL are both calculated as a 
  ! result of narrowing the PDF. 
  INTEGER, PARAMETER :: pc2eros_hybrid_allfaces = 2 
  ! Sink of qcl related to (qsat-qv), similarly to Tiedtke (1993) 
  ! eqn 30, but scaled according to surface area of cloud exposed 
  ! to clear sky. Change in CFL calculated assuming the change 
  ! in qcl was associated with a narrowing of PDF 
  ! (like Wilson et al. 2008). Surface area of exposed cloud is
  ! calculated considering the lateral sides and the top and bottom.
  INTEGER, PARAMETER :: pc2eros_hybrid_sidesonly = 3 
  ! As pc2_erosion_hybrid_method_all_faces but surface area of 
  ! exposed cloud is calculated considering the lateral sides only.

  ! Options for updating ice cloud fraction due to ice cloud fraction 
  ! falling in from layer above. 
  ! --------------------------------------------------------------------
  INTEGER, PARAMETER :: original_but_wrong = 1             
  ! Original code used in Wilson et al (2008a,b). It incorrectly 
  ! calculates the fall velocity of ice into the current layer from 
  ! the layer above, using the fall velocity in the current layer 
  ! (rather than the velocity in layer above). This allows ice to 
  ! fall into layer when there is no ice above! This option also uses 
  ! a hard-wired, globally constant wind_shear_factor, rather than 
  ! the shear derived from the model winds. 
  INTEGER, PARAMETER :: ignore_shear = 0 
  ! Correctly determines the fall velocity into the layer using the
  ! fall velocity of ice in the layer above. 
  INTEGER, PARAMETER :: real_shear = 2   
  ! Correctly determines the fall velocity into the layer using the 
  ! fall velocity of ice in the layer above. Adjusts the overhang by 
  ! using the true wind-shear calculated from the wind and also 
  ! translates the lateral displacement of the falling ice cloud due 
  ! to shear into a cloud fraction increment by considering the size 
  ! of the grid-box. 

! Constants migrated from c_pc2pdf.h

      ! Number of iterations in the initiation
      INTEGER,PARAMETER:: INIT_ITERATIONS=10

      ! Tolerance of cloud fraction for resetting of cloud
      ! Cloud_pc2_tol_2 should always be less than cloud_pc2_tol
      REAL,PARAMETER:: CLOUD_PC2_TOL   = 0.005
      REAL,PARAMETER:: CLOUD_PC2_TOL_2 = 0.001

      ! Tolerance of critical relative humidity for initiation of cloud
      REAL,PARAMETER:: RHCRIT_TOL=0.01

      ! Power that is used to weight the two values of G when they are
      ! merged together
      REAL,PARAMETER:: PDF_MERGE_POWER=0.5

      ! Power that describes the way G varies with s near the boundary
      ! of the cloud probability density function. For a "top-hat" it
      ! is equal to 0, for a "triangular" distribution it is equal to 1.
      REAL,PARAMETER:: PDF_POWER=0.0

      ! Parameters that govern the turbulent decrease in width of the
      ! PDF.  (dbs/dt)/bs = (DBSDTBS_TURB_0 + dQc/dt DBSDTBS_TURB_1)
      !                 * exp( - dbsdtbs_exp Q_c / (a_L qsat(T_L)))
      ! dbsdtbs_turb_0 is set in the UMUI
      REAL,PARAMETER:: DBSDTBS_TURB_1 = 0.0
      REAL,PARAMETER:: DBSDTBS_CONV   = 0.0
      REAL,PARAMETER:: dbsdtbs_exp    = 10.05

! Constant migrated from pc2_checks
  REAL, PARAMETER :: one_over_qcf0 = 1.0e4
! One_over_qcf0 is reciprocal of a reasonable in-cloud ice content.
  REAL, PARAMETER :: min_in_cloud_qcf = 1.0e-6
! Minium in-cloud ice water content ensured by reducing CFF
! Ensure value below is 1.0/value above
  REAL, PARAMETER :: one_over_min_in_cloud_qcf = 1.0e6
! One over MIN_IN_CLOUD_QCF.
  REAL, PARAMETER :: condensate_limit = 1.0e-10   
! Minimum value of condensate
  REAL, PARAMETER :: wcgrow                =   5.e-4 

END MODULE pc2_constants_mod
