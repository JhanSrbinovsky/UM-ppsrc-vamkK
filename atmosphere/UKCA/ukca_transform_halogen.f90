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
!  Description:
!    To transform one compound each for Br, Cl, N, and H to become
!    the tracer for the total elemental abundance, and to reverse that
!    transformation.  Contains SUBROUTINE transform.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE ukca_transform_halogen(tr_ukca, rows, row_length,      &
                              model_levels,                             &
                              off_x, off_y, tracers, halo_i, halo_j,    &
                              q, forward, timestep)

! This subroutine changes from representations of total chlorine in HCl
! and fractions of total chlorine in the other tracers, to proper
! mass mixing ratios for chlorine and bromine, or reverses this
! transformation.

      USE ASAD_MOD,       ONLY: advt
      USE UKCA_CONSTANTS, ONLY: c_hcl, c_cl, c_clo, c_cl2o2, c_hocl,    &
                                c_oclo, c_clono2, c_cfcl3,              &
                                c_cf2cl2, c_mecl, c_mecl, c_ccl4,       &
                                c_cf2clcfcl2, c_chf2cl,  c_hbr,         &
                                c_meccl3, c_brcl, c_cf2clbr, c_br,      &
                                c_hobr, c_brono2, c_mebr, c_cf3br,      &
                                c_ch2br2, c_n, c_no, c_no3, c_n2o5,     &
                                c_ho2no2, c_hono2, c_hono, c_meono2,    &
                                c_pan, c_ppan, c_ison, c_mpan, c_n2o,   &
                                c_h2o2, c_ch4, c_hcho, c_meooh,         &
                                c_h2, c_h, c_oh, c_ho2, c_meoo, c_h2o,  &
                                c_no2, c_hcl, c_bro
      USE ukca_option_mod,ONLY: jpctr
      USE parkind1,       ONLY: jprb, jpim
      USE yomhook,        ONLY: lhook, dr_hook
      USE ereport_mod,    ONLY: ereport
      USE Control_Max_Sizes
      IMPLICIT NONE

! Subroutine interface

      INTEGER, INTENT(IN) :: tr_ukca
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: off_x
      INTEGER, INTENT(IN) :: off_y
      INTEGER, INTENT(IN) :: halo_i
      INTEGER, INTENT(IN) :: halo_j
      INTEGER, INTENT(IN) :: timestep

      REAL, INTENT(INOUT) :: tracers(1-off_x:row_length+off_x,          &
                     1-off_y:rows+off_y, model_levels, tr_ukca)

      REAL, INTENT(INOUT) :: q(1-halo_i:row_length+halo_i,              &
                     1-halo_j:rows+halo_j,model_levels)

      LOGICAL, INTENT(IN) :: forward

! Local variables
      INTEGER :: errcode                ! Variable passed to ereport

      INTEGER, PARAMETER :: tr_max = 30

! Status of compounds
      INTEGER, PARAMETER :: unrestricted   = 0
      INTEGER, PARAMETER :: source_gas     = 1
      INTEGER, PARAMETER :: mixed_compound = 2

      LOGICAL, SAVE :: first=.TRUE.
      INTEGER :: i

! Positions in tracer array of selected tracers
      INTEGER, SAVE :: n_hcl =0
      INTEGER, SAVE :: n_bro =0
      INTEGER, SAVE :: n_no2 =0
      INTEGER, SAVE :: n_h2o=0
      INTEGER, SAVE :: n_h2os=0

! Number of compounds of conserved elements
      INTEGER, SAVE :: n_cl_tracers=0
      INTEGER, SAVE :: n_br_tracers=0
      INTEGER, SAVE :: n_n_tracers =0
      INTEGER, SAVE :: n_h_tracers =0

! Position in tracer array of compounds
      INTEGER, SAVE :: cl_tracers(tr_max)=0
      INTEGER, SAVE :: br_tracers(tr_max)=0
      INTEGER, SAVE :: n_tracers (tr_max)=0
      INTEGER, SAVE :: h_tracers (tr_max)=0

! Status of compounds
      INTEGER, SAVE :: s_cl_tracers(tr_max)=unrestricted
      INTEGER, SAVE :: s_br_tracers(tr_max)=unrestricted
      INTEGER, SAVE :: s_n_tracers (tr_max)=unrestricted
      INTEGER, SAVE :: s_h_tracers (tr_max)=unrestricted

! Conversion factors for VMR to MMR
      REAL, SAVE :: c_cl_tracers(tr_max)
      REAL, SAVE :: c_br_tracers(tr_max)
      REAL, SAVE :: c_n_tracers (tr_max)
      REAL, SAVE :: c_h_tracers (tr_max)

! Correction factor
      REAL, ALLOCATABLE :: corrfac(:,:,:)

      CHARACTER(LEN=72) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_TRANSFORM_HALOGEN',zhook_in,zhook_handle)
      IF (first) THEN
        DO i=1,jpctr
          SELECT CASE (advt(i))
            CASE('HCl       ')
              n_hcl  = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hcl
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('Cl        ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cl
              cl_tracers(n_cl_tracers) = i
            CASE('ClO       ','Clx       ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_clo
              cl_tracers(n_cl_tracers) = i
            CASE('Cl2O2     ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cl2o2/2.
              cl_tracers(n_cl_tracers) = i
            CASE('HOCl      ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_hocl
              cl_tracers(n_cl_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hocl
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('OClO      ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_oclo
              cl_tracers(n_cl_tracers) = i
            CASE('ClONO2    ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_clono2
              cl_tracers(n_cl_tracers) = i
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_clono2
              n_tracers(n_n_tracers) = i
              s_n_tracers(n_n_tracers) = mixed_compound
            CASE('CFCl3     ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cfcl3/3.
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
            CASE('CF2Cl2    ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cf2cl2/2.
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
            CASE('MeCl      ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_mecl
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_mecl ! SHOULD THIS BE /3?
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('CCl4      ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_ccl4/4.
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
            CASE('CF2ClCFCl2')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cf2clcfcl2/3.
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
            CASE('CHF2Cl    ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_chf2cl
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_chf2cl
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('MeCCl3    ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_meccl3/3.
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = source_gas
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_meccl3 ! SHOULD THIS BE /3?
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('BrCl      ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_brcl
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = mixed_compound
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_brcl
              br_tracers(n_br_tracers) = i
            CASE('CF2ClBr   ')
              n_cl_tracers = n_cl_tracers + 1
              c_cl_tracers(n_cl_tracers) = c_cf2clbr
              cl_tracers(n_cl_tracers) = i
              s_cl_tracers(n_cl_tracers) = mixed_compound
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_cf2clbr
              br_tracers(n_br_tracers) = i
              s_br_tracers(n_cl_tracers) = source_gas
! bromine tracers
            CASE('BrO       ','Brx       ')
              n_bro = i
            CASE('Br        ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_br
              br_tracers(n_br_tracers) = i
            CASE('HBr       ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_hbr
              br_tracers(n_br_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hbr
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('HOBr      ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_hobr
              br_tracers(n_br_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hobr
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('BrONO2    ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_brono2
              br_tracers(n_br_tracers) = i
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_brono2
              n_tracers(n_n_tracers) = i
              s_n_tracers(n_n_tracers) = mixed_compound
            CASE('MeBr      ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_mebr
              br_tracers(n_br_tracers) = i
              s_br_tracers(n_br_tracers) = source_gas
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_mebr ! SHOULD THIS BE /3?
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('CF3Br     ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_cf3br
              br_tracers(n_br_tracers) = i
              s_br_tracers(n_br_tracers) = source_gas
            CASE('CH2Br2    ')
              n_br_tracers = n_br_tracers + 1
              c_br_tracers(n_br_tracers) = c_ch2br2/2.
              br_tracers(n_br_tracers) = i
              s_br_tracers(n_br_tracers) = 1
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_ch2br2/2.
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = 2
! Nitrogen tracers
            CASE('N         ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_n
              n_tracers(n_n_tracers) = i
            CASE('NO        ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_no
              n_tracers(n_n_tracers) = i
            CASE('NO2       ','NOx       ')
              n_no2 = i
            CASE('NO3       ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_no3
              n_tracers(n_n_tracers) = i
            CASE('N2O5      ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_n2o5/2.
              n_tracers(n_n_tracers) = i
            CASE('HO2NO2    ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_ho2no2
              n_tracers(n_n_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_ho2no2
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('HONO2     ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_hono2
              n_tracers(n_n_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hono2
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('HONO      ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_hono
              n_tracers(n_n_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hono
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('MeONO2    ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_meono2
              n_tracers(n_n_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_meono2/3.
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('PAN       ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_pan
              n_tracers(n_n_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_pan/3.
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
           CASE('PPAN      ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_ppan
              n_tracers(n_n_tracers) = i
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_ppan/5.
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = mixed_compound
            CASE('ISON      ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_ison
              n_tracers(n_n_tracers) = i
            CASE('MPAN      ')
              n_n_tracers = n_n_tracers + 1
              c_n_tracers(n_n_tracers) = c_mpan
              n_tracers(n_n_tracers) = i
! mozart chemistry only commented out for now
!            CASE('ONIT      ')
!              n_n_tracers = n_n_tracers + 1
!              c_n_tracers(n_n_tracers) = c_onit
!              n_tracers(n_n_tracers) = i
!            CASE('ONITR     ')
!              n_n_tracers = n_n_tracers + 1
!              c_n_tracers(n_n_tracers) = c_onitr
!              n_tracers(n_n_tracers) = i
! Hydrogen tracers
            CASE('H2O2      ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_h2o2/2.
              h_tracers(n_h_tracers) = i
            CASE('CH4       ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_ch4/4.
              h_tracers(n_h_tracers) = i
              s_h_tracers(n_h_tracers) = source_gas
            CASE('HCHO      ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_hcho/2.
              h_tracers(n_h_tracers) = i
            CASE('MeOOH     ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_meooh/4.
              h_tracers(n_h_tracers) = i
            CASE('H2        ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_h2/2.
              h_tracers(n_h_tracers) = i
            CASE('H         ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_h
              h_tracers(n_h_tracers) = i
            CASE('OH        ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_oh
              h_tracers(n_h_tracers) = i
            CASE('HO2       ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_ho2
              h_tracers(n_h_tracers) = i
            CASE('MeOO      ')
              n_h_tracers = n_h_tracers + 1
              c_h_tracers(n_h_tracers) = c_meoo/3.
              h_tracers(n_h_tracers) = i
            CASE('H2O       ')
              n_h2o = i
            CASE('H2OS      ')
              n_h2os = i
            CASE DEFAULT
              cmessage='Species: '//advt(i)//' not included in CASE'
              errcode=-1*i
              CALL EREPORT('UKCA_TRANSFORM_HALOGEN',errcode,cmessage)
          END SELECT
        END DO
        first = .FALSE.
      END IF !first
      IF ((n_h2os > 0) .AND. (timestep == 1)) THEN
        tracers(:,:,:,n_h2os) = q(1-off_x:row_length+off_x,             &
                                    1-off_y:rows+off_y, :)
        CALL transform(row_length, rows, model_levels, off_x,           &
                   off_y, tr_ukca, tracers, n_h_tracers, h_tracers,     &
                   c_h_tracers, s_h_tracers, n_h2os, c_h2o/2.,          &
                   .FALSE.)
        WHERE (tracers(:,:,:,n_h2os) < 4.3491e-6)                       &
          tracers(:,:,:,n_h2os) = 4.3491e-6
        tracers(:,:,model_levels/2:model_levels,n_h2os) = 4.3491e-6
        CALL transform(row_length, rows, model_levels, off_x,           &
                   off_y, tr_ukca, tracers, n_h_tracers, h_tracers,     &
                   c_h_tracers, s_h_tracers, n_h2os, c_h2o/2.,          &
                   .TRUE.)
        q(1-off_x:row_length+off_x,1-off_y:rows+off_y,:) =              &
          tracers(:,:,:,n_h2os)
        tracers(:,:,:,n_h2o) = tracers(:,:,:,n_h2os)
        CALL transform(row_length, rows, model_levels, off_x,           &
                   off_y, tr_ukca, tracers, n_h_tracers, h_tracers,     &
                   c_h_tracers, s_h_tracers, n_h2os, c_h2o/2.,          &
                   .FALSE.)
      END IF

      IF (forward) THEN
        IF (n_h2os > 0) THEN
          CALL transform(row_length, rows, model_levels, off_x,         &
                     off_y, tr_ukca, tracers, n_h_tracers, h_tracers,   &
                     c_h_tracers, s_h_tracers, n_h2os, c_h2o/2.,        &
                     forward)

          ! NOTE: may cause issues with q feedback on!!
          q(1-off_x:row_length+off_x,1-off_y:rows+off_y,                &
            model_levels-20:model_levels) =                             &
            tracers(:,:,model_levels-20:model_levels,n_h2os)
          tracers(:,:,model_levels-20:model_levels,n_h2o) =             &
            tracers(:,:,model_levels-20:model_levels,n_h2os)
        END IF

        CALL transform(row_length, rows, model_levels, off_x,           &
                     off_y, tr_ukca, tracers, n_n_tracers, n_tracers,   &
                     c_n_tracers, s_n_tracers, n_no2, c_no2,            &
                     forward)

        CALL transform(row_length, rows, model_levels, off_x,           &
                     off_y, tr_ukca, tracers, n_cl_tracers, cl_tracers, &
                     c_cl_tracers, s_cl_tracers, n_hcl, c_hcl,          &
                     forward)

        CALL transform(row_length, rows, model_levels, off_x,           &
                     off_y, tr_ukca, tracers, n_br_tracers, br_tracers, &
                     c_br_tracers, s_br_tracers, n_bro, c_bro,          &
                     forward)

      ELSE

        CALL transform(row_length, rows, model_levels, off_x,           &
                     off_y, tr_ukca, tracers, n_br_tracers, br_tracers, &
                     c_br_tracers, s_br_tracers, n_bro, c_bro,          &
                     forward)

        CALL transform(row_length, rows, model_levels, off_x,           &
                     off_y, tr_ukca, tracers, n_cl_tracers, cl_tracers, &
                     c_cl_tracers, s_cl_tracers, n_hcl, c_hcl,          &
                     forward)

        CALL transform(row_length, rows, model_levels, off_x,           &
                     off_y, tr_ukca, tracers, n_n_tracers, n_tracers,   &
                     c_n_tracers, s_n_tracers, n_no2, c_no2,            &
                     forward)

        IF (n_h2os > 0) THEN
          tracers(:,:,:,n_h2os) = tracers(:,:,:,n_h2o)
          q(1-off_x:row_length+off_x,1-off_y:rows+off_y,:) =            &
            tracers(:,:,:,n_h2o)
          CALL transform(row_length, rows, model_levels, off_x,         &
                     off_y, tr_ukca, tracers, n_h_tracers, h_tracers,   &
                     c_h_tracers, s_h_tracers, n_h2os, c_h2o/2.,        &
                     forward)
        END IF
      END IF
      IF (lhook) CALL dr_hook('UKCA_TRANSFORM_HALOGEN',zhook_out,zhook_handle)
      RETURN

      CONTAINS

!====================================================================

      SUBROUTINE transform(row_length, rows, model_levels, off_x,       &
                              off_y, tr_ukca, tracers, n_tracers,       &
                              position,                                 &
                              c_tracers, s_tracers, major, c_major,     &
                              forward)

      USE ereport_mod, ONLY : ereport
      USE Control_Max_Sizes
      IMPLICIT NONE

! I/O variables
      INTEGER, INTENT(IN) :: row_length
      INTEGER, INTENT(IN) :: rows
      INTEGER, INTENT(IN) :: model_levels
      INTEGER, INTENT(IN) :: off_x
      INTEGER, INTENT(IN) :: off_y
      INTEGER, INTENT(IN) :: tr_ukca
      INTEGER, INTENT(IN) :: n_tracers

      INTEGER, INTENT(IN) :: major
      REAL, INTENT(IN) :: c_major
      LOGICAL, INTENT(IN) :: forward

      INTEGER, INTENT(IN) :: position(n_tracers)
      REAL, INTENT(IN) :: c_tracers(n_tracers)
      INTEGER, INTENT(IN) :: s_tracers(n_tracers)
      REAL, INTENT(INOUT) :: tracers(1-off_x:row_length+off_x,          &
                                     1-off_y:rows+off_y, model_levels,  &
                                     tr_ukca)

! Local variables

! Status of compounds
      INTEGER, PARAMETER :: unrestricted   = 0
      INTEGER, PARAMETER :: source_gas     = 1
      INTEGER, PARAMETER :: mixed_compound = 2

      REAL, PARAMETER :: rafeps = 1e-150

      INTEGER :: i

      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: total_all
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: total_source
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: total_nochange
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('TRANSFORM',zhook_in,zhook_handle)

      ALLOCATE(total_all(1-off_x:row_length+off_x,1-off_y:rows+off_y,   &
               model_levels))
      ALLOCATE(total_source(1-off_x:row_length+off_x,                   &
               1-off_y:rows+off_y,model_levels))
      ALLOCATE(total_nochange(1-off_x:row_length+off_x,                 &
               1-off_y:rows+off_y,model_levels))

! Sum up fractional contributions to total.
! Make sure they do not sum up to more than 1. If they do, rescale
! chlorine tracers.

      total_all = 0.
      total_source = 0.
      total_nochange = 0.

      DO i=1,n_tracers
        SELECT CASE (s_tracers(i))
          CASE (unrestricted)
            total_all = total_all +                                     &
              tracers(:,:,:,position(i)) / c_tracers(i)
          CASE (source_gas)
            total_source = total_source +                               &
              tracers(:,:,:,position(i)) / c_tracers(i)
          CASE (mixed_compound)
            total_nochange = total_nochange +                           &
              tracers(:,:,:,position(i)) / c_tracers(i)
        END SELECT
      END DO
      IF (.NOT. forward) THEN
        tracers(:,:,:,major) = tracers(:,:,:,major) +                   &
           c_major * (total_all + total_source + total_nochange)
      ELSE
! Subtract non-changing compounds and adjust where necessary.
        tracers(:,:,:,major) = tracers(:,:,:,major) -                   &
           c_major * total_nochange
        DO i=1,n_tracers
          IF (s_tracers(i) < mixed_compound) THEN
            WHERE (tracers(:,:,:,major) < rafeps)                       &
              tracers(:,:,:,position(i)) = rafeps
          END IF
        END DO
        WHERE (tracers(:,:,:,major) < rafeps)
          tracers(:,:,:,major) = rafeps
          total_source = rafeps
          total_all = rafeps
        ENDWHERE
! Subtract source gas compounds and rescale where necessary
        total_nochange = total_source * c_major /                       &
          MAX(tracers(:,:,:,major),rafeps)
        DO i=1,n_tracers
          SELECT CASE (s_tracers(i))
            CASE (unrestricted)
              WHERE (total_nochange > 1.)                               &
                tracers(:,:,:,position(i)) = rafeps
            CASE (source_gas)
              WHERE (total_nochange > 1.)                               &
                tracers(:,:,:,position(i)) =                            &
                  tracers(:,:,:,position(i)) / total_nochange
          END SELECT
        END DO
        WHERE (total_nochange > 1.) total_all = rafeps
        tracers(:,:,:,major) = MAX(tracers(:,:,:,major) -               &
          c_major * total_source,rafeps)
! Rescale remaining gases
        total_nochange = total_all * c_major /                          &
           MAX(tracers(:,:,:,major),rafeps)
        DO i=1,n_tracers
          IF (s_tracers(i) == unrestricted)                             &
            WHERE (total_nochange > 1.)                                 &
              tracers(:,:,:,position(i)) =                              &
              tracers(:,:,:,position(i)) / total_nochange
        END DO
        tracers(:,:,:,major) = MAX(tracers(:,:,:,major) -               &
          c_major * total_all,rafeps)
      END IF
      DEALLOCATE (total_all)
      DEALLOCATE (total_source)
      DEALLOCATE (total_nochange)

      IF (lhook) CALL dr_hook('TRANSFORM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE transform

      END SUBROUTINE ukca_transform_halogen
