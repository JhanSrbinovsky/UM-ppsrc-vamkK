! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
! *****************************COPYRIGHT*******************************
!
!  Description:  West scheme for cloud droplet activation.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office.
!  See: www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
!
!  West scheme for cloud droplet activation.
!  based on *activ_box* activation box model, but no longer a box model!
!
!  Purpose: calls aero_activ_abdulrazzak_ghan which calculates 
!  -------- number concentration of aerosol particles which become 
!           "activated" into cloud droplets
!
! THINGS AERO_ACTIV_ABDULRAZZAK_GHAN NEEDS:
!
! 
! kbdim, 
! klev, 
!  
! zcdncact(kbdim,klev,nmodes),number concentration of activated particles: INITIALLY 0
! zesw(kbdim,klev)         saturation water vapour pressure
! zrho(kbdim,klev)         air density[kg m-3]
! zn(kbdim,klev,nmodes)    aerosol number concentration for each mode [m-3]
! zxtm1(kbdim,klev,nmodes,ncp) tracer mass mixing ratio [kg kg-1]
! ztm1(kbdim,klev),        temperature[T]   
! zapm1(kbdim,klev),       pressure [K]
! zqm1(kbdim,klev),        specific humidity[kg kg-1]
! ztkem1(kbdim,klev),      turbulent kinetic energy [m2 s-2]: SET TO 0
! zvervel(kbdim,klev),     prescribe this
! zrdry(kbdim,klev,nmod),  dry count median radius (m)
! zsigma(nmod)             Geometric dispersion (sigma_g)
! --------------------------------------------------------------------------

      SUBROUTINE UKCA_ACTIVATE(row_length, rows, klev,                  &
                          wet_levels, bl_levels, kbdim,                 &
                          n_mode_tracers,                               &
                          n_mode_diags, n_chem_diags,                   &
                          tr_index,                                     &
                          mode_tracers,                                 &
                          mode_diags, chem_diags,                       &
                          p_theta_levels,                               &
                          t_theta_levels,                               &
                          q,                                            &
                          qsvp,                                         &
                          bl_tke,                                       &
                          vertvel,                                      &
                          liq_cloud_frac,                               &
                          cloud_liq_water,                              &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                          len_stashwork,STASHwork                       & 
                          )


 
      USE UKCA_CONSTANTS,   ONLY: zboltz,  &  ! Boltzmann's constant
                                  ra          ! gas constant for dry air

      USE UKCA_MODE_SETUP,  ONLY: nmodes, ncp, mode, component,         &
                                  modesol,                              &
                                  ii_nd,   &  ! tracer index for Number
                                  ii_md,   &  ! tracer index for MMR
                                  CP_SU,   &  ! Index to store SO4    cpt
                                  CP_BC,   &  ! Index to store BC     cpt
                                  CP_OC,   &  ! Index to store 1st OC cpt
                                  CP_CL,   &  ! Index to store NaCl   cpt
                                  CP_DU,   &  ! Index to store dust   cpt
                                  CP_SO  ! Index to store 2nd OC cpt


      USE UKCA_D1_DEFS,     ONLY: ukcaD1codes, item1_mode_diags,        &
                                  nmax_mode_diags, nm_spec,imode_first, &
                                  Nukca_D1items, mode_diag_sect,        &
                                  L_ukca_mode_diags, icd_cdnc, icd_cdnc3
      USE ukca_option_mod,  ONLY: L_ukca_sfix
      USE UKCA_ACTIV_MOD,   ONLY: activmklin, activmkskew
      USE ereport_mod,      ONLY: ereport
      USE PrintStatus_mod
      USE Control_Max_Sizes
      USE UM_ParVars,       ONLY: at_extremity
      USE yomhook,          ONLY: lhook, dr_hook
      USE parkind1,         ONLY: jprb, jpim
      USE Submodel_Mod
      IMPLICIT NONE
   
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end

! Interface:

      INTEGER, INTENT(IN) :: row_length        ! UM array dimension
      INTEGER, INTENT(IN) :: rows              ! UM array dimension
      INTEGER, INTENT(IN) :: klev              ! = model_levels
      INTEGER, INTENT(IN) :: wet_levels        ! UM array dimension
      INTEGER, INTENT(IN) :: bl_levels         ! UM array dimension
      INTEGER, INTENT(IN) :: kbdim             ! = theta_field_size=row_length*rows
      INTEGER, INTENT(IN) :: n_mode_tracers    ! # of mode tracers
      INTEGER, INTENT(IN) :: n_mode_diags      ! # of mode diagnostics
      INTEGER, INTENT(IN) :: n_chem_diags      ! # of mode diagnostics
      INTEGER, INTENT(IN) :: tr_index(n_mode_tracers) ! tracer item index

! aerosol tracers mass mixing ratio
      REAL, INTENT(IN) :: mode_tracers(row_length,rows,klev,            &
                                       n_mode_tracers)
! Pressure
      REAL, INTENT(IN) :: p_theta_levels(row_length,rows,klev)
! Temperature
      REAL, INTENT(IN) :: t_theta_levels(row_length,rows,klev)
! Specific humidity
      REAL, INTENT(IN) :: q(row_length,rows,wet_levels)
! Saturated
      REAL, INTENT(IN) :: qsvp(row_length,rows,wet_levels)
! Turbulent kinetic energy from BL scheme
      REAL, INTENT(IN) :: bl_tke(row_length,rows,bl_levels)
! w component of wind (from Section 0 Item 150)
      REAL, INTENT(IN) :: vertvel(row_length,rows,klev)
! liquid cloud fraction by volume (from Section 0 Item 267)
      REAL, INTENT(IN) :: liq_cloud_frac(row_length,rows,klev)
! cloud liquid water (from Section 4 Item 205)
      REAL, INTENT(IN) :: cloud_liq_water(row_length,rows,klev)
! MODE diagnostics array
      REAL, INTENT(INOUT) :: mode_diags(row_length,rows,klev,n_mode_diags)
! Chemistry diagnostics array (not transported)
      REAL, INTENT(INOUT) :: chem_diags(row_length,rows,klev,n_chem_diags)
! Diagnostic variables
      INTEGER, INTENT(IN) :: len_stashwork
      REAL, INTENT(INOUT) :: STASHwork(len_stashwork)

! Local variables:
      INTEGER, SAVE :: ldry(nmodes)  ! Address of dry radius in mode_diags
      INTEGER :: i_nuc=200           ! item number for nucl mode dry radius
      INTEGER :: itra
      INTEGER :: imode               ! loop index to modes
      INTEGER :: icp                 ! loop index to modes
      INTEGER :: k                   ! loop index
      INTEGER :: l                   ! loop index
      INTEGER :: i                   ! loop index
      INTEGER :: j                   ! loop index
      INTEGER :: m                   ! loop index
      INTEGER :: N                   ! loop index

 ! RW variables for ukca_abdulrazzak_ghan
            
      LOGICAL :: first=.TRUE.
                
      REAL :: ztm1(kbdim,klev)        ! temperature[K]
      REAL :: zapm1(kbdim,klev)       ! pressure [Pa=N m-2]
      REAL :: zrho(kbdim,klev)        ! air density (kg/m3)
      REAL :: zaird(kbdim, klev)      ! number concentration of air molecules (/m3)
      REAL :: zqm1(kbdim,klev)        ! specific humidity[kg kg-1]
      REAL :: zesw(kbdim,klev)        ! saturation water vapour pressure
      REAL :: ztkem1(kbdim,klev)      ! turbulent kinetic energy [m2 s-2]
      REAL :: zvervel(kbdim,klev)     ! large scale vertical velocity [m s-1]
      REAL :: zsmax(kbdim,klev)       ! maximum supersaturation [fraction]
      REAL :: zcdncact(kbdim,klev)    ! no. concentration of activated particles> [m-3]
      REAL :: zcdncact3(kbdim,klev)   ! <Number_activated**-1/3> [m]
      REAL :: zcdncactm(kbdim,klev,nmodes)  ! number concentration of activated
!                                           !  particles by mode [m-3]
      REAL :: zrdry(kbdim,klev,nmodes)      ! dry count median radius [m]
      REAL :: zn(kbdim,klev,nmodes)         ! aerosol no. concentration for modes [m-3]
      REAL :: zxtm1(kbdim,klev,nmodes,ncp)  ! component mass mixing ratio [kg kg-1]

      ! updraft velocity array for pdf
      REAL :: zsigmaw(kbdim,klev) ! std dev of gaussian distn of w
      REAL, PARAMETER :: sigwmin=0.1 ! define min value of sigmaw
      INTEGER, PARAMETER :: nwbins=20
       !   INTEGER, PARAMETER :: nwbins=1
      REAL, PARAMETER :: skewness=0.0   ! define fixed value of skewness of updraft
!                                       !  distribution
      REAL :: zvervel_min(kbdim,klev)   ! min limit of updraft vel pdf [m s-1]
      REAL :: zvervel_max(kbdim,klev)   ! max limit of updraft vel pdf [m s-1]
      REAL :: zwarr(kbdim,klev,nwbins)  ! lin array of vert vel [m s-1]
      REAL :: zwpdf(kbdim,klev,nwbins)  ! lin array of pdf of w
      REAL :: zwbin(kbdim,klev,nwbins)  ! w bin width [m s-1]
      REAL :: zwchar(kbdim,klev)        ! calculated characteristic updraught [m s-1]
      REAL :: zwbar(kbdim,klev)         ! reshaped large scale vertical velocity [m s-1]
      REAL :: zwalpha(kbdim,klev)       ! skewness of vertical velocity distribution 

! RW local reshaped fields for output

      REAL :: n_activated(row_length,rows,klev,nmodes) ! number concentration of
!                                                      ! activated particles by mode[m-3]
      REAL :: n_activ_sum(row_length,rows,klev)        ! total number concentration of
!                                                      !  activated particles [m-3]
      REAL :: n_activ_sum3(row_length,rows,klev)       ! <Number_activated**-1/3> [m]
      REAL :: smaxpc(row_length,rows,klev)             ! maximum supersaturation [%]
      REAL :: wchar(row_length,rows,klev)   ! calculated characteristic updraught [m s-1]
      REAL :: sigw(row_length,rows,klev)    ! spread of pdf of updraught [m s-1]
      REAL :: cdncwt(row_length,rows,klev)  ! weighted cdnc =
!                                           !         total cdnc * liq_cloud_frac [m-3]
      REAL :: cldflag(row_length,rows,klev) ! cloud flag=1 if cloud in grid box, else 0
      REAL :: cldbase(row_length,rows,klev) ! cloud base=1 and no cloud below, else 0
      REAL :: cdncflag(row_length,rows,klev) ! weighted cdnc = total cdnc * cldflg [m-3]
      REAL :: cdnc3flag(row_length,rows,klev) ! weighted <cdnc^-1/3> =
!                                             !  <cdnc^-1/3> * cldflg [m-3]
      REAL :: smaxflag(row_length,rows,klev)  ! weighted Smax = Smax * cldflg [%]
      REAL :: wcharflag(row_length,rows,klev) ! weighted wchar = wchar * cldflg [m s-1]
      REAL :: sigwflag(row_length,rows,klev)  ! weighted sigmaw = sigw * cldflg [m s-1]
      REAL :: tkeflag(row_length,rows,klev)   ! weighted BL TKE = sigw * cldflg [m s-1]

      CHARACTER(LEN=60) :: cmessage           ! error message

! Create an array of fixed supersaturations to run the code at:
      INTEGER,PARAMETER :: nsfix=19           ! number of elements in fixed-S array
      INTEGER,PARAMETER :: twonsfix=38        ! 19*2=38
      INTEGER :: isfix                        ! for looping

! array of values of fixed-s (*100 to get in %)
      REAL    :: zsfix(nsfix)=(/ 0.0002, 0.0004, 0.0006, 0.0008, & 
                                 0.001,  0.0016, 0.002,  0.0023, &
                                 0.003,  0.0033, 0.0038, 0.004,  &
                                 0.005,  0.006,  0.0075, 0.008,  &
                                 0.0085, 0.01,   0.012           /)        

      INTEGER :: icode              ! Return code  =0 Normal exit  >1 Error
      INTEGER :: sect               ! STASH section,
      INTEGER :: item               ! STASH item no.s
      INTEGER :: im_index           ! internal model index
      INTEGER :: im                 ! loop index
      REAL    :: n_ccn_2d(row_length,rows) !local
      INTEGER, PARAMETER :: ccnlev=1 ! sets the level on which to call fixeds
      REAL    :: zccn(kbdim,nsfix)
      REAL    :: n_ccn(row_length,rows,klev)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('UKCA_ACTIVATE',zhook_in,zhook_handle)

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)


      IF ( first .AND. PrintStatus >= PrStatus_Oper ) THEN
        WRITE(6,*) 'ROSALIND some INPUTS to ukca_activate:'
        WRITE(6,*) 'ROSALIND nmodes, ncp ', nmodes, ncp
        WRITE(6,*) 'ROSALIND mode ', mode
        WRITE(6,*) 'ROSALIND modesol',  modesol
      END IF

! Derived quantities
! ==================
      DO k=1,klev
         ztm1(:,k) = RESHAPE(t_theta_levels(:,:,k),(/kbdim/))
         zapm1(:,k) = RESHAPE(p_theta_levels(:,:,k),(/kbdim/))
         zrho(:,k) = zapm1(:,k)/(ztm1(:,k)*ra) ! now in kg/m3 
!(RA is specific gas constant for dry air in J kg^-1 K^-1
         zaird(:,k) = zapm1(:,k)/(ztm1(:,k)*zboltz) ! number air molecules /m3
         zqm1(:,k) = RESHAPE(q(:,:,k),(/kbdim/))
         zesw(:,k) = RESHAPE(qsvp(:,:,k),(/kbdim/))
      END DO
 
! Aerosol Tracers
! ===============
! Loop through the modes
      DO imode=1,nmodes
        IF(mode(imode)) THEN
         itra=ii_nd(imode)
! Set No Density (particles per m3) from aerosol tracer array
! = number particles per molecule air * number concentration of air molecules
         DO k=1,klev
            zn(:,k,imode)=RESHAPE(mode_tracers(:,:,k,itra),    &
            (/kbdim/))*zaird(:,k)
         END DO
      
! Loop through the components
         DO icp=1,ncp
            IF(component(imode,icp)) THEN
              itra=ii_md(imode,icp)
              IF(itra == 0) THEN
                write(6,*) '***** ERROR: MD_CP specified by COMPONENT'
                write(6,*) '*****        doesnt exist in II_MD'
                write(6,*) 'IMODE,ICP,II_MD=',imode,icp,ii_md(imode,icp)
                write(6,*) 'COMPONENT(IMODE,ICP)=',component(imode,icp)
                write(6,*) 'II_MD(IMODE,ICP)=',ii_md(imode,icp)
                cmessage=' MD_CP specified by COMPONENT not in II_MD'
                icode = 1
                CALL EREPORT('UKCA_ACTIVATE',icode,cmessage)
              END IF !itra
! Set mass mixing ratios from aerosol tracer array
              DO k=1,klev
                zxtm1(:,k,imode,icp) = RESHAPE(mode_tracers(:,:,k,itra) &
                                               ,(/kbdim/))
              END DO
            END IF  ! component
          END DO  ! loop over cpts
        END IF  ! mode(imode)
      END DO  ! loop over modes


      IF (first) THEN
! Check that the dry radius diagnostics have been requested and find
! the address in the mode_diags array.
        ldry(:)=IMDI
        k=0
        DO l=1,nmax_mode_diags
          IF (UkcaD1codes(imode_first+l-1)%item /= IMDI) THEN
            k=k+1
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc)   ldry(1)=k
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc+1) ldry(2)=k
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc+2) ldry(3)=k
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc+3) ldry(4)=k
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc+4) ldry(5)=k
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc+5) ldry(6)=k
            IF (UkcaD1codes(imode_first+l-1)%item ==                    &
                item1_mode_diags + i_nuc+6) ldry(7)=k
          END IF
        END DO
        DO imode=1,nmodes
           IF (mode(imode)) THEN
            IF (ldry(imode) == IMDI) THEN
               cmessage=' Address for dry radius in mode_diags not found'
               write(6,*) cmessage
               write(6,*) 'IMODE: ',imode
               write(6,*) 'Check that the dry radius diagnostics are set'
               icode = 1
               CALL EREPORT('UKCA_ACTIVATE',icode,cmessage)
            END IF
          END IF
        END DO
      END IF            ! first
      
      DO imode=1,nmodes
         IF (mode(imode)) THEN
            DO k=1,klev
               zrdry(:,k,imode)=RESHAPE(mode_diags(:,:,k,ldry(1)        &
               +imode-1),(/kbdim/))/2.0
            END DO
         END IF
      END DO


! Updraught velocity pdf
! =====================
      sigw(:,:,:)=0.0E0         ! initialise to 0
! Fix zsigmaw (std dev of updraught vel pdf) :
     ! zsigmaw(:,:)=(bl_tke(:,:,9))**0.5 ! m/s

! Define width of updraught velocity pdf based on BL TKE     
! (Factor of 2/3 due to assumption that TKE is isotropic)
      DO k=1,klev
         IF ( k <= bl_levels ) THEN
            sigw(:,:,k)=MAX(sigwmin, ((2.0/3.0)*bl_tke(:,:,k))**0.5) ! m/s
         ELSE
            sigw(:,:,k)=sigwmin
         END IF
!         sigw(:,:,k)=0.1 ! m/s [based on w_obs database]
         zsigmaw(:,k)=RESHAPE(sigw(:,:,k),(/kbdim/))
         zwbar(:,k)=RESHAPE(vertvel(:,:,k),(/kbdim/))
      END DO
        
      zvervel_min(:,:)=0.0 ! m/s
      zvervel_max(:,:)=4.0*zsigmaw(:,:) !m/s

! Set shape parameter (skewness) of skew-normal updraft distribution
      zwalpha(:,:)=skewness

! vertical velocity  
      IF (nwbins > 1) THEN
! generate pdf of vertical velocity  
        CALL activmklin(kbdim, klev, zvervel_min, zvervel_max, nwbins,  &
                        zwbin, zwarr)
        CALL activmkskew(kbdim, klev, zwarr, nwbins, zsigmaw(:,:),      &
                         zwbar(:,:), zwalpha(:,:), zwpdf)
      ELSE IF (nwbins == 1) THEN
! set fixed vertical velocity 
         zwbin(:,:,:)=1.0
         zwpdf(:,:,:)=1.0
         zwarr(:,:,1)=zvervel_max(:,:)
      ELSE 
        WRITE (6,*) 'RW: No smoke without fire. Specify nwbins>=1'
        cmessage = ' No smoke without fire. Specify nwbins>=1'
        icode = 1
        CALL EREPORT('UKCA_ACTIVATE',icode,cmessage)
      END IF

! Initialisation:
      zcdncactm(:,:,:)=0.0E0     ! initialise to 0
      zcdncact(:,:)=0.0E0        ! initialise to 0
      zcdncact3(:,:)=0.0E0       ! initialise to 0
      n_activ_sum(:,:,:)=0.0E0   ! initialise to 0
      n_activ_sum3(:,:,:)=0.0E0  ! initialise to 0
      n_activated(:,:,:,:)=0.0E0 ! initialise to 0
      n_ccn(:,:,:)=0.0E0         ! initialise to 0
      zccn(:,:)=0.0E0            ! initialise to 0
      zsmax(:,:)=0.0E0           ! initialise to 0
      smaxpc(:,:,:)=0.0E0        ! initialise to 0
      zwchar(:,:)=0.0E0          ! initialise to 0
      wchar(:,:,:)=0.0E0         ! initialise to 0
      cdncwt(:,:,:)=0.0E0        ! initialise to 0
      cldflag(:,:,:)=0.0E0       ! initialise to 0
      cldbase(:,:,:)=0.0E0       ! initialise to 0
      cdncflag(:,:,:)=0.0E0      ! initialise to 0
      cdnc3flag(:,:,:)=0.0E0     ! initialise to 0   
      smaxflag(:,:,:)=0.0E0      ! initialise to 0
      wcharflag(:,:,:)=0.0E0     ! initialise to 0
      sigwflag(:,:,:)=0.0E0      ! initialise to 0
      tkeflag(:,:,:)=0.0E0       ! initialise to 0
      
      IF ( PrintStatus == PrStatus_Diag ) THEN
        WRITE(6,*) 'ROSALIND: initialised values in ukca_activate:'
        WRITE(6,*) 'ROSALIND bl_tke',minval(bl_tke(:,:,:)),             &
                    maxval(bl_tke(:,:,:))
        WRITE(6,*) 'ROSALIND vertvel',minval(vertvel(:,:,:)),           &
                    maxval(vertvel(:,:,:))
        WRITE(6,*) 'ROSALIND liq_cloud_frac',                           &
                    minval(liq_cloud_frac(:,:,:)),                      &
                    maxval(liq_cloud_frac(:,:,:))
        WRITE(6,*) 'ROSALIND cloud_liq_water',                          &
                    minval(cloud_liq_water(:,:,:)),                     &
                    maxval(cloud_liq_water(:,:,:))
        WRITE(6,*) 'ROSALIND zcdncactm',minval(zcdncactm(:,:,:)),       &
                    maxval(zcdncactm(:,:,:))
        WRITE(6,*) 'ROSALIND zcdncact',minval(zcdncact(:,:)),           &
                    maxval(zcdncact(:,:))
        WRITE(6,*) 'ROSALIND zwchar',minval(zwchar(:,:)),               &
                    maxval(zwchar(:,:))
  !    WRITE(6,*) 'ROSALIND zesw', minval(zesw(:,:)), maxval(zesw(:,:))
  !    WRITE(6,*) 'ROSALIND zrho', minval(zrho(:,:)), maxval(zrho(:,:))

! write out min and max values of aerosol number concentration for each mode
        WRITE(6,*) 'ROSALIND zn(1)', minval(zn(:,:,1)), maxval(zn(:,:,1))
        WRITE(6,*) 'ROSALIND zn(2)', minval(zn(:,:,2)), maxval(zn(:,:,2))
        WRITE(6,*) 'ROSALIND zn(3)', minval(zn(:,:,3)), maxval(zn(:,:,3))
        WRITE(6,*) 'ROSALIND zn(4)', minval(zn(:,:,4)), maxval(zn(:,:,4))
        WRITE(6,*) 'ROSALIND zn(5)', minval(zn(:,:,5)), maxval(zn(:,:,5))
!      WRITE(6,*) 'ROSALIND zn(6)', minval(zn(:,:,6)), maxval(zn(:,:,6))
!      WRITE(6,*) 'ROSALIND zn(7)', minval(zn(:,:,7)), maxval(zn(:,:,7))
        WRITE(6,*) 'ROSALIND ztm1', minval(ztm1(:,:)), maxval(ztm1(:,:))      
        WRITE(6,*) 'ROSALIND zapm1', minval(zapm1(:,:)),maxval(zapm1(:,:))
        WRITE(6,*) 'ROSALIND zqm1', minval(zqm1(:,:)), maxval(zqm1(:,:))             
        WRITE(6,*) 'ROSALIND zvervel_max', minval(zvervel_max(:,:)),    &
                    maxval(zvervel_max(:,:))
        WRITE(6,*) 'ROSALIND zrdry', minval(zrdry(:,:,1:5)),            &
                    maxval(zrdry(:,:,1:5))
        WRITE(6,*) 'ROSALIND zwpdf', minval(zwpdf(:,:,:)),              &
                    maxval(zwpdf(:,:,:)) 
        WRITE(6,*) 'ROSALIND zwarr', minval(zwarr(:,:,:)),              &
                    maxval(zwarr(:,:,:))
        WRITE(6,*) 'ROSALIND zwbin', minval(zwbin(:,:,:)),              &
                    maxval(zwbin(:,:,:))
      END IF     ! Prstatus
      
! DEPENDS ON: ukca_abdulrazzak_ghan
      CALL UKCA_ABDULRAZZAK_GHAN(kbdim,   klev,                         &
                                 zesw,    zrho,                         &
                                 zn,      zxtm1,                        &
                                 ztm1,    zapm1,                        &
                                 zqm1,    zrdry,                        &
                                 nwbins,  zwarr,                        &
                                 zwpdf,   zwbin,                        &
                                 zsmax,   zwchar,                       &
                                 zcdncactm, zcdncact,                   &
                                 zcdncact3 )
      
      IF ( PrintStatus == PrStatus_Diag ) THEN
       WRITE(6,*) 'ROSALIND after ARG'
       WRITE(6,*) 'ROSALIND zcdncactm',minval(zcdncactm(:,:,:)),        &
                   maxval(zcdncactm(:,:,:))
       WRITE(6,*) 'ROSALIND zcdncact',minval(zcdncact(:,:)),            &
                   maxval(zcdncact(:,:))
       WRITE(6,*) 'ROSALIND zwchar',minval(zwchar(:,:)),                &
                   maxval(zwchar(:,:))
      END IF
                                               
! reshape zcdncactm to n_activated
      DO imode=1, nmodes
         IF (mode(imode)) THEN
            DO k=1,klev
               n_activated(:,:,k,imode)=RESHAPE(zcdncactm(:,k,imode),   &
                                         (/row_length,rows/))
            END DO
         END IF
      END DO

      IF ( PrintStatus == PrStatus_Diag )                               &
        WRITE(6,*) 'ROSALIND n_activated',minval(n_activated(:,:,:,:)), &
                    maxval(n_activated(:,:,:,:))


!reshape zcdncact to n_activ_sum
      DO k=1,klev
        n_activ_sum(:,:,k)=RESHAPE(zcdncact(:,k),(/row_length,rows/))  
        n_activ_sum3(:,:,k)=RESHAPE(zcdncact3(:,k),(/row_length,rows/))
!reshape zsmax to smaxpc and convert to %
        smaxpc(:,:,k)=100.0*(RESHAPE(zsmax(:,k),(/row_length,rows/)))
!reshape zwchar to wchar
        wchar(:,:,k)=RESHAPE(zwchar(:,k),(/row_length,rows/))
      END DO
      
      IF ( PrintStatus == PrStatus_Diag ) THEN
        WRITE(6,*) 'ROSALIND n_activ_sum',minval(n_activ_sum(:,:,:)),   &
                    maxval(n_activ_sum(:,:,:))
        WRITE(6,*) 'ROSALIND n_activ_sum3',minval(n_activ_sum3(:,:,:)), &
                    maxval(n_activ_sum3(:,:,:))
        WRITE(6,*) 'ROSALIND zsmax', minval(zsmax(:,:)),                &
                    maxval(zsmax(:,:))     
        WRITE(6,*) 'ROSALIND smaxpc', minval(smaxpc(:,:,:)),            &
                    maxval(smaxpc(:,:,:)) 
        WRITE(6,*) 'ROSALIND zsigmaw', minval(zsigmaw(:,:)),            &
                    maxval(zsigmaw(:,:))
      END IF

! calculate weighted CDNC (total CDNC weighted by liquid cloud fraction)
!  make a cloud flag diagnostic to tell us when and where there is > 1% cloud fraction
!  also make an unweighted CDNC diag, zero outside cloud.
      DO i=1, row_length
         DO j=1, rows
! For bottom level set CDNC if liquid cloud fraction > 1% and LWC > 1.0E-6 kg kg^-1
!           IF (liq_cloud_frac(i,j,1) > 0.01 .AND. cloud_liq_water(i,j,1) > 1.0E-06 ) THEN
! For bottom level set CDNC if liquid cloud fraction > 0 and LWC > 0 kg kg^-1
           IF (liq_cloud_frac(i,j,1) > 0.0 .AND.                        &
               cloud_liq_water(i,j,1) > 0.0 ) THEN
             cdncwt(i,j,1)=liq_cloud_frac(i,j,1)*n_activ_sum(i,j,1)
             cldflag(i,j,1)=1.0
             cldbase(i,j,1)=1.0
             cdncflag(i,j,1)=n_activ_sum(i,j,1)
             cdnc3flag(i,j,1)=n_activ_sum3(i,j,1)
             smaxflag(i,j,1)=smaxpc(i,j,1)
             wcharflag(i,j,1)=wchar(i,j,1)
             sigwflag(i,j,1)=sigw(i,j,1)
             tkeflag(i,j,1)=bl_tke(i,j,1) 
           END IF
           DO k=2, klev
! Define cloudy if liquid cloud fraction > 1% and LWC > 1.0E-6 kg kg^-1
!            IF (liq_cloud_frac(i,j,k) > 0.01 .AND. cloud_liq_water(i,j,k) > 1.0E-06 ) THEN
! Define cloudy if liquid cloud fraction > 0 and LWC > 0 kg kg^-1
             IF (liq_cloud_frac(i,j,k) > 0.0 .AND.                      &
                 cloud_liq_water(i,j,k) > 0.0 ) THEN
               cldflag(i,j,k)=1.0
! If layer below (i.e. k-1) is cloudy then set CDNC in this level (k) to number activated at cloud base 
               IF ( cldflag(i,j,k-1) == 1.0 ) THEN
                 cdncwt(i,j,k)=liq_cloud_frac(i,j,k)*cdncflag(i,j,k-1)
                 cdncflag(i,j,k)=cdncflag(i,j,k-1)
                 cdnc3flag(i,j,k)=cdnc3flag(i,j,k-1)
                 smaxflag(i,j,k)=smaxpc(i,j,k-1) 
                 wcharflag(i,j,k)=wchar(i,j,k-1)
                 sigwflag(i,j,k)=sigw(i,j,k-1)
! If layer below (k-1) is cloudy, then this layer(k) is not the cloud base, so
                 cldbase(i,j,k)=0.0
                 IF ( k <= bl_levels ) tkeflag(i,j,k)=bl_tke(i,j,k-1) 
               ELSE 
                 cdncwt(i,j,k)=liq_cloud_frac(i,j,k)*n_activ_sum(i,j,k)
                 cdncflag(i,j,k)=n_activ_sum(i,j,k)
                 cdnc3flag(i,j,k)=n_activ_sum3(i,j,k)
                 smaxflag(i,j,k)=smaxpc(i,j,k)
                 wcharflag(i,j,k)=wchar(i,j,k)
                 sigwflag(i,j,k)=sigw(i,j,k)
! But if layer below(k-1) is not cloudy, then this layer(k) is the cloud base 
                 cldbase(i,j,k)=1.0 
                 IF ( k <= bl_levels ) tkeflag(i,j,k)=bl_tke(i,j,k) 
               END IF
             END IF    ! liq_cloud_frac > 0 etc        
           END DO ! k
         ENDDO ! j
      ENDDO ! i 

! write CDNC and <CDNC^-1/3> to D1
      chem_diags(:,:,:,icd_cdnc)=cdncflag(:,:,:)
      chem_diags(:,:,:,icd_cdnc3)=cdnc3flag(:,:,:)

! Calculate CCN at fixed supersaturation:
      IF(L_ukca_sfix) THEN

! DEPENDS ON: ukca_fixeds
         CALL UKCA_FIXEDS(kbdim,                                        &
                          zccn,    zn(:,ccnlev,:),                      &
                          zxtm1(:,ccnlev,:,:),                          &
                          ztm1(:,ccnlev),                               &
                          zrdry(:,ccnlev,:),                            &
                          zsfix,   nsfix )

!reshape zccn to n_ccn: 
         DO isfix=1, nsfix
            n_ccn(:,:,isfix)=RESHAPE(zccn(:,isfix),                     &
                 (/row_length,rows/))
         END DO

      END IF !L_ukca_sfix

 
!  Write 3_D diagnostics to mode_diags array
!  N.B. L_ukca_mode_diags is set whenever a STASH request for a relevant
!  item is found, and the ukcaD1codes(N)%item is then set, otherwise it is IMDI    
      IF (L_ukca_mode_diags) THEN ! fill 3D array
         m=0
         DO N=1,Nukca_D1items
           IF (ukcaD1codes(N)%section == MODE_diag_sect .AND.           &
               ukcaD1codes(N)%item >= item1_mode_diags .AND.            &
               ukcaD1codes(N)%item <= item1_mode_diags+                 &
                                   nmax_mode_diags-1 .AND.              &
               ukcaD1codes(N)%item /= IMDI) THEN
             m=m+1

             SELECT CASE(UkcaD1codes(N)%item)
             CASE(item1_mode_diags+268) ! CDNC - NUC(sol)
               mode_diags(:,:,:,m)=n_activated(:,:,:,1)
             CASE(item1_mode_diags+269) ! CDNC - AIT(sol)
               mode_diags(:,:,:,m)=n_activated(:,:,:,2)
             CASE(item1_mode_diags+270)     ! CDNC - ACC(sol)
               mode_diags(:,:,:,m)=n_activated(:,:,:,3)
             CASE(item1_mode_diags+271)     ! CDNC - COR(sol)
               mode_diags(:,:,:,m)=n_activated(:,:,:,4)
             CASE(item1_mode_diags+272)     ! Smax %
               mode_diags(:,:,:,m)=smaxpc(:,:,:)
             CASE(item1_mode_diags+273)     ! cloud base
               mode_diags(:,:,:,m)=cldbase(:,:,:)
             CASE(item1_mode_diags+274)     ! sig_w [m/s]
               mode_diags(:,:,:,m)=sigw(:,:,:)
             CASE(item1_mode_diags+275)     ! liq cloud frac
               mode_diags(:,:,:,m)=liq_cloud_frac(:,:,:)
             CASE(item1_mode_diags+276)     ! total CDNC * liq cloud frac [m-3]
               mode_diags(:,:,:,m)=cdncwt(:,:,:)
             CASE(item1_mode_diags+277)     ! cloud flag
               mode_diags(:,:,:,m)=cldflag(:,:,:)
             CASE(item1_mode_diags+278)     ! total CDNC * cloud flag [m-3]
               mode_diags(:,:,:,m)=cdncflag(:,:,:)
             CASE(item1_mode_diags+279)     ! Smaxpc * cloud flag [%]
               mode_diags(:,:,:,m)=smaxflag(:,:,:)
             CASE(item1_mode_diags+280)     ! Wchar * cloud flag [m s-1]
               mode_diags(:,:,:,m)=wcharflag(:,:,:)
             CASE(item1_mode_diags+281)     ! Sigw * cloud flag [m s-1]
               mode_diags(:,:,:,m)=sigwflag(:,:,:)
             CASE(item1_mode_diags+282)     ! BL TKE * cloud flag [m2 s-2]
               mode_diags(:,:,:,m)=tkeflag(:,:,:)
             CASE(item1_mode_diags+283)     ! CCN at fixed S [m-3]
               IF(L_ukca_sfix) THEN
                  mode_diags(:,:,:,m)= n_ccn(:,:,:)
               ELSE 
                  mode_diags(:,:,:,m)=0.0
               END IF !L_ukca_sfix
            END SELECT

           END IF

         END DO
      END IF       ! L_ukca_mode_diags

      first=.FALSE.

      IF (lhook) CALL dr_hook('UKCA_ACTIVATE',zhook_out,zhook_handle)

      RETURN             
      END SUBROUTINE UKCA_ACTIVATE
