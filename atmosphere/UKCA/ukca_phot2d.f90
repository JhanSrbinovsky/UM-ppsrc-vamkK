! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module containing all routines relating to 2D photolysis
!  used in UKCA sub-model.
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      MODULE UKCA_PHOT2D

      USE ASAD_MOD,            ONLY: jpspj
      USE UKCA_CHEM_DEFS_MOD,  ONLY: ratj_t, ratj_defs
      USE yomhook,             ONLY: lhook, dr_hook
      USE parkind1,            ONLY: jprb, jpim
      
      USE conversions_mod,     ONLY: pi_over_180
      USE PrintStatus_mod,     ONLY: PrintStatus, PrStatus_Diag      
      USE ereport_mod,         ONLY: ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE
      SAVE
      PRIVATE


      INTEGER              :: myjppj  ! set = jppj below

      INTEGER,PARAMETER,PUBLIC        :: nolat=19   ! no of 2D lats
      INTEGER,PARAMETER,PUBLIC        :: nolev=17   ! no of 2D levs
      INTEGER,PARAMETER,PUBLIC        :: nlphot=51  ! no of 2D photol levs
      INTEGER,PARAMETER,PUBLIC        :: ntphot=3   ! no of times of day

!     Number of time intervals in data (74 corresponds to
!     a 5 day interval, with one extra set).

      INTEGER,PARAMETER               :: n_data_times=74  
      REAL,   PARAMETER               :: data_interval=5.0 ! interval in days

      LOGICAL, PARAMETER              :: L_optimised=.TRUE.  ! T for optimised read MONSooN        

!     Photolysis rates from 2-D model which have been interpolated onto 3-D

      REAL,ALLOCATABLE,DIMENSION(:,:,:,:),PUBLIC          :: pjin
 
! 2D model level pressures
      REAL,PUBLIC :: pr2d(nolev)                               

      PUBLIC UKCA_PHOTIN, UKCA_CURVE, UKCA_INPR2D

      CONTAINS

!-----------------------------------------------------------------------

      SUBROUTINE UKCA_PHOTIN(i_day_number,row_lengthda,tot_p_rows,     &
                      p_levelsda,p_fieldda,first_row, glon, glat,      & 
                      sinlat,pl,jppj) 

! Purpose: Subroutine to read in 2D photolysis rates, and interpolate on
!          3-d latitudes, heights. reconstruct daily curves (every 5 day
!          plus code to account for hour angle (longitude)
!
!          3 values are stored for each level and each latitude
!          symmetrical distribution plus zero values at dawn and
!          dusk gives 7 points altogether.
!          Based on PHOTIN.F from Cambridge TOMCAT model.

! dataset: 51(levels)x19(latitudes)x3(points)x74(every 5days)
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! ---------------------------------------------------------------------
!
        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE

      INTEGER, INTENT(IN) :: jppj
      INTEGER, INTENT(IN) :: i_day_number
      INTEGER, INTENT(IN) :: row_lengthda
      INTEGER, INTENT(IN) :: tot_p_rows
      INTEGER, INTENT(IN) :: p_levelsda
      INTEGER, INTENT(IN) :: p_fieldda
      INTEGER, INTENT(IN) :: first_row 
      INTEGER, INTENT(IN) :: glon 
      INTEGER, INTENT(IN) :: glat 

      REAL, INTENT(IN), DIMENSION(:)     :: SINLAT
      REAL, INTENT(IN), DIMENSION(:,:,:) :: PL

! Local variables

      INTEGER :: i,j 
      INTEGER :: idofy                         ! Day number
      INTEGER :: ipos
      LOGICAL,SAVE  :: first_call=.true.

      REAL :: fpos
      REAL :: pr2dj(nlphot)                    ! 2D photolysis level p
      REAL :: zmean_pl(tot_p_rows,p_levelsda)  ! zonal mean of 3D press
      REAL,dimension(:,:,:,:),allocatable    :: pjin2d

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_PHOTIN',zhook_in,zhook_handle)
      IF(first_call)  THEN
        myjppj = jppj
        call phot2d_allocate_memory (tot_P_ROWS,P_LEVELSDA,ntphot)
        first_call = .false.
      end if

      allocate(PJIN2D(NOLAT,NLPHOT,NTPHOT,MYJPPJ))

! Find nearest day in photolysis dataset - day 1 is 31. December

        idofy = i_day_number
        fpos  = idofy/data_interval + 1.0
        ipos  = nint(fpos)
        IF ((fpos-ipos*1.0) < 0.0) ipos = ipos-1

! Read in photolysis data

        IF (L_optimised) THEN
          CALL READ2D_OPT(ipos,fpos,pjin2d)
        ELSE
          CALL READ2D_ORIG(ipos,fpos,pjin2d)
        ENDIF

! Set up 2-D pressure arrays

        CALL UKCA_INPR2D(pr2d,pr2dj)

! Interpolate 2D photolysis rates onto 3-D levels and
! latitudes. Longitude comes later.

        CALL UKCA_CALC_ZMEAN_PRESS(row_lengthda, tot_p_rows,           &
                                 p_levelsda, glon, pl, zmean_pl) 
               
        CALL UKCA_INTERPJ(pjin2d,pr2dj,zmean_pl,row_lengthda,          &
                     tot_p_rows, p_levelsda,                           &
                     p_fieldda,sinlat,Pi_Over_180 )

        DEALLOCATE(PJIN2D)

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_PHOTIN',zhook_out,   &
                                zhook_handle)
        RETURN
      END SUBROUTINE UKCA_PHOTIN

!---------------------------------------------------------------------

        SUBROUTINE UKCA_CURVE(                                         &
                       pjinda,tloc,dayl,p_field,p_levels,              &
                       tot_p_rows,row_length,wks)
!
! Purpose: Subroutine to interpolate tropospheric photolysis rates
!          in time. Based on curve.F from Cambridge TOMCAT model.
!
!          Called from UKCA_CHEMISTRY_CTL.
!
! ---------------------------------------------------------------------
!
          USE ereport_mod, ONLY : ereport
          USE UM_ParVars
          USE Control_Max_Sizes
          IMPLICIT NONE

        INTEGER, INTENT(IN) :: p_field                   ! no of points
        INTEGER, INTENT(IN) :: p_levels                  ! no of vert
        INTEGER, INTENT(IN) :: tot_p_rows                ! no of rows
        INTEGER, INTENT(IN) :: row_length                ! no of cols

        REAL, DIMENSION(:), INTENT(IN) :: dayl           ! day length
        REAL, DIMENSION(:), INTENT(IN) :: tloc           ! local time
        REAL, DIMENSION(:,:,:), INTENT(IN) :: pjinda     ! 2D photolys

        REAL, DIMENSION(:,:), INTENT(OUT) :: wks         ! interpolated

! Local variables

        INTEGER :: i                                       ! loop variab
        INTEGER :: j                                       ! loop variables
        INTEGER :: k                                       ! loop variables
        INTEGER :: jr                                      ! loop variables

        REAL, PARAMETER :: tfrac1 = 0.04691008             ! determines
        REAL, PARAMETER :: tfrac2 = 0.23076534             ! 2D photolys

        REAL :: dawn                ! time of dawn
        REAL :: dusk                ! time of dusk
        REAL :: timel               ! local time
        REAL :: slope               ! slope used in linear interpolation
        REAL :: const               ! intercept used in linear interpola
        REAL :: fgmt(7)             ! times at which photol rates are va

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


! Initialise wks

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_CURVE',zhook_in,zhook_handle)
        wks = 0.0

! Calculate rates using a simple linear interpolation.

        DO j = 1,tot_p_rows
          DO i = 1,row_length
            k = i+(j-1)*row_length

! Non-Polar night

            IF (dayl(k) > 0.0) THEN
              dawn = 12.00 - (dayl(k)/2.0)
              dusk = 12.00 + (dayl(k)/2.0)

              fgmt(1) = dawn
              fgmt(2) = dawn + tfrac1*dayl(k)
              fgmt(3) = dawn + tfrac2*dayl(k)
              fgmt(4) = 12.00
              fgmt(5) = dawn + (1.0-tfrac2)*dayl(k)
              fgmt(6) = dawn + (1.0-tfrac1)*dayl(k)
              fgmt(7) = dusk

              timel = tloc(k)

              IF (timel > 24.0) timel = timel-24.0

              timel = min(timel,24.0)
              timel = max(timel, 0.0)

! Local Night-time

              IF (timel < dawn .OR. timel > dusk) THEN

                wks(k,1:MYjppj) = 0.0

! For the time between dawn and PJIN(1) or PJIN(5) and dusk

              ELSE IF ((timel >= dawn   .AND. timel < fgmt(2))         &
                   .OR.(timel > fgmt(6) .AND. timel <= dusk)) THEN

                IF (timel > fgmt(6)) timel = 24.00 - timel

                ! trap for -ve (timel-fgmt(1))
                IF ((fgmt(1) - timel) < 1.0E-6) timel = fgmt(1)
  
                DO jr = 1, myjppj
                  slope = pjinda(j,1,jr)/(fgmt(2) - fgmt(1))
                  wks(k,jr) = slope*(timel - fgmt(1))

                 IF (wks(k,jr) < 0.0) then
                   write(6,*) 'negative wks in ukca_curve 1',wks(k,jr)
                   write(6,*) i,j,k,jr,slope,fgmt(1),fgmt(2),timel
                   write(6,*) pjinda(j,1,jr),fgmt(1)-timel
                   write(6,*) fgmt,dawn,dusk
                   CALL EREPORT('UKCA_CURVE',jr,' Negative photolysis')
                 ENDIF
                ENDDO

! For the time between PJIN(1) and PJIN(2) or PJIN(4) and PJIN(5)

              ELSE IF ((timel >= fgmt(2) .AND. timel < fgmt(3))        &
                   .OR.(timel >  fgmt(5) .AND. timel <= fgmt(6))) THEN

                IF (timel > fgmt(5)) timel = 24.00 - timel

                DO jr = 1, myjppj
                  slope = (pjinda(j,2,jr)-pjinda(j,1,jr))/             &
                       (fgmt(3) - fgmt(2))
                  const = pjinda(j,1,jr)- slope* fgmt(2)
                  wks(k,jr)= slope*timel + const
                ENDDO

! For the time between PJIN(2), PJIN(3) and PJIN(4)

              ELSE IF (timel >= fgmt(3) .AND. timel <= fgmt(5)) THEN

                IF (timel > fgmt(4)) timel = 24.00 - timel

                DO jr = 1, myjppj
                  slope = (pjinda(j,3,jr)-pjinda(j,2,jr))/             &
                      (fgmt(4) - fgmt(3))
                  const = pjinda(j,2,jr)- slope* fgmt(3)
                  wks(k,jr)= slope*timel + const
                ENDDO

              ENDIF    ! end of IF (timel < dawn .OR. timel > dusk)

! End of the condition on the non polar night

            ELSE

              wks(k,1:myjppj) = 0.0

            ENDIF     ! end of IF (dayl(k) > 0.0)

          ENDDO       ! end of looping over row length
        ENDDO         ! end of looping over rows

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_CURVE',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_CURVE

! ----------------------------------------------------------------------

        SUBROUTINE READ2D_OPT(ipos,fpos,pjin2d)

! ----------------------------------------------------------------------
! Purpose: Subroutine to read 2D photolysis rates.
!          Based on read2d.F from Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ----------------------------------------------------------------------

          USE ukca_option_mod, ONLY: jppj, phot2d_dir
          USE ereport_mod, ONLY : ereport
          USE UM_ParVars
          USE Control_Max_Sizes
          USE chsunits_mod, ONLY : nunits
          IMPLICIT NONE

!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200


      INTEGER, INTENT(INOUT) :: ipos     ! integer position
      REAL,    INTENT(IN)    :: fpos     ! real position
      REAL, DIMENSION(:,:,:,:), INTENT(OUT) :: pjin2d   ! Photol rate

! Local variables
      INTEGER :: ifinx(myjppj)
      INTEGER :: j                             ! Loop variable
      INTEGER :: jr                            ! Loop variable
      INTEGER :: k                             ! Loop variable
      INTEGER :: kk                            ! Loop variable
      INTEGER :: in                            ! Index
      INTEGER :: ierror                        ! Error flag
      INTEGER :: info                          ! Tag for communication
      INTEGER :: errcode                       ! Variable passed to ereport

      REAL :: delpos
! Photol rates at all times
      REAL, allocatable, SAVE :: pjin2da(:,:,:,:,:)
      REAL :: pr(myjppj,3)
      REAL :: pfrac(1,3)

      LOGICAL :: L_exist                     ! Logical to check exist
      LOGICAL, SAVE :: L_first=.true.        ! Logical for firstcall 

      CHARACTER(len=256):: file2             ! Dir for photol files
      CHARACTER(len=256):: fnme(myjppj)      ! inc length of string 
      CHARACTER(len=10) :: csp(myjppj,jpspj)
      CHARACTER(len=10) :: cmnt(myjppj)      ! Comment line
      CHARACTER(len=72) :: cmessage          ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! 1.  Determine filenames containing specified photolysis rates


      IF (lhook) CALL dr_hook('UKCA_PHOT2D:READ2D_OPT',zhook_in,zhook_handle)
      pjin2d = 0.0
! Reset ipos if nearest day=365 or zero
      delpos = fpos - ipos*1.0
      IF (ipos == 74) ipos=1

! 1.1  Read filenames from photolysis ratefile on first call

      IF (L_first) THEN
        ALLOCATE(pjin2da(myjppj,n_data_times,nlphot,ntphot,nolat))

!       use module to get cmnt

        IF (SIZE(ratj_defs) /= jppj) THEN
          cmessage='size of ratj_defs is not equal to jppj'
          errcode=1
          CALL EREPORT('UKCA_PHOT2D.UKCA_READ2D',errcode,cmessage)
        ENDIF
        DO k=1,jppj
          cmnt(k)=TRIM(ADJUSTL(ratj_defs(k)%fname))
        ENDDO

!       1.2  Add '.bin' extension

        file2 = TRIM(phot2d_dir)//'/'
        DO jr = 1, MYjppj
          in = INDEX(cmnt(jr),' ') - 1
          IF ( in < 0 ) in = 10
          fnme(jr) = TRIM(ADJUSTL(file2))//                             &
               TRIM(ADJUSTL(cmnt(jr)(1:in)))//'.bin'
          IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag)             &
            WRITE(6,'(3A)') 'fnme =', fnme(jr), cmnt(jr)
        ENDDO

! 2.  Read specified photolysis rates

        IF (mype == 0) THEN
          DO jr = 1, MYjppj
            INQUIRE (FILE=FNME(JR), EXIST=L_exist)
            IF (.NOT. L_exist) THEN
              cmessage = 'File does:'//FNME(JR)//' not exist'
              CALL EREPORT('UKCA_PHOT2D.READ2D_OPT',jr,cmessage)
            ENDIF

          OPEN(ukca2pho_unit, FILE=fnme(jr), FORM='UNFORMATTED', IOSTAT=ierror)
            IF (ierror /= 0) THEN
              cmessage = ' Error opening file:'//FNME(JR)
              CALL EREPORT('UKCA_PHOT2D.READ2D_OPT',jr,cmessage)
            ENDIF

            READ(ukca2pho_unit) pjin2da(jr,:,:,:,:)
            CLOSE(ukca2pho_unit)

          ENDDO   ! jr
        ENDIF  ! mype
        l_first=.false.
      ENDIF  ! l_first

! Interpolate in time using the saved values
        IF (mype == 0) THEN
          DO jr = 1, MYjppj
            DO k = 1,nolat
              DO kk = 1,ntphot
                DO j = 1,nlphot
                  pjin2d(k,j,kk,jr) =                                  &
     &                  (pjin2da(jr,ipos+1,j,kk,k)-                    &
     &                   pjin2da(jr,ipos,j,kk,k))*delpos +             &
     &                   pjin2da(jr,ipos,j,kk,k)
              ENDDO
            ENDDO
          ENDDO
        ENDDO  ! jr
      ENDIF   ! mype

      CALL GC_RBCAST (1,nlphot*ntphot*nolat*MYjppj,0,nproc             &
                       ,info,pjin2d)

      IF (lhook) CALL dr_hook('UKCA_PHOT2D:READ2D_OPT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE READ2D_OPT

! ---------------------------------------------------------------------

      SUBROUTINE READ2D_ORIG(ipos,fpos,pjin2d)

! Purpose: Subroutine to read 2D photolysis rates.
!          Based on read2d.F from Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------

        USE ukca_option_mod, ONLY: jppj, phot2d_dir
        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        USE chsunits_mod, ONLY : nunits
        IMPLICIT NONE

!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200


      REAL, INTENT(IN) :: fpos

      INTEGER, INTENT(INOUT) :: ipos

      REAL, INTENT(OUT) :: pjin2d(nolat,nlphot,ntphot,myjppj) ! Photol rate

! Local variables

      INTEGER :: ifinx(myjppj)
      INTEGER :: ij                            ! Loop variable
      INTEGER :: j                             ! Loop variable
      INTEGER :: jr                            ! Loop variable
      INTEGER :: k                             ! Loop variable
      INTEGER :: kk                            ! Loop variable
      INTEGER :: in                            ! Index
      INTEGER :: ierror                        ! Error flag
      INTEGER :: info                          ! Tag for communication
      INTEGER :: errcode                       ! Variable passed to ereport

      REAL :: delpos
      REAL :: pjin2d1(nlphot,ntphot,nolat)! Photol rates straddling time
      REAL :: pjin2d2(nlphot,ntphot,nolat)! Photol rates straddling time
      REAL :: pr(myjppj,3)
      REAL :: pfrac(1,3)

      LOGICAL :: L_exist                       ! Logical used to check ex

      CHARACTER(len=80) :: file2             ! Dir for photol files
      CHARACTER(len=60) :: fnme(myjppj)      ! inc length of string 
      CHARACTER(len=10) :: csp(myjppj,jpspj)
      CHARACTER(len=10) :: cmnt(myjppj)      ! Comment line
      CHARACTER(len=72) :: cmessage          ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! 1.  Determine filenames from module info

      IF (lhook) CALL dr_hook('UKCA_PHOT2D:READ2D_ORIG',zhook_in,zhook_handle)
      IF (SIZE(ratj_defs) /= jppj) THEN
        cmessage='size of ratj_defs is not equal to jppj'
        errcode=1
        CALL EREPORT('UKCA_PHOT2D.UKCA_READ2D',errcode,cmessage)
      ENDIF 

      DO k=1,jppj
        cmnt(k)=TRIM(ADJUSTL(ratj_defs(k)%fname))
      ENDDO 


! 1.2  Add '.d' extension

      file2 = TRIM(phot2d_dir)//'/'

      DO jr = 1, myjppj
        in = index(cmnt(jr),' ') - 1
        IF ( in < 0 ) in = 10
        fnme(jr) = TRIM(ADJUSTL(file2))//TRIM(ADJUSTL(cmnt(jr)(1:in)))  &
             //'.dat'
        IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag)               &
          WRITE(6,'(3A)') 'fnme =', fnme(jr), cmnt(jr)
      ENDDO

! Reset ipos if nearest day=365 or zero

      delpos = fpos - ipos*1.0
      IF (ipos == 74) ipos=1

! 2.  Read specified photolysis rates

      pjin2d = 0.0
      IF (mype == 0) THEN
        DO jr = 1, myjppj
          INQUIRE (FILE=FNME(JR), EXIST=L_exist)
          IF (.NOT. L_exist) THEN
            cmessage = 'File does:'//FNME(JR)//' not exist'
            CALL EREPORT('UKCA_PHOT2D.READ2D_ORIG',jr,cmessage)
          ENDIF

          OPEN(ukca2pho_unit, FILE=fnme(jr), FORM='FORMATTED', IOSTAT=ierror)
          IF (ierror /= 0) THEN
            cmessage = 'Error reading file:'//FNME(JR)
            CALL EREPORT('UKCA_PHOT2D.READ2D_ORIG',jr,cmessage)
          ENDIF

          DO ij = 1,ipos
            read(ukca2pho_unit,*) (((pjin2d1(j,kk,k),j=1,51),kk=1,3)          &
                           ,k=nolat,1,-1)
          ENDDO

          READ(ukca2pho_unit,*) (((pjin2d2(j,kk,k),j=1,51),kk=1,3)            &
                        ,k=nolat,1,-1)

          DO k = 1,nolat
            DO kk = 1,ntphot
              DO j = 1,nlphot
                pjin2d(k,j,kk,jr) =                                    &
                         (pjin2d2(j,kk,k)-pjin2d1(j,kk,k))*            &
                          delpos + pjin2d1(j,kk,k)              
              ENDDO
            ENDDO
          ENDDO

          CLOSE(ukca2pho_unit)

        ENDDO
      ENDIF      ! End of IF mype statement

      CALL GC_RBCAST (1,nlphot*ntphot*nolat*myjppj,0,nproc             &
                       ,info,pjin2d)

      IF (lhook) CALL dr_hook('UKCA_PHOT2D:READ2D_ORIG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE READ2D_ORIG

!------------------------------------------------------------

      SUBROUTINE UKCA_INPR2D(pr2d,pr2dj)
!
! Purpose: Subroutine to calculate the pressure levels for the
!          2D photolysis rates. Original version taken from the
!          Cambridge TOMCAT model.
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
!
        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE

!       Local variables

        INTEGER, PARAMETER :: maxlev = 30

        INTEGER :: j                        ! Loop variable
        INTEGER :: ij                       ! Loop variable

        REAL, PARAMETER :: fac  = 1.284025417  ! Used to ensure that pre
        REAL, PARAMETER :: psur = 1.0e5        ! Surface pressure in Pas
        REAL, PARAMETER :: ares = 3.0          ! Factor related to verti
        REAL, PARAMETER :: eps  = 1.0e-10      ! Factor used to calculat

        REAL :: ee                       ! Temporary store
        REAL :: fj                       ! Factor related to vertical re
        REAL :: za                       ! Factor related to vertical re
        REAL :: pr2d(nolev)              ! Pressures of 2D model levels
        REAL :: pr2dj(nlphot)            ! Pressures of 2D photolysis le
        REAL :: pp(nlphot+1)
        REAL :: pres(maxlev)

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_INPR2D',zhook_in,zhook_handle)
        DO j = 1, nolev
          ee = exp((j-1)/2.0)
          pr2d(j) = psur/(ee * fac)
        ENDDO

!       2D pressure levels - normal 2-D order - up to level 30
!       for photolysis

        DO j = 1,maxlev-1
          ee = exp((j-1) / 2.0)
          pres(j) = psur / (ee * fac)
        ENDDO
        pres(maxlev) = pres(maxlev-1)

        fj    = 2.0/ares
        pp(1) = (1.0-fj)*alog(psur)+fj*alog(pres(1))
        pp(1) = exp(pp(1))

        DO ij = 2,nlphot+1
          za     = ij/(2.0*ares)
          fj     = 2.0*za+0.5+eps
          j      = int(fj)
          fj     = fj-j-eps
          j      = j+1
          pp(ij) = (1.0-fj)*alog(pres(j-1))+ fj*alog(pres(j))
          pp(ij) = exp(pp(ij))
        ENDDO

        pr2dj(1)=(psur+pp(1))*0.5

        DO ij = 2,nlphot
          pr2dj(ij) = (pp(ij)+pp(ij-1))*0.5
        ENDDO

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_INPR2D',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_INPR2D

!----------------------------------------------------------------

        SUBROUTINE UKCA_INTERPJ(pjin2d,pr2dj,zm_pl,lon,lat,lev,         & 
                                p_field,sinlat,degrad)
!
! Purpose: Subroutine to interpolate photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
          USE ereport_mod, ONLY : ereport
          USE UM_ParVars
          USE Control_Max_Sizes
          IMPLICIT NONE
        
        INTEGER :: lon                               ! No of longitudes
        INTEGER :: lat                               ! No of latitudes
        INTEGER :: lev                               ! No of levels
        INTEGER :: p_field                           ! No of spatial poi

        REAL, DIMENSION(:,:,:,:), INTENT(IN) :: pjin2d  ! 2D photo
        REAL, DIMENSION(:,:), INTENT(IN)     :: zm_pl   ! zmean 3D press 
        REAL, DIMENSION(:), INTENT(IN)       :: pr2dj   ! 2D photo
        REAL, DIMENSION(:), INTENT(IN)       :: sinlat  ! Sine (3D
        REAL, INTENT(IN)                     :: degrad  ! To conve
        
!       Local variables

        INTEGER :: i                     ! Loop variable
        INTEGER :: j                     ! Loop variable
        INTEGER :: jr                    ! Loop variable
        INTEGER :: k                     ! Loop variable
        INTEGER :: kk                    ! Loop variable
        INTEGER :: l                     ! Loop variable

        REAL, PARAMETER :: Npole = 90.0

        REAL :: p2d(nlphot)
        REAL :: lat2d(nolat)          ! 2D model latitudes
        REAL :: wks1(lat,nlphot)      ! Working array
        REAL :: wks2(nolat)           ! Working array
        REAL :: wks3(nlphot)          ! Working array
        REAL :: wks4(lat,lev)         ! Working array
        REAL :: lati
        REAL :: press
        REAL :: UKCA_FLUPJ            ! Function used for interpolation
        REAL :: delphi

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       Set up 2D latitudes. lat=1 pole nord
!       LAT2D() is the latitude in the centre of the 2D box in radians

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_INTERPJ',zhook_in,zhook_handle)
        delphi = Npole/nolat
        DO i = 1,nolat
          lat2d(i)=sin((90.0 - (2*i-1)*delphi)*degrad)
        ENDDO

        DO jr = 1,MYjppj              ! Loop over photolysis reactions

!         Interpolate linearly in Sin(lat) (KK is the point through the

          DO kk = 1,ntphot          ! Loop over times of day

            DO j = 1,nlphot
              DO i = 1,nolat
                wks2(i) = pjin2d(i,j,kk,jr)
              ENDDO

              DO k = 1,lat
                lati = sinlat((k-1)*lon+1)
                lati = MAX(lati, lat2d(nolat))
                lati = MIN(lati, lat2d(    1))
!DEPENDS ON: ukca_flupj
                wks1(k,j) = UKCA_FLUPJ(lati,lat2d,wks2,nolat)
                wks1(k,j) = MAX(wks1(k,j), 0.0)
              ENDDO
            ENDDO

!           Interpolate linearly in log(P)

            DO k = 1,lat
              DO j = 1,nlphot
                wks3(j) = wks1(k,j)
                p2d(j)  = LOG(pr2dj(j))
              ENDDO

              DO l = 1,lev
                press = LOG(zm_pl(k,l)) 
                press = MAX(press, p2d(nlphot))
                press = MIN(press, p2d(1))
!DEPENDS ON: ukca_flupj
                wks4(k,l) = UKCA_FLUPJ(press,p2d,wks3,nlphot)
                wks4(k,l) = MAX(wks4(k,l), 0.0)
              ENDDO
            ENDDO

            DO l = 1,lev
              DO k = 1,lat
                pjin(k,l,kk,jr) = wks4(k,l)
              ENDDO
            ENDDO


          ENDDO     ! End of loop over times of day
        ENDDO       ! End of loop over photolysis reactions

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_INTERPJ',zhook_out,zhook_handle)
        RETURN

      END SUBROUTINE UKCA_INTERPJ

!-----------------------------------------------------------------

      SUBROUTINE PHOT2D_ALLOCATE_MEMORY(lat,lev,ntphot)
!
! Purpose: Subroutine to allocate array to hold
!          2D photolysis rates
!
!          Called from UKCA_PHOTIN.
!
! ---------------------------------------------------------------------
        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE
        INTEGER,INTENT(IN)         :: lat,lev,ntphot

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

   
!       lat is the number of latitude on current PE, NOT nolat!
!       lev is numer of model level, NOT nolev

        IF (lhook) CALL dr_hook('UKCA_PHOT2D:PHOT2D_ALLOCATE_MEMORY',zhook_in,zhook_handle)
        ALLOCATE(pjin(lat,lev,ntphot,MYjppj))

      IF (lhook) CALL dr_hook('UKCA_PHOT2D:PHOT2D_ALLOCATE_MEMORY',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PHOT2D_ALLOCATE_MEMORY

!---------------------------------------------------------------------- 
         
      SUBROUTINE UKCA_CALC_ZMEAN_PRESS(lon, lat, lev,              & 
    &                                  glon,pl, zmean_pl) 
! 
! Purpose: Subroutine to calculate zonal mean pressure profile 
! 
!          Called from UKCA_PHOTIN. 
! 
! --------------------------------------------------------------------- 
      USE global_2d_sums_mod, ONLY: global_2d_sums
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE 


      INTEGER, INTENT(IN) :: lon        ! No of longitudes 
      INTEGER, INTENT(IN) :: lat        ! No of latitudes 
      INTEGER, INTENT(IN) :: lev        ! No of levels 
      INTEGER, INTENT(IN) :: glon       ! No of global lons 

      REAL, INTENT(IN)    :: pl(lon,lat,lev)   ! Pressure 

      REAL, INTENT(OUT)   :: zmean_pl (lat,lev) 

!     Local variables 

      INTEGER :: j,k                 ! Loop variables 
      INTEGER :: istat               ! Output status from GCOM routine 

      REAL    :: fac                 ! Factor used to calculate zmean 
      REAL    :: sumpl_row(lat,lev)  ! Global sum along rows 

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

         
      IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_CALC_ZMEAN_PRESS',zhook_in,zhook_handle)
      zmean_pl = 0.0 

      CALL global_2d_sums(pl, lon, 1, 0, 0, lev*lat,                    &
                          sumpl_row, gc_proc_row_group)

      fac      = 1.0/real(glon) 
      zmean_pl = sumpl_row*fac 
         
      IF (lhook) CALL dr_hook('UKCA_PHOT2D:UKCA_CALC_ZMEAN_PRESS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALC_ZMEAN_PRESS 

!---------------------------------------------------------------------- 

      END MODULE UKCA_PHOT2D
