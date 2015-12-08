! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_hyd
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: land surface

      Subroutine diagnostics_hyd(                                       &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length                       &
     &,                      at_extremity                               &
     &,                      land_points, dsm_levels                    &
!  Add inland basin outflow to arguments
     &,                      land_index,inlandout_atm                   &
     &,                      smc, surf_roff, sub_surf_roff              &
     &,                      snow_depth_land, snow_melt, canopy, t_soil &
     &,                      soil_layer_moisture                        &
     &,                      ntiles, snomlt_surf_htf, sthu, sthf        &
     &,                      tot_tfall, snow_tile, melt_tile            &
     &,                      rgrain, land_sea_mask                      &
     &,                      dun_roff, drain, qbase, qbase_zw           &
     &,                      fch4_wetl, fexp, gamtot, ti_mean, ti_sig   &
     &,                      fsat, fwetl, zw, sthzw                     &
     &,                      timestep                                   &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork                                                        &
     &     )

! Description:
!   Calculates hydrology-related diagnostics (held in STASH section 8).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   hydrology routines called previously. Each diagnostic is simply
!   copied into the STASHwork array to be passed on to STASH for
!   output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    208  Soil moisture content
!     23  Snow depth
!    201  Snow melt
!    234  Surface run-off
!    235  Sub-surface run-off
!    225  Soil temperature
!    223  Soil layer moisture
!    204  Surface run-off (accumulated)
!    205  Sub-surface run-off (accumulated)
!    209  Canopy water content
!    202  Land snow melt heat flux
!    231  Land snow melt rate
!    233  Canopy throughfall rate
!    236  Snow amount on tiles
!    237  Snow melt rate on tiles
!    238  Snow grain size on tiles
!    229  Unfrozen soil moisture fraction
!    230  Frozen soil moisture fraction
!    239  Baseflow
!    240  Dunne Runoff
!    241  Baseflow from water table layer
!    242  Wetland methane flux
!    252  Drainage out of "nshyd"th model layer
!    245  Inland basin outflow on atmos grid
!    252  Drainage out of bottom "nshyd"th soil layer (currently layer 4).

      USE mask_compression, ONLY:  expand_from_mask
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE Submodel_Mod
      IMPLICIT NONE

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain     ! indicator as to model type, ie global, lam

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, land_points                                                     &
                    ! No.of land points being processed, can be 0.
     &, dsm_levels                                                      &
     &, ntiles      ! No. of land-surface tiles ( MOSES II )

      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP


      Real                                                              &
     &  timestep

! Primary Arrays used in all models
      Integer                                                           &
     &  land_index(land_points)      ! set from land_sea_mask

      Real                                                              &
     &  snow_depth_land(land_points)                                    &
                                     !
     &, snow_melt(land_points)                                          &
                                ! snowmelt (kg/m2/s)
     &, canopy (land_points)                                            &
                              ! canopy water content (kg/m2)
     &, smc(land_points)                                                &
                            ! available soil moisture in the
!                                 rootzone (kg/m2).
     &, surf_roff(land_points)                                          &
                                ! surface runoff (kg/m2/s).
     &, sub_surf_roff(land_points)                                      &
                                    ! sub-surface runoff
! Declare inland basin outflow variable
     &, inlandout_atm(land_points)                                      &
                                    !inland basin outflow

     &, t_soil(land_points,dsm_levels)                                  &
     &, soil_layer_moisture(land_points,dsm_levels)                     &
     &, snomlt_surf_htf(row_length, rows)                               &
     &, sthu(land_points,dsm_levels)                                    &
                                       ! Unfrozen soil moisture
!                content of each layer as a fraction of saturation
     &, sthf(land_points,dsm_levels)                                    &
                                       ! Frozen soil moisture content
!                of each layer as a fraction of saturation
     &, tot_tfall(land_points)                                          &
                                       ! total throughfall (kg/m2/s)
     &, snow_tile(land_points,ntiles)                                   &
                                       ! Lying snow on tiles (kg/m2)
     &, melt_tile(land_points,ntiles)                                   &
                                       ! Snowmelt on tiles (kg/m2/s)
     &, rgrain(land_points,ntiles)                                      &
                                       ! Snow grain size (microns)

! Additional variables required for large-scale hydrology:
     &, qbase(land_points)                                              &
                                    ! Baseflow (kg/m2/s).
     &, dun_roff(land_points)                                           &
                                    ! Dunne runoff (kg/m2/s).
     &, qbase_zw(land_points)                                           &
                                    ! Baseflow from Zw layer (kg/m2/s).
     &, drain(land_points)                                              &
                                    ! Drainage out of "nshyd" later
!                                   ! (kg/m2/s).
     &, fch4_wetl(land_points)                                          &
                                    ! Wetland methane flux (kg C/m2/s).
     &, TI_MEAN(LAND_POINTS)                                            &
                                    ! Mean topographic index.
     &, TI_SIG(LAND_POINTS)                                             &
                                    ! Std. dev. of topographic index.
     &, FEXP(LAND_POINTS)                                               &
                                    ! Decay factor in Sat. Conductivity
!                                   !   in water table layer.
     &, GAMTOT(LAND_POINTS)                                             &
                                    ! Integrated complete Gamma
!                                   !   function.
     &, FSAT(LAND_POINTS)                                               &
                                    ! Surface saturation fraction.
     &, FWETL(LAND_POINTS)                                              &
                                    ! Wetland fraction.
     &, ZW(LAND_POINTS)                                                 &
                                    ! Water table depth (m).
     &, STHZW(LAND_POINTS)          ! Saturation fraction in water
!                                   !   table layer.


      Logical                                                           &
     &  land_sea_mask(row_length, rows)

      Integer                                                           &
     & PSLEVEL                                                          &
                     !  loop counter for pseudolevels
     &,PSLEVEL_OUT                                                      &
                     !  index for pseudolevels sent to STASH
     &,LEVEL                                                            &
                     !  loop counter for levels
     &,LEVEL_OUT     !  index for levels sent to STASH

      Logical                                                           &
     & PLLTILE(NTILES)                                                  &
                          ! pseudolevel list for surface types
     &,LIST(DSM_LEVELS)   ! level list for soil levels


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

! Local variables

      Integer                                                           &
     &  i, j, k, l                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      CHARACTER(LEN=80)  cmessage
      CHARACTER(LEN=*) RoutineName
      Parameter ( RoutineName='diagnostics_hyd')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     &  interp_data(row_length,rows)                                    &
     &, interp_data_3(row_length,rows,dsm_levels)

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

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1.  Initialisation.
! ----------------------------------------------------------------------


      IF (lhook) CALL dr_hook('DIAGNOSTICS_HYD',zhook_in,zhook_handle)
      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Soil moisture content
! ----------------------------------------------------------------------
! Item 8 208  smc

      If (sf(208,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = smc(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(208,8,im_index)),interp_data,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,208,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
! Snow depth
! ----------------------------------------------------------------------
! Item 8 023 snow_depth_land

      If (sf(023,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_depth_land(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(023,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,023,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 023)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Snow melt
! ----------------------------------------------------------------------
! Item 201

      If (sf(201,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_melt(l) * timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(201,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
! Surface run-off.
! ----------------------------------------------------------------------
! Item 234  surf_roff

      If (sf(234,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = surf_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(234,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,234,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 234)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Sub-surface run-off.
! ----------------------------------------------------------------------
! Item 235  sub_surf_roff

      If (sf(235,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sub_surf_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(235,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,235,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 235)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
!  Soil temperature
! ----------------------------------------------------------------------
! Item 8 225 t_soil

      If (sf(225,8)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               End Do
            End Do

            Do l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) = t_soil(l,k)
            End Do

         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,8,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,225,8,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,8,225,                                           &
     &        icode,cmessage)



         If (icode  >   0) then
            cmessage="Error in copydiag_3d( item 225)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Soil layer moisture
! ----------------------------------------------------------------------
! Item 8 223 soil_layer_moisture

      If (sf(223,8)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               End Do
            End Do
             Do l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) =                                   &
     &              soil_layer_moisture(l,k)
             End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,8,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,223,8,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,8,223,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 204  surf_roff

      If (sf(204,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = surf_roff(l) * timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204 )"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Sub-surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 205  sub_surf_roff

      If (sf(205,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sub_surf_roff(l)             &
     &     *  timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(205,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,205,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         endif
      endif

! ----------------------------------------------------------------------
! Canopy water content
! ----------------------------------------------------------------------
! Item 209

      If (sf(209,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = canopy(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(209,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,209,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 209)"
            goto 9999
         endif
      endif


! ----------------------------------------------------------------------
! Land snow melt heat flux (W/m2)
! ----------------------------------------------------------------------
! Item 202 SNOMLT_SURF_HTF

      IF (SF(202,8)) THEN
        CALL expand_from_mask (                                           &
            STASHWORK(SI(202,8,im_index)),                             &
            snomlt_surf_htf,                                           &
            land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Land snow melt heat rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 231 snow_melt

      IF (SF(231,8)) THEN
        CALL expand_from_mask (                                           &
            STASHWORK(SI(231,8,im_index)),                             &
            snow_melt,                                                 &
            land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Canopy throughfall rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 233 tot_tfall

      IF (SF(233,8)) THEN
        CALL expand_from_mask (                                           &
            STASHWORK(SI(233,8,im_index)),                             &
            tot_tfall,                                                 &
            land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Snow amount on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 236 snow_tile

      IF (SF(236,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,236,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 236 = snow_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                       &
                STASHWORK(SI(236,8,im_index)+(PSLEVEL_OUT-1)           &
                *row_length*rows),snow_tile(1,PSLEVEL_OUT),            &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow melt rate on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 237 melt_tile

      IF (SF(237,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,237,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 237 = melt_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                       &
                STASHWORK(SI(237,8,im_index)+(PSLEVEL_OUT-1)           &
                *row_length*rows),melt_tile(1,PSLEVEL_OUT),            &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow grain size on tiles
! ----------------------------------------------------------------------
! Item 238 rgrain

      IF (SF(238,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,238,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 238 = rgrain)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
            CALL expand_from_mask (                                       &
                STASHWORK(SI(238,8,im_index)+(PSLEVEL_OUT-1)           &
                *row_length*rows),rgrain(1,PSLEVEL_OUT),               &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Unfrozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 229 sthu

      IF (SF(229,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,229,8,im_index)),                       &
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_levels_list(item 229 = sthu)"
            goto 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
            CALL expand_from_mask (                                       &
                STASHWORK(SI(229,8,im_index)+(LEVEL_OUT-1)             &
                *row_length*rows),sthu(1,LEVEL_OUT),                   &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Frozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 230 sthf

      IF (SF(230,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,230,8,im_index)),                       &
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_levels_list(item 230 = sthf)"
            goto 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
            CALL expand_from_mask (                                       &
                STASHWORK(SI(230,8,im_index)+(LEVEL_OUT-1)             &
                *row_length*rows),sthf(1,LEVEL_OUT),                   &
                land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Baseflow
! ----------------------------------------------------------------------
! Item 239  baseflow

      If (sf(239,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(239,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,239,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 239)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Dunne Runoff
! ----------------------------------------------------------------------
! Item 240 Dunne Runoff

      If (sf(240,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = dun_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(240,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,240,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 240)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Baseflow from water table layer
! ----------------------------------------------------------------------
! Item 241  baseflow from water table layer

      If (sf(241,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase_zw(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(241,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,241,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Wetland methane flux
! ----------------------------------------------------------------------
! Item 242  Wetland methane flux

      If (sf(242,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(242,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,242,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if
      If (sf(243,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_mean(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(243,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,243,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if
      End if
      If (sf(244,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_sig(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(244,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,244,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if
      End if
      If (sf(251,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fexp(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(251,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,251,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 251)"
            goto 9999
         End if
      End if
      If (sf(246,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = gamtot(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(246,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,246,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if
      End if
      If (sf(247,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fsat(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(247,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,247,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if
      End if
      If (sf(248,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fwetl(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(248,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,248,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if
      End if
      If (sf(249,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = zw(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(249,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,249,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 249)"
            goto 9999
         End if
      End if
      If (sf(250,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sthzw(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(250,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,250,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 250)"
            goto 9999
         End if
      End if
! ----------------------------------------------------------------------
! Drainage out of  "nshyd"th model layer
! ----------------------------------------------------------------------
! Item 252  drainage

      If (sf(252,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = drain(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(252,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,252,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 252)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Output inland basin outflow on atmosphere grid

! --------------------------------------------------------------------
! Inland basin outflow (atmos grid)
! --------------------------------------------------------------------
! Item 245  Inland basin outflow

      If (sf(245,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = inlandout_atm(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(245,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,245,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if
! -------------------------------------------
 9999 continue
      If (icode /= 0) Then

        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      IF (lhook) CALL dr_hook('DIAGNOSTICS_HYD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE diagnostics_hyd
