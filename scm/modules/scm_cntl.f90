!
! MODULE scm_cntl---------------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT***************************************


MODULE scm_cntl_mod

!-------------------------------------------------------------------------------
! Description:
!   Declares variables and defines SCM namelist which is read from CNTLATM
!   
!-------------------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!-------------------------------------------------------------------------------

  USE scmoptype_defn

  IMPLICIT NONE



! Declare variables for diagnostics output streams
  CHARACTER(LEN=200) :: scm_nml

  INTEGER :: main_diag_switch
                     ! If set to zero the diagnostic system
                     ! is  off: no diagnostics will be calculated
                     ! or output. Default=1.

  INTEGER :: strm_switch(maxnstreams) = 0
                     ! If set to zero the respective
                     ! stream is off. Default=1 for stream 1,
                     ! 0 for all others.

  INTEGER :: strm_unit(maxnstreams) = 0 
                     ! The output unit to which the
                     ! stream will be written. Default=36+n.



  CHARACTER (len=lfilename) :: strm_filename(maxnstreams) = Default
                     ! The name of the file to which the stream
                     ! will be written. Default=compiler dependent.

  INTEGER :: strm_format(maxnstreams) = 3
                     ! Format of the output data:
                     ! 0 = can be read into PV-Wave with scmread2.pro
                     ! 1 = can be easily read with human eye
                     ! 2 = SSFM format, compatible with FSSI d'base
                     ! Default = 0.

  INTEGER :: strm_dumpstep(maxnstreams) = 1
                     ! Dumping period of the stream, i.e. number of timesteps
                     ! between dumps. If negative it's treated as a number
                     ! of seconds and is converted to the closest possible
                     ! non-zero number of timesteps. If the number of dumps
                     ! per day is `non-integer', a warning will be printed.
                     ! Default = 1.

  INTEGER :: strm_heed_hardwired(maxnstreams) = 1
                     ! If zero then diagnostics sent to this stream by the
                     ! hard-wired stream list provided in the call to
                     ! SCMoutput are ignored. Default = 1.

  INTEGER :: strm_heed_acceptlist(maxnstreams) = 1
                     ! If zero then the contents of
                     ! strm_acceptlist are ignored.

  INTEGER :: strm_heed_rejectlist(maxnstreams) = 1
                     ! If zero then the contents of
                     ! strm_rejectlist are ignored.

  CHARACTER (len=listlength) :: strm_acceptlist(maxnstreams) = ''
                     ! Comma-separated list of diagnostic (short) names
                     ! (as specifed in corresponding calls to SCMoutput)
                     ! which are to be sucked into the respective stream.
                     ! This provides a way of sending a diagnostic to a
                     ! stream even if that stream is not specified as a
                     ! destination in the hard-wired list of the
                     ! respective SCMoutput call. Default = ''

  CHARACTER (len=listlength) :: strm_rejectlist(maxnstreams) = ''
                     ! Like strm_acceptlist but this is a list of
                     ! diagnostics which are to be prevented from being
                     ! sent to the respective stream. If there is any
                     ! conflict between rejectlist and either acceptlist
                     ! and/or the "hard-wired" preferences, rejectlist
                     ! wins.  Default = ''

  INTEGER :: netcdf_chunksize
                     ! The ChunkSize input to NF90_Create(),
                     ! the routine that creates a NetCDF file.
                     ! Controls a space versus time trade-off:
                     ! memory allocated in the netcdf library
                     ! versus number of system calls.

                     ! Individual logicals for SCM diagnostics
                     ! packages:

  LOGICAL ::           &
    l_scmdiag_gen      &! General diagnostics  - default true
  , l_scmdiag_rad      &! Radiation            - default false
  , l_scmdiag_bl       &! Boundary layer       - default false
  , l_scmdiag_surf     &! Surface              - default false
  , l_scmdiag_land     &! Land points only     - default false
                        ! reset in Scm_Main if not land point
  , l_scmdiag_sea      &! Sea points only      - default false
                        ! reset in Scm_Main if not sea point
  , l_scmdiag_lsp      &! Large scale precip   - default false
  , l_scmdiag_conv     &! Convection           - default false
  , l_scmdiag_lscld    &! Large scale cloud    - default false
  , l_scmdiag_pc2      &! PC2                  - default false
                        ! set in Scm_Main to false if PC2 off
  , l_scmdiag_gwd      &! Gravity waves        - default false
  , l_scmdiag_forc     &! Forcing              - default false
  , l_scmdiag_incs      ! Increments           - default false


  NAMELIST/SCM_CNTL/                                                      &
    scm_nml, strm_filename, strm_switch

  !----------------------------------------------------------------------
  ! The DIAGS namelist
  ! This is read from the forcing namelist but may eventually
  ! move to being produced via the UMUI. As of UM?.?,
  ! strm_filename and strm_switch are taken from SCM_CNTL
  !----------------------------------------------------------------------
  NAMELIST/DIAGS/                                                         &
    main_diag_switch, strm_unit, strm_format, strm_dumpstep,              &
    strm_heed_hardwired, strm_heed_acceptlist, strm_heed_rejectlist,      &
    strm_acceptlist, strm_rejectlist, netcdf_chunksize,                   &
    l_scmdiag_gen,   l_scmdiag_rad, l_scmdiag_bl,   l_scmdiag_surf,       &
    l_scmdiag_land,  l_scmdiag_sea, l_scmdiag_lsp,  l_scmdiag_conv,       &
    l_scmdiag_lscld, l_scmdiag_pc2, l_scmdiag_gwd,  l_scmdiag_forc,       &
    l_scmdiag_incs

!=============================================================================
END MODULE scm_cntl_mod
