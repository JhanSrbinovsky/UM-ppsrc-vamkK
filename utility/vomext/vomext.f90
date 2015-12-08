! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! *********************************************************************
!+ Program VOM_Extract : Extracts Vertical Profiles from UM data.
!
PROGRAM VOM_Extract
  USE IO
    USE ereport_mod, ONLY : ereport, ereport_finalise
    
    USE PrintStatus_mod
    USE filenamelength_mod, ONLY :                              &
        filenamelength
    USE UM_Config, ONLY : &
        appInit, &
        exe_vomext
    IMPLICIT NONE
!
! Description :
!    Program to read in model data from UM Fieldsfiles and extrcat 
!    required data to run the 3DVOM model. 
!
! Method : 
!    1. Namelist is read in - this will contain a list of lat/long of
!       points and times for which data is to be extracted for.
!       Maximums allowed : 50 profiles, 10 times.
!    2. The lookup headers for all the input fieldsfiles are read in.
!       Maximum of 20 fieldsfiles allowed.
!    3. For each lat/long the 'nearest' row/col on the model grid is
!       determined. Method used is very basic. The algorithm could be
!       improved to determine the nearest point more accurately. It
!       has only been tested with Global Model Level fieldsfiles. For
!       rotated grids, new code would be needed to convert the namelist
!       lat/long to the rotated lat/long.
!    4. Determine pointers to the data in the fieldsfiles. currently,
!       6 pointers :
!         - Orography               Stash Code 33
!         - U/V                     Stash Code 2/3
!         - Theta                   Stash Code 4 
!         - Density                 Stash Code 253
!         - Pressure on Rho levels  Stash Code 407
!    5. Data is read in, extracted and written out. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! Language : FORTRAN 90
! This code is written to UMDP3 programming standards.
!
! Declarations :

  INTEGER :: icode            !  Error code
  INTEGER :: me_gc
  INTEGER :: nproc_gc

  CHARACTER(LEN=80) :: cmessage          ! Error Message

  CHARACTER(LEN=*), PARAMETER :: RoutineName = "VOM_Extract"

! DEPENDS ON: InitPrintStatus
  CALL InitPrintStatus
  CALL gc_init(' ',me_gc,nproc_gc)
  CALL appInit(exe_vomext)

! DEPENDS ON: Timer
  CALL Timer ( RoutineName, 1 )
  CALL ioInit()
  WRITE (6,*) ' ###################################### '
  WRITE (6,*) ' Running VOM Extract Utility to extract '
  WRITE (6,*) ' vertical profiles from UM data         '
  WRITE (6,*) ' ###################################### '
  WRITE (6,*) ' '

  icode = 0


  CALL UM_Profiles


  CALL Ereport_finalise( ) 

  CALL ioShutDown()
! DEPENDS ON: Timer
  CALL Timer ( RoutineName, 2 )

  WRITE (6,*) ' '
  WRITE (6,*) ' #######################################'
  WRITE (6,*) ' VOM Extract program completed normally.'
  WRITE (6,*) ' #######################################'

CONTAINS


!+ Subroutine UM_Profiles

! Subroutine Interface :
  SUBROUTINE UM_Profiles

    USE ereport_mod, ONLY : ereport
    USE PrintStatus_mod
    IMPLICIT NONE
!
! Description :
!     Routine for extract UM data for 3DVOM model
!
! Method :
!     See main program.
!
!
! Language : FORTRAN 90
! This code is written to UMDP3 programming standards.
!
! Declarations :
     
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

!   Local variables

    INTEGER :: icode                  !  Error code
    CHARACTER (LEN=80) :: cmessage    !  Error Message
    CHARACTER (Len=*), PARAMETER :: RoutineName = 'UM_Profiles'

    INTEGER :: len_env    = 6  !  Length of env. variable
    INTEGER :: env_var    = 0  !  Indicator that filename is in env var
    INTEGER :: read_only  = 0  !  Input Fieldsfiles   - read only
    INTEGER :: read_write = 1  !  Profile output file - read & write
    CHARACTER(LEN=6) env            !  Env Variable for input fieldsfiles

    INTEGER :: NL_unit_no   = 10  ! Unit No for Namelist File
    INTEGER :: prof_unit_no = 20  ! Unit No for Profile data
    INTEGER :: FF_unit_no         ! Unit No for FFs (31-)

    INTEGER :: Len_FixHd = 256
    !     Integer :: FixHd(Len_FixHd)  ! Compiler doesn't like this line ?
    INTEGER :: FixHd (256)

    INTEGER :: Len_IHead = 10
    INTEGER :: Len_RHead = 10
    INTEGER :: IntHead (10)
    REAL    :: RealHead(10)

    INTEGER :: ReturnCode
    INTEGER :: i,j,k,n

    INTEGER, PARAMETER :: max_n_times    = 10
    INTEGER, PARAMETER :: max_n_profiles = 50
    INTEGER, PARAMETER :: max_n_levs     = 70
    REAL :: orog  (max_n_profiles)
    REAL :: u     (max_n_levs,max_n_profiles)
    REAL :: v     (max_n_levs,max_n_profiles)
    REAL :: p     (max_n_levs,max_n_profiles)
    REAL :: rho   (max_n_levs,max_n_profiles)
    REAL :: theta (max_n_levs,max_n_profiles)

    CHARACTER (len=filenamelength) :: FileName   
    !  FileName for Namelist and Profile output file

    INTEGER :: Len1_Lookup 
    INTEGER :: Len2_Lookup
    INTEGER :: Len_Lookup
    INTEGER, ALLOCATABLE :: Lookup (:,:)
    INTEGER, ALLOCATABLE :: Lookup_FF_No (:)

    INTEGER :: n_lookups (0:20)
    INTEGER :: LookUp_Start_Address(20)
    INTEGER :: i_FF, n_FF
    INTEGER :: ihead

    REAL, ALLOCATABLE :: Data_In (:,:,:)
    REAL, ALLOCATABLE :: srce_lat (:)
    REAL, ALLOCATABLE :: srce_long(:)

    REAL :: a_io
    INTEGER :: len_io
    REAL    :: prof_lat
    REAL    :: prof_long
    INTEGER :: iprof
    INTEGER :: n_profiles
    INTEGER :: n_times
    INTEGER :: itime
    INTEGER :: irow (max_n_profiles)
    INTEGER :: icol (max_n_profiles)

    INTEGER :: ipt_orog
    INTEGER :: ipt_u
    INTEGER :: ipt_v
    INTEGER :: ipt_theta
    INTEGER :: ipt_density
    INTEGER :: ipt_rho_p

    INTEGER :: ix
    INTEGER :: iy

    !     Validity Time
    INTEGER :: vt_yy
    INTEGER :: vt_mm 
    INTEGER :: vt_dd
    INTEGER :: vt_hh

    INTEGER :: row_length
    INTEGER :: rows
    INTEGER :: nlevs

!     ---------------
!     Namelist UMPROF
!     ---------------

    REAL :: Latitude  (max_n_profiles) !  Latitude of profile
    REAL :: Longitude (max_n_profiles) !  Longitude of profile
    INTEGER :: FC_Time (max_n_times)   !  Forecast time data required
    NAMELIST /UMPROF/ Latitude, Longitude, FC_Time

!   -------------------------
!   Get filename for Namelist
!   -------------------------

    CALL Fort_Get_Env ('UNIT10', 6, FileName, filenamelength, ICode)

    WRITE (6,*) ' NameList File : ',FileName(1:LEN_TRIM(FileName))

!   ----------------------
!   Open the Namelist File
!   ----------------------

    OPEN ( Unit=NL_unit_no, File=FileName,       &
         Status='unknown', iostat= icode )

    WRITE (6,*) ' Namelist file opened : icode ',icode

    IF (icode /= 0) THEN
      CMessage = 'Error in opening namelist file'

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

!   -----------------
!   Read the Namelist
!   -----------------

    Latitude (:)  = RMDI
    Longitude (:) = RMDI
    FC_TIME (:)   = -1

    READ (NL_unit_no, NML=UMPROF)
    WRITE (6, NML=UMPROF)

!   -------------------
!   Close Namelist File
!   -------------------

    WRITE (6,*) ' closing namelist file'
    CLOSE ( NL_unit_no, iostat=icode )
    WRITE (6,*) ' icode ',icode
    IF (icode /= 0) THEN
      CMessage = 'Error in closing namelist file'

      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

!   ------------------------------------
!   Work out how many profiles requested
!   ------------------------------------

    n_profiles = 0
    DO i = 1, max_n_profiles
      IF ( Latitude(i) /= rmdi ) THEN
        n_profiles = n_profiles + 1
      END IF
    END DO
    WRITE (6,*) ' n_profiles ',n_profiles

!   ------------------------------------------
!   Work out how many forecast times requested
!   ------------------------------------------

    n_times = 0
    DO i = 1,max_n_times
      IF ( FC_time(i) /= -1 ) THEN
        n_times = n_times + 1
      END IF
    END DO
    WRITE (6,*) ' n_times ',n_times

!   ------------------------------------
!   Get filename for profile output file
!   ------------------------------------

    CALL Fort_Get_Env ('UNIT20', 6, FileName, filenamelength, ICode)

    WRITE (6,*) ' Output Profile File : ',           &
         FileName(1:LEN_TRIM(FileName))

!   ----------------------------
!   Open the profile output file
!   ----------------------------

!   Cannot use file_open, must use fortran open

    OPEN ( Unit=prof_unit_no, File=FileName,         &
         Status='unknown', iostat= icode )

    WRITE (6,*) ' output file opened : icode ',icode

    IF (icode /= 0) THEN
      CMessage = 'Error in opening profile output file'

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

!   ----------------------------------
!   Work out how many FFs to be opened
!   ----------------------------------

    n_FF = 0
    DO i_FF = 1, 20

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no

      FileName = 'NotSet'
      CALL Fort_Get_Env (env, 6, FileName, filenamelength, ICode)

      !       WRITE (6,*) ' ICode ',ICode
      !       WRITE (6,*) ' Filename ',Filename

      IF (Icode == -1) THEN

        ! Return Code is -1 if Env Var is not found
        WRITE (6,*) ' Env Var ',env,' not found.'
        CYCLE

      END IF

      IF (Icode == 0) THEN

        ! Env Var found
        WRITE (6,*) ' Env Var ',env,' found.'
        WRITE (6,*) ' Filename ',Filename
        n_FF = n_FF + 1

      END IF

    END DO
    WRITE (6,*) ' n_FF = ',n_FF

    IF (n_FF == 0) THEN
      ICode = 13 
      CMessage = 'No input FieldsFiles available?'

      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

!   --------------------
!   Open the fieldsfiles
!   --------------------

    n_lookups (:) = 0

    DO i_FF = 1, n_FF

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no
      WRITE (6,*) ' '

      !     Test env var exists

      CALL File_Open (ff_unit_no,env,len_env,read_only,env_var,icode)
      IF (icode /= 0) THEN
        CMessage = 'Error in opening fieldfiles'

        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     ------------------------
!     Read in the fixed header
!     ------------------------

! DEPENDS ON: Read_FLH
      CALL Read_FLH ( ff_unit_no,fixhd,len_fixhd,icode,cmessage )
      IF (icode > 0) THEN
        CMessage = 'Error in Read_FLH'

        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

      !     WRITE (6,*) ' Fixed Header for FFs '
      !     DO i=1,256 
      !     WRITE (6,*) i,fixhd(i)
      !     END DO

      IF (fixhd(5) == 3) THEN
        WRITE (6,*) 'This is a FieldsFile'
      ELSE
        CMessage= 'Invalid input file type'
        Returncode = 10

        CALL Ereport( RoutineName, ReturnCode, CMessage )
      END IF

!     -------------------------------
!     Move to start of Integer Header
!     -------------------------------

      CALL SetPos ( FF_unit_no, FixHd(100)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for Integer Header'

        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     ----------------------
!     Read in Integer Header
!     ----------------------

      CALL BuffIn ( FF_unit_no,                          &
           IntHead, Len_IHead, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_ihead) THEN      
! DEPENDS ON: IOERROR
        CALL IOERROR('buffer in of Integer header',A_IO,LEN_IO,  &
             len_lookup)
        CMESSAGE = 'I/O ERROR with buffin of Integer Header'
        ICode    = 11

        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     ---------------------------------------
!     Get grid dimensions from Integer Header
!     ---------------------------------------

      row_length = IntHead (6)
      rows       = IntHead (7)
      nlevs      = IntHead (8)
      WRITE (6,*) ' Model grid : ',row_length,' x ',rows
      WRITE (6,*) ' No of model levels in FFs ',nlevs

!     ----------------------------
!     Move to start of Real Header
!     ----------------------------

      CALL SetPos ( FF_unit_no, FixHd(105)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for Real Header'

        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     -------------------
!     Read in Real Header
!     -------------------

      CALL BuffIn ( FF_unit_no,                             &
           RealHead, Len_RHead, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_rhead) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR('buffer in of Real header',A_IO,LEN_IO,  &
             len_lookup)
        CMESSAGE = 'I/O ERROR with buffin of Real Header'
        ICode    = 11

        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     ---------------------------------------
!     Get lookup dimensions from fixed header
!     ---------------------------------------

      LookUp_Start_Address(i_FF) = FixHd(150)
      Len1_Lookup = FixHd(151)
      Len2_Lookup = FixHd(152)
      Len_Lookup  = Len1_Lookup * Len2_Lookup
      WRITE (6,*) ' len1_lookup ',len1_lookup,  &
           ' len2_lookup ',len2_lookup,  &
           ' len_lookup  ',len_lookup

!     -------------------------------
!     Allocate space for LookUp table
!     -------------------------------

      ALLOCATE ( LookUp (Len1_Lookup, Len2_Lookup) )

!     -----------------------------
!     Move to start of LookUp table
!     -----------------------------

      CALL SetPos ( FF_unit_no, FixHd(150)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for LookUp table'

        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     --------------------
!     Read in LookUP table
!     --------------------

      CALL BuffIn ( FF_unit_no,                           &
           LookUp, Len_Lookup, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_lookup) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR('buffer in of lookup header',A_IO,LEN_IO,  &
             len_lookup)
        CMESSAGE = 'I/O ERROR with buffin of Lookup Header'
        ICode    = 11

        CALL EReport (RoutineName, ICode, CMessage)
      END IF

!     -----------------------------
!     Find number of LookUp headers
!     -----------------------------

      DO i = 1, FixHd(152)
        IF ( Lookup (42,i) == -99 ) THEN
          EXIT
        ELSE
          n_lookups(i_FF) = n_lookups(i_FF) +1
        END IF
      END DO

!     Update total no of LookUp Entries
      n_lookups(0) = n_lookups(0) + n_lookups(i_FF)

      WRITE (6,*) ' n_lookups ',n_lookups(i_FF)

!     Deallocate LookUp Table for File i_FF
      DEALLOCATE (LookUp)

    END DO   !  Loop over i_FF

    WRITE (6,*) ' total n_lookups ',n_lookups(0)

!   ----------------------------------------
!   Total number of LookUp entries now known
!   Allocate space for all LookUp entries
!   ----------------------------------------

    ALLOCATE ( LookUp (64, n_lookups(0) ) )
    ALLOCATE ( LookUp_FF_No (n_lookups(0)) )

    ihead = 0

    DO i_FF = 1, n_FF

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no
      WRITE (6,*) ' '

!     -----------------------------
!     Move to start of LookUp table
!     -----------------------------

      CALL SetPos ( FF_unit_no, LookUp_Start_Address(i_FF)-1, ICode )

      IF (ICode /= 0) THEN
        CMessage = 'Error in SetPos for LookUp table'

        CALL Ereport ( RoutineName, icode, CMessage )
      END IF

!     --------------------
!     Read in LookUP table
!     --------------------

      Len_LookUp = 64 * n_lookups(i_FF)

      CALL BuffIn ( FF_unit_no,                                      &
           LookUp (1:,ihead+1:), Len_Lookup, Len_io, a_io )

      IF (A_IO /= -1.0 .OR. LEN_IO /= len_lookup) THEN
! DEPENDS ON: IOERROR
        CALL IOERROR('buffer in of lookup header',A_IO,LEN_IO,       &
             len_lookup)
        CMESSAGE = 'I/O ERROR with buffin of Lookup Header'
        ICode    = 11

        CALL EReport (RoutineName, ICode, CMessage)
      END IF

      DO i = 1, n_lookups (i_FF)
        LookUp_FF_No (ihead+i) = i_FF
      END DO

      ihead = ihead + n_lookups (i_FF)

    END DO

!   -----------------------------
!   Allocate space for input data
!   -----------------------------

    ALLOCATE ( Data_In  (row_length, rows, nlevs) )
    ALLOCATE ( Srce_Lat (rows) )
    ALLOCATE ( Srce_Long(row_length) )

!   -----------------------------------------------
!   Compute icol and jrow from lat/long in namelist
!   -----------------------------------------------

    WRITE (6,*) ' Source Grid details'
    WRITE (6,*) ' No of columns   ',IntHead(6)
    WRITE (6,*) ' No of rows      ',IntHead(7)
    WRITE (6,*) ' deltax          ',RealHead(1)
    WRITE (6,*) ' deltay          ',RealHead(2)
    WRITE (6,*) ' first latitude  ',RealHead(3)
    WRITE (6,*) ' first longitude ',RealHead(4)

    DO iprof = 1, n_profiles

      WRITE (6,*) ' latitude  in namelist ',latitude(iprof)
      WRITE (6,*) ' longitude in namelist ',longitude(iprof)

      prof_lat  = latitude (iprof)
      prof_long = longitude(iprof)

      !       WRITE (6,*) ' prof_lat  ',prof_lat
      !       WRITE (6,*) ' prof_long ',prof_long

      DO i = 1, Inthead(6)
        srce_long(i) = RealHead(4) + (i-1)*RealHead(1)
      END DO

      DO j = 1, IntHead(7)
        srce_lat(j) = RealHead(3) + (j-1)*RealHead(2)
      END DO

      ix = 0
      DO i = 1, IntHead(6)
        !         if (i < 10) WRITE (6,*) i,srce_long(i),prof_long
        IF ( srce_long(i) <= prof_long ) THEN
          ix = i
        ELSE
          EXIT
        END IF
      END DO
      !       WRITE (6,*) ' ix ',ix,' long ',srce_long(ix)

      iy = 0
      DO j = 1, IntHead(7)
        IF ( srce_lat(j) <= prof_lat ) THEN
          iy = j
        ELSE
          EXIT
        END IF
      END DO
      !       WRITE (6,*) ' iy ',iy,' lat ',srce_lat(iy)

      icol (iprof) = ix
      irow (iprof)  = iy
      WRITE (6,*) ' icol set to ',icol (iprof)
      WRITE (6,*) ' irow set to ',irow (iprof)

    END DO

!   ---------------
!   Loop over times
!   ---------------

    DO itime = 1, n_times

!     ----------------------------
!     Find pointers to data fields
!     ----------------------------

      IF (itime == 1) THEN
        ipt_orog = 0
        DO i = 1, n_lookups (0)
          IF ( Lookup(42,i) == 33 .AND.               &
               LookUp(14,i) == FC_Time(itime) ) THEN
            ipt_orog = i
            EXIT
          END IF
        END DO
      END IF

      ipt_u = 0
      DO i = 1, n_lookups (0)
        IF ( Lookup(42,i) == 2 .AND.                  & 
             LookUp(14,i) == FC_Time(itime) ) THEN
          ipt_u = i
          EXIT
        END IF
      END DO

      ipt_v = 0
      DO i = 1, n_lookups (0)
        IF ( Lookup(42,i) == 3 .AND.                  & 
             LookUp(14,i) == FC_Time(itime) ) THEN
          ipt_v = i
          EXIT
        END IF
      END DO

      ipt_theta = 0
      DO i = 1, n_lookups (0)
        IF ( Lookup(42,i) == 4 .AND.                  &
             LookUp(14,i) == FC_Time(itime) ) THEN
          ipt_theta = i
          EXIT
        END IF
      END DO

      ipt_density = 0
      DO i = 1, n_lookups (0)
        IF ( Lookup(42,i) == 253 .AND.                & 
             LookUp(14,i) == FC_Time(itime) ) THEN
          ipt_density = i
          EXIT
        END IF
      END DO

      ipt_rho_p = 0
      DO i = 1, n_lookups (0)
        IF ( Lookup(42,i) == 407 .AND.                &
             LookUp(14,i) == FC_Time(itime) ) THEN
          ipt_rho_p = i
          EXIT
        END IF
      END DO

      WRITE (6,*) ' ipt_orog    ',ipt_orog
      WRITE (6,*) ' ipt_u       ',ipt_u
      WRITE (6,*) ' ipt_v       ',ipt_v
      WRITE (6,*) ' ipt_theta   ',ipt_theta
      WRITE (6,*) ' ipt_density ',ipt_density
      WRITE (6,*) ' ipt_rho_p   ',ipt_rho_p

!     ----------------------------------------------
!     Extract validity time from LookUp header for U
!     ----------------------------------------------

      IF (ipt_u > 0) THEN
        vt_yy = LookUp (1, ipt_u)
        vt_mm = LookUp (2, ipt_u)
        vt_dd = LookUp (3, ipt_u)
        vt_hh = LookUp (4, ipt_u)
      END IF
      WRITE (6,*) ' VT : ',vt_yy,vt_mm,vt_dd,vt_hh

!     ---------------------
!     Get FieldFiles Number
!     ---------------------

      FF_Unit_No = 30 + LookUp_FF_No (ipt_u)
      WRITE (6,*) ' FF_Unit_No : ',FF_Unit_No

!     -----------------
!     Read in orography
!     -----------------

      IF (itime == 1) THEN   !   Read in orography for first time only

        WRITE (6,*) ' Reading in orography'

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,          &
                              1,                   &
                              ipt_orog,            &  
                              lookup,              &
                              64,                  &
                              data_in(1,1,1),      &
                              fixhd,               &
                              1,                   &
                              icode,               &
                              CMessage)        

        IF (icode /= 0) THEN
          WRITE (6,*) 'Error in ReadFlds_Serial ?'

          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

!     -----------------
!     Extract Orography
!     -----------------

        DO iprof = 1, n_profiles
          orog (iprof) = data_in ( icol(iprof), irow(iprof), 1 )
          WRITE (6,*) ' orog ',orog(iprof)
        END DO

      END IF

!     --------------
!     Read in U data
!     --------------

      !     ------------------------------------------
      !     ReadFlds_Serial can only read in one level
      !     at a time so must call for each level
      !     ------------------------------------------

      WRITE (6,*) ' Reading in U'

      DO k = 1, nlevs

        !       WRITE (6,*) ' Calling ReadFlds_Serial for level ',k

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,        &
                              1,                 &
                              ipt_u+k-1,         & 
                              lookup,            &
                              64,                &
                              data_in(1,1,k),    &
                              fixhd,             &
                              1,                 &
                              icode,             &
                              CMessage)

        IF (icode /= 0) THEN
          WRITE (6,*) 'Error in ReadFlds_Serial ?'

          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

      END DO

!     ----------------
!     Extract U data
!     ----------------

      DO iprof = 1, n_profiles
        DO k=1,nlevs
          u(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
          !       WRITE (6,*) ' u data ',k,u(k)
        END DO
      END DO

!     --------------
!     Read in V data
!     --------------

      WRITE (6,*) ' Reading in V'

      DO k = 1, nlevs

        !       WRITE (6,*) ' Calling ReadFlds_Serial for level ',k

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,        &
                              1,                 &
                              ipt_v+k-1,         &
                              lookup,            & 
                              64,                &
                              data_in(1,1,k),    &
                              fixhd,             &
                              1,                 &
                              icode,             &
                              CMessage)

        IF (icode /= 0) THEN
          WRITE (6,*) 'Error in ReadFlds_Serial ?'

          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

      END DO

!     ----------------
!     Extract V data
!     ----------------

      DO iprof = 1, n_profiles
        DO k=1,nlevs
          v(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
          !       WRITE (6,*) ' v data ',k,v(k)
        END DO
      END DO

!     --------------
!     Read in P data
!     --------------

      WRITE (6,*) ' Reading in P'

      DO k = 1, nlevs

        !       WRITE (6,*) ' Calling ReadFlds_Serial for level ',k

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,        &
                              1,                 &
                              ipt_rho_p+k-1,     &
                              lookup,            &
                              64,                &
                              data_in(1,1,k),    &
                              fixhd,             &
                              1,                 &
                              icode,             &
                              CMessage)

        IF (icode /= 0) THEN
          WRITE (6,*) 'Error in ReadFlds_Serial ?'

          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

      END DO

!     --------------
!     Extract P data
!     --------------

      DO iprof = 1, n_profiles
        DO k=1,nlevs
          p(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
          !       WRITE (6,*) ' p data ',k,p(k)
        END DO
      END DO

!     -------------
!     Read in Theta
!     -------------

      WRITE (6,*) ' Reading in theta'

      DO k = 1, nlevs

        !       WRITE (6,*) ' Calling ReadFlds_Serial for level ',k

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,         &
                              1,                  &
                              ipt_theta+k-1,      &  
                              lookup,             &
                              64,                 &
                              data_in(1,1,k),     &
                              fixhd,              &
                              1,                  &
                              icode,              &
                              CMessage)

        IF (icode /= 0) THEN
          WRITE (6,*) 'Error in ReadFlds_Serial ?'

          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

      END DO

!     ------------------
!     Extract Theta data
!     ------------------

      DO iprof = 1, n_profiles
        DO k=1,nlevs
          theta(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
          !       WRITE (6,*) ' theta data ',k,theta(k)
        END DO
      END DO

!     ---------------
!     Read in density
!     ---------------

      WRITE (6,*) ' Reading in density'

      DO k = 1, nlevs

        !       WRITE (6,*) ' Calling ReadFlds_Serial for level ',k

! DEPENDS ON: ReadFlds_Serial
        CALL ReadFlds_Serial (FF_unit_no,         &
                              1,                  &
                              ipt_density+k-1,    &
                              lookup,             & 
                              64,                 &
                              data_in(1,1,k),     &
                              fixhd,              &
                              1,                  &
                              icode,              &
                              CMessage)

        IF (icode /= 0) THEN
          WRITE (6,*) 'Error in ReadFlds_Serial ?'

          CALL Ereport ( RoutineName, icode, CMessage )
        END IF

      END DO

!     --------------------
!     Extract Density data
!     --------------------

      DO iprof = 1, n_profiles
        DO k=1,nlevs
          rho(k,iprof) = data_in ( icol(iprof), irow(iprof), k )
          !       WRITE (6,*) ' rho data ',k,rho(k)
        END DO
      END DO

!     --------------------------
!     Write out the profile data
!     --------------------------

      DO i = 1, n_profiles  !  Points

        WRITE (prof_unit_no,1)  &
             'YEAR MONTH DAY HOUR         FORECAST TIME (H)'
        WRITE (prof_unit_no,2) vt_yy,vt_mm,vt_dd,vt_hh,fc_time(itime)
        WRITE (prof_unit_no,3)  & 
             'LATITUDE LONGITUDE  ROW  COLUMN   HEIGHT  LEVELS'
        WRITE (prof_unit_no,4)  &
             latitude(i),longitude(i),irow(i),icol(i),orog(i),nlevs
        WRITE (prof_unit_no,5)  &
             ' Level    U         V        Theta      P        Rho'

        DO k = 1, nlevs     !  Levels

          !         WRITE (6,6)                            &
          !         z(k),u(k),v(k),theta(k),p(k),rho(k)    &
          !         k,u(k),v(k),theta(k),p(k),rho(k)       &
          !         k,u(k,i),v(k,i),theta(k,i),p(k,i),rho(k,i)

          WRITE (prof_unit_no,6)  &
               k,u(k,i),v(k,i),theta(k,i),p(k,i),rho(k,i)

        END DO

      END DO

    END DO    !  End of loop over forecast times

!   -------------------------
!   Deallocate space for data
!   -------------------------

    DEALLOCATE (Data_In)
    DEALLOCATE (Srce_Lat)
    DEALLOCATE (Srce_Long)

    DEALLOCATE (LookUp)
    DEALLOCATE (LookUp_FF_No)

!   ---------------------
!   Close the fieldsfiles
!   ---------------------

    WRITE (6,*) ' closing fieldsfiles'
    DO i_FF = 1, n_FF

      ff_unit_no = 30 + i_FF

      env = 'FILE  '
      WRITE (env(5:6),'(I2)') ff_unit_no
      WRITE (6,*) ' '

      CALL File_Close (ff_unit_no,env,len_env,env_var,0,icode)
      WRITE (6,*) ' icode ',icode
      IF (icode /= 0) THEN
        CMessage = 'Error in closing fieldsfiles'

        CALL Ereport ( RoutineName, Icode, Cmessage )
      END IF
    END DO

!   ---------------------
!   Close the output file
!   ---------------------

!   Cannot use file_close, use fortran close

    CLOSE ( prof_unit_no, iostat=icode )

    WRITE (6,*) ' icode ',icode
    IF (icode /= 0) THEN
      CMessage = 'Error in closing output file'

      CALL Ereport ( RoutineName, Icode, Cmessage )
    END IF

!   ----------------------------------------------------
!   These FORMAT statements must not be changed without
!   making corresponding changes in the 3DVOM model
!   ----------------------------------------------------

1   FORMAT (a45)
2   FORMAT (i4,1x,i2,4x,i2,2x,i2,10x,i2)
3   FORMAT (a48)
4   FORMAT (f7.3,2x,f7.3,2x,i5,2x,i5,2x,f8.3,2x,i3)
5   FORMAT (a52)
6   FORMAT (i3,1x,2(f10.4,1x),1x,f8.3,1x,f10.2,2x,f17.2)

    RETURN
  END SUBROUTINE UM_Profiles
END PROGRAM VOM_Extract
