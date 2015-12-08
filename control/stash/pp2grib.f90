! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE PP2GRIB------------------------------------------------
!
!  Purpose:
!   to code pp_header and un-packed data into grib
!
!  Programming standard: Unified Model Documentation Paper No 3
!
!  System component:
!
!  System task:
!
!  Documentation:
!
! Arguments:-------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
      SUBROUTINE PP2GRIB(FIELD,WORK_ARRAY,MAXDIM,NUM_CRAY_WORDS,        &
     &                   GRIB_PACKING,ILABEL,RLABEL,                    &
     &                   ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER                                                           &
     &     MAXDIM                                                       &
     &    ,NUM_CRAY_WORDS                                               &
     &    ,ILABEL(45)                                                   &
     &    ,WORK_ARRAY(MAXDIM)                                           &
     &    ,GRIB_PACKING                                                 &
                                 ! IN - profile for packing
     &    ,ICODE                 ! out - error code
      REAL                                                              &
     &     FIELD(MAXDIM)                                                &
     &    ,RLABEL(19)
      CHARACTER(LEN=80)                                                      &
     &     CMESSAGE              ! out - error message

!     LOCAL VARIABLES
      INTEGER IDIM                                                      &
     &       ,BLOCK0(4)                                                 &
     &       ,BLOCK1(21)                                                &
     &       ,BLOCK2(30)                                                &
     &       ,BLOCK3(2)                                                 &
     &       ,BLOCK4(2)                                                 &
     &       ,BITMAP(MAXDIM)                                            &
     &       ,WORK1(288)                                                &
     &       ,WORK2(500)                                                &
     &       ,ERROR                                                     &
     &       ,ERROR_UNIT                                                &
     &       ,QUASI(1)                                                  &
     &       ,LENVRT                                                    &
     &       ,WIDTH                                                     &
     &       ,WORDSZ                                                    &
     &       ,LENQ                                                      &
     &       ,LENGRB
      INTEGER                                                           &
     &        STASH_SECT_NO                                             &
     &       ,STASH_ITEM_NO                                             &
     &       ,TABLE2_VERSION                                            &
     &       ,TABLE2_ENTRY                                              &
     &       ,IB ,IC, BBB, ICENTURY                                     &
     &       ,D_TIME, T_UNITS                                           &
     &       ,MSG_LVL                                                   &
     &       ,IFLAG_MAX, IFLAG_MIN, IFLAG_VERTM                         &
                                                  ! flags for processing
     &       ,IFLAG_MEAN ,IFLAG_ZONAL                                   &
                                                  ! code
     &       ,IREM
      LOGICAL OROW                                                      &
     &       ,OBITMAP
      REAL VERTCO(4)                                                    &
     &    ,WORKR(288)                                                   &
     &    ,BLOCKR(20)                                                   &
     &    ,STORE(MAXDIM)                                                &
     &    ,FMAX ,FMIN, RANGE                                            &
     &    ,DLONZ

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      PARAMETER(MSG_LVL=2)  ! Errors only
      INTEGER :: icnt
      INTEGER :: ii, ij
      INTEGER :: len_bitmap


! --------------------------------------------------------------------
! initialise variables for call to coder

      IF (lhook) CALL dr_hook('PP2GRIB',zhook_in,zhook_handle)
      ICODE=0
      ERROR=0
      ERROR_UNIT=6
      LENQ   = 1
      WORDSZ = 64
      LENGRB = MAXDIM

      IDIM   = ILABEL(18)*ILABEL(19)      ! length of field

! --------------------------------------------------------------------
! Method of grib packing
! -----------------------
! PPXREF profiles 1-5 - use binary accuracy method (requires less space)
!                  6  - use width method, with simple packing to be
!                       the similar to the ECMWF MARS archive.
!
! Note in the case of -99 in a binary accuracy profile ie no packing,
!  use width =30, the maximum value which can safely be used with 32
!  bit numbers.
!
      IF (GRIB_PACKING == 6) THEN  ! Width method

        OROW = .FALSE.             ! Simple packing
        OBITMAP = .FALSE.          ! Do not use bitmap.
        BLOCK0(4)=0
        WIDTH = NINT(RLABEL(6))    ! required width
        IF (WIDTH >  30.or.WIDTH <= 0) THEN  ! check sensible value
          WIDTH=30
        ENDIF
      ELSE                         ! binary accuracy method

       OROW = .TRUE.               ! Row by row packing
       OBITMAP = .TRUE.            ! use bitmap for fields with
                                   ! missing data.
       BLOCK0(4)  = NINT(RLABEL(6)) ! Binary accuracy - power of 2
       WIDTH=0
       IF (BLOCK0(4) == -99) THEN  ! use maximum number of bits
          WIDTH = 30
       ENDIF
      ENDIF

! --------------------------------------------------------------------
!     SET UP VARIABLES FOR GRIB CODING ROUTINE
!
!  Section 0
! -------------
!
      BLOCK0(1)  = 1               ! GRIB Edition number
      STASH_SECT_NO = ILABEL(42)/1000
      STASH_ITEM_NO = MOD(ILABEL(42),1000)
! DEPENDS ON: stash_grib
      CALL STASH_GRIB(STASH_SECT_NO,STASH_ITEM_NO,                      &
     &                TABLE2_VERSION,TABLE2_ENTRY,ERROR)
      BLOCK0(2)  = TABLE2_VERSION  ! table 2 version number
      BLOCK0(3)  = 0               !length of message (OUTPUT ONLY)
!     BLOCK0(4)   set above

!  Section 1
! ------------

      BLOCK1(1)  = 74              ! ORIGINATING CENTRE
      BLOCK1(2)  = 45              ! MODEL IDENT NUMBER
      BLOCK1(3)  = 42              ! Grid ident number

! Bit map

      IF (OBITMAP) THEN
        ICNT=0
!       WRITE(6,*)' @@ MISSING DATA BIT-MAPPING  rmdi ',RLABEL(18)
        DO II=1,IDIM
          IF (FIELD(II)  /=  RLABEL(18)) THEN
            ICNT=ICNT+1
            STORE(ICNT)=FIELD(II)
            BITMAP(II)=1
          ELSE
            BITMAP(II)=0
          END IF
        ENDDO
        LEN_BITMAP  = IDIM
        IDIM        = ICNT
        IF (IDIM  /=  LEN_BITMAP) then   ! bitmap required
          BLOCK1(4)   = 192           ! Block ident flags
        ELSE                             ! no bitmap required
          LEN_BITMAP  = 1
          BLOCK1(4)   = 128           ! Block ident flags
        ENDIF
      ELSE
!
! Profile 6 - attempt to resemble grib in MARS archive
! replace missing data values by 0.0 - not an ideal solution
! but loose all accuaracy if keep UM missing data indicator.
!
        DO IJ=1,IDIM
          STORE(IJ)=FIELD(IJ)
          IF (field(ij) == RLABEL(18)) store(IJ)=0.0
        ENDDO
        LEN_BITMAP  = 1
        BLOCK1(4)   = 128           ! Block ident flags
      END IF

      BLOCK1(5)  = TABLE2_ENTRY      !parameter identification
! ------------------------------------------------------------------
!  Make use of LBPROC information  ILABEL(25)
!
!  Note not worked out a way of storing max and min information in
!  grib header

      IFLAG_MAX=0
      IFLAG_MIN=0
      IFLAG_VERTM=0
      IFLAG_MEAN=0
      IFLAG_ZONAL=0
      IREM=ILABEL(25)
      IF(IREM >= 8192) THEN     ! maximum value
        IFLAG_MAX=1
        IREM=IREM-8192
      ENDIF
      IF(IREM >= 4096) THEN     ! minimum value
        IFLAG_MIN=1
        IREM=IREM-4096
      ENDIF
      IF(IREM >= 2048) THEN     !  vertical mean
        IFLAG_VERTM=1
        IREM=IREM-2048
      ENDIF
! 1024 - difference between fields at 2 levels (not used in UM)
! 512 - Square root of field ( not used in UM)
! 256 - product of fields not used in UM output
      IF(IREM >= 128) THEN     ! time mean
        IFLAG_MEAN=1
        IREM=IREM-128
      ENDIF
      IF(IREM >= 64) THEN     ! Zonal mean
        IFLAG_ZONAL=1
      ENDIF
! 32, 16, 8, 4, 2 & 1 not used in UM output
! ------------------------------------------------------------------
!
! Level type information

      IF (ILABEL(26)  ==  9) THEN
        IF (IFLAG_VERTM == 1) THEN
          BLOCK1(6) = 110       ! vertical mean hybrid coordinates
                                ! Using code for layer information
        ELSE
          BLOCK1(6) = 109            ! Hybrid coordinates
        ENDIF
      ELSE IF (ILABEL(26)  ==  8) THEN
        IF (IFLAG_VERTM == 1) THEN
          BLOCK1(6) = 121       ! vertical mean pressure coordinates
                                ! Using code for layer information
        ELSE
          BLOCK1(6) = 100              ! pressure coordinates
        ENDIF
      ELSE IF (ILABEL(26)  ==  1) THEN
        BLOCK1(6) = 105              ! Height coordinates
      ELSE IF (ILABEL(26)  ==  128) THEN
        BLOCK1(6) = 102              ! Mean sea level pressure
      ELSE IF (ILABEL(26)  ==  129) THEN
        BLOCK1(6) = 1                ! surface
      ELSE IF (ILABEL(26)  ==  130) THEN
        BLOCK1(6) = 7                ! tropopause level
      ELSE IF (ILABEL(26)  ==  131) THEN
        BLOCK1(6) = 6                ! Max wind
      ELSE IF (ILABEL(26)  ==  132) THEN
        BLOCK1(6) = 4                ! Freezing level ?
      ELSE IF (ILABEL(26)  ==  10) THEN
        BLOCK1(6) = 107              ! Sigma coordinates
      ELSE IF (ILABEL(26)  ==  6) THEN
        BLOCK1(6) = 111              ! depth below land surface
                                     ! used for soil levels
      ELSE IF (ILABEL(26)  ==  133) THEN ! top of atmosphere
        BLOCK1(6) = 8                ! nominal top of atmosphere
      ELSE IF (ILABEL(26)  ==  275) THEN ! canopy height
        BLOCK1(6) = 1                ! At present redefine as surface
      ELSE IF (ILABEL(26)  ==  0) THEN
        BLOCK1(6) = 0                ! Unspecified
      ELSE
        CMESSAGE='PP2GRIB : unrecognised level coordinate'
      END IF

! Additional level information

      IF (IFLAG_VERTM == 1) THEN
        IF (ILABEL(26) == 9) THEN
! Note pp headers may not contain top level info as they should
          BLOCK1(7) = 19             ! fixed at present
          BLOCK1(8) =ILABEL(33)      ! bottom level number
        ELSE IF (ILABEL(26) == 8) THEN
          BLOCK1(7) = NINT(1100. - RLABEL(8)) ! Top pressure
          BLOCK1(8) = NINT(1100. - RLABEL(7)) ! bottom pressure
        ENDIF
      ELSE
        IF (ILABEL(26) == 9) THEN
          BLOCK1(7) =ILABEL(33)      ! model level number
          BLOCK1(8) =0
        ELSEIF (ILABEL(26) == 8) THEN
          BLOCK1(7) =ILABEL(33)      ! Pressure in hPa
          BLOCK1(8) =0
        ELSE
          BLOCK1(7)   = 0             ! Level descriptor
          BLOCK1(8)   = 0             !   "       "     (overflow)
        ENDIF
      ENDIF

! ------------------------------------------------------------------
! Time and date information
!
!  First use LBTIM to determine how to work out time date

      BBB=MOD(ILABEL(13),100)
      IC=MOD(ILABEL(13),10)
      IB=(BBB-IC)/10

! work out century using ilabel(7)
      ICENTURY=(ILABEL(7)-1)/100 + 1

! 1. Model time no year and month
! ---------------------------------
      IF (IC == 0) THEN
         CMESSAGE='PP2GRIB : cannot code date/time '
         ICODE=1
         WRITE(6,*)'PP2GRIB: not able to code at present'

! 2. Normal 365 day calendar
! --------------------------
! At present assumes all fields are forecasts not means

      ELSE IF(IC == 1) THEN

! a) normal forecasts less than 10 days

        IF (ILABEL(14) <  256) THEN
          BLOCK1(17) = 0               ! Time range indicator

! b) Copes with periods up to 65535 hours (2730 days or 7 years)

        ELSE IF (ILABEL(14) >= 256.and.ILABEL(14) <  65535) THEN
          BLOCK1(17) = 10              ! uses two octets for P1

! c) Periods longer than 7 years
        ELSE
         CMESSAGE='PP2GRIB : cannot code forecast period '
         ICODE=1
         WRITE(6,*)'PP2GRIB: cannot code forecast period ',ilabel(14)
        ENDIF

          BLOCK1(9)  = ILABEL(7)-(ICENTURY-1)*100     ! year
          BLOCK1(10) = ILABEL(8)       ! Month
          BLOCK1(11) = ILABEL(9)       ! Day
          BLOCK1(12) = ILABEL(10)      ! hour
          BLOCK1(13) = 0               ! minute
          BLOCK1(14) = 1               ! time unit
          BLOCK1(15) = ILABEL(14)      ! P1 (F/C period in hours)
          BLOCK1(16) = 0               ! P2
          BLOCK1(18) = 0               ! number of averages
          BLOCK1(19) = ICENTURY        ! Century of reference time

! 3. 360 day year - normal climate run calendar
! ---------------------------------------------

      ELSE IF (IC == 2) THEN

        IF (ib == 1) THEN     ! forecast fields
! As many climate runs exceed 7 years all forecast periods are
! recoded as analyses.

           BLOCK1(17) = 0               ! Time range indicator
!
           BLOCK1(9)  = ILABEL(1)-(ICENTURY-1)*100  ! Year in ref cent
           BLOCK1(10) = ILABEL(2)       ! Month
           BLOCK1(11) = ILABEL(3)       ! Day
           BLOCK1(12) = ILABEL(4)       ! hour
           BLOCK1(13) = ILABEL(5)       ! minute
           BLOCK1(14) = 1               ! time unit
           BLOCK1(15) = 0               ! P1
           BLOCK1(16) = 0               ! P2
           BLOCK1(18) = 0               ! number of averages
           BLOCK1(19) = ICENTURY        ! Century of reference time

        ELSE IF (ib ==  2) THEN       ! Time average

           BLOCK1(17) = 3               ! Time range indicator
           BLOCK1(9)  = ILABEL(1)-(ICENTURY-1)*100    ! year
           BLOCK1(10) = ILABEL(2)       ! Month
           BLOCK1(11) = ILABEL(3)       ! Day
           BLOCK1(12) = ILABEL(4)       ! hour
           BLOCK1(13) = ILABEL(5)       ! minute
! Work out p2 and appropriate units for p2
! DEPENDS ON: grib_time_int
           CALL GRIB_TIME_INT(ILABEL(1),ILABEL(2),ILABEL(3),ILABEL(4),  &
     &              ILABEL(5),ILABEL(7),ILABEL(8),ILABEL(9),ILABEL(10), &
     &              ILABEL(11),.TRUE.,D_TIME,T_UNITS)

           BLOCK1(14) = T_UNITS         ! time units
           BLOCK1(16) = D_TIME          ! P2
           BLOCK1(15) = 0               ! P1
           BLOCK1(18) = 1               ! number of averages
           BLOCK1(19) = ICENTURY        ! Century of reference time

         else    ! cannot code
           WRITE(6,*)'Not able to code for date type ib =',ib
           CMESSAGE='PP2GRIB: error for date time'
           ICODE=1
         ENDIF
      ELSE
          WRITE(6,*)'Not able to code for date type ic =',ic
          CMESSAGE='PP2GRIB: error for date time'
          ICODE=1
      ENDIF

      BLOCK1(20) = 0               ! Decimal scale factor
      BLOCK1(21) = 21              ! Length of BLOCK1

! --------------------------------------------------------------------
! Section 2 information
! ----------------------
      IF (BLOCK1(6)  ==  109) THEN
        BLOCK2(1) = 2                ! Number of vert coord parms
        BLOCK2(2) = 53               ! PV, PL or 255
        LENVRT   = 2
        VERTCO(1) = RLABEL(9)       ! A COORDINATE
        VERTCO(2) = RLABEL(7)        ! B COORDINATE
!  Can be uncommented when values for rlabel(1) & rlabel(2) are coded
!     ELSE IF (BLOCK1(6)  ==  110) THEN
!       BLOCK2(1) = 4                ! Number of vert coord parms
!       BLOCK2(2) = 53               ! PV, PL or 255
!       LENVRT   = 4
!       VERTCO(1) = RLABEL(10)      ! A COORDINATE top
!       VERTCO(2) = RLABEL(8)       ! B COORDINATE  top
!       VERTCO(3) = RLABEL(1)       ! A COORDINATE bottom (not coded)
!       VERTCO(4) = RLABEL(2)       ! B COORDINATE bottom (not coded)
      ELSE
        BLOCK2(1) = 0                ! Number of vert coord parms
        BLOCK2(2) = 255              ! PV, PL or 255
        LENVRT   = 1
        VERTCO(1) = 0.               ! A COORDINATE
        VERTCO(2) = 0.               ! B COORDINATE
      END IF

      BLOCK2(3) = 0                ! representation type ie lat longrid
      BLOCK2(4) = ILABEL(19)       ! Number of cols
      BLOCK2(5) = ILABEL(18)       ! Number of rows
      BLOCK2(6) = NINT((RLABEL(14)+RLABEL(15))*1000) ! Lat 1st pt.
      IF (BLOCK2(6) >  180000.0) THEN
        BLOCK2(6) = BLOCK2(6)-180000.0
      END IF
      IF (IFLAG_ZONAL == 1) THEN        ! zonal means
        BLOCK2(7) = NINT(RLABEL(16)*1000)              ! Lon 1st pt.
        BLOCK2(10) = NINT(RLABEL(17)*1000)             ! lon last pt.
        BLOCK2(11) = ABS(NINT((RLABEL(17)-RLABEL(16))*1000)) ! dlon
      ELSE
        BLOCK2(7) = NINT((RLABEL(16)+RLABEL(17))*1000) ! Lon 1st pt.
        BLOCK2(10) = NINT((RLABEL(16)+(ILABEL(19)*RLABEL(17)))*1000)
!                          ! Lon of extreme point
        BLOCK2(11) = ABS(NINT(RLABEL(17)*1000))
!                             ! Horizontal dirn increment
      ENDIF
      IF (BLOCK2(7) >  360000.0) THEN
        BLOCK2(7) = BLOCK2(7)-360000.0
      END IF
      IF (BLOCK2(10) >  360000.0) THEN
        BLOCK2(10) = BLOCK2(10)-360000.0
      END IF
      BLOCK2(8) = 128        ! resolution and component flags
      BLOCK2(9) = NINT((RLABEL(14)+(ILABEL(18)*RLABEL(15)))*1000)
!                                  Lat of extreme point
      IF (BLOCK2(9) >  180000.0) THEN
        BLOCK2(9) = BLOCK2(9)-180000.0
      END IF
      BLOCK2(12) = ABS(NINT(RLABEL(15)*1000))
!                          ! Vertical dirn increment
! Scanning mode flags
!  west to east is positive
!   If grid scans from west to east bit 1 is 0
!   if grid scans from east to west bit 1 is 1  ie add 128
!  south to north is positive
!   if grid scans from north to south bit 2 is 0
!   if grid scans from south to north bit 2 is 1 ie add 64

      BLOCK2(13) = 0                     ! Scanning mode flags
      IF (IFLAG_ZONAL == 1) THEN
        DLONZ=RLABEL(17)-RLABEL(16)
        IF (DLONZ  <   0.0) BLOCK2(13) = BLOCK2(13) +128
      ELSE
        IF (RLABEL(17)  <   0.0) BLOCK2(13) = BLOCK2(13) + 128
      ENDIF
      IF (RLABEL(15)  >   0.0) BLOCK2(13) = BLOCK2(13) + 64

      BLOCK2(14) = -NINT(RLABEL(11)*1000) ! Lat S Pole
      BLOCK2(15) = -NINT(RLABEL(12)*1000) ! Lon S Pole
      BLOCK2(16) = 0
      BLOCK2(17) = 0
      BLOCK2(18) = 0
      BLOCK2(19) = 0
      BLOCK2(20) = 0

! Section 3
! ----------

      BLOCK3(1)=0
      BLOCK3(2)=0

! Section 4
! ----------

      IF (OROW) THEN
        BLOCK4(1)  = 80  ! row by row packing
      ELSE
        BLOCK4(1)  = 0   ! simple packing
      END IF
      BLOCK4(2)  = 0

! -------------------------------------------------------------------

! Call grib encoder

      IF (ICODE == 0) THEN

! DEPENDS ON: coder

       CALL CODER(STORE,IDIM,VERTCO,LENVRT,BITMAP,LEN_BITMAP,QUASI,LENQ,&
     &           WIDTH,WORDSZ,BLOCK0,BLOCK1,BLOCK2,BLOCK3,BLOCK4,       &
     &           BLOCKR,WORK_ARRAY,LENGRB,NUM_CRAY_WORDS,               &
     &           ERROR,WORK1,WORK2,WORKR,ERROR_UNIT,MSG_LVL)


        IF (WIDTH >  30) then
          CMESSAGE='PP2GRIB: trying to use more than 30 bits for grib'
          ICODE=0             ! don't enforce failure at the moment
      WRITE(6,*)'WARNING: grib requires more than 30 bits for accuracy',&
     &            WIDTH,' stash code ',ILABEL(42)
        ENDIF
      ELSE
        WRITE(6,*)'PP2GRIB: CODER not called for field ',ilabel(42)
        WRITE(6,*)CMESSAGE
      ENDIF

! set output length in header

      ILABEL(15) = NUM_CRAY_WORDS

      IF (lhook) CALL dr_hook('PP2GRIB',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PP2GRIB
