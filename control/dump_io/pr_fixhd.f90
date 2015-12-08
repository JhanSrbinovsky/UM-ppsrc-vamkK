! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_FIXHD---------------------------------------
!
!    Purpose: Prints out fixed length header record and checks
!             validity of information.
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Dump I/O

SUBROUTINE pr_fixhd                                               &
(fixhd,len_fixhd,len_inthd,len_realhd,len1_levdepc                &
,len2_levdepc,len1_rowdepc,len2_rowdepc,len1_coldepc,len2_coldepc &
,len1_flddepc,len2_flddepc,len_extcnst,len_dumphist,len_cfi1      &
,len_cfi2,len_cfi3,len1_lookup,len2_lookup,len_data               &
,icode,cmessage)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE PrintStatus_mod
IMPLICIT NONE

INTEGER                                                           &
 len_fixhd                                                        &
               !IN Length of fixed length header
,len_inthd                                                        &
               !IN Length of integer header
,len_realhd                                                       &
               !IN Length of real header
,len1_levdepc                                                     &
               !IN 1st dim of level dep consts
,len2_levdepc                                                     &
               !IN 2nd dim of level dep consts
,len1_rowdepc                                                     &
               !IN 1st dim of row dep consts
,len2_rowdepc                                                     &
               !IN 2nd dim of row dep consts
,len1_coldepc                                                     &
               !IN 1st dim of column dep consts
,len2_coldepc                                                     &
               !IN 2nd dim of column dep consts
,len1_flddepc                                                     &
               !IN 1st dim of field dep consts
,len2_flddepc                                                     &
               !IN 2nd dim of field dep consts
,len_extcnst                                                      &
               !IN Length of extra constants
,len_dumphist                                                     &
               !IN Length of history block
,len_cfi1                                                         &
               !IN Length of comp field index 1
,len_cfi2                                                         &
               !IN Length of comp field index 2
,len_cfi3                                                         &
               !IN Length of comp field index 3
,len1_lookup                                                      &
               !IN 1st dim of lookup
,len2_lookup   !IN 2nd dim of lookup

INTEGER                                                           &
 fixhd(len_fixhd)                                                 &
                  !IN Fixed length header
,len_data                                                         &
                  !IN Length of real data
,icode          !OUT Return code; successful=0
                !                 error > 0

CHARACTER(len=80)                                                 &
 cmessage       !OUT Error message if ICODE > 0

! -------------------------------------------------------------
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
! Local variables:---------------------------------------------
INTEGER i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!--------------------------------------------------------------

IF (lhook) CALL dr_hook('PR_FIXHD',zhook_in,zhook_handle)
icode=0
cmessage=' '

IF ( printstatus >= prstatus_oper ) THEN
WRITE(6,'('' '')')
WRITE(6,'('' FIXED LENGTH HEADER'')')
WRITE(6,'('' -------------------'')')

WRITE(6,'('' Dump format version'',I6)')fixhd(1)
WRITE(6,'('' UM Version No      '',I6)')fixhd(12)
END IF

IF(fixhd(2) == 1)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Atmospheric data'')')
  END IF
ELSE IF(fixhd(2) == 2)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Oceanic data'')')
  END IF
ELSE IF (fixhd(2) == 4) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Wave sub-model data'')')
  END IF
ELSE
  IF ( printstatus >= prstatus_min ) THEN
    WRITE(6,'(A,i9)')                                             &
    ' ***FATAL ERROR*** Invalid data type: FIXHD(2)=',fixhd(2)
  END IF
  icode=4
  cmessage='PR_FIXHD: Consistency check'
  IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
  RETURN
END IF

IF(fixhd(3) == 1)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' On hybrid levels'')')
  END IF
ELSE IF(fixhd(3) == 2)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' On sigma levels'')')
  END IF
ELSE IF(fixhd(3) == 3)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' On pressure levels'')')
  END IF
ELSE IF(fixhd(3) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' On depth levels'')')
  END IF
ELSE IF(fixhd(3) == 5)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Charney-Phillips on radius levels'')')
  END IF
ELSE IF (fixhd(3) == 6) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'(A)')                                                &  
    ' Wave model direction levels and frequency pseudo-levels'
  END IF
ELSE IF(fixhd(3) == imdi.AND.fixhd(5) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'(A)')                                                & 
    ' Missing data indicator used for vert coord type'
  END IF
ELSE
  IF ( printstatus >= prstatus_min ) THEN
    WRITE(6,'(A,I9)')                                                &
    ' ***FATAL ERROR*** Invalid level type: FIXHD(3)=',fixhd(3)
  END IF
  icode=4
  cmessage='PR_FIXHD: Consistency check'
  IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
  RETURN
END IF

IF(fixhd(4) == 0)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Over global domain'')')
  END IF
ELSE IF(fixhd(4) == 1)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Over N. Hemispheric domain'')')
  END IF
ELSE IF(fixhd(4) == 2)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Over S. Hemispheric domain'')')
  END IF
ELSE IF(fixhd(4) == 3)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Over LAM domain with no wrap around'')')
  END IF
ELSE IF(fixhd(4) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Over LAM domain with wrap around'')')
  END IF
ELSE IF(fixhd(4) == 103)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Over rotated LAM domain'')')
  END IF
ELSE
  IF ( printstatus >= prstatus_min ) THEN
    WRITE(6,'(A,I9)')                                             &
    ' ***FATAL ERROR*** Invalid domain: FIXHD(4)=',fixhd(4)
  END IF
  icode=4
  cmessage='PR_FIXHD: Consistency check'
  IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
  RETURN
END IF

IF(fixhd(5) == 1)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Instantaneous dump'')')
  END IF
ELSE IF(fixhd(5) == 2)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Meaned dump'')')
  END IF
ELSE IF(fixhd(5) == 3)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' FIELDS file'')')
  END IF
ELSE IF(fixhd(5) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Ancillary dataset'')')
  END IF
ELSE IF(fixhd(5) == 5)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Boundary dataset'')')
  END IF
ELSE IF(fixhd(5) == 6)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' AC Observation File'')')
  END IF
ELSE IF(fixhd(5) == 7)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Var Observation File'')')
  END IF
ELSE IF(fixhd(5) == 8)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Cx file (model columns at ob locations)'')')
  END IF
ELSE IF(fixhd(5) == 9)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Covariance File'')')
  END IF
ELSE IF (fixhd(5) == 10) THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE (6, '('' OPS Obstore file '')')
  END IF
ELSE
  IF ( printstatus >= prstatus_min ) THEN
    WRITE(6,'(A,I9)')                                             &
    ' ***FATAL ERROR*** Invalid dump type: FIXHD(5)=',fixhd(5)
  END IF
  icode=4
  cmessage='PR_FIXHD: Consistency check'
  IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
  RETURN
END IF

IF ( printstatus >= prstatus_oper ) THEN
WRITE(6,'('' Exp No ='',I6,'' Run Id ='',I6)') fixhd(7),fixhd(6)
END IF

IF(fixhd(8) == 1)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Gregorian calendar'')')
  END IF
ELSE IF(fixhd(8) == 2)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' 360-day calendar'')')
  END IF
ELSE IF(fixhd(8) == imdi.AND.fixhd(5) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Missing data indcator used as calendar'')')
  END IF
ELSE
  IF ( printstatus >= prstatus_min ) THEN
    WRITE(6,'(A,I9)')                                             &
    ' *** FATAL ERROR *** Invalid calendar type: FIXHD(8) =',fixhd(8)
  END IF
  icode=4
  cmessage='PR_FIXHD: Consistency check'
  IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
  RETURN
END IF

IF(fixhd(9) == 1)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Arakawa A grid'')')
  END IF
ELSE IF(fixhd(9) == 2)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Arakawa B grid'')')
  END IF
ELSE IF(fixhd(9) == 3)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Arakawa C grid'')')
  END IF
ELSE IF(fixhd(9) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Arakawa D grid'')')
  END IF
ELSE IF(fixhd(9) == 5)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Arakawa E grid'')')
  END IF
ELSE IF(fixhd(9) == 6)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' V-AT-POLES/ENDGAME grid'')')
  END IF
ELSE IF(fixhd(9) == imdi.AND.fixhd(5) == 4)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Missing data indicator used for grid type'')')
  END IF
ELSE IF(fixhd(9) == imdi.AND.fixhd(5) == 5)THEN
  IF ( printstatus >= prstatus_oper ) THEN
    WRITE(6,'('' Missing data indicator used for grid type'')')
  END IF
ELSE
  IF ( printstatus >= prstatus_min ) THEN
    WRITE(6,'(A,I9)')                                             &
    ' *** FATAL ERROR *** Invalid grid type: FIXHD(9) =',fixhd(9)
  END IF
  icode=4
  cmessage='PR_FIXHD: Consistency check'
  IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
  RETURN
END IF

IF ( printstatus >= prstatus_oper ) THEN

  IF ( fixhd(5) == 5 ) THEN    !  Boundary dataset
WRITE(6,'(A)')                                                    & 
'                       Year  Month Day Hour Min  Sec  DayNo  '
    WRITE(6,'('' First Validity time ='',6I5,I7)')(fixhd(i),i=21,27)
    WRITE(6,'('' Last  Validity time ='',6I5,I7)')(fixhd(i),i=28,34)
    WRITE(6,'('' Interval            ='',6I5,I7)')(fixhd(i),i=35,41)

  ELSE
WRITE(6,'(A)')                                                    & 
'                 Year  Month Day Hour Min  Sec  DayNo  '
    WRITE(6,'('' Data time     ='',6I5,I7)')(fixhd(i),i=21,27)
    WRITE(6,'('' Validity time ='',6I5,I7)')(fixhd(i),i=28,34)
    WRITE(6,'('' Creation time ='',6I5,I7)')(fixhd(i),i=35,41)

  END IF

  WRITE(6,'(A)') &
      '                        Start     1st dim    2nd dim   1st parm    2nd parm'
  WRITE(6,'('' Integer Consts   '',2I11,11X,I11)')fixhd(100),     &
  fixhd(101),len_inthd
  WRITE(6,'('' Real Consts      '',2I11,11X,I11)')fixhd(105),     &
  fixhd(106),len_realhd
  WRITE(6,'('' Level Dep Consts '',5I11)')fixhd(110),             &
  fixhd(111),fixhd(112),len1_levdepc,len2_levdepc
  WRITE(6,'('' Row Dep Consts   '',5I11)')fixhd(115),             &
  fixhd(116),fixhd(117),len1_rowdepc,len2_rowdepc
  WRITE(6,'('' Column Dep Consts'',5I11)')fixhd(120),             &
  fixhd(121),fixhd(122),len1_coldepc,len2_coldepc
  WRITE(6,'('' Fields of Consts '',5I11)')fixhd(125),             &
  fixhd(126),fixhd(127),len1_flddepc,len2_flddepc
  WRITE(6,'('' Extra Consts     '',2I11,11X,I11)')fixhd(130),     &
  fixhd(131),len_extcnst
  WRITE(6,'('' History Block    '',2I11,11X,I11)')fixhd(135),     &
  fixhd(136),len_dumphist
  WRITE(6,'('' CFI No 1         '',2I11,11X,I11)')fixhd(140),     &
  fixhd(141),len_cfi1
  WRITE(6,'('' CFI No 2         '',2I11,11X,I11)')fixhd(142),     &
  fixhd(143),len_cfi2
  WRITE(6,'('' CFI No 3         '',2I11,11X,I11)')fixhd(144),     &
  fixhd(145),len_cfi3
  WRITE(6,'('' Lookup Tables    '',5I11)')fixhd(150),             &
  fixhd(151),fixhd(152),len1_lookup,len2_lookup
  WRITE(6,'('' Model Data       '',2I11,11X,I11)')fixhd(160),     &
  fixhd(161),len_data
END IF

! Check model parameters against header record entries

IF (fixhd(101) >  0) THEN
  IF (len_inthd /= fixhd(101)) THEN
    IF ( printstatus >= prstatus_min ) THEN
      WRITE(6,'('' *ERROR* Integer Consts'')')
      WRITE(6,'('' Parameter and header values dont match'')')
    END IF
    icode=4
    cmessage='PR_FIXHD: Consistency check'
    IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF (fixhd(105) >  0) THEN
  IF (len_realhd /= fixhd(106)) THEN
    IF ( printstatus >= prstatus_min ) THEN
      WRITE(6,'('' *ERROR* Real Consts'')')
      WRITE(6,'('' Parameter and header values dont match'')')
    END IF
    icode=4
    cmessage='PR_FIXHD: Consistency check'
    IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
    RETURN
  END IF
END IF

IF (fixhd(110) >  0) THEN
  IF (len1_levdepc /= 0) THEN
    IF (len1_levdepc/=fixhd(111).OR.len2_levdepc/=fixhd(112)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* Level Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(115) >  0) THEN
  IF (len1_rowdepc /= 0) THEN
    IF (len1_rowdepc/=fixhd(116).OR.len2_rowdepc/=fixhd(117)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* Row Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(120) >  0) THEN
  IF (len1_coldepc /= 0) THEN
    IF (len1_coldepc/=fixhd(121).OR.len2_coldepc/=fixhd(122) )THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* Column Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(125) >  0) THEN
  IF (len1_flddepc /= 0) THEN
    IF (len1_flddepc/=fixhd(126).OR.len2_flddepc/=fixhd(127) )THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* Fields of Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(130) >  0) THEN
  IF (len_extcnst /= 0) THEN
    IF (len_extcnst /= fixhd(131)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* Extra Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(135) >  0) THEN
  IF (len_dumphist /= 0) THEN
    IF (len_dumphist /= fixhd(136)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* History File'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(140) >  0) THEN
  IF (len_cfi1 /= 0) THEN
    IF (len_cfi1 /= fixhd(141)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* CFI No 1'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(142) >  0) THEN
  IF (len_cfi2 /= 0) THEN
    IF (len_cfi2 /= fixhd(143)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* CFI No 2'')')
      END IF
      IF ( printstatus >= prstatus_oper ) THEN
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(144) >  0) THEN
  IF (len_cfi3 /= 0) THEN
    IF (len_cfi3 /= fixhd(145) )THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* CFI No 3'')')
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (fixhd(150) >  0) THEN
  IF (len2_lookup /= imdi) THEN
    IF (len1_lookup/=fixhd(151).OR.len2_lookup/=fixhd(152)) THEN
      IF ( printstatus >= prstatus_min ) THEN
        WRITE(6,'('' *ERROR* Lookup Table'')')
      END IF
      IF ( printstatus >= prstatus_oper ) THEN
        WRITE(6,'('' Parameter and header values dont match'')')
      END IF
      icode=4
      cmessage='PR_FIXHD: Consistency check'
      IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
      RETURN
    END IF
  END IF
END IF

IF (lhook) CALL dr_hook('PR_FIXHD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_fixhd
