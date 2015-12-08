! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Calculates addressing of the output dump

Module Rcf_Address_Mod

!  Subroutine Rcf_Address  - Loops over items and calls Rcf_Primary
!  Subroutine Rcf_Primary  - Calculated data lengths and addresses
!  Function Rcf_Disct_Lev  - Tests if level is discrete/continuous
!  Subroutine Rcf_PSLevCod - Decodes pseudo-level sizes
!
! Description:
!   Calculates the lengths, levels and addresses for the output
!   dump.
!
! Method:
!   Rcf_Address loops over STASHmaster records, and for those which
!   may be primary calls Rcf_Primary. This determines whether or not
!   the field is requied (via tstmsk) and all the related sizes and
!   addressing.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.

! Some variables (previously in Cstash) that are used only
! for addressing purposes (mostly in Address and Primary)
! They are all *current* values from a STASHmaster file.

Use Rcf_Ppx_Info_Mod, Only :    &
    STM_Record_Type,            &
    Ppxptr,                     &
    Ppxref_Sections,            &
    Ppxref_Items,               &
    STM_OptCodeLen

Use Rcf_Address_Vars_Mod, Only : &
    ISPACE,                & ! Space code
    IGP,                   & ! Grid of data code
    ILEV,                  & ! Level type code
    IBOT,                  & ! First level code
    ITOP,                  & ! Last level code
    IFLAG,                 & ! Level compression flag
    IOPN,                  & ! Sectional option code
    VMSK,                  & ! Integer equiv of bin vers mask
    IPSEUDO,               & ! Pseudo dimension type
    IPFIRST,               & ! First pseudo dim code
    IPLAST,                & ! Last pseudo dim code
    HALO                     ! Halo type code

Contains

!-------------------------------------------------------------------
! This routine determines the Addressing of the output dump for
! reconfiguration - includes calls to primary etc to determine
! which fields are included.
!--------------------------------------------------------------------

SUBROUTINE Rcf_Address(input_lookup)

Use Rcf_Exppx_Mod, Only :      &
    Rcf_Exppx

Use Submodel_Mod, Only :       &
    N_Internal_Model,          &
    Submodel_Ident,            &
    Submodel_for_IM,           &
    Internal_Model_List,       &
    N_Submodel_Partition_Max

Use Rcf_Model_Mod, Only :      &
    Lextra

Use Rcf_NRecon_Mod             ! Use all of this...

Use Rcf_Recon_Mod, Only :      &
    Var_Recon

USE ereport_mod, ONLY : ereport
USE lookup_addresses, ONLY : item_code
IMPLICIT NONE

! Arguments
INTEGER, OPTIONAL, INTENT(IN) :: input_lookup(:,:)

! Local scalars:
Integer    :: Im_ident  !Internal model identifier (absolute - CSMID)
Integer    :: Im_index  !Internal model index (expt. dependent)
Integer    :: Sm_ident  !Submodel identifier (absolute)
Integer    :: ISEC
Integer    :: IITM
Integer    :: RLEVS
Integer    :: RADDRESS
Integer    :: PIrow
Integer    :: I
Integer    :: ErrorStatus
Integer    :: item_code_val
Logical    :: LADDR
Logical    :: LMASK
Logical    :: input_only

Character (Len=80) :: Cmessage
! Local arrays:
!    Submodel definitions array: stores list of Im_index's
!     for each submodel partition
Integer SM_def(N_SUBMODEL_PARTITION_MAX,N_INTERNAL_MODEL_MAX)

! Local structures
Type (STM_Record_Type) :: STM_Rec

! 1.  Set STASHIN addresses and input lengths for primary fields
!   The address loop for primary fields is performed for each
!   internal model in turn. Hence, each internal model's primary
!   data occupies a contiguous block in D1. The order of these blocks
!   is the same as the order of the internal models given in the
!   array INTERNAL_MODEL_LIST.
!   User-defined prognostics are included in this primary addressing
!   routine, since they are incorporated into the ppxref lookup
!   arrays PPXI, PPXC in routine GETPPX.

!   Initialisation

PrimDataLen(:)   = 0
SM_def( :, : )   = 0

IF (PRESENT(input_lookup)) THEN
  input_only = .TRUE.
ELSE
  input_only = .FALSE.
END IF

!   Obtain submodel definitions and store in SMdef array
Do Im_index = 1,N_INTERNAL_MODEL
   !   Submodel ident.
   Sm_ident =   SUBMODEL_FOR_IM(Im_index)
   !   Internal model index
   SM_def(Sm_ident,Im_index) = Im_index
End Do

!   Primary address loop

!     Loop over submodel partitions
Do  Sm_ident = 1,N_SUBMODEL_PARTITION_MAX
  !       Initialise LEXTRA
  LEXTRA(Sm_ident)=0

  !     Initialise address for reconfiguration
  RADDRESS = 1

  !      Loop over internal models for each SM partition
  Do Im_index = 1,N_INTERNAL_MODEL
    
    !       Test whether current SM contains this IM
    If (SM_def(Sm_ident,Im_index) > 0) Then

      !        Obtain internal model identifier
      Im_ident   = INTERNAL_MODEL_LIST(Im_index)

      If(Im_ident == SUBMODEL_IDENT) Then
        PIrow  = 0
        ! Now that reconfiguration can take place over any section for VAR loop
        ! over all sections
        Do ISEC = 0,ppxref_sections
          If (ISEC == 0  .OR.    & ! Prognostic variables
              ISEC == 33 .OR.    & ! Free tracers
              ISEC == 34 .OR.    & ! UKCA tracers
              Var_Recon  .OR.    & ! VAR reconfiguration
              input_only) Then     ! Input fields only.
            !       Loop over wanted section items
            recondat_node => RecondatList( im_index, isec ) 
            Do IITM   = 1,PPXREF_ITEMS
              !   Check whether there is a primary field corresponding
              !   to this item number
              If ( PPXPTR( Im_ident, ISEC, IITM) /= 0) Then
                STM_Rec = Rcf_Exppx( Im_ident, Isec, Iitm )
                VMSK    = STM_Rec % version_mask
                ISPACE  = STM_Rec % space_code
                IGP     = STM_Rec % grid_type
                ILEV    = STM_Rec % lv_code
                IBOT    = STM_Rec % lb_code
                ITOP    = STM_Rec % lt_code
                HALO    = STM_Rec % halo_type
                Do I = 1, STM_OptCodeLen / 5
                  IOPN(I) = STM_Rec % opt_code(I)
                End Do
                IFLAG   = STM_Rec % lev_flag
                IPSEUDO = STM_Rec % pt_code
                IPFIRST = STM_Rec % pf_code
                IPLAST  = STM_Rec % pl_code

                IF (input_only) THEN
                  DO i = LBOUND(input_lookup,2), UBOUND(input_lookup,2)
                    item_code_val = input_lookup(item_code,i)
                    IF (isec == item_code_val/1000 .AND.  &
                        iitm == MOD(item_code_val,1000)) THEN
                      CALL rcf_primary(isec,iitm,im_index,im_ident, &
                                       sm_ident,raddress,pirow)
                      ! Now exit loop.
                      EXIT
                    END IF
                  END DO
                ELSE
                  If ( (ISPACE == 2) .OR. (ISPACE == 3) .OR. (ISPACE == 9) &
                      .OR.(ISPACE == 5)  .OR.(ISPACE == 10) &
                      .OR.(ISPACE == 8) .OR. Var_Recon) Then !Primary variable
! VAR reconfiguration in 4D-VAR needs to include any section in the
! reconfiguation if specified by option code. 
! Find out whether the primary is included for this version
! DEPENDS ON: tstmsk
                    CALL TSTMSK(Im_ident,ISEC,LMASK,LADDR,ErrorStatus,CMESSAGE)

                    IF (laddr) THEN
                      CALL Rcf_Primary( ISEC, IITM, Im_index, Im_ident, &
                                        Sm_ident, RADDRESS, PIrow )
                    END IF
                  END IF
                End If
              End If    !  PPXPTR(m,s,i) .ne. 0
            End Do    !  Loop over items
            recondat_node => Null()
          End If ! ISEC or VAR
        End Do    !  Loop over 'primary' sections
      End If
    End If    !  test whether SM contains IM
  End Do     !  Loop over Im_index
End Do      !  Loop over SM partitions

Return
End Subroutine Rcf_Address

!- End of subroutine code -------------------------------------------

!===================================================================

!+Compute data lengths and addresses for primary fields
! Subroutine Interface:

SUBROUTINE Rcf_Primary( ISEC, IITM, Im_index, Im_ident, &
                        Sm_ident, RADDRESS, PIrow )

Use Rcf_NRecon_Mod, Only :  &
    Recondat_Node,           &
    RecondatList,           &
    DumpProgLevs,       &
    PrimDataLen

Use Rcf_Recon_Mod, Only : &
    var_recon,            &
    l_interp_input_only

USE nlsizes_namelist_mod, ONLY : &
    ntiles,             &
    nice,               &
    tr_vars

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Address_Length_Mod, Only : &
    Rcf_Address_Length

Use Ereport_Mod, Only : &
    Ereport
IMPLICIT NONE

! Subroutine arguments:
!   Scalar arguments with intent(in):
Integer, Intent(In)    :: ISEC     ! Current section number
Integer, Intent(In)    :: IITM     ! Current section 0 item number
Integer, Intent(In)    :: Im_ident ! Current internal model number
Integer, Intent(In)    :: Im_index ! pos'n in internal model list
Integer, Intent(In)    :: Sm_ident ! Submodel identifier (absolute)
Integer, Intent(InOut) :: RADDRESS   ! Address for reconfiguration
Integer, Intent(InOut) :: PIrow    ! Counter for ProgItems array

! Local variables
Logical            :: MODEL_LEV
Integer            :: RLEVS      ! No. of levels for reconfiguration
Integer            :: RPLEVS     ! no. of pseudo-levels
Integer            :: IL1,IL2
Integer            :: IPL1,IPL2
Integer            :: LEN        ! Data length for primary field
Integer            :: ErrorStatus
Character (Len=80) :: Cmessage

!- End of Header ---------------------------------------------------

  If (ISPACE /= 10) Then

    ! Find address length per level
    CALL Rcf_Address_Length( IGP, HALO, LEN )

! If tracer, multiply length by number of tracer vars for
! LBC calculation
    IF ( ( IGP .EQ. 26 ) .AND. ( ITOP .EQ. 11 ) ) THEN
      LEN = LEN * TR_VARS
    END IF

! If level type is 0, is LBC so multiply size by levels
    IF ( ILEV .EQ. 0 ) THEN
      CALL Rcf_Level_Code( IBOT, IL1, Output_Grid )
      CALL Rcf_Level_Code( ITOP, IL2, Output_Grid )
      LEN = LEN * (IL2 - IL1 + 1)
    END IF

    MODEL_LEV = Rcf_Disct_Lev( ILEV )
    If (MODEL_LEV .OR.(ILEV == 5 .AND. IPSEUDO /= 0)) Then
      ! Field has model levels - decode level codes
      If (ILEV  /=  5) Then
        CALL Rcf_Level_Code( IBOT, IL1, Output_Grid )
        CALL Rcf_Level_Code( ITOP, IL2, Output_Grid )
      Else
        IL1=1
        IL2=1
      End If

      ! No. of model levels (for reconfiguration)
      RLEVS=IL2-IL1+1
      ! Initialise first & last pseudo level indices
      IPL1 =0
      IPL2 =0

      If (Iflag == 0 .AND. Ipseudo /= 0) Then
        ! Primary with input on all available pseudo levels -
        !   decode pseudo level codes
        CALL Rcf_PSLevCod( IPFIRST, NTILES, NICE, IPL1, 'F' )
        CALL Rcf_PSLevCod( IPLAST , NTILES, NICE, IPL2, 'L' )
      End If

      RPLEVS=IPL2-IPL1+1
      ! Multiply length per level by no. of levels
      LEN=LEN*(IL2-IL1+1)*(IPL2-IPL1+1)

      If (ISPACE /= 4 .AND. ISPACE /= 9) Then
        ! Increment no. of headers
        DumpProgLevs   (Im_ident)=   DumpProgLevs(Im_ident) &
                                      +(IL2-IL1+1)*(IPL2-IPL1+1)
      End If

    Else  ! Not model levels
      RLEVS=1
      RPLEVS=1
      If (ISPACE /= 4 .AND. ISPACE /= 9) Then
        DumpProgLevs   (Im_ident)=DumpProgLevs   (Im_ident)+1
      End If
    End If

    ! Addresses are set relative to the beginning of the primary data,
    !  since the primary data starts at the beginning of D1.
    If (ISPACE /= 9) Then
      ! Start address for this primary field
      ! Increment LPRIM by LEN (=data length for this primary field)
      PrimDataLen(Im_ident)      =PrimDataLen(Im_ident)+LEN
    EndIf

    ! Store levels, lengths and addresses required for reconfiguration
    !                                                in array Recondat
    If (ISPACE /= 4 .AND. ISPACE /= 9) Then
      If (ISEC == 0 .OR. ISEC == 33 .OR. ISEC == 34 .OR. &
          VAR_RECON .OR. l_interp_input_only) Then
        Allocate(recondat_node % recondat_info)
        recondat_node % recondat_info % sec_item=IITM + 1000*ISEC
        recondat_node % recondat_info % rlevs=RLEVS
        recondat_node % recondat_info % len=LEN
        recondat_node % recondat_info % raddress=RADDRESS
        recondat_node % recondat_info % rplevs=RPLEVS
        
        Allocate(recondat_node % next)
        recondat_node => recondat_node % next
      Else
        ErrorStatus=1
        Cmessage='Rcf_Primary : Invalid value for ISEC, section no.'
        Call Ereport( 'Rcf_Primary', ErrorStatus, Cmessage )
      End If
      RADDRESS = RADDRESS+LEN
    End If
  End If  ! ISPACE .ne. 10

Return

End Subroutine Rcf_Primary

!====================================================================

!+Test whether level type is discrete (model) or continuous (non-model)
! Function Interface:
Logical FUNCTION Rcf_Disct_Lev( LEV_CODE )

Use Ereport_Mod, Only : &
    Ereport

IMPLICIT NONE

! Function arguments:
Integer, Intent(In)          :: LEV_CODE !Level code from STASHmaster

! ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'Disct_Lev'
Integer                      :: ErrorStatus
Character (Len=80)           :: CMESSAGE

!- End of Header ----------------------------------------------

If (LEV_CODE == 1 .OR. LEV_CODE == 2 .OR. LEV_CODE == 6 .OR. &
                       LEV_CODE == 10) Then
  Rcf_Disct_Lev=.TRUE.
Else If (LEV_CODE .GE. 0 .AND. LEV_CODE .LE. 10) Then
  Rcf_Disct_Lev=.FALSE.
Else
  Rcf_Disct_Lev=.FALSE.
  ErrorStatus=1
  CMESSAGE='Rcf_DISCT_LEV : Invalid level type in STASHmaster'

  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

End Function Rcf_Disct_Lev
!- End of Function code --------------------------------------------

!===================================================================

!+Decode the STASH pseudo level code
! Subroutine Interface:
SUBROUTINE Rcf_PSLevCod( ILIN, NTILES, NICE, ILOUT, SWTCH )

Use Ereport_Mod, Only : &
    Ereport

USE rad_input_mod, ONLY : &
    h_swbands,            &
    h_lwbands

! JULES
USE ancil_info, ONLY: nsmax

  USE nstypes
  USE clmchfcg_scenario_mod, ONLY: nsulpat

IMPLICIT NONE
! Description:
!   Sets ILOUT to an appropriate pseudo level size according
!   to the value of IL
!
! Subroutine arguments:
Integer, Intent(In)           :: ILIN    ! Model pseudo level code
Integer, Intent(In)           :: NTILES  ! Number of surface tiles
Integer, Intent(In)           :: NICE    ! Number of seaice catagories
Character (Len=1), Intent(In) :: SWTCH
Integer, Intent(Out)          :: ILOUT   ! An actual pseudo level

! Local variables
Integer              :: ErrorStatus
Character (Len=80)   :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_PSLevCod'


!- End of Header --------------------------------------------------

If (SWTCH == 'F') Then
  If (ILIN == 1) Then
    ILOUT=1
! Ocean assimilation groups - removed
  Else
    Write(6,*) &
     &   'MSG FROM RCF_PSLEVCOD: ', &
     &   'INAPPROPRIATE FIRST PSEUDO LEVEL CODE FOUND ',ILIN
    ErrorStatus=2
    Cmessage = 'Error from PSLevCod!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

  End If
Else If (SWTCH == 'L') Then
  If (ILIN == 1) Then
    ILOUT=H_SWBANDS
  Else If (ILIN == 2) Then
    ILOUT=H_LWBANDS
  Else If (ILIN == 4) Then
    ! Last frequency (wave model)
! Wave model not used in RCF
!    ILOUT=NFRE
     ILOUT=0
  Else If (ILIN == 5) Then
    ! Last wave train (wave model)
!   ILOUT=NWTRAIN
    ILOUT=0
  Else If ( ILIN  ==  6 ) Then
    ! Last index for HadCM2 sulphate loading patterns.
    ILOUT = NSULPAT
  Else If ( ILIN  ==  7 ) Then
    ! All surface types
    ILOUT = NTYPE
  Else If ( ILIN  ==  8 ) Then
    ! Plant functional types only
    ILOUT = NPFT
  Else If ( ILIN  ==  9 ) Then
    ! All tiles
    ILOUT = NTILES
  Else If ( ILIN  ==  10 ) Then
    ! All seaice catagories
    ILOUT = NICE
  Else If ( ILIN  ==  11 ) Then
    !New snow scheme: all tiles <times> max no. of snow levels
    ILOUT = NTILES*NSMAX
  Else
    Write(6,*) &
     &   'MSG FROM RCF_PSLEVCOD: ', &
     &   'INAPPROPRIATE LAST PSEUDO LEVEL CODE FOUND ',ILIN
    ErrorStatus=2
    Cmessage = 'Error from PSLevCod'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

  End If

End If

Return

End Subroutine Rcf_PSLevCod

End Module Rcf_Address_Mod
