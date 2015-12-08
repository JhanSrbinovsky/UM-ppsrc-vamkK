! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! + Gathers any type of field from many processors to one processor

Module Rcf_General_Gather_Field_Mod

!  Subroutine Rcf_General_Gather_Field - generalised gather field
!
! Description:
!   GAthers any type of field from many processors to one processor by
!   consideration of type - only copes with full fields
!
! Method:
!   Land compressed fields are uncompressed, gathered and recompressed
!   LBCs are dealt with by their own routine
!   Zonal fields have their own routine
!   All other fields are "normal" and use rcf_gather_field
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

Contains

Subroutine Rcf_General_Gather_Field( LOCAL_FIELD,  GLOBAL_FIELD , &
                                     LOCAL_SIZE,   GLOBAL_SIZE ,  &
                                     STASH_RECORD, GATHER_PE )

Use Rcf_Lsm_Mod, Only  :     &
    Local_Land_Field,        &
    Local_Atmos_Landmask,    &
    Glob_Atmos_Landmask

USE UM_ParVars, Only :                    &
    blsizep,                blsizeu,           &
    blsizev,                glsizeu,           &
    glsizev,                glsizep,           &
    blsizer,                glsizer,           &
    lasize,                 gc_all_proc_group, &
    fld_type_p,             fld_type_u,        &
    fld_type_v,             fld_type_r,        &
    mype

Use Ereport_Mod, Only : &
     &    Ereport

Use Rcf_Gather_Atmos_LBCs_Mod, Only : &
     &    Rcf_Gather_Atmos_LBCs

Use Rcf_Gather_Ocean_LBCs_Mod, Only : &
     &    Rcf_Gather_Ocean_LBCs

Use Rcf_Gather_Zonal_Field_Mod, Only : &
     &    Rcf_Gather_Zonal_Field

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field

Use Rcf_Global_To_Local_Mod, Only : &
    Rcf_Get_Fld_Type

Use Rcf_Ppx_Info_Mod, Only : &
    STM_Record_Type

USE mask_compression, ONLY: compress_to_mask, expand_from_mask

USE UM_ParParams

USE cppxref_mod

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)   :: Global_Size   ! size of global_field
Integer, Intent(In)   :: Gather_PE     ! PE to collect global data on
Integer, Intent(Out)  :: Local_Size    ! size of local field

Real, Intent(In)      :: Local_Field(*)      ! local part of field
Real, Intent(Out)     :: Global_Field( Global_Size )

Type( STM_Record_Type ), Intent(In) :: Stash_Record

! Local variables
Integer       :: dummy            ! ignored argument
Integer       :: fld_type         ! P or U or V
Integer       :: local_x          ! local x size
Integer       :: local_y          ! local y size
Integer       :: global_x         ! global x size
Integer       :: global_y         ! global y size
Integer       :: ErrorStatus      ! error code

Character (Len=*), Parameter :: RoutineName = 'Rcf_General_Gather_Field'
Character (Len=80)           :: Cmessage

Real          :: buf_expand( glsizep(1) * glsizep(2) )
Real          :: buf_expand_local(          &
    lasize(1,fld_type_p,halo_type_single) * &
    lasize(2,fld_type_p,halo_type_single) )

!===================================================================

! Choose action depending on grid code
Select Case( Stash_Record % grid_type )

!------------------------------------------------------------------
! Land compressed fields
!------------------------------------------------------------------

  Case( ppx_atm_compressed )

! Unpack the local field out to full (local) field size and
! put this into the array buf_expand_local

    CALL expand_from_mask(buf_expand_local, LOCAL_FIELD, &
                          local_atmos_landmask,          &
                          lasize(1,fld_type_p,halo_type_single)* &
                          lasize(2,fld_type_p,halo_type_single), &
                          dummy)

! Now gather in all the processors local fields into the global
! field (array buf_expand)

    Call Rcf_Gather_Field(buf_expand_local, buf_expand,          &
                          lasize(1,fld_type_p,halo_type_single), &
                          lasize(2,fld_type_p,halo_type_single), &
                          glsizep(1),        glsizep(2),         &
                          GATHER_PE,        GC_ALL_PROC_GROUP  )

! And now pack the global field (buf_expand) back to land points
! and put into the array GLOBAL_FIELD.

    IF (mype .EQ. 0) THEN
      CALL compress_to_mask(buf_expand, GLOBAL_FIELD, &
                          glob_atmos_landmask,      &
                          glsizep(1)*glsizep(2),dummy)
    ENDIF

    Local_Size = local_land_field

!-------------------------------------------------------------------
! Atmosphere Lateral boundary fields
!-------------------------------------------------------------------
  Case( ppx_atm_lbc_theta, ppx_atm_lbc_u, ppx_atm_lbc_v )

    Call Rcf_Gather_Atmos_LBCs( Global_field, global_size, &
                                Local_field,  Local_size,  &
                                Stash_Record, Gather_PE )


!-------------------------------------------------------------------
! Ocean Lateral boundary fields
!-------------------------------------------------------------------
  Case ( ppx_ocn_rim )

    Call Rcf_Gather_Ocean_LBCs( Global_field, Global_size, &
                                Local_field,  Local_size,  &
                                Stash_Record, Gather_PE )

!-------------------------------------------------------------------
! Zonal fields
!-------------------------------------------------------------------
  Case ( ppx_atm_tzonal, ppx_atm_uzonal, &
         ppx_ocn_tzonal, ppx_ocn_uzonal )

    If ( Stash_Record % grid_type == ppx_atm_tzonal .OR. &
         Stash_Record % grid_type == ppx_ocn_tzonal ) Then   ! P grid
      global_y = glsizep(2)
      local_y  = blsizep(2)

    Else
      global_y = glsizeu(2)
      local_y  = blsizeu(2)
    End If

    Call Rcf_Gather_Zonal_Field( Local_Field, Global_Field,            &
                                 local_y,     global_y,                &
                                 1,           Stash_Record % grid_type,&
                                 Gather_PE )

    Local_Size = local_y
!-------------------------------------------------------------------
! "Normal" fields
!-------------------------------------------------------------------
  Case &
!     atmosphere grids
       ( ppx_atm_tall,    &! Atmos T points
         ppx_atm_tland,   &! Atmos T land points
         ppx_atm_tsea,    &! Atmos T sea points
         ppx_atm_uall,    &! Atmos U points
         ppx_atm_uland,   &! Atmos U land points
         ppx_atm_usea,    &! Atmos U sea points
         ppx_atm_cuall,   &! Atmos C grid U pts
         ppx_atm_cvall,   &! Atmos C grid V pts
         ppx_atm_ozone,   &! Atmos ozone field
         ppx_atm_river,   &! Atmos river-routing field
!     ocean grids
         ppx_ocn_tcomp,   &! Ocean "Compressed" T
         ppx_ocn_ucomp,   &! Ocean "Compressed" u
         ppx_ocn_tall,    &! Ocean T points (cyc)
         ppx_ocn_uall,    &! Ocean U points (cyc)
         ppx_ocn_cuall,   &! Ocean C grid U pts
         ppx_ocn_cvall,   &! Ocean C grid V pts
         ppx_ocn_tfield,  &! Ocean T points
         ppx_ocn_ufield) ! Ocean U points

    fld_type = Rcf_Get_Fld_Type(Stash_Record % grid_type)

    If (fld_type == fld_type_p) Then
      global_x = glsizep(1)
      global_y = glsizep(2)
      local_x  = blsizep(1)
      local_y  = blsizep(2)
    Else If (fld_type == fld_type_u) Then
      global_x = glsizeu(1)
      global_y = glsizeu(2)
      local_x  = blsizeu(1)
      local_y  = blsizeu(2)
    Else If (fld_type == fld_type_v) Then
      global_x = glsizev(1)
      global_y = glsizev(2)
      local_x  = blsizev(1)
      local_y  = blsizev(2)
    Else If (fld_type == fld_type_r) Then
      global_x = glsizer(1)
      global_y = glsizer(2)
      local_x  = blsizer(1)
      local_y  = blsizer(2)
    End If

    Call Rcf_Gather_Field( Local_Field, Global_Field,      &
                           local_x,     local_y,           &
                           global_x,    global_y,          &
                           Gather_PE,   GC_ALL_PROC_GROUP )

    Local_Size = local_x * local_y

!-------------------------------------------------------------------
! Any other type of field
!-------------------------------------------------------------------
   Case Default

    ErrorStatus = 10
    Cmessage = 'Field type not recognised for Gather'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select


Return
End Subroutine Rcf_General_Gather_Field
End Module Rcf_General_Gather_Field_Mod



