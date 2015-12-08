! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Misc

! Purpose: Holds PARAMETER definitions to describe the structure of
!          each STASHmaster file record plus some valid entries.      

MODULE cppxref_mod

IMPLICIT NONE

! Primary file record definition
      ! length of ID in a record
      INTEGER, PARAMETER :: ppxref_idlen = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      INTEGER, PARAMETER :: ppxref_charlen = 36

      ! number of packing profiles
      INTEGER, PARAMETER :: ppxref_pack_profs = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      INTEGER, PARAMETER :: ppxref_codelen = 33 + ppxref_pack_profs

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      INTEGER, PARAMETER :: ppx_charword = ((ppxref_charlen+3)/4)

      ! read buffer record length
      INTEGER, PARAMETER :: ppx_recordlen = ppx_charword+ppxref_codelen
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      INTEGER, PARAMETER ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      INTEGER, PARAMETER ::  ppx_section_number = 2  ! Section number
                                                     ! address
      INTEGER, PARAMETER ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      INTEGER, PARAMETER ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      INTEGER, PARAMETER ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      INTEGER, PARAMETER ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      INTEGER, PARAMETER ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      INTEGER, PARAMETER ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      INTEGER, PARAMETER ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      INTEGER, PARAMETER ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      INTEGER, PARAMETER ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      INTEGER, PARAMETER ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      INTEGER, PARAMETER ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_user_code      =26  ! User code address
      INTEGER, PARAMETER ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_cf_levelcode   =27
      INTEGER, PARAMETER ::  ppx_cf_fieldcode   =28
      INTEGER, PARAMETER ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      INTEGER, PARAMETER ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      INTEGER, PARAMETER ::  ppx_halo_type      =33
      INTEGER, PARAMETER ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      INTEGER, PARAMETER ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      INTEGER, PARAMETER :: ppx_atm_tall=1        ! All T points (atmos)
      INTEGER, PARAMETER :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      INTEGER, PARAMETER :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      INTEGER, PARAMETER :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      INTEGER, PARAMETER :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      INTEGER, PARAMETER :: ppx_atm_uall=11       ! All u points (atmos)
      INTEGER, PARAMETER :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      INTEGER, PARAMETER :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      INTEGER, PARAMETER :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      INTEGER, PARAMETER :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      INTEGER, PARAMETER :: ppx_atm_scalar=17     ! Scalar (atmos)
      INTEGER, PARAMETER :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      INTEGER, PARAMETER :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      INTEGER, PARAMETER :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      INTEGER, PARAMETER :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      INTEGER, PARAMETER :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      INTEGER, PARAMETER :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      INTEGER, PARAMETER :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      INTEGER, PARAMETER :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      INTEGER, PARAMETER :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      INTEGER, PARAMETER :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      INTEGER, PARAMETER :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      INTEGER, PARAMETER :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      INTEGER, PARAMETER :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      INTEGER, PARAMETER :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      INTEGER, PARAMETER :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      INTEGER, PARAMETER :: ppx_ocn_scalar=47     ! Scalar (ocean)
      INTEGER, PARAMETER :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      INTEGER, PARAMETER :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      INTEGER, PARAMETER :: ppx_ocn_lbc_u=53      ! on T & U grids
      INTEGER, PARAMETER :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      INTEGER, PARAMETER :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      INTEGER, PARAMETER :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      INTEGER, PARAMETER :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_full_level=1      ! Model full level
      INTEGER, PARAMETER :: ppx_half_level=2      ! Model half level
      INTEGER, PARAMETER :: ppx_rho_level=1       ! Model rho level
      INTEGER, PARAMETER :: ppx_theta_level=2     ! Model theta level
      INTEGER, PARAMETER :: ppx_single_level=5    ! Model single level
      INTEGER, PARAMETER :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_type_real=1       ! Real data type
      INTEGER, PARAMETER :: ppx_type_int=2        ! INTEGER data type
      INTEGER, PARAMETER :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      INTEGER, PARAMETER :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      INTEGER, PARAMETER :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      INTEGER, PARAMETER :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      INTEGER, PARAMETER :: ppx_lbvc_unset   =  0 ! unset
      INTEGER, PARAMETER :: ppx_lbvc_height  =  1 ! height
      INTEGER, PARAMETER :: ppx_lbvc_depth   =  2 ! depth (ocean)
      INTEGER, PARAMETER :: ppx_lbvc_soil_level =  6 ! level (soil)
      INTEGER, PARAMETER :: ppx_lbvc_pressure=  8 ! pressure
      INTEGER, PARAMETER :: ppx_lbvc_hyb_pres=  9 ! hybrid pressure
      INTEGER, PARAMETER :: ppx_lbvc_theta   = 19 ! potential T
      INTEGER, PARAMETER :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      INTEGER, PARAMETER :: ppx_lbvc_PV      = 82 ! potential vorticity
      INTEGER, PARAMETER :: ppx_lbvc_msl     =128 ! mean sea level
      INTEGER, PARAMETER :: ppx_lbvc_surface =129 ! surface
      INTEGER, PARAMETER :: ppx_lbvc_freezing =132 ! freezing 
      INTEGER, PARAMETER :: ppx_lbvc_atmtop  =133 ! top of atmosphere
      INTEGER, PARAMETER :: ppx_lbvc_canopy_water =275 ! canopy water


END MODULE cppxref_mod
