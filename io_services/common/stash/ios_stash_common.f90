! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Description of the contents of an async stash control buffer:
!

! Each item packed into the dispatch buffer is accompanied by a control record
! of 64 bit integer data that describes that payload (and what to do with it). 
! The client routine will insert the  items automatically, so
! that the buffer looks like (after N packs):
!
! <auto1>[user1]<auto2>[user2]...<autoN>[userN]
!
! <auto> consists of <ios_stash_control_auto_len> words, the length 
! of user is specified in the <auto> part, but at most, 
! ios_async_control_sz_max words. 
!
! loc_* variables describe the location in one of these metadata 
! blocks of a particular item. Occasionally for items with only a small 
! possible number range, we will use a single 64 bit integer to host
! up to 4 variables using the simplistic pack4/unpack4 routines
!
! Please see UMDP C10 for further detail.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE IOS_Stash_Common
  
  USE IOS_Common
  USE um_types,  ONLY: integer64
  IMPLICIT NONE
 

! Size of the <auto> record 



  INTEGER, PARAMETER  :: ios_stash_control_auto_len  = 2


! Location of fields in the 
! ***** <auto> *****
! record for various items
  INTEGER, PARAMETER  :: loc_record_start            = 1
  INTEGER, PARAMETER  :: loc_record_len_control      = 2

! Only used in debug build
  INTEGER, PARAMETER  :: loc_record_handle           = 3
  INTEGER, PARAMETER  :: loc_record_record           = 4
  INTEGER, PARAMETER  :: loc_record_len_data         = 5
  INTEGER, PARAMETER  :: loc_record_slot             = 6

! Maximum size of the [user] record 
  INTEGER, PARAMETER  :: ios_async_control_sz_max    = 14

! Location of fields in the 
! ***** [user] *****
! record for various items
  INTEGER, PARAMETER  :: loc_record_type             = 1
  INTEGER, PARAMETER  :: loc_fld_type                = 2
  INTEGER, PARAMETER  :: loc_preprocess_flag         = 3
  INTEGER, PARAMETER  :: loc_subdomain_flag          = 4
  INTEGER, PARAMETER  :: loc_packing_flag            = 5
  INTEGER, PARAMETER  :: loc_boundary                = 6
  INTEGER, PARAMETER  :: loc_pack_type               = 7

! Used by stash only
  INTEGER, PARAMETER  :: loc_comp_accry              = 8
  INTEGER, PARAMETER  :: loc_dmi                     = 9

! Used differently by stash/dumps
  INTEGER, PARAMETER  :: loc_disk_block              = 10

! Used by dumps only
  INTEGER, PARAMETER  :: loc_disk_block_start        = 11

  INTEGER, PARAMETER  :: loc_seek_target             = 12
  INTEGER, PARAMETER  :: loc_landmaskcompress        = 13
  INTEGER, PARAMETER  :: loc_id                      = 14

! Maximum total per field for metadata
  INTEGER, PARAMETER  :: ios_max_cntrlwords_per_field &
       = ios_async_control_sz_max + ios_stash_control_auto_len


! Parameterised possible values of metadata variables. 
  INTEGER, PARAMETER  :: ios_stash_record_start      = -99999

  INTEGER, PARAMETER  :: ios_stash_distributed_field = 12
  INTEGER, PARAMETER  :: ios_repeat_record           = 385
  INTEGER, PARAMETER  :: ios_stash_data              = 123

  INTEGER, PARAMETER  :: ios_stash_nopreprocess      = 434
  INTEGER, PARAMETER  :: ios_stash_preprocess        = 777

  INTEGER, PARAMETER  :: ios_partial_field           = 221
  INTEGER, PARAMETER  :: ios_full_field              = 121

  INTEGER, PARAMETER  :: ios_packing                 = 7
  INTEGER, PARAMETER  :: ios_no_packing              = 21


  INTEGER, PARAMETER  :: ios_packing_type_wgdos      = 734
  INTEGER, PARAMETER  :: ios_packing_type_pack21     = 437

  INTEGER, PARAMETER  :: ios_packing_type_landmask   = 87
  INTEGER, PARAMETER  :: ios_packing_type_nolandmask = 98

  INTEGER, PARAMETER  :: ios_unset_int               = 0
  REAL   , PARAMETER  :: ios_unset_real              = 0.0

! tag space for metadata
  INTEGER, PARAMETER  :: ios_async_control_tag_base  = 1000000
  INTEGER, PARAMETER  :: ios_async_control_tag_incr  = 1

  INTEGER             :: mpl_threading
  INTEGER             :: ios_async_comm
  INTEGER             :: ios_local_comm

  LOGICAL, SAVE       :: async_stash_on              = .FALSE.
  LOGICAL, SAVE       :: async_dump_on               = .FALSE.

! Single constant for packing/unpacking data:
  INTEGER(KIND=integer64), PARAMETER, PRIVATE  :: imax = 10000

  CHARACTER (LEN=132), PRIVATE   :: ios_stsh_cmn_message



! Options from IOSCNTL:
!
! Defaults cope with a large model of N512 or UKV
! type sizes whilst keeping memory use sane
! consequently this will not give good performance
  INTEGER             :: IOS_AsyncNumSlots           = 10
  INTEGER             :: IOS_AsyncMaxFieldsInPack    = 76
  LOGICAL             :: IOS_AsyncSendNull           = .FALSE.
  LOGICAL             :: IOS_AsyncDoStats            = .TRUE.

  TYPE IOS_async_buffer
    SEQUENCE
    INTEGER              :: tag
    INTEGER              :: control_offset
    INTEGER              :: buffer_offset
    INTEGER              :: num_records_in_pack
    INTEGER              :: state
    INTEGER              :: unit
    INTEGER              :: handle 
    INTEGER              :: bufferType
    REAL                 :: time_create
    REAL                 :: time_update
    REAL                 :: time_complete






    INTEGER(KIND=integer64), POINTER     :: control(:)
    REAL, POINTER        :: buffer(:)

  END TYPE IOS_async_buffer





  TYPE(IOS_async_buffer), &
      POINTER            :: slot(:)
  INTEGER, POINTER       :: astash_mpi_requests(:)


CONTAINS

  LOGICAL FUNCTION isUsingAsyncStash() RESULT(r)
    IMPLICIT NONE
    r=async_stash_on
  END FUNCTION isUsingAsyncStash

  LOGICAL FUNCTION isUsingAsyncDumps() RESULT(r)
    IMPLICIT NONE
    r=async_dump_on
  END FUNCTION isUsingAsyncDumps

  SUBROUTINE useAsyncStash(active)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: active
    async_stash_on=active
  END SUBROUTINE useAsyncStash

  SUBROUTINE useAsyncDump(active)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: active
    async_dump_on=active
  END SUBROUTINE useAsyncDump

 SUBROUTINE IOS_pack4(word,a,b,c,d)
     USE IOS_Common
     IMPLICIT NONE
     INTEGER, INTENT(IN)                  :: a
     INTEGER, INTENT(IN)                  :: b
     INTEGER, INTENT(IN)                  :: c
     INTEGER, INTENT(IN)                  :: d
     INTEGER(KIND=integer64), INTENT(OUT) :: word
  
 
     IF (a>=imax .OR. &
         b>=imax .OR. &
         c>=imax .OR. &
         d>=imax .OR. &
         a<0     .OR. &
         b<0     .OR. &
         c<0     .OR. &
         d<0           ) THEN
 
       CALL IOS_Ereport('ios_Stash_common:pack4',99, &
           'Number out of range')
     END IF
 
     word=                      &
         d +                    &
         c * imax +             &
         b * imax * imax +      &
         a * imax * imax * imax
  
   END SUBROUTINE IOS_pack4
 
   SUBROUTINE IOS_unpack4(word,a,b,c,d)
     IMPLICIT NONE
     INTEGER(KIND=integer64), INTENT(IN)  :: word
     INTEGER, INTENT(OUT)                 :: a
     INTEGER, INTENT(OUT)                 :: b
     INTEGER, INTENT(OUT)                 :: c
     INTEGER, INTENT(OUT)                 :: d
! These (and imax) set to integer64 to prevent an overflow during 32-bit builds:
     INTEGER(KIND=integer64)              :: tmp
     INTEGER(KIND=integer64), PARAMETER   :: imax2=imax*imax
     INTEGER(KIND=integer64), PARAMETER   :: imax3=imax*imax2
 
     tmp=word
     a=tmp/(imax3)
     tmp=tmp-a*(imax3)
     b=tmp/(imax2)
     tmp=tmp-b*(imax2)
     c=tmp/(imax)
     tmp=tmp-c*(imax)
     d=tmp
  
   END SUBROUTINE IOS_unpack4

END MODULE IOS_Stash_Common

