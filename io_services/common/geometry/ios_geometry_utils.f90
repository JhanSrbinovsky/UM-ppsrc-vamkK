! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Provides utility functionality for manipulating decomposition data

MODULE IOS_Geometry_utils

  USE IOS_Model_geometry
  USE yomhook, ONLY: lhook, dr_hook  
  USE parkind1, ONLY: jprb, jpim

! A simple derived type for conveying a boundary.
  TYPE box
     INTEGER :: N
     INTEGER :: S
     INTEGER :: E
     INTEGER :: W
  END TYPE box 

  CHARACTER (LEN=132), PRIVATE :: geom_util_message

  ! Classification values
  INTEGER, PARAMETER  :: not_computed=-3
  INTEGER, PARAMETER  :: no_intersection=-2
  INTEGER, PARAMETER  :: complete_intersection=-1
  ! >0 implies the number of 'patches' from the local domain
  ! =0 implies 'not relevent'

  INTEGER, PARAMETER  :: WesterlyZone=1
  INTEGER, PARAMETER  :: EasterlyZone=2

! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS


  SUBROUTINE printBox(b)
    TYPE(box), INTENT(IN) :: b
    WRITE(6,'(A,i5,A,i5,A,i5,A,i5,A,i5,A)')          &
        'Box: S=',b%s,' N=',b%n,' W=',b%w,' E=',b%e, &
        '(data_size=',(b%n-b%s+1)*(b%e-b%w+1),')' 
  END SUBROUTINE printBox

  ! provide coordinates in the local domain on processor 'cpu', corresponding 
  ! to the bit of the local domain included in the global region
  ! defined by 'box'. This routine ASSUMES that cpu intersects the box_in
  ! so you need to call classify_subdomain first
  SUBROUTINE getLocalSubdomainBounds(cpu,fld,box_in,localSubdomain)

    IMPLICIT NONE
    TYPE(box), INTENT(IN)        :: box_in   ! A global box
    TYPE(box), INTENT(OUT)       :: localSubdomain  ! A local box
    INTEGER, INTENT(IN)          :: fld
    INTEGER, INTENT(IN)          :: cpu
    TYPE(box)                    :: globalSubdomain
    TYPE(box)                    :: processorDomain
    TYPE(box)                    :: modelDomain
    LOGICAL                      :: ew_wrap
    INTEGER                      :: East
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook                            &
        ('IOS_GEOMETRY_UTILS:GETLOCALSUBDOMAINBOUNDS', &
        zhook_in,zhook_handle)

    ! Normalize the box
    globalSubdomain=box_in
    CALL Normalize_Box(globalSubdomain,fld,ew_wrap)

    ! Set the whole domain in local coords first
    localSubdomain%S=1
    localSubdomain%N=size_map(2,fld,cpu)
    localSubdomain%W=1
    localSubdomain%E=size_map(1,fld,cpu)

    ! my domain in global coords
    processorDomain%S=offset_map(2,fld,cpu)
    processorDomain%N=offset_map(2,fld,cpu)+size_map(2,fld,cpu)-1
    processorDomain%W=offset_map(1,fld,cpu)
    processorDomain%E=offset_map(1,fld,cpu)+size_map(1,fld,cpu)-1

    modelDomain%S=1
    modelDomain%N=atm_global_points(2,fld)
    modelDomain%W=1
    modelDomain%E=atm_global_points(1,fld)
    
    ! and then pull the sides in

    ! South
    IF (globalSubdomain%S >= processorDomain%S) THEN
       IF (globalSubdomain%S <= processorDomain%N) THEN
          localSubdomain%S=globalSubdomain%S-processorDomain%S+1
       ELSE
          CALL IOS_Ereport('getLocalSubdomainBounds',99, &
              'South boundary is North of my domain')
       END IF
    END IF

    !North
    IF (globalSubdomain%N <= processorDomain%N) THEN
       IF (globalSubdomain%N >= processorDomain%S ) THEN
          localSubdomain%N=globalSubdomain%N-processorDomain%S+1
       ELSE
          CALL IOS_Ereport('getLocalSubdomainBounds',99, &
              'North boundary is South of my domain')
       END IF
    END IF


    IF (ew_wrap) THEN
      East=globalSubdomain%E-modelDomain%E
    
    !East
      IF (East >= processorDomain%W) THEN
        IF (East <= processorDomain%E) THEN
          localSubdomain%E=East-processorDomain%W+1
        END IF
      END IF
    !West
    IF (globalSubdomain%W <= processorDomain%E) THEN
      IF (globalSubdomain%W >= processorDomain%W) THEN
        localSubdomain%W=globalSubdomain%W-processorDomain%W+1
      END IF
    END IF
    ELSE
      East=globalSubdomain%E
      IF (East <= processorDomain%E) THEN
        IF (East >= processorDomain%W) THEN
          localSubdomain%E=East-processorDomain%W+1
        ELSE
          CALL printBox(box_in)
          CALL printBox(processorDomain)
          CALL IOS_Ereport('getLocalSubdomainBounds',99, &
              'East boundary is west of my domain')        
        END IF
      END IF
    !West
    IF (globalSubdomain%W >= processorDomain%W) THEN
      IF (globalSubdomain%W <= processorDomain%E) THEN
        localSubdomain%W=globalSubdomain%W-processorDomain%W+1
      ELSE
        CALL printBox(box_in)
        CALL printBox(processorDomain)
        CALL IOS_Ereport('getLocalSubdomainBounds',99,&
            'West boundary is east of my domain')
      END IF
    END IF
      
    END IF
    
    IF (lhook) CALL dr_hook('IOS_GEOMETRY_UTILS:GETLOCALSUBDOMAINBOUNDS', &
        zhook_out,zhook_handle)

  END SUBROUTINE getLocalSubdomainBounds

  SUBROUTINE getLocalSubdomainBoundsCyclic(cpu,fld,box_in,localSubdomain, &
      region)

    IMPLICIT NONE
    TYPE(box), INTENT(IN)        :: box_in   ! A global box
    TYPE(box), INTENT(OUT)       :: localSubdomain  ! A local box
    INTEGER, INTENT(IN)          :: fld
    INTEGER, INTENT(IN)          :: cpu
    INTEGER, INTENT(IN)          :: region
    TYPE(box)                    :: globalSubdomain
    TYPE(box)                    :: processorDomain
    TYPE(box)                    :: modelDomain
    LOGICAL                      :: ew_wrap
    INTEGER                      :: East
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook &
        ('IOS_GEOMETRY_UTILS:GETLOCALSUBDOMAINBOUNDSCYCLIC', &
        zhook_in,zhook_handle)

    ! Normalize the box
    globalSubdomain=box_in
    CALL Normalize_Box(globalSubdomain,fld,ew_wrap)

    ! Set the whole domain in local coords first
    localSubdomain%S=1
    localSubdomain%N=size_map(2,fld,cpu)
    localSubdomain%W=1
    localSubdomain%E=size_map(1,fld,cpu)

    ! my domain in global coords
    processorDomain%S=offset_map(2,fld,cpu)
    processorDomain%N=offset_map(2,fld,cpu)+size_map(2,fld,cpu)-1
    processorDomain%W=offset_map(1,fld,cpu)
    processorDomain%E=offset_map(1,fld,cpu)+size_map(1,fld,cpu)-1

    modelDomain%S=1
    modelDomain%N=atm_global_points(2,fld)
    modelDomain%W=1
    modelDomain%E=atm_global_points(1,fld)
    
    ! and then pull the sides in

    ! South
    IF (globalSubdomain%S >= processorDomain%S) THEN
       IF (globalSubdomain%S <= processorDomain%N) THEN
          localSubdomain%S=globalSubdomain%S-processorDomain%S+1
       ELSE
          CALL IOS_Ereport('getLocalSubdomainBoundsCyclic',99, &
              'South boundary is North of my domain')
       END IF
    END IF

    !North
    IF (globalSubdomain%N <= processorDomain%N) THEN
       IF (globalSubdomain%N >= processorDomain%S ) THEN
          localSubdomain%N=globalSubdomain%N-processorDomain%S+1
       ELSE
          CALL IOS_Ereport('getLocalSubdomainBoundsCyclic',99,&
              'North boundary is South of my domain')
       END IF
    END IF

    ! East west wraping is a given for the cyclic case
    
    IF (region==WesterlyZone)THEN
      ! The westerly zone is the easterly bit of the local cpu's
      ! data

      ! East : The  boundary must not move, otherwise
      !        we aren't cyclic at all!

      ! West :
      localSubdomain%W=globalSubdomain%W-processorDomain%W+1 
    ELSE IF (region==EasterlyZone)THEN
      ! The easterly zone is the westerly bit of the local cpu's
      ! data

      ! West : The  boundary must not move, otherwise
      !        we aren't cyclic at all!

      ! East : compute the on processor column of the easterly edge       
      East=globalSubdomain%E-modelDomain%E
      localSubdomain%E=East-processorDomain%W+1       
    ELSE
      CALL IOS_Ereport('getLocalSubdomainBoundsCyclic',99,&
          'Bad region: must provide WesterlyZone or EasterlyZone')
      
    END IF
    
    IF (lhook) CALL dr_hook &
        ('IOS_GEOMETRY_UTILS:GETLOCALSUBDOMAINBOUNDSCYCLIC', &
        zhook_out,zhook_handle)

  END SUBROUTINE getLocalSubdomainBoundsCyclic


  ! Describe how the box, box_in intersects with my local domain,
  ! the options are 'it doesnt','completely', or some number 
  ! of overlaps (with 1 processor, a global model the most overlaps 
  ! is 4. For simplicity we will assume that only EW wrapping will occur. 
  FUNCTION classify_subdomain(cpu,fld,box_in) RESULT(R)
    TYPE(box), INTENT(IN)  :: box_in   ! The domain
    INTEGER, INTENT(IN)    :: fld
    INTEGER, INTENT(IN)    :: cpu
    INTEGER                :: R
    LOGICAL                :: ew_wrap  ! Does domain wraps the meridian
    TYPE(box)              :: bx       ! Local copy of box_in
    TYPE(box)              :: me       ! My processor domain
    INTEGER                :: bxWrapEast
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook                       &
        ('IOS_GEOMETRY_UTILS:CLASSIFY_SUBDOMAIN', &
        zhook_in,zhook_handle)

    R=not_computed
    bx=box_in
    
    me%N=offset_map(2,fld,cpu)+size_map(2,fld,cpu)
    me%S=offset_map(2,fld,cpu)
    me%E=offset_map(1,fld,cpu)+size_map(1,fld,cpu)
    me%W=offset_map(1,fld,cpu)
    
    CALL Normalize_Box(bx,fld,ew_wrap)
    
    IF (ew_wrap) THEN
      bxWrapEast=bx%E-atm_global_points(1,fld)
     
      IF (                                        &
          bx%S <= me%S .AND. bx%N >= me%N .AND.   &
          (bx%W <= me%W .OR. bxWrapEast  >= me%E)) THEN
        ! The box encloses this processor
        R=complete_intersection
      ELSE IF ( &
          bx%S > me%N .OR. bx%N < me%S .OR. &
          (bx%W > me%E .AND. bxWrapEast < me%W ) ) THEN
        R=no_intersection
      ELSE IF ( &
          bx%W <= me%E .AND. bxWrapEast >= me%W ) THEN
        R=2
      ELSE 
        R=1
      END IF
    ELSE ! No wrapping
      IF ( &
          bx%S <= me%S .AND. bx%N >= me%N .AND. &
          bx%W <= me%W .AND. bx%E >= me%E) THEN
        ! The box encloses this processor
        R=complete_intersection
      ELSE IF ( &
          bx%S > me%N .OR. bx%N < me%S .OR. &
          bx%W > me%E .OR. bx%E < me%W ) THEN
        R=no_intersection         
      ELSE ! Without wrapping there can be at most one patch
        R=1
      END IF
    END IF

    IF (lhook) CALL dr_hook &
        ('IOS_GEOMETRY_UTILS:CLASSIFY_SUBDOMAIN', &
        zhook_out,zhook_handle)

  END FUNCTION classify_subdomain

  ! Put the coordinates in a standard form such that 
  ! the east bound is off-map and east of west and check for
  ! errors
  SUBROUTINE Normalize_Box(bx,fld,ew_wrap)

    IMPLICIT NONE
    TYPE(box), INTENT(INOUT) :: bx
    INTEGER, INTENT(IN)      :: fld
    LOGICAL, INTENT(OUT)     :: ew_wrap
    REAL(KIND=jprb)          :: zhook_handle

    IF (lhook) CALL dr_hook &
        ('IOS_GEOMETRY_UTILS:NORMALIZE_BOX', &
        zhook_in,zhook_handle)

    ew_wrap=.FALSE.

    IF (bx%W > bx%E) THEN
       bx%E=bx%E+atm_global_points(1,fld)      
    END IF

! If east is off-map, set the flag
    IF (bx%E > atm_global_points(1,fld))&
         ew_wrap=.TRUE.

! Check for other wierd conditions (we don't handle these yet)

    IF (bx%S > bx%N) THEN
      WRITE(geom_util_message,'(A,I6,A,I6,A)') &
          'Normalize Box: South(',bx%S,        &
          ') greater than North(',bx%N,        &
          ') condition not supported yet'
      CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
    END IF
    
    IF (bx%S < 1) THEN
      WRITE(geom_util_message,'(A,I6,A)') &
          'Normalize Box: South(',bx%S,   &
          ') is off map'
      CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
    END IF
    
    IF (bx%N > atm_global_points(2,fld)) THEN
      WRITE(geom_util_message,'(A,I6,A,I6,A)') &
          'Normalize Box: North(',bx%N,        &
          ') is off map (',atm_global_points(2,fld),')'
      CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )
    END IF
    
    IF (bx%W > atm_global_points(1,fld)) THEN
      WRITE(geom_util_message,'(A,I6,A,I6,A)') &
          'Normalize Box: West(',bx%W,         &
          ') is off map (',atm_global_points(1,fld),')'
      CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )  
    END IF

    IF (bx%W < 1) THEN
      WRITE(geom_util_message,'(A,I6,A,I6,A)') &
          'Normalize Box: West(',bx%W,         &
          ') is off map (',atm_global_points(1,fld),')'
      CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )  
    END IF
          
    IF (ew_wrap) THEN
      IF (bx%E <= atm_global_points(1,fld)) THEN
        WRITE(geom_util_message,'(A,I6,A,I6,A)') &
            'Normalize Box: East WRAP(',bx%E,    &
            ') is off map (',atm_global_points(1,fld),')'
        CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )  
      END IF
      IF (bx%E > 2*atm_global_points(1,fld)) THEN
        WRITE(geom_util_message,'(A,I6,A,I6,A)') &
            'Normalize Box: East WRAP(',bx%E,    &
            ') is off map (',atm_global_points(1,fld),')'
        CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )  
      END IF
      
    ELSE
      IF (bx%E < 1) THEN
        WRITE(geom_util_message,'(A,I6,A,I6,A)') &
            'Normalize Box: East(',bx%E,         &
            ') is off map (',atm_global_points(1,fld),')'
        CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )  
      END IF
      IF (bx%E > atm_global_points(1,fld)) THEN
        WRITE(geom_util_message,'(A,I6,A,I6,A)') &
            'Normalize Box: East(',bx%E,         &
            ') is off map (',atm_global_points(1,fld),')'
        CALL IOS_Ereport('IOS_Geometry_utils', 99, geom_util_message )  
      END IF
      
    END IF
    IF (lhook) CALL dr_hook &
        ('IOS_GEOMETRY_UTILS:NORMALIZE_BOX', &
        zhook_out,zhook_handle)
  END SUBROUTINE Normalize_Box

END MODULE IOS_Geometry_utils

