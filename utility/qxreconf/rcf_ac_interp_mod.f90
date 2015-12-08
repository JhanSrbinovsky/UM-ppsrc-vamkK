! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Interpolations between A and C grids

MODULE Rcf_AC_interp_Mod

!  Subroutine Rcf_U_to_P
!  Subroutine Rcf_V_to_P
!  Subroutine Rcf_P_to_U
!  Subroutine Rcf_P_to_V
!
! Description:
!   Interpolates full-level (ie no horizontal decomposition) fields
!   between A and C grids.
!
! Method:
!   Linear interpolation in relevant (x or y) direction.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------
! Routine to interpolate U field onto P points
!------------------------------------------------------------------
SUBROUTINE Rcf_U_to_P( u_level, u_at_p, grid)

USE Rcf_Grid_Type_Mod, ONLY : &
    grid_type

USE UM_ParVars, ONLY : &
    bound,                  &
    bc_cyclic
    
USE Rcf_UMhead_Mod, ONLY : &
    um_header_type
   
USE Rcf_Headaddress_Mod, ONLY :       &
    FH_GridStagger,                   &
    FH_GridStagger_Endgame

IMPLICIT NONE

! Arguments
TYPE( grid_type ), INTENT(IN)     :: grid
REAL,              INTENT(IN)     :: u_level( grid % glob_u_row_length,&
                                              grid % glob_u_rows )
REAL,              INTENT(OUT)    :: u_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )

! Local variables
INTEGER     :: i
INTEGER     :: j
LOGICAL                           :: L_endgame

L_endgame  = grid % grid_stagger == FH_GridStagger_Endgame


IF (L_endgame) THEN

  DO j = 1, grid % glob_u_rows
    DO i = 1,grid % glob_p_row_length - 1
      u_at_p(i,j) = ( u_level(i+1,j) +  u_level(i,j) ) * 0.5
    END DO
  END DO
  
  ! assume this module only employed for LAMs (rotated grids)
  IF (bound(1) == bc_cyclic) THEN  
    DO j = 1, grid % glob_u_rows
       i = grid % glob_p_row_length
       u_at_p(i,j) = ( u_level(1,j) +  u_level(i,j) ) * 0.5    
    END DO 
    
  ELSE
    ! copy for end P pt.
    DO j = 1, grid % glob_u_rows
       i = grid % glob_p_row_length
       u_at_p(i,j) = u_level(i,j) 
    END DO  
  
  END IF
  
  
ELSE
  
  DO j = 1, grid % glob_u_rows
    DO i = 2,grid % glob_u_row_length
      u_at_p(i,j) = ( u_level(i-1,j) +  u_level(i,j) ) * 0.5
    END DO

!  If ( grid % glob_u_row_length == grid % glob_p_row_length ) Then
  ! Wrapping LAM

!   u_at_p(1,j) =(u_level(1,j) + u_level( grid % glob_u_row_length, j))&
!                 * 0.5

!  Else    ! non-wrapping

    ! Copy for end P points.
    u_at_p(1,j) = u_level(1,j)
!  End If
  END DO
END IF  

RETURN
END SUBROUTINE Rcf_U_to_P

!------------------------------------------------------------------
! Routine to interpolate V field onto P points
!------------------------------------------------------------------
SUBROUTINE Rcf_V_to_P( v_level, v_at_p, grid )

USE Rcf_Grid_Type_Mod, ONLY : &
    grid_type

USE UM_ParVars, ONLY : &
    bound,                  &
    bc_cyclic

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type
   
USE Rcf_Headaddress_Mod, ONLY :       &
    FH_GridStagger,                   &
    FH_GridStagger_Endgame


IMPLICIT NONE

! Arguments
TYPE( grid_type ), INTENT(IN)     :: grid
REAL,              INTENT(IN)     :: v_level( grid % glob_v_row_length,&
                                              grid % glob_v_rows )
REAL,              INTENT(OUT)    :: v_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )

! Local variables
INTEGER                           :: i
INTEGER                           :: j
LOGICAL                           :: L_endgame

L_endgame  = grid % grid_stagger == FH_GridStagger_Endgame

IF (L_endgame) THEN

  DO i = 1, grid % glob_v_row_length
    DO j = 1,grid % glob_p_rows-1
      v_at_p(i,j) = (v_level(i, j+1) + v_level(i,j) ) * 0.5
    END DO
  END DO

  IF (bound(2) == bc_cyclic) THEN
    DO i = 1, grid % glob_v_row_length
     j = grid % glob_p_rows
     v_at_p(i,j) = (v_level(i, 1) + v_level(i,j) ) * 0.5
    END DO  
  ELSE
    DO i = 1, grid % glob_v_row_length
       j = grid % glob_p_rows
      v_at_p(i,j) = (v_level(i,j+1) + v_level(i,j) ) * 0.5
    END DO  
  END IF

ELSE

  DO i = 1, grid % glob_v_row_length
    DO j = 2,grid % glob_v_rows
      v_at_p(i,j) = (v_level(i, j-1) + v_level(i,j) ) * 0.5
    END DO

    ! Copy for end points

    IF (bound(2) == bc_cyclic) THEN
      v_at_p(i,1)=(v_level(i,1)+v_level(i,grid % glob_v_rows))* 0.5
    ELSE
      v_at_p(i,1) = v_level(i,1)
      v_at_p(i,grid % glob_p_rows) = v_level(i,grid % glob_v_rows)
    END IF
  END DO
END IF
RETURN
END SUBROUTINE Rcf_V_to_P

!------------------------------------------------------------------
! Routine to interpolate P field onto U points
!------------------------------------------------------------------
SUBROUTINE Rcf_P_to_U( u_at_p, u_level, grid )

USE Rcf_Grid_Type_Mod, ONLY : &
    grid_type

USE UM_ParVars, ONLY : &
    bound,                  &
    bc_cyclic

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type
   
USE Rcf_Headaddress_Mod, ONLY :       &
    FH_GridStagger,                   &
    FH_GridStagger_Endgame

IMPLICIT NONE

! Arguments
TYPE( grid_type ), INTENT(IN)     :: grid
REAL,              INTENT(IN)     :: u_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )
REAL,              INTENT(OUT)    :: u_level( grid % glob_u_row_length,&
                                              grid % glob_u_rows )

! Local variables
INTEGER                           :: i
INTEGER                           :: j
LOGICAL                           :: L_endgame

L_endgame  = grid  % grid_stagger == FH_GridStagger_Endgame

IF (L_endgame) THEN
  DO j = 1, grid % glob_u_rows
    DO i = 2, grid % glob_p_row_length 
      u_level(i,j) = ( u_at_p(i,j) +  u_at_p(i-1,j) ) * 0.5
    END DO
  END DO
  
  IF(bound(1) == bc_cyclic)THEN
    DO j = 1, grid % glob_u_rows
      u_level(1,j) = ( u_at_p(1,j) +               &
                       u_at_p(grid % glob_p_row_length,j) ) * 0.5
    END DO
  ELSE  
    DO j = 1, grid % glob_u_rows
      ! COPY end points
      u_level(1,j) = u_at_p(1,j)
    END DO
  END IF
ELSE

  DO j = 1, grid % glob_u_rows
    DO i = 1, grid % glob_p_row_length - 1
      u_level(i,j) = ( u_at_p(i,j) +  u_at_p(i+1,j) ) * 0.5
    END DO

    IF (bound(1) == bc_cyclic) THEN
      u_level( grid % glob_u_row_length, j ) =                        &
        (u_at_p( grid % glob_p_row_length, j ) + u_at_p(1, j) ) * 0.5
    END IF
  END DO
  
  ! copy data into last u pt.
  DO j = 1, grid % glob_u_rows
      ! COPY end points
      u_level(grid % glob_u_row_length,j) =        &
                       u_at_p(grid % glob_p_row_length,j)
  END DO  
END IF

RETURN
END SUBROUTINE Rcf_P_to_U

!------------------------------------------------------------------
! Routine to interpolate P field onto V points
!------------------------------------------------------------------
SUBROUTINE Rcf_P_to_V( v_at_p, v_level, grid )

USE Rcf_Grid_Type_Mod, ONLY : &
    grid_type

USE UM_ParVars, ONLY : &
    bound,                  &
    bc_cyclic

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type
   
USE Rcf_Headaddress_Mod, ONLY :       &
    FH_GridStagger,                   &
    FH_GridStagger_Endgame

IMPLICIT NONE

! Arguments
TYPE( grid_type ), INTENT(IN)     :: grid
REAL,              INTENT(IN)     :: v_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )
REAL,              INTENT(OUT)    :: v_level( grid % glob_v_row_length,&
                                              grid % glob_v_rows )

! Local variables
INTEGER                           :: i
INTEGER                           :: j
LOGICAL                           :: L_endgame

L_endgame  = grid  % grid_stagger == FH_GridStagger_Endgame

IF (L_endgame) THEN
  DO i = 1, grid % glob_v_row_length
    DO j = 2, grid % glob_p_rows
    v_level(i,j) = ( v_at_p(i,j) + v_at_p(i,j-1) ) * 0.5
    END DO
  END DO
  
  IF(bound(2) == bc_cyclic)THEN
    DO i = 1, grid % glob_v_row_length
      v_level(i,1) = ( v_at_p(i,1) + v_at_p(i,grid % glob_p_rows) ) * 0.5
    END DO
  ELSE
  
    DO i = 1, grid % glob_v_row_length      
      ! Copy for end points
      v_level(i,1) = v_at_p(i,1) 
      v_level(i,grid % glob_v_rows)=v_at_p(i,grid % glob_p_rows)      
    END DO
  END IF

ELSE

  DO i = 1, grid % glob_v_row_length
    DO j = 1, grid % glob_v_rows - 1
    v_level(i,j) = ( v_at_p(i,j) + v_at_p(i,j+1) ) * 0.5
    END DO
  END DO
  IF(bound(2) == bc_cyclic)THEN
    DO i = 1, grid % glob_v_row_length
      v_level(i,grid % glob_v_rows)=(v_at_p(i,grid % glob_v_rows) + &
                             v_at_p(i,1) ) *0.5
    END DO
  ELSE
    DO i = 1, grid % glob_v_row_length
      j = grid % glob_v_rows
      v_level(i,j) = ( v_at_p(i,j) + v_at_p(i,j+1) ) * 0.5
    END DO
  END IF

END IF
RETURN
END SUBROUTINE Rcf_P_to_V

END MODULE Rcf_AC_interp_Mod
