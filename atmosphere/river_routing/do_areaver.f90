! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: River Routing
MODULE do_areaver_mod

IMPLICIT NONE

CONTAINS


      SUBROUTINE DO_AREAVER(gaps_lambda_srce,gaps_phi_srce             &
      ,data_srce,gaps_lambda_targ,gaps_phi_targ,count_targ,base_targ   &
      ,lrow_targ,want,mask_targ,index_srce,weight,data_targ            &
      ,global_src_lambda_gaps,global_src_phi_gaps, grid, recv_concern, &
      recv_size)
!LL   Subroutine DO_AREAVER_PAR ----------------------------------------
!LL
!LL Purpose:
!LL
!LL   
!LL
!LL   Performs parallel area-averaging to transform data from the source 
!LL   grid to the target grid, or adjust the values on the source grid to
!LL   have the area-averages supplied on the target grid. The latter 
!LL   mode is intended for adjusting values obtained by interpolating 
!LL   from "target" to "source" in order to conserve the area-averages.
!LL   This mode should be used ONLY if each source box belongs in
!LL   exactly one target box. ADJUST=0 selects normal area-averaging,
!LL   ADJUST=1 selects adjustment by addition (use this mode for fields
!LL   which may have either sign), ADJUST=2 selects adjustment by
!LL   multiplication (for fields which are positive-definite or
!LL   negative-definite).
!LL
!LL   For two-way conservative coupling, ADJUST=3 makes an adjustment
!LL   field for fields which may have either sign, ADJUST=4 makes an
!LL   adjustment field for fields which are positive-definite or
!LL   negative-definite, ADJUST=5 performs conservative adjustment
!LL   by addition (use this mode for fields which may have either sign)
!LL   and ADJUST=6 selects conservative adjustment by multiplication
!LL   (for fields which are positive-definite or negative-definite).
!LL
!LL   The shape of the source and target grids are specified by their
!LL   dimensions GAPS_aa_bb, which give the number of gaps in the
!LL   aa=LAMBDA,PHI coordinate in the bb=SRCE,TARG grid. (The product
!LL   of GAPS_LAMBDA_bb and GAPS_PHI_bb is the number of boxes in the
!LL   bb grid.)
!LL
!LL   The input and output data are supplied as 2D arrays DATA_SRCE and
!LL   DATA_TARG, whose first dimensions should also be supplied. Speci-
!LL   fying these sizes separately from the actual dimensions of the
!LL   grids allows for columns and rows in the arrays to be ignored.
!LL   A target land/sea mask should be supplied in MASK_TARG, with the
!LL   value indicating wanted points specified in WANT. Points which
!LL   are unwanted or which lie outside the source grid are not altered
!LL   in DATA_TARG. DATA_SRCE can optionally be supplied with its rows
!LL   in reverse order (i.e. with the first row corresponding to
!LL   minimum LAMBDA).
!LL
!LL   The arrays COUNT_TARG, BASE_TARG, INDEX_SRCE and WEIGHT should be
!LL   supplied as returned by PRE_AREAVER q.v.
!LL
!LL   Programming Standard, UMDP3 vn8.3
!LL
!LLEND -----------------------------------------------------------------
!

        USE regrid_utils, ONLY: global_to_local_gridpt, find_value,      &
             gridpt_outside_proc_domain
        USE regrid_types
        USE UM_ParVars, ONLY: mype, g_blsize, fld_type_p
        
        
        IMPLICIT NONE
        
      INTEGER, INTENT(IN) ::                                            &
       gaps_lambda_srce                                                 &
                               !   number lambda gaps in source grid
      ,gaps_phi_srce                                                    &
                               !   number phi gaps in source grid
      ,gaps_lambda_targ                                                 &
                               !   number lambda gaps in target grid
      ,gaps_phi_targ                                                    &
                               !   number phi gaps in target grid
      ,lrow_targ                                                        &
                               !   first dimension of target arrays
      ,count_targ(gaps_lambda_targ,gaps_phi_targ)                       &
!                              !   no. of source boxes in target box
      ,base_targ(gaps_lambda_targ,gaps_phi_targ)                        &
!                              !   first index in list for target box
      ,index_srce(*)                                                    &
                               !   list of source box indices
      ,global_src_lambda_gaps                                           &
                               !   lambda gaps in global src grid
      ,global_src_phi_gaps                                              &
                               !   phi gaps in global src grid
      ,grid                                                             &
                               ! grid type (e.g. atmos, river ...)
      ,recv_size 

      LOGICAL, INTENT(IN) ::                                            &
       want                                                             &
                               !   indicator of wanted points in mask
      ,mask_targ(gaps_lambda_targ, gaps_phi_targ)  
      !   land/sea mask for target grid
    
      REAL, INTENT(IN) ::                                               &
       data_srce(gaps_lambda_srce,gaps_phi_srce)                                          
                               !    data on source grid

      REAL, INTENT(IN) ::                                               &
       weight(*)                                                        
                               !   list of weights for source boxes

      REAL, INTENT(OUT) ::                                              &
       data_targ(gaps_lambda_targ, gaps_phi_targ)  
                               !    data on target grid

      TYPE(CONCERN), INTENT(IN) :: recv_concern(recv_size)

! local variables 

      INTEGER                                                           &
       ip                                                               &
                               ! pointer into lists
      ,i                                                                &
                               ! loop index
      ,ix1(gaps_lambda_srce*gaps_phi_srce)                              &
!                              ! working srce lambda indices
      ,iy1(gaps_lambda_srce*gaps_phi_srce)                              &
!                              ! working srce phi indices
      ,ix2,iy2                                                          &
                               ! working targ lambda/phi indices         
      ,src_phi_gap                                                      &
                               ! phi gap in global domain src grid
      ,src_lambda_gap
                               ! lambda gap in global domain src grid

      INTEGER                                                           &
       lambda_gap0, lambda_gapf, phi_gap0, phi_gapf, xSrc, ySrc                     
            ! the local start and end points in the src grid 

      REAL                                                              &
       temp_targ, data_src                                                      

      LOGICAL found, not_within 

!
!     loop over all target boxes and calculate values as required.
!
!     the weights and source box indices are recorded in continuous
!     lists. count_targ indicates how many consecutive entries in these
!     lists apply to each target box.
!
      DO iy2=1,gaps_phi_targ
         DO ix2=1,gaps_lambda_targ
            IF (mask_targ(ix2,iy2).EQV.want) THEN
               IF (count_targ(ix2,iy2) /= 0) THEN
                  temp_targ=0.
                  DO i=1,count_targ(ix2,iy2)
                     ip=base_targ(ix2,iy2)+i
                     
                     !! determine if src requested is within boundary 
                     !! in global src grid 
                     src_lambda_gap = MOD(index_srce(ip)-1,             &
                          global_src_lambda_gaps)+1
                     src_phi_gap = (index_srce(ip)-1)/                  &
                          global_src_lambda_gaps+1
                     

                     !! check if src point is in proc domain
                     not_within = GRIDPT_OUTSIDE_PROC_DOMAIN(           &
                     src_lambda_gap, src_phi_gap, grid)

                     data_src = 0
                     
                     !! if not check in your recv Concerns
                     IF(not_within) THEN
                      
                        found = .FALSE.
                        CALL FIND_VALUE(src_lambda_gap, src_phi_gap,    &
                             recv_size, recv_concern, data_src, found)
                        IF(.NOT. found) print*, "Error, Src Not Found!!"

                     ELSE
                        xsrc = src_lambda_gap
                        ysrc = src_phi_gap
                        CALL GLOBAL_TO_LOCAL_GRIDPT(xsrc, ysrc, grid)
                        data_src = data_srce(xsrc, ysrc)        
                    
                     END IF
                     !! apply weighting 
                     temp_targ=temp_targ+weight(ip)*data_src                
                  END DO
               ELSE
                  !! one to one 
                  temp_targ=data_targ(ix2,iy2)
               END IF
               data_targ(ix2,iy2)=temp_targ
            END IF
         END DO
      END DO

      RETURN
    END SUBROUTINE DO_AREAVER
END MODULE do_areaver_mod
