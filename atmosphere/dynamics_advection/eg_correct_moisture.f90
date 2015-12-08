! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.  
! For further details please refer to the file COPYRIGHT.txt  
! which you should have received as part of this distribution.  
! *****************************COPYRIGHT*******************************
!  
  MODULE eg_conserv_moist_mod  
  IMPLICIT NONE

! Description: 
!            This routine enforces the mass conservation for
!            the 6 moisture variables. It should work with
!            either mixing ratios and dry density or with
!            specific humdities and wet density. 
!            If L_conserve_mass=.FALSE., then the routine make 
!            sure that the moisture variables are above the 
!            minimum values required only (doesn't impose
!            any mass conservation constraint). 
!  
! Method: ENDGame formulation version 3.02
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Dynamics Advection
!  
! Code description:
! Language: Fortran 90.  
! This code is written to UMDP3 standards.  

  CONTAINS  
  SUBROUTINE eg_correct_moisture(rho_n, rho_np1,       &
                 q1_n, q2_n, q3_n, q4_n, q5_n, q6_n,   &
     q1_np1, q2_np1, q3_np1, q4_np1, q5_np1, q6_np1,   &
                q1_s, q2_s, q3_s, q4_s, q5_s, q6_s,    &
      q1min, L_q4, L_q5, L_q6, L_conserve_mass, mype,  &
                        L_conserv_smooth_lap           ) 
  
       USE eg_check_conserv_mod
       USE eg_mass_conserv_mod
       USE PrintStatus_mod
       USE Field_Types
       USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims,        &
                                        tdims_s, array_dims

       USE parkind1,              ONLY: jpim, jprb       !DrHook
       USE yomhook,               ONLY: lhook, dr_hook   !DrHook  
       USE eg_swap_bounds_mod

       IMPLICIT NONE      
       INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
       INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
       REAL(KIND=jprb)               :: zhook_handle
            
       INTEGER, INTENT(IN) ::  mype 
       LOGICAL, INTENT(IN) :: L_conserv_smooth_lap  
  
       REAL, INTENT(IN)   ::   rho_n(pdims_s%i_start:pdims_s%i_end,   &
                                     pdims_s%j_start:pdims_s%j_end,   &
                                     pdims_s%k_start:pdims_s%k_end)
       REAL, INTENT(IN)   :: rho_np1(pdims_s%i_start:pdims_s%i_end,   &
                                     pdims_s%j_start:pdims_s%j_end,   &
                                     pdims_s%k_start:pdims_s%k_end)     

       REAL, INTENT(IN)   ::    q1_n(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN)   ::    q2_n(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN)   ::    q3_n(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN)   ::    q4_n(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)

       REAL, INTENT(IN)   ::    q5_n(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::      q6_n(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::      q1_s(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::      q2_s(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::      q3_s(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::      q4_s(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::      q5_s(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(IN) ::     q6_s (tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end) 
       
       REAL, INTENT(INOUT) :: q1_np1(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(INOUT) :: q2_np1(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(INOUT) :: q3_np1(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(INOUT) :: q4_np1(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(INOUT) :: q5_np1(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       REAL, INTENT(INOUT) :: q6_np1(tdims_s%i_start:tdims_s%i_end,   &
                                     tdims_s%j_start:tdims_s%j_end,   &
                                     tdims_s%k_start:tdims_s%k_end)
       
       REAL,INTENT(IN)     :: q1min
       LOGICAL, INTENT(IN) ::  L_q4, L_q5, L_q6, L_conserve_mass 
      
       ! locals
    
       INTEGER :: number_qs, k

       REAL,    ALLOCATABLE ::   qs_n(:,:,:,:)
       REAL,    ALLOCATABLE :: qs_np1(:,:,:,:)
       REAL,    ALLOCATABLE ::   qs_s(:,:,:,:)

       REAL,    ALLOCATABLE ::     qsmin(:)  
       INTEGER, ALLOCATABLE :: qswitches(:)  

       TYPE (array_dims) :: superdim

       IF (lhook) CALL dr_hook('EG_CORRECT_MOISTURE',zhook_in,zhook_handle)
        
       IF ( L_conserve_mass ) THEN

         number_qs = 6
         ALLOCATE(    qs_n(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,number_qs))

         ALLOCATE(  qs_np1(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,number_qs))

         ALLOCATE(    qs_s(tdims_s%i_start:tdims_s%i_end,             &
                           tdims_s%j_start:tdims_s%j_end,             &
                           tdims_s%k_start:tdims_s%k_end,number_qs))

         ALLOCATE(     qsmin(number_qs) )     
         ALLOCATE( qswitches(number_qs) )
   
! copy the 6 moisture variable fields and sources into
! superarrays of dimension 6
         
         DO k = tdims_s%k_start, tdims_s%k_end 
           qs_n(:,:,k,1) = q1_n(:,:,k)
           qs_n(:,:,k,2) = q2_n(:,:,k)
           qs_n(:,:,k,3) = q3_n(:,:,k)
           qs_n(:,:,k,4) = q4_n(:,:,k)
           qs_n(:,:,k,5) = q5_n(:,:,k)
           qs_n(:,:,k,6) = q6_n(:,:,k)
   
           qs_s(:,:,k,1) = q1_s(:,:,k)
           qs_s(:,:,k,2) = q2_s(:,:,k)
           qs_s(:,:,k,3) = q3_s(:,:,k)
           qs_s(:,:,k,4) = q4_s(:,:,k)
           qs_s(:,:,k,5) = q5_s(:,:,k)
           qs_s(:,:,k,6) = q6_s(:,:,k)
      
           qs_np1(:,:,k,1) = q1_np1(:,:,k)
           qs_np1(:,:,k,2) = q2_np1(:,:,k)
           qs_np1(:,:,k,3) = q3_np1(:,:,k)
           qs_np1(:,:,k,4) = q4_np1(:,:,k)
           qs_np1(:,:,k,5) = q5_np1(:,:,k)
           qs_np1(:,:,k,6) = q6_np1(:,:,k)                        
         END DO
   
! make sure the boundary condition filed(surface)=field(level 1) is imposed.
! This may be not necessary here because it may already be done higher up
! but for safety we keep it.

         qs_n(:,:,tdims_s%k_start,:) =   qs_n(:,:,tdims_s%k_start+1,:) 
         qs_s(:,:,tdims_s%k_start,:) =   qs_s(:,:,tdims_s%k_start+1,:) 
         qs_np1(:,:,tdims_s%k_start,:) = qs_np1(:,:,tdims_s%k_start+1,:) 

! set switches to 1 (means impose mass conservation). usuallay the first
! 3 moisture variables are always there but the last 3 are optional
   
         qswitches(1:3) = 1
         qswitches(4:5) = 0 

         IF ( L_q4 )  qswitches(4) = 1
         IF ( L_q5 )  qswitches(5) = 1
         IF ( L_q6 )  qswitches(6) = 1
   
         qsmin(1)           = q1min
         qsmin(2:number_qs) = 0.0
 
! print mass conservation error before correction if needed 
   
         IF (printstatus > prstatus_normal) THEN
           IF( mype == 0) THEN
             WRITE(6,fmt='(A)') 'Error in mass conservation for',     &
                                ' moisture before correction'
           END IF            
           CALL eg_check_mass_conservation(rho_n, rho_np1, qs_n,      &
                                           qs_np1, qs_s,              & 
                      qswitches, number_qs, mype,'eg_correct_moisture')
         END IF

! enfore mass conservation
    
         CALL eg_mass_conservation(rho_n, rho_np1, qs_n, qs_np1, qs_s, &
                                   qsmin, qswitches, number_qs,        &
                                   L_conserv_smooth_lap                )
      
         ! may be there is no need for this swap-bounds if the moisture 
         ! fileds are swap-bounded later on -- but it's safer to leave it 
 
         superdim         = tdims_s
         superdim%k_start = 1
         superdim%k_end   = number_qs*(tdims%k_end-tdims%k_start+1)

         CALL eg_swap_bounds(qs_np1, superdim, fld_type_p, .FALSE.) 

! put the corrected moisture variables back into their original arrays
! unpack them from the superarray

         DO k = tdims_s%k_start, tdims_s%k_end 
           q1_np1(:,:,k) = qs_np1(:,:,k,1) 
           q2_np1(:,:,k) = qs_np1(:,:,k,2)
           q3_np1(:,:,k) = qs_np1(:,:,k,3) 
           q4_np1(:,:,k) = qs_np1(:,:,k,4) 
           q5_np1(:,:,k) = qs_np1(:,:,k,5)
           q6_np1(:,:,k) = qs_np1(:,:,k,6)     
         END DO

! print mass conservation error after correction if needed 
   
         IF (printstatus > prstatus_normal) THEN
           IF( mype == 0) THEN
             WRITE(6,fmt='(A)')'Error in mass conservation for',      &
                               ' moisture after correction'
           END IF     
           CALL eg_check_mass_conservation(rho_n, rho_np1, qs_n,      &
                                           qs_np1, qs_s,              & 
                      qswitches, number_qs, mype,'eg_correct_moisture')
         END IF
    

         DEALLOCATE (qswitches)
         DEALLOCATE (qsmin)
         DEALLOCATE (qs_s)
         DEALLOCATE (qs_np1)
         DEALLOCATE (qs_n)
      
       ELSE

! if L_conserve_mass==.false. then simply make sure that the moisture
! variables are above the minimum values required

         q1_np1 = MAX(q1_np1, q1min)
         q2_np1 = MAX(q2_np1, 0.0  )
         q3_np1 = MAX(q3_np1, 0.0  )
         q4_np1 = MAX(q4_np1, 0.0  )
         q5_np1 = MAX(q5_np1, 0.0  )
         q6_np1 = MAX(q6_np1, 0.0  )

       END IF  ! L_conserve_mass
  
       IF (lhook) CALL dr_hook('EG_CORRECT_MOISTURE',zhook_out,zhook_handle)
       
       RETURN
       END SUBROUTINE eg_correct_moisture

END MODULE eg_conserv_moist_mod
