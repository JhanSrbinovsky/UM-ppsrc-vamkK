! *****************************COPYRIGHT******************************* 
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT******************************* 
!                                                                       
!  Purpose: To assign dry deposition rates in s-1 to array dpd in       
!           ASAD module asad_mod. This is passed for use in UKCA chemistry
!                                                                       
!           Original version from Cambridge TOMCAT model  
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!                                                                       
!           Called from ASAD_CDRIVE.                                    
!                                                                       
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!                                                                       
! Code description:                                                     
!   Language: FORTRAN 90                                                
!   This code is written to UMDP3 v6 programming standards.             
!                                                                       
! --------------------------------------------------------------------- 
!                                                                       
        SUBROUTINE UKCA_DRYDEP(nlev, dryrt, n_points)                   

        USE ASAD_MOD,       ONLY: ndepd, nldepd, dpd
        USE ukca_option_mod,ONLY: L_ukca_intdd, jpdd
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE Control_Max_Sizes
        IMPLICIT NONE                                                   
                                                                        
        INTEGER, INTENT(IN) :: nlev                 ! Level number      
        INTEGER, INTENT(IN) :: n_points             ! No of spatial poin
                                                                        
        REAL, INTENT(IN) :: dryrt(n_points,jpdd)    ! Dry deposition rat
                                                                        
!       Local variables                                                 
                                                                        
        INTEGER :: i                                ! Loop variable     
        INTEGER :: js                               ! Loop variable     
        INTEGER :: nspec                            ! Pointer for specie

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('UKCA_DRYDEP',zhook_in,zhook_handle)
                                                                        
        IF (nlev == 1 .OR. L_ukca_intdd) THEN
! Bottom level or interactive dry deposition      
                                                                        
          DO js = 1,ndepd                                               
            nspec = nldepd(js)                                          
            DO i = 1,n_points                                           
              dpd(i,nspec) = dryrt(i,js)                                
            ENDDO                                                       
          ENDDO                                                         
                                                                        
        ELSE
! All other levels when not interactive 
                                                                        
          DO js = 1,ndepd                                               
            nspec = nldepd(js)                                          
            DO i = 1,n_points                                           
              dpd(i,nspec) = 0.0                                        
            ENDDO                                                       
          ENDDO                                                         
                                                                        
        ENDIF                                                           
                                                                        
        IF (lhook) CALL dr_hook('UKCA_DRYDEP',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE UKCA_DRYDEP    
