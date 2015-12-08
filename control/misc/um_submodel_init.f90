! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Initialise model for submodel and internal model coupling
!
! Subroutine Interface:
      SUBROUTINE UM_Submodel_Init(ErrorStatus)

      USE Submodel_mod
      USE check_iostat_mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!   UM_Submodel_Init initialises the model with information specifying
!   internal model and submodel partitions for the run, which is
!   required for control of coupling when more than one internal model
!   is present.
!
! Method:
!   The routine reads information from the user interface, providing
!   lists of internal models and their associated submodel data
!   partitions. This is required in both the reconfiguration and the
!   model as a prior step to calculating addressing in STASH_PROC.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
! Declarations:
!
!
! An alternative common block required by TYPD1
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=1

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! Subroutine arguments

!   ErrorStatus
      INTEGER      ErrorStatus          ! Error flag (0 = OK)

! Local parameters:

! Local scalars:
      INTEGER                                                           &
     & s                                                                &
                       ! submodel loop
     &,i                                                                &
                       ! internal model loop
     &,sm                                                               &
                       ! submodel identifier
     &,im                                                               &
                       ! internal model identifier
     &,sm_prev                                                          &
                       ! previous submodel identifier
     &,im_prev         ! previous internal model identifier

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays:

! Function & Subroutine calls: None

!- End of header

!      IF (lhook) CALL dr_hook('UM_SUBMODEL_INIT',zhook_in,zhook_handle)

!
! 1. Initialise lists before obtaining values for this experiment.
!
      do i=1,N_INTERNAL_MODEL_MAX
         INTERNAL_MODEL_LIST(i)      = 0
         SUBMODEL_FOR_IM(i)          = 0
      enddo   ! i over internal model list

      do im=1,INTERNAL_ID_MAX
         SUBMODEL_PARTITION_INDEX(im) = 0
         INTERNAL_MODEL_INDEX(im) = 0
         LAST_IM_IN_SM(im)=.false.
      enddo   ! im over internal model ids

      do s=1,N_SUBMODEL_PARTITION_MAX
         SUBMODEL_PARTITION_LIST(s)= 0
         SUBMODEL_FOR_SM(s)=0
      enddo  ! s over submodel list

      do sm=1,SUBMODEL_ID_MAX
         N_INTERNAL_FOR_SM(sm)      = 0
      enddo  ! sm over submodel ids

!
! 2. Obtain internal model and submodel identifiers from umui
!    generated namelist.
!

      READ (UNIT=5, NML=NSUBMODL, IOSTAT=ErrorStatus) 
      CALL check_iostat(errorstatus, "namelist NSUBMODL")      
!
!
! 3. Check umui supplied values.
!
!
! 3.1 Check for umui supplied dimensions against parameter maxima.

      if(N_INTERNAL_MODEL >  N_INTERNAL_MODEL_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Too many internal ',  &
     &   'models =',N_INTERNAL_MODEL,                                   &
     &   ' :You need to increase N_INTERNAL_MODEL_MAX'
         ErrorStatus=1       ! Set error flag
      endif
!
! 3.2 Check umui suppiled values are valid
!
      do i=1,N_INTERNAL_MODEL ! loop over internal models

        im = INTERNAL_MODEL_LIST(i) ! internal model identifier
        if(im <= 0.or.im >  INTERNAL_ID_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Illegal internal ',   &
     &   'model identifier=',im,                                        &
     &   ' :Check values in namelist NSUBMODL supplied by umui'
         ErrorStatus=1       ! Set error flag
        endif

        sm = SUBMODEL_FOR_IM(i)     ! submodel for this internal model
        if(sm <= 0.or.sm >  SUBMODEL_ID_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Illegal submodel ',   &
     &   'dump identifier=',sm,                                         &
     &   ' :Check values in namelist NSUBMODL supplied by umui'
         ErrorStatus=1       ! Set error flag
        endif

      enddo ! i=1,N_INTERNAL_MODEL
!
! 4. Form internal model and submodel description arrays.
!
      sm_prev = 0             ! Null value of submodel identifier
      N_SUBMODEL_PARTITION=0  ! Count no. of submodel partitions

      do i=1,N_INTERNAL_MODEL ! loop over internal models

        im = INTERNAL_MODEL_LIST(i) ! internal model identifier
        sm = SUBMODEL_FOR_IM(i)     ! submodel for this internal model
        INTERNAL_MODEL_INDEX(im)=i  ! sequence no. for STASH arrays

        if(sm /= sm_prev) then  ! new submodel

           N_SUBMODEL_PARTITION = N_SUBMODEL_PARTITION+1
           SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) = sm

!   Since this is a new submodel, the previous internal model must be
!   the last internal model in its submodel partition.
           IF(N_SUBMODEL_PARTITION >  1) THEN ! Not first dump
              LAST_IM_IN_SM(im_prev) = .true.
           ENDIF

        endif                   ! test on new submodel
        SUBMODEL_FOR_SM(IM) = N_SUBMODEL_PARTITION

        SUBMODEL_PARTITION_INDEX(im)=sm
        N_INTERNAL_FOR_SM(sm)=N_INTERNAL_FOR_SM(sm)+1

        im_prev=im
        sm_prev=sm

      enddo ! i=1,N_INTERNAL_MODEL

      LAST_IM_IN_SM(im) = .true.  ! last im in list is last im in sm

!
! 5. Check calculated dimensions against parameter maxima.

      if(N_SUBMODEL_PARTITION >  N_SUBMODEL_PARTITION_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Too many submodels =',&
     &   N_SUBMODEL_PARTITION,                                          &
     &   ' You need to increase N_SUBMODEL_PARTITION_MAX'
         ErrorStatus=1       ! Set error flag
      endif
!
!     Need a copy of No of submodels for use by TYPD1.
      ALT_N_SUBMODEL_PARTITION=N_SUBMODEL_PARTITION

      if (ALT_N_SUBMODEL_PARTITION_MAX /= N_SUBMODEL_PARTITION_MAX)THEN
        write(6,*)'UM_Submodel_In: Mismatch in parameters '
        WRITE(6,*)'N_SUBMODEL_PARTITION_MAX and '
        WRITE(6,*)'ALT_N_SUBMODEL_PARTITION_MAX. '
        WRITE(6,*)'They should be identical '
        ErrorStatus=1
        endif
!        IF (lhook) CALL dr_hook('UM_SUBMODEL_INIT',zhook_out,zhook_handle)
        RETURN
      END SUBROUTINE UM_Submodel_Init
