! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP.

MODULE decomp_db
  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE field_types
  USE decomp_params
  USE domain_params
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParCore
  USE UM_ParParams
  USE Atmos_Max_Sizes, ONLY : &
      rows_max,               &
      row_length_max,         &
      max_halo_size
  USE gcom_mod
  
  IMPLICIT NONE
  
  INTERFACE decompose
    MODULE PROCEDURE      &
        decompose_smexe,  &
        decompose_full
  END INTERFACE
  
  ! Global data describing decompositions.
  TYPE Decomp_DB_type
    INTEGER :: bound    (ndim_max)                     =0
    INTEGER :: glsize   (ndim_max, nfld_max)           =0
    INTEGER :: gridsize (ndim_max)                     =0
    INTEGER :: halosize (ndim_max, nhalo_max)          =0
    INTEGER :: neighbour(4)                            =0

    INTEGER :: g_pe_index_EW &
        (1-Max_Halo_Size:row_length_max+Max_Halo_Size) =0
    INTEGER :: g_pe_index_NS &
        (1-Max_Halo_Size:rows_max+Max_Halo_Size)       =0

    INTEGER :: first_comp_pe       
    INTEGER :: last_comp_pe
    INTEGER :: nproc
    INTEGER :: gc_proc_row_group
    INTEGER :: gc_proc_col_group
    INTEGER :: gc_all_proc_group
    INTEGER :: sb_model_domain
    LOGICAL :: set = .FALSE.
    
    ! Shortcut aliases, used by reconfig code
    INTEGER, POINTER  :: glsizep(:)             =>NULL()
    INTEGER, POINTER  :: glsizeu(:)             =>NULL()
    INTEGER, POINTER  :: glsizev(:)             =>NULL()
    INTEGER, POINTER  :: glsizer(:)             =>NULL()
    
    ! Allocatable structures, due to 'maxproc' index
    INTEGER, POINTER :: g_lasize      (:,:,:,:) =>NULL()
    INTEGER, POINTER :: g_blsize      (:,:,:)   =>NULL()
    INTEGER, POINTER :: g_datastart   (:,:)     =>NULL()
    INTEGER, POINTER :: g_datastart_f (:,:,:)   =>NULL()
    INTEGER, POINTER :: g_datastartr  (:,:)     =>NULL()
    INTEGER, POINTER :: g_gridpos     (:,:)     =>NULL()
    
  END TYPE Decomp_DB_type

  TYPE (Decomp_DB_type), SAVE, TARGET :: DecompDB(Max_Decomps )  
  
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1
  
CONTAINS

  SUBROUTINE Allocate_Decomposition(deco)

    ! Initialise shortcuts and allocate structures. 

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: deco

    ALLOCATE(DecompDB(deco)%g_lasize      &
        (ndim_max,nfld_max,nhalo_max,0:nproc_max-1))
    ALLOCATE(DecompDB(deco)%g_blsize      (ndim_max,nfld_max,0:nproc_max-1))
    ALLOCATE(DecompDB(deco)%g_datastart   (ndim_max,0:nproc_max-1))
    ALLOCATE(DecompDB(deco)%g_datastart_f (ndim_max,nfld_max,0:nproc_max-1))
    ALLOCATE(DecompDB(deco)%g_datastartr  (ndim_max,0:nproc_max-1))
    ALLOCATE(DecompDB(deco)%g_gridpos     (ndim_max,0:nproc_max-1))

    ! Zero entries in case of unused field types.
    DecompDB(deco)%g_lasize     (:,:,:,0:nproc_max-1)=0
    DecompDB(deco)%g_blsize     (:,:,  0:nproc_max-1)=0
    DecompDB(deco)%g_datastart  (:,    0:nproc_max-1)=0
    DecompDB(deco)%g_datastart_f(:,:,  0:nproc_max-1)=0
    DecompDB(deco)%g_datastartr (:,    0:nproc_max-1)=0
    DecompDB(deco)%g_gridpos    (:,    0:nproc_max-1)=0

! Set up aliases
    DecompDB (deco) % glsizeu=>&
        DecompDB (deco) % glsize(1:ndim_max,fld_type_u)
    DecompDB (deco) % glsizev=>&
        DecompDB (deco) % glsize(1:ndim_max,fld_type_v)
    DecompDB (deco) % glsizep=>&
        DecompDB (deco) % glsize(1:ndim_max,fld_type_p)
    DecompDB (deco) % glsizer=>&
        DecompDB (deco) % glsize(1:ndim_max,fld_type_r)
  END SUBROUTINE Allocate_Decomposition

! Wrapper for small execs  
  SUBROUTINE decompose_smexe(global_row_len, global_n_rows,         &
      extended_halo_EW,extended_halo_NS,     &
      tot_levels,noSmallHalos)

    USE rimtypes, ONLY : nrima_max
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: &      
        global_row_len,    & ! number of E-W points of entire model
        global_n_rows,     & ! number of P rows of entire model
        extended_halo_EW,  & ! EW halo size
        extended_halo_NS,  & ! NS halo size
        tot_levels           ! total number of levels
    LOGICAL, OPTIONAL   :: noSmallHalos

    INTEGER ::  local_row_len
    INTEGER ::  local_n_rows      

    ! These are not a factor for decompose_smexe, so we will 
    ! just use a dummy zero array to avoid risk that the module 
    ! version has already been sensibly initialised. 
    INTEGER ::  dummy_rimwidth(nrima_max)
    dummy_rimwidth(:)=0
    
    CALL decompose(decomp_smexe,           &
        global_row_len, global_n_rows,     &
        tot_levels,                        &
        0,0,                               &
        mt_smexe,                          &
        1,1,                               &
        extended_halo_EW,extended_halo_NS, &
        dummy_rimwidth,nrima_max,          &
        local_row_len, local_n_rows,       &
        noSmallHalos)

  END SUBROUTINE decompose_smexe


! Parallel UM: Perform 2D data decomposition
  SUBROUTINE decompose_full(deco,            &
      global_row_len, global_n_rows,         &
      tot_levels,                            &
      river_rows, river_row_length,          &
      model_type,                            &
      nproc_ew, nproc_ns,                    &
      extended_halo_ew,                      &
      extended_halo_ns,                      &
      rimwidth, nrima_max,                   &
      local_row_len, local_n_rows,           &
      noSmallHalos)
    
    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE

! Description:
! This routine performs a 2D decomposition - taking the global X
! (global_row_len) and Y (global_n_rows) data sizes and decomposing
! across nproc_EW processors in the X direction and nproc_NS processors
! in the Y direction.
! The local data size is returned via local_row_len and local_n_rows.
! These values will include a data halo for boundary updates.

! Method:
! The local data sizes are calculated and stored in the COMMON block
! DECOMPDB. The boundary conditions are set (cyclic in East/West
! direction if *DEF,GLOBAL).
! Subroutine Arguments:

    INTEGER, INTENT(IN)  :: deco
    INTEGER, INTENT(IN)  ::  global_row_len 
                       ! number of E-W points of entire model
    INTEGER, INTENT(IN)  ::  global_n_rows
                       ! number of P rows of entire model
    INTEGER, INTENT(IN)  ::  tot_levels 
                       ! total number of levels
    INTEGER, INTENT(IN)  ::  river_rows
                       ! number of rows in river routing model
    INTEGER, INTENT(IN)  ::  river_row_length
                       ! number of E-W points for river model
    INTEGER, INTENT(IN)  ::  model_type
                       ! type (Global,LAM etc) of model
    INTEGER, INTENT(IN)  ::  nproc_ew
                       ! number of processors East-West
    INTEGER, INTENT(IN)  ::  nproc_ns
                       ! number of processors North-South
    INTEGER, INTENT(IN)  ::  extended_halo_ew
                       ! size of extended EW halo
    INTEGER, INTENT(IN)  ::  extended_halo_ns
                       ! size of extended NS halo
    INTEGER, INTENT(IN)  ::  nrima_max
                       ! size of rimwidth
    INTEGER, INTENT(IN)  ::  rimwidth(nrima_max)
                       ! size of blending region in the bcs    
    INTEGER, INTENT(OUT) ::  local_row_len
                       ! number of E-W points of this process
    INTEGER, INTENT(OUT) ::  local_n_rows      
                       ! number of rows of this process

    LOGICAL, OPTIONAL    :: noSmallHalos

! Local variables
    INTEGER                       :: smallHaloSize
    INTEGER                       :: iproc
    INTEGER                       :: iproc_x
    INTEGER                       :: iproc_y
    INTEGER                       :: ifld
    INTEGER                       :: ihalo
    INTEGER                       :: idim
    INTEGER                       :: ipt
    INTEGER                       :: irest
    INTEGER                       :: start
    INTEGER                       :: size1
    INTEGER                       :: size
    INTEGER                       :: prow_n
    INTEGER                       :: prow_s
    INTEGER                       :: info
    INTEGER                       :: in_atm_decomp
    INTEGER                       :: max_ns
    INTEGER                       :: max_ew
    INTEGER                       :: my_comm       ! current communicator

    LOGICAL                       :: at_north
    LOGICAL                       :: small_size


    INTEGER                       :: size_x(0:nproc_ew-1)
    INTEGER                       :: size_y(0:nproc_ns-1)
    INTEGER                       :: start_x(0:nproc_ew-1)
    INTEGER                       :: start_y(0:nproc_ns-1)

    ! For river routing
    INTEGER                       :: rsize_x(0:nproc_ew-1)
    INTEGER                       :: rsize_y(0:nproc_ns-1)
    INTEGER                       :: rstart_x(0:nproc_ew-1)
    INTEGER                       :: rstart_y(0:nproc_ns-1)

    INTEGER                       :: errorstatus   ! =0 normal exit 
                                                   ! >0 error exit
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    CHARACTER(LEN=256)            :: cmessage      ! Error message
    CHARACTER(LEN=*), PARAMETER   :: routinename='DECOMP_DB:DECOMPOSE'


    IF (lhook) CALL dr_hook(routinename,zhook_in,zhook_handle)

! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

    IF (mype  ==  0) THEN

      IF (deco /= decomp_smexe) THEN
        IF ((nproc_ew  /=  1) .AND. (MOD(nproc_ew,2)  /=  0)) THEN
        errorstatus=2
        WRITE(cmessage,'(A,I3,A)')                                    &
            'Cannot run with an odd (',nproc_ew,                      &
            ') number of processors in the East-West direction.'        
        CALL ereport(routinename,errorstatus,cmessage)
      END IF

      IF (global_row_len /= 1 .AND. &
          MOD(global_row_len,2)  /=  0) THEN
        errorstatus=3
        WRITE(cmessage,'(A,I3,A)')                                    &
            'Cannot run with an odd (',global_row_len,                &
            ') number of points in the East-West direction.'       
        CALL ereport(routinename,errorstatus,cmessage)
      END IF
    END IF

    IF (extended_halo_ew  >   max_halo_size) THEN
      errorstatus=4
      WRITE(cmessage,'(A,I2,A,I2)')                                        &
          'East-West extended halo size (',extended_halo_ew,               &
          ') is too large. The maximum permitted size is Max_Halo_Size=',  &
          max_halo_size            
      CALL ereport(routinename,errorstatus,cmessage)
    END IF

    IF (extended_halo_ns  >   max_halo_size) THEN
      errorstatus=4
      WRITE(cmessage,'(A,I2,A,I2)')                                        &
          'North-South extended halo size (',extended_halo_ns,             &
          ') is too large. The maximum permitted size is Max_Halo_Size=',  &
          max_halo_size      
      CALL ereport(routinename,errorstatus,cmessage)
      END IF
    END IF
    
    CALL Allocate_Decomposition(deco)
    
! ------------------------------------------------------------------

    smallHaloSize=1
    IF (PRESENT(noSmallHalos)) THEN
      IF (noSmallHalos) THEN
        smallHaloSize=0
      END IF
    END IF
        
    DecompDB(deco)%sb_model_domain=model_type
    DecompDB(deco)%halosize(1,halo_type_single)   = smallHaloSize
    DecompDB(deco)%halosize(2,halo_type_single)   = smallHaloSize
    DecompDB(deco)%halosize(3,halo_type_single)   = 0

    DecompDB(deco)%halosize(1,halo_type_extended) = &
        extended_halo_ew
    DecompDB(deco)%halosize(2,halo_type_extended) = &
        extended_halo_ns
    DecompDB(deco)%halosize(3,halo_type_extended) = 0


    DecompDB(deco)%halosize(1,halo_type_no_halo) = 0
    DecompDB(deco)%halosize(2,halo_type_no_halo) = 0
    DecompDB(deco)%halosize(3,halo_type_no_halo) = 0


! ------------------------------------------------------------------
! 1.0 Set up global data size
! ------------------------------------------------------------------

    DecompDB(deco)%glsize(1,fld_type_p) =            &
        global_row_len
    DecompDB(deco)%glsize(2,fld_type_p) =            &
        global_n_rows
    DecompDB(deco)%glsize(3,fld_type_p) =            &
        tot_levels

    DecompDB(deco)%glsize(1,fld_type_u) =            &
        global_row_len
    DecompDB(deco)%glsize(2,fld_type_u) =            &
        global_n_rows
    DecompDB(deco)%glsize(3,fld_type_u) =            &
        tot_levels

    DecompDB(deco)%glsize(1,fld_type_v) =            &
        global_row_len
    IF (model_type  /=  mt_bi_cyclic_lam ) THEN
      IF (l_vatpoles) THEN
!       V-AT-POLES/ENDGAME : assume that global_n_rows = number of u/p rows
        DecompDB(deco)%glsize(2,fld_type_v) =          &
            global_n_rows+1
      ELSE 
        DecompDB(deco)%glsize(2,fld_type_v) =          &
            global_n_rows-1
      END IF ! vatpoles
    ELSE
      DecompDB(deco)%glsize(2,fld_type_v) =          &
          global_n_rows
    END IF
    DecompDB(deco)%glsize(3,fld_type_v) =            &
        tot_levels

    DecompDB(deco)%glsize(1,fld_type_r) =            &
        river_row_length
    DecompDB(deco)%glsize(2,fld_type_r) =            &
        river_rows
    DecompDB(deco)%glsize(3,fld_type_r) =            &
        1

        ! ------------------------------------------------------------------
        ! 2.0 Calculate decomposition
        ! ------------------------------------------------------------------

    ! select processors to use for the data decomposition

    DecompDB(deco)%nproc=nproc_ew*nproc_ns
    DecompDB(deco)%first_comp_pe = 0
    DecompDB(deco)%last_comp_pe = DecompDB(deco)%nproc-1

    !     Set the grid size

    DecompDB(deco)%gridsize(1) = nproc_ew
    DecompDB(deco)%gridsize(2) = nproc_ns
    DecompDB(deco)%gridsize(3) = 1

! Work out the decomposition in the East-West direction. As far as
! possible each processor has the same local row length. However, if
! this is not possible, the extra points are distributed symetrically
! such that each processor has the same number of points as the
! processor on the opposite side of the globe.

    start=1
    size=global_row_len/nproc_ew  ! local data size on each processor
                                  ! assuming nproc_EW divides exactly
                                  ! into global_row_len.
    irest=global_row_len-(size*nproc_ew)
                                  ! If it doesn't divide exactly then
                                  ! irest contains the number of left
                                  ! over points that need to be
                                  ! allocated to processors

! Check the domains are big enough for the extended halos
    IF ((mype  ==  0) .AND. (size  <   extended_halo_ew)) THEN
      errorstatus=5
      WRITE(cmessage,'(A,I3,A,I3,A,I3,A)')                            &
          'Too many processors in the East-West direction (',nproc_ew,    &
          ') to support the extended halo size (',extended_halo_ew,       &
          '). Try running with ',(global_row_len/extended_halo_ew),       &
          ' processors.'  

      CALL ereport(routinename,errorstatus,cmessage)
    END IF
    
    
    DO iproc=1,nproc_ew
      start_x(iproc-1)=start

      IF (iproc  <=  nproc_ew/2) THEN
        ! processor in the first half of the row
        IF (iproc  <=  (irest/2)) THEN
          size_x(iproc-1)=size+1 ! gets one of the extra points
        ELSE
          size_x(iproc-1)=size   ! gets the standard row length
        END IF
      ELSE
        ! processor in the second half of the row
        IF (iproc-(nproc_ew/2)  <=  (irest/2)) THEN
          size_x(iproc-1)=size+1 ! gets one of the extra points
        ELSE
          size_x(iproc-1)=size   ! gets the standard row length
        END IF
      END IF

      start=start+size_x(iproc-1)

    END DO

    IF (deco/=decomp_smexe) THEN

! Work out the decomposition in the East-West direction for the
! river routing grid. As far as possible all processors are given
! the same number of points. If this is not possible the extra "n"
! points are given to the first "n" processors in the LPG.
      start=1
      size=river_row_length/nproc_ew  !local data size on each processor
                                    !assuming nproc_EW divides exactly
                                    !into global_row_len.
      irest=river_row_length-(size*nproc_ew)
                                    !If it doesn't divide exactly then
                                    !irest contains the number of left
                                    !over points that need to be
                                    !allocated to processors

      DO iproc=1,nproc_ew
        rstart_x(iproc-1)=start

        IF (iproc  <=  nproc_ew/2) THEN
          ! processor in the first half of the row
          IF (iproc  <=  (irest/2)) THEN
             rsize_x(iproc-1)=size+1 ! gets one of the extra points
          ELSE
             rsize_x(iproc-1)=size   ! gets the standard row length
          END IF
        ELSE
          ! processor in the second half of the row
          IF (iproc-(nproc_ew/2)  <=  (irest/2)) THEN
            rsize_x(iproc-1)=size+1 ! gets one of the extra points
          ELSE
           rsize_x(iproc-1)=size   ! gets the standard row length
          END IF
        END IF
     
        start=start+rsize_x(iproc-1)     
      END DO
    END IF

! Work out the decomposition in the North-South direction. As far as
! possible each processor has the same number of rows. However, if this
! is not possible, the extra rows are distributed thus:
! - an extra row is given to the Northern most processor
! - the remaining extra rows are distributed symetrically around the
!   equator, starting at the processor(s) closest to the equator.

    start=1
    size=global_n_rows/nproc_ns  ! local data size on each processor
                                 ! assuming nproc_NS divides exactly
                                 ! into global_n_rows
    irest=global_n_rows-(size*nproc_ns)
                                 ! If it doesn't divide exactly then
                                 ! irest contains the number of left
                                 ! over points that need to be
                                 ! allocated to processors

! Check the domains are big enough for the extended halos
    IF (mype  ==  0) THEN
      small_size=.FALSE.
      IF (irest  >  0) THEN
        IF (size  <  extended_halo_ns) small_size=.TRUE.
      ELSE
        IF (size  <=  extended_halo_ns) small_size=.TRUE.
      END IF    ! irest > 0
      IF ( small_size ) THEN
        errorstatus=5

        WRITE(cmessage,'(A,I3,A,I3,A,I3,A)')                            &
            'Too many processors in the North-South direction (',nproc_ns,  &
            ') to support the extended halo size (',extended_halo_ns,       &
            '). Try running with ',((global_n_rows-1)/extended_halo_ns),    &
            ' processors.'


        CALL ereport(routinename,errorstatus,cmessage)
      END IF ! small_size
    END IF ! mype == 0


    DO iproc=1,nproc_ns
      size_y(iproc-1)=size
    END DO

    IF (irest  >=  1) THEN
     IF (l_vatpoles) THEN
      ! give Southern most processors an extra row
      size_y(0)=size+1
      irest=irest-1
     ELSE
      ! give Northern most processors an extra row
      size_y(nproc_ns-1)=size+1
      irest=irest-1
     END IF ! vatpoles
    END IF

! Set up pointers to processor rows to which we will add extra rows
! to. These start around the equator, and will work out towards
! the poles.
    
    IF (MOD(nproc_ns,2)  ==  0) THEN  ! Even number of NS processors
      prow_s=nproc_ns/2
      prow_n=prow_s+1
    ELSE  ! Odd number of NS processors
      prow_s=(nproc_ns/2)+1
      prow_n=prow_s
    END IF

    DO WHILE (irest  >=  1)

      IF (prow_n  ==  prow_s) THEN
        size_y(prow_n-1)=size+1
        irest=irest-1
      ELSE
        size_y(prow_s-1)=size+1
        irest=irest-1
        IF (irest  >=  1) THEN
          size_y(prow_n-1)=size+1
          irest=irest-1
        END IF
      END IF

      prow_s=MAX(1,prow_s-1)
      prow_n=MIN(nproc_ns,prow_n+1)

    END DO

    DO iproc=1,nproc_ns
      start_y(iproc-1)=start
      start=start+size_y(iproc-1)
    END DO

! Work out the decomposition in the North-South direction for the
! river routing grid. As far as possible all processors are given
! the same number of rows. If this is not possible the extra "n" rows
! are given to the first "n" processors in the LPG.
    start=1
    size=river_rows/nproc_ns  ! local data size on each processor
                              ! assuming nproc_NS divides exactly
                              ! into global_row_len.
    irest=river_rows - (size*nproc_ns)
                              ! If it doesn't divide exactly then
                              ! irest contains the number of left
                              ! over points that need to be
                              ! allocated to processors

    DO iproc=0, nproc_ns-1
      rstart_y( iproc ) = start

      IF (iproc < irest) THEN
        rsize_y( iproc ) = size + 1
      ELSE
        rsize_y( iproc ) = size
      END IF

      start = start + rsize_y( iproc )
    END DO

! Set the local data shape and offsets of each processor

    DO iproc_y=0,nproc_ns-1

      IF (iproc_y  ==  (nproc_ns-1)) THEN
        at_north=.TRUE.  ! Doing the nothernmost processors
      ELSE
        at_north=.FALSE. ! Not doing the northernmost processors
      END IF

      DO iproc_x=0,nproc_ew-1

        iproc=DecompDB(deco)%first_comp_pe+         &
            iproc_x+(iproc_y*nproc_ew)

! Set the position in the Logical Processor Grid

        DecompDB(deco)%g_gridpos(1,iproc)=iproc_x
        DecompDB(deco)%g_gridpos(2,iproc)=iproc_y
        DecompDB(deco)%g_gridpos(3,iproc)=0

! Set the number of local datapoints (blsize) on the processor

!  Fields on P grid:
        DecompDB(deco)%g_blsize(1,fld_type_p,iproc)=                    &
            size_x(iproc_x)
        DecompDB(deco)%g_blsize(2,fld_type_p,iproc)=                    &
            size_y(iproc_y)
        DecompDB(deco)%g_blsize(3,fld_type_p,iproc)=                    &
            tot_levels

! Fields on U grid:
        DecompDB(deco)%g_blsize(1,fld_type_u,iproc)=                    &
            size_x(iproc_x)
        DecompDB(deco)%g_blsize(2,fld_type_u,iproc)=                    &
            size_y(iproc_y)
        DecompDB(deco)%g_blsize(3,fld_type_u,iproc)=                    &
            tot_levels

! Fields on V grid:
        DecompDB(deco)%g_blsize(1,fld_type_v,iproc)=                    &
            size_x(iproc_x)
        IF (at_north) THEN
          IF (l_vatpoles) THEN
            DecompDB(deco)%g_blsize(2,fld_type_v,iproc)=                  &
                size_y(iproc_y)+1
          ELSE
            DecompDB(deco)%g_blsize(2,fld_type_v,iproc)=                  &
                size_y(iproc_y)-1
          END IF ! vatpoles
          IF (model_type  ==  mt_bi_cyclic_lam ) THEN
            DecompDB(deco)%g_blsize(2,fld_type_v,iproc)=                &
                size_y(iproc_y)
          END IF
        ELSE
          DecompDB(deco)%g_blsize(2,fld_type_v,iproc)=                  &
              size_y(iproc_y)
        END IF
        DecompDB(deco)%g_blsize(3,fld_type_v,iproc)=                    &
            tot_levels
            
    IF (deco/=decomp_smexe) THEN
!  Fields on R (river) grid:
      DecompDB(deco)%g_blsize(1,fld_type_r,iproc)=                    &
          rsize_x(iproc_x)
      DecompDB(deco)%g_blsize(2,fld_type_r,iproc)=                    &
          rsize_y(iproc_y)
      DecompDB(deco)%g_blsize(3,fld_type_r,iproc)=1
    ELSE
!  Fields on R (river) grid smexe:
      DecompDB(deco)%g_blsize(1,fld_type_r,iproc)=0
      DecompDB(deco)%g_blsize(2,fld_type_r,iproc)=0
      DecompDB(deco)%g_blsize(3,fld_type_r,iproc)=0
    END IF

! Set the number of points including the halos on the processor

        DO ihalo=1,nhalo_max
          DO ifld=1,nfld_max
            DO idim=1,ndim_max
              DecompDB(deco)%g_lasize(idim,ifld,ihalo,iproc)=              &
                  DecompDB(deco)%g_blsize(idim,ifld,iproc)+            &
                  2*DecompDB(deco)%halosize(idim,ihalo)
            END DO  ! idim
          END DO  ! ifld
        END DO  ! ihalo
        
! Set the starting point in the global domain

        DecompDB(deco)%g_datastart(1,iproc)=         &
            start_x(iproc_x)
        DecompDB(deco)%g_datastart(2,iproc)=         &
            start_y(iproc_y)
        DecompDB(deco)%g_datastart(3,iproc)=1

        DecompDB(deco)%g_datastartr(1,iproc)=         &
             rstart_x(iproc_x)
        DecompDB(deco)%g_datastartr(2,iproc)=         &
             rstart_y(iproc_y)
        
        DecompDB(deco)%g_datastart_f(1,fld_type_p,iproc)=start_x(iproc_x)
        DecompDB(deco)%g_datastart_f(2,fld_type_p,iproc)=start_y(iproc_y)
        DecompDB(deco)%g_datastart_f(3,fld_type_p,iproc)=1
        
        DecompDB(deco)%g_datastart_f(1,fld_type_u,iproc)=start_x(iproc_x)
        DecompDB(deco)%g_datastart_f(2,fld_type_u,iproc)=start_y(iproc_y)
        DecompDB(deco)%g_datastart_f(3,fld_type_u,iproc)=1
        
        DecompDB(deco)%g_datastart_f(1,fld_type_v,iproc)=start_x(iproc_x)
        DecompDB(deco)%g_datastart_f(2,fld_type_v,iproc)=start_y(iproc_y)
        DecompDB(deco)%g_datastart_f(3,fld_type_v,iproc)=1

        
        IF (deco/=decomp_smexe) THEN
          DecompDB(deco)%g_datastart_f(1,fld_type_r,iproc)=rstart_x(iproc_x)
          DecompDB(deco)%g_datastart_f(2,fld_type_r,iproc)=rstart_y(iproc_y)
          DecompDB(deco)%g_datastart_f(3,fld_type_r,iproc)=1
        ELSE      
          DecompDB(deco)%g_datastart_f(1,fld_type_r,iproc)=0
          DecompDB(deco)%g_datastart_f(2,fld_type_r,iproc)=0
          DecompDB(deco)%g_datastart_f(3,fld_type_r,iproc)=1
        END IF

      END DO ! iproc_x
    END DO ! iproc_y
    
! Set up the pe_index_EW array - for each point along a global row
! it indicates the PE index (along the processor row) which
! contains that point

    DO iproc_x=0,nproc_ew-1
      DO ipt=DecompDB(deco)%g_datastart(1,iproc_x),  &
          DecompDB(deco)%g_datastart(1,iproc_x)+  &
          size_x(iproc_x)
        DecompDB(deco)%g_pe_index_ew(ipt)=iproc_x
      END DO
    END DO
    
! And fill in the halos at either end

    DO ipt=1-extended_halo_ew,0
      DecompDB(deco)%g_pe_index_ew(ipt)=0
    END DO
    
    DO ipt=global_row_len+1,global_row_len+extended_halo_ew
      DecompDB(deco)%g_pe_index_ew(ipt)=nproc_ew-1
    END DO
    
! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

    DO iproc_y=0,nproc_ns-1
      DO ipt=DecompDB(deco)%g_datastart(2,iproc_y*nproc_ew),            &
          DecompDB(deco)%g_datastart(2,iproc_y*nproc_ew)+            &
          size_y(iproc_y)
        DecompDB(deco)%g_pe_index_ns(ipt)=iproc_y
      END DO
    END DO
    
! And fill in the halos at either end

    DO ipt=1-extended_halo_ns,0
      DecompDB(deco)%g_pe_index_ns(ipt)=0
    END DO
    
    DO ipt=global_n_rows+1,global_n_rows+extended_halo_ns
      DecompDB(deco)%g_pe_index_ns(ipt)=nproc_ns-1
    END DO
    
! ------------------------------------------------------------------
! 3.0 Set boundary conditions
! ------------------------------------------------------------------
        
    IF (model_type  ==  mt_smexe) THEN
      ! These aren't really used
      DecompDB(deco)%bound(1) = bc_static
      DecompDB(deco)%bound(2) = bc_static
    ELSE IF (model_type  ==  mt_global) THEN
      DecompDB(deco)%bound(1) = bc_cyclic
      DecompDB(deco)%bound(2) = bc_overpole
    ELSE IF (model_type  ==  mt_lam) THEN
      DecompDB(deco)%bound(1) = bc_static
      DecompDB(deco)%bound(2) = bc_static
    ELSE IF (model_type  ==  mt_cyclic_lam) THEN
      DecompDB(deco)%bound(1) = bc_cyclic
      DecompDB(deco)%bound(2) = bc_static
    ELSE IF (model_type  ==  mt_bi_cyclic_lam) THEN
      DecompDB(deco)%bound(1) = bc_cyclic
      DecompDB(deco)%bound(2) = bc_cyclic
    ELSE
      errorstatus=10
      WRITE(cmessage,'(A,I3)')                                     &
          'Unrecognised model_type: ',model_type


      CALL ereport(routinename,errorstatus,cmessage)
    END IF
    DecompDB(deco)%bound(3) = bc_static

    CALL set_neighbour(deco)

! ------------------------------------------------------------------
! 4.0 Return the new data sizes 
! ------------------------------------------------------------------

! Set up the GCOM groups (effectively communicators in MPI):
    CALL gc_get_communicator(my_comm, info) 

! 1) Group of all processors on my row

    IF ( DecompDB(deco)%gridsize(2)  ==  1)          &
        THEN
      DecompDB(deco)%gc_proc_row_group=my_comm
    ELSE
      CALL gcg_split(mype,nproc_max,                                  &
          DecompDB(deco)%g_gridpos(2,mype),            &
          info,                                                         &
          DecompDB(deco)%gc_proc_row_group)
    END IF

! 2) Group of all processors on my column

    IF ( DecompDB(deco)%gridsize(1)  ==  1)          &
        THEN
      DecompDB(deco)%gc_proc_col_group=my_comm
    ELSE
      CALL gcg_split(mype,nproc_max,                                  &
          DecompDB(deco)%g_gridpos(1,mype),            &
          info,                                                         &
          DecompDB(deco)%gc_proc_col_group)
    END IF

! 3) Group of all processors in the atmosphere model
    IF (DecompDB(deco)%nproc  ==  nproc_max)        &
        THEN
      DecompDB(deco)%gc_all_proc_group=my_comm
    ELSE
      IF ((mype  >=  DecompDB(deco)%first_comp_pe)  &
          .AND.                                                         &
          (mype  <=  DecompDB(deco)%last_comp_pe) )  &
          THEN
        in_atm_decomp=1
      ELSE
        in_atm_decomp=0
      END IF

      CALL gcg_split(mype,nproc_max,in_atm_decomp,info,               &
          DecompDB(deco)%gc_all_proc_group)
    END IF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

    DecompDB(deco)%set=.TRUE.

! And return the new horizontal dimensions

    local_row_len=DecompDB(deco)%g_blsize(1,fld_type_p,mype)
    local_n_rows=DecompDB(deco)%g_blsize(2,fld_type_p,mype)
    
! ------------------------------------------------------------------
! 5.0 Check decomposition valid for LAM's           
! ------------------------------------------------------------------

!Check that the decomposition is valid for the mes.  The halos and the
!rimwidth must be contained within a single processor - they cannot
!cross over into another one as the indexing for the lbcs will go
!awry.
    
    IF (model_type  ==  mt_lam) THEN
      IF ((mype  ==  nproc_ew*(nproc_ns-1)) .OR.                      &
          (mype  ==  nproc_ew-1) ) THEN   !ie nw corner or se corner
        DO idim = 1, nrima_max
          IF(nproc_ns  == 1)THEN
            size=MAX(2*rimwidth(idim) +1,                               &
                extended_halo_ns + rimwidth(idim) +2 )
          ELSE
            size=extended_halo_ns + rimwidth(idim) +2
          END IF
          size1=extended_halo_ns + rimwidth(idim) +2
          IF(size  >    local_n_rows ) THEN
            max_ns = global_n_rows /size1
            errorstatus = 4
            WRITE(cmessage, '(A,A,I2)')                                 &
                'Too many processors in the North-South direction.',        &
                'The maximum permitted is ',max_ns

            CALL ereport(routinename,errorstatus,cmessage)
          END IF
        END DO
      END IF
      IF (mype  ==  nproc_ew*nproc_ns-1 .OR.                            &
          mype  ==  0 ) THEN        ! ie sw or ne corner
        DO idim = 1, nrima_max
          IF(nproc_ew  == 1)THEN
            size=MAX(2*rimwidth(idim) +1,                               &
                extended_halo_ew + rimwidth(idim) +2 )
          ELSE
            size=extended_halo_ew + rimwidth(idim) +2
          END IF
          size1=extended_halo_ew + rimwidth(idim) +2
          IF(size  >   local_row_len ) THEN
            max_ew = global_row_len / size1
            errorstatus = 4
            WRITE(cmessage, '(A,A,I2)')                                 &
                'Too many processors in the East-West direction.',      &
                'The maximum permitted is ',2* INT(max_ew/2)


            CALL ereport(routinename,errorstatus,cmessage)
          END IF
        END DO
      END IF
    END IF
    IF (lhook) CALL dr_hook(routinename,zhook_out,zhook_handle)
    RETURN

  END SUBROUTINE decompose_full

! Parallel UM: Sets up array of tids of adjacent processes.

! Subroutine interface:

  SUBROUTINE set_neighbour(decomp_type)

    IMPLICIT NONE

! Description:
! This routine finds the tids of the North, South, East and West
! neighbouring processes. It takes account of the boundary
! condition in each dimension (X and Y) which can be either:
! cyclic : wrap around
! overpole : data moves from processor on opposite side of pole
! static : no wrap around

! Method:
! The tid of each neighbouring process is calculated (taking into
! account the relevant boundary conditions) and placed in the
! neighbour array.

    INTEGER, INTENT(IN)           :: decomp_type  ! Decomposition type 
    ! to update neighbours  

    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook('SET_NEIGHBOUR',zhook_in,zhook_handle)

! Set Northen neighbour

    IF ( DecompDB(decomp_type)%g_gridpos(2,mype)  <        &
        (DecompDB(decomp_type)%gridsize(2)-1) ) THEN
      ! This processor is not at the North of the LPG

      DecompDB(decomp_type)%neighbour(pnorth) =            &
          mype + DecompDB(decomp_type)%gridsize(1)
      
    ELSE IF (DecompDB(decomp_type)%bound(2)  ==  bc_cyclic) THEN
      ! This processor at the top of LPG, and has cyclic BCs.

      DecompDB(decomp_type)%neighbour(pnorth) =                      &
          mype - DecompDB(decomp_type)%nproc +                       &
          DecompDB(decomp_type)%gridsize(1)
      
    ELSE IF (DecompDB(decomp_type)%bound(2)  ==  bc_overpole) THEN
      ! This processor passes data over the pole
      
      IF ((DecompDB(decomp_type)%g_gridpos(1,mype)+1)  <=            &
          (DecompDB(decomp_type)%gridsize(1)/2)) THEN
        DecompDB(decomp_type)%neighbour(pnorth) =                    &
            mype+DecompDB(decomp_type)%gridsize(1)/2
      ELSE
        DecompDB(decomp_type)%neighbour(pnorth) =                    &
            mype-DecompDB(decomp_type)%gridsize(1)/2
      END IF
      
    ELSE
      
      ! This processor at top of LPG and has static BCs
      DecompDB(decomp_type)%neighbour(pnorth) =nodomain
      
    END IF
    
    ! Set Southern neighbour
    
    IF ( DecompDB(decomp_type)%g_gridpos(2,mype)  >   0) THEN
      
      ! This processor is not at the South of the LPG
      DecompDB(decomp_type)%neighbour(psouth) =                      &
          mype - DecompDB(decomp_type)%gridsize(1)
      
    ELSE IF (DecompDB(decomp_type)%bound(2)  ==  bc_cyclic) THEN
      
      ! This processor at the bottom of LPG, and has cyclic BCs.
      DecompDB(decomp_type)%neighbour(psouth) =                      &
          mype + DecompDB(decomp_type)%nproc -                       &
          DecompDB(decomp_type)%gridsize(1)
      
    ELSE IF (DecompDB(decomp_type)%bound(2)  ==  bc_overpole) THEN
      
      ! This processor passes data over the pole
      IF ((DecompDB(decomp_type)%g_gridpos(1,mype)+1)  <=            &
          (DecompDB(decomp_type)%gridsize(1)/2)) THEN
        DecompDB(decomp_type)%neighbour(psouth) =                    &
            mype+DecompDB(decomp_type)%gridsize(1)/2
      ELSE
        DecompDB(decomp_type)%neighbour(psouth) =                    &
            mype-DecompDB(decomp_type)%gridsize(1)/2
      END IF
      
    ELSE
      
      ! This processor at top of LPG and has static BCs
      DecompDB(decomp_type)%neighbour(psouth) =nodomain
      
    END IF
    
    ! Set Western neighbour
    
    IF ( DecompDB(decomp_type)%g_gridpos(1,mype)  >   0) THEN
      
      ! This processor is not at the left of the LPG
      DecompDB(decomp_type)%neighbour(pwest) = mype - 1
      
    ELSE IF (DecompDB(decomp_type)%bound(1)  ==  bc_cyclic) THEN
      
      ! This processor at the left of the LPG, and has cyclic BCs.
      DecompDB(decomp_type)%neighbour(pwest) =                        &
          mype + DecompDB(decomp_type)%gridsize(1) - 1
      
    ELSE
      
      ! This processor at top of LPG and has static BCs
      DecompDB(decomp_type)%neighbour(pwest) =nodomain
      
    END IF
    
    ! Set Eastern neighbour
    IF ( DecompDB(decomp_type)%g_gridpos(1,mype)  <                   &
        (DecompDB(decomp_type)%gridsize(1)-1) ) THEN
      
      ! This processor is not at the right of the LPG
      DecompDB(decomp_type)%neighbour(peast) = mype + 1
      
    ELSE IF (DecompDB(decomp_type)%bound(1)  ==  bc_cyclic) THEN
      
      !       This processor at the left of the LPG, and has cyclic BCs.
      DecompDB(decomp_type)%neighbour(peast) =                        &
          mype - DecompDB(decomp_type)%gridsize(1) + 1
      
    ELSE
      
      !       This processor at top of LPG and has static BCs
      DecompDB(decomp_type)%neighbour(peast) =nodomain
      
    END IF
    
    IF (lhook) CALL dr_hook('SET_NEIGHBOUR',zhook_out,zhook_handle)
    
  END SUBROUTINE set_neighbour
  
END MODULE decomp_db
