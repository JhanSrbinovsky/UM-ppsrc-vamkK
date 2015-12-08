! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module: REGRID_ALLOC_CALC ------------------------------------------
!
!  Purpose: Location of routines which initialise the 
!  persistent data structures of regridding routines. This has been 
!  written in a semi-generic manner for regridding from a source grid to 
!  a target grid. But primarily for river routing regridding
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: MPP

MODULE REGRID_ALLOC_CALC

  USE mpl
  USE regrid_types
  
  IMPLICIT NONE     
   
CONTAINS

  ! This routine is provided to concern information required by other 
  ! processes and receive concerns required by this process
  SUBROUTINE SWAP_CONCERNS(send_concern, recv_concern, error, cmessage) 

    USE regrid_utils, ONLY: error_check_mpl
    IMPLICIT NONE
    
    TYPE(CONCERN), DIMENSION(:), INTENT(INOUT) :: send_concern
    ! this is not initialised apart from the proc data member 
    ! and size of grid point data 

    TYPE(CONCERN), DIMENSION(:), INTENT(INOUT) :: recv_concern
    ! this should be initialised with the concern information 
    ! needed for field data from this process's subdomain

    CHARACTER(len=*), INTENT(OUT) :: cmessage 
    ! if an error occurs this argument is set with an explanation 
    ! of the error 

    INTEGER, INTENT(OUT) :: error 
    ! set to -1 if an error occurs, 0 otherwise

    ! locals 

    INTEGER i, my_comm
    ! loop counter 

    INTEGER X_TAG, Y_TAG
    ! mpl tags 

    INTEGER s_concern_size, r_concern_size, request_size
    ! extent of send and recev concerns passed

    INTEGER, ALLOCATABLE :: requests(:), statuses(:,:), ierror(:)
    ! for receiving requests, mpl status, and mpl error codes 
    
    cmessage = ""
    error = 0
    X_TAG = 1
    Y_TAG = 2

    CALL gc_get_communicator(my_comm, error)

    s_concern_size = SIZE(send_concern,1)
    r_concern_size = SIZE(recv_concern,1)
    request_size = s_concern_size*2 + r_concern_size*2

    ! allocate requests for sends and recvs
    ALLOCATE(requests(request_size))
    ALLOCATE(ierror(request_size))

    requests = 0
    ierror = 0

    ! send lambda/phi point components telling appropriate 
    ! proc what grid points you need
    DO i=1, r_concern_size
       
      CALL MPL_ISEND(recv_concern(i)%x, recv_concern(i)%size,           &
           MPL_INTEGER, recv_concern(i)%proc_num, X_TAG,                &
           my_comm, requests(i), ierror(i))
      
      
      CALL MPL_ISEND(recv_concern(i)%y, recv_concern(i)%size,           &
           MPL_INTEGER, recv_concern(i)%proc_num, Y_TAG,                &
           my_comm, requests(i+r_concern_size), ierror(i))
      
    END DO

    ! recv lambda/phi point components which tells you
    ! what grid points other procs need
    DO i=1, s_concern_size
      
      CALL MPL_IRECV(send_concern(i)%x, send_concern(i)%size,           &
           MPL_INTEGER, send_concern(i)%proc_num, X_TAG,                &
           my_comm, requests(i+(r_concern_size*2)),                     &
           ierror(i+(r_concern_size*2)))
      
      CALL MPL_IRECV(send_concern(i)%y, send_concern(i)%size,           &
           MPL_INTEGER, send_concern(i)%proc_num, Y_TAG,                &
           my_comm, requests(i+(r_concern_size*2) +                     &
           s_concern_size), ierror(i+(r_concern_size*2))) 
      
    END DO
    
    ! check error for send and receives 
    CALL ERROR_CHECK_MPL(ierror, request_size, error)
    
    ALLOCATE(statuses(MPL_STATUS_SIZE, request_size))
    statuses = 0
    
    CALL MPL_WAITALL(request_size, requests, statuses, i) 

    ! account for case 
    IF(request_size /= 0) THEN
      ierror(1) = i
    ! check error for waitall 
      CALL ERROR_CHECK_MPL(ierror, 1, i)
    END IF
    
    IF(i /= 0) error = 1
    
    ! release request mem resource
    DEALLOCATE(requests)
    DEALLOCATE(statuses)
    DEALLOCATE(ierror)

    RETURN
  END SUBROUTINE SWAP_CONCERNS


  ! convenience function which calculates proc send and 
  ! recv list, the information concerns are based on global_index 
  ! array passed in
  SUBROUTINE CALC_PROC_LIST(global_index, index_size, send_list,          &
        recv_list, send_size, recv_size, p_count_send, p_count_recv,      &
        g_row_length)

    USE UM_ParVars, ONLY: g_pe_index_EW, g_pe_index_NS, nproc, nproc_x
    USE UM_ParCore, ONLY: mype
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: global_index(index_size), index_size,          &
         g_row_length 
    ! encodes the 2d grid points needed by this process, which can be 
    ! decoded with the grid points global row length.

    INTEGER, INTENT(OUT) :: send_list(nproc), recv_list(nproc)
    ! send_list, a list of processes that will send regridding 
    ! information to this process and vice versa for recv_list
    
    INTEGER, INTENT(OUT) :: recv_size(nproc), send_size(nproc)
    ! the number of grid points that processes list in send and recv list
    ! will be sending receiving respectively

    INTEGER, INTENT(OUT) :: p_count_send, p_count_recv
    ! actual number processes in send_list and recv_list (nproc is just
    ! the max amount of send and recv processes)

    ! locals     
    INTEGER n_grid_proc_list(nproc), n_grid_proc_list_all(nproc*nproc)
    INTEGER pe_col, pe_row, proc, ierror
    INTEGER xpt, ypt, i, j, my_comm

    send_list = 0
    recv_list = 0
    recv_size = 0
    send_size = 0
    ierror = 0
    n_grid_proc_list = 0
    n_grid_proc_list_all = 0

    DO i=1, index_size 
      xpt = MOD(global_index(i)-1, g_row_length)+1
      ypt = (global_index(i)-1)/g_row_length+1
      
      pe_col = g_pe_index_EW(xpt)
      pe_row = g_pe_index_NS(ypt)
      proc = pe_row*nproc_x + pe_col 

      
      ! record number of grid pts process needs to send       
      IF(proc /= mype) THEN
        




        n_grid_proc_list(proc+1) = n_grid_proc_list(proc+1)             &
             + 1
      END IF
    END DO

    CALL GC_GET_COMMUNICATOR(my_comm, ierror)
    
    CALL MPL_ALLGATHER(n_grid_proc_list, nproc, MPL_INTEGER,            &
         n_grid_proc_list_all, nproc, MPL_INTEGER,                      &
         my_comm, ierror)
    
    
    p_count_send = 0 
    p_count_recv = 0
    
    DO i=0, nproc-1
      DO j=0, nproc-1
        !! look at own list and determine procs you will be 
        !! receving from and size 
        IF(i==mype) THEN 
          
          IF(n_grid_proc_list_all(i*nproc+j+1)>0) THEN
            !! number of recv concerns
            p_count_recv = p_count_recv + 1
            !! size for proc
            recv_size(p_count_recv) = n_grid_proc_list_all(         &
                 i*nproc+j+1)
            !! proc to recv from
            recv_list(p_count_recv) = j
            
          END IF
          
        ELSE 
          !! look at other procs list and determine who wants
          !! you to send to them
          
          IF(j==mype) THEN
            IF(n_grid_proc_list_all(i*nproc+j+1)>0) THEN
              !! number of send concerns
              p_count_send = p_count_send + 1
              !! size of data to send proc
              send_size(p_count_send) = n_grid_proc_list_all(      &
                   i*nproc+j+1)
              !! proc to send to
              send_list(p_count_send) = i
            END IF
          END IF
        END IF
        
      END DO
    END DO
    
  END SUBROUTINE CALC_PROC_LIST
  
  
  
  ! Calculate and allocates space for send and receive concerns 
  ! in order to regrid from one grid to the other
  SUBROUTINE CALC_ALLOC_CONCERNS(send_concern, recv_concern,            &
       index_size, global_index, g_row_length, g_rows, grid, error,     &
       cmessage)
    
    USE UM_ParVars, ONLY: nproc, g_pe_index_EW, g_pe_index_NS, nproc_x
    USE UM_ParCore, ONLY: mype 
    USE ereport_mod, ONLY: ereport 
    IMPLICIT NONE

    TYPE(CONCERN), POINTER :: send_concern(:)
    TYPE(CONCERN), POINTER :: recv_concern(:)
    ! these allocated within this routine and set with coordination 
    ! information needed to regrid this processes subdomain from 
    ! src to target grid

    INTEGER, INTENT(IN) :: index_size, global_index(index_size)
    ! lat and long start grid location in global indexing
    ! contains information on the grid to calculate 'regrid to'
    ! information

    INTEGER, INTENT(IN) :: grid    
    ! grid one is regridding from 

    INTEGER, INTENT(IN) :: g_row_length, g_rows 
    ! this the row length and rows of the src grid 
    ! across all processes 

    CHARACTER(len=*), INTENT(OUT) :: cmessage
    ! set if error occurs, otherwise empty

    INTEGER, INTENT(OUT) :: error
    ! set to -1 if an error occurs otherwise 0 

    ! local variables

    INTEGER num_procs_send, num_procs_recv
    ! number of procs for which this processor is a send concern 
    ! and receive concern respectively
    
    INTEGER send_list(nproc), recv_list(nproc)
    ! a list to possibly contain all processes you need to send to 
    ! or recv from 

    INTEGER recv_size(nproc), send_size(nproc)
    ! number of grid points each process in list will be sending 
    ! to help allocate space on respective processes


    INTEGER p_count_send, p_count_recv, pe_col, pe_row
    LOGICAL found_proc
    INTEGER lam, phi, proc, ierror, i, j, k, procCount
    INTEGER :: xpt, ypt
    

    error = 0
    cmessage = ""
    p_count_recv = 0
    p_count_send = 0

    ! retrieve send and recv list as well as sizes 
    CALL CALC_PROC_LIST(global_index, index_size, send_list, recv_list, &
         send_size, recv_size, p_count_send, p_count_recv, g_row_length)

    ! allocate resources for recv and send concerns

    ALLOCATE(recv_concern(p_count_recv))
    ALLOCATE(send_concern(p_count_send)) 

    ! allocate for recv and send concern space needed
    ! for send and recv
    
    DO i=1, p_count_recv
      
      ALLOCATE(recv_concern(i)%x(recv_size(i)))
      ALLOCATE(recv_concern(i)%y(recv_size(i)))
      ALLOCATE(recv_concern(i)%field(recv_size(i)))
      recv_concern(i)%proc_num = recv_list(i) 
      
      !! init size to zero in order to
      !! use for grid initialisation
      recv_concern(i)%size = 0 
      
    END DO
    
    
    DO i=1, p_count_send
      send_concern(i)%proc_num = send_list(i) 
      send_concern(i)%size = send_size(i)
      ALLOCATE(send_concern(i)%x(send_size(i)))
      ALLOCATE(send_concern(i)%y(send_size(i)))
      ALLOCATE(send_concern(i)%field(send_size(i)))
    END DO
    
    ! populate recv_concern with the data points you want
    
    DO i=1, index_size
      xpt = MOD(global_index(i)-1, g_row_length)+1
      ypt = (global_index(i)-1)/g_row_length+1
      
      pe_col = g_pe_index_EW(xpt)
      pe_row = g_pe_index_NS(ypt)
      proc = pe_row*nproc_x + pe_col 
      
      !! find relevant proc number and store grid point
      !! location, if relevant proc can't be found this 
      !! is a fatal error
      IF(proc /= mype) THEN
        
        found_proc = .FALSE.
        DO j=1, p_count_recv
          
          IF(proc==(recv_concern(j)%proc_num)) THEN
            
            recv_concern(j)%size = recv_concern(j)%size + 1
            
            !! update to next grid point
            k = recv_concern(j)%size
            recv_concern(j)%x(k) = xpt
            recv_concern(j)%y(k) = ypt
            found_proc = .TRUE.
            
          END IF
        END DO
        
        IF(.NOT. found_proc) THEN
          
          PRINT*, "proc, xpt, ypt, :", proc, xpt, ypt
          CALL FLUSH(6)
          
          ! proc has not been found in receive list
          CALL EREPORT("CALC_ALLOC_CONCERNS", 1,                        &
               "Concerns not calclated successfully")
          
        END IF
      END IF
      
    END DO
    
    ! swap concerns with other processors
    CALL SWAP_CONCERNS(send_concern, recv_concern, error, cmessage) 
    
    IF (error /= 0) THEN
      CALL EREPORT("CALC_ALLOC_CONCERNS", 1, cmessage)
    END IF
    
    RETURN
  END SUBROUTINE CALC_ALLOC_CONCERNS
  
  
  ! calculates and allocates space for regridding via map max method 
  SUBROUTINE CALC_ALLOC_MAX_CONCERNS(send_concern_max, recv_concern_max,&
       contribution, base_tr, count_tr, row_length, rows,               &
       weight, src_grid, global_index, index_size, global_row_length,   &
       global_rows, local_targ_row_length, local_targ_rows, targ_grid,  &
       trip_mask, want, global_src_row_length, error, cmessage)
    
    USE UM_ParVars , ONLY: nproc, g_pe_index_EW, g_pe_index_NS, nproc_x,&
         nproc_y
    USE UM_ParCore , ONLY: mype 
    USE regrid_utils, ONLY: local_to_global_gridpt, get_proc_for_gridpt,&
         sort_contributors
    USE ereport_mod, ONLY: ereport 
    IMPLICIT NONE

    TYPE(CONCERN_MAX), POINTER :: send_concern_max(:)
    TYPE(CONCERN_MAX), POINTER :: recv_concern_max(:)

    ! the send and recv concern max to allocated and initialised with 
    ! coordination information need for map_max regridding

    TYPE(CONTRIBUTION_INFO), POINTER :: contribution(:)
    ! more convenient form of recv_concern_max

    INTEGER, INTENT(IN) :: rows, row_length
    ! the local row and row_length of the target grid

    INTEGER, INTENT(IN) :: global_index(index_size), index_size 
    ! stores indexes of contributing target pixels to src 

    INTEGER, INTENT(IN) :: global_src_row_length
    ! global row length of src grid

    REAL, INTENT(IN) :: weight(index_size)
    ! contirbuting weight of src field points pointed to by global_index 

    LOGICAL, INTENT(IN) :: trip_mask(row_length, rows)
    ! trip mask used to indicate wanted point
    
    LOGICAL, INTENT(IN) :: want
    ! mask indicator to note wanted points

    INTEGER, INTENT(IN) :: base_tr(row_length, rows)
    ! starting index of contributors to target field pixel (x,y) in 
    ! global_index 
    INTEGER, INTENT(IN) :: count_tr(row_length, rows) 
    ! the number of targets that contributes in for a given pixel
    ! hence, global_index(base_tr(x,y):base_tr(x,y)+count_tr(x,y)) contribute 
    ! to target field pixel x,y

    INTEGER, INTENT(IN) ::  global_rows, global_row_length
    ! row and row length of the targ grid across processes 

    INTEGER, INTENT(IN) :: src_grid, targ_grid
    ! src grid, the grid being regridded from and target grid; the grid
    ! being regridded to

    INTEGER, INTENT(IN) :: local_targ_row_length, local_targ_rows
    ! number of rows, row length of the target grid on this process

    INTEGER, INTENT(OUT) :: error 
    ! set to -1 if an error occurs and 0 otherwise

    CHARACTER(len=*), INTENT(OUT) :: cmessage
    ! set if an error occurs explanating error
    
    ! locals

    INTEGER send_list(nproc), recv_list(nproc)
    ! a list to possibly contain all processes you need to send to 
    ! or recv from 

    INTEGER recv_size_a(nproc), send_size_a(nproc)
    ! number of grid points each process in list will be sending 
    ! to help allocate space on respective processes

    INTEGER p_count_send, p_count_recv
    ! actual size of process and send/recv size lists
    
    INTEGER i, j, k, l, m, ip, targx, targy, x, y, contrib_size
    INTEGER sc_size, dat_size, recv_size, rc_size
    REAL weight_max
    INTEGER proc, max_ip 
    LOGICAL found_proc


    error = 0
    cmessage = ""
    
    ! calculate who you would list of processors you would send to and 
    ! receive from if you were going to regrid via area average
    CALL CALC_PROC_LIST(global_index, index_size, send_list, recv_list, &
         send_size_a, recv_size_a, p_count_send, p_count_recv,          &
         global_row_length)
    
    ! for the processors you would receive from this 
    ! represent src grid which have contributed to 
    ! target points in this subdomain 
    recv_size = p_count_recv 
    ALLOCATE(send_concern_max(recv_size))
    
    DO i=1, recv_size 
      dat_size = recv_size_a(i) 
      
      send_concern_max(i)%proc_num = recv_list(i)
      ALLOCATE(send_concern_max(i)%x(dat_size))
      ALLOCATE(send_concern_max(i)%y(dat_size))
      
      ALLOCATE(send_concern_max(i)%xtarg(dat_size))
      ALLOCATE(send_concern_max(i)%ytarg(dat_size))
      
      ALLOCATE(send_concern_max(i)%weight(dat_size))
      ALLOCATE(send_concern_max(i)%field(dat_size))
      ALLOCATE(send_concern_max(i)%contribute(dat_size))
      
      send_concern_max(i)%size = 0
    END DO
    
    
    ! for each src point that exist in this process's subdomain 
    ! determine which target point outwith this subdomain it contributes
    ! to and store this in send_concern_max 
    sc_size = recv_size
    
    DO j=1, rows  
      DO i=1, row_length  


          IF(count_tr(i,j) /= 0) THEN
          
            weight_max = 0.
            
            ! determine which src pixel contributes to this
            ! this target the most
            DO k=1, count_tr(i,j) 
              ip = base_tr(i,j) + k
              IF(weight(ip) > weight_max) THEN
                weight_max = weight(ip)
                max_ip = ip 
              END IF
            END DO
            
            DO k=1, count_tr(i,j) 
              
              ip = base_tr(i,j) + k
              x = MOD(global_index(ip)-1, global_row_length)+1
              y = (global_index(ip)-1)/global_row_length+1
              
              proc =  GET_PROC_FOR_GRIDPT(x, y, targ_grid)
              
              ! consider all processes except calling process
              IF(proc /= mype) THEN
                
                found_proc = .FALSE.
                ! look for process this src point contributes 
                ! to its subdomain
                DO l=1, sc_size 
                  IF(proc==send_concern_max(l)%proc_num) THEN 
                    send_concern_max(l)%size = send_concern_max(l)%size + 1
                    m = send_concern_max(l)%size
                    
                    ! store src grid point (global indexing)  
                    ! this target contributes to
                    
                    send_concern_max(l)%xtarg(m) = x
                    send_concern_max(l)%ytarg(m) = y
                    
                    ! the target pixel currently being considered
                    targx = i
                    targy = j
                    
                    ! change local grid point to a 
                    ! global one
                    CALL LOCAL_TO_GLOBAL_GRIDPT(targx, targy, src_grid)
                    
                    ! record the global 2d indices the contributing 
                    ! target
                    send_concern_max(l)%x(m) = targx 
                    send_concern_max(l)%y(m) = targy
                  
                    ! store target's weight contribution to src
                    send_concern_max(l)%weight(m) = weight(ip) 
                    
                    ! indicate if target actually contributed 
                    ! to the most by this source 
                    ! (max contributor based regridder)
                    send_concern_max(l)%contribute(m) = (ip == max_ip) &
                         .AND. (trip_mask(i,j) .EQV. want)      
                    
                    found_proc = .TRUE.
                    
                  END IF ! proc = concern%proc_num
                END DO ! 1, sc_size
                
                !! proc should be found if not fatal error!
                IF(.NOT. found_proc) THEN
                  print*, "Fatal error! proc not found: ", mype
                END IF
                
              END IF ! proc /= mype
            END DO ! 1 to cout_targ
          END IF ! count_tr_targ /= 0
!        END IF ! trip_mask .EQV. want
      END DO ! i
    END DO ! j
      
    
    ! the data needed by this process equally reflect the send
    ! concerns needed for area average regridding. all target points 
    ! which contributed to determining current source will now have 
    ! the source value and their weights returned to them in 
    ! recv concern max 

    ! recv_concern_max resource should mirror send_concern
    recv_size = p_count_send 
    
    ALLOCATE(recv_concern_max(recv_size)) 

    DO i=1, recv_size 
      dat_size = send_size_a(i)
      recv_concern_max(i)%size = dat_size
      recv_concern_max(i)%proc_num = send_list(i)
      
      ALLOCATE(recv_concern_max(i)%x(dat_size))
      ALLOCATE(recv_concern_max(i)%y(dat_size))
      
      ALLOCATE(recv_concern_max(i)%xtarg(dat_size))
      ALLOCATE(recv_concern_max(i)%ytarg(dat_size))
      
      ALLOCATE(recv_concern_max(i)%weight(dat_size))
      ALLOCATE(recv_concern_max(i)%field(dat_size))
      ALLOCATE(recv_concern_max(i)%contribute(dat_size))
      
    END DO
    
    ! swap concern data that does not change per regrid 
    ! iteration for recv concern
    CALL SWAP_CONCERNS_MAX(send_concern_max, recv_concern_max, error,   &
         cmessage)
    
    IF (error /= 0) CALL EREPORT("CALC_ALLOC_MAX_CONCERNS", -1,         &
         "swap concern max failed")
    
    contrib_size = local_targ_rows*local_targ_row_length
    recv_size = p_count_send
    
    ! allocate space to contain all local target field points
    ALLOCATE(contribution(local_targ_rows*local_targ_row_length))
    
    
    ! calculate and intitalise contribution array 
    CALL CALC_CONTRIBUTORS(base_tr, count_tr, weight, row_length, rows, &
         global_index, index_size, global_row_length, targ_grid,        &
         src_grid, local_targ_row_length, recv_concern_max, recv_size,  &
         contribution, contrib_size, error)

    IF (error /= 0) CALL EREPORT("CALC_ALLOC_MAX_CONCERNS", -1,         &
         "Failed calc contribtutions")

    CALL SORT_CONTRIBUTORS(contribution, contrib_size,                  &
         global_src_row_length, error)

    IF (error /= 0) CALL EREPORT("CALC_ALLOC_MAX_CONCERNS", -1,         &
         "Failed sort_contributions")

    RETURN
  END SUBROUTINE CALC_ALLOC_MAX_CONCERNS
  

 
  ! the purpose of this routine is to recv grid points so 
  ! this processor know what field values it needs to send
  ! to other processors
  SUBROUTINE SWAP_CONCERNS_MAX(send_concern_max, recv_concern_max,      &
       error, cmessage)
    
    USE regrid_utils, ONLY: error_check_mpl
    IMPLICIT NONE

    TYPE(CONCERN_MAX), DIMENSION(:), INTENT(INOUT) :: send_concern_max
    ! contains what this process needs to send to allow other procs
    ! to regrid their field 
    ! the size of data and procs to send to, must be set before entry
    ! to this routine

    TYPE(CONCERN_MAX), DIMENSION(:), INTENT(INOUT) :: recv_concern_max 
    ! set to info needed by other procs to regrid their subdomain
    ! the size of data and procs to recv from must be set before entry
    ! to this routine
    
    INTEGER, INTENT(OUT) :: error 
    ! set to -1 if an error occurs, otherwise 0

    CHARACTER(len=*), INTENT(OUT) :: cmessage 
    ! set with info explaining error should one occurs, otherise empty

    ! locals 

    INTEGER i, my_comm
    ! loop index
    INTEGER X_TAG, Y_TAG, X_SRC_TAG, Y_SRC_TAG, W_TAG, C_TAG
    ! mpl tags to make sure procs receive the right information
    INTEGER s_concern_size, r_concern_size, request_size, req_offset,   &
         req_pt, error_i(1)
    ! extent of send and recev concerns passed

    INTEGER, ALLOCATABLE :: requests(:), statuses(:,:), ierror(:)
    
    error = 0 
    cmessage = ""
    CALL gc_get_communicator(my_comm, error)
    s_concern_size = SIZE(send_concern_max,1)
    r_concern_size = SIZE(recv_concern_max,1)
    request_size = s_concern_size*6 + r_concern_size*6

    ! allocate requests for sends and recvs
    ALLOCATE(requests(request_size))
    ALLOCATE(ierror(request_size))

    requests = 0
    req_pt = 0

    W_TAG = 1
    X_TAG = 2
    Y_TAG = 3
    X_SRC_TAG = 4
    Y_SRC_TAG = 5
    C_TAG = 6
    

    ! send lambda/phi gird point which tell what target maps 
    ! to what source and the weighting to apply to the associated
    ! field point value
    DO i=1, s_concern_size

      req_pt =  (i-1)*6+1
      CALL MPL_ISEND(send_concern_max(i)%weight,                        &
           send_concern_max(i)%size, MPL_REAL,                          &
           send_concern_max(i)%proc_num, W_TAG, my_comm,                &
           requests(req_pt), ierror(req_pt))
      
      req_pt =  (i-1)*6+2
      CALL MPL_ISEND(send_concern_max(i)%x,                             &
           send_concern_max(i)%size, MPL_INTEGER,                       &
           send_concern_max(i)%proc_num, X_TAG, my_comm,                &
           requests(req_pt), ierror(req_pt))
      
      req_pt =  (i-1)*6+3
      CALL MPL_ISEND(send_concern_max(i)%y,                             &
           send_concern_max(i)%size, MPL_INTEGER,                       &
           send_concern_max(i)%proc_num, Y_TAG, my_comm,                &  
           requests(req_pt), ierror(req_pt))
      
      req_pt =  (i-1)*6+4
      CALL MPL_ISEND(send_concern_max(i)%xtarg,                         &
           send_concern_max(i)%size, MPL_INTEGER,                       &
           send_concern_max(i)%proc_num, X_SRC_TAG, my_comm,            &
           requests(req_pt), ierror(req_pt))
      
      req_pt =  (i-1)*6+5
      CALL MPL_ISEND(send_concern_max(i)%ytarg,                         &
           send_concern_max(i)%size, MPL_INTEGER,                       &
           send_concern_max(i)%proc_num, Y_SRC_TAG, my_comm,            &
           requests(req_pt), ierror(req_pt))
      
      req_pt =  (i-1)*6+6
      CALL MPL_ISEND(send_concern_max(i)%contribute,                    &
           send_concern_max(i)%size, MPL_LOGICAL,                       &
           send_concern_max(i)%proc_num, C_TAG, my_comm,                &
           requests(req_pt), ierror(req_pt))
      
    END DO
    
    req_offset = s_concern_size*6

    ! recv lambda/phi gird point which tell what target maps 
    ! to what source and the weighting to apply to the associated
    ! field point value
    DO i=1, r_concern_size 
      
      req_pt = (i-1)*6+1 + req_offset
      CALL MPL_IRECV(recv_concern_max(i)%weight,                        &
           recv_concern_max(i)%size, MPL_REAL,                          &
           recv_concern_max(i)%proc_num, W_TAG, my_comm,                &   
           requests(req_pt), ierror(req_pt))
      
      req_pt = (i-1)*6+2 + req_offset
      CALL MPL_IRECV(recv_concern_max(i)%x,                             &  
           recv_concern_max(i)%size, MPL_INTEGER,                       &
           recv_concern_max(i)%proc_num, X_TAG, my_comm,                &
           requests(req_pt), ierror(req_pt))
      
      req_pt = (i-1)*6+3 + req_offset
      CALL MPL_IRECV(recv_concern_max(i)%y,                             &
           recv_concern_max(i)%size, MPL_INTEGER,                       &
           recv_concern_max(i)%proc_num, Y_TAG, my_comm,                &
           requests(req_pt), ierror(req_pt))
      
      req_pt = (i-1)*6+4 + req_offset
      CALL MPL_IRECV(recv_concern_max(i)%xtarg,                         &
           recv_concern_max(i)%size, MPL_INTEGER,                       &
           recv_concern_max(i)%proc_num, X_SRC_TAG, my_comm,            &
           requests(req_pt), ierror(req_pt))
      
      req_pt = (i-1)*6+5 + req_offset
      CALL MPL_IRECV(recv_concern_max(i)%ytarg,                         &
           recv_concern_max(i)%size, MPL_INTEGER,                       &
           recv_concern_max(i)%proc_num, Y_SRC_TAG, my_comm,            &
           requests(req_pt), ierror(req_pt))
      
      req_pt = (i-1)*6+6 + req_offset
      CALL MPL_IRECV(recv_concern_max(i)%contribute,                    &
           recv_concern_max(i)%size, MPL_LOGICAL,                       &
           recv_concern_max(i)%proc_num, C_TAG, my_comm,                &
           requests(req_pt), ierror(req_pt))
      
    END DO

    ! check mpl errors of recv and send
    CALL ERROR_CHECK_MPL(ierror, request_size,     &
         error)
    
    ALLOCATE(statuses(MPL_STATUS_SIZE, request_size))
    statuses = 0
    
    CALL MPL_WAITALL(request_size, requests, statuses, i) 
    
    ! check mpl errors of wait all

    IF(request_size /= 0) THEN
      error_i(1) = i
      CALL ERROR_CHECK_MPL(error_i, 1, error)
    END IF

    ! release request mem resource
    DEALLOCATE(requests)
    DEALLOCATE(statuses)
    DEALLOCATE(ierror)

    RETURN
  END SUBROUTINE SWAP_CONCERNS_MAX
  
  ! The routine gathers up target grid points which 
  ! contribute to this process's subdomain source
  ! points and also the target point's weight 
  SUBROUTINE CALC_CONTRIBUTORS(base_tr, count_tr, weight, row_length,   &
       rows, g_index_src, lenl, g_targ_row_length, targ_grid, src_grid, &
       l_targ_row_length, recv_concern_max, recv_size, contribution,    &
       contrib_size, error)

    USE regrid_utils, ONLY: local_to_global_gridpt,                     &
         get_proc_for_gridpt, global_to_local_gridpt,                   &
         gridpt_outside_proc_domain
    USE UM_ParVars, ONLY: nproc
    USE UM_ParCore, ONLY: mype
    USE ereport_mod, ONLY: ereport
    IMPLICIT NONE
 
    INTEGER, INTENT(IN) :: contrib_size, recv_size
    ! the size of contribution and recv_concern_max arrays respectively
 
    TYPE(CONCERN_MAX), DIMENSION(recv_size), INTENT(IN) :: recv_concern_max 
    ! contains information on src grid points this process needs to receive 
    ! to regrid 

    TYPE(CONTRIBUTION_INFO), DIMENSION(contrib_size), INTENT(INOUT) ::  &
         contribution
    ! this is set here using recv_concern_max, the contribution array is needed
    ! to compress and simplify information in recv concern max (to a more 
    ! usable form,  which for instance allows preserving reproducability 
    ! despite different decompositions). 

    INTEGER, INTENT(IN) :: g_index_src(lenl), lenl
    ! stores indexes of contributing src pixels to target

    REAL, INTENT(IN) :: weight(lenl)
    ! contirbuting weight of src field points pointed to by global_index    

    INTEGER, INTENT(IN) :: l_targ_row_length
    ! the local target row length

    INTEGER, INTENT(IN) :: rows, row_length
    ! the local src row and row lengt

    INTEGER, INTENT(IN) :: base_tr(row_length, rows)
    ! the starting point ("base") of contributing pixels in 
    ! global index

    INTEGER, INTENT(IN) :: count_tr(row_length, rows)
    ! the number of src points that contributes in for a given pixel
    ! hence, elements of global_index(base_tr(x,y):base_tr(x,y)+
    ! count_tr(x,y)) contribute to src field pixel x,y

    INTEGER, INTENT(IN) :: src_grid, targ_grid, g_targ_row_length
    ! the src grid we're regridding from and target grid we're regrdding 
    ! to (e.g. atmos, river)

    INTEGER, INTENT(OUT) :: error 
    ! set to -1 if error occurs, 0 otherwise
    
    ! local variables

    INTEGER i, j, k, ip 
    ! loop counter variables 

    INTEGER xpt, ypt, y_max, x_max
    ! grid point references 

    INTEGER index, dat_size
    ! used to track 1-D global src grid point index 

    REAL weight_max
    ! tracks the max weight contributed to a src field point by the target

    INTEGER, ALLOCATABLE :: incr(:)
    ! allows dynamic incrementing of index sizes 

    ALLOCATE(incr(contrib_size))

    DO i=1, contrib_size 
      contribution(i)%size = 0
    END DO

    ! a bit of a hack but fortran doesn't have native linked lists, so...
    ! so first we retrieve the space required of contributors for the target
    ! field subdomain which can then be used to allocate the space required 
    ! to store contribution information
    
    ! firstly account for space required by src contributors within 
    ! prcoess's subdomain
    DO j=1, rows
      DO i=1, row_length 
        
        weight_max = 0
        
        IF(count_tr(i,j) /= 0) THEN 
          
          ! get indices for target which src contributes 
          ! to the most
          DO k=1, count_tr(i,j)
            ip = base_tr(i,j) + k 
            IF(weight(ip) > weight_max) THEN
              weight_max = weight(ip) 
              x_max = MOD(g_index_src(ip)-1,g_targ_row_length)+1
              y_max = (g_index_src(ip)-1)/g_targ_row_length+1
            END IF
          END DO
          
          ! do not consider points contributed by src 
          ! outside proc domain
          IF(.NOT. GRIDPT_OUTSIDE_PROC_DOMAIN(x_max, y_max, targ_grid)) &
               THEN 
            
            xpt = x_max 
            ypt = y_max 
            CALL GLOBAL_TO_LOCAL_GRIDPT(xpt, ypt, targ_grid)
            
            index = ((ypt-1)*l_targ_row_length)+xpt
            
            !! increment to make space for this 
            !! contributing src
            contribution(index)%size = contribution(index)%size + 1
          END IF  ! gridpt_outside_proc_domain
        END IF ! count_tr /= 0
        
      END DO
    END DO
    
    
    ! now account for space needed by contribution 
    ! by src points outside this process's subdomain
    ! which contribute to this process's target 
    ! subdomain
    DO i=1, recv_size 
      
      dat_size = recv_concern_max(i)%size  
      DO j=1, dat_size 
        
        !! account src that actually contributes
        !! i.e. it is the maximum contributor
        IF(recv_concern_max(i)%contribute(j)) THEN
          
          !! target gridpt that src contributes to 
          x_max = recv_concern_max(i)%xtarg(j) 
          y_max = recv_concern_max(i)%ytarg(j) 
          
          xpt = x_max 
          ypt = y_max 
          CALL GLOBAL_TO_LOCAL_GRIDPT(xpt, ypt, targ_grid)
          
          index = ((ypt-1)*l_targ_row_length)+xpt 
          
          IF(GRIDPT_OUTSIDE_PROC_DOMAIN(x_max, y_max, targ_grid)) THEN
            error = -1 
          END IF
          
          !! increment to make space for this contributing src
          contribution(index)%size = contribution(index)%size + 1
        END IF
        
      END DO
    END DO
    
    ! allocate space to actually store src
    ! contribution info
    DO i=1, contrib_size  
      dat_size = contribution(i)%size 
      IF(dat_size>0) THEN
        ALLOCATE(contribution(i)%x(dat_size))
        ALLOCATE(contribution(i)%y(dat_size))
        ALLOCATE(contribution(i)%weight(dat_size))
        ALLOCATE(contribution(i)%contrib_proc(dat_size)) !! debug info
      END IF
    END DO
    
    incr = 0
    
    ! now do everthing above again except allocation, and also now 
    ! initialising the newly allocated space with target contribution 
    ! info
    DO j=1, rows
      DO i=1, row_length  
        
        weight_max = 0
        
        IF(count_tr(i,j) /= 0) THEN 
          
          ! get indices for src which target contributes 
          ! to the most
          DO k=1, count_tr(i,j)
            ip = base_tr(i,j) + k 
            IF(weight(ip) > weight_max) THEN
              weight_max = weight(ip) 
              x_max = MOD(g_index_src(ip)-1,g_targ_row_length)+1
              y_max = (g_index_src(ip)-1)/g_targ_row_length+1
            END IF
          END DO
          
          ! do not consider target points outside this process's domain
          IF(.NOT.GRIDPT_OUTSIDE_PROC_DOMAIN(x_max, y_max, targ_grid))  &
               THEN 
            
            xpt = x_max 
            ypt = y_max 
            CALL GLOBAL_TO_LOCAL_GRIDPT(xpt, ypt, targ_grid)
            
            ! get src 1-d index 
            index = ((ypt-1)*l_targ_row_length)+xpt 
            
            incr(index) = incr(index) + 1
            xpt = i 
            ypt = j 
            CALL LOCAL_TO_GLOBAL_GRIDPT(xpt, ypt, src_grid)
            contribution(index)%x(incr(index)) = xpt
            contribution(index)%y(incr(index)) = ypt 
            contribution(index)%weight(incr(index)) = weight_max 
            contribution(index)%contrib_proc(incr(index)) = mype
            
          END IF
          
        END IF
      END DO
    END DO
    
    ! and now contributing srcs outside the process's 
    ! subdomain
    DO i=1, recv_size  
      
      dat_size = recv_concern_max(i)%size 
      DO j=1, dat_size  
        
        ! add only if target contributes
        IF(recv_concern_max(i)%contribute(j)) THEN
          
          x_max = recv_concern_max(i)%xtarg(j) 
          y_max = recv_concern_max(i)%ytarg(j) 
          
          xpt = x_max 
          ypt = y_max 
          CALL GLOBAL_TO_LOCAL_GRIDPT(xpt, ypt, targ_grid)
          
          index = ((ypt-1)*l_targ_row_length)+xpt  
          incr(index) = incr(index) + 1
          
          contribution(index)%x(incr(index)) = recv_concern_max(i)%x(j)
          contribution(index)%y(incr(index)) = recv_concern_max(i)%y(j)
          contribution(index)%weight(incr(index)) = recv_concern_max(i  &
               )%weight(j)
          contribution(index)%contrib_proc(incr(index)) =               &
               recv_concern_max(i)%proc_num
        END IF
      END DO
    END DO
    
    ! release incrementer resource
    DEALLOCATE(incr)
    
    RETURN
  END SUBROUTINE CALC_CONTRIBUTORS
  
END MODULE REGRID_ALLOC_CALC
