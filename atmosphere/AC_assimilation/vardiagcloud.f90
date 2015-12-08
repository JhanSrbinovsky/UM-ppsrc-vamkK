! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Diagnostic liquid and frozen cloud scheme for deriving IAU increments
!+ Calculate updated cloud variables as documented in VSDP 61

MODULE var_diagcloud_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE Var_DiagCloud(field_size,     & ! in
                         p_theta_levels, & ! in
                         RHc,            & ! in
                         IncrementIce,   & ! in
                         qT,             & ! inout
                         CMessage,       & ! inout
                         ICode,          & ! inout
                         CL,             & ! out   (optional)
                         qCL,            & ! inout (optional)
                         dCLdp,          & ! out   (optional)
                         dCLdqT,         & ! out   (optional)
                         dCLdT,          & ! out   (optional)
                         dqCLdp,         & ! out   (optional)
                         dqCLdqT,        & ! out   (optional)
                         dqCLdT,         & ! out   (optional)
                         CF,             & ! out   (optional)
                         qCF,            & ! inout (optional)
                         qCFmax,         & ! inout (optional)
                         dCFdp,          & ! out   (optional)
                         dCFdqcf,        & ! out   (optional)
                         dCFdT,          & ! out   (optional)
                         dqCFdqCL,       & ! out   (optional)
                         dqCFdT,         & ! out   (optional)
                         T,              & ! in    (optional)
                         BGqcl,          & ! in    (optional)
                         BGqcf,          & ! in    (optional)
                         BGT,            & ! in    (optional)
                         TL)               ! inout (optional)


! Description:
!
!   Documentation located in VAR Scientific Documentation Paper 61
!
!   NOTE 1: Only request gradients when calling from VAR or OPS and
!           the temperature input is T; not TL.
!
!   NOTE 2: This code is intended to mirror the subroutine of the same
!           name held in VarMod_InterpColumns. Differences between UM and
!           GEN versions are commented in-line. These are primarily
!           differences in CALLs to subroutines for calculating
!           saturated vapour pressure and a commenting-out of linear
!           code in the UM-only version. At some future time, it may be
!           possible for the UM, VAR & OPS all to access this code from
!           a single source file.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: FORTRAN 90
!
! Declarations:

USE IAU_mod, ONLY :                  &
    DiagCloud_ApplyCompLimits,       &
    DiagCloud_QN_CompRegimeLimit,    &
    DiagCloud_NumMaxLoops,           &
    DiagCloud_Tol_q,                 & ! tolerarance for qcl==0 and qcf==0 tests
    DiagCloud_Tol_FM                   ! tolerarance for FM==0 tests

USE water_constants_mod, ONLY: lc
USE atmos_constants_mod, ONLY : cp

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

! Common blocks:


! Subroutine arguments:

INTEGER,       INTENT(IN)              ::                         &
    field_size

REAL,          INTENT(IN)              ::                         &
    p_theta_levels (field_size),                                  &
    RHc            (field_size)

LOGICAL,       INTENT(IN)              ::                         &
    IncrementIce                          ! =T for ice incrementing

REAL,          INTENT(IN),    OPTIONAL ::                         &
    T         (:),                                                &
    BGqcl     (:),                                                &
    BGqcf     (:),                                                &
    BGT       (:)

INTEGER,       INTENT(INOUT)           ::                         &
    ICode

REAL,          INTENT(INOUT)           ::                         &
    qT         (field_size)

REAL,          INTENT(INOUT), OPTIONAL ::                         &
    qcl            (:),                                           &
    qcf            (:),                                           &
    qcfMax         (:),                                           &
    TL             (:)

CHARACTER(len=*), INTENT(INOUT)        ::                         &
    CMessage

REAL, INTENT(OUT),   OPTIONAL          ::                         &
    CL             (:),                                           &
    CF             (:)

! These are not marked as OUT purely for the UM
! to stop some compilers from complaining.
! For use in VAR INTENT(OUT) can be used.
REAL, OPTIONAL                         ::                         &
    dCLdp          (:),                                           &
    dCLdqT         (:),                                           &
    dCLdT          (:),                                           &
    dqCLdp         (:),                                           &
    dqCLdqT        (:),                                           &
    dqCLdT         (:),                                           &
    dCFdp          (:),                                           &
    dCFdqcf        (:),                                           &
    dCFdT          (:),                                           &
    dqCFdqCL       (:),                                           &
    dqCFdT         (:)

! Local constants:

CHARACTER(len=*), PARAMETER :: RoutineName='Var_DiagCloud'

! Local variables:

INTEGER           ::              &
    calc,                         &
    loop,                         &
    maxLoops,                     &
    statL,                        &
    statF,                        &
    x

LOGICAL           ::              &
    oscillate(field_size),        &
    converged(field_size),        &
    OutsideCompRegime(field_size)

REAL              ::              &
    A,                            &
    accuracy0,                    &
    accuracy    (field_size),     &
    aL          (field_size),     &
    dqcf        (field_size),     &
    delta       (field_size),     &
    p_temp      (field_size),     &
    qCL_calc    (field_size),     &
    qCL_1       (field_size),     &
    qCL_2       (field_size),     &
    qCL_3       (field_size),     &
    qSatW       (field_size),     &
    qSatW_minus (field_size),     &
    qSatW_plus  (field_size),     &
    qSat_temp   (field_size),     &
    Qc          (field_size),     &
    QN          (field_size),     &
    Tmod        (field_size),     &
    Tmod_temp   (field_size)

REAL, ALLOCATABLE ::              &
    alphaL        (:),            &
    dPsatdT_minus (:),            &
    dPsatdT_plus  (:),            &
    f_of_T        (:),            & ! parameterised F(T)
    f_of_T_plus   (:),            & ! parameterised F(T+T')
    FM            (:),            & ! Fg * M
    Fplus         (:),            & ! incremented Fg
    FplusM        (:),            & ! F^+ * M
    qCL_4         (:),            &
    qcfMin        (:),            &
    qcfMax0       (:),            &
    qcfin         (:),            &
    qcf_1         (:),            &
    qcf_2         (:),            &
    qcf_3         (:),            &
    qcl2_1        (:),            &
    qcl2_2        (:),            &
    qcl2_3        (:),            &
    B             (:),            &
    C             (:),            &
    D             (:),            &
    E             (:),            &
    F             (:),            &
    Fg            (:),            & ! initial guess for F
    G             (:),            &
    GCx           (:),            &
    H             (:),            &
    I             (:),            &
    J             (:),            &
    K             (:),            &
    L             (:),            &
    M             (:),            & ! mult factor for M's
    N             (:),            &
    O             (:),            &
    V             (:),            &
    W             (:),            &
    Y             (:),            &
    Z             (:),            &
    multiplier    (:)

REAL :: test(field_size)

LOGICAL :: assignPseudoBG
LOGICAL :: assign_qCFmax
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

REAL :: huge_a
REAL :: huge_log

!- End of header ------------------------------------------------------

huge_a   = HUGE(a)
huge_log = LOG(huge_a)*0.99

IF (lhook) CALL dr_hook('Var_DiagCloud',zhook_in,zhook_handle)

ICode = 0
statL = 0
statF = 0

IF (DiagCloud_NumMaxLoops < 150) THEN

  WRITE (6,*)'(DiagCloud_NumMaxLoops cannot be set to less than 150.)'
  WRITE (6,*)'(Resetting DiagCloud_NumMaxLoops to 150.)'
  DiagCloud_NumMaxLoops = 150

END IF

! Check that one temperature vector is passed in
IF (.NOT. PRESENT(T) .AND. .NOT. PRESENT(TL)) THEN

  ICode    = 1
  CMessage = 'Error from Var_DiagCloud: No temperature data'

ELSE IF (PRESENT(T) .AND. PRESENT(TL)) THEN

  ICode    = 1
  CMessage = 'Error from Var_DiagCloud: T & TL supplied'

ELSE

!----------------------------------------------------------------------
! [1]: Set fixed and initial values.
!----------------------------------------------------------------------

  ! Fractional accuracy to which cloud water calculation
  accuracy  = 0.001
  accuracy0 = accuracy(1)

  ! Initial qcf increment
  dqcf      = 0.0

  ! Determine whether initial/final values are needed
  !  assignPseudoBG=T when x^s (simplified field) = x^g (bg field)
  IF (PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.                   &
      PRESENT(qcl)   .AND. PRESENT(qcf)) THEN
    assignPseudoBG = (ALL(BGqcl == qcl)                   .AND.   &
                      ALL(BGqcf == qcf)                   .AND.   &
                      IncrementIce) ! bg values have already been assigned
  ELSE
    assignPseudoBG = .FALSE. !x^s /= x^g
  END IF

  ! Allocate extra array if qcfMax is to be calculated
  IF (PRESENT(qcfMax)) THEN
    assign_qCFmax = (ALL(qT == qcfMax)) ! check to see if qcmfmax's redundant
    IF (assign_qCFmax) ALLOCATE ( qCL_4 (field_size) ) ! allocate if redundant
  ELSE
    assign_qCFmax = .FALSE.  ! no  qcfMax to be assigned
  END IF

  ! Setup maximum number of iterations for cloud water calculation
  ! (i) If calculating incremented qcf, allocate arrays, setup initial
  !   values and calculate FM & FplusM
  ! (ii) otherwise just set maxloops
  ! eq numbers refer to VSDP61
  IF (.NOT. AssignPseudoBG .AND. IncrementIce   .AND.             &
            PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.             &
            PRESENT(qcl)   .AND. PRESENT(qcf)   .AND.             &
            .NOT. assign_qCFmax) THEN

    ! [1.1] allocation
    ALLOCATE ( Fg          (field_size) )! ice partioning factor
    ALLOCATE ( f_of_T_plus (field_size) ) ! parameterised F(T+T')
    ALLOCATE ( FM          (field_size) ) ! Fg*M
    ALLOCATE ( FplusM      (field_size) ) ! F+*M
    ALLOCATE ( M           (field_size) ) ! mult factor
    ALLOCATE ( qcfin       (field_size) )
    ALLOCATE ( qcfMin      (field_size) )
    ALLOCATE ( qcfMax0     (field_size) )
    ALLOCATE ( qcf_1       (field_size) )
    ALLOCATE ( qcf_2       (field_size) )
    ALLOCATE ( qcf_3       (field_size) )
    ALLOCATE ( qcl2_1      (field_size) )
    ALLOCATE ( qcl2_2      (field_size) )
    ALLOCATE ( qcl2_3      (field_size) )

    ! Allocate late, deallocate early to reduce memory footprint
    ALLOCATE ( f_of_T      (field_size) ) ! parameterised F(T)
    ALLOCATE ( Fplus       (field_size) ) ! incremented Fg

    ! [1.2] define increments to Fg (ice partitioning factor)
    ! case  qcf^g  qcl^g      F^g
    ! i     >0     >0         eq.1
    ! ii    >0     =0         1
    ! iii   0      >0         0
    ! iv    0      0          0
    f_of_T      = 0.5 * (1.0 - TANH((BGT - 260.65) * 0.18)) ! f(T_bg)
    f_of_T_plus = 0.5 * (1.0 - TANH((T   - 260.65) * 0.18)) ! f(T_bg +T')

    loop1: WHERE (BGqcf > DiagCloud_Tol_q)

      loop1_1: WHERE (BGqcl > DiagCloud_Tol_q) ! case 1: BGqcf /= 0, BGqcl /= 0
        Fg    = BGqcf / (BGqcf + BGqcl)        !   eq. (18), using F^g:=F
        M     = (BGqcf / (BGqcf + qcl)) / Fg   !   eq. (19) using qcf=BGqcf
        Fplus = Fg + (f_of_T_plus - f_of_T) &  ! eq. 22
                   * MIN((MIN(1.0, 1.0 / M) - Fg), Fg)
      ELSEWHERE                                ! case ii: BGqcf /= 0, BGqcl = 0
        Fg     = 1.0
        Fplus = 1.0
        M= (BGqcf / (BGqcf + qcl))             ! no tests needed
      ENDWHERE loop1_1

    ELSEWHERE ! case iii,iv: BGqcf = 0

      Fg     = 0.0 !
      loop1_2: WHERE (f_of_T_plus - f_of_T > 0.0)  ! don't need threshold here.
        Fplus = f_of_T_plus - f_of_T  ! use +ve increment,  T'<0
        M     = 1.0                   ! use M=1 to enable incrementing
      ELSEWHERE
        Fplus = 0.0                   ! don't increment
        M     = 0.0                   ! keep M=0
      ENDWHERE loop1_2

    ENDWHERE loop1


    ! [1.3] define values for iteration
    maxLoops = DiagCloud_NumMaxLoops
    FM     = MIN(MAX(Fg * M, 0.0), 1.0) ! Fg*M needed to solve for qcf (eq. 23)
    FplusM = MIN(MAX(Fplus * M, 0.0), 1.0) ! updated version of F*M

    qcfMin    = 0.0
    qcfin     = qcf
    qcfMax0   = qcfMax
    loop2: WHERE (qcfMax > DiagCloud_Tol_q)
      ! Ensure that qcf is in the allowed range
      qcf       = MAX(MIN(qcfMax,qcfin),qcfMin)
      dqcf      = qcf - qcfin  ! increment wrt initial value
    ELSEWHERE
      ! qcf incrementing not possible
      dqcf = 0.0 ! can't increment when qcfMax=0
    ENDWHERE loop2

    test      = 1.0
    qcf_1     = qcf
    qcf_2     = qcf
    qcf_3     = qcf
    qcl2_1    = 0.0
    qcl2_2    = 0.0
    qcl2_3    = 0.0

    ! Tests for qcf calculation
    ! i. Skip for (1-FM) small. Want to avoide divide-by-zero in
    !   calculation of qcf.
    ! ii. Skip if qcl is already small.
    ! iii. Skip if qcl=BG qcf (to some tolerance)
    loop3: WHERE ( (1.0 - FM) * (1.0 - FM) <= DiagCloud_Tol_FM .OR. &
                   qcl <= DiagCloud_Tol_q                      .OR. &
                   (1.0 - FplusM) * (1.0 - FplusM) <= DiagCloud_Tol_FM .OR.  &
                   BGqcf == BGqcf + accuracy * qcl * (1.0-accuracy) )

                   ! these points don't need to be iterated
                   test = huge_a   ! set test flag to divide by (1-FM)
                   qcf=qcfin
                   dqcf = 0.0
    ENDWHERE loop3


  ELSE

    maxLoops = 50 ! default value

  END IF

  ! [1.4] deallocate
  IF  (ALLOCATED(FPlus))  DEALLOCATE (FPlus)
  IF  (ALLOCATED(f_of_T)) DEALLOCATE (f_of_T)


  ! [1.5] calculate qsat, sat specific humidity wrt water
  IF (PRESENT(TL)) THEN

! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW, TL, p_theta_levels, field_size)

    Tmod = TL

  ELSE IF (PRESENT(T)) THEN

! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW_plus,     T+0.05,                        &
                    p_theta_levels, field_size)
! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW_minus,    T-0.05,                        &
                    p_theta_levels, field_size)

    aL = 1.0 /                                                    &
        (1.0 + Lc * ((qSatW_plus - qSatW_minus) / 0.1) / Cp)

    Tmod = T

  END IF

!----------------------------------------------------------------------
! [2]: Calculate qcl & derived quantities
!----------------------------------------------------------------------

  OutsideCompRegime = .FALSE.

  converged = .FALSE.
  loop      = 0
  oscillate = .FALSE.
  qCL_calc  = 0.0
  qCL_1     = 0.0
  qCL_2     = 0.0
  qCL_3     = 0.0
  IF (assign_qCFmax) qCL_4 = 0.0

  DO WHILE (COUNT(converged) + COUNT(oscillate) < field_size .AND.    &
             loop < maxLoops)

    loop = loop + 1
    calc = field_size - COUNT(converged)

    Tmod_temp(1:calc) = PACK(Tmod, .NOT. converged)

    IF (PRESENT(T)) THEN

      p_temp   (1:calc) = PACK(p_theta_levels, .NOT. converged)
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (qSat_temp(1:calc),                           &
                      Tmod_temp(1:calc),                          &
                      p_temp   (1:calc),                          &
                      calc)
      qSatW = UNPACK(qSat_temp(1:calc), .NOT. converged, qSatW)


    ELSE IF (PRESENT(TL)) THEN

      p_temp   (1:calc) = PACK(p_theta_levels, .NOT. converged)
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (qSat_temp(1:calc),                           &
                      Tmod_temp(1:calc)+0.05,                     &
                      p_temp   (1:calc),                          &
                      calc)
      qSatW_plus = UNPACK(qSat_temp(1:calc),                      &
            .NOT. converged, qSatW_plus)
! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (qSat_temp(1:calc),                           &
                      Tmod_temp(1:calc)-0.05,                     &
                      p_temp   (1:calc),                          &
                      calc)
      qSatW_minus = UNPACK(qSat_temp(1:calc),                     &
          .NOT. converged, qSatW_minus)

      WHERE (.NOT. converged)

        aL = 1.0 /                                                &
            (1.0 + Lc * ((qSatW_plus - qSatW_minus) * 10.0) / Cp)

      ENDWHERE

    END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x)

!$OMP DO SCHEDULE(STATIC)
    DO x=1, field_size 
    
      oscillate(x) = .FALSE.

      ! When assigning qCFmax, set qT such that we find the point
      !   where qT = qcl
      IF (assign_qCFmax) THEN
        IF (loop == 0) THEN
          qT(x)  = 0.0
        ELSE 
          IF (qCL_calc(x) > DiagCloud_Tol_q) THEN
            qT(x) = qCL_calc(x)
          ELSE 
            qT(x) = 0.0
          END IF
        END IF
      END IF

      IF (.NOT. converged(x)) THEN

        Qc(x) = aL(x) * ((qT(x)-dqcf(x)) - qSatW(x))
        delta(x)    = qSatW(x) * aL(x) * (1.0 - RHc(x)) * 0.25 !RHc not changing
        QN(x)       = Qc(x) / delta(x)
        QN(x)       = MAX(QN(x), -huge_log)
        qCL_calc(x) = Qc(x) + delta(x) * LOG( EXP(-QN(x)) + 1.0 )

      END IF

      IF (DiagCloud_ApplyCompLimits) THEN

        ! Avoid tiny deviations of qCL_calc from its bounding values that may
        ! lead to machine precision issues later on:
        IF (.NOT. converged(x) .AND.                                    &
            QN(x) < -DiagCloud_QN_CompRegimeLimit) THEN
          qCL_calc(x) = 0.0
          OutsideCompRegime(x) = .TRUE.
        ELSE IF (.NOT. converged(x) .AND.                               &
                 QN(x) >  DiagCloud_QN_CompRegimeLimit) THEN
          qCL_calc(x) = Qc(x)
          OutsideCompRegime(x) = .TRUE.
        ELSE 
          OutsideCompRegime(x) = .FALSE.
        END IF

      END IF

      IF (IncrementIce) THEN
        
        IF (.NOT. converged(x)) THEN
          
          qCL_calc(x) = MAX(qCL_calc(x), 0.0)
          
        END IF
        
      END IF

    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO x=1, field_size 
      
      IF (PRESENT(TL)) THEN
        
        IF (.NOT. converged(x)) THEN
          Tmod(x) = TL(x) + Lc * qCL_calc(x) / Cp
        END IF
        
      ELSE IF (PRESENT(T)) THEN
        
        IF(.NOT. converged(x)) THEN
          Tmod(x) = T(x)  - Lc * qCL_calc(x) / Cp
        END IF
        
      END IF
    END DO
!$OMP END DO NOWAIT
    
    ! Apply convergence criteria
    IF (IncrementIce .AND. .NOT. assign_qCFmax) THEN

!$OMP DO SCHEDULE(STATIC)
      DO x=1, field_size      
        
        IF (.NOT. converged(x)   .AND.                                  &
            (ABS(qCL_1(x) - qCL_calc(x)) <                              &
             accuracy(x) * MIN(qCL_calc(x), qCL_1(x)) .OR.              &  
             qCL_calc(x) <= DiagCloud_Tol_q )) THEN
          
          converged(x) = .TRUE.
          
        ELSE IF (.NOT. converged(x)   .AND.                             &
                 qCL_2(x) == qCL_calc(x) .AND.                          &
                 qCL_3(x) == qCL_1(x)) THEN
          
          oscillate(x) = .TRUE.
          
        END IF
        
        IF (oscillate(x)) THEN
          
          qCL_calc(x)  = (qCL_calc(x) + qCL_1(x)) * 0.5
          oscillate(x) = .FALSE.
          converged(x) = .TRUE.
          
        END IF

      END DO
!$OMP END DO NOWAIT
    ELSE IF (assign_qCFmax) THEN

!$OMP DO SCHEDULE(STATIC)
      DO x=1, field_size      
        
        IF (.NOT. converged(x) .AND.                                    &
             ((qCL_calc(x) == qCL_1(x) .AND.                            &
               qCL_1(x)    == qCL_2(x) .AND.                            &
               qCL_2(x)    == qCL_3(x)         ) .OR.                   &
              (qCL_2(x)    == qCL_calc(x) .AND.                         &
               qCL_3(x)    == qCL_1(x)         ) .OR.                   &
              (qCL_3(x)    == qCL_calc(x)      ) .OR.                   &
              (qCL_4(x)    == qCL_calc(x)      ))) THEN
          
          converged(x) = .TRUE.
          IF (qcfMax(x) < qCL_calc(x)) THEN
            
            qcfMax(x) = 0.0
            
          ELSE 
            
            qcfMax(x) = MAX(qcfMax(x) - qT(x), BGqcf(x)) * (1.0 + accuracy(x))
            
          END IF
          
        END IF
        
        IF (loop == maxLoops) THEN
          
          IF (.NOT. converged(x)) THEN
            
            converged(x) = .TRUE.
            IF (qcfMax(x) < qCL_calc(x)) THEN
              qcfMax(x) = 0.0
            ELSE 
              qcfMax(x) = MAX(qcfMax(x) - qT(x), BGqcf(x)) *            &
              (1.0 + accuracy(x))
            END IF
            
          END IF
          
        END IF
        
      END DO
!$OMP END DO NOWAIT

    ELSE

!$OMP DO SCHEDULE(STATIC)
      DO x=1, field_size      

        IF (.NOT. converged(x) .AND.                                    &
             ABS(qCL_1(x) - qCL_calc(x)) <                              &
             accuracy(x) * MIN(qCL_calc(x), qCL_1(x))) THEN
          
          converged(x) = .TRUE.
          
        ELSE IF (.NOT. converged(x)      .AND.                          &
                 qCL_2(x) == qCL_calc(x) .AND.                          &
                 qCL_3(x) == qCL_1(x)) THEN
          
          oscillate(x) = .TRUE.
          
        END IF
      END DO
!$OMP END DO NOWAIT

    END IF

!$OMP DO SCHEDULE(STATIC)
    DO x=1, field_size 
      
      IF (assign_qCFmax) THEN
        
        IF(.NOT. converged(x)) THEN
          
          qCL_4(x) = qCL_3(x)
          
        END IF
        
      END IF
      
      IF(.NOT. converged(x)) THEN
        
        qCL_3(x)    = qCL_2(x)
        qCL_2(x)    = qCL_1(x)
        qCL_1(x)    = qCL_calc(x)
        
      END IF

    END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL       

    ! Control iteration for calculating incremented cloud water
    IF (.NOT. AssignPseudoBG .AND. IncrementIce   .AND.           &
              PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.           &
              PRESENT(qcl)   .AND. PRESENT(qcf)   .AND.           &
        .NOT. (assign_qCFmax .OR. loop == maxloops)) THEN

      ! test == huge_a serves as an indicator that the qcf
      !   calculation is converged
      IF (ANY(converged)) THEN
        
        DO x=1,field_size

          IF (converged(x) .AND. loop /= maxloops) THEN
            IF(test(x) /= huge_a) THEN
              IF (qcfMax(x) <= DiagCloud_Tol_q ) THEN

                qcf(x)  = qcfin(x)
                test(x) = huge_a
                converged(x) = .FALSE.

              ELSE IF (OutsideCompRegime(x) .AND. &
                       qCL_calc(x) <= DiagCloud_Tol_q ) THEN

                qcf(x)  = qcfin(x)
                test(x) = huge_a
                converged(x) = .FALSE.

              ELSE

                ! Set test = the next value of qcf
                test(x) = qCL_calc(x) * (FplusM(x) / (1.0 - FplusM(x)))

              END IF

            ELSE IF(test(x) == huge_a .AND.                           &
                    qCL_calc(x) <= DiagCloud_Tol_q .AND.              &
                    qcf(x)      <= DiagCloud_Tol_q) THEN

              ! Ensure that qcf increments are zero if qCL_calc is zero
              qcf(x) = qcfin(x)
              IF (qcfin(x) > DiagCloud_Tol_q) converged(x) = .FALSE.

            END IF

            IF (test(x) == huge_a) THEN

             ! Do nothing

            ELSE IF ((test(x) == 0.0 .AND. qCF(x) <= DiagCloud_Tol_q) .OR.  &
                (test(x) == test(x) + accuracy0 *          &
                 qCL_calc(x) * (1.0 - accuracy0) .AND.     &
                 qcf(x)  == qcf(x) + accuracy0 *           &
                 qCL_calc(x) * (1.0 - accuracy0) .AND.     &
                 test(x) /= 0.0 .AND. qcf(x) > DiagCloud_Tol_q)) THEN

              ! Ensure that qcf increments are zero where qcl is tiny
              IF (qcfin(x) > DiagCloud_Tol_q ) converged(x) = .FALSE.

              qcf(x)  = qcfin(x)
              test(x) = huge_a

            ELSE IF ((qCL_calc(x) == qCL_calc(x) +                  &
                      accuracy0 * test(x) * (1.0 - accuracy0) .AND. &
                      qcl2_1(x)  == qcl2_1(x)  +                    &
                      accuracy0 * test(x) * (1.0 - accuracy0) .AND. &
                      qCL_calc(x) > DiagCloud_Tol_q           .AND. &
                      qcl2_1(x)   > DiagCloud_Tol_q)) THEN

              ! qcf calculation is converged if qcf in cannot be seen
              !   by qcl
              test(x) = huge_a

            ELSE

              qcf(x) = qCL_calc(x) * (FplusM(x) / (1.0 - FplusM(x)))

              IF((qcf(x) == qcf_1(x) .AND. qcf(x) > DiagCloud_Tol_q .AND.    &
                  qcf(x) /= qcfMax0(x))               .OR.        &
                  ((ABS(qCL_calc(x) - qcl2_1(x)) < accuracy0 *    &
                  MIN(qcl2_1(x), qCL_calc(x)) .OR.                &
                  qCL_calc(x) == qcl2_1(x))          .AND.        &
                  (ABS(qcf(x) - qcf_1(x)) < accuracy0 *           &
                  MIN(qcf_1(x), qcf(x)) .OR. qcf(x) == qcf_1(x)))) THEN

                ! qcf calculation is converged if qcf and qcl changes
                !   are within predefined fractional accuracy
                test(x)=huge_a

              ELSE IF (qCL_calc(x) > DiagCloud_Tol_q  .AND.     &
                       qcl2_1(x) > DiagCloud_Tol_q    .AND.     &
                      (qCL_calc(x) == qcl2_2(x) .OR. &
                       qCL_calc(x) == qcl2_1(x))) THEN

                ! qcl calculation is oscillating
                IF(qCL_calc(x) == qCL_1(x)) THEN

                  ! qcf is converged if oscillation is null
                  test(x)=huge_a
                  qCL_calc(x) = 0.5*(qCL_calc(x) + qcl2_1(x))
                  qcf(x) = qCL_calc(x) * (FplusM(x) / (1.0 - FplusM(x)))

                ELSE

                  ! tighten accuracy requirement
                  accuracy(x)  = accuracy(x) * 0.1
                  converged(x) = .FALSE.
                  qcf_3(x)  = qcf_2(x)
                  qcf_2(x)  = qcf_1(x)
                  qcf_1(x)  = qcf(x)
                  qcl2_3(x) = qcl2_2(x)
                  qcl2_2(x) = qcl2_1(x)
                  qcl2_1(x) = qCL_calc(x)

                END IF

              ELSE IF((ABS(qcfMax(x) - qcfMin(x)) <                         &
                       accuracy(x) * MIN(qcfMax(x), qcfMin(x))   .OR.       &
                       qcfMax(x) == qcfMin(x)) .OR.                         &
                      (qcfMax(x) < qcfMax0(x) .AND.                         &
                       qcfMin(x) > DiagCloud_Tol_q .AND.                    &
                       ABS(qcfMax(x) / (1.0 + Accuracy(x)) -                &
                       qcfMin(x) * (1.0 + Accuracy(x))) <                   &
                       accuracy(x) * (qcfMax(x) / (1.0 + Accuracy(x)) +     &
                       qcfMin(x) * (1.0 + Accuracy(x))))) THEN

                ! Bounds on qcf calculation are within accuracy
                IF (qcf(x) > qcfMax(x) .OR. qcf(x) < qcfMin(x)) &
                  qcf(x) = 0.5*(qcfMax(x)+qcfMin(x))

                qcf_3(x)     = qcf(x)
                qcf_2(x)     = qcf(x)
                qcf_1(x)     = qcf(x)
                qcl2_3(x)    = qcl2_2(x)
                qcl2_2(x)    = qcl2_1(x)
                qcl2_1(x)    = qCL_calc(x)
                converged(x) = .FALSE.
                accuracy(x)  = 0.1 * accuracy(x)

              ELSE IF(((qcf(x) < qcf_1(x) .AND. qcf_1(x) > qcf_2(x)) .OR.   &
                       (qcf(x) > qcf_1(x) .AND. qcf_1(x) < qcf_2(x)) ).AND. &
                      ABS(qcf(x) - qcf_1(x)) > ABS(qcf_1(x) - qcf_2(x))) THEN

                ! qcf solution diverging
                qcf(x) = MAX(MIN(qcf(x),qcfMax(x)),qcfMin(x))
                qcf(x) = qcf_2(x) + (qcf_1(x) - qcf_2(x)) * (qcf_1(x) -  &
                         qcf_2(x)) / ((qcf_1(x) - qcf_2(x)) + &
                         (qcf_1(x) - qcf(x)))

                IF (qcf(x) < qcf_1(x)) THEN
                  qcfMax(x)=MIN(qcfMax(x),qcf_1(x)*(1.0+Accuracy(x)))
                  qcfMin(x)=MAX(qcfMin(x),qcf_2(x)/(1.0+Accuracy(x)))
                END IF

                IF (qcf(x) > qcf_1(x)) THEN
                  qcfMax(x)=MIN(qcfMax(x),qcf_2(x)*(1.0+Accuracy(x)))
                  qcfMin(x)=MAX(qcfMin(x),qcf_1(x)/(1.0+Accuracy(x)))
                END IF

                qcf_3(x)     = qcf(x)
                qcf_2(x)     = qcf(x)
                qcf_1(x)     = qcf(x)
                qcl2_3(x)    = qcl2_2(x)
                qcl2_2(x)    = qcl2_1(x)
                qcl2_1(x)    = qCL_calc(x)
                converged(x) = .FALSE.

              ELSE IF (((qcf(x) < qcf_1(x) .AND. qcf_1(x) > qcf_2(x)) .OR. &
                        (qcf(x) > qcf_1(x) .AND. qcf_1(x) < qcf_2(x)))) THEN

                ! qcf solution converging
                IF (qcf(x) < qcf_1(x)) THEN
                  qcfMax(x)=MIN(qcfMax(x),qcf_1(x)*(1.0+Accuracy(x)))
                  qcfMin(x)=MAX(qcfMin(x),qcf_2(x)/(1.0+Accuracy(x)))
                END IF

                IF (qcf(x) > qcf_1(x)) THEN
                  qcfMax(x)=MIN(qcfMax(x),qcf_2(x)*(1.0+Accuracy(x)))
                  qcfMin(x)=MAX(qcfMin(x),qcf_1(x)/(1.0+Accuracy(x)))
                END IF

                qcf(x) = qcf_2(x) + (qcf_1(x) - qcf_2(x)) * &
                         (qcf_1(x) - qcf_2(x)) /  &
                         ((qcf_1(x) - qcf_2(x)) + (qcf_1(x) - qcf(x) ))

                qcf_3(x)     = qcf(x)
                qcf_2(x)     = qcf(x)
                qcf_1(x)     = qcf(x)
                qcl2_3(x)    = qcl2_2(x)
                qcl2_2(x)    = qcl2_1(x)
                qcl2_1(x)    = qCL_calc(x)
                converged(x) = .FALSE.

              ELSE IF(((qcf(x) < qcf_1(x) .AND. qcf_1(x) < qcf_2(x)).OR.  &
                      (qcf(x) > qcf_1(x) .AND. qcf_1(x) > qcf_2(x)))) THEN

                ! qcf solution increasing or decreasing
                IF (qcf(x) < qcf_1(x)) THEN
                  qcfMax(x)=MIN(qcfMax(x),qcf_1(x)*(1.0+Accuracy(x)))
                ELSE IF (qcf(x) > qcf_1(x)) THEN
                  qcfMin(x)=MAX(qcfMin(x),qcf_1(x)/(1.0+Accuracy(x)))
                END IF

                qcf_3(x)     = qcf_2(x)
                qcf_2(x)     = qcf_1(x)
                qcf_1(x)     = qcf(x)
                qcl2_3(x)    = qcl2_2(x)
                qcl2_2(x)    = qcl2_1(x)
                qcl2_1(x)    = qCL_calc(x)
                converged(x) = .FALSE.

              ELSE IF(qcl2_1(x) == qcl2_3(x)   .AND. &
                      qCL_calc(x) == qcl2_2(x) .AND. &
                      (qcfMax(x) == qcfMin(x)  .AND. &
                      qcfMax(x) /= qcfMax0(x)  .AND. &
                      qcfMax(x) > DiagCloud_Tol_q .AND.         &
                      qcf(x) == qcfMax(x))) THEN

                ! qcf solution converged to machine accuracy
                test(x)   = huge_a
                qcf_3(x)  = qcf(x)
                qcf_2(x)  = qcf(x)
                qcf_1(x)  = qcf(x)
                qcl2_3(x) = qcl2_2(x)
                qcl2_2(x) = qcl2_1(x)
                qcl2_1(x) = qCL_calc(x)

              ELSE

                IF(.NOT. (FplusM(x) <= SQRT(DiagCloud_Tol_FM) )  .AND. &
                   .NOT. (qcf(x) /= qcf_1(x) .AND. &
                    qcf_1(x) == qcf_2(x))) THEN

                  ! Preceding logic should make this section unused
                  test(x)      = huge_a
                  qcf(x)       = qcfin(x)
                  converged(x) = .FALSE.

                ELSE

                  converged(x) = .FALSE.

                  IF (qcf_1(x) < qcfMin(x) .OR. qcf_1(x) > qcfMax(x)) THEN

                    ! qcf solution ouside known bounds
                    qcf(x)    = 0.5*(qcfMin(x)+qcfMax(x))
                    qcf_3(x)  = qcf(x)
                    qcf_2(x)  = qcf(x)
                    qcf_1(x)  = qcf(x)
                    qcfMax(x) = qcfMax0(x)
                    qcfMin(x) = 0.0
                    qcf_1(x)  = qcf(x)
                    qcf_2(x)  = qcf(x)
                    qcf_3(x)  = qcf(x)
                    qcl2_1(x) = qcl(x)
                    qcl2_2(x) = qcl(x)
                    qcl2_3(x) = qcl(x)

                  ELSE IF (qcfMax(x) <= qcfMin(x)) THEN

                    ! Computational precision issue or convergence
                    qcfMax(x) = qcfMax0(x)
                    qcfMin(x) = 0.0
                    qcf_3(x)  = qcf(x)
                    qcf_2(x)  = qcf(x)
                    qcf_1(x)  = qcf(x)

                  ELSE IF (qcf(x) < qcfMin(x)) THEN

                    qcfMax(x)=MIN(qcfMax(x),qcf_1(x)*(1.0+Accuracy(x)))

                    IF(qcf_3(x) < qcf_2(x) .AND. qcf_2(x) == qcf_1(x) .AND. &
                      ABS(qcf_1(x) - qcf_3(x)) < accuracy(x) * &
                      (qcf_1(x)+qcf_3(x))) THEN

                      ! Need greater accuracy in order to converge
                      accuracy(x) = accuracy(x)*0.1
                      qcf(x)      = 0.5*(qcf_3(x)+qcf_1(x))
                      qcf_3(x)    = qcf(x)
                      qcf_2(x)    = qcf(x)
                      qcf_1(x)    = qcf(x)

                    ELSE IF (qcf_3(x) < qcf_2(x) .AND. &
                             qcf_2(x) == qcf_1(x)) THEN

                      ! Reset qcf to know mid-point
                      qcf(x)   = 0.5*(qcf_3(x)+qcf_1(x))
                      qcf_3(x) = qcf(x)
                      qcf_2(x) = qcf(x)
                      qcf_1(x) = qcf(x)

                    ELSE

                      ! Reset qcf to know mid-point
                      qcf(x)   = 0.5*(qcfMin(x)+qcf_1(x))
                      qcf_3(x) = qcf_1(x)
                      qcf_2(x) = qcf(x)
                      qcf_1(x) = qcf(x)

                    END IF

                  ELSE IF (qcf(x) > qcfMax(x)) THEN

                    qcfMin(x)=MAX(qcfMin(x),qcf_1(x)/(1.0+Accuracy(x)))

                    IF (qcf_3(x) > qcf_2(x) .AND. qcf_2(x) == qcf_1(x)  &
                       .AND. ABS(qcf_1(x) - qcf_3(x)) < accuracy(x) *   &
                       (qcf_1(x)+qcf_3(x))) THEN

                      ! Need greater accuracy in order to converge
                      accuracy(x) = accuracy(x)*0.1
                      qcf(x)      = 0.5*(qcf_3(x)+qcf_1(x))
                      qcf_3(x)    = qcf(x)
                      qcf_2(x)    = qcf(x)
                      qcf_1(x)    = qcf(x)

                    ELSE IF(qcf_3(x) > qcf_2(x) .AND. &
                            qcf_2(x) == qcf_1(x)) THEN

                      ! Reset qcf to know mid-point
                      qcf(x)=0.5*(qcf_3(x)+qcf_1(x))
                      qcf_3(x) = qcf(x)
                      qcf_2(x) = qcf(x)
                      qcf_1(x) = qcf(x)

                    ELSE

                      ! Reset qcf to know mid-point
                      qcf(x)=0.5*(qcfMax(x)+qcf_1(x))
                      qcf_3(x) = qcf_1(x)
                      qcf_2(x) = qcf(x)
                      qcf_1(x) = qcf(x)

                    END IF

                  ELSE

                    ! Set up succession of previous qcf values
                    qcf_3(x) = qcf_2(x)
                    qcf_2(x) = qcf_1(x)
                    qcf_1(x) = qcf(x)

                  END IF
                  ! Set up succession of previous qCL_calc values
                  qcl2_3(x) = qcl2_2(x)
                  qcl2_2(x) = qcl2_1(x)
                  qcl2_1(x) = qCL_calc(x)

                END IF
              END IF
            END IF

            ! Assign qcf increment to dqcf
            dqcf(x) = qcf(x) - qcfin(x)

          END IF

        END DO

        IF (loop > maxLoops - 50) THEN

          WHERE(.NOT. converged)

            ! Bail out and just do liquid incrementing if necessary
            qcf   = qcfin
            test  = huge_a
            dqcf  = 0.0

          ENDWHERE

          IF (statF == 0) THEN

            statF = COUNT(.NOT. converged)

          END IF

        END IF
      END IF
    END IF

 END DO

  IF (loop == maxLoops .OR. statF /= 0) THEN

    statL    = COUNT(.NOT. converged)
    ICode    = -1

    IF (PRESENT(TL))THEN
      IF (COUNT(TL < 183.10 .AND. Tmod_temp > 183.10 .AND. &
          .NOT. converged) == COUNT (.NOT. converged)) THEN
        WRITE(CMessage,'(A,I8,A)') 'Var_DiagCloud: all ', &
            COUNT(.NOT.(converged)),                      &
            ' unconverged points at first call due to T~183.15K'
      ELSE
        WRITE(CMessage, '(A,I3,A,I8,A,I8,A,I8)') 'Var_DiagCloud:', loop,  &
        ' loops unconverged, liquid=', statL, ' ice=', statF, ' / ', field_size
      END IF
    ELSE
      WRITE(CMessage, '(A,I3,A,I8,A,I8,A,I8)') 'Var_DiagCloud:', loop,  &
      ' loops unconverged, liquid=', statL, ' ice=', statF, ' / ', field_size
    END IF

  END IF

  !----------------------------------------------------------------------
  ! [3]: Assign output values
  !----------------------------------------------------------------------

  IF (PRESENT(TL)) THEN
    TL = Tmod
  END IF

  IF (PRESENT(qCL)) THEN
    qCL = qCL_calc
  END IF

  ! Calculate qcf, given calculated qcl and inputs of qcl & qcf
  IF (assignPseudoBG .AND. .NOT. assign_qCFmax) THEN

    ALLOCATE ( Fg (field_size)  )
    ALLOCATE ( M (field_size)  )
    ALLOCATE ( FM (field_size) )

    ! [3.1] Calculate F and M
    ! definition of cases follows section [1.2]
    loop4: WHERE (BGqcf > DiagCloud_Tol_q)

      loop4_1: WHERE (BGqcl > DiagCloud_Tol_q) ! case 1: BGqcf /= 0, BGqcl /= 0
        Fg =  BGqcf / (BGqcf + BGqcl)     ! eq. (18), F^g
        M  = (BGqcf / (BGqcf + qcl)) / Fg ! eq. (19) using qcf=BGqcf

      ELSEWHERE !   case ii: BGqcf /= 0, BGqcl = 0
        Fg = 1.0
        M  = (BGqcf / (BGqcf + qcl))  ! no tests needed
      ENDWHERE loop4_1

    ELSEWHERE ! case iii,iv: BGqcf = 0

      Fg    = 0.0
      M     = 0.0

    ENDWHERE loop4

    ! [3.2] calculate qcf
    FM  = MIN(MAX(Fg * M, 0.0), 1.0)
    loop5: WHERE ( (1.0 - FM) * (1.0 - FM) <= DiagCloud_Tol_FM .OR. &
                   qcl <= DiagCloud_Tol_q )
      qcf = BGqcf
    ELSEWHERE
      qcf = qcl * ( FM / (1.0 - FM))
    ENDWHERE loop5

  ELSE IF (.NOT. AssignPseudoBG .AND. IncrementIce   .AND.        &
                 PRESENT(BGqcl) .AND. PRESENT(BGqcf) .AND.        &
                 PRESENT(qcl)   .AND. PRESENT(qcf)) THEN
    ! NULL
  END IF

  ! Calculate CL according to scheme in use
  IF (PRESENT(CL) .AND. .NOT. IncrementIce) THEN

    QN = Qc / (2.0 * delta)
    QN = MIN(MAX(QN,                                              &
                 -huge_log),                                      &
             huge_log)
    CL = ( 1.0 + (exp(QN) - exp(-QN)) /                           &
                 (exp(QN) + exp(-QN))  ) / 2.0

    IF (PRESENT(qCL)) THEN
      WHERE(qCL <= DiagCloud_Tol_q)
        CL = 0.0
      ENDWHERE
      
      WHERE(CL == 0.0)
        qCL = 0.0
      ENDWHERE
    END IF

  ELSE IF (PRESENT(CL) .AND. IncrementIce) THEN

    CL = 1.0 - EXP(-qCL / delta)

    IF (PRESENT(qCL)) THEN
      WHERE(qCL <= DiagCloud_Tol_q)
        CL = 0.0
      ENDWHERE
      
      WHERE(CL == 0.0)
        qCL = 0.0
      ENDWHERE
    END IF

  END IF

  ! Calculate CF
  IF (PRESENT(CF) .AND. IncrementIce) THEN

! DEPENDS ON: qsat_wat
    CALL QSAT_WAT (qSatW_plus,     T,                             &
                    p_theta_levels, field_size)


    CF = 1.0 - EXP(-4.0*qCF/((1.0-RHc)*qSatW_plus))

    IF (PRESENT(qCF)) THEN
      WHERE(qCF <= DiagCloud_Tol_q)
        CF = 0.0
      ENDWHERE
      
      WHERE(CF == 0.0)
        qCF = 0.0
      ENDWHERE
    END IF

  END IF


!----------------------------------------------------------------------
! [4]: Deallocate arrays
!----------------------------------------------------------------------

  IF  (ALLOCATED(alphaL)       ) DEALLOCATE ( alphaL        )
  IF  (ALLOCATED(dPsatdT_minus)) DEALLOCATE ( dPsatdT_minus )
  IF  (ALLOCATED(dPsatdT_plus) ) DEALLOCATE ( dPsatdT_plus  )
  IF  (ALLOCATED(FM)           ) DEALLOCATE ( FM            )
  IF  (ALLOCATED(Fg)           ) DEALLOCATE ( Fg            )
  IF  (ALLOCATED(Fplus)        ) DEALLOCATE ( Fplus         )
  IF  (ALLOCATED(FplusM)       ) DEALLOCATE ( FplusM        )
  IF  (ALLOCATED(f_of_T)       ) DEALLOCATE ( f_of_T        )
  IF  (ALLOCATED(f_of_T_Plus)  ) DEALLOCATE ( f_of_T_plus   )
  IF  (ALLOCATED(multiplier)   ) DEALLOCATE ( multiplier    )
  IF  (ALLOCATED(qCL_4)        ) DEALLOCATE ( qCL_4         )
  IF  (ALLOCATED(qCFmax0)      ) DEALLOCATE ( qCFmax0       )
  IF  (ALLOCATED(qCFmin)       ) DEALLOCATE ( qCFmin        )
  IF  (ALLOCATED(qCFin)        ) DEALLOCATE ( qCFin         )
  IF  (ALLOCATED(qCF_1)        ) DEALLOCATE ( qCF_1         )
  IF  (ALLOCATED(qCF_2)        ) DEALLOCATE ( qCF_2         )
  IF  (ALLOCATED(qCF_3)        ) DEALLOCATE ( qCF_3         )
  IF  (ALLOCATED(qCL2_1)       ) DEALLOCATE ( qCL2_1        )
  IF  (ALLOCATED(qCL2_2)       ) DEALLOCATE ( qCL2_2        )
  IF  (ALLOCATED(qCL2_3)       ) DEALLOCATE ( qCL2_3        )
  IF  (ALLOCATED(B)            ) DEALLOCATE ( B             )
  IF  (ALLOCATED(C)            ) DEALLOCATE ( C             )
  IF  (ALLOCATED(D)            ) DEALLOCATE ( D             )
  IF  (ALLOCATED(E)            ) DEALLOCATE ( E             )
  IF  (ALLOCATED(F)            ) DEALLOCATE ( F             )
  IF  (ALLOCATED(G)            ) DEALLOCATE ( G             )
  IF  (ALLOCATED(GCx)          ) DEALLOCATE ( GCx           )
  IF  (ALLOCATED(H)            ) DEALLOCATE ( H             )
  IF  (ALLOCATED(I)            ) DEALLOCATE ( I             )
  IF  (ALLOCATED(J)            ) DEALLOCATE ( J             )
  IF  (ALLOCATED(K)            ) DEALLOCATE ( K             )
  IF  (ALLOCATED(L)            ) DEALLOCATE ( L             )
  IF  (ALLOCATED(M)            ) DEALLOCATE ( M             )
  IF  (ALLOCATED(N)            ) DEALLOCATE ( N             )
  IF  (ALLOCATED(O)            ) DEALLOCATE ( O             )
  IF  (ALLOCATED(V)            ) DEALLOCATE ( V             )
  IF  (ALLOCATED(W)            ) DEALLOCATE ( W             )
  IF  (ALLOCATED(Y)            ) DEALLOCATE ( Y             )
  IF  (ALLOCATED(Z)            ) DEALLOCATE ( Z             )

END IF

IF (lhook) CALL dr_hook('Var_DiagCloud',zhook_out,zhook_handle)
! When copying this subroutine into VarMod_InterpColumns,
! replace this IF code with the following:
!IF (UseTrace) CALL Gen_TraceExit (RoutineName)

RETURN
END SUBROUTINE Var_DiagCloud

END MODULE var_diagcloud_mod
