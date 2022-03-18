!123456
      ! !*************************************************************************
      !                                                                        !
      !                      ============================                      !
      !                         PROGRAM FUNC_DECOM_ALS                         !
      !                      ============================                      !
      !                                                                        !
      !                                                                        !      
      ! Purpose                                                                !
      ! =======                                                                !
      !                                                                        !
      ! Written by    : Subhajit Banerjee, UC Davis                            !
      ! Date          : April 2015                                             !
      !                                                                        !
      ! Use:  gfortran -fopenmp -g -fcheck=all -Wall -o func_decom_als
      !                                             func_decom_als.f90
      !*************************************************************************
      !
      PROGRAM FUNC_DECOM_ALS
        USE NRTYPE
        USE OMP_LIB
        USE, INTRINSIC :: ISO_Fortran_env
        IMPLICIT NONE

        ! INTERFACE
        INTERFACE
          FUNCTION SEMI_NORM(U, N)
            USE NRTYPE
            REAL(DP)                        :: SEMI_NORM
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(IN)            :: U(1:N)
          END FUNCTION SEMI_NORM
        END INTERFACE

        INTERFACE
          SUBROUTINE PRINT_MATRIX(M, N, A)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)          :: M, N
            REAL(DP), INTENT(IN)              :: A(M,N)
          END SUBROUTINE PRINT_MATRIX
        END INTERFACE

        INTERFACE       
          FUNCTION FUNC1(X, D)
            USE NRTYPE
            REAL(DP)                    :: FUNC1
            INTEGER(I4B), INTENT(IN)    :: D
            REAL(DP), INTENT(IN)        :: X(1:D)
          END FUNCTION FUNC1
        END INTERFACE

        INTERFACE       
          FUNCTION FUNC2(X, D, C, W)
            USE NRTYPE
            REAL(DP)                            :: FUNC2
            INTEGER(I4B), INTENT(IN)            :: D
            REAL(DP), INTENT(IN)                :: C(1:D), W(1:D)
            REAL(DP), INTENT(IN)                :: X(1:D)
          END FUNCTION FUNC2
        END INTERFACE

        INTERFACE       
          FUNCTION FUNC3(X, D)
            USE NRTYPE
            REAL(DP)                    :: FUNC3
            INTEGER(I4B), INTENT(IN)    :: D
            REAL(DP), INTENT(IN)        :: X(1:D)
          END FUNCTION FUNC3
        END INTERFACE

        INTERFACE       
          FUNCTION FUNC4(X, D, T, C, X0)
            USE NRTYPE
            REAL(DP)                    :: FUNC4
            INTEGER(I4B), INTENT(IN)    :: D
            REAL(DP), INTENT(IN)        :: X(1:D)
            REAL(DP), INTENT(IN)        :: C(1:D), X0(1:D)
            REAL(DP), INTENT(IN)        :: T
          END FUNCTION FUNC4
        END INTERFACE

        INTERFACE       
          FUNCTION FUNC5(X, D)
            USE NRTYPE
            REAL(DP)                    :: FUNC5
            INTEGER(I4B), INTENT(IN)    :: D
            REAL(DP), INTENT(IN)        :: X(1:D)
          END FUNCTION FUNC5
        END INTERFACE

        INTERFACE       
          FUNCTION FUNC6(X, D)
            USE NRTYPE
            REAL(DP)                    :: FUNC6
            INTEGER(I4B), INTENT(IN)    :: D
            REAL(DP), INTENT(IN)        :: X(1:D)
          END FUNCTION FUNC6
        END INTERFACE

        INTERFACE
          FUNCTION FUNC7(X, D, K, Q, DIR_NUM)
            USE NRTYPE
            REAL(DP)                    :: FUNC7
            INTEGER(I4B), INTENT(IN)    :: D, Q, DIR_NUM
            REAL(DP), INTENT(IN)        :: K
            REAL(DP), INTENT(IN)        :: X(1:D)
          END FUNCTION FUNC7
        END INTERFACE

        INTERFACE
          FUNCTION FUNC8(X, D, K, Q, DIR_NUM)
            USE NRTYPE
            REAL(DP)                    :: FUNC8
            INTEGER(I4B), INTENT(IN)    :: D, Q, DIR_NUM
            REAL(DP), INTENT(IN)        :: K
            REAL(DP), INTENT(IN)        :: X(1:D)
          END FUNCTION FUNC8
        END INTERFACE

        INTERFACE
          FUNCTION FUNC9(L, M, X0, ALPHA, X)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)    :: L, M
            REAL(DP), INTENT(IN)        :: X(1:3), X0(1:3)
            REAL(DP), INTENT(IN)        :: ALPHA
            REAL(DP)                    :: FUNC9
          END FUNCTION FUNC9
        END INTERFACE


        INTERFACE
          SUBROUTINE ALS_ORTH(ERR_DECAY, R_OUT, C, S, D, M, R_MAX,      &
     &      C_MIN, C_MAX, POLY_TYPE, U, SAMPLE_SIZE, COORD, EPS,        &
     &      RUNDATAFILE, SOLVER, PRBLM_TYPE, FLAG_CONV,                 &
     &      MAX_ITER, T_REGU_STAT, LAMBDA, LOW_LIM, UP_LIM)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)             :: D, M, R_MAX
            INTEGER(I4B), INTENT(IN)             :: SAMPLE_SIZE
            INTEGER(I4B), INTENT(IN)             :: MAX_ITER
            REAL(DP), ALLOCATABLE, INTENT(OUT)   :: C(:,:,:)
            INTEGER(I4B), INTENT(OUT)            :: R_OUT
            LOGICAL, INTENT(OUT)                 :: FLAG_CONV
            REAL(DP), INTENT(IN)                 :: C_MIN, C_MAX
            REAL(DP), INTENT(OUT)                :: ERR_DECAY
            CHARACTER(80), INTENT(IN)            :: POLY_TYPE
            REAL(DP), INTENT(IN)                 :: U(SAMPLE_SIZE)
            REAL(DP), INTENT(IN)                 :: COORD(SAMPLE_SIZE,D)
            REAL(DP), INTENT(IN)                 :: EPS, LAMBDA
            REAL(DP), ALLOCATABLE, INTENT(OUT)   :: S(:)
            CHARACTER(80), INTENT(IN)            :: RUNDATAFILE, SOLVER
            CHARACTER(80), INTENT(IN)            :: PRBLM_TYPE
            CHARACTER(80), INTENT(IN)            :: T_REGU_STAT
            REAL(DP), OPTIONAL, INTENT(IN)       :: LOW_LIM, UP_LIM
          END SUBROUTINE ALS_ORTH
        END INTERFACE

        INTERFACE
          SUBROUTINE ALS_BRN(ERR_DECAY, R_OUT, C, S, D, M, R_MAX, C_MIN,&
     &                   C_MAX, U, SAMPLE_SIZE, COORD, EPS, LOW_LIM,    &
     &                   UP_LIM, RUNDATAFILE, SOLVER, PRBLM_TYPE,       &
     &                   FLAG_CONV, MAX_ITER, T_REGU_STAT, LAMBDA)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)             :: D, M, R_MAX
            INTEGER(I4B), INTENT(IN)             :: MAX_ITER
            INTEGER(I4B), INTENT(IN)             :: SAMPLE_SIZE
            REAL(DP), ALLOCATABLE, INTENT(OUT)   :: C(:,:,:)
            INTEGER(I4B), INTENT(OUT)            :: R_OUT
            LOGICAL, INTENT(OUT)                 :: FLAG_CONV
            REAL(DP), INTENT(IN)                 :: C_MIN, C_MAX
            REAL(DP), INTENT(OUT)                :: ERR_DECAY
            REAL(DP), INTENT(IN)                 :: U(SAMPLE_SIZE)
            REAL(DP), INTENT(IN)                 :: COORD(SAMPLE_SIZE,D)
            REAL(DP), INTENT(IN)                 :: EPS, LAMBDA
            REAL(DP), INTENT(IN)                 :: LOW_LIM, UP_LIM
            REAL(DP), ALLOCATABLE, INTENT(OUT)   :: S(:)
            CHARACTER(80), INTENT(IN)            :: RUNDATAFILE, SOLVER
            CHARACTER(80), INTENT(IN)            :: PRBLM_TYPE
            CHARACTER(80), INTENT(IN)            :: T_REGU_STAT
          END SUBROUTINE ALS_BRN
        END INTERFACE        

        INTERFACE
          SUBROUTINE NUM_INTEGRATION(INTEGRAL_OUT, C, S, POLY_TYPE, M,  &
     &                           D, R, INT_ORDER, A, B)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)            :: M, D, INT_ORDER, R
            REAL(DP), INTENT(IN)                :: A, B
            REAL(DP), INTENT(IN)                :: C(1:R, 1:(M+1),1:D)
            REAL(DP), INTENT(IN)                :: S(1:R)
            CHARACTER(80)                       :: POLY_TYPE
            REAL(DP), INTENT(OUT)               :: INTEGRAL_OUT
          END SUBROUTINE NUM_INTEGRATION
        END INTERFACE

        INTERFACE
          SUBROUTINE NUM_INTEGRATION_BRN(INTEGRAL_OUT, C, S, M, D, R,   &
     &                               INT_ORDER, A, B)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)            :: M, D, INT_ORDER, R
            REAL(DP), INTENT(IN)                :: A, B
            REAL(DP), INTENT(IN)                :: C(1:R, 1:(M+1),1:D)
            REAL(DP), INTENT(IN)                :: S(1:R)
            REAL(DP), INTENT(OUT)               :: INTEGRAL_OUT
          END SUBROUTINE NUM_INTEGRATION_BRN
        END INTERFACE

        INTERFACE
          subroutine timestamp()
          
          END subroutine timestamp
        END INTERFACE

        INTERFACE
          SUBROUTINE POSTPROCESSOR(PLOTFILE, NFINEDIV, FINEGRID, C, S,  &
     &                             R, D, M,POLY_TYPE, LOW_LIM, UP_LIM)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)    :: D, M, NFINEDIV, R
            CHARACTER(80), INTENT(IN)   :: POLY_TYPE
            REAL(DP), INTENT(IN)        :: FINEGRID(1:(NFINEDIV+1)**D,D)
            REAL(DP), INTENT(IN)        :: C(1:R, 1:(M+1), 1:D)
            REAL(DP), INTENT(IN)        :: S(1:R)
            REAL(DP), INTENT(IN)        :: LOW_LIM, UP_LIM
            CHARACTER(80), INTENT(OUT)  :: PLOTFILE
          END SUBROUTINE POSTPROCESSOR
        END INTERFACE

        INTERFACE
          SUBROUTINE COEFF_FILE_HELM(FUNC_TYPE_NUM,  WAVE_NUM, Q,       &
     &                               DIR_NUM, C, R, S, D, M)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: FUNC_TYPE_NUM, Q
            INTEGER(I4B), INTENT(IN)            :: R, M, D, DIR_NUM
            REAL(DP), INTENT(IN)                :: WAVE_NUM
            REAL(DP), INTENT(IN)                :: C(1:R, 1:(M+1), 1:D)
            REAL(DP), INTENT(IN)                :: S(1:R)
          END SUBROUTINE COEFF_FILE_HELM
        END INTERFACE

        INTERFACE
          SUBROUTINE CARTPROD2D(XY, X, Y, M, N)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: M, N
            REAL(DP), INTENT(IN)            :: X(1:M), Y(1:N)
            REAL(DP), INTENT(OUT)           :: XY(M*N,2)
          END SUBROUTINE CARTPROD2D
        END INTERFACE

        INTERFACE
          SUBROUTINE CARTPROD3D(XYZ, X, Y, Z, M, N, P)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: M, N, P
            REAL(DP), INTENT(IN)            :: X(1:M), Y(1:N), Z(1:P)
            REAL(DP), INTENT(OUT)           :: XYZ(M*N*P,3)
          END SUBROUTINE CARTPROD3D
        END INTERFACE

        ! END INTERFACE

        ! Local Varibales:
        CHARACTER(80)               :: GRID_TYPE, POLY_TYPE, PRBLM_TYPE
        CHARACTER(80)               :: XFILE, GXFILE
        INTEGER(I4B)                :: SAMPLE_SIZE, MAX_ITER
        INTEGER(I4B)                :: R_OUT
        REAL(DP), ALLOCATABLE       :: COORD(:,:)
        REAL(DP), ALLOCATABLE       :: U(:)
        INTEGER(I4B)                :: NUM_THREADS = 4
        REAL(DP)                    :: ID
        INTEGER(I4B)                :: I, J, K
        INTEGER(I4B)                :: D, R_MAX, INT_ORDER, M
        INTEGER(I4B), ALLOCATABLE   :: M_ARR(:)
        INTEGER(I4B)                :: LEN_M_ARR 
        REAL(DP)                    :: EPS, TOL, WAVE_NUM, LAMBDA
        REAL(DP)                    :: LOW_LIM, UP_LIM, INTEGRAL_OUT
        INTEGER(I4B)                :: Q, DIR_NUM
        REAL(DP)                    :: ALPHA
        INTEGER(I4B)                :: ELL, EMM
        REAL(DP), ALLOCATABLE       :: X0(:), CFACT(:)
        REAL(DP), ALLOCATABLE       :: C(:,:,:)
        REAL(DP)                    :: ERR_DECAY
        REAL(DP), ALLOCATABLE       :: S(:)
        CHARACTER(80)               :: PROBLEMDATAFILE, OUTFILENAME
        CHARACTER(80)               :: RUNDATAFILE, COORDFILE
        CHARACTER(80)               :: PLOTFILE
        INTEGER(I4B)                :: NFINEDIV
        REAL(DP), ALLOCATABLE       :: FINEGRID(:,:)
        REAL(DP), ALLOCATABLE       :: XFINE(:), YFINE(:), ZFINE(:)
        REAL(DP)                    :: HFINEX, HFINEY, HFINEZ
        CHARACTER(80)               :: T_REGU_STAT
        INTEGER(I4B)                :: FUNC_TYPE_NUM
        REAL(DP), ALLOCATABLE       :: C_PARAM(:), W_PARAM(:)
        REAL(DP)                    :: C_MIN, C_MAX
        INTEGER(I4B)                :: COUNTER = 0
        REAL(DP), ALLOCATABLE       :: INT_MAT(:), ERR_MAT(:)
        INTEGER(I4B), ALLOCATABLE   :: R_MAT(:)
        CHARACTER(80)               :: SOLVER
        INTEGER(I4B)                :: T1, T2, CLOCK_RATE
        INTEGER(I4B)                :: CLOCK_MAX
        LOGICAL                     :: FLAG_CONV
        ! YET TO FIGURE OUT HOW TO USE THIS TIMING SCHEME
        !real ( kind = 4 ) t_array(2)
        !real ( kind = 8 ) time1
        !real ( kind = 8 ) time2
        !integer ( kind = 4 ) rep
        !integer ( kind = 4 ) l_log

        ! Main program begins:
        ! Defining the Problem Domain:
        WRITE(*,*)"Begining Execution..."
        WRITE(*,*)" "
        WRITE(*,'(A)',ADVANCE="NO")"Specify the INPUT data fie:  "
        READ(*,*) PROBLEMDATAFILE
        WRITE(*,*)""
        WRITE(*,*)"Processing problem data..."
        WRITE(*,*)""
        OPEN(UNIT = 11, FILE = PROBLEMDATAFILE, STATUS = 'OLD',         &
     &       ACTION = 'READ')                   ! read problem data
        READ(11, *, END = 110)                  ! header: GRID TYPE
        READ(11, *, END = 110) GRID_TYPE        ! GRID TYPE                        
        SELECT CASE (GRID_TYPE)
          CASE ('STRUCTURED')
            READ(11, *, END = 110)                  ! header: # OF NODES
            READ(11, *, END = 110) SAMPLE_SIZE      ! # OF NODES
            ! Number of scattered data-sites (have to find a reasonable data-site
            ! density to get accurate results) \equiv to # of realizations
            ALLOCATE( COORD(1:SAMPLE_SIZE, 1:D) )
            COORD = 0.D0
            READ(11, *, END = 110)                  ! header: DOMAIN
            READ(11, *, END = 110) UP_LIM, LOW_LIM  ! 1-D DOMAIN BOUNDS
            ! Nodal Coordinates:
            DO I = 1, SAMPLE_SIZE
              READ(11, *, END = 110) COORD(I,:)     ! REALIZATION PTS
            ENDDO
          CASE ('UNSTRUCTURED')
            READ(11, *, END = 110)                  ! header: DATA SITES
            READ(11, *, END = 110) XFILE            ! GRIDFILE NAME
            OPEN(UNIT = 101, FILE = XFILE, STATUS = 'OLD',              &
     &           ACTION = 'READ')
            READ(101, *, END = 510)                  ! header: # OF NODES
            READ(101, *, END = 510) SAMPLE_SIZE      ! # OF NODES
            READ(101, *, END = 510)                  ! header: DIMENSION
            READ(101, *, END = 510) D                ! dimension
            ! Number of scattered data-sites (have to find a reasonable data-site
            ! density to get accurate results) \equiv to # of realizations
            ALLOCATE( COORD(1:SAMPLE_SIZE, 1:D) )
            COORD = 0.D0
            ! Nodal Coordinates:
            READ(101, *, END = 510)                  ! header: COORDS
            DO I = 1, SAMPLE_SIZE
              READ(101, *, END = 510) COORD(I,:)     ! REALIZATION PTS
            ENDDO
  510       CLOSE(UNIT = 101)
          CASE DEFAULT
            PRINT *, 'DATA SITES NOT SPECIFIED!'
            STOP
        END SELECT
        ALLOCATE( X0(1:D), CFACT(1:D) )
        READ(11, *, END = 110)                  ! header: ALS MAXITER
        READ(11, *, END = 110) MAX_ITER         ! # MAX # ALS ITER
        READ(11, *, END = 110)                  ! header: FUNCTION
        READ(11, *, END = 110) FUNC_TYPE_NUM    ! FUNCTION TYPE #
        IF (FUNC_TYPE_NUM .EQ. 4) THEN
           READ(11, *, END = 110)                ! header: center of Gaussian
           DO I = 1, D
             READ(11, *, END = 110) X0(I)        ! CENTER OF GAUSSIAN
           ENDDO
           READ(11, *, END = 110)                ! VECTOR DECAY FACT
           DO I = 1, D
             READ(11, *, END = 110) CFACT(I)     ! DECAY FACTOR - I
           ENDDO
           READ(11, *, END = 110)                ! DECAY FACTOR - II
           READ(11, *, END = 110) ALPHA
        ELSEIF (FUNC_TYPE_NUM .EQ. 7) THEN
          READ(11, *, END = 110)                ! header: WAVE NUMBER
          READ(11, *, END = 110) WAVE_NUM
          READ(11, *, END = 110)                ! header: # OF
                                                ! PLANE-WAVE ENRICHMENTS
          READ(11, *, END = 110) Q
          READ(11, *, END = 110)                ! header: DIR_NUM
          READ(11, *, END = 110) DIR_NUM
        ELSEIF (FUNC_TYPE_NUM .EQ. 8) THEN
          READ(11, *, END = 110)                ! header: WAVE NUMBER
          READ(11, *, END = 110) WAVE_NUM
          READ(11, *, END = 110)                ! header: # OF
                                                ! PLANE-WAVE ENRICHMENTS
          READ(11, *, END = 110) Q
          READ(11, *, END = 110)                ! header: DIR_NUM
          READ(11, *, END = 110) DIR_NUM
        ELSEIF (FUNC_TYPE_NUM .EQ. 9) THEN
          READ(11, *, END = 110)                ! header: DEG of LEG_POLY
          READ(11, *, END = 110) ELL
          READ(11, *, END = 110)                ! header: ORD of LEG_POLY
          READ(11, *, END = 110) EMM
          READ(11, *, END = 110)                ! header: center of Gaussian
          DO I = 1, D
            READ(11, *, END = 110) X0(I)        ! CENTER OF GAUSSIAN
          ENDDO
          READ(11, *, END = 110)                ! header: DECAY FACTOR
          READ(11, *, END = 110) ALPHA
        ELSEIF (FUNC_TYPE_NUM .EQ. 10) THEN
          READ(11, *, END = 110)                ! header: FUNC FILE NAME
          READ(11, *, END = 110) GXFILE         
        ENDIF
        READ(11, *, END = 110)                  ! header: ERROR TOL
        READ(11, *, END = 110) TOL              ! TOLERANCE
        READ(11, *, END = 110)                  ! header: LEN OF POLY ORDER ARRAY
        READ(11, *, END = 110) LEN_M_ARR        ! LENGTH POLY ORDER ARRAY
        READ(11, *, END = 110)                  ! header: POLY ORDER ARRAY
        ALLOCATE( M_ARR(1:LEN_M_ARR) )
        M_ARR = 0
        DO I = 1, LEN_M_ARR
          READ(11, *, END = 110) M_ARR(I)       ! POLY ORDER ARRAY
        ENDDO           
        READ(11, *, END = 110)                  ! header: MAX RANK
        READ(11, *, END = 110) R_MAX            ! MAX RANK
        READ(11, *, END = 110)                  ! header: POLY TYPE
        READ(11, *, END = 110) POLY_TYPE        ! ORTH POLY TYPE
        READ(11, *, END = 110)                  ! header: ORDER OF QUADRATURE
        READ(11, *, END = 110) INT_ORDER        ! ORDER OF QUADRATURE
        READ(11, *, END = 110)                  ! header: CMIN
        READ(11, *, END = 110) C_MIN            ! C_MIN
        READ(11, *, END = 110)                  ! header: CMAX
        READ(11, *, END = 110) C_MAX            ! C_MAX
        READ(11, *, END = 110)                  ! header: LAPACK SOLVER
        READ(11, *, END = 110) SOLVER           ! SOLVER
        READ(11, *, END = 110)                  ! header: PROBLEM TYPE
        READ(11, *, END = 110) PRBLM_TYPE       ! PROBLEM TYPE (LS/NR.EQ.)
        READ(11, *, END = 110)                  ! header: REGULARIZATION STATE
        READ(11, *, END = 110) T_REGU_STAT      ! 'ON' FOR NRML_EQN, 'OFF' FOR LLS
        READ(11, *, END = 110)                  ! header: REGULARIZATION PARAM
        READ(11, *, END = 110) LAMBDA           ! PARAMETER VALUE
        READ(11, *, END = 110)                  ! header: OUT-FILE NAME
        READ(11, *, END = 110) OUTFILENAME      ! OUT-FILE NAME
        READ(11, *, END = 110)                  ! header: RUNTIME DATA FILE
        READ(11, *, END = 110) RUNDATAFILE      ! RUNTIME DATA FILE NAME
        READ(11, *, END = 110)                  ! header: GRIDFILE
        READ(11, *, END = 110) COORDFILE        ! GRIDFILE NAME
        READ(11, *, END = 110)                  ! header: PLOTFILE
        READ(11, *, END = 110) PLOTFILE         ! PLOTFILE NAME
        READ(11, *, END = 110)                  ! header: RESOLUTION
        READ(11, *, END = 110) NFINEDIV
  110   CLOSE(UNIT = 11)

        WRITE(*,*) " "
        WRITE(*,*) "Read problem data file successfully!"
        WRITE(*,*) " "
        ! YET TO FIGURE OUT HOW THIS WORKS:
        !time1 = etime ( t_array )
        call system_clock ( t1, clock_rate, clock_max )

        !WRITE(*,*) "Now writing THE COORDINATES to external data file: "
        !WRITE(*,*) "..."
        !OPEN(UNIT = 12, FILE = COORDFILE, STATUS = 'NEW')
        ! Nodes:
        !WRITE(12, *)"THE COORDINATES: "
        !DO I = 1, SIZE(COORD,1)
        !  WRITE(12, *)"Node = ", I, " Coordinates = ", COORD(I,:)      
        !ENDDO
        !WRITE(12, *) " "
        !WRITE(12, *)"THE POLYNOMIAL ORDER ARRAY: "
        !DO I = 1, LEN_M_ARR
        !  WRITE(12, *)"M_ARR(", I, ") = ", M_ARR(I)
        !ENDDO
        !WRITE(12, *) " "
        !CLOSE(UNIT = 12)
        !WRITE(*,*) "DONE writing THE COORDINATES! "


        !$OMP PARALLEL
        !$ WRITE ( *, '(A)' ) '  Use OpenMP FOR Parallel Execution.'
        !$ WRITE ( *, '(A,I4)' ) &
        !$ 'The number of processors available is ', OMP_GET_NUM_PROCS( )
        !$ ID = OMP_GET_THEREAD_NUM ( )
        !$ IF ( ID .EQ. 0 ) THEN
        !$   WRITE ( *, '(A,I4)' ) &
        !$   'The number of threads available is ', OMP_GET_NUM_THREADS ( )
        !$ ENDIF

        !$ WRITE(*, '(A)', ADVANCE = "NO")'Specify the number of threads: '
        !$ READ(*,*) NUM_THREADS

        ! Allocating number of threads:
        !$ CALL OMP_SET_NUM_THREADS(NUM_THREADS)

        ALLOCATE( U(1:SAMPLE_SIZE) )
        U = 0.D0
        SELECT CASE (FUNC_TYPE_NUM)
          CASE (1)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC1((/ COORD(J,:) /), D)
            ENDDO
           !$OMP ENDDO
          CASE (2)
            ALLOCATE( C_PARAM(1:D), W_PARAM(1:D) )
            C_PARAM = 9.D0/D
            W_PARAM = 1.D0/D
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC2((/ COORD(J,:) /), D, C_PARAM, W_PARAM)
            ENDDO
           !$OMP ENDDO
          CASE (3)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC3((/ COORD(J,:) /), D)
            ENDDO
           !$OMP ENDDO
          CASE (4)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC4((/ COORD(J,:) /), D, ALPHA, CFACT, X0)
            ENDDO
           !$OMP ENDDO
          CASE (5)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC5((/ COORD(J,:) /), D)
            ENDDO
           !$OMP ENDDO
          CASE (6)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC6((/ COORD(J,:) /), D)
            ENDDO
           !$OMP ENDDO
          CASE (7)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC7((/ COORD(J,:) /), D, WAVE_NUM, Q, DIR_NUM)
            ENDDO
           !$OMP ENDDO
          CASE (8)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC8((/ COORD(J,:) /), D, WAVE_NUM, Q, DIR_NUM)
            ENDDO
           !$OMP ENDDO
          CASE (9)
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              U(J) = FUNC9(ELL, EMM, X0, ALPHA, (/ COORD(J,:) /))
            ENDDO
            !$OMP ENDDO
          CASE (10)
            OPEN(UNIT = 102, FILE = GXFILE, STATUS = 'OLD',             &
     &           ACTION = 'READ')
            !$OMP DO
            DO J = 1, SAMPLE_SIZE
              READ(102, *, END = 610) U(J)
            ENDDO
            !$OMP ENDDO
  610       CLOSE(UNIT = 102)
          CASE DEFAULT
            PRINT *, 'NOT CODED YET!'
            STOP
        END SELECT
        !$OMP END PARALLEL

        ! Maximum error between the actual function and the separated rep.
        ! TOL is the user specified value which scales the seminorm
        ! (discrete L2-norm) of the analytical function. So, EPS is, in
        ! a way, a problem (function) specific tolerence which makes the 
        ! code function independent as EPS (the tuned tolerance) gets
        ! used as iput in the ALS subroutines. 
        
        EPS = TOL*SEMI_NORM(U, SAMPLE_SIZE)

        ALLOCATE( ERR_MAT(LEN_M_ARR), R_MAT(R_MAX),                 &
     &            INT_MAT(LEN_M_ARR) )

        ! *********************** IMPORTANT NOTE **********************
        ! This is the the outermost loop of the entire code. The user
        ! inputs a vector M which indicates the allowed polynomial basis
        ! orders to be tried on. Starting from the first order (M(1))
        ! the ALS routine is called to seek the desired separated
        ! representation. If the requirements are met (in the form of
        ! TOL and EPS) for some polynomial order M(J) where 
        ! M(1) <= J <= M(end) then no more polynomials after M(J) is
        ! tried out. 
        ! *************************************************************
        
        I = 0
        DO WHILE (I .LT. LEN_M_ARR)
          COUNTER = COUNTER + 1
          I = I + 1
          M = M_ARR(I)

          WRITE(*,*)' '
          WRITE(*,*)'Subroutine ALS evoked for polynomial order: ', M
          WRITE(*,*)' '

          SELECT CASE (POLY_TYPE)
            CASE ('ORTH_LEGENDRE')
              CALL ALS_ORTH(ERR_DECAY, R_OUT, C, S, D, M, R_MAX, C_MIN, &
     &             C_MAX, POLY_TYPE, U, SAMPLE_SIZE, COORD, EPS,        &
     &             RUNDATAFILE, SOLVER, PRBLM_TYPE, FLAG_CONV,          &
     &             MAX_ITER, T_REGU_STAT, LAMBDA)

             CASE ('ORTH_CHEBYSHEV')
              CALL ALS_ORTH(ERR_DECAY, R_OUT, C, S, D, M, R_MAX, C_MIN, &
     &             C_MAX, POLY_TYPE, U, SAMPLE_SIZE, COORD, EPS,        &
     &             RUNDATAFILE, SOLVER, PRBLM_TYPE, FLAG_CONV,          &
     &             MAX_ITER, T_REGU_STAT, LAMBDA, LOW_LIM, UP_LIM)

            CASE ('BERNSTEIN')
              CALL ALS_BRN(ERR_DECAY, R_OUT, C, S, D, M, R_MAX, C_MIN,  &
     &             C_MAX, U, SAMPLE_SIZE, COORD, EPS, LOW_LIM, UP_LIM,  &
     &             RUNDATAFILE, SOLVER, PRBLM_TYPE, FLAG_CONV,          &
     &             MAX_ITER, T_REGU_STAT, LAMBDA)
            CASE DEFAULT
              PRINT *, ' '
              PRINT *, 'NO OTHER POLYNOMIALS ARE CODED TO SOLVE THE ALS'
              PRINT *, ' '
              STOP
            END SELECT
          ERR_MAT(COUNTER) = ERR_DECAY
          R_MAT(COUNTER) = R_OUT

          IF (FLAG_CONV .EQV. .TRUE.) THEN
            WRITE(*,*)' '
            WRITE(*,*)'The solution converged with Separation rank = ', &
     &                 R_OUT
            WRITE(*,*)' '

            SELECT CASE (POLY_TYPE)
              CASE ('ORTH_LEGENDRE')
                CALL NUM_INTEGRATION(INTEGRAL_OUT, C, S, POLY_TYPE, M,  &
     &               D, R_OUT, INT_ORDER, LOW_LIM, UP_LIM)
              CASE ('ORTH_CHEBYSHEV')
                CALL NUM_INTEGRATION(INTEGRAL_OUT, C, S, POLY_TYPE, M,  &
     &               D, R_OUT, INT_ORDER, LOW_LIM, UP_LIM)
              CASE ('BERNSTEIN')
                CALL NUM_INTEGRATION_BRN(INTEGRAL_OUT, C, S, M, D,      &
     &               R_OUT, INT_ORDER, LOW_LIM, UP_LIM)
              CASE DEFAULT
                PRINT *, ' '
                PRINT *, 'NO OTHER POLYNOMIALS ARE CODED FOR QUADRATURE'
                PRINT *, ' '
                STOP
            END SELECT

            WRITE(*,*)' '
            WRITE(*,*)'Numerical Integration Yields: ', INTEGRAL_OUT
            WRITE(*,*)' '
            
            INT_MAT(COUNTER) = INTEGRAL_OUT

            GO TO 19991

          ELSEIF (FLAG_CONV .EQV. .FALSE.) THEN
            WRITE(*,*)' '
            WRITE(*,*)'No solution is possible for polynomial'
            WRITE(*,*)' approximation of order = ', M
            WRITE(*,*)' '
            
            SELECT CASE (POLY_TYPE)
              CASE ('ORTH_LEGENDRE')
                CALL NUM_INTEGRATION(INTEGRAL_OUT, C, S, POLY_TYPE, M,  &
     &               D, R_OUT, INT_ORDER, LOW_LIM, UP_LIM)
              CASE ('ORTH_CHEBYSHEV')
                CALL NUM_INTEGRATION(INTEGRAL_OUT, C, S, POLY_TYPE, M,  &
     &               D, R_OUT, INT_ORDER, LOW_LIM, UP_LIM)
              CASE ('BERNSTEIN')
                CALL NUM_INTEGRATION_BRN(INTEGRAL_OUT, C, S, M, D,      &
     &               R_OUT, INT_ORDER, LOW_LIM, UP_LIM)
              CASE DEFAULT
                PRINT *, ' '
                PRINT *, 'NO OTHER POLYNOMIALS ARE CODED FOR QUADRATURE'
                PRINT *, ' '
                STOP
            END SELECT

            WRITE(*,*)'Numerical Integration Yields: ', INTEGRAL_OUT
            WRITE(*,*)'RESULT MAY BE INACCURATE DUE TO NON-CONVERGENCE.'
            WRITE(*,*)' '

            INT_MAT(COUNTER) = INTEGRAL_OUT

            IF (I .LE. LEN_M_ARR) THEN
              WRITE(*,*)' '
              WRITE(*,*)'Moving on to the next poly. order in the range'
              WRITE(*,*)' '
            ELSE
              WRITE(*,*)' '
              WRITE(*,*)'NO GOOD SOLUTION FOUND FOR ANY OF THE          &
     &                   POLYNOMIALS IN THE GIVEN RANGE'
              WRITE(*,*)' '
            ENDIF
          ENDIF

        ENDDO
        !       !$OMP ENDDO
        !       !$OMP END PARALLEL
       
19991   CONTINUE

        WRITE(*,*)' '
        WRITE(*,*) "Now writing the output to external data file: "
        WRITE(*,*) "..."
        WRITE(*,*)' '
        OPEN(UNIT = 13, FILE = OUTFILENAME, STATUS = 'NEW')
        ! Nodes:
        WRITE(13, *)"ERROR NORM VALUES: "
        DO I = 1, COUNTER
          WRITE(13, *) "POLYNOMIAL ORDER = ", M_ARR(I), " , ERROR = ",  &
     &          ERR_MAT(I)      
        ENDDO
        WRITE(13, *) " "

        WRITE(13, *)"RANK INFORMATION: "
        DO I = 1, COUNTER
          WRITE(13, *) "POLYNOMIAL ORDER = ", M_ARR(I),                 &
     &        " , RANK = ", R_MAT(I)
        ENDDO

        WRITE(13, *) " "

        DO I = 1, COUNTER
          DO K = 1, D
            WRITE(13, *) "FOR DIMENSION # ", K, ", VALUES OF ALL THE    &
     &                   COEFFS --> ", C(1:R_MAT(I), 1:M_ARR(I), K)
          ENDDO
        ENDDO
        WRITE(13, *) " "
        WRITE(13, *)"NUMERICAL INTEGRATION: "
        DO I = 1, COUNTER
          WRITE(13, *) "POLYNOMIAL ORDER = ", M_ARR(I),                 &
     &      " , CONVERGED/UNCONVERGED RANK = ", R_MAT(I),               &
     &      ", NUMERICAL INTEGRATION YIELDS: ", INT_MAT(I)
        ENDDO
        WRITE(13, *) " "
        CLOSE(UNIT = 13)

        !************************************************************
        !************************************************************

        ! Writing another MATLAB input-file to plot the 
        ! solution if dimension D == 2 or 3. In case of 2-dimensions we
        ! use triangulation in MATLAB. For 3-dimensions iso surface is
        ! used. The format of the output in 30dimensions shouldbe
        ! consitent with the meshgrid function of matlab. 
        
        IF (D .EQ. 2) THEN
          HFINEX = (UP_LIM - LOW_LIM)/NFINEDIV
          HFINEY = HFINEX
          ALLOCATE( XFINE(NFINEDIV+1), YFINE(NFINEDIV+1) )
          ALLOCATE( FINEGRID((NFINEDIV+1)**2,2) )
          XFINE = LOW_LIM + HFINEX*(/ (J, J = 0, NFINEDIV) /)
          YFINE = XFINE

          CALL CARTPROD2D(FINEGRID, XFINE, YFINE, (NFINEDIV+1),         &
     &                    (NFINEDIV+1))

          CALL POSTPROCESSOR(PLOTFILE, NFINEDIV, FINEGRID, C, S,        &
     &                       R_OUT, D, M, POLY_TYPE, LOW_LIM, UP_LIM)

        ELSEIF (D .EQ. 3) THEN
          HFINEX = (UP_LIM - LOW_LIM)/NFINEDIV
          HFINEY = HFINEX
          HFINEZ = HFINEX
          ALLOCATE( XFINE(NFINEDIV+1), YFINE(NFINEDIV+1),               &
     &              ZFINE(NFINEDIV+1) )
          ALLOCATE( FINEGRID((NFINEDIV+1)**3,3) )
          XFINE = LOW_LIM + HFINEX*(/ (J, J = 0, NFINEDIV) /)
          YFINE = XFINE
          ZFINE = XFINE

          CALL CARTPROD3D(FINEGRID, XFINE, YFINE, ZFINE, (NFINEDIV+1),  &
     &                    (NFINEDIV+1), (NFINEDIV+1))

          CALL POSTPROCESSOR(PLOTFILE, NFINEDIV, FINEGRID, C, S,        &
     &                       R_OUT, D, M, POLY_TYPE, LOW_LIM, UP_LIM)

        ENDIF

        IF ((FUNC_TYPE_NUM .EQ. 7) .OR. (FUNC_TYPE_NUM .EQ. 8)) THEN
          CALL COEFF_FILE_HELM(FUNC_TYPE_NUM, WAVE_NUM, Q, DIR_NUM, C, &
     &                    R_OUT, S, D, M)
        ENDIF
          
        WRITE(*,*)' '
        WRITE(*,*)'...'
        WRITE(*,*)' '
        WRITE(*,*) "Done writing the output data files. "
        WRITE(*,*)' '
        WRITE(*,*)'PROGRAM EXECUTION ENDS.'
        WRITE(*,*)' '
        WRITE(*,*)'GOOD BYE!'
        WRITE(*,*)' '

        call system_clock ( t2, clock_rate, clock_max )
        write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) /       &
     &                                      real ( clock_rate )

      END PROGRAM FUNC_DECOM_ALS

      !==================================================================
      !==================================================================

      SUBROUTINE POSTPROCESSOR(PLOTFILE, NFINEDIV, FINEGRID, C, S, R,   &
     &                         D, M, POLY_TYPE, LOW_LIM, UP_LIM)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE SEP_RESP(F_SEP, C, S, Y, POLY_TYPE, R, M, D, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)          :: R, M, D
            CHARACTER(80), INTENT(IN)         :: POLY_TYPE
            REAL(DP), INTENT(IN)              :: C(1:R, 1:(M+1), 1:D)
            REAL(DP), INTENT(IN)              :: S(1:R), Y(1:D)
            REAL(DP), INTENT(IN)              :: A, B
            REAL(DP), INTENT(OUT)             :: F_SEP          
          END SUBROUTINE SEP_RESP
        END INTERFACE
        
        INTERFACE
          SUBROUTINE SEP_RESP_BRN(F_SEP, C, S, Y, R, M, D, LOW_LIM,     &
     &                            UP_LIM)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)          :: R, M, D
            REAL(DP), INTENT(IN)              :: C(1:R, 1:(M+1), 1:D)
            REAL(DP), INTENT(IN)              :: S(1:R), Y(1:D)
            REAL(DP), INTENT(IN)              :: LOW_LIM, UP_LIM
            REAL(DP), INTENT(OUT)             :: F_SEP
          END SUBROUTINE SEP_RESP_BRN
        END INTERFACE

        INTEGER(I4B), INTENT(IN)       :: D, M, NFINEDIV, R
        CHARACTER(80), INTENT(IN)      :: POLY_TYPE
        REAL(DP), INTENT(IN)           :: FINEGRID(1:(NFINEDIV+1)**D,D)
        REAL(DP), INTENT(IN)           :: C(1:R, 1:(M+1), 1:D)
        REAL(DP), INTENT(IN)           :: S(1:R)
        REAL(DP), INTENT(IN)           :: LOW_LIM, UP_LIM
        CHARACTER(80), INTENT(OUT)     :: PLOTFILE

        ! LOCAL VARIABLES:
        INTEGER(I4B)            :: J
        REAL(DP)                :: F_SEP
        REAL(DP)                :: SAMPLE_POINT(1:D)
        
        OPEN(UNIT = 31, FILE = PLOTFILE, STATUS = 'NEW')          

        !$OMP PARALLEL
        !$OMP DO
        DO J = 1, (NFINEDIV+1)**D

          SAMPLE_POINT = FINEGRID(J,:)

          IF (POLY_TYPE .EQ. 'ORTH_LEGENDRE') THEN
              CALL SEP_RESP(F_SEP, C, S, SAMPLE_POINT, POLY_TYPE, R,    &
     &                      M, D, LOW_LIM, UP_LIM)
            ELSEIF (POLY_TYPE .EQ. 'ORTH_CHEBYSHEV') THEN
              CALL SEP_RESP(F_SEP, C, S, SAMPLE_POINT, POLY_TYPE, R,    &
     &                      M, D, LOW_LIM, UP_LIM)
            ELSEIF (POLY_TYPE .EQ. 'BERNSTEIN') THEN
              CALL SEP_RESP_BRN(F_SEP, C, S, SAMPLE_POINT, R,           &
     &                          M, D, LOW_LIM, UP_LIM)
          ENDIF
            
          WRITE(31,*) SAMPLE_POINT, F_SEP

        ENDDO
        !$OMP ENDDO        
        !$OMP END PARALLEL
        CLOSE(UNIT = 31)

      END SUBROUTINE POSTPROCESSOR

      !==================================================================
      !==================================================================
      
      SUBROUTINE COEFF_FILE_HELM(FUNC_TYPE_NUM,  WAVE_NUM, Q, DIR_NUM,  &
     &                           C, R, S, D, M)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          CHARACTER(LEN=1024) FUNCTION STRINT(K)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)          :: K
          END FUNCTION STRINT
        END INTERFACE

        INTERFACE
          CHARACTER(LEN=1024) FUNCTION STRFLT(K)
            USE NRTYPE
            REAL(DP), INTENT(IN)          :: K
          END FUNCTION STRFLT
        END INTERFACE

        INTEGER(I4B), INTENT(IN)              :: FUNC_TYPE_NUM, Q, R, M
        INTEGER(I4B), INTENT(IN)              :: DIR_NUM, D
        REAL(DP), INTENT(IN)                  :: WAVE_NUM
        REAL(DP), INTENT(IN)                  :: C(1:R, 1:(M+1), 1:D)
        REAL(DP), INTENT(IN)                  :: S(1:R)

        ! LOCAL VARIABLES:
        INTEGER(I4B)             :: I, J, L, DIMEN, FILE_ID
        CHARACTER(1024)          :: FILENAME

        IF (FUNC_TYPE_NUM .EQ. 7) THEN
          FILENAME = 'COS_COEFFS_Q_'//STRINT(Q)//'_K_'//                &
     &               STRFLT(WAVE_NUM)//'_DIR_NUM_'//STRINT(DIR_NUM)//   &
     &               '.dat'
        ELSEIF (FUNC_TYPE_NUM .EQ. 8) THEN
          FILENAME = 'SIN_COEFFS_Q_'//STRINT(Q)//'_K_'//                &
     &               STRFLT(WAVE_NUM)//'_DIR_NUM_'//STRINT(DIR_NUM)//   &
     &               '.dat'
        ENDIF

        !FILENAME = ADJUSTL(FILENAME)
        FILE_ID = 100*DIR_NUM+101

        OPEN(UNIT = FILE_ID, FILE = FILENAME, STATUS = 'NEW')
        WRITE(FILE_ID, *) 'RANK:'
        WRITE(FILE_ID, *) R
        WRITE(FILE_ID, *) 'MAX_ORDER_OF_POLYNOMIAL:'
        WRITE(FILE_ID, *) M
        WRITE(FILE_ID, *) 'DIMENSION:'
        WRITE(FILE_ID, *) D
        WRITE(FILE_ID, *) 'FACTORS:'
        !$OMP PARALLEL
        !$OMP DO
        DO I = 1, R
          WRITE(FILE_ID, *) S(I)
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL

        !$OMP PARALLEL
        !$OMP DO
        DO DIMEN = 1, D
          !$OMP DO
          DO L = 1, R
            WRITE(FILE_ID, *) ( C(L, J, DIMEN), J = 1, M+1)
          ENDDO
          !$OMP ENDDO
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL

        CLOSE(UNIT = FILE_ID)

       END SUBROUTINE COEFF_FILE_HELM

      !==================================================================
      !==================================================================

      character(len=1024) function strint(k)
      !   "Convert an integer to string."
      use nrtype
      implicit none
      integer(i4b), intent(in) :: k
      write (strint, *) k
      strint = adjustl(strint)
      end function strint

      !==================================================================
      !==================================================================

      character(len=1024) function strflt(k)
      !   "Convert a float to string."
      use nrtype
      implicit none
      real(dp), intent(in)     :: k
      write (strflt, *) k
      strflt = adjustl(strflt)
      end function strflt


      !==================================================================
      !==================================================================

      SUBROUTINE CARTPROD2D(XY, X, Y, M, N)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN)        :: M, N
        REAL(DP), INTENT(IN)            :: X(1:M), Y(1:N)
        REAL(DP), INTENT(OUT)           :: XY(M*N,2)
        
        ! LOCAL VARIABLES
        INTEGER(I4B)        :: I, J, COUNTER

        COUNTER = 0

        !$OMP PARALLEL SHARED ( I, J )
        !$OMP DO 
        DO J = 1, N
          !$OMP DO 
          DO I = 1, M
            COUNTER = COUNTER + 1
            XY(COUNTER,:) = (/ X(I), Y(J) /)
          ENDDO
          !$OMP ENDDO
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL
      END SUBROUTINE CARTPROD2D

      !==================================================================
      !==================================================================

      SUBROUTINE CARTPROD3D(XYZ, X, Y, Z, M, N, P)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN)        :: M, N, P
        REAL(DP), INTENT(IN)            :: X(1:M), Y(1:N), Z(1:P)
        REAL(DP), INTENT(OUT)           :: XYZ(M*N*P,3)
        ! LOCAL VARIABLES
        INTEGER(I4B)        :: I, J, K, COUNTER
        
        COUNTER = 0

        !$OMP PARALLEL SHARED ( I, J, K )
        !$OMP DO
        DO K = 1, P
          !$OMP DO 
          DO J = 1, N
            !$OMP DO 
            DO I = 1, M
              COUNTER = COUNTER + 1
              XYZ(COUNTER,:) = (/ X(I), Y(J), Z(K) /)
            ENDDO
            !$OMP ENDDO
          ENDDO
          !$OMP ENDDO
        ENDDO
        !$OMP ENDDO

        !$OMP END PARALLEL
      END SUBROUTINE CARTPROD3D

      !==================================================================
      !==================================================================
      
      SUBROUTINE DIAG(DIAG_MAT, N)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN)        :: N
        REAL(DP), INTENT(OUT)           :: DIAG_MAT(N,N)
        ! LOCAL VARIABLES:
        INTEGER(I4B)            :: I, J
        DIAG_MAT = 0.D0
        
        DO J = 1, N
          DO I = 1, N
            DIAG_MAT(I,J) = (I/J)*(J/I)
          ENDDO
        ENDDO
      
      END SUBROUTINE DIAG

      !==================================================================
      !==================================================================

      SUBROUTINE KRON(A_KRON_B, A, M, N, B, P, Q)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN)        :: M, N, P, Q
        REAL(DP), INTENT(IN)            :: A(M,N), B(P,Q)
        REAL(DP), INTENT(OUT)           :: A_KRON_B(M*P,N*Q)
        ! LOCAL VARIABLES:
        INTEGER(I4B)                    :: ROW_BLOCK, COL_BLOCK
        !$OMP PARALLEL
        !$OMP DO
        DO ROW_BLOCK = 1, M
          !$OMP DO
          DO COL_BLOCK = 1, N
            A_KRON_B((ROW_BLOCK-1)*P+1:ROW_BLOCK*P, (COL_BLOCK-1)*Q+1:  &
     &      COL_BLOCK*Q) = A(ROW_BLOCK, COL_BLOCK)*B
          ENDDO
          !$OMP ENDDO
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL
      END SUBROUTINE KRON

      !==================================================================
      !==================================================================

      SUBROUTINE SEP_RESP(F_SEP, C, S, Y, POLY_TYPE, R, M, D, A, B)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE U_HAT_1D(U_HAT, C, M, POLY_TYPE, X, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: M
            REAL(DP), INTENT(IN)                :: C(1:(M+1))
            CHARACTER(80), INTENT(IN)           :: POLY_TYPE
            REAL(DP), INTENT(IN)                :: X(1)
            REAL(DP), INTENT(IN)                :: A, B
            REAL(DP), INTENT(OUT)               :: U_HAT
          END SUBROUTINE U_HAT_1D
        END INTERFACE

        INTEGER(I4B), INTENT(IN)                :: R, M, D
        CHARACTER(80), INTENT(IN)               :: POLY_TYPE
        REAL(DP), INTENT(IN)                    :: C(1:R, 1:(M+1), 1:D)
        REAL(DP), INTENT(IN)                    :: S(1:R), Y(1:D)
        REAL(DP), INTENT(IN)                    :: A, B
        REAL(DP), INTENT(OUT)                   :: F_SEP
        ! LOCAL VARIABLES:
        REAL(DP), ALLOCATABLE :: U_SEP(:,:)
        REAL(DP), ALLOCATABLE :: C_TEMP(:)
        INTEGER(I4B)          :: K, L

        ALLOCATE(  U_SEP(R,D) )
        U_SEP = 0.D0

        !$OMP PARALLEL
        !$OMP DO 
        DO K = 1, D
          !$OMP DO 
          DO L = 1, R
            ALLOCATE( C_TEMP(1:(M+1)) )
            C_TEMP = 0.D0
            C_TEMP = RESHAPE(C(L,:,K), (/ (M+1) /))
            CALL U_HAT_1D(U_SEP(L,K), C_TEMP, M, POLY_TYPE, (/ Y(K) /), &
     &                    A, B)
            DEALLOCATE( C_TEMP )
          ENDDO
          !$OMP ENDDO
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL

        F_SEP = 0.D0
        DO L = 1,R
          F_SEP = F_SEP + S(L)*PRODUCT(U_SEP(L,:))
        ENDDO

      END SUBROUTINE SEP_RESP

      !==================================================================
      !==================================================================

      SUBROUTINE SEP_RESP_BRN(F_SEP, C, S, Y, R, M, D, LOW_LIM, UP_LIM)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE U_HAT_1D_BRN(U_HAT, C, N, X, LOW_LIM, UP_LIM) 
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: N        ! ORDER OF POLYNOMIAL
            REAL(DP), INTENT(IN)                :: C(1:(N+1))
            REAL(DP), INTENT(IN)                :: X
            REAL(DP), INTENT(IN)                :: LOW_LIM, UP_LIM
            REAL(DP), INTENT(OUT)               :: U_HAT
          END SUBROUTINE U_HAT_1D_BRN
        END INTERFACE

        INTEGER(I4B), INTENT(IN)                :: R, M, D
        REAL(DP), INTENT(IN)                    :: C(1:R, 1:(M+1), 1:D)
        REAL(DP), INTENT(IN)                    :: S(1:R), Y(1:D)
        REAL(DP), INTENT(IN)                    :: LOW_LIM, UP_LIM
        REAL(DP), INTENT(OUT)                   :: F_SEP
        
        ! LOCAL VARIABLES:
        REAL(DP), ALLOCATABLE :: U_SEP(:,:)
        REAL(DP), ALLOCATABLE :: C_TEMP(:)
        INTEGER(I4B)          :: K, L

        ALLOCATE(  U_SEP(R,D) )
        U_SEP = 0.D0

        !$OMP PARALLEL
        !$OMP DO 
        DO K = 1, D
          !$OMP DO 
          DO L = 1, R
            ALLOCATE( C_TEMP(1:(M+1)) )
            C_TEMP = RESHAPE(C(L,:,K), (/ (M+1) /))
            CALL U_HAT_1D_BRN(U_SEP(L,K), C_TEMP, M, Y(K), LOW_LIM,     &
     &                        UP_LIM)
            DEALLOCATE( C_TEMP )
          ENDDO
          !$OMP ENDDO
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL

        F_SEP = 0.D0
        DO L = 1,R
          F_SEP = F_SEP + S(L)*PRODUCT(U_SEP(L,:))
        ENDDO

      END SUBROUTINE SEP_RESP_BRN
      
      !==================================================================
      !==================================================================

      SUBROUTINE LINSPACE(X, A, B, N)
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP), INTENT(IN)                    :: A, B 
        INTEGER(I4B), INTENT(IN)                :: N 
        REAL(DP), INTENT(OUT), DIMENSION(1:N)   :: X
        ! LOCAL VARIBALE:
        REAL(DP)        :: DELT
        INTEGER(I4B)    :: I 
        DELT = (B - A)/(N - 1)
        X = A + DELT*(/ (I, I=0, N-1) /) 
      END SUBROUTINE LINSPACE

      !==================================================================
      !==================================================================

      SUBROUTINE TNSR_UHAT_D_DIM(PROD_U_HAT, NODES, M, D, C_RED,        &
     &                           POLY_TYPE, A, B)
        ! C_red[(M + 1), (d - 1)] --> contains all the coefficient at a given rank
        ! except the coeffiecients which go to the opt-loop
        ! dim_arr[d - 1] --> contains all the dimension except the one in opt-loop
        ! nodes[d -1] --> At a given sample point, it contains all the
        ! coordinate component of the sample point except the one direction
        ! that goes to the opt-loop.
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE U_HAT_1D(U_HAT, C, M, POLY_TYPE, X, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: M
            REAL(DP), INTENT(IN)                :: C(1:(M+1))
            CHARACTER(80), INTENT(IN)           :: POLY_TYPE
            REAL(DP), INTENT(IN)                :: X(1)
            REAL(DP), INTENT(OUT)               :: U_HAT
            REAL(DP), INTENT(IN)                :: A, B
          END SUBROUTINE U_HAT_1D
        END INTERFACE
        INTEGER(I4B), INTENT(IN)        :: M, D
        REAL(DP), INTENT(IN)            :: C_RED(1:(M+1),1:(D-1))
        REAL(DP), INTENT(IN)            :: NODES(1:(D-1))
        CHARACTER(80), INTENT(IN)       :: POLY_TYPE
        REAL(DP), INTENT(IN)            :: A, B
        REAL(DP), INTENT(OUT)           :: PROD_U_HAT

        ! LOCAL VARIABLES:
        INTEGER(I4B)        :: I
        REAL(DP)            :: PROD_U_HAT_VEC(1:(D-1))

        !$OMP PARALLEL
        !$OMP DO 
        DO I = 1,(D-1)
          CALL U_HAT_1D(PROD_U_HAT_VEC(I), C_RED(:,I), M, POLY_TYPE,    &
     &                 (/ NODES(I) /), A, B)
        ENDDO
        !$OMP ENDDO
        PROD_U_HAT = PRODUCT(PROD_U_HAT_VEC)

        !$OMP END PARALLEL

      END SUBROUTINE TNSR_UHAT_D_DIM

      !==================================================================
      !==================================================================

      SUBROUTINE U_HAT_1D(U_HAT, C, M, POLY_TYPE, X, A, B)
        ! c[M+1,1]
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTERFACE
          SUBROUTINE ORTH_POLY(PNX, POLY_TYPE, N, X, SIZE_OF_X, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
            INTEGER(I4B), INTENT(IN)        :: SIZE_OF_X       ! # OF PNTS TO EVAL       
            REAL(DP), INTENT(IN)            :: X(1:SIZE_OF_X)  ! PNTS TO EVAL
            CHARACTER(80), INTENT(IN)       :: POLY_TYPE
            REAL(DP), INTENT(OUT)           :: PNX(1:SIZE_OF_X)
            REAL(DP), OPTIONAL, INTENT(IN)  :: A
            REAL(DP), OPTIONAL, INTENT(IN)  :: B
          END SUBROUTINE ORTH_POLY
        END INTERFACE
 
        INTEGER(I4B), INTENT(IN)            :: M
        REAL(DP), INTENT(IN)                :: C(1:(M+1))
        CHARACTER(80), INTENT(IN)           :: POLY_TYPE
        REAL(DP), INTENT(IN)                :: X(1)
        REAL(DP), INTENT(IN)                :: A, B
        REAL(DP), INTENT(OUT)               :: U_HAT
        ! LOCAL VARIABLES:
        INTEGER(I4B)        :: I
        REAL(DP)            :: U_HAT_VEC(1:(M+1))
        !$OMP PARALLEL
        !$OMP DO 
        DO I = 1, (M+1)
          SELECT CASE (POLY_TYPE)
            CASE ('ORTH_LEGENDRE')
              CALL ORTH_POLY(U_HAT_VEC(I), POLY_TYPE, I-1, X(1), 1)
            CASE ('ORTH_CHEBYSHEV')
              CALL ORTH_POLY(U_HAT_VEC(I), POLY_TYPE, I-1, X(1), 1, A,  &
     &                        B)
          END SELECT
              
          !PRINT *, 'U_HAT FOR ORDER = '
          !PRINT *, I
          !PRINT *, 'IS = '
          !PRINT *, U_HAT_VEC(I)
          !PAUSE
        ENDDO
        !$OMP ENDDO
        U_HAT = DOT_PRODUCT(C, U_HAT_VEC)

        !$OMP END PARALLEL
      END SUBROUTINE U_HAT_1D

      !==================================================================
      !==================================================================
      
      SUBROUTINE ORTH_POLY(PNX, POLY_TYPE, N, X, SIZE_OF_X, A, B)
        ! Evaluates Legendre or Hermite Polynomial
        ! for a set of SIZE_OF_X points input as the array X, 
        ! this routine spits-out the value of the orthogonal 
        ! polynomials of type POLY_TYPE and stores the result in
        ! the array PNX of size SIZE_OF_X. 
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! The recursive method is too slow
        !INTERFACE 
        !  RECURSIVE FUNCTION LEGENDRE(X, N)  RESULT(L_NX)
        !    USE NRTYPE
        !    REAL(DP)                    :: L_NX
        !    INTEGER(I4B), INTENT(IN)    :: N
        !    REAL(DP), INTENT(IN)        :: X
        !  END FUNCTION LEGENDRE
        !END INTERFACE

        ! The recursive method is too slow
        !INTERFACE
        !  RECURSIVE FUNCTION HERMITE(X, N)  RESULT(H_NX)
        !    USE NRTYPE
        !    REAL(DP)                    :: H_NX
        !    INTEGER(I4B), INTENT(IN)    :: N
        !    REAL(DP), INTENT(IN)        :: X
        !  END FUNCTION HERMITE
        !END INTERFACE

        INTERFACE
          RECURSIVE FUNCTION FACTORIAL(N)  RESULT(FACT)
            USE NRTYPE
            INTEGER(I4B)                :: FACT
            INTEGER(I4B), INTENT(IN)    :: N
          END FUNCTION FACTORIAL
        END INTERFACE

        INTERFACE
          FUNCTION NCHOOSEK(N, K)
            USE NRTYPE
            INTEGER(I4B)                    :: NCHOOSEK
            INTEGER(I4B), INTENT(IN)        :: N, K
          END FUNCTION NCHOOSEK
        END INTERFACE

        INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
        INTEGER(I4B), INTENT(IN)        :: SIZE_OF_X       ! # OF PNTS TO EVAL       
        REAL(DP), INTENT(IN)            :: X(1:SIZE_OF_X)  ! PNTS TO EVAL
        CHARACTER(80), INTENT(IN)       :: POLY_TYPE
        REAL(DP), INTENT(OUT)           :: PNX(1:SIZE_OF_X)
        REAL(DP), OPTIONAL, INTENT(IN)  :: A
        REAL(DP), OPTIONAL, INTENT(IN)  :: B
        ! LOCAL VARIABLES
        REAL(DP), ALLOCATABLE       :: SUMMAND(:,:)
        INTEGER(I4B)                :: I, M, K
        REAL(DP)                    :: Y(1:SIZE_OF_X)

        SELECT CASE (POLY_TYPE)
          CASE ('ORTH_HERMITE')
            ! PHYSICIST'S HERMITE POLYNOMIAL
            ALLOCATE( SUMMAND(1:SIZE_OF_X, 1:(FLOOR(N/2.D0) + 1) ))
            !$OMP PARALLEL
            !$OMP DO
            DO I = 1, SIZE_OF_X
              !$OMP DO 
              DO M = 0, FLOOR(N/2.D0)
                SUMMAND(I,M+1) = SUMMAND(I,M+1) + (-1)**M*(2.D0*X(i))** &
     &          (N - 2*M)/(FACTORIAL(M)*FACTORIAL(N - 2*M))
              ENDDO
              !$OMP ENDDO
              PNX(I) = SUM(SUMMAND(I,:))
              ! The recursive method is too slow
              !PNX(I) = HERMITE(X(I), N)
            ENDDO
            !$OMP ENDDO
            !$OMP END PARALLEL
            PNX = FACTORIAL(N)*PNX
        
          CASE ('ORTH_LEGENDRE')
            ALLOCATE( SUMMAND(SIZE_OF_X, (N + 1)) )
            !$OMP PARALLEL
            !$OMP DO
            DO I = 1, SIZE_OF_X
              !$OMP DO 
              DO K = 0, N
              !  ! FORMULA -- 1:
                SUMMAND(I,K+1) = ((NCHOOSEK(N, K))**2)*                 &
     &                        ((X(I) - 1.D0)**(N - K))*(X(I) + 1.D0)**K
              ENDDO
              !$OMP ENDDO
              PNX(I) = SUM(SUMMAND(I,:))
              ! The recursive method is too slow
              !PNX(I) = LEGENDRE(X(I), N)
            ENDDO
            !$OMP ENDDO
            !$OMP END PARALLEL
            PNX = 2.D0**(-N)*PNX

          CASE ('ORTH_CHEBYSHEV')
            ! SCALED INTERVAL:
            Y = (X - (A + B)/2.D0)/(05.D0*(B - A))
            !$OMP PARALLEL
            !$OMP DO
            DO I = 1, SIZE_OF_X
              PNX(I) = COS(N*ACOS(Y(I)))
            ENDDO
            !$OMP ENDDO
            !$OMP END PARALLEL

          CASE DEFAULT
            PRINT *, 'NOT CODED YET!'
            STOP
        END SELECT
        RETURN
      END SUBROUTINE ORTH_POLY

      !==================================================================
      !==================================================================

      RECURSIVE FUNCTION LEGENDRE(X, N)  RESULT(L_NX)
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: L_NX
        INTEGER(I4B), INTENT(IN)    :: N
        REAL(DP), INTENT(IN)        :: X
        IF (N .EQ. 0) THEN
          L_NX = 1.D0
        ELSEIF (N .EQ. 1) THEN
          L_NX = X
        ELSE
          L_NX = (2.D0*N - 1.D0)/N * X * LEGENDRE(X, N-1) -             &
     &           (N - 1.D0)/N *  LEGENDRE(X, N-2)
        ENDIF
      END FUNCTION LEGENDRE

      !==================================================================
      !==================================================================

      RECURSIVE FUNCTION HERMITE(X, N)  RESULT(H_NX)
        ! PHYSICIST'S DEFINITION:
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: H_NX
        INTEGER(I4B), INTENT(IN)    :: N
        REAL(DP), INTENT(IN)        :: X
        IF (N .EQ. 0) THEN
          H_NX = 1.D0
        ELSEIF (N .EQ. 1) THEN
          H_NX = 2.D0 * X
        ELSE
          H_NX = 2.D0 * X * HERMITE(X, N-1) -                           &
     &           2.D0 * (N - 1.D0) *  HERMITE(X, N-2)
        ENDIF
      END FUNCTION HERMITE

      !==================================================================
      !==================================================================


      SUBROUTINE TNSR_UHAT_D_DIM_BRN(PROD_U_HAT, NODES, M, D, C_RED,    &
     &                               LOW_LIM, UP_LIM)
        ! C_red[(M + 1), (d - 1)] --> contains all the coefficient at a given rank
        ! except the coeffiecients which go to the opt-loop
        ! dim_arr[d - 1] --> contains all the dimension except the one in opt-loop
        ! nodes[d -1] --> At a given sample point, it contains all the
        ! coordinate component of the sample point except the one direction
        ! that goes to the opt-loop.
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE U_HAT_1D_BRN(U_HAT, C, N, X, LOW_LIM, UP_LIM) 
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: N        ! ORDER OF POLYNOMIAL
            REAL(DP), INTENT(IN)                :: C(1:(N+1))
            REAL(DP), INTENT(IN)                :: X
            REAL(DP), INTENT(IN)                :: LOW_LIM, UP_LIM
            REAL(DP), INTENT(OUT)               :: U_HAT
          END SUBROUTINE U_HAT_1D_BRN
        END INTERFACE

        INTEGER(I4B), INTENT(IN)        :: M, D
        REAL(DP), INTENT(IN)            :: C_RED(1:(M+1),1:(D-1))
        REAL(DP), INTENT(IN)            :: NODES(1:(D-1))
        REAL(DP), INTENT(IN)            :: LOW_LIM, UP_LIM
        REAL(DP), INTENT(OUT)           :: PROD_U_HAT

        ! LOCAL VARIABLES:
        INTEGER(I4B)        :: I
        REAL(DP)            :: PROD_U_HAT_VEC(1:(D-1))

        !$OMP PARALLEL
        !$OMP DO 
        DO I = 1,(D-1)
          CALL U_HAT_1D_BRN(PROD_U_HAT_VEC(I), C_RED(:,I), M, NODES(I), &
     &                      LOW_LIM, UP_LIM)
        ENDDO
        !$OMP ENDDO
        PROD_U_HAT = PRODUCT(PROD_U_HAT_VEC)

        !$OMP END PARALLEL

      END SUBROUTINE TNSR_UHAT_D_DIM_BRN

      !==================================================================
      !==================================================================

      SUBROUTINE U_HAT_1D_BRN(U_HAT, C, N, X, LOW_LIM, UP_LIM)
        ! c[M+1,1]
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTERFACE
          SUBROUTINE BERNSTEIN(BRN_BASIS_POLY, N, A, B, X)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
            REAL(DP), INTENT(IN)            :: X               ! POINT TO EVAL
            REAL(DP), INTENT(IN)            :: A, B            ! DOMAIN OF THE PROBLEM
                                                               ! X \IN [A, B]
            REAL(DP), INTENT(OUT)           :: BRN_BASIS_POLY(0:N)
          END SUBROUTINE BERNSTEIN
        END INTERFACE
 
        INTEGER(I4B), INTENT(IN)            :: N        ! ORDER OF POLYNOMIAL
        REAL(DP), INTENT(IN)                :: C(1:(N+1))
        REAL(DP), INTENT(IN)                :: X
        REAL(DP), INTENT(IN)                :: LOW_LIM, UP_LIM
        REAL(DP), INTENT(OUT)               :: U_HAT

        ! LOCAL VARIABLES:
        REAL(DP)            :: BRN_BASIS_POLY(0:N)

        CALL BERNSTEIN(BRN_BASIS_POLY, N, LOW_LIM, UP_LIM, X)

        U_HAT = DOT_PRODUCT(C, BRN_BASIS_POLY)

      END SUBROUTINE U_HAT_1D_BRN

      !==================================================================
      !==================================================================

        SUBROUTINE BERNSTEIN(BRN_BASIS_POLY, N, A, B, X)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        INTERFACE
          FUNCTION NCHOOSEK(N, K)
            USE NRTYPE
            INTEGER(I4B)                    :: NCHOOSEK
            INTEGER(I4B), INTENT(IN)        :: N, K
          END FUNCTION NCHOOSEK
        END INTERFACE

        INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
        REAL(DP), INTENT(IN)            :: X               ! POINT TO EVAL
        REAL(DP), INTENT(IN)            :: A, B            ! DOMAIN OF THE PROBLEM
                                                           ! X \IN [A, B]
        REAL(DP), INTENT(OUT)           :: BRN_BASIS_POLY(0:N)

        ! LOCAL VARIABLES
        INTEGER(I4B)                :: J

        !$OMP PARALLEL
        !$OMP DO
        DO J = 0, N
          BRN_BASIS_POLY(J) = NCHOOSEK(N, J)*(X - A)**J*                &
     &                        (B - X)**(N - J)/(B - A)**N
        ENDDO
        !$OMP ENDDO
        !$OMP END PARALLEL

        END SUBROUTINE BERNSTEIN
        
      !==================================================================
      !==================================================================
     
      RECURSIVE FUNCTION FACTORIAL(N)  RESULT(FACT)
        USE NRTYPE
        IMPLICIT NONE
        INTEGER(I4B)                :: FACT
        INTEGER(I4B), INTENT(IN)    :: N
        IF (N .EQ. 0) THEN
          FACT = 1
        ELSE
          FACT = N*FACTORIAL(N - 1)
        ENDIF
      END FUNCTION FACTORIAL

      !==================================================================
      !==================================================================

      FUNCTION NCHOOSEK(N, K)
        USE NRTYPE
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          RECURSIVE FUNCTION FACTORIAL(N)  RESULT(FACT)
            USE NRTYPE
            INTEGER(I4B)                :: FACT
            INTEGER(I4B), INTENT(IN)    :: N
          END FUNCTION FACTORIAL
        END INTERFACE
        INTEGER(I4B)                    :: NCHOOSEK
        INTEGER(I4B), INTENT(IN)        :: N, K
        NCHOOSEK = FACTORIAL(N)/FACTORIAL(K)/FACTORIAL(N - K)
      END FUNCTION NCHOOSEK

      !==================================================================
      !==================================================================

      FUNCTION SEMI_NORM(U, N)
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                        :: SEMI_NORM
        INTEGER(I4B), INTENT(IN)        :: N
        REAL(DP), INTENT(IN)            :: U(1:N)
        SEMI_NORM = SQRT((DOT_PRODUCT(U,U))/N)
      END FUNCTION SEMI_NORM

      !==================================================================
      !==================================================================
      
      FUNCTION FUNC1(X, D)
        ! FUNC1 = SIN(X**2 + Y**2 + ...)
        ! ACCEPTS ANY VECTOR OF DIMENSION D
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC1
        INTEGER(I4B), INTENT(IN)    :: D
        REAL(DP), INTENT(IN)        :: X(1:D)
        ! TEMP VARIABLES:
        REAL(DP)        :: ARG(1:D)
        INTEGER(I4B)    :: I
        ARG = (/ (X(I)**2, I = 1, D) /)
        FUNC1 = SIN(SUM(ARG))
      END FUNCTION FUNC1

      !==================================================================
      !==================================================================

      ! Ref: High dimensional polynomial interpolation on sparse grids
      ! Volker Barthelmann, Erich Novak, and Klaus Ritter
      ! Advances in Computational Mathematics 12 (2000) 273–288
      ! Page 284, function # 1:
      !
      ! Oscillatory
      FUNCTION FUNC2(X, D, C, W)
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                            :: FUNC2
        INTEGER(I4B), INTENT(IN)            :: D
        REAL(DP), INTENT(IN)                :: C(1:D), W(1:D)
        REAL(DP), INTENT(IN)                :: X(1:D)
        ! TEMP VARIABLES:
        FUNC2 = COS(2.D0*PI_D*W(1) + DOT_PRODUCT(C, X))

      END FUNCTION FUNC2

      !==================================================================
      !==================================================================
            
      FUNCTION FUNC3(X, D)
        ! FUNC3 = COS(X**2 + Y**2 + ...)
        ! ACCEPTS ANY VECTOR OF DIMENSION D
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC3
        INTEGER(I4B), INTENT(IN)    :: D
        REAL(DP), INTENT(IN)        :: X(1:D)
        ! TEMP VARIABLES:
        REAL(DP)        :: ARG(1:D)
        INTEGER(I4B)    :: I
        ARG = (/ (X(I)**2, I = 1, D) /)
        FUNC3 = COS(SUM(ARG))
      END FUNCTION FUNC3

      !==================================================================
      !==================================================================

      ! Ref: High dimensional polynomial interpolation on sparse grids
      ! Volker Barthelmann, Erich Novak, and Klaus Ritter
      ! Advances in Computational Mathematics 12 (2000) 273–288
      ! Page 284, function # 4:
      
      ! Gaussian:
      ! FUNC4 = exp(-\sum(c_i**2*t*(x_i-(x0)_i))) summed over all the dimensions
      ! The higher the 
      ! The parameter x0 acts as a shift parameter, and the difficulty 
      ! of the functions is monotonically increasing with the ci > 0.
      FUNCTION FUNC4(X, D, T, C, X0)     
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC4
        INTEGER(I4B), INTENT(IN)    :: D
        REAL(DP), INTENT(IN)        :: X(1:D)
        REAL(DP), INTENT(IN)        :: C(1:D), X0(1:D)
        REAL(DP), INTENT(IN)        :: T

        ! TEMP VARIABLES:
        REAL(DP)        :: ARG(1:D)
        INTEGER(I4B)    :: I
!        REAL(DP)        :: C(1:D), W(1:D)
!        REAL(DP)        :: T = 0.5D0
!        C = SQRT(2.D0)
!        W = 1.D0/8.D0
        ARG = (/ ( C(I)**2*T*(X(I) - X0(I))**2, I = 1, D) /)
        FUNC4 = EXP(-SUM(ARG))
      END FUNCTION FUNC4

      !==================================================================
      !==================================================================

      FUNCTION FUNC5(X, D)
        ! FUNC5 := PRODUCT PEAK
        ! ACCEPTS ANY VECTOR OF DIMENSION D
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC5
        INTEGER(I4B), INTENT(IN)    :: D
        REAL(DP), INTENT(IN)        :: X(1:D)
        ! TEMP VARIABLES:
        REAL(DP)        :: ARG(1:D)
        INTEGER(I4B)    :: I
        REAL(DP)        :: C(1:D), W(1:D)
        C = 7.25D0/D
        W = 1.D0/D
        ARG = (/ ( (C(I)**(-2) + (X(I) - W(I))**2)**(-1), I = 1, D) /)
        FUNC5 = PRODUCT(ARG)
      END FUNCTION FUNC5

      !==================================================================
      !==================================================================

      FUNCTION FUNC6(X, D)
        ! Corner Peak:
        ! FUNC6 = (1 + \sum(c_i*x_i))**(-d-1); summed over all the
        ! dimensions
        ! ACCEPTS ANY VECTOR OF DIMENSION D
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC6
        INTEGER(I4B), INTENT(IN)    :: D
        REAL(DP), INTENT(IN)        :: X(1:D)
        ! TEMP VARIABLES:
        REAL(DP)        :: ARG(1:D)
        INTEGER(I4B)    :: I
        REAL(DP)        :: C(1:D)
        C = 1.85D0/D
        ARG = (/ ( C(I)*X(I), I = 1, D ) /)
        FUNC6 = (1.D0 + SUM(ARG))**(- D - 1)
      END FUNCTION FUNC6

      !==================================================================
      !==================================================================
      
      FUNCTION FUNC7(X, D, K, Q, DIR_NUM)
        ! Helmholtz stifness:
        ! Ref: http://www.bettess.co.uk/pete_web/pete_prof/IGRfell03_files/reportsrf02.pdf
        ! FUNC8 SIN(K*(\SUM C_I*X_I)) summed over all the dimensions
        ! ACCEPTS ANY VECTOR OF DIMENSION D
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC7
        INTEGER(I4B), INTENT(IN)    :: D, Q, DIR_NUM
        REAL(DP), INTENT(IN)        :: K
        REAL(DP), INTENT(IN)        :: X(1:D)
        
        ! TEMP VARIABLES:
        REAL(DP)        :: XI(1:2)
        REAL(DP)        :: THETA_N, ARG
        ! AS PER PUFEM, IN 2-DIMENSIONS XI(1) = COS(2*PI*DIR_NUM/Q), 
        !                               XI(2) = SIN(2*PI*DIR_NUM/Q)
        ! WITH DIR_NUM = 0, 1, ..., Q-1; Q = # OF ENRICHMENTS = # OF EQUISPACED 
        ! ENRICHMENT DIRECTIONS. 

        THETA_N = 2.D0*PI_D*DIR_NUM/Q + PI_D/4.D0
        XI(1) = COS(THETA_N)
        XI(2) = SIN(THETA_N)
        IF (D .EQ. 2) THEN
          ARG = K*DOT_PRODUCT(XI, X)
          FUNC7 = COS(ARG)
        ELSE
          WRITE(*,*)'ERROR!'
          WRITE(*,*)'The Code is Valid for only 2-dimensions as of now.'
          STOP
        ENDIF

      END FUNCTION FUNC7

      !==================================================================
      !==================================================================

!      FUNCTION FUNC7(X, D, K)
!        ! Helmholtz stifness:
!        ! Ref: http://www.bettess.co.uk/pete_web/pete_prof/IGRfell03_files/reportsrf02.pdf
!        ! FUNC7 COS(K*(\SUM C_I*X_I)) summed over all the dimensions
!        ! ACCEPTS ANY VECTOR OF DIMENSION D
!        USE NRTYPE
!        IMPLICIT NONE
!        REAL(DP)                    :: FUNC7
!        INTEGER(I4B), INTENT(IN)    :: D
!        REAL(DP), INTENT(IN)        :: K
!        REAL(DP), INTENT(IN)        :: X(1:D)
!        
!        ! TEMP VARIABLES:
!        REAL(DP)        :: ARG(1:D)
!        INTEGER(I4B)    :: I
!        REAL(DP)        :: C(1:D)
!        ! THETA = 45 DEGREES
!        ! RESTRICTIONS ON CVECTOR C: 
!        ! 1. C(1)**2 + (2)**2 + ... + C(D)**2 = 1
!        ! 2. THE ABOVE IMPLIES -1<= C(I) <= 1 FOR ALL I = 1, ..., D.
!        ! 3. AS PER PUFEM, IN 2-DIMENSIONS C(1) = COS(2*PI*N/Q), 
!        !                                  C(2) = SIN(2*PI*N/Q)
!        ! WITH N = 0, 1, ..., Q-1; Q = # OF ENRICHMENTS = # OF EQUISPACED 
!        ! ENRICHMENT DIRECTIONS. 
!        !C = 1.D0/SQRT(2.D0)
!        C = -1.D0
!        ARG = (/ ( C(I)*X(I), I = 1, D ) /)
!        ! The following is one of the factors of a typical stiffness
!        ! matrix term arising in Helmholtz equation.
!        FUNC7 = COS(K*SUM(ARG))
!      END FUNCTION FUNC7
!
      !==================================================================
      !==================================================================

      FUNCTION FUNC8(X, D, K, Q, DIR_NUM)
        ! Helmholtz stifness:
        ! Ref: http://www.bettess.co.uk/pete_web/pete_prof/IGRfell03_files/reportsrf02.pdf
        ! FUNC8 SIN(K*(\SUM C_I*X_I)) summed over all the dimensions
        ! ACCEPTS ANY VECTOR OF DIMENSION D
        USE NRTYPE
        IMPLICIT NONE
        REAL(DP)                    :: FUNC8
        INTEGER(I4B), INTENT(IN)    :: D, Q, DIR_NUM
        REAL(DP), INTENT(IN)        :: K
        REAL(DP), INTENT(IN)        :: X(1:D)
        
        ! TEMP VARIABLES:
        REAL(DP)        :: XI(1:2)
        REAL(DP)        :: THETA_N, ARG
        ! AS PER PUFEM, IN 2-DIMENSIONS XI(1) = COS(2*PI*DIR_NUM/Q), 
        !                               XI(2) = SIN(2*PI*DIR_NUM/Q)
        ! WITH DIR_NUM = 0, 1, ..., Q-1; Q = # OF ENRICHMENTS = # OF EQUISPACED 
        ! ENRICHMENT DIRECTIONS. 

        THETA_N = 2.D0*PI_D*DIR_NUM/Q + PI_D/4.D0
        XI(1) = COS(THETA_N)
        XI(2) = SIN(THETA_N)
        IF (D .EQ. 2) THEN
          ARG = K*DOT_PRODUCT(XI, X)
          FUNC8 = SIN(ARG)
        ELSE
          WRITE(*,*)'ERROR!'
          WRITE(*,*)'The Code is Valid for only 2-dimensions as of now.'
          STOP
        ENDIF

      END FUNCTION FUNC8

      !==================================================================
      !==================================================================

      !FUNCTION FUNC7(X, D, Q, N, WV_NUM)
      !  ! Plane Wave Enrichments:
      !  ! FUNC7 = COS(K*XI.X), XI = [COS(THETA_N), SIN(THETA_N)]
      !  ! THETA_N = 2*PI*N/Q, N = 0, 1, ..., (Q - 1)
      !  ! ACCEPTS ANY VECTOR X OF DIMENSION D
      !  USE NRTYPE
      !  IMPLICIT NONE
      !  REAL(DP)                    :: FUNC7
      !  INTEGER(I4B), INTENT(IN)    :: D, Q, N
      !  REAL(DP), INTENT(IN)        :: X(1:D)
      !  REAL(DP), INTENT(IN)        :: WV_NUM       ! THE WAVE NUMBER
      !  ! TEMP VARIABLES:
      !  REAL(DP)        :: XI(1:2)
      !  INTEGER(I4B)    :: I
      !  REAL(DP)        :: THETA_N, ARG
      !  THETA_N = 2.D0*PI_D*N/Q
      !  XI(1) = COS(THETA_N)
      !  XI(2) = SIN(THETA_N)
      !  IF (D .EQ. 2) THEN
      !    ARG = WV_NUM*DOT_PRODUCT(XI, X)
      !    FUNC7 = COS(ARG)
      !  ELSE
      !    WRITE(*,*)'ERROR!'
      !    WRITE(*,*)'The Code is Valid for only 2-dimensions as of now.'
      !    STOP
      !  ENDIF
      !END FUNCTION FUNC7

      !==================================================================
      !==================================================================

      !FUNCTION FUNC8(X, D, Q, N, WV_NUM)
      !  ! Plane Wave Enrichments:
      !  ! FUNC8 = SIN(K*XI.X), XI = [COS(THETA_N), SIN(THETA_N)]
      !  ! THETA_N = 2*PI*N/Q, N = 0, 1, ..., (Q - 1)
      !  ! ACCEPTS ANY VECTOR X OF DIMENSION D
      !  USE NRTYPE
      !  IMPLICIT NONE
      !  REAL(DP)                    :: FUNC8
      !  INTEGER(I4B), INTENT(IN)    :: D, Q, N
      !  REAL(DP), INTENT(IN)        :: X(1:D)
      !  REAL(DP), INTENT(IN)        :: WV_NUM       ! THE WAVE NUMBER
      !  ! TEMP VARIABLES:
      !  REAL(DP)        :: XI(1:2)
      !  INTEGER(I4B)    :: I
      !  REAL(DP)        :: THETA_N, ARG
      !  THETA_N = 2.D0*PI_D*N/Q
      !  XI(1) = COS(THETA_N)
      !  XI(2) = SIN(THETA_N)
      !  IF (D .EQ. 2) THEN
      !    ARG = WV_NUM*DOT_PRODUCT(XI, X)
      !    FUNC7 = SIN(ARG)
      !  ELSE
      !    WRITE(*,*)'ERROR!'
      !    WRITE(*,*)'The Code is Valid for only 2-dimensions as of now.'
      !    STOP
      !  ENDIF
      !END FUNCTION FUNC8

      !==================================================================
      !==================================================================

      FUNCTION FUNC9(L, M, X0, ALPHA, X)
        ! R_{nl}(r)*Y^m_l(\theta,\phi) => Product of radial solution of 
        ! Schrodinger equation and the spherical harmonics. The radial
        ! solution is taken to be the model pseduso-potential which
        ! mimics a Gaussian. 
        ! 
        ! The Spherical harmonic (the quantum mechanics version)
        ! incorporating Condon–Shortley phase
        ! Works only for 3-dimensions
        ! Ref: specialFunctions.f90 of pufee package
        ! =====
        ! NOTE:
        ! =====
        ! real spherical harmonic function
        ! Real orthonormal set defined in terms of the complex spherical
        ! harmonics Y_l^m:
        !           Ylm,                        m = 0                                
        !  val =   (-1)**m sqrt(2) Re(Ylm),     m > 0                                
        !          (-1)**|m| sqrt(2) Im(Yl|m|), m < 0 
          
        USE NRTYPE
        IMPLICIT NONE
        
        INTERFACE
          RECURSIVE FUNCTION FACTORIAL(N)  RESULT(FACT)
            USE NRTYPE
            INTEGER(I4B)                :: FACT
            INTEGER(I4B), INTENT(IN)    :: N
          END FUNCTION FACTORIAL
        END INTERFACE
            
        INTERFACE
          FUNCTION PLGNDR(L, M, X)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)    :: L, M
            REAL(DP), INTENT(IN)        :: X
            REAL(DP)                    :: PLGNDR
          END FUNCTION PLGNDR
        END INTERFACE          

        INTEGER(I4B), INTENT(IN)    :: L, M
        REAL(DP), INTENT(IN)        :: X(1:3), X0(1:3)
        REAL(DP), INTENT(IN)        :: ALPHA
        REAL(DP)                    :: FUNC9
        
        ! TEMP VARIABLES:
        REAL(DP)        :: XX, YY, ZZ
        REAL(DP)        :: R, PHI, THETA, NORM_CONST
        REAL(DP)        :: ARG(1:3)
        INTEGER(I4B)    :: I, M_ABS
        
        XX = X(1)
        YY = X(2)
        ZZ = X(3)

        ARG = (/ ( (X(I) - X0(I))**2, I = 1, 3) /)

        R = SQRT(XX*XX + YY*YY +ZZ*ZZ)
        THETA = ACOS(ZZ/R)
        PHI = ATAN(YY/XX)
               
        IF (M .EQ. 0) THEN

          NORM_CONST = (-1)**M*SQRT((2*L + 1)/4.D0/PI)

          FUNC9 =  NORM_CONST*PLGNDR(L, M, COS(THETA))

        ELSEIF (M .LT. 0) THEN

          M_ABS = ABS(M)

          NORM_CONST = (-1)**M_ABS*SQRT(2.D0)*SQRT((2*L + 1)            &
     &                 *FACTORIAL(L - M_ABS)/4.D0/PI                    &
     &                 /FACTORIAL(L + M_ABS))

          FUNC9 = NORM_CONST*PLGNDR(L, M_ABS, COS(THETA))*SIN(M_ABS*PHI)

        ELSEIF (M .GT. 0) THEN

          NORM_CONST = (-1)**M*SQRT(2.D0)*SQRT((2*L + 1)*               &
     &                 FACTORIAL(L - M)/4.D0/PI/FACTORIAL(L + M))
             
          FUNC9 = NORM_CONST*PLGNDR(L, M, COS(THETA))*COS(M*PHI)

        ENDIF

        FUNC9 = EXP(-ALPHA*SUM(ARG))*FUNC9

        RETURN
      END FUNCTION FUNC9

      !==================================================================
      !==================================================================

      FUNCTION PLGNDR(L, M, X)
        ! Computes the associated Legendre polynomial P^M_L(X).
        ! Here M and L are integers satisfying
        ! 0 ≤ M ≤ L; while X lies in the range [−1, 1].
        ! REF: http://perso.ens-lyon.fr/christophe.winisdoerffer/INTRO_NUM/NumericalRecipesinF77.pdf


        USE NRTYPE
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN)    :: L, M
        REAL(DP), INTENT(IN)        :: X
        REAL(DP)                    :: PLGNDR

        ! TEMP VARIABLES:        
        INTEGER(I4B)                :: I, LL
        REAL(DP)                    :: FACT, PLL, PMM, PMMP1, SOMX2
        IF (M.LT.0.OR.M.gt.L.OR.abs(X).GT.1.) THEN
            WRITE(*,*) 'Bad arguments in PLGNDR!'
            WRITE(*,*) 'EXITING'
            STOP
        ENDIF
        
        ! Compute P^M_M
        PMM = 1.D0

        IF (M.gt.0) THEN
            SOMX2 = sqrt((1.D0 - X)*(1.D0 + X))
            FACT = 1.D0

            DO I = 1, M
                PMM = -PMM*FACT*SOMX2
                FACT = FACT + 2.D0
            ENDDO
        ENDIF

        IF (L.EQ.M) THEN
            PLGNDR = PMM
        ELSE
            PMMP1 = X*(2*M+1)*PMM   !Compute P^M_{M+1}
            IF (L.EQ.M+1) THEN
                PLGNDR = PMMP1
            ELSE 
                !Compute P^M_L; L > M + 1
                DO LL = M+2, L
                    PLL = (X*(2*LL-1)*PMMP1 - (LL + M - 1)*PMM)/(LL - M)
                    PMM = PMMP1
                    PMMP1 = PLL
                ENDDO

                PLGNDR = PLL
            ENDIF
        ENDIF
        RETURN
      END FUNCTION PLGNDR

      !==================================================================
      !==================================================================

      SUBROUTINE NUM_INTEGRATION(INTEGRAL_OUT, C, S, POLY_TYPE, M, D,   &
     &                           R, INT_ORDER, A, B)
        ! Numerical Integration routine for 1-D, with inputs as:
        ! i. poly_type
        ! ii. c0, c1, c2, etc, i.e., all the columns of C.
        USE NRTYPE
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE NUM_INTEGRATION_1D(INT_OUT_1D, C_VEC,              &
     &                              POLY_TYPE, INT_ORDER, A, B, M)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: INT_ORDER, M
            CHARACTER(80)                       :: POLY_TYPE
            REAL(DP), INTENT(IN)                :: C_VEC(M+1)
            REAL(DP), INTENT(IN)                :: A, B
            REAL(DP), INTENT(OUT)               :: INT_OUT_1D
          END SUBROUTINE NUM_INTEGRATION_1D
        END INTERFACE

        ! Interfacing part of the dummy variables
        INTEGER(I4B), INTENT(IN)            :: M, D, INT_ORDER, R
        REAL(DP), INTENT(IN)                :: A, B
        REAL(DP), INTENT(IN)                :: C(1:R, 1:(M+1), 1:D)
        REAL(DP), INTENT(IN)                :: S(1:R)
        CHARACTER(80)                       :: POLY_TYPE
        REAL(DP), INTENT(OUT)               :: INTEGRAL_OUT

        ! LOCAL VARIABLES
        INTEGER(I4B)                :: L, K
        REAL(DP), ALLOCATABLE       :: INT_OUT_1DV(:), C_VECTR(:)
        REAL(DP)                    :: INT_OUT_VEC(1:R)
        
        DO L = 1, R
          ALLOCATE( INT_OUT_1DV(1:D) )
          DO K = 1, D
            ALLOCATE( C_VECTR(1:M+1) )
            C_VECTR = RESHAPE(C(L,:,K), (/(M+1) /))
            CALL NUM_INTEGRATION_1D(INT_OUT_1DV(K), C_VECTR,            &
     &                              POLY_TYPE, INT_ORDER, A, B, M)
            DEALLOCATE( C_VECTR )
          ENDDO
          INT_OUT_VEC(L) = S(L)*PRODUCT(INT_OUT_1DV)
          DEALLOCATE( INT_OUT_1DV )
        ENDDO

        INTEGRAL_OUT = SUM(INT_OUT_VEC)
        RETURN

      END SUBROUTINE NUM_INTEGRATION

      !==================================================================
      !==================================================================

      SUBROUTINE NUM_INTEGRATION_1D(INT_OUT_1D, C_VEC, POLY_TYPE,       &
     &                              INT_ORDER, A, B, M)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE GAUSSPOINTS(X, N)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(OUT)           :: X(1:N)
          END SUBROUTINE GAUSSPOINTS
        END INTERFACE

        INTERFACE
          SUBROUTINE GAUSSWEIGHTS(W, N)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(OUT)           :: W(1:N)
          END SUBROUTINE GAUSSWEIGHTS
        END INTERFACE

        INTERFACE
          SUBROUTINE ORTH_POLY(PNX, POLY_TYPE, N, X, SIZE_OF_X, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
            INTEGER(I4B), INTENT(IN)        :: SIZE_OF_X       ! # OF PNTS TO EVAL       
            REAL(DP), INTENT(IN)            :: X(1:SIZE_OF_X)  ! PNTS TO EVAL
            CHARACTER(80), INTENT(IN)       :: POLY_TYPE
            REAL(DP), INTENT(OUT)           :: PNX(1:SIZE_OF_X)
            REAL(DP), OPTIONAL, INTENT(IN)  :: A, B
          END SUBROUTINE ORTH_POLY
        END INTERFACE

        ! INTERFACING VARIABLES:
        INTEGER(I4B), INTENT(IN)            :: INT_ORDER, M
        CHARACTER(80)                       :: POLY_TYPE
        REAL(DP), INTENT(IN)                :: C_VEC(M+1)
        REAL(DP), INTENT(IN)                :: A, B
        REAL(DP), INTENT(OUT)               :: INT_OUT_1D
        ! LOCAL VARIABLES:
        INTEGER(I4B)                :: I, J
        REAL(DP)                    :: PSI(1:(M+1), 1:INT_ORDER)
        REAL(DP)                    :: GP_ARRAY(INT_ORDER, 2)
        REAL(DP)                    :: TERM(1:INT_ORDER)

        CALL GAUSSPOINTS(GP_ARRAY(:, 1), INT_ORDER)
        CALL GAUSSWEIGHTS(GP_ARRAY(:, 2), INT_ORDER)
        
        SELECT CASE (POLY_TYPE)
          CASE ('ORTH_LEGENDRE')
            !$OMP PARALLEL
            !$OMP DO
            DO I = 0, M
              !$OMP DO
              DO J = 1, INT_ORDER
                CALL ORTH_POLY(PSI((I+1),J), POLY_TYPE, I,              &
     &            (/ (0.5D0*((B - A)*GP_ARRAY(J,1) + (B + A))) /), 1)
              ENDDO
              !$OMP ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP END PARALLEL
          CASE ('ORTH_CHEBYSHEV')
            !$OMP PARALLEL
            !$OMP DO
            DO I = 0, M
              !$OMP DO
              DO J = 1, INT_ORDER
                CALL ORTH_POLY(PSI((I+1),J), POLY_TYPE, I,              &
     &            (/ (0.5D0*((B - A)*GP_ARRAY(J,1) + (B + A))) /), 1,   &
     &            A, B)
              ENDDO
              !$OMP ENDDO
            ENDDO
            !$OMP ENDDO
            !$OMP END PARALLEL
        END SELECT


        SELECT CASE (POLY_TYPE)
          CASE ('ORTH_LEGENDRE')
            !$OMP PARALLEL
            !$OMP DO
            DO I = 1, INT_ORDER
              TERM(I) = GP_ARRAY(I,2)*(DOT_PRODUCT(C_VEC, PSI(:,I)))
            ENDDO
            !OMP END DO
            !$OMP END PARALLEL
          CASE ('ORTH_CHEBYSHEV')
            !$OMP PARALLEL
            !$OMP DO
            DO I = 1, INT_ORDER
              TERM(I) = GP_ARRAY(I,2)*(DOT_PRODUCT(C_VEC, PSI(:,I)))
            ENDDO
            !OMP END DO
            !$OMP END PARALLEL
          CASE DEFAULT
            PRINT *, 'HERMITE AND OTHER WEIGHTS NOT CODED YET!'
            STOP
        END SELECT
        INT_OUT_1D = 0.5D0*(B - A)*SUM(TERM)
        RETURN

      END SUBROUTINE NUM_INTEGRATION_1D

      !==================================================================
      !==================================================================

      SUBROUTINE NUM_INTEGRATION_BRN(INTEGRAL_OUT, C, S, M, D, R,       &
     &                           INT_ORDER, A, B)
        ! Numerical Integration routine for 1-D, with inputs as:
        ! i. poly_type
        ! ii. c0, c1, c2, etc, i.e., all the columns of C.
        USE NRTYPE
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE NUM_INTEGRATION_1D_BRN(INT_OUT_1D, C_VEC,          &
     &                                  INT_ORDER, A, B, M)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)            :: INT_ORDER, M
            REAL(DP), INTENT(IN)                :: C_VEC(M+1)
            REAL(DP), INTENT(IN)                :: A, B
            REAL(DP), INTENT(OUT)               :: INT_OUT_1D
          END SUBROUTINE NUM_INTEGRATION_1D_BRN
        END INTERFACE

        ! Interfacing part of the dummy variables
        INTEGER(I4B), INTENT(IN)            :: M, D, INT_ORDER, R
        REAL(DP), INTENT(IN)                :: A, B
        REAL(DP), INTENT(IN)                :: C(1:R, 1:(M+1), 1:D)
        REAL(DP), INTENT(IN)                :: S(1:R)
        REAL(DP), INTENT(OUT)               :: INTEGRAL_OUT

        ! LOCAL VARIABLES
        INTEGER(I4B)                :: L, K
        REAL(DP), ALLOCATABLE       :: INT_OUT_1DV(:), C_VECTR(:)
        REAL(DP)                    :: INT_OUT_VEC(1:R)
        
        DO L = 1, R
          ALLOCATE( INT_OUT_1DV(1:D) )
          DO K = 1, D
            ALLOCATE( C_VECTR(1:M+1) )
            C_VECTR = RESHAPE(C(L,:,K), (/(M+1) /))

            CALL NUM_INTEGRATION_1D_BRN(INT_OUT_1DV(K), C_VECTR,        &
     &                                  INT_ORDER, A, B, M)

            DEALLOCATE( C_VECTR )
          ENDDO

          INT_OUT_VEC(L) = S(L)*PRODUCT(INT_OUT_1DV)
          DEALLOCATE( INT_OUT_1DV )
        ENDDO

        INTEGRAL_OUT = SUM(INT_OUT_VEC)

      END SUBROUTINE NUM_INTEGRATION_BRN

      !==================================================================
      !==================================================================

      SUBROUTINE NUM_INTEGRATION_1D_BRN(INT_OUT_1D, C_VEC, INT_ORDER,   &
     &                                  A, B, M)
        USE NRTYPE
        USE OMP_LIB
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE GAUSSPOINTS(X, N)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(OUT)           :: X(1:N)
          END SUBROUTINE GAUSSPOINTS
        END INTERFACE

        INTERFACE
          SUBROUTINE GAUSSWEIGHTS(W, N)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(OUT)           :: W(1:N)
          END SUBROUTINE GAUSSWEIGHTS
        END INTERFACE

        INTERFACE
          SUBROUTINE BERNSTEIN(BRN_BASIS_POLY, N, A, B, X)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
            REAL(DP), INTENT(IN)            :: X               ! POINT TO EVAL
            REAL(DP), INTENT(IN)            :: A, B            ! DOMAIN OF THE PROBLEM
                                                               ! X \IN [A, B]
            REAL(DP), INTENT(OUT)           :: BRN_BASIS_POLY(0:N)
          END SUBROUTINE BERNSTEIN
        END INTERFACE

        ! INTERFACING VARIABLES:
        INTEGER(I4B), INTENT(IN)            :: INT_ORDER, M
        REAL(DP), INTENT(IN)                :: C_VEC(M+1)
        REAL(DP), INTENT(IN)                :: A, B
        REAL(DP), INTENT(OUT)               :: INT_OUT_1D
        ! LOCAL VARIABLES:
        INTEGER(I4B)                :: I, J
        REAL(DP)                    :: PSI(1:(M+1), 1:INT_ORDER)
        REAL(DP)                    :: GP_ARRAY(INT_ORDER, 2)
        REAL(DP)                    :: TERM(1:INT_ORDER)

        CALL GAUSSPOINTS(GP_ARRAY(:, 1), INT_ORDER)
        CALL GAUSSWEIGHTS(GP_ARRAY(:, 2), INT_ORDER)
        
        !$OMP PARALLEL
        !!$OMP DO
        !DO I = 0, M
          !$OMP DO
          DO J = 1, INT_ORDER
            CALL BERNSTEIN(PSI(:,J), M, A, B,                           &
     &           (0.5D0*((B - A)*GP_ARRAY(J,1) + (B + A))))


          ENDDO
          !$OMP ENDDO
        !ENDDO
        !!$OMP ENDDO
        !$OMP END PARALLEL

        !$OMP PARALLEL
        !$OMP DO
        DO I = 1, INT_ORDER
          TERM(I) = GP_ARRAY(I,2)*(DOT_PRODUCT(C_VEC, PSI(:,I)))
        ENDDO
        !OMP END DO
        !$OMP END PARALLEL

        INT_OUT_1D = 0.5D0*(B - A)*SUM(TERM)

      END SUBROUTINE NUM_INTEGRATION_1D_BRN

      !==================================================================
      !==================================================================

      !SUBROUTINE ALS_OBJ_FUNC(OBJ_FVAL, C_VEC, C, S, R, M, N, U,        &
     !&                        DIM_OPT, D, POLY_TYPE, COORD)
        ! C_VEC[r, (M + 1)] & C[r, (M + 1), d]
        ! N = SAMPLE_SIZE
        ! the argument input, "C_VEC" must be a vector
        ! I.e., the array: c_arr[r, (M + 1)]
        ! c_arr = [c_0^1 c_1^1 ... c_M^1;
        !          c_0^2 c_1^2 ... c_M^2;
        !           .          ...  .
        !           .          ...  .
        !           .          ...  .
        !          c_0^r c_1^r ... c_M^r]
        !  should be rearranged to --
        !  c_vec = [c_0^1 c_0^2 ... c_0^r ... c_M^1 ... c_M^r]
        !USE NRTYPE
        !USE OMP_LIB
        !IMPLICIT NONE
        ! INTERFACE
        !INTERFACE
        !  SUBROUTINE ORTH_POLY(PNX, POLY_TYPE, N, X, SIZE_OF_X, A, B)
        !    USE NRTYPE
        !    USE OMP_LIB
        !    INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
        !    INTEGER(I4B), INTENT(IN)        :: SIZE_OF_X       ! # OF PNTS TO EVAL       
        !    REAL(DP), INTENT(IN)            :: X(1:SIZE_OF_X)  ! PNTS TO EVAL
        !    CHARACTER(80), INTENT(IN)       :: POLY_TYPE
        !    REAL(DP), INTENT(OUT)           :: PNX(1:SIZE_OF_X)
        !    REAL(DP), OPTIONAL, INTENT(IN)  :: A, B
        !  END SUBROUTINE ORTH_POLY
        !END INTERFACE

        !INTERFACE
        !  SUBROUTINE TNSR_UHAT_D_DIM(PROD_U_HAT, NODES, M, D, C_RED,    &
     !&                               POLY_TYPE)
        !    USE NRTYPE
        !    USE OMP_LIB
        !    INTEGER(I4B), INTENT(IN)        :: M, D
        !    REAL(DP), INTENT(IN)            :: C_RED(1:(M+1),1:(D-1))
        !    REAL(DP), INTENT(IN)            :: NODES(1:(D-1))
        !    CHARACTER(80), INTENT(IN)       :: POLY_TYPE
        !    REAL(DP), INTENT(OUT)           :: PROD_U_HAT
        !  END SUBROUTINE TNSR_UHAT_D_DIM
        !END INTERFACE

        !INTERFACE
        !  FUNCTION SEMI_NORM(U, N)
        !    USE NRTYPE
        !    REAL(DP)                        :: SEMI_NORM
        !    INTEGER(I4B), INTENT(IN)        :: N
        !    REAL(DP), INTENT(IN)            :: U(1:N)
        !  END FUNCTION SEMI_NORM
        !END INTERFACE

      !  ! Interfacing part of the dummy variables:
      !  INTEGER(I4B), INTENT(IN)            :: R, M, N, DIM_OPT, D
      !  REAL(DP), INTENT(IN)                :: C_VEC(R*(M+1)), U(1:N)
      !  REAL(DP), INTENT(IN)                :: C(R,(M+1),D), S(1:R)
      !  REAL(DP), INTENT(IN)                :: COORD(N,D)
      !  CHARACTER(80), INTENT(IN)           :: POLY_TYPE
      !  REAL(DP), INTENT(OUT)               :: OBJ_FVAL

      !  ! Local Variables:
      !  REAL(DP)                            :: SUMM, PROD_U_HAT
      !  REAL(DP), ALLOCATABLE               :: C_ARR(:,:)
      !  REAL(DP)                            :: U_HAT(1:N)
      !  REAL(DP), ALLOCATABLE               :: NODES(:)
      !  REAL(DP), ALLOCATABLE               :: C_RED(:,:), PNX(:)
      !  INTEGER(I4B), ALLOCATABLE           :: DIM_ARR(:)
      !  INTEGER(I4B)                        :: J, J1, J2, L, K, ALPHA
      !  INTEGER(I4B)                        :: SIZE_OF_X = 1
      !  ALLOCATE( C_ARR(R, (M+1)) )
      !  C_ARR = 0.D0
      !  C_ARR = RESHAPE(C_VEC, (/ R, (M + 1) /))
      !  U_HAT = 0.D0

      !  ! DIM_ARR(D - 1) -- Contains all the dimension except the one in
      !  ! opt-loop, i.e., DIM_OPT
      !  ALLOCATE( DIM_ARR(D-1) )                   ! **** need to deallocate
      !  DIM_ARR = 0
      !  IF (DIM_OPT .EQ. 1) THEN
      !    DIM_ARR = (/ (J, J = 2, D) /)
      !  ELSEIF (DIM_OPT .EQ. D) THEN
      !    DIM_ARR = (/ (J, J = 1, (D-1)) /)
      !  ELSE
      !    DIM_ARR(1:(DIM_OPT-1)) = (/ (J1, J1 = 1, (DIM_OPT-1)) /)
      !    DIM_ARR(DIM_OPT:(D-1)) = (/ (J2, J2 = (DIM_OPT+1), D) /)
      !  ENDIF

      !  DO J = 1, N
      !    ALLOCATE( NODES(D-1) )                ! **** need to deallocate
      !    NODES = 0.D0
      !    IF (DIM_OPT .EQ. 1) THEN
      !      NODES = COORD(J,2:D)
      !    ELSEIF (DIM_OPT .EQ. D) THEN
      !      NODES = COORD(J,1:D-1)
      !    ELSE
      !      NODES(1:DIM_OPT-1) = COORD(J,1:(DIM_OPT-1))
      !      NODES(DIM_OPT:D-1) = COORD(J,(DIM_OPT+1):D)
      !    ENDIF

      !    DO L = 1, R
      !      ! C_RED((M + 1), (d - 1)) -- Contains all the coefficient at a
      !      ! given rank except the coeffiecients which go to the opt-loop
      !      ! Computed/Selected for each rank L
      !      ALLOCATE( C_RED(1:(M+1), 1:(D - 1)) )   ! **** need to deallocate
      !      C_RED = 0.D0
      !      IF (DIM_OPT .EQ. 1) THEN
      !        C_RED = RESHAPE(C(L,:,2:D), (/ (M+1), (D-1) /))
      !      ELSEIF (DIM_OPT .EQ. D) THEN
      !        C_RED = RESHAPE(C(L,:,1:(D-1)), (/ (M+1), (D-1) /))
      !      ELSE
      !        C_RED(:,1:(DIM_OPT-1)) = C(L,:,1:(DIM_OPT-1))
      !        C_RED(:,DIM_OPT:D-1) = C(L,:,(DIM_OPT+1):D)
      !      ENDIF
            
      !      SUMM = 0.D0
      !      DO K = 1, D
      !        IF (K .EQ. DIM_OPT) THEN
      !          !$OMP DO 
      !          DO ALPHA = 0, M
      !            ! Product over all the dimensions except "dim_opt"
      !            ! and then summed over l for a given node (j)
      !            ALLOCATE( PNX(1:SIZE_OF_X) )
      !            PNX = 0.D0
      !           CALL ORTH_POLY(PNX, POLY_TYPE, ALPHA,                 &
     !&                       (/ COORD(J,DIM_OPT) /), SIZE_OF_X)

      !            SUMM = SUMM + C_ARR(L,ALPHA+1)*PNX(1)
      !            DEALLOCATE( PNX )
      !          ENDDO
      !        ENDIF
      !      ENDDO

      !     CALL TNSR_UHAT_D_DIM(PROD_U_HAT, NODES, M, D, C_RED,        &
     !&                           POLY_TYPE)
      !      U_HAT(J) = U_HAT(J) + PROD_U_HAT*SUMM*S(L)

      !      DEALLOCATE( C_RED )

      !    ENDDO
          
      !    DEALLOCATE( NODES )

      !  ENDDO

      !  OBJ_FVAL = (SEMI_NORM((U - U_HAT), N))**2

      !END SUBROUTINE ALS_OBJ_FUNC

      !==================================================================
      !==================================================================

!      Auxiliary routine: printing a matrix.

      SUBROUTINE PRINT_MATRIX(M, N, A)
        USE NRTYPE
        IMPLICIT NONE
        INTEGER(I4B), INTENT(IN)          :: M, N
        REAL(DP), INTENT(IN)              :: A(M,N)

        ! LOCAL VARIABLES:
        INTEGER          I, J
        DO I = 1, M
           WRITE(*,*) ( A( I, J ), J = 1, N )
        END DO
!
        RETURN
      END SUBROUTINE PRINT_MATRIX

      !==================================================================
      !==================================================================

      ! Uses the ALS scheme assuming the polynomial basis of type
      ! POLY_TYPE (Legendre/Tchebyshev). The algorithm assumes the
      ! polynomial in each variable (x, y, z, ...) to be of order M 
      ! (same for all the dimensions). The ALS scheme starts with rank 
      ! R = 1; and increments the rank by 1 if convergence criterion is
      ! not met or MAXITER (= 20) is reached. This goes on till 
      ! R = R_MAX (in case convergence is not achieved) or 
      ! R_OUT (<= R_MAX) ==> R_OUT in that case becomes the output. The
      ! corresponding coeeficients of expansion (output of the LAPACK
      ! solver) are reshaped and stored in  C(1:R_OUT, 1:(M + 1), 1:D)


      SUBROUTINE ALS_ORTH(ERR_DECAY, R_OUT, C, S, D, M, R_MAX, C_MIN,   &
     &               C_MAX, POLY_TYPE, U, SAMPLE_SIZE, COORD, EPS,      &
     &               RUNDATAFILE, SOLVER, PRBLM_TYPE, FLAG_CONV,        &
     &               MAX_ITER, T_REGU_STAT, LAMBDA, LOW_LIM, UP_LIM)
        USE NRTYPE
        USE OMP_LIB
        USE, INTRINSIC :: ISO_Fortran_env
        IMPLICIT NONE
        ! INTERFACE
        INTERFACE
          SUBROUTINE PRINT_MATRIX(M, N, A)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)          :: M, N
            REAL(DP), INTENT(IN)              :: A(M,N)
          END SUBROUTINE PRINT_MATRIX
        END INTERFACE

        INTERFACE
          SUBROUTINE TNSR_UHAT_D_DIM(PROD_U_HAT, NODES, M, D, C_RED,    &
     &                               POLY_TYPE, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: M, D
            REAL(DP), INTENT(IN)            :: C_RED(1:(M+1),1:(D-1))
            REAL(DP), INTENT(IN)            :: NODES(1:(D-1))
            CHARACTER(80), INTENT(IN)       :: POLY_TYPE
            REAL(DP), INTENT(IN)            :: A, B
            REAL(DP), INTENT(OUT)           :: PROD_U_HAT
          END SUBROUTINE TNSR_UHAT_D_DIM
        END INTERFACE
 
        INTERFACE
          SUBROUTINE ORTH_POLY(PNX, POLY_TYPE, N, X, SIZE_OF_X, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N               ! ORDER OF POLYNOMIAL
            INTEGER(I4B), INTENT(IN)        :: SIZE_OF_X       ! # OF PNTS TO EVAL       
            REAL(DP), INTENT(IN)            :: X(1:SIZE_OF_X)  ! PNTS TO EVAL
            CHARACTER(80), INTENT(IN)       :: POLY_TYPE
            REAL(DP), INTENT(OUT)           :: PNX(1:SIZE_OF_X)
            REAL(DP), OPTIONAL, INTENT(IN)  :: A, B
          END SUBROUTINE ORTH_POLY
        END INTERFACE

        INTERFACE
          FUNCTION SEMI_NORM(U, N)
            USE NRTYPE
            REAL(DP)                        :: SEMI_NORM
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(IN)            :: U(1:N)
          END FUNCTION SEMI_NORM
        END INTERFACE

        INTERFACE
          SUBROUTINE SEP_RESP(F_SEP, C, S, Y, POLY_TYPE, R, M, D, A, B)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)         :: R, M, D
            CHARACTER(80), INTENT(IN)        :: POLY_TYPE
            REAL(DP), INTENT(IN)             :: C(1:R, 1:(M+1), 1:D)
            REAL(DP), INTENT(IN)             :: S(1:R), Y(1:D)
            REAL(DP), INTENT(IN)             :: A, B
            REAL(DP), INTENT(OUT)            :: F_SEP
          END SUBROUTINE SEP_RESP
        END INTERFACE

        INTERFACE
          SUBROUTINE DIAG(DIAG_MAT, N)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(OUT)           :: DIAG_MAT(N,N)
          END SUBROUTINE DIAG
        END INTERFACE

        INTERFACE
          SUBROUTINE KRON(A_KRON_B, A, M, N, B, P, Q)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: M, N, P, Q
            REAL(DP), INTENT(IN)            :: A(M,N), B(P,Q)
            REAL(DP), INTENT(OUT)           :: A_KRON_B(M*P,N*Q)
          END SUBROUTINE KRON
        END INTERFACE

        ! Interfacing part of the dummy variables
        INTEGER(I4B), INTENT(IN)               :: D, M, R_MAX
        INTEGER(I4B), INTENT(IN)               :: SAMPLE_SIZE, MAX_ITER
        REAL(DP), ALLOCATABLE, INTENT(OUT)     :: C(:,:,:)
        INTEGER(I4B), INTENT(OUT)              :: R_OUT
        LOGICAL, INTENT(OUT)                   :: FLAG_CONV
        REAL(DP), INTENT(IN)                   :: C_MIN, C_MAX
        REAL(DP), ALLOCATABLE, INTENT(OUT)     :: S(:)
        REAL(DP), INTENT(OUT)                  :: ERR_DECAY
        CHARACTER(80), INTENT(IN)              :: POLY_TYPE, PRBLM_TYPE
        REAL(DP), INTENT(IN)                   :: U(SAMPLE_SIZE)
        REAL(DP), INTENT(IN)                   :: COORD(SAMPLE_SIZE,D)
        REAL(DP), INTENT(IN)                   :: EPS, LAMBDA
        CHARACTER(80), INTENT(IN)              :: RUNDATAFILE, SOLVER
        CHARACTER(80), INTENT(IN)              :: T_REGU_STAT
        REAL(DP), OPTIONAL, INTENT(IN)         :: LOW_LIM
        REAL(DP), OPTIONAL, INTENT(IN)         :: UP_LIM
        
        ! LOCAL VARIABLES:
        INTEGER(I4B)            :: R
        INTEGER(I4B)            :: K, L, J, J1, J2, ALPHA, JJ
        INTEGER(I4B)            :: DIM_OPT
        INTEGER(I4B)            :: SIZE_OF_X = 1
        REAL(DP), ALLOCATABLE   :: AMAT(:,:)
        REAL(DP), ALLOCATABLE   :: C_SOL_ARR(:,:), C_RED(:,:)
        REAL(DP), ALLOCATABLE   :: U_H(:,:), DIM_ARR(:),  NODES(:)
        REAL(DP), ALLOCATABLE   :: NORM_UJL(:), F_SEP_APPRX(:)
        REAL(DP)                :: PROD_U_HAT
        REAL(DP), ALLOCATABLE   :: PNX(:), RHS(:), LHS(:,:)
        REAL(DP)                :: SEMI_NRM
        REAL(DP)                :: ERR_OLD, ERR_NEW, ERR_TOL, TOL
        INTEGER(I4B)            :: ITER
        REAL(DP), ALLOCATABLE   :: F_SEP_ITER(:)
        REAL(DP), ALLOCATABLE   :: T_REG_MAT(:,:), S_MAT(:,:)

        ! Parameters used by LAPACK equation solver:                                                
        INTEGER(I4B)                            :: LDA, LDB, INFO
        INTEGER(I4B)                            :: MMAX, NMAX, RANK
        INTEGER(I4B), ALLOCATABLE               :: PIVOT(:)
        INTEGER(I4B)                            :: NRHS, NN, MM
        INTEGER(I4B)                            :: LWMAX, LWORK
        REAL(DP), ALLOCATABLE                   :: WORK(:), SS(:)
        REAL(DP), ALLOCATABLE                   :: CC(:), XX(:,:), RR(:)
        REAL(DP), ALLOCATABLE                   :: BERR(:), FERR(:)
        REAL(DP), ALLOCATABLE                   :: AF(:,:)
        REAL(DP), ALLOCATABLE                   :: IDENTITY(:,:)
        INTEGER(I4B), ALLOCATABLE               :: IWORK(:)
        REAL(DP)                                :: RCOND
        INTEGER(I4B)                            :: LD_LHS, LD_RHS
        CHARACTER(1)                            :: EQUED
        INTEGER(I4B)                            :: LDAF, LDX

        FLAG_CONV = .TRUE.
        !MAX_ITER = 10

        ! Difference among EPS(user specified), TOL, and ERR_TOL:
        ! EPS     = user specified tolerance, basically similar to a
        !           discrtete L2-norm (L2-disc) of the difference between the
        !           analytical solution and the separated reps. 
        !
        ! TOL     = EPS*1.D-06 (some fraction of EPS); TOL is < EPS as
        !           it is compared against a "DROP in L2-disc error" 
        !           (the ERR_TOL) rather than just error value (L2-disc)
        !
        ! ERR_TOL = ABS(ERR_NEW - ERR_OLD); where ERR_NEW = L2-disc at
        !           iteration # (I+1) and 
        !           ERR_OLD = L2-disc at iteration # I.
        !
        !           Computes the absolute value of change in error (decrease or increase)
        !           at ITER. We need to have a certain "DROP in error" at
        !           every iteration. If the "error drop" is not sufficiently
        !           small (governed by TOL above) then we move on to the next 
        !           ITER <= MAX_ITER for the same given rank.

        TOL = EPS*1.D-06

        ! The big loop for R:
        R = 1       ! Start with Rank 1 apprximation and then move up if
                    ! the solution is non-convergent

        DO WHILE (R .LE. R_MAX)
          !
          IF (FLAG_CONV .EQV. .FALSE.) THEN
            DEALLOCATE( F_SEP_APPRX, C, S )
          ENDIF

          ALLOCATE( F_SEP_APPRX(SAMPLE_SIZE) )
          
          ! Re-allocating C and S
          ALLOCATE( C(1:R, 1:(M + 1), 1:D) )
          ALLOCATE( S(1:R) )

          ! The following set of data remains the same for all dimensions in one
          ! iteration, i.e., for a given rank. When the rank is incremented by 1,
          ! these values again change/re-initialized, irrespective of the values from
          ! last rank-iteration. 
          ! So, we actually relinquish all the effort put in
          ! to calculate the "optimum solution" for the lower rank. 
          ! That's UNDESIRABLE!

          ! EX: j = n + FLOOR((m+1-n)*u)  ! We want to choose one from m-n+1 integers
          ! The intrinsic random_number(u) returns a real number u (or an array of such) 
          ! from the uniform distribution over the interval [0,1). 
          ! [That is, it includes 0 but not 1.]
          ! To have a discrete uniform distribution on the integers {n, n+1, ..., m-1, m}
          ! carve the continuous distribution up into m+1-n equal sized chunks, mapping
          ! each chunk to an integer. One way could be:
          ! call random_number(u)
          ! As you can see, for the initial question for {0, 1, 2, 3, 4, 5} this reduces to
          ! call random_number(u)
          ! j = FLOOR(6*u)            ! n=0 and m=5


          ! Randomly initialize R and S:
          CALL RANDOM_NUMBER(C)
          C = C_MIN + (C_MAX + 1.D0 - C_MIN)*C   ! This is in [C_MIN, C_MAX]

          CALL RANDOM_NUMBER(S)

          ! The following two needs to be re-initialized for every rank
          ! R as we are doing everything from scratch as the rank gets
          ! incremented

          ERR_TOL = 10.D0*TOL
          ITER = 0
          
          DO WHILE (ERR_TOL .GT. TOL) 
            IF (ITER .GE. MAX_ITER) THEN
              WRITE(*,*)' '
              WRITE(*,*)'EXCEEDED MAXIMUM ALLOWED', MAX_ITER,           &
     &                  'ITERATIONS FOR RANK: ', R
              WRITE(*,*)'MOVING ON TO THE NEXT RANK'
              WRITE(*,*)' ' 
              GO TO 998         ! SHOULD JUMP TO NEXT RANK
            ELSE
              ITER = ITER + 1
              WRITE(*,*)' '
              WRITE(*,*)'STARTING ITERATION #', ITER, 'OF RANK: ', R
              WRITE(*,*)' '

              ! Alternating direction loop:
              DO K = 1, D
                DIM_OPT = K
                ! AMAT is initialized for every dimension. The calculated value of
                ! C ans S are included in the already initialized C and
                ! S for every new AMAT

                ! This is a rectangular system; overdetermined system
                ! solved by method of least squares.
                ALLOCATE( AMAT(SAMPLE_SIZE, R*(M + 1)) )   ! **** need to deallocate

                ! DIM_ARR(D - 1) -- Contains all the dimension except the one in
                ! opt-loop, i.e., DIM_OPT
                ALLOCATE( DIM_ARR(D-1) )                   ! **** need to deallocate
                IF (DIM_OPT .EQ. 1) THEN
                  DIM_ARR = (/ (J, J = 2, D) /)
                ELSEIF (DIM_OPT .EQ. D) THEN
                  DIM_ARR = (/ (J, J = 1, (D-1)) /)
                ELSE
                  DIM_ARR(1:(DIM_OPT-1)) = (/ (J1, J1 = 1, (DIM_OPT-1)) &
     &              /)
                  DIM_ARR(DIM_OPT:(D-1)) = (/ (J2, J2 = (DIM_OPT+1), D) &
     &             /)
                ENDIF

                !$OMP DO
                DO L = 1, R
                  ! C_RED((M + 1), (d - 1)) -- Contains all the coefficient at a
                  ! given rank except the coeffiecients which go to the opt-loop
                  ! Computed/Selected for each rank L
                  ALLOCATE( C_RED(1:(M+1), 1:(D - 1)) )   ! **** need to deallocate
                  IF (DIM_OPT .EQ. 1) THEN
                    C_RED = RESHAPE(C(L,:,2:D), (/ (M+1), (D-1) /))
                  ELSEIF (DIM_OPT .EQ. D) THEN
                    C_RED = RESHAPE(C(L,:,1:(D-1)), (/ (M+1), (D-1) /))
                  ELSE
                    C_RED(:,1:(DIM_OPT-1)) = C(L,:,1:(DIM_OPT-1))
                    C_RED(:,DIM_OPT:D-1) = C(L,:,(DIM_OPT+1):D)
                  ENDIF

                  !$OMP DO
                  DO J = 1, SAMPLE_SIZE
                    ! NODES(d -1) -- At a given sample point, it contains the all the
                    ! coordinate component of the sample point except the one, whose direction
                    ! goes to the opt-loop.
                    ALLOCATE( NODES(D-1) )                ! **** need to deallocate
                    IF (DIM_OPT .EQ. 1) THEN
                      NODES = COORD(J,2:D)
                    ELSEIF (DIM_OPT .EQ. D) THEN
                      NODES = COORD(J,1:D-1)
                    ELSE
                      NODES(1:DIM_OPT-1) = COORD(J,1:(DIM_OPT-1))
                      NODES(DIM_OPT:D-1) = COORD(J,(DIM_OPT+1):D)
                    ENDIF
                    
                    CALL TNSR_UHAT_D_DIM(PROD_U_HAT, NODES, M, D,       &
     &                                   C_RED, POLY_TYPE, LOW_LIM,     &
     &                                   UP_LIM)
                    !$OMP DO    
                    DO ALPHA = 0, M
                      ! Product over all the dimensions except "dim_opt"
                      ! and then summed over l for a given node (j)
                    
                      ALLOCATE( PNX(1:SIZE_OF_X) )

                      SELECT CASE (POLY_TYPE)
                        CASE ('ORTH_LEGENDRE')
                          CALL ORTH_POLY(PNX, POLY_TYPE, ALPHA,         &
     &                             (/ COORD(J, DIM_OPT) /), SIZE_OF_X)
                        CASE ('ORTH_CHEBYSHEV')
                          CALL ORTH_POLY(PNX, POLY_TYPE, ALPHA,         &
     &                             (/ COORD(J, DIM_OPT) /), SIZE_OF_X,  &
     &                              LOW_LIM, UP_LIM)
                      END SELECT
                      
                      AMAT(J,((L-1)*(M+1)+(1+ALPHA))) =                 &
     &                                        S(L)*PNX(1)*PROD_U_HAT
                      DEALLOCATE( PNX )
 
                    ENDDO       ! ALPHA LOOP ENDS
                    !$OMP ENDDO

                    DEALLOCATE( NODES )

                  ENDDO         ! J LOOP OVER ALL NODES ENDS
                  !$OMP ENDDO

                  DEALLOCATE( C_RED )

                ENDDO       ! L LOOP ENDS

                SELECT CASE (PRBLM_TYPE)

                  CASE ('NRML_EQTN')
                    !  ##############################################################  !
                    !  ####################### BAD TECHNIQUES #######################  !
                    !  ##############################################################  !

                    IF (SOLVER .EQ. 'DGELS') THEN
                      PRINT *, 'ERROR! ERROR!'
                      WRITE(*,*)'THE NORMAL EQUATIONS CAN NOT BE        &
     &                           SOLVED WITH LEAST SQUARE SOLVER DGELS.'
                      PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                      PRINT *, 'THE FOLLOWING: '
                      PRINT *, '1. DSYSV'
                      PRINT *, '2. DPOSV'
                      PRINT *, '3. DGESV'
                      PRINT *, '4. DGESVX '
                      PRINT *, 'EXITING'
                      STOP
                    ELSEIF (SOLVER .EQ. 'DGELSS') THEN
                      PRINT *, 'ERROR! ERROR!'
                      WRITE(*,*) 'THE NORMAL EQUATIONS CAN NOT BE SOLVED&
     &                 WITH LEAST SQUARE SOLVERS DGELSS'
                      PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                      PRINT *, 'THE FOLLOWING: '
                      PRINT *, '1. DSYSV'
                      PRINT *, '2. DPOSV'
                      PRINT *, '3. DGESV'
                      PRINT *, '4. DGESVX '
                      PRINT *, 'EXITING'
                      STOP
                    ENDIF

                    ALLOCATE( LHS(R*(M+1), R*(M+1)) )
                    ALLOCATE( RHS(R*(M+1)) )
                    RHS = MATMUL(TRANSPOSE(AMAT),U)
                    
                    IF (T_REGU_STAT .EQ. 'OFF') THEN
                      ! BADLY CONDITIONED (A**T*A) -- DON'T USE THIS:
                      ! *********** DON'T USE THIS ************
                      LHS = MATMUL(TRANSPOSE(AMAT),AMAT)      ! *** VERY HIGH CONDITION NUMBER
                    ELSEIF (T_REGU_STAT .EQ. 'ON') THEN
                      ! TIKHONOV REGULARIZATION WITH UNIT REGULARIZATION
                      ! MATRIX
                      ALLOCATE( IDENTITY((M+1), (M+1)) )
                      CALL DIAG(IDENTITY, (M+1))
                      ALLOCATE( S_MAT(R,1) )
                      S_MAT = RESHAPE(S, (/ R, 1 /))
                      ALLOCATE( T_REG_MAT(R*(M+1), R*(M+1)) )
                      CALL KRON(T_REG_MAT, S_MAT, R, 1, IDENTITY,       &
     &                      (M+1), (M+1))
                      LHS = MATMUL(TRANSPOSE(AMAT),AMAT) +              &
     &                      LAMBDA**2*T_REG_MAT
                      DEALLOCATE( IDENTITY, T_REG_MAT, S_MAT )
                    ENDIF


                    ! ## -------------------------------------------------------   ##
                    ! ## Available SIMPLE and DIVIDE AND CONQUER DRIVER routines:  ##
                    ! ## -------------------------------------------------------   ##
                    SELECT CASE (SOLVER)
                      CASE ('DSYSV')
                        ! ******************************************************************
                        !
                        ! **************************  SOLVER # 1   *************************
                        !
                        ! ******************************************************************
                        ! ========
                        ! DSYSV.F:
                        ! ========
                        ! The program computes the solution to the system of linear equations
                        ! with a real symmetric matrix A and multiple right-hand sides B,
                        ! where A is the coefficient matrix:
                        ! REFERENCE: 
                        ! https://software.intel.com/sites/products/documentation/doclib
                        ! /mkl_sa/11/mkl_lapack_examples/dsysv_ex.f.htm
                        NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                        NRHS = 1
                        LD_LHS = NN
                        LD_RHS = NN
                        LWMAX = 100*NN
                        ALLOCATE( PIVOT(NN) )
                        ALLOCATE( WORK(LWMAX) )
                        ! Query the optimal workspace.
                        LWORK = -1
                        CALL DSYSV('Lower', NN, NRHS, LHS, LD_LHS,      &
      &                           PIVOT, RHS, LD_RHS, WORK, LWORK, INFO)
                        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

                        ! Solve the SQUARE SYSTEM OF equations A*X = B.
                        CALL DSYSV('Lower', NN, NRHS, LHS, LD_LHS,      &
      &                           PIVOT, RHS, LD_RHS, WORK, LWORK, INFO )

                        ! Check for the exact singularity.
                        IF( INFO.GT.0 ) THEN
                          WRITE(*,*)'The element of the diagonal  '
                          WRITE(*,*)'factor D(',INFO,',',INFO,') is  '
                          WRITE(*,*)'zero, so that D is singular; the'
                          WRITE(*,*)'solution could not be computed.'
                          STOP
                        ENDIF
                        DEALLOCATE( PIVOT )

                      CASE ('DPOSV')
                        ! ******************************************************************
                        !
                        ! **************************  SOLVER # 2   *************************
                        !
                        ! ******************************************************************
                        !
                        ! ========
                        ! DPOSV.F
                        ! ========
                        ! The program computes the solution to the system of linear
                        ! equations with a symmetric positive-definite matrix A and multiple
                        ! right-hand sides B, where A is the coefficient matrix:
                        ! REFERENCE:
                        ! https://software.intel.com/sites/products/documentation/
                        ! doclib/mkl_sa/11/mkl_lapack_examples/dposv_ex.f.htm
                        NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                        NRHS = 1
                        LD_LHS = NN
                        LD_RHS = NN
                        ALLOCATE( PIVOT(NN) )
                        ALLOCATE( WORK(1) )     ! THIS IS A HACK
                        CALL DPOSV('Upper', NN, NRHS, LHS, LD_LHS, RHS, &
     &                             LD_RHS, INFO)
                        IF( INFO.GT.0 ) THEN
                          WRITE(*,*)'The leading minor of order ',INFO
                          WRITE(*,*)' is not positive definite; '
                          WRITE(*,*)'the solution could not be '
                          WRITE(*,*)'computed.'
                          STOP
                        ENDIF
                        DEALLOCATE( PIVOT )
                      CASE ('DGESV')
                        ! ******************************************************************
                        !
                        ! **************************  SOLVER # 3  *************************
                        !
                        ! ******************************************************************
                        !
                        ! ========
                        ! DGESV.F:
                        ! ========
                        ! The program computes the solution to the system of linear
                        !  equations with a square matrix A and multiple
                        !  right-hand sides B, where A is the coefficient matrix:
                        ! REFERENCE:
                        ! 1. https://software.intel.com/sites/products/documentation/doclib/mkl_sa
                        ! /11/mkl_lapack_examples/dgesv_ex.f.htm
                        ! 
                        ! 2. http://faculty.washington.edu/rjl/uwamath583s11/sphinx/notes
                        ! /html/lapack_examples.html
                        !
                        NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                        NRHS = 1
                        LD_LHS = NN
                        LD_RHS = NN
                        ALLOCATE( PIVOT(NN) )
                        ALLOCATE( WORK(1) )     ! THIS IS A HACK

                        CALL DGESV(NN, NRHS, LHS, LD_LHS, PIVOT, RHS,   &
     &                             LD_RHS, INFO)

                        ! Note: The solution is returned in vector RHS and the array
                        ! LHS has been changed.

                        IF (INFO .GT. 0) THEN
                          WRITE(*,*)'ERROR!'
                          WRITE(*,*)'The diagonal element of triangular'
                          WRITE(*,*)'factor of LHS; '
                          WRITE(*,*)'LHS(',INFO,',',INFO,') is zero,'
                          WRITE(*,*)'so, that LHS is singular;'
                          WRITE(*,*)'the solution  could not be '
                          WRITE(*,*)'computed.'
                          STOP
                        ENDIF
                      DEALLOCATE( PIVOT )

                    CASE ('DGESVX')
                      ! ## -----------------------------------------   ##
                      ! ## Available EXPERT and RRR DRIVER routines:   ##
                      ! ## -----------------------------------------   ##
                      !
                      ! ******************************************************************
                      !
                      ! **************************  SOLVER # 4  *************************
                      !
                      ! ******************************************************************
                      !
                      ! ========
                      ! DGESVX.F
                      ! ========
                      ! The program computes the solution to the system of linear
                      ! equations with a square matrix A and multiple
                      ! right-hand sides B, where A is the coefficient matrix.
                      ! Additionally, we obtain error estimates for the solutions, information on
                      ! scaling, an estimate of the reciprocal of the
                      ! condition number of the scaled matrix A and an
                      ! estimate of the reciprocal of the pivot growth factor for 
                      ! the factorization of A are also output.
                      NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                      NMAX = NN
                      NRHS = 1
                      ALLOCATE( CC(NMAX), RR(NMAX) )
                      LDAF = NMAX
                      LDX = NMAX
                      LD_LHS = NN
                      LD_RHS = NN
                      ALLOCATE( XX(LDX, NRHS) )
                      ALLOCATE( BERR(NRHS), FERR(NRHS) )
                      ALLOCATE( AF(LDAF, NMAX) )
                      ALLOCATE( IWORK(NMAX) )
                      ALLOCATE( WORK(4*NMAX) )
                      ALLOCATE( PIVOT(NN) )
                
                      CALL DGESVX('E', 'N', NN,                         &
     &                       NRHS, LHS, LD_LHS, AF, LDAF, PIVOT, EQUED, &
     &                       RR, CC, RHS, LD_RHS, XX, LDX, RCOND, FERR, &
     &                       BERR, WORK, IWORK, INFO)
                      RHS = XX(:,1)

                      IF ((INFO .EQ. 0) .OR. (INFO .EQ. NN+1)) THEN
                        ! Print solution, error bounds, condition number, the form
                        ! of equilibration and the pivot growth factor
                        !IFAIL = 0
                        ! Ref:
                        ! X04CAF is an easy-to-use routine to print a real matrix
                        ! stored in a two-dimensional array.
                        ! http://www.nag.co.uk/numeric/fl/manual/xhtml/X04/x04caf.xml
                        !CALL X04CAF('General', ' ', NN, NRHS, XX, LDX,        &
!     &                        'Solution(s)', IFAIL)
                        ! This is not found in LAPACK/ATLAS

                        WRITE (*,*)
                        WRITE (*,*) 'Backward errors (m/c-dependent)'
                        WRITE (*,99996) (BERR(JJ),JJ =1, NRHS)
                        WRITE (*,*)
                        WRITE (*,*) 'Estimated forward error bounds'
                        WRITE (*,*) ' (m/C-dependent)'
                        WRITE (*,99996) (FERR(JJ),JJ=1, NRHS)
                        WRITE (*,*)

                        IF (EQUED .EQ. 'N') THEN
                          WRITE (*,*) 'A has not been equilibrated'
                        ELSEIF (EQUED .EQ. 'R') THEN
                          WRITE(*,*)'A has been row scaled as'
                          WRITE(*,*) 'diag(R)*A'
                        ELSEIF (EQUED .EQ. 'C') THEN
                          WRITE(*,*)'A has been column scaled:'
                          WRITE(*,*) 'A*diag(C)'
                        ELSEIF (EQUED .EQ. 'B') THEN
                          WRITE (*,*)'A has been row and column         &
     &                                scaled as diag(R)*A*diag(C)'
                        END IF
                        WRITE (*,*)
                        WRITE (*,*)'Reciprocal condition # estimate of'
                        WRITE (*,*)'the of scaled matrix: '
                        WRITE (*,*) RCOND
                        WRITE (*,*)
                        WRITE (*,*) 'Reciprocal pivot growth fact.'
                        WRITE (*,99996) WORK(1)

                        IF (INFO .EQ. NN+1) THEN
                          WRITE (*,*)
                          WRITE (*,*)'The matrix A is singular to'
                          WRITE (*,*)'working precision.'
                          WRITE(*,*)'The Results maybe inaccurate'
                          WRITE(*,*)
                        ENDIF
                      ELSE
                        WRITE (*,99997) 'The (', INFO, ',', INFO, ')',  &
     &                                ' element of the factor U is zero'
                        STOP
                      ENDIF

99996                 FORMAT ((3X,1P,7E11.1))
99997                 FORMAT (1X,A,I3,A,I3,A,A)
                      DEALLOCATE( PIVOT, CC, XX, RR )
                      DEALLOCATE( BERR, FERR, AF, IWORK )

                    CASE DEFAULT
                      WRITE(*,*)'SOLVERS OTHER THAN DSYSV, DPOSV, AND   &
     &                           DGESV ARE NOT USED TO SOLVE THE        &
     &                           NORMAL EQUATION. EXITING.'
                      STOP
                  END SELECT
                  
                CASE ('LLS')      ! LINEAR LEAST SQUARE -- SOLVERS ARE BETTER
                  IF (SOLVER .EQ. 'DSYSV') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DSYSV.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ELSEIF (SOLVER .EQ. 'DPOSV') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DPOSV.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ELSEIF (SOLVER .EQ. 'DGESV') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DGESV.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ELSEIF (SOLVER .EQ. 'DGESVX') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DGESVX.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ENDIF

                  !  ##############################################################  !
                  !  ####################### GOOD TECHNIQUES ######################  !
                  !  ##############################################################  !
                               

                  ! ## -----------------------------------------   ##
                  ! ##       EXPERT LEAST-SQUARE TYPE SOLVERS      ##
                  ! ## -----------------------------------------   ##

                  ! ******************************************************************
                  !
                  ! **************************  SOLVER # 1   *************************
                  !
                  ! ******************************************************************
                  !

                  ! The solver - 'DGELSS' doesn't converge for Legendre and
                  ! Bernstein; hence, the solver 'DGELS' is most
                  ! frequently used for all types of polynomial basis
                  ! functions
                  
                  SELECT CASE (SOLVER)
                    CASE ('DGELS')
                      ! ========
                      ! DGELS.F:
                      ! ========
                      ! Program computes the least squares solution to the overdetermined linear
                      ! system A*X = B with full rank matrix A using QR factorization,
                      ! where A is the coefficient matrix.//
                      ! DGELS solves overdetermined or underdetermined real linear systems
                      ! involving an M-by-N matrix A, or its transpose, using a QR or LQ
                      ! factorization of A.  It is assumed that A has full rank.
                      ! REFERENCE:
                      ! https://software.intel.com/sites/products/documentation/doclib/mkl_sa
                      ! /11/mkl_lapack_examples/dgels_ex.f.htm
                      ! Parameters

                      NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                      LDA = SAMPLE_SIZE
                      LDB = LDA
                      MM = SAMPLE_SIZE     ! The number of rows of the matrix A. MM >= 0.
                      LWMAX = 100*NN
                      ALLOCATE( PIVOT(NN) )
                      ALLOCATE( WORK(LWMAX) )
                      ALLOCATE( RHS(SAMPLE_SIZE) )
                      RHS = U         
                      NRHS = 1
                      ALLOCATE( LHS(1,1) )            ! THIS IS A HACK

                      ! AMAT == (input/output) DOUBLE PRECISION array, dimension (LDA,N)
                      ! On entry, the M-by-N matrix A.
                      ! On exit,
                      ! if M >= N, A is overwritten by details of its QR
                      ! factorization as returned by DGEQRF;
                      ! if M <  N, A is overwritten by details of its LQ
                      ! factorization as returned by DGELQF.

                      ! Query the optimal workspace.
                      LWORK = -1      ! If LWORK = -1, then a workspace query is assumed; the routine
                                      ! only calculates the optimal size of the WORK array, returns
                                      ! this value as the first entry of the WORK array, and no error
                                      ! message related to LWORK is issued by XERBLA.
                      CALL DGELS('N', MM, NN, NRHS, AMAT, LDA, RHS, LDB,&
     &                           WORK, LWORK, INFO)

                      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

                      ! Solve the equations A*X = B.
                      ! OR in other words -- 
                      ! Solve the least squares problem min( norm2(b - Ax) ) for x

                      CALL DGELS('N', MM, NN, NRHS, AMAT, LDA, RHS, LDB,&
     &                           WORK, LWORK, INFO)
                      
                      IF (INFO .GT. 0) THEN
                        WRITE(*,*)'The diagonal element ',INFO,' of '
                        WRITE(*,*)'the triangular  factor of A is zero,'
                        WRITE(*,*)' so that A does not have full rank;'
                        WRITE(*,*)' the least-squares solution could' 
                        WRITE(*,*)' not be computed.'
                        STOP
                      ENDIF

                      DEALLOCATE( PIVOT )

                    CASE ('DGELSS')
                      ! ******************************************************************
                      !
                      ! **************************  SOLVER # 2  *************************
                      !
                      ! ******************************************************************
                      ! PROF. BAI: "In addition, if you anticipate A might be extremely 
                      ! ill-conditioned, say near 10^16, then safer (and more
                      ! expensive) solver is DGELSS."
                      ! =========
                      ! DGELSS.F
                      ! =========
                      ! DGELSS - computes the minimum norm solution to a real linear
                      ! least squares problem:
                      ! Minimize 2-norm(| b - A*x |).
                      ! using the singular value decomposition (SVD) of A. A is an
                      ! M-by-N matrix which may be rank-deficient.
                      ! Several right hand side vectors b and solution vectors x can
                      ! be handled in a single call; they are stored as the columns
                      ! of the M-by-NRHS right hand side matrix B and the N-by-NRHS
                      ! solution matrix X.

                      ! The effective rank of A is determined by treating as zero
                      ! those singular values which are less than RCOND times the
                      ! largest singular value.
                      ! 
                      ! Parameters
                      NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                      NMAX = NN
                      LDA = SAMPLE_SIZE
                      LDB = LDA
                      MM = SAMPLE_SIZE     ! The number of rows of the matrix A. MM >= 0.
                      MMAX = MM
                      LWORK = 3*NMAX + 4*MMAX*(MMAX + NMAX)
                      ALLOCATE( SS(NMAX), WORK(LWORK) )
                      ALLOCATE( RHS(SAMPLE_SIZE) )
                      RHS = U         
                      NRHS = 1
                      ALLOCATE( LHS(1,1) )            ! THIS IS A HACK

                      ! Choose RCOND to reflect the relative accuracy of the input data
                      RCOND = 1.0D-16     ! Keep it as small as possible
                
                      ! Solve the least squares problem min( norm2(b - Ax) ) for the x
                      ! of minimum norm.
                      !
                      CALL DGELSS(MM, NN, NRHS, AMAT, LDA, RHS, LDB, SS,&
     &                      RCOND, RANK, WORK, LWORK, INFO)

                      IF (INFO .EQ. 0) THEN
                        WRITE (*,*)
                        WRITE (*,*) 'Tolerance used to estimate the rank&
     &                              of A'
                        WRITE (*,99998) RCOND
                        WRITE (*,*) 'Estimated rank of A'
                        WRITE (*,99994) RANK
                      ELSE
                        WRITE (*,*) 'SVD algorithm failed to converge'
                        STOP
                      ENDIF
                      DEALLOCATE( SS )

99998                 FORMAT (3X,1P,E11.2)
99994                 FORMAT (1X,I6)

                    CASE DEFAULT
                    PRINT *, 'NO OTHER LEAST-SQUARE SOLVERS CAN BE '
                    PRINT *, 'USED TO SOLVE THIS PROBLEM.'
                    STOP
                  END SELECT

                CASE DEFAULT
                  PRINT *, 'NO OTHER PROBLEM FORMULATION EXISTS'
                  PRINT *, 'EXCEPT THE NORM. EQN. AND LLS SCHEMES.'
                  PRINT *, 'EXITING!'
                  STOP
                END SELECT

                ! POST-PROCESSING STARTS
                ALLOCATE( C_SOL_ARR((M+1),R) )          ! need to deallocate
                ! the first NN elements of RHS are overwritten by the
                ! solution vector (LS solution)
                C_SOL_ARR = RESHAPE(RHS(1:NN), (/ (M+1), R /) )
                ALLOCATE( NORM_UJL(R) )                 ! need to deallocate
                ALLOCATE( U_H(SAMPLE_SIZE, R) )         ! need to deallocate
        
                !$OMP DO
                DO L = 1, R
                  !$OMP DO
                  DO J = 1, SAMPLE_SIZE
                    ! For the gien dimension k == dim_opt
                    CALL  U_HAT_1D(U_H(J,L), C_SOL_ARR(:,L), M,         &
     &                           POLY_TYPE, (/ COORD(J,DIM_OPT) /),     &
     &                           LOW_LIM, UP_LIM)
                  ENDDO
                  !$OMP ENDDO
                ENDDO
                !$OMP ENDDO

                !$OMP DO
                DO L = 1, R
                  NORM_UJL(L) = SEMI_NORM(U_H(:,L), SAMPLE_SIZE)
                  C_SOL_ARR(:,L) = C_SOL_ARR(:,L)/NORM_UJL(L)
                  S(L) = S(L)*NORM_UJL(L)               ! C changed to use
                                                        ! in next dimension K
                ENDDO
                !$OMP ENDDO
        
                C(:,:,DIM_OPT) = TRANSPOSE(C_SOL_ARR)   ! C changed to use
                                                        ! in next dimension K
                ! POST-PROCESSING ENDS
                !$OMP END PARALLEL
                DEALLOCATE( AMAT, C_SOL_ARR, NORM_UJL, U_H )
                DEALLOCATE( RHS, LHS ) ! DEALLOCATE( SS )
                DEALLOCATE( DIM_ARR )
                DEALLOCATE( WORK )

              !************************************************************!  
              ! Closes the alternating direction DO LOOP K = 1, D          !
              !************************************************************!

              ENDDO

              WRITE(*,*)' '
              WRITE(*,*)'ITERATION #', ITER, 'OF RANK: ', R, 'ENDS'
              WRITE(*,*)'FOR ALL THE DIMENSIONS'
              WRITE(*,*)' '

              ALLOCATE( F_SEP_ITER(1:SAMPLE_SIZE) )
              F_SEP_ITER = 0.D0
              !$OMP PARALLEL
              !$OMP DO 
              DO J = 1, SAMPLE_SIZE
                CALL SEP_RESP(F_SEP_ITER(J), C, S, (/ COORD(J,:) /),    &
     &                        POLY_TYPE, R, M, D, LOW_LIM, UP_LIM)
              ENDDO
              !$OMP ENDDO
              !$OMP END PARALLEL
                
              ! Computes the absolute value of change in error (decrease or increase)
              ! at ITER. We need to have a certain "DROP in error" at
              ! every iteration. If the "error drop" is not sufficiently
              ! small (governed by TOL) then we move on to the next 
              ! ITER <= MAX_ITER for the same given rank.

              IF (ITER .EQ. 1) THEN
                ERR_NEW = (SEMI_NORM((U - F_SEP_ITER), SAMPLE_SIZE))**2
                ERR_TOL = ABS(ERR_NEW)
              ELSE
                ERR_OLD = ERR_NEW
                ERR_NEW = (SEMI_NORM((U - F_SEP_ITER), SAMPLE_SIZE))**2
                ERR_TOL = ABS(ERR_NEW - ERR_OLD)
              ENDIF

              DEALLOCATE( F_SEP_ITER )
              
              
              IF (ERR_TOL .LE. TOL) THEN
                WRITE(*,*)' '
                WRITE(*,*)'SUCCESS! ERROR AT ITER #', ITER, 'DECREASED  &
     &                     SUFFICIENTLY W.R.T. ITER #', (ITER - 1)
                WRITE(*,*)' '
              ELSE
                WRITE(*,*)' '
                WRITE(*,*)'ERROR IN ITER # ',ITER, 'IS: ', ERR_TOL
                WRITE(*,*)' '
                WRITE(*,*)'NEED TO ACHIEVE AT LEAST: ', TOL
                WRITE(*,*)'HENCE MOVING ON TO THE NEXT ITERATION...'
              ENDIF

            !************************************************************!  
            !     CLOSES THE INNER IF LOOP (ITER .GT. MAX_ITER)          !
            !************************************************************!

            ENDIF

          !************************************************************!  
          !     CLOSES THE INNER DO-WHILE LOOP (ERR_TOL .GT. TOL)      !
          !************************************************************!

          ENDDO

  998     CONTINUE
          !$OMP PARALLEL
          !$OMP DO 
          DO J = 1, SAMPLE_SIZE
            CALL SEP_RESP(F_SEP_APPRX(J), C, S, (/ COORD(J,:) /),       &
     &                     POLY_TYPE, R, M, D, LOW_LIM, UP_LIM)
          ENDDO
          !$OMP ENDDO
          !$OMP END PARALLEL

          SEMI_NRM = (SEMI_NORM((U - F_SEP_APPRX), SAMPLE_SIZE))**2
                  
          IF (SEMI_NRM .LE. EPS) THEN
            WRITE(*,*)' '
            PRINT *, 'Solution converged at rank: ', R
            WRITE(*,*)' '
            FLAG_CONV = .TRUE.
            R_OUT = R
            ERR_DECAY = SEMI_NRM
            !EXIT
            RETURN
          ELSE
            FLAG_CONV = .FALSE.
            WRITE(*,*)' '
            WRITE(*,*)'Solution DID NOT CONVERGE at rank: ', R
            WRITE(*,*)' '

            ! Incrementing R by 1 if not convergered in the last iteration
            R = R + 1
          IF (R .LE. R_MAX) THEN
            WRITE(*,*)' '
            WRITE(*,*)'Moving on to the next rank:', R
            WRITE(*,*)' '
          ENDIF
        ENDIF

        !************************************************************!  
        !        ENDS THE MAIN DO-WHILE LOOP (R .LE. R_MAX)          !
        !************************************************************!
          
        ENDDO

        IF (FLAG_CONV .EQV. .FALSE.) THEN
          R_OUT = R_MAX
          ERR_DECAY = SEMI_NRM
        ENDIF

      END SUBROUTINE ALS_ORTH

      !==================================================================
      !==================================================================

      SUBROUTINE ALS_BRN(ERR_DECAY, R_OUT, C, S, D, M, R_MAX, C_MIN,    &
     &               C_MAX, U, SAMPLE_SIZE, COORD, EPS, LOW_LIM,        &
     &               UP_LIM, RUNDATAFILE, SOLVER, PRBLM_TYPE,           &
     &               FLAG_CONV, MAX_ITER, T_REGU_STAT, LAMBDA)
        USE NRTYPE
        USE OMP_LIB
        USE, INTRINSIC :: ISO_Fortran_env
        IMPLICIT NONE

        ! INTERFACE
        INTERFACE
          SUBROUTINE PRINT_MATRIX(M, N, A)
            USE NRTYPE
            INTEGER(I4B), INTENT(IN)          :: M, N
            REAL(DP), INTENT(IN)              :: A(M,N)
          END SUBROUTINE PRINT_MATRIX
        END INTERFACE

        INTERFACE
          SUBROUTINE TNSR_UHAT_D_DIM_BRN(PROD_U_HAT, NODES, M, D,       &
     &                                   C_RED,LOW_LIM, UP_LIM)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: M, D
            REAL(DP), INTENT(IN)            :: C_RED(1:(M+1),1:(D-1))
            REAL(DP), INTENT(IN)            :: NODES(1:(D-1))
            REAL(DP), INTENT(IN)            :: LOW_LIM, UP_LIM
            REAL(DP), INTENT(OUT)           :: PROD_U_HAT
          END SUBROUTINE TNSR_UHAT_D_DIM_BRN
        END INTERFACE

        INTERFACE
          FUNCTION SEMI_NORM(U, N)
            USE NRTYPE
            REAL(DP)                        :: SEMI_NORM
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(IN)            :: U(1:N)
          END FUNCTION SEMI_NORM
        END INTERFACE

        INTERFACE
          SUBROUTINE SEP_RESP_BRN(F_SEP, C, S, Y, R, M, D, LOW_LIM,     &
     &                            UP_LIM)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)         :: R, M, D
            REAL(DP), INTENT(IN)             :: C(1:R, 1:(M+1), 1:D)
            REAL(DP), INTENT(IN)             :: S(1:R), Y(1:D)
            REAL(DP), INTENT(IN)             :: LOW_LIM, UP_LIM
            REAL(DP), INTENT(OUT)            :: F_SEP
          END SUBROUTINE SEP_RESP_BRN
        END INTERFACE

        INTERFACE
          SUBROUTINE DIAG(DIAG_MAT, N)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: N
            REAL(DP), INTENT(OUT)           :: DIAG_MAT(N,N)
          END SUBROUTINE DIAG
        END INTERFACE

        INTERFACE
          SUBROUTINE KRON(A_KRON_B, A, M, N, B, P, Q)
            USE NRTYPE
            USE OMP_LIB
            INTEGER(I4B), INTENT(IN)        :: M, N, P, Q
            REAL(DP), INTENT(IN)            :: A(M,N), B(P,Q)
            REAL(DP), INTENT(OUT)           :: A_KRON_B(M*P,N*Q)
          END SUBROUTINE KRON
        END INTERFACE

        ! Interfacing part of the dummy variables
        INTEGER(I4B), INTENT(IN)               :: D, M, R_MAX
        INTEGER(I4B), INTENT(IN)               :: SAMPLE_SIZE, MAX_ITER
        REAL(DP), ALLOCATABLE, INTENT(OUT)     :: C(:,:,:)
        INTEGER(I4B), INTENT(OUT)              :: R_OUT
        LOGICAL, INTENT(OUT)                   :: FLAG_CONV
        REAL(DP), INTENT(IN)                   :: C_MIN, C_MAX
        REAL(DP), ALLOCATABLE, INTENT(OUT)     :: S(:)
        REAL(DP), INTENT(OUT)                  :: ERR_DECAY
        REAL(DP), INTENT(IN)                   :: U(SAMPLE_SIZE)
        REAL(DP), INTENT(IN)                   :: COORD(SAMPLE_SIZE,D)
        REAL(DP), INTENT(IN)                   :: EPS, LAMBDA
        REAL(DP), INTENT(IN)                   :: LOW_LIM, UP_LIM
        CHARACTER(80), INTENT(IN)              :: RUNDATAFILE, SOLVER
        CHARACTER(80), INTENT(IN)              :: PRBLM_TYPE
        CHARACTER(80), INTENT(IN)              :: T_REGU_STAT
        
        ! LOCAL VARIABLES:
        INTEGER(I4B)            :: R
        INTEGER(I4B)            :: K, L, J, J1, J2, ALPHA, JJ
        INTEGER(I4B)            :: DIM_OPT
        REAL(DP), ALLOCATABLE   :: AMAT(:,:)
        REAL(DP), ALLOCATABLE   :: C_SOL_ARR(:,:), C_RED(:,:)
        REAL(DP), ALLOCATABLE   :: U_H(:,:), DIM_ARR(:),  NODES(:)
        REAL(DP), ALLOCATABLE   :: NORM_UJL(:), F_SEP_APPRX(:)
        REAL(DP)                :: PROD_U_HAT
        REAL(DP), ALLOCATABLE   :: RHS(:), PNX(:), LHS(:,:)
        REAL(DP)                :: SEMI_NRM
        REAL(DP)                :: ERR_OLD, ERR_NEW, ERR_TOL, TOL
        INTEGER(I4B)            :: ITER
        REAL(DP), ALLOCATABLE   :: F_SEP_ITER(:)
        REAL(DP), ALLOCATABLE   :: T_REG_MAT(:,:), S_MAT(:,:)
        ! Parameters used by LAPACK equation solver:                                                
        INTEGER(I4B)                            :: LDA, LDB, INFO
        INTEGER(I4B)                            :: MMAX, NMAX, RANK
        INTEGER(I4B), ALLOCATABLE               :: PIVOT(:)
        INTEGER(I4B)                            :: NRHS, NN, MM
        INTEGER(I4B)                            :: LWMAX, LWORK
        REAL(DP), ALLOCATABLE                   :: WORK(:), SS(:)
        REAL(DP), ALLOCATABLE                   :: AF(:,:)
        REAL(DP), ALLOCATABLE                   :: IDENTITY(:,:)
        REAL(DP), ALLOCATABLE                   :: BERR(:), FERR(:)
        REAL(DP), ALLOCATABLE                   :: CC(:), RR(:), XX(:,:)
        INTEGER(I4B), ALLOCATABLE               :: IWORK(:)
        REAL(DP)                                :: RCOND
        INTEGER(I4B)                            :: LD_LHS, LD_RHS
        CHARACTER(1)                            :: EQUED
        INTEGER(I4B)                            :: LDAF, LDX

        ! Keep this bracket reasonably big to get rid-off singular
        ! matrices while solving the liner system
        !C_MIN = -10.D0
        !C_MAX = 10.D0
        R = 1
        FLAG_CONV = .TRUE.
        !MAX_ITER = 20

        DO WHILE (R .LE. R_MAX)
          !
          IF (FLAG_CONV .EQV. .FALSE.) THEN
            DEALLOCATE( F_SEP_APPRX, C, S )
          ENDIF

          ALLOCATE( F_SEP_APPRX(SAMPLE_SIZE) )
           ! Re-initializing C and S
          ALLOCATE( C(1:R, 1:(M + 1), 1:D) )
          ALLOCATE( S(1:R) )

          ! The following set of data remains the same for all dimensions in one
          ! iteration, i.e., for a given rank. When the rank is incremented by 1,
          ! these values again change/re-initialized, irrespective of the values from
          ! last rank-iteration. 
          ! So, we actually relinquish all the effort put in
          ! to calculate the "optimum solution" for the lower rank. That's bad!

          ! EX: j = n + FLOOR((m+1-n)*u)  ! We want to choose one from m-n+1 integers
          ! The intrinsic random_number(u) returns a real number u (or an array of such) 
          ! from the uniform distribution over the interval [0,1). 
          ! [That is, it includes 0 but not 1.]
          ! To have a discrete uniform distribution on the integers {n, n+1, ..., m-1, m}
          ! carve the continuous distribution up into m+1-n equal sized chunks, mapping
          ! each chunk to an integer. One way could be:
          ! call random_number(u)
          ! As you can see, for the initial question for {0, 1, 2, 3, 4, 5} this reduces to
          ! call random_number(u)
          ! j = FLOOR(6*u)            ! n=0 and m=5

          CALL RANDOM_NUMBER(C)
          C = C_MIN + (C_MAX + 1.D0 - C_MIN)*C   ! This is in [C_MIN, C_MAX]
          !C = C_MIN + (C_MAX - C_MIN)*C   ! This is in [C_MIN, C_MAX)

          CALL RANDOM_NUMBER(S)

          TOL = EPS*1.D-06
          ERR_TOL = 10.D0*TOL
          ITER = 0

          !WRITE(*,*)' '
          !WRITE(*,*)'THE SOLUTION PROCESS STARTS FOR RANK = ', R
          !WRITE(*,*)' '
          
          DO WHILE (ERR_TOL .GT. TOL) 
            IF (ITER .GT. MAX_ITER) THEN
              !WRITE(*,*)' '
              !WRITE(*,*)'EXCEEDED MAXIMUM ALLOWED', MAX_ITER,           &
     !&                  'ITERATIONS FOR RANK: ', R
              !WRITE(*,*)' '
              GO TO 999
            ELSE
              ITER = ITER + 1
              !WRITE(*,*)' '
              !WRITE(*,*)'ITERATION #', ITER, 'OF RANK: ', R, 'STARTS...'
              !WRITE(*,*)' '

              DO K = 1, D
                DIM_OPT = K
                ! AMAT is initialized for every dimension. The calculated value of
                ! C ans S are included in the already initialized C and
                ! S for every new AMAT

                ! This is a rectangular system; overdetermined system
                ! solved by method of least squares.
                ALLOCATE( AMAT(SAMPLE_SIZE, R*(M + 1)) )   ! **** need to deallocate

                ! DIM_ARR(D - 1) -- Contains all the dimension except the one in
                ! opt-loop, i.e., DIM_OPT
                ALLOCATE( DIM_ARR(D-1) )                   ! **** need to deallocate
                IF (DIM_OPT .EQ. 1) THEN
                  DIM_ARR = (/ (J, J = 2, D) /)
                ELSEIF (DIM_OPT .EQ. D) THEN
                  DIM_ARR = (/ (J, J = 1, (D-1)) /)
                ELSE
                  DIM_ARR(1:(DIM_OPT-1)) = (/ (J1, J1 = 1, (DIM_OPT-1)) &
     &              /)
                  DIM_ARR(DIM_OPT:(D-1)) = (/ (J2, J2 = (DIM_OPT+1), D) &
     &             /)
                ENDIF

                !$OMP DO
                DO L = 1, R
                  ! C_RED((M + 1), (d - 1)) -- Contains all the coefficient at a
                  ! given rank except the coeffiecients which go to the opt-loop
                  ! Computed/Selected for each rank L
                  ALLOCATE( C_RED(1:(M+1), 1:(D - 1)) )   ! **** need to deallocate
                  IF (DIM_OPT .EQ. 1) THEN
                    C_RED = RESHAPE(C(L,:,2:D), (/ (M+1), (D-1) /))
                  ELSEIF (DIM_OPT .EQ. D) THEN
                    C_RED = RESHAPE(C(L,:,1:(D-1)), (/ (M+1), (D-1) /))
                  ELSE
                    C_RED(:,1:(DIM_OPT-1)) = C(L,:,1:(DIM_OPT-1))
                    C_RED(:,DIM_OPT:D-1) = C(L,:,(DIM_OPT+1):D)
                  ENDIF

                  !$OMP DO
                  DO J = 1, SAMPLE_SIZE
                    ! NODES(d -1) -- At a given sample point, it contains the all the
                    ! coordinate component of the sample point except the one, whose direction
                    ! goes to the opt-loop.
                    ALLOCATE( NODES(D-1) )                ! **** need to deallocate
                    IF (DIM_OPT .EQ. 1) THEN
                      NODES = COORD(J,2:D)
                    ELSEIF (DIM_OPT .EQ. D) THEN
                      NODES = COORD(J,1:D-1)
                    ELSE
                      NODES(1:DIM_OPT-1) = COORD(J,1:(DIM_OPT-1))
                      NODES(DIM_OPT:D-1) = COORD(J,(DIM_OPT+1):D)
                    ENDIF
                    
                    CALL TNSR_UHAT_D_DIM_BRN(PROD_U_HAT, NODES, M, D,   &
     &                                       C_RED, LOW_LIM, UP_LIM)

                    ALLOCATE( PNX(0:M))

                    CALL BERNSTEIN(PNX, M, LOW_LIM, UP_LIM,             &
     &                             COORD(J, DIM_OPT))
   
                    !$OMP DO    
                    DO ALPHA = 0, M
                      ! Product over all the dimensions except "dim_opt"
                      ! and then summed over l for a given node (j)                   
                      AMAT(J,((L-1)*(M+1)+(1+ALPHA))) =                 &
     &                                       S(L)*PNX(ALPHA)*PROD_U_HAT
                    ENDDO       ! ALPHA LOOP ENDS
                    !$OMP ENDDO

                    DEALLOCATE( NODES, PNX )

                  ENDDO         ! J LOOP OVER ALL NODES ENDS
                  !$OMP ENDDO

                  DEALLOCATE( C_RED )

                ENDDO       ! L LOOP ENDS

                SELECT CASE (PRBLM_TYPE)
                  CASE ('NRML_EQTN')
                    IF (SOLVER .EQ. 'DGELS') THEN
                      PRINT *, 'ERROR! ERROR!'
                      WRITE(*,*) 'THE NORMAL EQUATIONS CAN NOT BE SOLVED&
     &                 WITH LEAST SQUARE SOLVERS DGELS.'
                      PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                      PRINT *, 'THE FOLLOWING: '
                      PRINT *, '1. DSYSV'
                      PRINT *, '2. DPOSV'
                      PRINT *, '3. DGESV'
                      PRINT *, '4. DGESVX '
                      PRINT *, 'EXITING'
                      STOP
                    ELSEIF (SOLVER .EQ. 'DGELSS') THEN
                      PRINT *, 'ERROR! ERROR!'
                      WRITE(*,*) 'THE NORMAL EQUATIONS CAN NOT BE SOLVED&
     &                 WITH LEAST SQUARE SOLVERS DGELSS'
                      PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                      PRINT *, 'THE FOLLOWING: '
                      PRINT *, '1. DSYSV'
                      PRINT *, '2. DPOSV'
                      PRINT *, '3. DGESV'
                      PRINT *, '4. DGESVX '
                      PRINT *, 'EXITING'
                      STOP
                    ENDIF

                    !  ##############################################################  !
                    !  ####################### BAD TECHNIQUES #######################  !
                    !  ##############################################################  !
    
                    ! BADLY CONDITIONED (A**T*A) -- DON'T USE THIS:
                    ALLOCATE( LHS(R*(M+1), R*(M+1)) )
                    ALLOCATE( RHS(R*(M+1)) )
                    RHS = MATMUL(TRANSPOSE(AMAT),U)

                    IF (T_REGU_STAT .EQ. 'OFF') THEN
                      ! BADLY CONDITIONED (A**T*A) -- DON'T USE THIS:
                      ! *********** DON'T USE THIS ************
                      LHS = MATMUL(TRANSPOSE(AMAT),AMAT)      ! *** VERY HIGH CONDITION NUMBER
                    ELSEIF (T_REGU_STAT .EQ. 'ON') THEN
                      ! TIKHONOV REGULARIZATION WITH UNIT REGULARIZATION
                      ! MATRIX
                      ALLOCATE( IDENTITY((M+1), (M+1)) )
                      CALL DIAG(IDENTITY, (M+1))
                      ALLOCATE( S_MAT(R,1) )
                      S_MAT = RESHAPE(S, (/ R, 1 /))
                      ALLOCATE( T_REG_MAT(R*(M+1), R*(M+1)) )
                      CALL KRON(T_REG_MAT, S_MAT, R, 1, IDENTITY,       &
     &                      (M+1), (M+1))
                      LHS = MATMUL(TRANSPOSE(AMAT),AMAT) +              &
     &                      LAMBDA**2*T_REG_MAT
                      DEALLOCATE( IDENTITY, T_REG_MAT, S_MAT )
                    ENDIF                   

                    ! ## -------------------------------------------------------   ##
                    ! ## Available SIMPLE and DIVIDE AND CONQUER DRIVER routines:  ##
                    ! ## -------------------------------------------------------   ##
                    SELECT CASE (SOLVER)

                      CASE ('DSYSV')
                        ! ******************************************************************
                        !
                        ! **************************  SOLVER # 1   *************************
                        !
                        ! ******************************************************************
                        ! ========
                        ! DSYSV.F:
                        ! ========
                        ! The program computes the solution to the system of linear equations
                        ! with a real symmetric matrix A and multiple right-hand sides B,
                        ! where A is the coefficient matrix:
                        ! REFERENCE: 
                        ! https://software.intel.com/sites/products/documentation/doclib
                        ! /mkl_sa/11/mkl_lapack_examples/dsysv_ex.f.htm
                        NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                        NRHS = 1
                        LD_LHS = NN
                        LD_RHS = NN
                        LWMAX = 100*NN
                        ALLOCATE( PIVOT(NN) )
                        ALLOCATE( WORK(LWMAX) )
                        ! Query the optimal workspace.
                        LWORK = -1
                        CALL DSYSV('Lower', NN, NRHS, LHS, LD_LHS,      &
      &                           PIVOT, RHS, LD_RHS, WORK, LWORK, INFO)
                        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

                        ! Solve the SQUARE SYSTEM OF equations A*X = B.
                        CALL DSYSV('Lower', NN, NRHS, LHS, LD_LHS,      &
      &                           PIVOT, RHS, LD_RHS, WORK, LWORK, INFO )

                        ! Check for the exact singularity.
                        IF( INFO.GT.0 ) THEN
                          WRITE(*,*)'The element of the diagonal  '
                          WRITE(*,*)'factor D(',INFO,',',INFO,') is  '
                          WRITE(*,*)'zero, so that D is singular; the'
                          WRITE(*,*)'solution could not be computed.'
                          STOP
                        ENDIF
                        DEALLOCATE( PIVOT )

                      CASE ('DPOSV')
                        ! ******************************************************************
                        !
                        ! **************************  SOLVER # 2   *************************
                        !
                        ! ******************************************************************
                        !
                        ! ========
                        ! DPOSV.F
                        ! ========
                        ! The program computes the solution to the system of linear
                        ! equations with a symmetric positive-definite matrix A and multiple
                        ! right-hand sides B, where A is the coefficient matrix:
                        ! REFERENCE:
                        ! https://software.intel.com/sites/products/documentation/
                        ! doclib/mkl_sa/11/mkl_lapack_examples/dposv_ex.f.htm
                        NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                        NRHS = 1
                        LD_LHS = NN
                        LD_RHS = NN
                        ALLOCATE( PIVOT(NN) )
                        ALLOCATE( WORK(1) )     ! THIS IS A HACK
                        CALL DPOSV('Upper', NN, NRHS, LHS, LD_LHS, RHS, &
     &                             LD_RHS, INFO)
                        IF( INFO.GT.0 ) THEN
                          WRITE(*,*)'The leading minor of order ',INFO
                          WRITE(*,*)' is not positive definite; '
                          WRITE(*,*)'the solution could not be '
                          WRITE(*,*)'computed.'
                          STOP
                        ENDIF
                        DEALLOCATE( PIVOT )
                      CASE ('DGESV')
                        ! ******************************************************************
                        !
                        ! **************************  SOLVER # 3  *************************
                        !
                        ! ******************************************************************
                        !
                        ! ========
                        ! DGESV.F:
                        ! ========
                        ! The program computes the solution to the system of linear
                        !  equations with a square matrix A and multiple
                        !  right-hand sides B, where A is the coefficient matrix:
                        ! REFERENCE:
                        ! 1. https://software.intel.com/sites/products/documentation/doclib/mkl_sa
                        ! /11/mkl_lapack_examples/dgesv_ex.f.htm
                        ! 
                        ! 2. http://faculty.washington.edu/rjl/uwamath583s11/sphinx/notes
                        ! /html/lapack_examples.html
                        !
                        NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                        NRHS = 1
                        LD_LHS = NN
                        LD_RHS = NN
                        ALLOCATE( PIVOT(NN) )
                        ALLOCATE( WORK(1) )     ! THIS IS A HACK
                        CALL DGESV(NN, NRHS, LHS, LD_LHS, PIVOT, RHS,   &
     &                             LD_RHS, INFO)

                        ! Note: The solution is returned in vector RHS and the array
                        ! LHS has been changed.

                        IF (INFO .GT. 0) THEN
                          WRITE(*,*)'ERROR!'
                          WRITE(*,*)'The diagonal element of triangular'
                          WRITE(*,*)'factor of LHS; '
                          WRITE(*,*)'LHS(',INFO,',',INFO,') is zero,'
                          WRITE(*,*)'so, that LHS is singular;'
                          WRITE(*,*)'the solution  could not be '
                          WRITE(*,*)'computed.'
                          STOP
                        ENDIF
                      DEALLOCATE( PIVOT )

                    CASE ('DGESVX')
                      ! ## -----------------------------------------   ##
                      ! ## Available EXPERT and RRR DRIVER routines:   ##
                      ! ## -----------------------------------------   ##
                      !
                      ! ******************************************************************
                      !
                      ! **************************  SOLVER # 4  *************************
                      !
                      ! ******************************************************************
                      !
                      ! ========
                      ! DGESVX.F
                      ! ========
                      ! The program computes the solution to the system of linear
                      ! equations with a square matrix A and multiple
                      ! right-hand sides B, where A is the coefficient matrix.
                      ! Additionally, we obtain error estimates for the solutions, information on
                      ! scaling, an estimate of the reciprocal of the
                      ! condition number of the scaled matrix A and an
                      ! estimate of the reciprocal of the pivot growth factor for 
                      ! the factorization of A are also output.
                      NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                      NMAX = NN
                      NRHS = 1
                      ALLOCATE( CC(NMAX), RR(NMAX) )
                      LDAF = NMAX
                      LDX = NMAX
                      LD_LHS = NN
                      LD_RHS = NN
                      ALLOCATE( XX(LDX, NRHS) )
                      ALLOCATE( BERR(NRHS), FERR(NRHS) )
                      ALLOCATE( AF(LDAF, NMAX) )
                      ALLOCATE( IWORK(NMAX) )
                      ALLOCATE( WORK(4*NMAX) )
                      ALLOCATE( PIVOT(NN) )
                
                      CALL DGESVX('E', 'N', NN,                         &
     &                       NRHS, LHS, LD_LHS, AF,  LDAF, PIVOT, EQUED,&
     &                       RR, CC, RHS, LD_RHS, XX, LDX, RCOND, FERR, &
     &                       BERR, WORK, IWORK, INFO)

                      RHS = XX(:,1)

                      IF ((INFO .EQ. 0) .OR. (INFO .EQ. NN+1)) THEN
                        ! Print solution, error bounds, condition number, the form
                        ! of equilibration and the pivot growth factor
                        !IFAIL = 0
                        ! Ref:
                        ! X04CAF is an easy-to-use routine to print a real matrix
                        ! stored in a two-dimensional array.
                        ! http://www.nag.co.uk/numeric/fl/manual/xhtml/X04/x04caf.xml
                        !CALL X04CAF('General', ' ', NN, NRHS, XX, LDX,        &
!     &                        'Solution(s)', IFAIL)
                        ! This is not found in LAPACK/ATLAS
                        WRITE (*,*)
                        WRITE (*,*) 'Backward errors (m/c-dependent)'
                        WRITE (*,99999) (BERR(JJ),JJ =1, NRHS)
                        WRITE (*,*)
                        WRITE (*,*) 'Estimated forward error bounds'
                        WRITE (*,*) ' (m/C-dependent)'
                        WRITE (*,99999) (FERR(JJ),JJ=1, NRHS)
                        WRITE (*,*)

                        IF (EQUED .EQ. 'N') THEN
                          WRITE (*,*) 'A has not been equilibrated'
                        ELSEIF (EQUED .EQ. 'R') THEN
                          WRITE(*,*)'A has been row scaled as'
                          WRITE(*,*) 'diag(R)*A'
                        ELSEIF (EQUED .EQ. 'C') THEN
                          WRITE(*,*)'A has been column scaled:'
                          WRITE(*,*) 'A*diag(C)'
                        ELSEIF (EQUED .EQ. 'B') THEN
                          WRITE (*,*)'A has been row and column         &
     &                                scaled as diag(R)*A*diag(C)'
                        END IF
                        WRITE (*,*)
                        WRITE (*,*)'Reciprocal condition # estimate of'
                        WRITE (*,*)'the of scaled matrix: '
                        WRITE (*,*) RCOND
                        WRITE (*,*)
                        WRITE (*,*) 'Reciprocal pivot growth fact.'
                        WRITE (*,99999) WORK(1)

                        IF (INFO .EQ. NN+1) THEN
                          WRITE (*,*)
                          WRITE (*,*)'The matrix A is singular to'
                          WRITE (*,*)'working precision.'
                          WRITE(*,*)'The Results maybe inaccurate'
                          WRITE(*,*)
                        ENDIF
                      ELSE
                        WRITE (*,99991) 'The (', INFO, ',', INFO, ')',  &
     &                                ' element of the factor U is zero'
                        STOP
                      ENDIF

99999                 FORMAT ((3X,1P,7E11.1))
99991                 FORMAT (1X,A,I3,A,I3,A,A)
                      DEALLOCATE( PIVOT, CC, XX, RR )
                      DEALLOCATE( BERR, FERR, AF, IWORK )

                    CASE DEFAULT
                      WRITE(*,*)'SOLVERS OTHER THAN DSYSV, DPOSV, AND   &
     &                           DGESV ARE NOT USED TO SOLVE THE        &
     &                           NORMAL EQUATION. EXITING.'
                      STOP
                  END SELECT
                  
                CASE ('LLS')      ! LINEAR LEAST SQUARE -- SOLVERS ARE BETTER
                  IF (SOLVER .EQ. 'DSYSV') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DSYSV.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ELSEIF (SOLVER .EQ. 'DPOSV') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DPOSV.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ELSEIF (SOLVER .EQ. 'DGESV') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DGESV.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ELSEIF (SOLVER .EQ. 'DGESVX') THEN
                    PRINT *, 'ERROR! ERROR!'
                    WRITE(*,*) 'THE LEAST SQUARE PROBLEM CAN NOT BE     &
     &              SOLVED WITH SOLVERS: DGESVX.'
                    PRINT *, 'CONSIDER CHANGING THE SOLVER TO ONE OF'
                    PRINT *, 'THE FOLLOWING: '
                    PRINT *, '1. DGELS'
                    PRINT *, '2. DGELSS'
                    PRINT *, 'EXITING'
                    STOP
                  ENDIF

                  !  ##############################################################  !
                  !  ####################### GOOD TECHNIQUES ######################  !
                  !  ##############################################################  !
                               

                  ! ## -----------------------------------------   ##
                  ! ##       EXPERT LEAST-SQUARE TYPE SOLVERS      ##
                  ! ## -----------------------------------------   ##

                  ! ******************************************************************
                  !
                  ! **************************  SOLVER # 1   *************************
                  !
                  ! ******************************************************************
                  !

                  SELECT CASE (SOLVER)
                    CASE ('DGELS')
                      ! ========
                      ! DGELS.F:
                      ! ========
                      ! Program computes the least squares solution to the overdetermined linear
                      ! system A*X = B with full rank matrix A using QR factorization,
                      ! where A is the coefficient matrix.//
                      ! DGELS solves overdetermined or underdetermined real linear systems
                      ! involving an M-by-N matrix A, or its transpose, using a QR or LQ
                      ! factorization of A.  It is assumed that A has full rank.
                      ! REFERENCE:
                      ! https://software.intel.com/sites/products/documentation/doclib/mkl_sa
                      ! /11/mkl_lapack_examples/dgels_ex.f.htm
                      ! Parameters
                      NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                      LDA = SAMPLE_SIZE
                      LDB = LDA
                      MM = SAMPLE_SIZE     ! The number of rows of the matrix A. MM >= 0.
                      LWMAX = 100*NN
                      ALLOCATE( PIVOT(NN) )
                      ALLOCATE( WORK(LWMAX) )
                      ALLOCATE( RHS(SAMPLE_SIZE) )
                      RHS = U         
                      NRHS = 1
                      ALLOCATE( LHS(1,1) )            ! THIS IS A HACK

                      ! AMAT == (input/output) DOUBLE PRECISION array, dimension (LDA,N)
                      ! On entry, the M-by-N matrix A.
                      ! On exit,
                      ! if M >= N, A is overwritten by details of its QR
                      ! factorization as returned by DGEQRF;
                      ! if M <  N, A is overwritten by details of its LQ
                      ! factorization as returned by DGELQF.

                      ! Query the optimal workspace.
                      LWORK = -1      ! If LWORK = -1, then a workspace query is assumed; the routine
                                      ! only calculates the optimal size of the WORK array, returns
                                      ! this value as the first entry of the WORK array, and no error
                                      ! message related to LWORK is issued by XERBLA.
                      CALL DGELS('N', MM, NN, NRHS, AMAT, LDA, RHS, LDB,&
     &                           WORK, LWORK, INFO)

                      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

                      ! Solve the equations A*X = B.
                      ! OR in other words -- 
                      ! Solve the least squares problem min( norm2(b - Ax) ) for x

                      CALL DGELS('N', MM, NN, NRHS, AMAT, LDA, RHS, LDB,&
     &                           WORK, LWORK, INFO)
                      IF (INFO .GT. 0) THEN
                        WRITE(*,*)'The diagonal element ',INFO,' of '
                        WRITE(*,*)'the triangular  factor of A is zero,'
                        WRITE(*,*)' so that A does not have full rank;'
                        WRITE(*,*)' the least-squares solution could' 
                        WRITE(*,*)' not be computed.'
                        STOP
                      ENDIF

                      DEALLOCATE( PIVOT )

                    CASE ('DGELSS')
                      ! ******************************************************************
                      !
                      ! **************************  SOLVER # 2  *************************
                      !
                      ! ******************************************************************
                      ! PROF. BAI: "In addition, if you anticipate A might be extremely 
                      ! ill-conditioned, say near 10^16, then safer (and more
                      ! expensive) solver is DGELSS."
                      ! =========
                      ! DGELSS.F
                      ! =========
                      ! DGELSS - computes the minimum norm solution to a real linear
                      ! least squares problem:
                      ! Minimize 2-norm(| b - A*x |).
                      ! using the singular value decomposition (SVD) of A. A is an
                      ! M-by-N matrix which may be rank-deficient.
                      ! Several right hand side vectors b and solution vectors x can
                      ! be handled in a single call; they are stored as the columns
                      ! of the M-by-NRHS right hand side matrix B and the N-by-NRHS
                      ! solution matrix X.

                      ! The effective rank of A is determined by treating as zero
                      ! those singular values which are less than RCOND times the
                      ! largest singular value.
                      ! 
                      ! Parameters
                      NN = R*(M+1)         ! The number of columns of the matrix A. NN >= 0.
                      NMAX = NN
                      LDA = SAMPLE_SIZE
                      LDB = LDA
                      MM = SAMPLE_SIZE     ! The number of rows of the matrix A. MM >= 0.
                      MMAX = MM
                      LWORK = 3*NMAX + 4*MMAX*(MMAX + NMAX)
                      ALLOCATE( SS(NMAX), WORK(LWORK) )
                      ALLOCATE( RHS(SAMPLE_SIZE) )
                      RHS = U         
                      NRHS = 1
                      ALLOCATE( LHS(1,1) )            ! THIS IS A HACK

                      ! Choose RCOND to reflect the relative accuracy of the input data
                      RCOND = 1.0D-16     ! Keep it as small as possible
                
                      ! Solve the least squares problem min( norm2(b - Ax) ) for the x
                      ! of minimum norm.
                      !
                      CALL DGELSS(MM, NN, NRHS, AMAT, LDA, RHS, LDB, SS,&
     &                      RCOND, RANK, WORK, LWORK, INFO)

                      IF (INFO .EQ. 0) THEN
                        WRITE (*,*)
                        WRITE (*,*) 'Tolerance used to estimate the rank&
     &                              of A'
                        WRITE (*,99998) RCOND
                        WRITE (*,*) 'Estimated rank of A'
                        WRITE (*,99997) RANK
                      ELSE
                        WRITE (*,*) 'SVD algorithm failed to converge'
                        STOP
                      ENDIF
                      DEALLOCATE( SS )

99998                 FORMAT (3X,1P,E11.2)
99997                 FORMAT (1X,I6)

                    CASE DEFAULT
                    PRINT *, 'NO OTHER LEAST-SQUARE SOLVERS CAN BE '
                    PRINT *, 'USED TO SOLVE THIS PROBLEM.'
                    STOP
                  END SELECT

                CASE DEFAULT
                  PRINT *, 'NO OTHER PROBLEM FORMULATION EXISTS'
                  PRINT *, 'EXCEPT THE NORM. EQN. AND LLS SCHEMES.'
                  PRINT *, 'EXITING!'
                  STOP
                END SELECT

                ! POST-PROCESSING STARTS
                ALLOCATE( C_SOL_ARR((M+1),R) )          ! need to deallocate
                ! the first NN elements of RHS are overwritten by the
                ! solution vector (LS solution)
                C_SOL_ARR = RESHAPE(RHS(1:NN), (/ (M+1), R /) )
                ALLOCATE( NORM_UJL(R) )                 ! need to deallocate
                ALLOCATE( U_H(SAMPLE_SIZE, R) )         ! need to deallocate
        
                !$OMP DO
                DO L = 1, R
                  !$OMP DO
                  DO J = 1, SAMPLE_SIZE
                    ! For the gien dimension k == dim_opt
                    CALL U_HAT_1D_BRN(U_H(J,L),  C_SOL_ARR(:,L), M,     &
     &                                COORD(J,DIM_OPT), LOW_LIM, UP_LIM)
                  ENDDO
                  !$OMP ENDDO
                ENDDO
                !$OMP ENDDO

                !$OMP DO
                DO L = 1, R
                  NORM_UJL(L) = SEMI_NORM(U_H(:,L), SAMPLE_SIZE)
                  C_SOL_ARR(:,L) = C_SOL_ARR(:,L)/NORM_UJL(L)
                  S(L) = S(L)*NORM_UJL(L)               ! C changed to use
                                                        ! in next dimension K
                ENDDO
                !$OMP ENDDO
        
                C(:,:,DIM_OPT) = TRANSPOSE(C_SOL_ARR)   ! C changed to use
                                                    ! in next dimension K
                ! POST-PROCESSING ENDS
                !$OMP END PARALLEL
                DEALLOCATE( AMAT, C_SOL_ARR, NORM_UJL, U_H )
                DEALLOCATE( RHS, LHS ) ! DEALLOCATE( SS )
                DEALLOCATE( DIM_ARR )
                DEALLOCATE( WORK )
              ENDDO         ! K LOOP ENDS

              !WRITE(*,*)' '
              !WRITE(*,*)'ITERATION #', ITER, 'OF RANK: ', R, 'ENDS'
              !WRITE(*,*)'FOR ALL THE DIMENSIONS'
              !WRITE(*,*)' '


              ALLOCATE( F_SEP_ITER(1:SAMPLE_SIZE) )
              !$OMP PARALLEL
              !$OMP DO 
              DO J = 1, SAMPLE_SIZE
                CALL SEP_RESP_BRN(F_SEP_ITER(J), C, S,  COORD(J,:), R,  &
     &                        M, D, LOW_LIM, UP_LIM)
              ENDDO
              !$OMP ENDDO
              !$OMP END PARALLEL

              IF (ITER .EQ. 1) THEN
                ERR_NEW = (SEMI_NORM((U - F_SEP_ITER), SAMPLE_SIZE))**2
                ERR_OLD = 0.D0
              ELSE
                ERR_OLD = ERR_NEW
                ERR_NEW = (SEMI_NORM((U - F_SEP_ITER), SAMPLE_SIZE))**2
              ENDIF

              DEALLOCATE( F_SEP_ITER )
              ERR_TOL = ABS(ERR_NEW - ERR_OLD)
              !WRITE(*,*)' '
              !WRITE(*,*)'ERROR IN ITER # ',ITER, 'IS: ', ERR_TOL
              !WRITE(*,*)' '
              !WRITE(*,*)'NEED TO ACHIEVE AT LEAST: ', TOL
              !WRITE(*,*)' '
              IF (ERR_TOL .LE. TOL) THEN
                !WRITE(*,*)' '
                !WRITE(*,*)'SUCCESS! ITERATIONS END.'
                !WRITE(*,*)' '
              ENDIF
            ENDIF
          ENDDO         ! CLOSES THE INNER WHILE LOOP

  999     CONTINUE
          !$OMP PARALLEL
          !$OMP DO 
          DO J = 1, SAMPLE_SIZE
            CALL SEP_RESP_BRN(F_SEP_APPRX(J), C, S,  COORD(J,:), R, M,  &
     &                        D, LOW_LIM, UP_LIM)
          ENDDO
          !$OMP ENDDO
          !$OMP END PARALLEL

          SEMI_NRM = (SEMI_NORM((U - F_SEP_APPRX), SAMPLE_SIZE))**2
                  
          IF (SEMI_NRM .LE. EPS) THEN
            !WRITE(*,*)' '
            !PRINT *, 'Solution converged at rank ', R
            !WRITE(*,*)' '
            FLAG_CONV = .TRUE.
            R_OUT = R
            ERR_DECAY = SEMI_NRM
            !EXIT
            RETURN
          ELSE
            FLAG_CONV = .FALSE.
            !WRITE(*,*)' '
            !WRITE(*,*)'Solution DID NOT CONVERGE'
            !WRITE(*,*)' '
            ! Incrementing R by 1 if not convergered in the last iteration
            R = R + 1
          IF (R .LE. R_MAX) THEN
            !WRITE(*,*)' '
            !WRITE(*,*)'Moving on to the next rank:', R
            !WRITE(*,*)' '
          ENDIF
        ENDIF
          
        ENDDO           ! ENDS THE MAIN DO-WHILE LOOP

        IF (FLAG_CONV .EQV. .FALSE.) THEN
          R_OUT = R_MAX
          ERR_DECAY = SEMI_NRM
        ENDIF

      END SUBROUTINE ALS_BRN

      !==================================================================
      !==================================================================
      

      subroutine timestamp ( )
      !*****************************************************************************80
      !
      !! TIMESTAMP prints the current YMDHMS date as a time stamp.
      !
      !  Example:
      !
      !    31 May 2001   9:45:54.872 AM
      !
      !  Licensing:
      !
      !    This code is distributed under the GNU LGPL license.
      !
      !  Modified:
      !
      !    18 May 2013
      !
      !  Author:
      !
      !    John Burkardt
      !
      !  Parameters:
      !
      !    None
      !
      implicit none

      character ( len = 8 ) ampm
      integer ( kind = 4 ) d
      integer ( kind = 4 ) h
      integer ( kind = 4 ) m
      integer ( kind = 4 ) mm
      character ( len = 9 ), parameter, dimension(12) :: month = (/     &
     &   'January  ', 'February ', 'March    ', 'April    ',            &
     &   'May      ', 'June     ', 'July     ', 'August   ',            &
     &   'September', 'October  ', 'November ', 'December ' /)
      integer ( kind = 4 ) n
      integer ( kind = 4 ) s
      integer ( kind = 4 ) values(8)
      integer ( kind = 4 ) y

      call date_and_time ( values = values )

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if ( h < 12 ) then
        ampm = 'AM'
      else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h < 12 ) then
          ampm = 'PM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
     &   d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim      &
     &   ( ampm )

      return
     end subroutine timestamp

      !==================================================================
      !==================================================================

      !SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      !CHARACTER*(*)    DESC
      !INTEGER          M, N, LDA
      !DOUBLE PRECISION A( LDA, * )
      !INTEGER          I, J
      !WRITE(*,*)
      !WRITE(*,*) DESC
      !DO I = 1, M
      !   WRITE(*,9998) ( A( I, J ), J = 1, N )
      !END DO
! 9998 FORMAT( 11(:,1X,F6.2) )
      !RETURN
      !END
!     Auxiliary routine: printing norms of matrix columns.
      !SUBROUTINE PRINT_VECTOR_NORM( DESC, M, N, A, LDA )
      !CHARACTER*(*)    DESC
      !INTEGER          M, N, LDA
      !DOUBLE PRECISION A( LDA, * )

      !DOUBLE PRECISION TEMP
      !INTEGER          I, J

      !WRITE(*,*)
      !WRITE(*,*) DESC
      !DO J = 1, N
      !   TEMP = 0.0
      !   DO I = 1, M
      !      TEMP = TEMP + A( I, J )*A( I, J )
      !   END DO
      !   WRITE(*,9998,ADVANCE='NO') TEMP
      !END DO
      !WRITE(*,*)
! 9998 FORMAT( 11(:,1X,F6.2) )
      !RETURN
      !END

      !==================================================================
      !==================================================================

      SUBROUTINE GAUSSPOINTS(X, N)
        USE NRTYPE
        IMPLICIT NONE
        ! Interfacing part of the dummy variables
        INTEGER(I4B), INTENT(IN)        :: N
        REAL(DP), INTENT(OUT)           :: X(1:N)
        !  For 1 <= n <= 20, returns the abscissas x of an n point
        !  Gauss-Legendre quadrature rule over the interval [-1,1].

        if ( n == 1 ) then
             x(1) = 0.d0
        elseif ( n == 2 ) then
            x(1) = -0.577350269189625764509148780502
            x(2) =  0.577350269189625764509148780502
        elseif ( n == 3 ) then
            x(1) = -0.774596669241483377035853079956
            x(2) =  0.0000000000000000000000000000000
            x(3) =  0.774596669241483377035853079956
        elseif ( n == 4 ) then
            x(1) = -0.861136311594052575223946488893
            x(2) = -0.339981043584856264802665759103
            x(3) =  0.339981043584856264802665759103
            x(4) =  0.861136311594052575223946488893
        elseif ( n == 5 ) then
            x(1) = -0.906179845938663992797626878299
            x(2) = -0.538469310105683091036314420700
            x(3) =  0.000000000000000000000000000000
            x(4) =  0.538469310105683091036314420700
            x(5) =  0.906179845938663992797626878299
        elseif ( n == 6 ) then
            x(1) = -0.932469514203152027812301554494
            x(2) = -0.661209386466264513661399595020
            x(3) = -0.238619186083196908630501721681
            x(4) =  0.238619186083196908630501721681
            x(5) =  0.661209386466264513661399595020
            x(6) =  0.932469514203152027812301554494
        elseif ( n == 7 ) then
            x(1) = -0.949107912342758524526189684048
            x(2) = -0.741531185599394439863864773281
            x(3) = -0.405845151377397166906606412077
            x(4) =  0.000000000000000000000000000000
            x(5) =  0.405845151377397166906606412077
            x(6) =  0.741531185599394439863864773281
            x(7) =  0.949107912342758524526189684048
        elseif ( n == 8 ) then
            x(1) = -0.960289856497536231683560868569
            x(2) = -0.796666477413626739591553936476
            x(3) = -0.525532409916328985817739049189
            x(4) = -0.183434642495649804939476142360
            x(5) =  0.183434642495649804939476142360
            x(6) =  0.525532409916328985817739049189
            x(7) =  0.796666477413626739591553936476
            x(8) =  0.960289856497536231683560868569
        elseif ( n == 9 ) then
            x(1) = -0.968160239507626089835576202904
            x(2) = -0.836031107326635794299429788070
            x(3) = -0.613371432700590397308702039341
            x(4) = -0.324253423403808929038538014643
            x(5) =  0.000000000000000000000000000000
            x(6) =  0.324253423403808929038538014643
            x(7) =  0.613371432700590397308702039341
            x(8) =  0.836031107326635794299429788070
            x(9) =  0.968160239507626089835576202904
        elseif ( n == 10 ) then
            x(1) =  -0.973906528517171720077964012084
            x(2) =  -0.865063366688984510732096688423
            x(3) =  -0.679409568299024406234327365115
            x(4) =  -0.433395394129247290799265943166
            x(5) =  -0.148874338981631210884826001130
            x(6) =   0.148874338981631210884826001130
            x(7) =   0.433395394129247290799265943166
            x(8) =   0.679409568299024406234327365115
            x(9) =   0.865063366688984510732096688423
            x(10) =  0.973906528517171720077964012084
        elseif ( n == 11 ) then
            x(1)= -0.978228658146056992803938001123
            x(2)= -0.887062599768095299075157769304
            x(3)= -0.730152005574049324093416252031
            x(4)= -0.519096129206811815925725669459
            x(5)= -0.269543155952344972331531985401
            x(6)=  0.000000000000000000000000000000
            x(7)=  0.269543155952344972331531985401
            x(8)=  0.519096129206811815925725669459
            x(9)=  0.730152005574049324093416252031
            x(10)= 0.887062599768095299075157769304
            x(11)= 0.978228658146056992803938001123
        elseif ( n == 12 ) then
            x(1)= -0.981560634246719250690549090149
            x(2)= -0.904117256370474856678465866119
            x(3)= -0.769902674194304687036893833213
            x(4)= -0.587317954286617447296702418941
            x(5)= -0.367831498998180193752691536644
            x(6)= -0.125233408511468915472441369464
            x(7)=  0.125233408511468915472441369464
            x(8)=  0.367831498998180193752691536644
            x(9)=  0.587317954286617447296702418941
            x(10)= 0.769902674194304687036893833213
            x(11)= 0.904117256370474856678465866119
            x(12)= 0.981560634246719250690549090149
        elseif ( n == 13 ) then
            x(1)= -0.984183054718588149472829448807
            x(2)= -0.917598399222977965206547836501
            x(3)= -0.801578090733309912794206489583
            x(4)= -0.642349339440340220643984606996
            x(5)= -0.448492751036446852877912852128
            x(6)= -0.230458315955134794065528121098
            x(7)=  0.000000000000000000000000000000
            x(8)=  0.230458315955134794065528121098
            x(9)=  0.448492751036446852877912852128
            x(10)= 0.642349339440340220643984606996
            x(11)= 0.801578090733309912794206489583
            x(12)= 0.917598399222977965206547836501
            x(13)= 0.984183054718588149472829448807
        elseif ( n == 14 ) then
            x(1)= -0.986283808696812338841597266704
            x(2)= -0.928434883663573517336391139378
            x(3)= -0.827201315069764993189794742650
            x(4)= -0.687292904811685470148019803019
            x(5)= -0.515248636358154091965290718551
            x(6)= -0.319112368927889760435671824168
            x(7)= -0.108054948707343662066244650220
            x(8)=  0.108054948707343662066244650220
            x(9)=  0.319112368927889760435671824168
            x(10)= 0.515248636358154091965290718551
            x(11)= 0.687292904811685470148019803019
            x(12)= 0.827201315069764993189794742650
            x(13)= 0.928434883663573517336391139378
            x(14)= 0.986283808696812338841597266704
        elseif ( n == 15 ) then
            x(1)= -0.987992518020485428489565718587
            x(2)= -0.937273392400705904307758947710
            x(3)= -0.848206583410427216200648320774
            x(4)= -0.724417731360170047416186054614
            x(5)= -0.570972172608538847537226737254
            x(6)= -0.394151347077563369897207370981
            x(7)= -0.201194093997434522300628303395
            x(8)=  0.000000000000000000000000000000
            x(9)=  0.201194093997434522300628303395
            x(10)= 0.394151347077563369897207370981
            x(11)= 0.570972172608538847537226737254
            x(12)= 0.724417731360170047416186054614
            x(13)= 0.848206583410427216200648320774
            x(14)= 0.937273392400705904307758947710
            x(15)= 0.987992518020485428489565718587
        elseif (n == 16) then
            x(1)= -0.989400934991649932596154173450
            x(2)= -0.944575023073232576077988415535
            x(3)= -0.865631202387831743880467897712
            x(4)= -0.755404408355003033895101194847
            x(5)= -0.617876244402643748446671764049
            x(6)= -0.458016777657227386342419442984
            x(7)= -0.281603550779258913230460501460
            x(8)= -0.0950125098376374401853193354250
            x(9)=  0.0950125098376374401853193354250
            x(10)= 0.281603550779258913230460501460
            x(11)= 0.458016777657227386342419442984
            x(12)= 0.617876244402643748446671764049
            x(13)= 0.755404408355003033895101194847
            x(14)= 0.865631202387831743880467897712
            x(15)= 0.944575023073232576077988415535
            x(16)= 0.989400934991649932596154173450
        elseif (n == 17) then
            x(1)= -0.990575475314417335675434019941
            x(2)= -0.950675521768767761222716957896
            x(3)= -0.880239153726985902122955694488
            x(4)= -0.781514003896801406925230055520
            x(5)= -0.657671159216690765850302216643
            x(6)= -0.512690537086476967886246568630
            x(7)= -0.351231763453876315297185517095
            x(8)= -0.178484181495847855850677493654
            x(9)=  0.000000000000000000000000000000
            x(10)= 0.178484181495847855850677493654
            x(11)= 0.351231763453876315297185517095
            x(12)= 0.512690537086476967886246568630
            x(13)= 0.657671159216690765850302216643
            x(14)= 0.781514003896801406925230055520
            x(15)= 0.880239153726985902122955694488
            x(16)= 0.950675521768767761222716957896
            x(17)= 0.990575475314417335675434019941
        elseif (n == 18) then
            x(1)= -0.991565168420930946730016004706
            x(2)= -0.955823949571397755181195892930
            x(3)= -0.892602466497555739206060591127
            x(4)= -0.803704958972523115682417455015
            x(5)= -0.691687043060353207874891081289
            x(6)= -0.559770831073947534607871548525
            x(7)= -0.411751161462842646035931793833
            x(8)= -0.251886225691505509588972854878
            x(9)= -0.0847750130417353012422618529358
            x(10)= 0.084775013041735301242261852936
            x(11)= 0.251886225691505509588972854878
            x(12)= 0.411751161462842646035931793833
            x(13)= 0.559770831073947534607871548525
            x(14)= 0.691687043060353207874891081289
            x(15)= 0.803704958972523115682417455015
            x(16)= 0.892602466497555739206060591127
            x(17)= 0.955823949571397755181195892930
            x(18)= 0.991565168420930946730016004706
        elseif (n == 19) then
            x(1)= -0.992406843843584403189017670253
            x(2)= -0.960208152134830030852778840688
            x(3)= -0.903155903614817901642660928532
            x(4)= -0.822714656537142824978922486713
            x(5)= -0.720966177335229378617095860824
            x(6)= -0.600545304661681023469638164946
            x(7)= -0.464570741375960945717267148104
            x(8)= -0.316564099963629831990117328850
            x(9)= -0.160358645640225375868096115741
            x(10)= 0.00000000000000000000000000000
            x(11)= 0.160358645640225375868096115741
            x(12)= 0.316564099963629831990117328850
            x(13)= 0.464570741375960945717267148104
            x(14)= 0.600545304661681023469638164946
            x(15)= 0.720966177335229378617095860824
            x(16)= 0.822714656537142824978922486713
            x(17)= 0.903155903614817901642660928532
            x(18)= 0.960208152134830030852778840688
            x(19)= 0.992406843843584403189017670253
        elseif (n == 20) then
            x(1)= -0.993128599185094924786122388471
            x(2)= -0.963971927277913791267666131197
            x(3)= -0.912234428251325905867752441203
            x(4)= -0.839116971822218823394529061702
            x(5)= -0.746331906460150792614305070356
            x(6)= -0.636053680726515025452836696226
            x(7)= -0.510867001950827098004364050955
            x(8)= -0.373706088715419560672548177025
            x(9)= -0.227785851141645078080496195369
            x(10)= -0.0765265211334973337546404093988
            x(11)= 0.0765265211334973337546404093988
            x(12)= 0.227785851141645078080496195369
            x(13)= 0.373706088715419560672548177025
            x(14)= 0.510867001950827098004364050955
            x(15)= 0.636053680726515025452836696226
            x(16)= 0.746331906460150792614305070356
            x(17)= 0.839116971822218823394529061702
            x(18)= 0.912234428251325905867752441203
            x(19)= 0.963971927277913791267666131197
            x(20)= 0.993128599185094924786122388471
        else
          STOP 'GAUSSPOINTS: Fatal error! Illegal value of n.'
        endif

      END SUBROUTINE GAUSSPOINTS

      !==================================================================
      !==================================================================

      SUBROUTINE GAUSSWEIGHTS(W, N)
        USE NRTYPE
        IMPLICIT NONE
        ! Interfacing part of the dummy variables
        INTEGER(I4B), INTENT(IN)        :: N
        REAL(DP), INTENT(OUT)           :: W(1:N)
        !  For 1 <= N <= 20, returns the weights W of an
        !  N point Gauss-Legendre quadrature rule over the interval [-1,1].

        if ( n == 1 ) then
            w(1) = 2.d0
        elseif ( n == 2 ) then
            w(1) =  1.d0 
            w(2) =  w(1)
        elseif ( n == 3 ) then
            w(1) =  0.5555555555555555555555555555565
            w(2) =  0.8888888888888888888888888888889
            w(3) =  0.5555555555555555555555555555565
        elseif ( n == 4 ) then
            w(1) = 0.347854845137453857373063949222
            w(2) = 0.652145154862546142626936050778
            w(3) = 0.652145154862546142626936050778
            w(4) = 0.347854845137453857373063949222
        elseif ( n == 5 ) then
            w(1) = 0.236926885056189087514264040720
            w(2) = 0.478628670499366468041291514836
            w(3) = 0.568888888888888888888888888889
            w(4) = 0.478628670499366468041291514836
            w(5) = 0.236926885056189087514264040720
        elseif ( n == 6 ) then
            w(1) = 0.171324492379170345040296142173
            w(2) = 0.360761573048138607569833513838
            w(3) = 0.467913934572691047389870343990
            w(4) = 0.467913934572691047389870343990
            w(5) = 0.360761573048138607569833513838
            w(6) = 0.171324492379170345040296142173
        elseif ( n == 7 ) then
            w(1) = 0.129484966168869693270611432679
            w(2) = 0.279705391489276667901467771424
            w(3) = 0.381830050505118944950369775489
            w(4) = 0.417959183673469387755102040816
            w(5) = 0.381830050505118944950369775489
            w(6) = 0.279705391489276667901467771424
            w(7) = 0.129484966168869693270611432679
        elseif ( n == 8 ) then
            w(1) = 0.101228536290376259152531354310
            w(2) = 0.222381034453374470544355994426
            w(3) = 0.313706645877887287337962201987
            w(4) = 0.362683783378361982965150449277
            w(5) = 0.362683783378361982965150449277
            w(6) = 0.313706645877887287337962201987
            w(7) = 0.222381034453374470544355994426
            w(8) = 0.101228536290376259152531354310
        elseif ( n == 9 ) then
            w(1) = 0.0812743883615744119718921581105
            w(2) = 0.180648160694857404058472031243
            w(3) = 0.260610696402935462318742869419
            w(4) = 0.312347077040002840068630406584
            w(5) = 0.330239355001259763164525069287
            w(6) = 0.312347077040002840068630406584
            w(7) = 0.260610696402935462318742869419
            w(8) = 0.180648160694857404058472031243
            w(9) = 0.0812743883615744119718921581105
        elseif ( n == 10 ) then
            w(1) =  0.0666713443086881375935688098933
            w(2) =  0.149451349150580593145776339658
            w(3) =  0.219086362515982043995534934228
            w(4) =  0.269266719309996355091226921569
            w(5) =  0.295524224714752870173892994651
            w(6) =  0.295524224714752870173892994651
            w(7) =  0.269266719309996355091226921569
            w(8) =  0.219086362515982043995534934228
            w(9) =  0.149451349150580593145776339658
            w(10) = 0.0666713443086881375935688098933
        elseif ( n == 11 ) then
            w(1)= 0.0556685671161736664827537204425
            w(2)= 0.125580369464904624634694299224
            w(3)= 0.186290210927734251426097641432
            w(4)= 0.233193764591990479918523704843
            w(5)= 0.262804544510246662180688869891
            w(6)= 0.272925086777900630714483528336
            w(7)= 0.262804544510246662180688869891
            w(8)= 0.233193764591990479918523704843
            w(9)= 0.186290210927734251426097641432
            w(10)= 0.125580369464904624634694299224
            w(11)= 0.0556685671161736664827537204425
        elseif ( n == 12 ) then
            w(1)= 0.0471753363865118271946159614850
            w(2)= 0.106939325995318430960254718194
            w(3)= 0.160078328543346226334652529543
            w(4)= 0.203167426723065921749064455810
            w(5)= 0.233492536538354808760849898925
            w(6)= 0.249147045813402785000562436043
            w(7)= 0.249147045813402785000562436043
            w(8)= 0.233492536538354808760849898925
            w(9)= 0.203167426723065921749064455810
            w(10)= 0.160078328543346226334652529543
            w(11)= 0.106939325995318430960254718194
            w(12)= 0.0471753363865118271946159614850
        elseif ( n == 13 ) then
            w(1)= 0.0404840047653158795200215922010
            w(2)= 0.0921214998377284479144217759538
            w(3)= 0.138873510219787238463601776869
            w(4)= 0.178145980761945738280046691996
            w(5)= 0.207816047536888502312523219306
            w(6)= 0.226283180262897238412090186040
            w(7)= 0.232551553230873910194589515269
            w(8)= 0.226283180262897238412090186040
            w(9)= 0.207816047536888502312523219306
            w(10)= 0.178145980761945738280046691996
            w(11)= 0.138873510219787238463601776869
            w(12)= 0.0921214998377284479144217759538
            w(13)= 0.0404840047653158795200215922010
        elseif ( n == 14 ) then
            w(1)= 0.0351194603317518630318328761382
            w(2)= 0.0801580871597602098056332770629
            w(3)= 0.121518570687903184689414809072
            w(4)= 0.157203167158193534569601938624
            w(5)= 0.185538397477937813741716590125
            w(6)= 0.205198463721295603965924065661
            w(7)= 0.215263853463157790195876443316
            w(8)= 0.215263853463157790195876443316
            w(9)= 0.205198463721295603965924065661
            w(10)= 0.185538397477937813741716590125
            w(11)= 0.157203167158193534569601938624
            w(12)= 0.121518570687903184689414809072
            w(13)= 0.0801580871597602098056332770629
            w(14)= 0.0351194603317518630318328761382
        elseif ( n == 15 ) then
            w(1)= 0.0307532419961172683546283935772
            w(2)= 0.0703660474881081247092674164507
            w(3)= 0.107159220467171935011869546686
            w(4)= 0.139570677926154314447804794511
            w(5)= 0.166269205816993933553200860481
            w(6)= 0.186161000015562211026800561866
            w(7)= 0.198431485327111576456118326444
            w(8)= 0.202578241925561272880620199968
            w(9)= 0.198431485327111576456118326444
            w(10)= 0.186161000015562211026800561866
            w(11)= 0.166269205816993933553200860481
            w(12)= 0.139570677926154314447804794511
            w(13)= 0.107159220467171935011869546686
            w(14)= 0.0703660474881081247092674164507
            w(15)= 0.0307532419961172683546283935772
        elseif ( n == 16 ) then
            w(1)= 0.0271524594117540948517805724560
            w(2)= 0.0622535239386478928628438369944
            w(3)= 0.0951585116824927848099251076022
            w(4)= 0.124628971255533872052476282192
            w(5)= 0.149595988816576732081501730547
            w(6)= 0.169156519395002538189312079030
            w(7)= 0.182603415044923588866763667969
            w(8)= 0.189450610455068496285396723208
            w(9)= 0.189450610455068496285396723208
            w(10)= 0.182603415044923588866763667969
            w(11)= 0.169156519395002538189312079030
            w(12)= 0.149595988816576732081501730547
            w(13)= 0.124628971255533872052476282192
            w(14)= 0.0951585116824927848099251076022
            w(15)= 0.0622535239386478928628438369944
            w(16)= 0.0271524594117540948517805724560
        elseif ( n == 17 ) then
            w(1)= 0.0241483028685479319601100262876
            w(2)= 0.0554595293739872011294401653582
            w(3)= 0.0850361483171791808835353701911
            w(4)= 0.111883847193403971094788385626
            w(5)= 0.135136368468525473286319981702
            w(6)= 0.154045761076810288081431594802
            w(7)= 0.168004102156450044509970663788
            w(8)= 0.176562705366992646325270990113
            w(9)= 0.179446470356206525458265644262
            w(10)= 0.176562705366992646325270990113
            w(11)= 0.168004102156450044509970663788
            w(12)= 0.154045761076810288081431594802
            w(13)= 0.135136368468525473286319981702
            w(14)= 0.111883847193403971094788385626
            w(15)= 0.0850361483171791808835353701911
            w(16)= 0.0554595293739872011294401653582
            w(17)= 0.0241483028685479319601100262876
        elseif ( n == 18 ) then
            w(1)= 0.0216160135264833103133427102665
            w(2)= 0.0497145488949697964533349462026
            w(3)= 0.0764257302548890565291296776166
            w(4)= 0.100942044106287165562813984925
            w(5)= 0.122555206711478460184519126800
            w(6)= 0.140642914670650651204731303752
            w(7)= 0.154684675126265244925418003836
            w(8)= 0.164276483745832722986053776466
            w(9)= 0.169142382963143591840656470135
            w(10)= 0.169142382963143591840656470135
            w(11)= 0.164276483745832722986053776466
            w(12)= 0.154684675126265244925418003836
            w(13)= 0.140642914670650651204731303752
            w(14)= 0.122555206711478460184519126800
            w(15)= 0.100942044106287165562813984925
            w(16)= 0.0764257302548890565291296776166
            w(17)= 0.0497145488949697964533349462026
            w(18)= 0.0216160135264833103133427102665
        elseif ( n == 19 ) then
            w(1)= 0.0194617882297264770363120414644
            w(2)= 0.0448142267656996003328381574020
            w(3)= 0.0690445427376412265807082580060
            w(4)= 0.0914900216224499994644620941238
            w(5)= 0.111566645547333994716023901682
            w(6)= 0.128753962539336227675515784857
            w(7)= 0.142606702173606611775746109442
            w(8)= 0.152766042065859666778855400898
            w(9)= 0.158968843393954347649956439465
            w(10)= 0.161054449848783695979163625321
            w(11)= 0.158968843393954347649956439465
            w(12)= 0.152766042065859666778855400898
            w(13)= 0.142606702173606611775746109442
            w(14)= 0.128753962539336227675515784857
            w(15)= 0.111566645547333994716023901682
            w(16)= 0.0914900216224499994644620941238
            w(17)= 0.0690445427376412265807082580060
            w(18)= 0.0448142267656996003328381574020
            w(19)= 0.0194617882297264770363120414644
        elseif ( n == 20 ) then
            w(1)= 0.0176140071391521183118619623519
            w(2)= 0.0406014298003869413310399522749
            w(3)= 0.0626720483341090635695065351870
            w(4)= 0.0832767415767047487247581432220
            w(5)= 0.101930119817240435036750135480
            w(6)= 0.118194531961518417312377377711
            w(7)= 0.131688638449176626898494499748
            w(8)= 0.142096109318382051329298325067
            w(9)= 0.149172986472603746787828737002
            w(10)= 0.152753387130725850698084331955
            w(11)= 0.152753387130725850698084331955
            w(12)= 0.149172986472603746787828737002
            w(13)= 0.142096109318382051329298325067
            w(14)= 0.131688638449176626898494499748
            w(15)= 0.118194531961518417312377377711
            w(16)= 0.101930119817240435036750135480
            w(17)= 0.0832767415767047487247581432220
            w(18)= 0.0626720483341090635695065351870
            w(19)= 0.0406014298003869413310399522749
            w(20)= 0.0176140071391521183118619623519
        else
            STOP 'GAUSSWEIGHTS: Fatal error! Illegal value of n.'
        endif

      END SUBROUTINE GAUSSWEIGHTS
subroutine imtqlx ( n, d, e, z )

!*****************************************************************************80
!
!! IMTQLX diagonalizes a symmetric tridiagonal matrix.
!
!  Discussion:
!
!    This routine is a slightly modified version of the EISPACK routine to 
!    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
!
!    The authors thank the authors of EISPACK for permission to use this
!    routine. 
!
!    It has been modified to produce the product Q' * Z, where Z is an input 
!    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
!    The changes consist (essentially) of applying the orthogonal 
!    transformations directly to Z as they are generated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    27 December 2009
!
!  Author:
!
!    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Sylvan Elhay, Jaroslav Kautsky,
!    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
!    Interpolatory Quadrature,
!    ACM Transactions on Mathematical Software,
!    Volume 13, Number 4, December 1987, pages 399-415.
!
!    Roger Martin, James Wilkinson,
!    The Implicit QL Algorithm,
!    Numerische Mathematik,
!    Volume 12, Number 5, December 1968, pages 377-383.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) D(N), the diagonal entries of the matrix.
!    On output, the information in D has been overwritten.
!
!    Input/output, real ( kind = 8 ) E(N), the subdiagonal entries of the 
!    matrix, in entries E(1) through E(N-1).  On output, the information in
!    E has been overwritten.
!
!    Input/output, real ( kind = 8 ) Z(N).  On input, a vector.  On output,
!    the value of Q' * Z, where Q is the matrix that diagonalizes the
!    input symmetric tridiagonal matrix.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ), parameter :: itn = 30
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mml
  real ( kind = 8 ) p
  real ( kind = 8 ) prec
  real ( kind = 8 ) r
  real ( kind = 8 ) s
  real ( kind = 8 ) z(n)

  prec = epsilon ( prec )

  if ( n == 1 ) then
    return
  end if

  e(n) = 0.0D+00

  do l = 1, n

    j = 0

    do

      do m = l, n

        if ( m == n ) then
          exit
        end if

        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) ) then
          exit
        end if

      end do

      p = d(l)

      if ( m == l ) then
        exit
      end if

      if ( itn <= j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IMTQLX - Fatal error!'
        write ( *, '(a)' ) '  Iteration limit exceeded.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  L = ', l
        write ( *, '(a,i8)' ) '  M = ', m
        write ( *, '(a,i8)' ) '  N = ', n
        stop 1
      end if

      j = j + 1
      g = ( d(l+1) - p ) / ( 2.0D+00 * e(l) )
      r =  sqrt ( g * g + 1.0D+00 )
      g = d(m) - p + e(l) / ( g + sign ( r, g ) )
      s = 1.0D+00
      c = 1.0D+00
      p = 0.0D+00
      mml = m - l

      do ii = 1, mml

        i = m - ii
        f = s * e(i)
        b = c * e(i)

        if ( abs ( g ) <= abs ( f ) ) then
          c = g / f
          r =  sqrt ( c * c + 1.0D+00 )
          e(i+1) = f * r
          s = 1.0D+00 / r
          c = c * s
        else
          s = f / g
          r =  sqrt ( s * s + 1.0D+00 )
          e(i+1) = g * r
          c = 1.0D+00 / r
          s = s * c
        end if

        g = d(i+1) - p
        r = ( d(i) - g ) * s + 2.0D+00 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        f = z(i+1)
        z(i+1) = s * z(i) + c * f
        z(i) = c * z(i) - s * f

      end do

      d(l) = d(l) - p
      e(l) = g
      e(m) = 0.0D+00

    end do

  end do
!
!  Sorting.
!
  do ii = 2, n

    i = ii - 1
    k = i
    p = d(i)

    do j = ii, n
      if ( d(j) < p ) then
        k = j
        p = d(j)
      end if
    end do

    if ( k /= i ) then
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end if

  end do

  return
end
subroutine p_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! P_EXPONENTIAL_PRODUCT: exponential products for P(n,x).
!
!  Discussion:
!
!    Let P(n,x) represent the Legendre polynomial of degree n.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -1.0 <= X <= +1.0 ) exp(B*X) * P(I,X) * P(J,X) dx
!
!    We will estimate these integrals using Gauss-Legendre quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x(1)
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call p_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x(1) = x_table(k)
    call p_polynomial_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x(1) ) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine p_integral ( n, value )

!*****************************************************************************80
!
!! P_INTEGRAL evaluates a monomial integral associated with P(n,x).
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x < +1 ) x^n dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the exponent.
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) value

  if ( mod ( n, 2 ) == 1 ) then
    value = 0.0D+00
  else
    value = 2.0D+00 / real ( n + 1, kind = 8 )
  end if

  return
end
subroutine p_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
!
!  Discussion:
!
!     1
!     0     1
!    -1/2   0      3/2
!     0    -3/2    0     5/2
!     3/8   0    -30/8   0     35/8
!     0    15/8    0   -70/8    0     63/8
!    -5/16  0    105/16  0   -315/16   0    231/16
!     0   -35/16   0   315/16   0   -693/16   0    429/16
!
!     1.00000
!     0.00000  1.00000
!    -0.50000  0.00000  1.50000
!     0.00000 -1.50000  0.00000  2.5000
!     0.37500  0.00000 -3.75000  0.00000  4.37500
!     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
!    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
!     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 February 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the 
!    Legendre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n <= 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0:i-2) =          real (   - i + 1, kind = 8 ) * c(i-2,0:i-2) &
                        / real (     i,     kind = 8 )
    c(i,1:i) = c(i,1:i) + real ( i + i - 1, kind = 8 ) * c(i-1,0:i-1) &
                        / real (     i,     kind = 8 )
  end do
 
  return
end
subroutine p_polynomial_prime ( m, n, x, vp )

!*****************************************************************************80
!
!! P_POLYNOMIAL_PRIME evaluates the derivative of Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) VP(M,0:N), the values of the derivatives of the
!    Legendre polynomials of order 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) vp(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00
  vp(1:m,0) = 0.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
  vp(1:m,1) = 1.0D+00
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
    vp(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * ( v(1:m,i-1) &
                                                   + x(1:m) * vp(1:m,i-1) ) &
                - real (     i - 1, kind = 8 ) *   vp(1:m,i-2)               ) &
                / real (     i,     kind = 8 )
 
  end do
 
  return
end
subroutine p_polynomial_prime2 ( m, n, x, vpp )

!*****************************************************************************80
!
!! P_POLYNOMIAL_PRIME2: second derivative of Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(0,X) = 1
!    P(1,X) = X
!    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
!
!    P'(0,X) = 0
!    P'(1,X) = 1
!    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
!
!    P"(0,X) = 0
!    P"(1,X) = 0
!    P"(N,X) = ( (2*N-1)*(2*P(N-1,X)+X*P"(N-1,X)-(N-1)*P"(N-2,X) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) VPP(M,0:N), the second derivative of the
!    Legendre polynomials of order 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) vp(m,0:n)
  real ( kind = 8 ) vpp(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00
  vp(1:m,0) = 0.0D+00
  vpp(1:m,0) = 0.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
  vp(1:m,1) = 1.0D+00
  vpp(1:m,1) = 0.0D+00
 
  do i = 2, n
 
    v(1:m,i) = &
      ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
      - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
      / real (     i,     kind = 8 )
 
    vp(1:m,i) = &
      ( real ( 2 * i - 1, kind = 8 ) * ( v(1:m,i-1) + x(1:m) * vp(1:m,i-1) ) &
      - real (     i - 1, kind = 8 ) *   vp(1:m,i-2)               ) &
      / real (     i,     kind = 8 )

    vpp(1:m,i) = &
      ( real ( 2 * i - 1, kind = 8 ) * ( 2.0D+00 * vp(1:m,i-1) &
                                         + x(1:m) * vpp(1:m,i-1) ) &
      - real (     i - 1, kind = 8 ) *   vpp(1:m,i-2)               ) &
      / real (     i,     kind = 8 )

  end do
 
  return
end


subroutine p_polynomial_value ( m, n, x, v )

!*****************************************************************************80
!
!! P_POLYNOMIAL_VALUE evaluates the Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(n,1) = 1.
!    P(n,-1) = (-1)^N.
!    | P(n,x) | <= 1 in [-1,1].
!
!    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
!
!    The formula is:
!
!      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,x) =      1
!    P( 1,x) =      1 X
!    P( 2,x) = (    3 X^2 -       1)/2
!    P( 3,x) = (    5 X^3 -     3 X)/2
!    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
!
!  Recursion:
!
!    P(0,x) = 1
!    P(1,x) = x
!    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
!
!    P'(0,x) = 0
!    P'(1,x) = 1
!    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
  end do
 
  return
end
subroutine p_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! P_POLYNOMIAL_VALUES returns values of the Legendre polynomials P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the function.
!
!    Output, real ( kind = 8 ) X, the point where the function is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 22

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1000000000000000D+01, &
     0.2500000000000000D+00, &
    -0.4062500000000000D+00, &
    -0.3359375000000000D+00, &
     0.1577148437500000D+00, &
     0.3397216796875000D+00, &
     0.2427673339843750D-01, &
    -0.2799186706542969D+00, &
    -0.1524540185928345D+00, &
     0.1768244206905365D+00, &
     0.2212002165615559D+00, &
     0.0000000000000000D+00, &
    -0.1475000000000000D+00, &
    -0.2800000000000000D+00, &
    -0.3825000000000000D+00, &
    -0.4400000000000000D+00, &
    -0.4375000000000000D+00, &
    -0.3600000000000000D+00, &
    -0.1925000000000000D+00, &
     0.8000000000000000D-01, &
     0.4725000000000000D+00, &
     0.1000000000000000D+01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3,  3,  3, &
     3 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.00D+00, &
    0.10D+00, &
    0.20D+00, &
    0.30D+00, &
    0.40D+00, &
    0.50D+00, &
    0.60D+00, &
    0.70D+00, &
    0.80D+00, &
    0.90D+00, &
    1.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine p_polynomial_zeros ( nt, t )

!*****************************************************************************80
!
!! P_POLYNOMIAL_ZEROS: zeros of Legendre function P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), the zeros.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  
  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:nt) = sqrt ( bj(1:nt) )

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( 2.0D+00 )

  call imtqlx ( nt, t, bj, wts )

  return
end
subroutine p_power_product ( p, e, table )

!*****************************************************************************80
!
!! P_POWER_PRODUCT: power products for Legendre polynomial P(n,x).
!
!  Discussion:
!
!    Let P(n,x) represent the Legendre polynomial of degree n.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -1.0 <= X <= +1.0 ) X^E * P(i,x) * P(j,x) dx
!
!    We will estimate these integrals using Gauss-Legendre quadrature.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, integer ( kind = 4 ) E, the exponent of X in the integrand.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) e
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x(1)
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( ( e + 1 ) / 2 )

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call p_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x(1) = x_table(k)
    call p_polynomial_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    if ( e == 0 ) then
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
        end do
      end do
    else
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * x(1) ** e * h_table(i) * h_table(j)
        end do
      end do
    end if

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine p_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! P_QUADRATURE_RULE: quadrature for Legendre function P(n,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), WTS(NT), the points and weights
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)
  
  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = real ( i * i, kind = 8 ) / real ( 4 * i * i - 1, kind = 8 )
  end do
  bj(1:nt) = sqrt ( bj(1:nt) )

  wts(1) = sqrt ( 2.0D+00 )
  wts(2:nt) = 0.0D+00

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2

  return
end
subroutine pm_polynomial_value ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).
!
!  Differential equation:
!
!    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0
!
!  First terms:
!
!    M = 0  ( = Legendre polynomials of first kind P(N,X) )
!
!    Pm(0,0,x) =    1
!    Pm(1,0,x) =    1 X
!    Pm(2,0,x) = (  3 X^2 -   1)/2
!    Pm(3,0,x) = (  5 X^3 -   3 X)/2
!    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
!    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
!    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
!    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
!
!    M = 1
!
!    Pm(0,1,x) =   0
!    Pm(1,1,x) =   1 * SQRT(1-X^2)
!    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
!    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
!    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)
!
!    M = 2
!
!    Pm(0,2,x) =   0
!    Pm(1,2,x) =   0
!    Pm(2,2,x) =   3 * (1-X^2)
!    Pm(3,2,x) =  15 * (1-X^2) * X
!    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)
!
!    M = 3
!
!    Pm(0,3,x) =   0
!    Pm(1,3,x) =   0
!    Pm(2,3,x) =   0
!    Pm(3,3,x) =  15 * (1-X^2)^1.5
!    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X
!
!    M = 4
!
!    Pm(0,4,x) =   0
!    Pm(1,4,x) =   0
!    Pm(2,4,x) =   0
!    Pm(3,4,x) =   0
!    Pm(4,4,x) = 105 * (1-X^2)^2
!
!  Recursion:
!
!    if N < M:
!      Pm(N,M,x) = 0
!    if N = M:
!      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
!      all the odd integers less than or equal to N.
!    if N = M+1:
!      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
!    if M+1 < N:
!      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 May 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the point at which the function is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(mm)

  cx(1:mm,0:n) = 0.0D+00
!
!  J = M is the first nonzero function.
!
  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00

    fact = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * fact * sqrt ( 1.0D+00 - x(1:mm)**2 )
      fact = fact + 2.0D+00
    end do

  end if
!
!  J = M + 1 is the second nonzero function.
!
  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if
!
!  Now we use a three term recurrence.
!
  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
                 / real (     j - m,     kind = 8 )
  end do

  return
end
subroutine pm_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! PM_POLYNOMIAL_VALUES returns values of Legendre polynomials Pm(n,m,x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.0000000000000000D+00, &
    -0.5000000000000000D+00, &
     0.0000000000000000D+00, &
     0.3750000000000000D+00, &
     0.0000000000000000D+00, &
    -0.8660254037844386D+00, &
    -0.1299038105676658D+01, &
    -0.3247595264191645D+00, &
     0.1353164693413185D+01, &
    -0.2800000000000000D+00, &
     0.1175755076535925D+01, &
     0.2880000000000000D+01, &
    -0.1410906091843111D+02, &
    -0.3955078125000000D+01, &
    -0.9997558593750000D+01, &
     0.8265311444100484D+02, &
     0.2024442836815152D+02, &
    -0.4237997531890869D+03, &
     0.1638320624828339D+04, &
    -0.2025687389227225D+05 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 0, 0, &
    0, 1, 1, 1, &
    1, 0, 1, 2, &
    3, 2, 2, 3, &
    3, 4, 4, 5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    1,  2,  3,  4, &
    5,  1,  2,  3, &
    4,  3,  3,  3, &
    3,  4,  5,  6, &
    7,  8,  9, 10 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.00D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.20D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00, &
    0.25D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine pmn_polynomial_value ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PMN_POLYNOMIAL_VALUE: normalized Legendre polynomial Pmn(n,m,x).
!
!  Discussion:
!
!    The unnormalized associated Legendre functions P_N^M(X) have
!    the property that
!
!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
!
!    By dividing the function by the square root of this term,
!    the normalized associated Legendre functions have norm 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 March 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) factor
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) x(mm)

  if ( m < 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M is ', m
    write ( *, '(a)' ) '  but M must be nonnegative.'
    stop 1
  end if
 
  if ( n < m ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PMN_POLYNOMIAL_VALUE - Fatal error!'
    write ( *, '(a,i8)' ) '  Input value of M = ', m
    write ( *, '(a,i8)' ) '  Input value of N = ', n
    write ( *, '(a)' ) '  but M must be less than or equal to N.'
    stop 1
  end if

  cx(1:mm,0:n) = 0.0D+00

  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00
    factor = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0D+00 - x(1:mm)**2 )
      factor = factor + 2.0D+00
    end do
  end if

  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if

  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
                 / real (     j - m,     kind = 8 )
  end do
!
!  Normalization.
!
  do j = m, n
    factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * r8_factorial ( j - m ) ) &
      / ( 2.0D+00 * r8_factorial ( j + m ) ) )
    cx(1:mm,j) = cx(1:mm,j) * factor
  end do

  return
end
subroutine pmn_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! PMN_POLYNOMIAL_VALUES: normalized Legendre polynomial Pmn(n,m,x).
!
!  Discussion:
!
!    In Mathematica, the unnormalized function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by the factor:
!
!      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7071067811865475D+00, &
    0.6123724356957945D+00, &
   -0.7500000000000000D+00, &
   -0.1976423537605237D+00, &
   -0.8385254915624211D+00, &
    0.7261843774138907D+00, &
   -0.8184875533567997D+00, &
   -0.1753901900050285D+00, &
    0.9606516343087123D+00, &
   -0.6792832849776299D+00, &
   -0.6131941618102092D+00, &
    0.6418623720763665D+00, &
    0.4716705890038619D+00, &
   -0.1018924927466445D+01, &
    0.6239615396237876D+00, &
    0.2107022704608181D+00, &
    0.8256314721961969D+00, &
   -0.3982651281554632D+00, &
   -0.7040399320721435D+00, &
    0.1034723155272289D+01, &
   -0.5667412129155530D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine pmns_polynomial_value ( mm, n, m, x, cx )

!*****************************************************************************80
!
!! PMNS_POLYNOMIAL_VALUE: sphere-normalized Legendre polynomial Pmns(n,m,x).
!
!  Discussion:
!
!    The unnormalized associated Legendre functions P_N^M(X) have
!    the property that
!
!      Integral ( -1 <= X <= 1 ) ( P_N^M(X) )^2 dX 
!      = 2 * ( N + M )! / ( ( 2 * N + 1 ) * ( N - M )! )
!
!    By dividing the function by the square root of this term,
!    the normalized associated Legendre functions have norm 1.
!
!    However, we plan to use these functions to build spherical
!    harmonics, so we use a slightly different normalization factor of
!
!      sqrt ( ( ( 2 * N + 1 ) * ( N - M )! ) / ( 4 * pi * ( N + M )! ) ) 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MM, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the maximum first index of the Legendre
!    function, which must be at least 0.
!
!    Input, integer ( kind = 4 ) M, the second index of the Legendre function,
!    which must be at least 0, and no greater than N.
!
!    Input, real ( kind = 8 ) X(MM), the evaluation points.
!
!    Output, real ( kind = 8 ) CX(MM,0:N), the function values.
!
  implicit none

  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n

  real ( kind = 8 ) cx(mm,0:n)
  real ( kind = 8 ) factor
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x(mm)

  cx(1:mm,0:n) = 0.0D+00

  if ( m <= n ) then
    cx(1:mm,m) = 1.0D+00
    factor = 1.0D+00
    do j = 1, m
      cx(1:mm,m) = - cx(1:mm,m) * factor * sqrt ( 1.0D+00 - x(1:mm) ** 2 )
      factor = factor + 2.0D+00
    end do
  end if

  if ( m + 1 <= n ) then
    cx(1:mm,m+1) = x(1:mm) * real ( 2 * m + 1, kind = 8 ) * cx(1:mm,m)
  end if

  do j = m + 2, n
    cx(1:mm,j) = ( real ( 2 * j     - 1, kind = 8 ) * x(1:mm) * cx(1:mm,j-1) &
                 + real (   - j - m + 1, kind = 8 ) *           cx(1:mm,j-2) ) &
                 / real (     j - m,     kind = 8 )
  end do
!
!  Normalization.
!
  do j = m, n
    factor = sqrt ( ( real ( 2 * j + 1, kind = 8 ) * r8_factorial ( j - m ) ) &
      / ( 4.0D+00 * r8_pi * r8_factorial ( j + m ) ) )
    cx(1:mm,j) = cx(1:mm,j) * factor
  end do

  return
end
subroutine pmns_polynomial_values ( n_data, n, m, x, fx )

!*****************************************************************************80
!
!! PMNS_POLYNOMIAL_VALUES: sphere-normalized Legendre polynomial Pmns(n,m,x).
!
!  Discussion:
!
!    In Mathematica, the unnormalized function can be evaluated by:
!
!      LegendreP [ n, m, x ]
!
!    The function is normalized by dividing by the factor:
!
!      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, integer ( kind = 4 ) M, 
!    real ( kind = 8 ) X, the arguments of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.2820947917738781D+00, &
     0.2443012559514600D+00, &
    -0.2992067103010745D+00, &
    -0.07884789131313000D+00, &
    -0.3345232717786446D+00, &
     0.2897056515173922D+00, &
    -0.3265292910163510D+00, &
    -0.06997056236064664D+00, &
     0.3832445536624809D+00, &
    -0.2709948227475519D+00, &
    -0.2446290772414100D+00, &
     0.2560660384200185D+00, &
     0.1881693403754876D+00, &
    -0.4064922341213279D+00, &
     0.2489246395003027D+00, &
     0.08405804426339821D+00, &
     0.3293793022891428D+00, &
    -0.1588847984307093D+00, &
    -0.2808712959945307D+00, &
     0.4127948151484925D+00, &
    -0.2260970318780046D+00 /)
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save, dimension ( n_max ) :: m_vec = (/ &
    0, 0, 1, 0, &
    1, 2, 0, 1, &
    2, 3, 0, 1, &
    2, 3, 4, 0, &
    1, 2, 3, 4, &
    5 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  1,  2, &
    2,  2,  3,  3, &
    3,  3,  4,  4, &
    4,  4,  4,  5, &
    5,  5,  5,  5, &
    5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    m = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    m = m_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine pn_pair_product ( p, table )

!*****************************************************************************80
!
!! PN_PAIR_PRODUCT: pair products for normalized Legendre polynomial Pn(n,x).
!
!  Discussion:
!
!    Let Pn(n,x) represent the normalized Legendre polynomial of degree n.  
!
!    To check orthonormality, we compute
!
!      Tij = Integral ( -1.0 <= X <= +1.0 ) Pn(i,x) * Pn(j,x) dx
!
!    We will estimate these integrals using Gauss-Legendre quadrature.
!
!    The computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x(1)
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call p_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x(1) = x_table(k)
    call pn_polynomial_value ( 1, p, x, h_table )

    do i = 0, p
      do j = 0, p
        table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine pn_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre Pn(n,x).
!
!  Discussion:
!
!    Pn(n,x) = P(n,x) * sqrt ( (2n+1)/2 )
!
!          1       x       x^2     x^3     x^4      x^5    x^6     x^7
!
!    0   0.707
!    1   0.000   1.224
!    2  -0.790   0.000   2.371
!    3   0.000  -2.806   0.000   4.677
!    4   0.795   0.000  -7.954   0.000   9.280
!    5   0.000   4.397   0.000 -20.520   0.000   18.468
!    6  -0.796   0.000  16.731   0.000 -50.193    0.000  36.808
!    7   0.000  -5.990   0.000  53.916   0.000 -118.616   0.000  73.429 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the 
!    normalized Legendre polynomials of degree 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) t

  if ( n < 0 ) then
    return
  end if
!
!  Compute P(i,x) coefficients.
!
  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( 0 < n ) then
    c(1,1) = 1.0D+00
  end if

  do i = 2, n
    c(i,0:i-2) =          real (   - i + 1, kind = 8 ) * c(i-2,0:i-2) &
                        / real (     i,     kind = 8 )
    c(i,1:i) = c(i,1:i) + real ( i + i - 1, kind = 8 ) * c(i-1,0:i-1) &
                        / real (     i,     kind = 8 )
  end do
!
!  Normalize them.
!
  do i = 0, n
    t = sqrt ( real ( 2 * i + 1, kind = 8 ) / 2.0D+00 )
    c(i,0:i) = c(i,0:i) * t
  end do
 
  return
end
subroutine pn_polynomial_value ( m, n, x, v )

!*****************************************************************************80
!
!! PN_POLYNOMIAL_VALUE evaluates the normalized Legendre polynomials Pn(n,x).
!
!  Discussion:
!
!    The normalized Legendre polynomials are orthonormal under the inner product
!    defined as integration from -1 to 1:
!
!      Integral ( -1 <= x <= +1 ) Pn(i,x) * Pn(j,x) dx = delta(i,j)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) norm
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  call p_polynomial_value ( m, n, x, v )

  do j = 0, n
    norm = sqrt ( 2.0D+00 / real ( 2 * j + 1, kind = 8 ) )
    v(1:m,j) = v(1:m,j) / norm
  end do
 
  return
end
function r8_factorial ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL computes the factorial of N.
!
!  Discussion:
!
!    factorial ( N ) = product ( 1 <= I <= N ) I
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 January 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the factorial function.
!    If N is less than 1, the function value is returned as 1.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL, the factorial of N.
!
  implicit none

  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n

  r8_factorial = 1.0D+00

  do i = 1, n
    r8_factorial = r8_factorial * real ( i, kind = 8 )
  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_PRINT prints an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_PRINT_SOME prints some of an R8MAT.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) A_FIRST, A_LAST, the first and last entries.
!
!    Output, real ( kind = 8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = 8 ) * a   &
             + real (     i - 1, kind = 8 ) * b ) &
             / real ( n     - 1, kind = 8 )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
  end do

  return
end
subroutine r8vec2_print ( n, a1, a2, title )

!*****************************************************************************80
!
!! R8VEC2_PRINT prints an R8VEC2.
!
!  Discussion:
!
!    An R8VEC2 is a dataset consisting of N pairs of R8's, stored
!    as two separate vectors A1 and A2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of components of the vector.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), the vectors to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  integer ( kind = 4 ) i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
  end do

  return
end

     ! Hermite library:
!subroutine h_integral ( n, value )

!*****************************************************************************80
!
!! H_INTEGRAL evaluates the integral of H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!    The integral computed is:
!
!      integral ( -oo < x < +oo ) H(i,x) exp(-x^2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the integral.  
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
!  implicit none

!  integer ( kind = 4 ) n
!  real ( kind = 8 ) r8_factorial2
!  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
!  real ( kind = 8 ) value

!  if ( mod ( n, 2 ) == 1 ) then

!    value = 0.0D+00

!  else

!    value = r8_factorial2 ( n - 1 ) * sqrt ( r8_pi ) / 2.0D+00 ** ( n / 2 )

!  end if

!  return
!end
subroutine h_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! H_POLYNOMIAL_COEFFICIENTS: coefficients of H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     2
!     2     -2     0      4
!     3      0   -12      0      8
!     4     12     0    -48      0      16
!     5      0   120      0   -160       0    32
!     6   -120     0    720      0    -480     0     64
!     7      0 -1680      0   3360       0 -1344      0   128
!     8   1680     0 -13440      0   13440     0  -3584     0    256
!     9      0 30240      0 -80640       0 48384      0 -9216      0 512
!    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the polynomials
!    of orders 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 2.0D+00
 
  do i = 2, n
    c(i,0)     =  -2.0D+00 * real ( i - 1, kind = 8 ) * c(i-2,0)
    c(i,1:i-2) =   2.0D+00                            * c(i-1,0:i-3)  &
                  -2.0D+00 * real ( i - 1, kind = 8 ) * c(i-2,1:i-2)
    c(i,  i-1) =   2.0D+00                            * c(i-1,  i-2)
    c(i,  i  ) =   2.0D+00                            * c(i-1,  i-1)
  end do
 
  return
end
subroutine h_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! H_POLYNOMIAL_VALUE evaluates H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!    H(N,X) = (-1)^N * exp ( X^2 ) * dn/dXn ( exp(-X^2 ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 October 2002
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Larry Andrews,
!    Special Functions of Mathematics for Engineers,
!    Second Edition, 
!    Oxford University Press, 1998.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) P(M,0:N), the values of the first N+1 Hermite
!    polynomials at the point X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = 2.0D+00 * x(1:m)
 
  do j = 2, n
    p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
             - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
 
  return
end
subroutine h_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! H_POLYNOMIAL_VALUES: tabulated values of H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!    In Mathematica, the function can be evaluated by:
!
!      HermiteH[n,x]
!
!  Differential equation:
!
!    Y'' - 2 X Y' + 2 N Y = 0
!
!  First terms:
!
!      1
!      2 X
!      4 X^2     -  2
!      8 X^3     - 12 X
!     16 X^4     - 48 X^2     + 12
!     32 X^5    - 160 X^3    + 120 X
!     64 X^6    - 480 X^4    + 720 X^2    - 120
!    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
!    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
!    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
!   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
!
!  Recursion:
!
!    H(0,X) = 1,
!    H(1,X) = 2*X,
!    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
!    = sqrt ( PI ) * 2^N * N!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
      0.1000000000000000D+01, &
      0.1000000000000000D+02, &
      0.9800000000000000D+02, &
      0.9400000000000000D+03, &
      0.8812000000000000D+04, &
      0.8060000000000000D+05, &
      0.7178800000000000D+06, &
      0.6211600000000000D+07, &
      0.5206568000000000D+08, &
      0.4212712000000000D+09, &
      0.3275529760000000D+10, &
      0.2432987360000000D+11, &
      0.1712370812800000D+12, &
      0.0000000000000000D+00, &
      0.4100000000000000D+02, &
     -0.8000000000000000D+01, &
      0.3816000000000000D+04, &
      0.3041200000000000D+07 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5, 5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine h_polynomial_zeros ( nt, z )

!*****************************************************************************80
!
!! H_POLYNOMIAL_ZEROS: zeros of H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) Z(NT), the zeros of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) z(nt)

  z(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, z, bj, wts )

  return
end
subroutine h_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! H_QUADRATURE_RULE: quadrature for H(i,x).
!
!  Discussion:
!
!    H(i,x) is the physicist's Hermite polynomial of degree I.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), WTS(NT), the points and weights 
!    of the rule.
!
  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2

  return
end
function he_double_product_integral ( i, j )

!*****************************************************************************80
!
!! HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
!    Princeton, 2010,
!    ISBN13: 978-0-691-14212-8,
!    LC: QA274.23.X58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, the polynomial indices.
!
!    Output, real ( kind = 8 ) HE_DOUBLE_PRODUCT_INTEGRAL, the value of 
!    the integral.
!
  implicit none

  real ( kind = 8 ) he_double_product_integral
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r8_factorial
  real ( kind = 8 ) value

  if ( i == j ) then
    value = r8_factorial ( i )
  else
    value = 0.0D+00
  end if

  he_double_product_integral = value

  return
end

!subroutine he_integral ( n, value )

!*****************************************************************************80
!
!! HE_INTEGRAL evaluates the integral of He(i,x).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!    The integral computed is:
!
!      integral ( -oo < x < +oo ) He(i,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the integral.  
!    0 <= N.
!
!    Output, real ( kind = 8 ) VALUE, the value of the integral.
!
!  implicit none

!  integer ( kind = 4 ) n
!  real ( kind = 8 ) r8_factorial2
!  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
!  real ( kind = 8 ) value

!  if ( mod ( n, 2 ) == 1 ) then

!    value = 0.0D+00

!  else

!    value = r8_factorial2 ( n - 1 ) * sqrt ( 2.0D+00 * r8_pi )

!  end if

!  return
!end
subroutine he_polynomial_coefficients ( n, c )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_COEFFICIENTS: coefficients of He(i,x).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!  First terms:
!
!    N/K     0     1      2      3       4     5      6    7      8    9   10
!
!     0      1
!     1      0     1
!     2     -1     0      1
!     3      0    -3      0      1
!     4      3     0     -6      0       1
!     5      0    15      0    -10       0     1
!     6    -15     0     45      0     -15     0      1
!     7      0  -105      0    105       0   -21      0     1
!     8    105     0   -420      0     210     0    -28     0      1
!     9      0   945      0  -1260       0   378      0   -36      0   1
!    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Output, real ( kind = 8 ) C(0:N,0:N), the coefficients of the polynomials
!    of orders 0 through N.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) c(0:n,0:n)
  integer ( kind = 4 ) i

  if ( n < 0 ) then
    return
  end if

  c(0:n,0:n) = 0.0D+00

  c(0,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  c(1,1) = 1.0D+00
 
  do i = 2, n
    c(i,0)     =              - real ( i - 1, kind = 8 ) * c(i-2,0)
    c(i,1:i-2) = c(i-1,0:i-3) - real ( i - 1, kind = 8 ) * c(i-2,1:i-2)
    c(i,  i-1) = c(i-1,  i-2)
    c(i,  i  ) = c(i-1,  i-1)
  end do
 
  return
end
subroutine he_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_VALUE evaluates He(i,x).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!  Differential equation:
!
!    ( exp ( - 0.5 * x^2 ) * He(n,x)' )' + n * exp ( - 0.5 * x^2 ) * He(n,x) = 0
!
!  First terms:
!
!   1
!   X
!   X^2  -  1
!   X^3  -  3 X
!   X^4  -  6 X^2 +   3
!   X^5  - 10 X^3 +  15 X
!   X^6  - 15 X^4 +  45 X^2 -   15
!   X^7  - 21 X^5 + 105 X^3 -  105 X
!   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
!   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
!   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
!
!  Recursion:
!
!    He(0,X) = 1,
!    He(1,X) = X,
!    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
!
!  Orthogonality:
!
!    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
!    = sqrt ( 2 * pi ) * N! * delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real X(M), the evaluation points.
!
!    Output, real P(M,0:N), the values of the probabilist's Hermite polynomials 
!    of index 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  p(1:m,0) = 1.0D+00

  if ( 0 < n ) then
    p(1:m,1) = x(1:m)
  end if

  do j = 2, n
    p(1:m,j) = x(1:m) * p(1:m,j-1) - real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
 
  return
end
subroutine he_polynomial_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_VALUES: tabulated values of He(i,x).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!    In Mathematica, the function can be evaluated by:
!
!      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
!
!  First terms:
!
!   1
!   X
!   X^2  -  1
!   X^3  -  3 X
!   X^4  -  6 X^2 +   3
!   X^5  - 10 X^3 +  15 X
!   X^6  - 15 X^4 +  45 X^2 -   15
!   X^7  - 21 X^5 + 105 X^3 -  105 X
!   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
!   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
!   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
!
!  Recursion:
!
!    He(0,X) = 1,
!    He(1,X) = X,
!    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
!
!  Norm:
!
!    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
!    = sqrt ( 2 * pi ) * N! * delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 18

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.000000000000000D+00, &
    5.000000000000000D+00, &
    24.00000000000000D+00, &
    110.0000000000000D+00, &
    478.0000000000000D+00, &
    1950.000000000000D+00, &
    7360.000000000000D+00, &
    25100.00000000000D+00, &
    73980.00000000000D+00, &
    169100.0000000000D+00, &
    179680.0000000000D+00, &
   -792600.0000000000D+00, &
   -5939480.000000000D+00, &
    0.000000000000000D+00, &
    6.281250000000000D+00, &
    6.000000000000000D+00, &
    18.00000000000000D+00, &
    90150.00000000000D+00 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
     0,  1,  2, &
     3,  4,  5, &
     6,  7,  8, &
     9, 10, 11, &
    12,  5,  5, &
     5,  5,  5 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    5.0D+00, &
    0.0D+00, &
    0.5D+00, &
    1.0D+00, &
    3.0D+00, &
    1.0D+01 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine he_polynomial_zeros ( nt, z )

!*****************************************************************************80
!
!! HE_POLYNOMIAL_ZEROS: zeros of He(i,x).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the degree of the polynomial.
!
!    Output, real ( kind = 8 ) Z(NT), the zeros of the polynomial.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) wts(nt)
  real ( kind = 8 ) z(nt)

  z(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, z, bj, wts )

  z(1:nt) = z(1:nt) * sqrt ( 2.0D+00 )

  return
end
subroutine he_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! HE_QUADRATURE_RULE: quadrature for He(i,x).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), WTS(NT), the points and weights 
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, t, bj, wts )

  t(1:nt) = t(1:nt) * sqrt ( 2.0D+00 )
  wts(1:nt) = wts(1:nt) ** 2 * sqrt ( 2.0D+00 )

  return
end
function he_triple_product_integral ( i, j, k )

!*****************************************************************************80
!
!! HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
!
!  Discussion:
!
!    He(i,x) represents the probabilist's Hermite polynomial.
!
!    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
!    Princeton, 2010,
!    ISBN13: 978-0-691-14212-8,
!    LC: QA274.23.X58.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, J, K, the polynomial indices.
!
!    Output, real ( kind = 8 ) HE_TRIPLE_PRODUCT_INTEGRAL, the value 
!    of the integral.
!
  implicit none

  real ( kind = 8 ) he_triple_product_integral
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_factorial
  integer ( kind = 4 ) s
  real ( kind = 8 ) value

  s = ( i + j + k ) / 2

  if ( s < max ( i, j, k ) ) then
    value = 0.0D+00
  else if ( mod ( i + j + k, 2 ) /= 0 ) then
    value = 0.0D+00
  else
    value = r8_factorial ( i ) / r8_factorial ( s - i ) &
          * r8_factorial ( j ) / r8_factorial ( s - j ) &
          * r8_factorial ( k ) / r8_factorial ( s - k )
  end if

  he_triple_product_integral = value

  return
end
subroutine hen_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! HEN_EXPONENTIAL_PRODUCT: exponential product exp(b*x)*Hen(i,x)*Hen(j,x).
!
!  Discussion:
!
!    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) 
!        exp(B*X) * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact, and moreoever, the
!    table will be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the 
!    polyonomial factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the weighted integral of 
!    exp(B*X) * Hen(I,X) * Hen(J,X).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call he_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k);
    call hen_polynomial_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x ) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine hen_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! HEN_POLYNOMIAL_VALUE evaluates Hen(i,x).
!
!  Discussion:
!
!    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
!
!    These polynomials satisfy the orthonormality condition:
!
!      Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * Hen(M,X) Hen(N,X) dX 
!      = delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) P(M,0:N), the values of the polynomials of 
!    index 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m)

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = x(1:m)
 
  do j = 2, n
    p(1:m,j) = x(1:m) * p(1:m,j-1) - real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
!
!  Normalize.
!
  fact = 1.0D+00
  do j = 0, n
    p(1:m,j) = p(1:m,j) / sqrt ( fact * sqrt ( 2.0D+00 * r8_pi ) )
    fact = fact * real ( j + 1, kind = 8 )
  end do

  return
end
subroutine hen_power_product ( p, e, table )

!*****************************************************************************80
!
!! HEN_POWER_PRODUCT: power products, x^e*Hen(i,x)*Hen(j,x).
!
!  Discussion:
!
!    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) 
!        X^E * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!
!    When E is 0, the computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, integer ( kind = 4 ) E, the exponent of X in the integrand.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the weighted integral of 
!    X^E * Hen(I,X) * Hen(J,X).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  integer ( kind = 4 ) e
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( ( e + 1 ) / 2 )

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call he_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hen_polynomial_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    if ( e == 0 ) then
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
        end do
      end do
    else
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * x ** e * h_table(i) * h_table(j)
        end do
      end do
    end if

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine hf_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! HF_EXPONENTIAL_PRODUCT: exponential products, exp(b*x)*Hf(i,x)*Hf(j,x).
!
!  Discussion:
!
!    Hf(i,x) represents the Hermite function of "degree" I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) exp(B*X) * Hf(I,X) * Hf(J,X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact, and moreoever, the
!    table will be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the integral of exp(B*X) * Hf(I,X) * Hf(J,X).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call hf_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hf_function_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x ) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine hf_function_value ( m, n, x, f )

!*****************************************************************************80
!
!! HF_FUNCTION_VALUE evaluates Hf(i,x).
!
!  Discussion:
!
!    Hf(i,x) represents the Hermite function of "degree" I.  
!
!    The Hermite function of degree n is related to the physicist's
!    Hermite polynomial H(n,x):
!
!      Hf(n,x) = H(n,x) * exp ( - 0.5 * x^2 ) / sqrt ( 2^n n! sqrt ( pi ) )
!
!    The Hermite functions are orthonormal:
!
!      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) F(M,0:N), the values of the Hermite functions 
!    of index 0 through N at the evaluation points.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) f(m,0:n)
  integer ( kind = 4 ) j
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) x(m)

  f(1:m,0) = exp ( - 0.5D+00 * x(1:m)**2 ) / sqrt ( sqrt ( r8_pi ) )

  if ( n == 0 ) then
    return
  end if

  f(1:m,1) = 2.0D+00 * exp ( - 0.5D+00 * x(1:m)**2 ) * x(1:m) &
    / sqrt ( 2.0D+00 * sqrt ( r8_pi ) )

  do j = 2, n
    f(1:m,j) = ( sqrt ( 2.0D+00 ) * x(1:m) * f(1:m,j-1) &
      - sqrt ( real ( j - 1, kind = 8 ) ) * f(1:m,j-2) ) &
      / sqrt ( real ( j, kind = 8 ) )
  end do
 
  return
end
subroutine hf_function_values ( n_data, n, x, fx )

!*****************************************************************************80
!
!! HF_FUNCTION_VALUES: tabulated values of Hf(i,x).
!
!  Discussion:
!
!    Hf(i,x) represents the Hermite function of "degree" I.  
!
!    In Mathematica, the function can be evaluated by:
!
!      Hf(n,x) = HermiteH[n,x] 
!        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
!
!    The Hermite functions are orthonormal:
!
!      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) N_DATA.  The user sets N_DATA to 0
!    before the first call.  On each call, the routine increments N_DATA by 1,
!    and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, integer ( kind = 4 ) N, the order of the polynomial.
!
!    Output, real ( kind = 8 ) X, the point where the polynomial is evaluated.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 23

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.7511255444649425D+00,  0.0000000000000000D+00, -0.5311259660135985D+00, &
    0.0000000000000000D+00,  0.4599685791773266D+00,  0.0000000000000000D+00, &
    0.4555806720113325D+00,  0.6442883651134752D+00,  0.3221441825567376D+00, &
   -0.2630296236233334D+00, -0.4649750762925110D+00, -0.5881521185179581D-01, &
    0.3905052515434106D+00,  0.2631861423064045D+00, -0.2336911435996523D+00, &
   -0.3582973361472840D+00,  0.6146344487883041D-01,  0.3678312067984882D+00, &
    0.9131969309166278D-01,  0.4385750950032321D+00, -0.2624689527931006D-01, &
    0.5138426125477819D+00,  0.9355563118061758D-01 /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ), save, dimension ( n_max ) :: n_vec = (/ &
    0,  1,  2,  &
    3,  4,  5,  &
    0,  1,  2,  &
    3,  4,  5,  &
    6,  7,  8,  &
    9, 10, 11,  &
   12,  5,  5,  &
    5,  5  /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 0.5D+00, 2.0D+00, &
    3.0D+00, 4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    n = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    n = n_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine hf_power_product ( p, e, table )

!*****************************************************************************80
!
!! HF_POWER_PRODUCT: power products x^e*Hf(i,x)*Hf(j,x).
!
!  Discussion:
!
!    Hf(i,x) represents the Hermite function of "degree" I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) X^E * Hf(I,X) * Hf(J,X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!
!    When E is 0, the computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, integer ( kind = 4 ) E, the exponent of X in the integrand.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the integral of X^E * Hf(I,X) * Hf(J,X).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  integer ( kind = 4 ) e
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( ( e + 1 ) / 2 )

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call hf_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hf_function_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    if ( e == 0 ) then
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
        end do
      end do
    else
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * x ** e * h_table(i) * h_table(j)
        end do
      end do
    end if

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine hf_quadrature_rule ( nt, t, wts )

!*****************************************************************************80
!
!! HF_QUADRATURE_RULE: quadrature for Hf(i,x).
!
!  Discussion:
!
!    Hf(i,x) represents the Hermite function of "degree" I.  
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NT, the order of the rule.
!
!    Output, real ( kind = 8 ) T(NT), WTS(NT), the points and weights
!    of the rule.
!
  implicit none

  integer ( kind = 4 ) nt

  real ( kind = 8 ) bj(nt)
  integer ( kind = 4 ) i
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) t(nt)
  real ( kind = 8 ) wts(nt)

  t(1:nt) = 0.0D+00

  do i = 1, nt
    bj(i) = sqrt ( real ( i, kind = 8 ) / 2.0D+00 )
  end do

  wts(1:nt) = 0.0D+00
  wts(1) = sqrt ( sqrt ( r8_pi ) )

  call imtqlx ( nt, t, bj, wts )

  wts(1:nt) = wts(1:nt) ** 2 * exp ( t(1:nt) ** 2 )

  return
end
subroutine hn_exponential_product ( p, b, table )

!*****************************************************************************80
!
!! HN_EXPONENTIAL_PRODUCT: exponential products exp(b*x)*Hn(i,x)*Hn(j,x).
!
!  Discussion:
!
!    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of exp(B*X) with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) 
!        exp(B*X) * Hn(I,X) * Hn(J,X) exp(-X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!    Because of the exponential factor exp(B*X), the quadrature will not 
!    be exact.
!
!    However, when B = 0, the quadrature is exact, and moreoever, the
!    table will be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, real ( kind = 8 ) B, the coefficient of X in the exponential factor.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the weighted integral of 
!    exp(B*X) * Hn(I,X) * Hn(J,X).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = ( 3 * p + 4 ) / 2

  allocate ( x_table(1:order) )
  allocate ( w_table(1:order) )

  call h_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hn_polynomial_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    do j = 0, p
      do i = 0, p
        table(i,j) = table(i,j) &
          + w_table(k) * exp ( b * x ) * h_table(i) * h_table(j)
      end do
    end do

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end
subroutine hn_polynomial_value ( m, n, x, p )

!*****************************************************************************80
!
!! HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
!
!  Discussion:
!
!    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.
!
!    These polynomials satisfy the orthonormality condition:
!
!      Integral ( -oo < X < +oo ) 
!        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
!    NIST Handbook of Mathematical Functions,
!    Cambridge University Press, 2010,
!    ISBN: 978-0521192255,
!    LC: QA331.N57.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to compute.
!    Note that polynomials 0 through N will be computed.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) P(M,0:N), the values of the polynomials of 
!    index 0 through N.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) fact
  integer ( kind = 4 ) j
  real ( kind = 8 ) p(m,0:n)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ) two_power
  real ( kind = 8 ) x(m)

  p(1:m,0) = 1.0D+00

  if ( n == 0 ) then
    return
  end if

  p(1:m,1) = 2.0D+00 * x(1:m)
 
  do j = 2, n
    p(1:m,j) = 2.0D+00 * x(1:m) * p(1:m,j-1) &
      - 2.0D+00 * real ( j - 1, kind = 8 ) * p(1:m,j-2)
  end do
!
!  Normalize.
!
  fact = 1.0D+00
  two_power = 1.0D+00
  do j = 0, n
    p(1:m,j) = p(1:m,j) / sqrt ( fact * two_power * sqrt ( r8_pi ) )
    fact = fact * real ( j + 1, kind = 8 )
    two_power = two_power * 2.0D+00
  end do

  return
end
subroutine hn_power_product ( p, e, table )

!*****************************************************************************80
!
!! HN_POWER_PRODUCT: power products x^e*Hn(i,x)*Hn(j,x).
!
!  Discussion:
!
!    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.
!
!    For polynomial chaos applications, it is of interest to know the
!    value of the integrals of products of X with every possible pair
!    of basis functions.  That is, we'd like to form
!
!      Tij = Integral ( -oo < X < +oo ) X^E * Hn(I,X) * Hn(J,X) exp(-X*X) dx
!
!    We will estimate these integrals using Gauss-Hermite quadrature.
!
!    When E is 0, the computed table should be the identity matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) P, the maximum degree of the polyonomial 
!    factors.  0 <= P.
!
!    Input, integer ( kind = 4 ) E, the exponent of X in the integrand.
!    0 <= E.
!
!    Output, real ( kind = 8 ) TABLE(0:P,0:P), the table of integrals.  
!    TABLE(I,J) represents the weighted integral of X^E * Hn(I,X) * Hn(J,X).
!
  implicit none

  integer ( kind = 4 ) p

  real ( kind = 8 ) b
  integer ( kind = 4 ) e
  real ( kind = 8 ) h_table(0:p)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) order
  real ( kind = 8 ) table(0:p,0:p)
  real ( kind = 8 ), allocatable :: w_table(:)
  real ( kind = 8 ) x
  real ( kind = 8 ), allocatable :: x_table(:)

  table(0:p,0:p) = 0.0D+00

  order = p + 1 + ( ( e + 1 ) / 2 )

  allocate ( x_table(order) )
  allocate ( w_table(order) )

  call h_quadrature_rule ( order, x_table, w_table )

  do k = 1, order

    x = x_table(k)
    call hn_polynomial_value ( 1, p, x, h_table )
!
!  The following formula is an outer product in H_TABLE.
!
    if ( e == 0 ) then
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) + w_table(k) * h_table(i) * h_table(j)
        end do
      end do
    else
      do i = 0, p
        do j = 0, p
          table(i,j) = table(i,j) &
            + w_table(k) * x ** e * h_table(i) * h_table(j)
        end do
      end do
    end if

  end do

  deallocate ( w_table )
  deallocate ( x_table )

  return
end


!function r8_factorial2 ( n )

!*****************************************************************************80
!
!! R8_FACTORIAL2 computes the double factorial function.
!
!  Discussion:
!
!    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
!                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
!
!  Example:
!
!     N    Factorial2(N)
!
!     0     1
!     1     1
!     2     2
!     3     3
!     4     8
!     5    15
!     6    48
!     7   105
!     8   384
!     9   945
!    10  3840
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the argument of the double factorial
!    function.  If N is less than 1, the value is returned as 1.0.
!
!    Output, real ( kind = 8 ) R8_FACTORIAL2, the value.
!
!  implicit none

!  integer ( kind = 4 ) n
!  real ( kind = 8 ) r8_factorial2
!  real ( kind = 8 ) r8_n

!  if ( n < 1 ) then
!    r8_factorial2 = 1.0D+00
!    return
!  end if

!  r88n = real ( n, kind = 8 )
!  r8_factorial2 = 1.0D+00

!  do while ( 1.0D+00 < r8_n )
!    r8_factorial2 = r8_factorial2 * r8_n
!    r8_n = r8_n - 2.0D+00
!  end do

!  return
!end

