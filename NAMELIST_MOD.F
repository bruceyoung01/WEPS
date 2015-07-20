! $ID: NAMELIST_MOD.F V1.0 06/21/2012 10:10 BRUCE EXP$
!
!******************************************************************************
!  MODULE NAMELIST_MOD READS THE NAMELIST.WEPS FILE IN THE 'RUN' DIRECTORY AT 
!  THE START RUN AND PASSES THE INFO TO SEVERAL OTHER MODULES, SUCH AS 
!  'NEI_MOD', AND 'FLAMBE_MOD'.
!
!  MODULE VARIABLES:
!  ============================================================================
!  (1 ) 
!
!  MODULE ROUTINES:
!  ============================================================================
!  (1 ) READ_NAMELIST              : MAIN SUBROUTINE TO READ THE NAMELIST.WEPS
!  (2 ) READ_ONE_LINE              : READ ONE LINE OF THE NAMELIST.WEPS
!  (3 ) SPLIT_ONE_LINE             : SPLIT ONE LINE INTO A SET
!  (4 ) READ_TIME_MENU             : READ TIME MENU
!  (5 ) READ_DOMAIN_MENU           : READ DOMAIN MENU
!  (6 ) READ_EMISSION_MENU         : READ SWITCH OPTIONS OF EMISSION
!  (7 ) READ_INJETION_HEIGHT_MENU  : READ INJECTION HEIGHT MENU
!  (8 ) READ_SMOKE_SPECIES_MENU    : READ SMOKE SPECIES MENU
!  (9 ) READ_DIR_MENU              : READ DIRECTORY MENU
!  (10) READ_DIAGNOSTIC_MENU       : READ DIAGNOSTIC MENU
!  (11) READ_DEBUG_MENU            : READ DEBUG MENU
!
!  NOTES:
!  ============================================================================
!  (1 ) ORGINALLY WRITTEN BY BRUCE. (06/21/2012)
!  (2 ) MODIFIED FROM RICHARD'S MODIFICATION FROM GEOS-CHEM MODEL. (06/21/2012)
!  (3 ) MODIFIED READ_EMISSION_MENU BY BRUCE. (07/10/2012)
!  (4 ) ADD BASE SPECIES TO CALCULATE EMISSION RATIO BY BRUCE. (01/13/2014)
!******************************************************************************
!
      MODULE NAMELIST_MOD 

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES
      ! AND ROUTINES FROM BEING SEEN OUTSIDE 'NAMELIST_MOD.F'
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ...EXCEPT THESE ROUTINES
      PUBLIC :: READ_NAMELIST
      LOGICAL, PUBLIC :: F_PM25     ! if .true., output PM2.5
                                    ! otherwise, OC and BC


      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL               :: VERBOSE  = .FALSE.
      INTEGER, PARAMETER    :: FIRSTCOL = 31
      INTEGER, PARAMETER    :: MAXDIM   = 255
      INTEGER, PARAMETER    :: IU_NAME  = 11
      CHARACTER (LEN = 255) :: FILENAME = '../run/namelist.weps'

      !=================================================================
      ! MODULE ROUTINES -- FOLLOW BELOW THE 'CONTAINS' STATEMENT
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

!  $ID: READ_NAMELIST V01 06/21/2012 10:41 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_NAMELIST READS INPUT FILE 'namelist.weps'
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL. (06/21/2012)
!******************************************************************************
!

      SUBROUTINE READ_NAMELIST

      !REFERENCES TO F90 MODULES
      USE CHARPAK_MOD, ONLY : STRREPL
      USE ERROR_MOD,   ONLY : ERROR_STOP

      ! LOCAL VARIABLES
      LOGICAL               :: EOF
      INTEGER               :: IOS
      CHARACTER (LEN = 255) :: TOPTITLE
      CHARACTER (LEN = 255) :: LINE
      CHARACTER (LEN = 255) :: MESSAGE
      CHARACTER (LEN = 1)   :: TAB   = ACHAR(9)
      CHARACTER (LEN = 1)   :: SPACE = ACHAR(32)

      !=================================================================
      ! READ_NAMELIST BEGINS HERE!
      !=================================================================

      ! SCREEN PRINT
      WRITE(6, '(A  )') REPEAT('=', 79)
      WRITE(6, 100    ) TRIM(FILENAME)
100   FORMAT('READ_NAMELIST: Reading ', A)

      ! OPEN FILE
      OPEN(IU_NAME, FILE = TRIM(FILENAME), STATUS = 'OLD', IOSTAT = IOS)
      IF (IOS /= 0) THEN
       MESSAGE = 'Fail to open file: '//TRIM(FILENAME)
       CALL ERROR_STOP(MESSAGE, 'NAMELIST_MOD.F')
      ENDIF

      ! READ TOPTITLE FOR BINARY PUNCH FILE
      TOPTITLE = READ_ONE_LINE(EOF)
      IF (EOF) RETURN

      ! LOOP UNTIL EOF
      DO

       ! READ A LINE FROM THE FILE namelist.weps, EXIT IF EOF
       LINE = READ_ONE_LINE(EOF)
       IF (EOF) EXIT

        ! REPLACE TAB CHARACTERS IN LINE (IF ANY) W/ SPACES
        CALL STRREPL(LINE, TAB, SPACE)

       IF      (INDEX(LINE, 'TIME MENU'            ) >0) THEN
        CALL READ_TIME_MENU

       ELSE IF (INDEX(LINE, 'DOMAIN MENU'          ) >0) THEN
        CALL READ_DOMAIN_MENU

       ELSE IF (INDEX(LINE, 'EMISSION MENU'        ) >0) THEN
        CALL READ_EMISSION_MENU

       ELSE IF (INDEX(LINE, 'INJECTION HEIGHT MENU') >0) THEN
        CALL READ_INJECTION_HEIGHT_MENU

       ELSE IF (INDEX(LINE, 'SMOKE SPECIES MENU'   ) >0) THEN
        CALL READ_SMOKE_SPECIES_MENU

       ELSE IF (INDEX(LINE, 'DIRECTORY MENU'       ) >0) THEN
        CALL READ_DIR_MENU

       ELSE IF (INDEX(LINE, 'DIAGNOSTIC MENU'      ) >0) THEN
        CALL READ_DIAGNOSTIC_MENU

       ELSE IF (INDEX(LINE, 'DEBUG MENU'           ) >0) THEN
        CALL READ_DEBUG_MENU

       ELSE IF (INDEX(LINE, 'END OF FILE'          ) >0) THEN
        EXIT

       ENDIF

      ENDDO

      ! CLOSE INPUT FILE
      CLOSE (IU_NAME)

      END SUBROUTINE READ_NAMELIST

!------------------------------------------------------------------------------
!
!  $ID: READ_ONE_LINE V01 06/21/2012 12:01 BRUCE EXP$
!
!******************************************************************************
!  FUNCTION READS A LINE FROM THE INPUT FILE. IF THE GLOBAL VARIABLE
!  VERBOSE IS SET, THE LINE WILL BE PRINTED TO STDOUT. READ_ONE_LINE CAN
!  TRAP AN UNEXPECTED EOF IF LOCATION IS PASSED. OTHERWISE, IT WILL PASS
!  A LOGICAL FLAG BACK TO THE CALLING ROUTINE, WHERE THE ERROR TRAPPING
!  WILL BE DONE.
!
!  ARGUMENTS AS OUTPUT:
!  ============================================================================
!  (1 ) EOF       (CHARACTER) : LOGICAL FLAG DENOTING EOF CONDITION
!  (2 ) LOCATION  (CHARACTER) : NAME OF CALLING ROUTINE; TRAPS PREMATURE EOF
!
!  FUNCTION VALUE:
!  ============================================================================
!  (1 ) LINE      (CHARACTER) : A LINE OF TEXT AS READ FROM THE FILE
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM GEOS-CHEM MODEL BY BRUCE. (06/21/2012)
!******************************************************************************
!

      FUNCTION READ_ONE_LINE(EOF, LOCATION) RESULT(LINE)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! ARGUMENTS
      LOGICAL,             INTENT(OUT)          :: EOF
      CHARACTER (LEN = *), INTENT(IN), OPTIONAL :: LOCATION

      ! LOCAL VARIABLES
      INTEGER                                   :: IOS
      CHARACTER (LEN = 255)                     :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE BEGINS HERE!
      !=================================================================

      ! INITIALIZE
      EOF = .FALSE.

      ! READ A LINE FROM THE FILE
      READ(IU_NAME, '(A)', IOSTAT = IOS) LINE

      ! IO STATUS < 0: EOF CONDITION
      IF (IOS < 0) THEN
       EOF = .TRUE.

       ! TRAP UNEXPECTED EOF -- STOP W/ ERROR MSG IF LOCATION IS PASSED
       ! OTHERWISE, RETURN EOF TO THE CALLING PROGRAM
       IF (PRESENT(LOCATION)) THEN
        MSG = 'READ_ONE_LINE: error at: '//TRIM(LOCATION)
        WRITE(6, '(A)') MSG
        WRITE(6, '(A)') 'Unexpected end of file encountered!'
        WRITE(6, '(A)') 'Stop in READ_ONE_LINE (NAMELIST_MOD.F)'
        WRITE(6, '(A)') REPEAT('=', 79)
        STOP
       ELSE
        RETURN
       ENDIF
      ENDIF

      ! IO STATUS > 0: TRUE I/O ERROR CONDITION
      IF (IOS .GT. 0) CALL ERROR_STOP('Error when reading', 
     &                                'READ_ONE_LINE: 1')

      ! PRINT THE LINE (IF NECESSARY)
      IF (VERBOSE) WRITE(6, '(A)') TRIM(LINE)

      ! RETURN TO CALLING PROGRAM
      END FUNCTION READ_ONE_LINE

!------------------------------------------------------------------------------
!
!  $ID: SPLIT_ONE_LINE V01 06/21/2012 15:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE SPLIT_ONE_LINE READS A LINE FROM THE INPUT FILE (VIA ROUTINE 
!  READ_ONE_LINE), AND SEPARATES IT INTO SUBSTRINGS.
!  SPLIT_ONE_LINE ALSO CHECKS TO SEE IF THE NUMBER OF SUBSTRINGS FOUND
!  IS EQUAL TO THE NUMBER OF SUBSTRINGS THAT WE EXPECTED TO FIND.
!  HOWEVER, IF YOU DON'T KNOW A-PRIORI HOW MANY SUBSTRINGS TO EXPECT
!  A-PRIORI, YOU CAN SKIP THE ERROR CHECK.
!
!  ARGUMENTS AS INPUT:
!  ============================================================================
!  (3 ) N_EXP      (INTEGER  ) : NUMBER OF SUBSTRINGS WE EXPECT TO FIND
!                                (N_EXP < 0 WILL SKIP THE ERROR CHECK!)
!  (4 ) LOCATION   (CHARACTER) : NAME OF ROUTINE THAT CALLED SPLIT_ONE_LINE
!
!  ARGUMENTS AS OUTPUT:
!  ============================================================================
!  (1 ) SUBSTRS    (CHARACTER) : ARRAY OF SUBSTRINGS (SEPARATED BY " ")
!  (2 ) N_SUBSTRS  (INTEGER  ) : NUMBER OF SUBSTRINGS ACTUALLY FOUND
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFED FROM GEOS-CHEM MODEL BY BRUCE. (06/21/2012)
!******************************************************************************
!

      SUBROUTINE SPLIT_ONE_LINE(SUBSTRS, N_SUBSTRS, N_EXP, LOCATION)

      ! REFERENCES TO F90 MODULES
      USE CHARPAK_MOD, ONLY : STRSPLIT

      ! ARGUMENTS
      CHARACTER (LEN = 255), INTENT(OUT) :: SUBSTRS(MAXDIM)
      INTEGER,               INTENT(OUT) :: N_SUBSTRS
      INTEGER,               INTENT(IN)  :: N_EXP
      CHARACTER (LEN = *  ), INTENT(IN)  :: LOCATION

      ! LOCAL VARIABLES
      LOGICAL                            :: EOF
      CHARACTER (LEN = 255)              :: LINE, MSG

      !=================================================================
      ! SPLIT_ONE_LINE BEGINS HERE!
      !=================================================================

      ! CREATE ERROR MSG
      MSG = 'SPLIT_ONE_LINE: error at '//TRIM(LOCATION)

      !=================================================================
      ! READ A LINE FROM DISK
      !=================================================================
      LINE = READ_ONE_LINE(EOF)

      ! STOP ON END-OF-FILE W/ ERROR MSG
      IF (EOF) THEN
       WRITE(6, '(A)') TRIM(MSG)
       WRITE(6, '(A)') 'End of file encountered!'
       WRITE(6, '(A)') 'Stop in SPLIT_ONE_LINE (NAMELIST_MOD.F)'
       WRITE(6, '(A)') REPEAT('=', 79)
      ENDIF

      !=================================================================
      !  SPLIT THE LINES BETWEEN SPACES -- START AT COLUMN FIRSTCOL
      !=================================================================
      CALL STRSPLIT(LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS)

      ! SOMETIMES WE DON'T KNOW HOW MANY SUBSTRINGS TO EXPECT,
      ! IF N_EXP IS GREATER THAN MAXDIM, THEN SKIP THE ERROR CHECK
      IF (N_EXP < 0) RETURN

      ! RETURN TO CALLING PROGRAM
      END SUBROUTINE SPLIT_ONE_LINE

!------------------------------------------------------------------------------
!  $ID: READ_TIME_MENU V01 06/22/2012 09:23 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_TIME_MENU READS TIME VARIABLES OF TIME MENU IN THE
!  FILE NAMELIST.WEPS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) MAX_DOM     (INTEGER) : # OF DOMAINS                           [-----]
!  (2 ) START_YEAR  (INTEGER) : START YEAR OF EMISSION DATA            [year ]
!  (3 ) START_MONTH (INTEGER) : START MONTH OF EMISSION DATA           [month]
!  (4 ) START_DAY   (INTEGER) : START DAY OF EMISSION DATA             [day  ]
!  (5 ) START_HOUR  (INTEGER) : START HOUR OF EMISSION DATA            [hour ]
!  (6 ) END_YEAR    (INTEGER) : END YEAR OF EMISSION DATA              [year ]
!  (7 ) END_MONTH   (INTEGER) : END MONTH OF EMISSION DATA             [month]
!  (8 ) END_DAY     (INTEGER) : END DAY OF EMISSION DATA               [day  ]
!  (9 ) END_HOUR    (INTEGER) : END HOUR OF EMISSION DATA              [hour ]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL BY BRUCE. (06/22/2012)
!******************************************************************************
!
      SUBROUTINE READ_TIME_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : MAX_DOM
      USE NAMELIST_ARRAY_MOD, ONLY : START_YEAR, START_MONTH
      USE NAMELIST_ARRAY_MOD, ONLY : START_DAY,  START_HOUR
      USE NAMELIST_ARRAY_MOD, ONLY : END_YEAR,   END_MONTH
      USE NAMELIST_ARRAY_MOD, ONLY : END_DAY,    END_HOUR
      USE ERROR_MOD,          ONLY : ERROR_STOP

      ! LOCAL VARIABLES
      INTEGER               :: N, I
      CHARACTER (LEN = 255) :: SUBSTRS(MAXDIM)
      CHARACTER (LEN = 255) :: ERRMSG

      !=================================================================
      ! READ_TIME_MENU BEGINS HERE!
      !=================================================================

      ! # OF DOMAIN
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_TIME_MENU: 1')
      READ(SUBSTRS(1:N), *) MAX_DOM

      ! START YEAR
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 2')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) START_YEAR(I)
      ENDDO

      ! START MONTH
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 3')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) START_MONTH(I)
      ENDDO

      ! START DAY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 4')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) START_DAY(I)
      ENDDO

      ! START HOUR
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 5')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) START_HOUR(I)
      ENDDO

      ! END YEAR
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 6')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) END_YEAR(I)
      ENDDO

      ! END MONTH
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 7')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) END_MONTH(I)
      ENDDO

      ! END DAY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 8')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) END_DAY(I)
      ENDDO

      ! END HOUR
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_TIME_MENU: 9')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) END_HOUR(I)
      ENDDO

      !SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_TIME_MENU: 10')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% TIME MENU %%%'
      WRITE(6, 110     ) '# OF DOMAIN  : ', MAX_DOM
      WRITE(6, 100     ) 'START YEAR   : ', START_YEAR(1: MAX_DOM)
      WRITE(6, 100     ) 'START MONTH  : ', START_MONTH(1: MAX_DOM)
      WRITE(6, 100     ) 'START DAY    : ', START_DAY(1: MAX_DOM)
      WRITE(6, 100     ) 'START HOUR   : ', START_HOUR(1: MAX_DOM)
      WRITE(6, 100     ) 'END YEAR     : ', END_YEAR(1: MAX_DOM)
      WRITE(6, 100     ) 'END MONTH    : ', END_MONTH(1: MAX_DOM)
      WRITE(6, 100     ) 'END DAY      : ', END_DAY(1: MAX_DOM)
      WRITE(6, 100     ) 'END HOUR     : ', END_HOUR(1: MAX_DOM)
100   FORMAT( A, 999I11)
110   FORMAT( A, I11   )

      !=================================================================
      ! CHECK THE INPUTS
      !=================================================================

      ! # OF DOMAINS
      IF (MAX_DOM .LE. 0) THEN
       WRITE(6, '(A)') 'Warning: MAX_DOM is less than 1'
       ERRMSG = 'Stopped, please check MAX_DOM in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

      DO I = 1, MAX_DOM
      ! GIVE A WARNING IF YEAR, MONTH, DAY, HOUR ARE OUT OF RANGE
       IF (START_YEAR(I)  .LE. 0 .OR. START_YEAR(I)  .GE. 9999) 
     &  WRITE(6, '(A)') 'Warning: START YEAR is out of range: 0-9999'

       IF (START_MONTH(I) .LT. 1 .OR. START_MONTH(I) .GT. 12  ) 
     &  WRITE(6, '(A)') 'Warning: START MONTH is out of range: 1-12'

       IF (START_DAY(I)   .LT. 1 .OR. START_DAY(I)   .GT. 31  ) 
     &  WRITE(6, '(A)') 'Warning: START DAY is out of range: 1-31'

       IF (START_HOUR(I)  .LT. 0 .OR. START_HOUR(I)  .GT. 24  ) 
     &  WRITE(6, '(A)') 'Warning: START HOUR is out of range: 0-24'

       IF (END_YEAR(I)    .LE. 0 .OR. END_YEAR(I)    .GE. 9999) 
     &  WRITE(6, '(A)') 'Warning: END YEAR is out of range: 0-9999'
      
       IF (END_MONTH(I)   .LT. 1 .OR. END_MONTH(I)   .GT.  12 ) 
     &  WRITE(6, '(A)') 'Warning: END MONTH is out of range: 1-12'

       IF (END_DAY(I)     .LT. 1 .OR. END_DAY(I)     .GT.  31 ) 
     &  WRITE(6, '(A)') 'Warning: END DAY is out of range: 1-31'

       IF (END_HOUR(I)    .LT. 0 .OR. END_HOUR(I)    .GT. 24  ) 
     &  WRITE(6, '(A)') 'Warning: END HOUR is out of range: 0-24'
      ENDDO

      END SUBROUTINE READ_TIME_MENU


!------------------------------------------------------------------------------
!
!  $ID: READ_DOMAIN_MENU V01 06/22/2012 10:34 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_DOMAIN_MENU READS DOMAIN MENU IN THE FILE
!  namelist.weps.
!
!  SUPPORTED PROJECTIONS
!  ---------------------
!  Cylindrical Lat/Lon (proj_code = PROJ_LATLON)
!  Mercator (proj_code = PROJ_MERC)
!  Lambert Conformal (proj_code = PROJ_LC)
!  Gaussian (proj_code = PROJ_GAUSS)
!  Polar Stereographic (proj_code = PROJ_PS)
!  Rotated Lat/Lon (proj_code = PROJ_ROTLL)
!
!  VARIABLES:
!  ============================================================================
!  (1 ) PROJ_CODE   (INTEGER) : MAP PROJECTION                           [---]
!  (2 ) DOM_ID      (INTEGER) : DOMAIN ORDER                             [---]
!  (3 ) DX          (REAL*8 ) : GRID SPACING IN METERS IN THE X-DIRECTION, USED 
!                               ONLY FOR PS, LC, AND MERC PROJECTIONS    [m  ]
!  (4 ) DY          (REAL*8)  : GRID SPACING IN METERS IN THE Y-DIRECTION, USED 
!                               ONLY FOR PS, LC, AND MERC PROJECTIONS    [m  ]
!  (5 ) E_WE        (INTEGER) : # OF GRID IN THE X-DIRECTION (FROM WPS)  [---]
!  (6 ) E_SN        (INTEGER) : # OF GRID IN THE Y-DIRECTION (FROM WPS)  [---]
!  (7 ) E_VERT      (INTEGER) : # OF VERTICAL LEVELS OF EMISSION         [---]
!  (8 ) CORNER_LAT  (REAL*8)  : LATITUDE OF SOUTHWEST CORNER             [deg]
!  (9 ) CORNER_LON  (REAL*8)  : LONGITUDE OF SOUTHWEST CORNER            [deg]
!  (10) STAND_LON   (REAL*8)  : STANDARD LONGITUDE (SEE WPS FOR DETAIL)  [deg]
!  (11) TRUELAT1    (REAL*8)  : FIRST TRUE LATITUDE (FOR ALL PROJECTIONS)[deg]
!  (12) TRUELAT2    (REAL*8)  : SECOND TRUE LATITUDE (FOR LC ONLY)       [deg]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM BY BRUCE. (06/22/2012)
!******************************************************************************
!
      SUBROUTINE READ_DOMAIN_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : MAX_DOM, PROJ_CODE
      USE NAMELIST_ARRAY_MOD, ONLY : DOM_ID, DX, DY
      USE NAMELIST_ARRAY_MOD, ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD, ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD, ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD, ONLY : TRUELAT1, TRUELAT2
      USE ERROR_MOD,          ONLY : ERROR_STOP


      ! LOCAL VARIABLES
      INTEGER               :: N, I
      CHARACTER (LEN = 255) :: SUBSTRS1, SUBSTRS(MAXDIM)
      CHARACTER (LEN = 255) :: ERRMSG

      !=================================================================
      ! READ_DOMAIN_MENU BEGINS HERE!
      !=================================================================

      ! PROJECTION CODE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_DOMAIN_MENU: 1 ')
      IF      (SUBSTRS(1) .EQ. 'PROJ_LATLON' ) THEN
       PROJ_CODE = 0
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_LC'     ) THEN
       PROJ_CODE = 1
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_PS'     ) THEN
       PROJ_CODE = 2
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_MERC'   ) THEN
       PROJ_CODE = 3
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_GAUSS'  ) THEN
       PROJ_CODE = 4
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_CYL'    ) THEN
       PROJ_CODE = 5
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_CASSINI') THEN
       PROJ_CODE = 6
      ELSE IF (SUBSTRS(1) .EQ. 'PROJ_ROTLL'  ) THEN
       PROJ_CODE = 203
      ENDIF

      ! DOMAIN ORDER
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 2 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) DOM_ID(I)
      ENDDO

      ! DX
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 3 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) DX(I)
      ENDDO

      ! DY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 4 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) DY(I)
      ENDDO

      ! E_WE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 5 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) E_WE(I)
      ENDDO

      ! E_SN
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 6 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) E_SN(I)
      ENDDO

      ! E_VERT
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_DOMAIN_MENU: 7 ')
      READ(SUBSTRS(1:N), *) E_VERT

      ! CORNER_LAT
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 8 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) CORNER_LAT(I)
      ENDDO

      ! CORNER_LON
      CALL SPLIT_ONE_LINE(SUBSTRS, N, MAX_DOM, 'READ_DOMAIN_MENU: 9 ')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) CORNER_LON(I)
      ENDDO

      ! STAND_LON
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_DOMAIN_MENU: 10')
      READ(SUBSTRS(1:N), *) STAND_LON

      ! TRUELAT1
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_DOMAIN_MENU: 11')
      READ(SUBSTRS(1:N), *) TRUELAT1

      ! TRUELAT2
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_DOMAIN_MENU: 12')
      READ(SUBSTRS(1:N), *) TRUELAT2

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1,       'READ_DOMAIN_MENU: 13')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% DOMAIN MENU %%%:'
      WRITE(6, 100     ) 'PROJ_CODE    : ', SUBSTRS(1)
      WRITE(6, 120     ) 'DOMAIN ORDER : ', DOM_ID(1: MAX_DOM)
      WRITE(6, 110     ) 'DX           : ', DX(1: MAX_DOM)
      WRITE(6, 110     ) 'DY           : ', DY(1: MAX_DOM)
      WRITE(6, 120     ) 'E_WE         : ', E_WE(1: MAX_DOM)
      WRITE(6, 120     ) 'E_SN         : ', E_SN(1: MAX_DOM)
      WRITE(6, 160     ) 'E_VERT       : ', E_VERT
      WRITE(6, 130     ) 'CORNER_LAT   : ', CORNER_LAT(1: MAX_DOM)
      WRITE(6, 130     ) 'CORNER_LON   : ', CORNER_LON(1: MAX_DOM)
      WRITE(6, 140     ) 'STAND_LON    : ', STAND_LON
      WRITE(6, 140     ) 'TRUELAT1     : ', TRUELAT1
      WRITE(6, 140     ) 'TRUELAT2     : ', TRUELAT2

100   FORMAT(A, 1X, A15 )
110   FORMAT(A, 999F11.1)
120   FORMAT(A, 999I11  )
130   FORMAT(A, 999F11.5)
140   FORMAT(A, F11.5   )
160   FORMAT(A, I11     )

      !=================================================================
      ! CHECK THE INPUTS
      !=================================================================

      !GIVE A WARNING IF THE INPUTS ARE OUT OF RANGE
      DO I = 1, MAX_DOM
       ! DOMAIN ORDER
       IF (DOM_ID(I) .NE. I) THEN
        WRITE(6, 150) 'Warning: DOM_ID of domain ', I, 
     &                ' does not match domain order'
        ERRMSG = 'Stopped, please check DOM_ID in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

       ! DX
       IF (DX(I) .LE. 0.0) THEN
        WRITE(6, 150) 'Warning: DX of domain ', I, 
     &                ' is less than 0.0'
        ERRMSG = 'Stopped, please check DX in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

       ! DY
       IF (DY(I) .LE. 0.0) THEN
        WRITE(6, 150) 'Warning: DY of domain ', I, 
     &                ' is less than 0.0'
        ERRMSG = 'Stopped, please check DY in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

       ! E_WE
       IF (E_WE(I) .LE. 0) THEN
        WRITE(6, 150) 'Warning: E_WE of domain ', I, 
     &                ' is less than 0'
        ERRMSG = 'Stopped, please check E_WE in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

       ! E_SN
       IF (E_SN(I) .LE. 0) THEN
        WRITE(6, 150) 'Warning: E_SN of domain ', I, 
     &                ' is less than 0'
        ERRMSG = 'Stopped, please check E_SN in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

       ! CORNER_LAT
       IF (CORNER_LAT(I) .LT. -90.0 .OR. CORNER_LAT(I) .GT. 90.0) THEN
        WRITE(6, 150) 'Warning: CORNER_LAT of domain ', I, 
     &                ' is out of range: -90.0-90.0 deg'
        ERRMSG = 'Stopped, please check CORNER_LAT in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

       ! CORNER_LON
       IF (CORNER_LON(I) .LT. -180.0 .OR. CORNER_LON(I) .GT. 180.0) THEN
        WRITE(6, 150) 'Warning: CORNER_LON of domain ', I, 
     &                ' is out of range: -180.0-180.0 deg'
        ERRMSG = 'Stopped, please check CORNER_LON in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

      ENDDO

      ! E_VERT
      IF (E_VERT .LT. 0) THEN
       WRITE(6, '(A)') 'Warning: E_VERT is less than 0'
       ERRMSG = 'Stopped, please check E_VERT in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

      ! STAND_LON
      IF (STAND_LON .LT. -180.0 .OR. STAND_LON .GT. 180.0) THEN
       WRITE(6, '(A)') 'Warning: STAND_LON is out of range: ' 
     &               //' -180.0-180.0 deg'
       ERRMSG = 'Stopped, please check STAND_LON in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

      ! TRUELAT1
      IF (TRUELAT1 .LT. -90.0 .OR. TRUELAT1 .GT. 90.0) THEN
       WRITE(6, '(A)') 'Warning: TRUELAT1 is out of range: ' 
     &               //' -90.0-90.0 deg'
       ERRMSG = 'Stopped, please check TRUELAT1 in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

      ! TRUELAT2
      IF (TRUELAT2 .LT. -90.0 .OR. TRUELAT2 .GT. 90.0) THEN
       WRITE(6, '(A)') 'Warning: TRUELAT2 is out of range: ' 
     &               //' -90.0-90.0 deg'
       ERRMSG = 'Stopped, please check TRUELAT2 in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

150   FORMAT(A, I2, A)

      END SUBROUTINE READ_DOMAIN_MENU


!------------------------------------------------------------------------------
!
!  $ID: READ_EMISSION_MENU V01 07/04/2012 09:58 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_EMISSION_MENU READS THE SWITCHES OF NEI OPTIONS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) LNEI              (LOGICAL) : SWITCH ON/OFF FOR NEI              [---]
!  (2 ) LNEIPOINT_DAILY   (LOGICAL) : SWITCH ON/OFF FOR NEI DAILY POINT  [---]
!  (3 ) LNEIPOINT_HOURLY  (LOGICAL) : SWITCH ON/OFF FOR NEI HOURLY POINT [---]
!  (4 ) LNEIAREA_DAILY    (LOGICAL) : SWITCH ON/OFF FOR NEI DAILY AREA   [---]
!  (5 ) LNEIAREA_HOURLY   (LOGICAL) : SWITCH ON/OFF FOR NEI HOURLY AREA  [---]
!  (6 ) LINTEX            (LOGICAL) : SWITCH ON/OFF FOR INTEX            [---]
!  (7 ) LFLAMBE           (LOGICAL) : SWITCH ON/OFF FOR FLAMBE           [---]
!  (8 ) LFINN             (LOGICAL) : SWITCH ON/OFF FOR FINN             [---]
!  (9 ) LGBBEP            (LOGICAL) : SWITCH ON/OFF FOR GBBEP            [---]
!  (10) LGFED             (LOGICAL) : SWITCH ON/OFF FOR GFED             [---]
!  (11) LSEVIRI           (LOGICAL) : SWITCH ON/OFF FOR SEVIRI           [---]
!  (12) LGFAS             (LOGICAL) : SWITCH ON/OFF FOR GFAS             [---]
!  (13) LQFED             (LOGICAL) : SWITCH ON/OFF FOR QFED             [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!  (2 ) ADD SEVIRI FIRE EMISSION LOGIC BY BRUCE. (07/17/2012)
!******************************************************************************
!
      SUBROUTINE READ_EMISSION_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : MAX_DOM
      USE NAMELIST_ARRAY_MOD, ONLY : LNEI
      USE NAMELIST_ARRAY_MOD, ONLY : LNEIPOINT_DAILY, LNEIPOINT_HOURLY
      USE NAMELIST_ARRAY_MOD, ONLY : LNEIAREA_DAILY,  LNEIAREA_HOURLY
      USE NAMELIST_ARRAY_MOD, ONLY : LINTEX
      USE NAMELIST_ARRAY_MOD, ONLY : LFLAMBE
      USE NAMELIST_ARRAY_MOD, ONLY : LFINN
      USE NAMELIST_ARRAY_MOD, ONLY : LGBBEP
      USE NAMELIST_ARRAY_MOD, ONLY : LGFED
      USE NAMELIST_ARRAY_MOD, ONLY : LSEVIRI
      USE NAMELIST_ARRAY_MOD, ONLY : LGFAS
      USE NAMELIST_ARRAY_MOD, ONLY : LQFED

      ! LOCAL VARIABLES
      INTEGER               :: N, I
      CHARACTER (LEN = 255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_EMISSION__MENU BEGINS HERE!
      !=================================================================

      ! LNEI
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 1')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LNEI(I)
      ENDDO

      ! LNEIPOINT_DAILY
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 2')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LNEIPOINT_DAILY(I)
      ENDDO

      ! LNEIPOINT_HOURLY
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 3')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LNEIPOINT_HOURLY(I)
      ENDDO

      ! LNEIAREA_DAILY
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 4')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LNEIAREA_DAILY(I)
      ENDDO

      ! LNEIAREA_HOURLY
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 5')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LNEIAREA_HOURLY(I)
      ENDDO

      ! LINTEX
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 6')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LINTEX(I)
      ENDDO

      ! LFLAMBE
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 7')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LFLAMBE(I)
      ENDDO

      ! LFINN
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 8')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LFINN(I)
      ENDDO

      ! LGBBEP
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU : 9')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LGBBEP(I)
      ENDDO

      ! LGFED
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU :10')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LGFED(I)
      ENDDO

      ! LSEVIRI
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU :11')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LSEVIRI(I)
      ENDDO

      ! LGFAS
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU :12')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LGFAS(I)
      ENDDO

      ! LQFED
      CALL SPLIT_ONE_LINE(SUBSTRS,N,MAX_DOM,'READ_EMISSION_MENU :13')
      DO I = 1, MAX_DOM
       READ(SUBSTRS(I), *) LQFED(I)
      ENDDO

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS,N,1,      'READ_EMISSION_MENU :14')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% EMISSION MENU %%%       : '
      WRITE(6, 110     ) 'TURN ON NEI                 : ',  
     &                    LNEI(1: MAX_DOM)
      WRITE(6, 110     ) ' - TURN ON NEI POINT DAILY  : ',  
     &                    LNEIPOINT_DAILY(1: MAX_DOM)
      WRITE(6, 110     ) ' - TURN ON NEI POINT HOURLY : ',  
     &                    LNEIPOINT_HOURLY(1: MAX_DOM)
      WRITE(6, 110     ) ' - TURN ON NEI AREA DAILY   : ',  
     &                    LNEIAREA_DAILY(1: MAX_DOM)
      WRITE(6, 110     ) ' - TURN ON NEI AREA HOURLY  : ',  
     &                    LNEIAREA_HOURLY(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON INTEX               : ',  
     &                    LINTEX(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON FLAMBE              : ',  
     &                    LFLAMBE(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON FINN                : ',  
     &                    LFINN(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON GBBEP               : ',  
     &                    LGBBEP(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON GFED                : ',  
     &                    LGFED(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON SEVIRI              : ',  
     &                    LSEVIRI(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON GFAS                : ',  
     &                    LGFAS(1: MAX_DOM)
      WRITE(6, 110     ) 'TURN ON QFED                : ',  
     &                    LQFED(1: MAX_DOM)

110   FORMAT(A, 999L5)

      END SUBROUTINE READ_EMISSION_MENU

!------------------------------------------------------------------------------
!
!  $ID: READ_INJECTION_HEIGHT_MENU V01 06/23/2012 11:15 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_INJECTION_HEIGHT_MENU READS INJECTION HEIGHT MENU FOR
!  SMOKE EMISSION.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) INJ_HEIGHT   (REAL*8) : INJECTION HEIGHT OF SMOKE EMISSION       [m  ]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL BY BRUCE. (06/23/2012)
!******************************************************************************
!
      SUBROUTINE READ_INJECTION_HEIGHT_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : INJ_HEIGHT
      USE ERROR_MOD,          ONLY : ERROR_STOP


      ! LOCAL VARIABLES
      INTEGER               :: N
      CHARACTER (LEN = 255) :: SUBSTRS(MAXDIM)
      CHARACTER (LEN = 255) :: ERRMSG

      !=================================================================
      ! READ_INJECTION_HEIGHT_MENU BEGINS HERE!
      !=================================================================

      ! INJECTION HEIGHT
      CALL SPLIT_ONE_LINE(SUBSTRS,N,1,'READ_INJECTION_HEIGHT_MENU: 1')
      READ(SUBSTRS(1:N), *) INJ_HEIGHT

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS,N,1,'READ_INJECTION_HEIGHT_MENU: 2')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% INJECTION HEIGHT MENU %%%: '
      WRITE(6, 100) 'INJ_HEIGHT   : ', INJ_HEIGHT
100   FORMAT(A, F10.2)

      !=================================================================
      ! CHECK THE INPUTS
      !=================================================================

      IF (INJ_HEIGHT .LT. 0.0) THEN
       WRITE(6, '(A)') 'Warning: INJ_HEIGHT is less than 0.0 m'
       ERRMSG = 'Stopped, please check INJ_HEIGHT in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

      END SUBROUTINE READ_INJECTION_HEIGHT_MENU

!------------------------------------------------------------------------------
!
!  $ID: READ_SMOKE_SPECIES_MENU V01 06/23/2012 11:55 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_SMOKE_SPECIES_MANU READS SPECIES WHICH SMOKE IS ADDED
!  IN THE BACKGROUND EMISSION.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) N_SPECIES    (INTEGER) : # OF SPECIES WHICH SMOKE IS ADDED IN
!                                THE BACKGROUND EMISSION                 [---]
!  (2 ) SPECIES      (INTEGER) : SPECIES' INDICES WHICH SMOKE IS ADDED 
!                                IN THE BACKGROUND EMISSION              [---]
!  (3 ) P_SPECIES    (REAL*8)  : PERCENTAGE OF SPECIES WHICH SMOKE IS
!                                ADDED IN THE BACKGROUND EMISSION        [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL BY BRUCE. (06/23/2012)
!******************************************************************************
!
      SUBROUTINE READ_SMOKE_SPECIES_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : N_SPECIES, SPECIES
      USE NAMELIST_ARRAY_MOD, ONLY : P_SPECIES
      USE NAMELIST_ARRAY_MOD, ONLY : NBASE_SPECIES, BASE_SPECIES
      USE ERROR_MOD,          ONLY : ERROR_STOP


      ! LOCAL VARIABLES
      INTEGER               :: N, I
      CHARACTER (LEN = 255) :: SUBSTRS(MAXDIM)
      CHARACTER (LEN = 255) :: ERRMSG

      !=================================================================
      ! READ_SMOKE_SPECIES_MENU BEGINS HERE!
      !=================================================================

      ! N_SPECIES
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'Read_SMOKE_SPECIES_MENU: 1')
      READ(SUBSTRS(1:N), *) N_SPECIES

      ! SPECIES
      CALL SPLIT_ONE_LINE(SUBSTRS, N, N_SPECIES, 
     &                    'Read_SMOKE_SPECIES_MENU: 2')
      DO I = 1, N_SPECIES
       READ(SUBSTRS(I), *) SPECIES(I)
      ENDDO

      ! % OF SPECIES
      CALL SPLIT_ONE_LINE(SUBSTRS, N, N_SPECIES, 
     &                    'Read_SMOKE_SPECIES_MENU: 3')
      DO I = 1, N_SPECIES
       READ(SUBSTRS(I), *) P_SPECIES(I)
      ENDDO

      ! NBASE_SPECIES
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'Read_SMOKE_SPECIES_MENU: 4')
      READ(SUBSTRS(1:N), *) NBASE_SPECIES

      ! BASE SPECIES FROM GFED
      CALL SPLIT_ONE_LINE(SUBSTRS, N, NBASE_SPECIES,
     &                    'Read_SMOKE_SPECIES_MENU: 5')
      DO I = 1, NBASE_SPECIES
         READ(SUBSTRS(1:N), *) BASE_SPECIES(I)
      ENDDO

      F_PM25 = .FALSE.
      IF (SPECIES(1) == 'PM2p5') THEN
         F_PM25 = .TRUE.
      ENDIF

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'Read_SMOKE_SPECIES_MENU: 6')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% SMOKE SPECIES MENU %%%:'
      WRITE(6, 110     ) 'N_SPECIES    : ', N_SPECIES
      WRITE(6, 120     ) 'SPECIES      : ', SPECIES(1:N_SPECIES)
      WRITE(6, 130     ) '% OF SPECIES : ', P_SPECIES(1:N_SPECIES)
      WRITE(6, 120     ) 'BASE SPECIES : ', BASE_SPECIES

110   FORMAT(A, I8)
120   FORMAT(A, 6X, 999A8)
130   FORMAT(A, 999F8.3)
      !=================================================================
      ! CHECK THE INPUTS
      !=================================================================

      ! N_SPECIES
      IF (N_SPECIES .LT. 0) THEN
       WRITE(6, '(A)') 'Warning: N_SPECIES is less than 0'
       ERRMSG = 'Stopped, please check N_SPECIES in namelist.weps'
       CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
      ENDIF

      DO I = 1, N_SPECIES
       ! % OF SPECIES
       IF ( P_SPECIES(I) .LT. 0 .OR. P_SPECIES(I) .GT. 100) THEN
        WRITE(6, 130) 'Warning: SPECIES ', I, 
     &                ' is out of range: '//'0-100'
        ERRMSG = 'Stopped, please check % OF SPECIES in namelist.weps'
        CALL ERROR_STOP(ERRMSG, 'NAMELIST_MOD.F')
       ENDIF

      ENDDO

      END SUBROUTINE READ_SMOKE_SPECIES_MENU

!------------------------------------------------------------------------------
!
!  $ID: READ_DIR_MENU V01 06/23/2012 15:21 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_DIR_MENU READS DIRECTORY OF EACH EMISSION, INCLUDING 
!  BACKGROUND AND FIRE EMISSIONS.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) NEIDIR     (CHARACTER) : NEI 2005 DIRECTORY                      [---]
!  (2 ) INTEXDIR   (CHARACTER) : INTEX-B DIRECTORY                       [---]
!  (3 ) FLAMBEDIR  (CHARACTER) : FLAMBE DIRECTORY                        [---]
!  (4 ) FINNDIR    (CHARACTER) : FINN DIRECTORY                          [---]
!  (5 ) GBBEPDIR   (CHARACTER) : GBBEP DIRECTORY                         [---]
!  (6 ) GFEDDIR    (CHARACTER) : GFED DIRECTORY                          [---]
!  (7 ) SEVIRIDIR  (CHARACTER) : SEVIRI DIRECTORY                        [---]
!  (8 ) GFASDIR    (CHARACTER) : GFAS DIRECTORY                          [---]
!  (9 ) QFEDDIR    (CHARACTER) : QFED DIRECTORY                          [---]
!  (10) OUTPUTDIR  (CHARACTER) : OUTPUT DIRECTORY                        [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL BY BRUCE. (06/23/2012)
!  (2 ) ADD SEVIRIDIR BY BRUCE. (07/16/2012)
!******************************************************************************
!
      SUBROUTINE READ_DIR_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : NEIDIR,    INTEXDIR, 
     &                               FLAMBEDIR, FINNDIR,  
     &                               GBBEPDIR,  GFEDDIR,  
     &                               SEVIRIDIR, GFASDIR,  
     &                               QFEDDIR,             
     &                               ERDIR,     OUTPUTDIR

      ! LOCAL VARIABLES
      INTEGER               :: N
      CHARACTER (LEN = 255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DIR_MENU BEGINS HERE!
      !=================================================================

      ! NEI DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 1')
      READ(SUBSTRS(1:N), '(A)') NEIDIR

      ! INTEX DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 2')
      READ(SUBSTRS(1:N), '(A)') INTEXDIR

      ! FLAMBE DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 3')
      READ(SUBSTRS(1:N), '(A)') FLAMBEDIR

      ! FINN DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 4')
      READ(SUBSTRS(1:N), '(A)') FINNDIR

      ! GBBEP DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 5')
      READ(SUBSTRS(1:N), '(A)') GBBEPDIR

      ! GFED DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 6')
      READ(SUBSTRS(1:N), '(A)') GFEDDIR

      ! SEVIRI DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 7')
      READ(SUBSTRS(1:N), '(A)') SEVIRIDIR

      ! GFAS DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 8')
      READ(SUBSTRS(1:N), '(A)') GFASDIR

      ! QFED DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 9')
      READ(SUBSTRS(1:N), '(A)') QFEDDIR

      ! ER DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 10')
      READ(SUBSTRS(1:N), '(A)') ERDIR

      ! OUTPUT DIRECTORY
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 11')
      READ(SUBSTRS(1:N), '(A)') OUTPUTDIR

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIR_MENU: 12')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% DIRECTORY MENU %%%: '
      WRITE(6, 120     ) 'NEIDIR       : ', TRIM(NEIDIR)
      WRITE(6, 120     ) 'INTEXDIR     : ', TRIM(INTEXDIR)
      WRITE(6, 120     ) 'FLAMBEDIR    : ', TRIM(FLAMBEDIR)
      WRITE(6, 120     ) 'FINNDIR      : ', TRIM(FINNDIR)
      WRITE(6, 120     ) 'GBBEPDIR     : ', TRIM(GBBEPDIR)
      WRITE(6, 120     ) 'GFEDDIR      : ', TRIM(GFEDDIR)
      WRITE(6, 120     ) 'SEVIRIDIR    : ', TRIM(SEVIRIDIR)
      WRITE(6, 120     ) 'GFASDIR      : ', TRIM(GFASDIR)
      WRITE(6, 120     ) 'QFEDDIR      : ', TRIM(QFEDDIR)
      WRITE(6, 120     ) 'OUTPUTDIR    : ', TRIM(OUTPUTDIR)
120   FORMAT(A, A)
      
      END SUBROUTINE READ_DIR_MENU

!------------------------------------------------------------------------------
!
!  $ID: READ_DIAGNOSTIC_MENU V01 06/23/2012 15:44 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_DIAGNOSTIC_MENU READS DIAGNOSTIC MENU IN THE FILE 
!  namelist.weps
!
!  VARIABLES:
!  ============================================================================
!  (1 ) LDIAG     (LOGICAL) : SWITCH ON/OFF FOR DIAGNOSTIC               [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RMT MODEL BY BRUCE. (06/23/2012)
!******************************************************************************
!
      SUBROUTINE READ_DIAGNOSTIC_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : LDIAG

      ! LOCAL VARIABLES
      INTEGER               :: N
      CHARACTER (LEN = 255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DIAGNOSTIC_MENU BEGINS HERE!
      !=================================================================

      ! LDIAG
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIAGNOSTIC_MENU: 1')
      READ(SUBSTRS(1:N), *) LDIAG

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DIAGNOSTIC_MENU: 2')

      ! SCREEN PRINT
      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% DIAGNOSTIC MENU %%%: '
      WRITE(6, 110     ) 'TURN ON DIAGNOSTIC?  : ', LDIAG
110   FORMAT(A, L5)

      END SUBROUTINE READ_DIAGNOSTIC_MENU

!------------------------------------------------------------------------------
!
!  $ID: READ_DEBUG_MENU V01 06/23/2012 15:57 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_DEBUG_MENU READS DEBUG SWITCHES IN THE DEBUG MENU OF
!  namelist.weps
!
!  VARIBALES:
!  ============================================================================
!  (1 ) LPRT          (LOGICAL) : SWITCH ON/OFF FOR SCREEN PRINTING      [---]
!  (2 ) LDEBUG_NEI    (LOGICAL) : SWITCH ON/OFF FOR NEI_MOD DEBUG        [---]
!  (3 ) LDEBUG_INTEX  (LOGICAL) : SWITCH ON/OFF FOR INTEX_MOD DEBUG      [---]
!  (4 ) LDEBUG_FLAMBE (LOGICAL) : SWITCH ON/OFF FOR FLAMBE_MOD DEBUG     [---]
!  (5 ) LDEBUG_FINN   (LOGICAL) : SWITCH ON/OFF FOR FINN_MOD_DEBUG       [---]
!  (6 ) LDEBUG_GBBEP  (LOGICAL) : SWITCH ON/OFF FOR GBBEP_MOD DEBUG      [---]
!  (7 ) LDEBUG_GFED   (LOGICAL) : SWITCH ON/OFF FOR GFED_MOD DEBUG       [---]
!  (8 ) LDEBUG_SEVIRI (LOGICAL) : SWITCH ON/OFF FOR SEVIRI_MOD DEBUG     [---]
!  (9 ) LDEBUG_GFAS   (LOGICAL) : SWITCH ON/OFF FOR GFAS_MOD DEBUG       [---]
!  (10) LDEBUG_QFED   (LOGICAL) : SWITCH ON/OFF FOR QFED_MOD DEBUG       [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM RICHARD'S ART-RTM MODEL BY BRUCE. (06/23/2012)
!  (2 ) ADD SEVIRI DEBUG LOGIC BY BRUCE. (07/17/2012)
!******************************************************************************
!
      SUBROUTINE READ_DEBUG_MENU

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : LPRT
      USE NAMELIST_ARRAY_MOD, ONLY : LDEBUG_NEI,    LDEBUG_INTEX
      USE NAMELIST_ARRAY_MOD, ONLY : LDEBUG_FLAMBE, LDEBUG_FINN
      USE NAMELIST_ARRAY_MOD, ONLY : LDEBUG_GBBEP,  LDEBUG_GFED
      USE NAMELIST_ARRAY_MOD, ONLY : LDEBUG_SEVIRI, LDEBUG_GFAS
      USE NAMELIST_ARRAY_MOD, ONLY : LDEBUG_QFED

      ! LOCAL VARIABLES
      INTEGER               :: N
      CHARACTER (LEN = 255) ::SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_DEBUG_MENU BEGINS HERE!
      !=================================================================

      ! LPRT
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  1')
      READ(SUBSTRS(1:N), *) LPRT

      ! LDEBUG_NEI
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  2')
      READ(SUBSTRS(1:N), *) LDEBUG_NEI

      ! LDEBUG_INTEX
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  3')
      READ(SUBSTRS(1:N), *) LDEBUG_INTEX

      ! LDEBUG_FLAMBE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  4')
      READ(SUBSTRS(1:N), *) LDEBUG_FLAMBE

      ! LDEBUG_FINN
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  5')
      READ(SUBSTRS(1:N), *) LDEBUG_FINN

      ! LDEBUG_GBBEP
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  6')
      READ(SUBSTRS(1:N), *) LDEBUG_GBBEP

      ! LDEBUG_GFED
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  7')
      READ(SUBSTRS(1:N), *) LDEBUG_GFED

      ! LDEBUG_SEVIRI
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  8')
      READ(SUBSTRS(1:N), *) LDEBUG_SEVIRI

      ! LDEBUG_GFAS
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU:  9')
      READ(SUBSTRS(1:N), *) LDEBUG_GFAS

      ! LDEBUG_QFED
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU: 10')
      READ(SUBSTRS(1:N), *) LDEBUG_QFED

      ! SEPARATOR LINE
      CALL SPLIT_ONE_LINE(SUBSTRS, N, 1, 'READ_DEBUG_MENU: 11')

      !=================================================================
      ! PRINT TO SCREEN
      !=================================================================

      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, '(   A)') '%%% DEBUG MENU %%%: '
      WRITE(6, 100     ) 'TURN ON SCREEN PRINT : ', LPRT
      WRITE(6, 100     ) 'TURN ON NEI DEBUG    : ', LDEBUG_NEI
      WRITE(6, 100     ) 'TURN ON INTEX DEBUG  : ', LDEBUG_INTEX
      WRITE(6, 100     ) 'TURN ON FLAMBE DEBUG : ', LDEBUG_FLAMBE
      WRITE(6, 100     ) 'TURN ON FINN DEBUG   : ', LDEBUG_FINN
      WRITE(6, 100     ) 'TURN ON GBBEP DEBUG  : ', LDEBUG_GBBEP
      WRITE(6, 100     ) 'TURN ON GFED DEBUG   : ', LDEBUG_GFED
      WRITE(6, 100     ) 'TURN ON SEVIRI DEBUG : ', LDEBUG_SEVIRI
      WRITE(6, 100     ) 'TURN ON GFAS DEBUG   : ', LDEBUG_GFAS
      WRITE(6, 100     ) 'TURN ON QFED DEBUG   : ', LDEBUG_QFED

100   FORMAT(A, L5)

      END SUBROUTINE READ_DEBUG_MENU

!------------------------------------------------------------------------------

      ! END OF MODULE
      END MODULE NAMELIST_MOD
