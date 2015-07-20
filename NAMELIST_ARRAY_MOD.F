!  $ID: NAMELIST_ARRAY_MOD V01 06/23/2012 16:21 BRUCE EXP$
!
!******************************************************************************
!  MODULE NAMELIST_ARRAY_MOD CONTAINS ALL THE VARIABLES THAT READ FROM
!  namelist.weps BY NAMELIST_MOD.F
!
!  MODULE VARIABLES:
!  ============================================================================
!  TIME MENU:
!  (1 ) MAX_DOM     (INTEGER) : # OF DOMAINS                           [-----]
!  (2 ) START_YEAR  (INTEGER) : START YEAR OF EMISSION DATA            [year ]
!  (3 ) START_MONTH (INTEGER) : START MONTH OF EMISSION DATA           [month]
!  (4 ) START_DAY   (INTEGER) : START DAY OF EMISSION DATA             [day  ]
!  (5 ) START_HOUR  (INTEGER) : START HOUR OF EMISSION DATA            [hour ]
!  (6 ) END_YEAR    (INTEGER) : END YEAR OF EMISSION DATA              [year ]
!  (7 ) END_MONTH   (INTEGER) : END MONTH OF EMISSION DATA             [month]
!  (8 ) END_DAY     (INTEGER) : END DAY OF EMISSION DATA               [day  ]
!  (9 ) END_HOUR    (INTEGER) : END HOUR OF EMISSION DATA              [hour ]
!------------------------------------------------------------------------------
!  DOMAIN MENU:
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
!------------------------------------------------------------------------------
!  EMISSION MENU:
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
!------------------------------------------------------------------------------
!  INJECTION HEIGHT MENU:
!  (1 ) INJ_HEIGHT   (REAL*8)       : INJECTION HEIGHT OF SMOKE EMISSION [m  ]
!------------------------------------------------------------------------------
!  SMOKE SPECIES MENU:
!  (1 ) N_SPECIES    (INTEGER)      : # OF SPECIES WHICH SMOKE IS ADDED IN
!                                     THE BACKGROUND EMISSION            [---]
!  (2 ) SPECIES    (CHARACTER)      : SPECIES WHICH SMOKE IS ADDED 
!                                     IN THE BACKGROUND EMISSION         [---]
!  (3 ) P_SPECIES  (REAL*8)         : PERCENTAGE OF SPECIES WHICH SMOKE IS
!                                     ADDED IN THE BACKGROUND EMISSION   [---]
!  (4 ) NBASE_SPECIES(INTEGER)      : # OF BASE SPECIES FROM GFED        [---]
!  (5 ) BASE_SPECIES (CHARACTER)    : BASE SPECIES(PM2p5, TPM) FROM GFED [---]
!------------------------------------------------------------------------------
!  DIRECTORY MENU:
!  (1 ) NEIDIR     (CHARACTER) : NEI 2005 DIRECTORY                      [---]
!  (2 ) INTEXDIR   (CHARACTER) : INTEX-B DIRECTORY                       [---]
!  (3 ) FLAMBEDIR  (CHARACTER) : FLAMBE DIRECTORY                        [---]
!  (4 ) FINNDIR    (CHARACTER) : FINN DIRECTORY                          [---]
!  (5 ) GBBEPDIR   (CHARACTER) : GBBEP DIRECTORY                         [---]
!  (6 ) GFEDDIR    (CHARACTER) : GFED DIRECTORY                          [---]
!  (7 ) SEVIRIDIR  (CHARACTER) : SEVIRI DIRECTORY                        [---]
!  (8 ) GFASDIR    (CHARACTER) : GFAS DIRECTORY                          [---]
!  (9 ) QFEDDIR    (CHARACTER) : QFED DIRECTORY                          [---]
!  (10) ERDIR      (CHARACTER) : EMISSION RATIO DIRECTORY                [---]
!  (11) OUTPUTDIR  (CHARACTER) : OUTPUT DIRECTORY                        [---]
!------------------------------------------------------------------------------
!  DIAGNOSTIC MANU
!  (1 ) LDIAG         (LOGICAL) : SWITCH ON/OFF FOR DIAGNOSTIC           [---]
!------------------------------------------------------------------------------
!  DEBUG MENU:
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
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/23/2012)
!  (2 ) ADD SWITCHES FOR EMISSIONS BY BRUCE.  (07/10/2012)
!  (3 ) ADD BASE SPECIES TO CALCULATE EMISSION RATIO BY BRUCE. (01/13/2014)
!******************************************************************************
!
      MODULE NAMELIST_ARRAY_MOD

      ! MAKE EVERYTHING PUBLIC ...
      PUBLIC

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! DECLARE A MAXIMUM VALUE TO DEFINE ARRAY
      INTEGER, PARAMETER        :: MAX_N = 5
      ! DECLARE A MAXIMUM VALUE TO DEFINE SMOKE SPECIES ARRAY
      INTEGER, PARAMETER        :: MAX_S = 50

      ! TIME MENU
      INTEGER                   :: MAX_DOM
      INTEGER, DIMENSION(MAX_N) :: START_YEAR
      INTEGER, DIMENSION(MAX_N) :: START_MONTH
      INTEGER, DIMENSION(MAX_N) :: START_DAY
      INTEGER, DIMENSION(MAX_N) :: START_HOUR
      INTEGER, DIMENSION(MAX_N) :: END_YEAR
      INTEGER, DIMENSION(MAX_N) :: END_MONTH
      INTEGER, DIMENSION(MAX_N) :: END_DAY
      INTEGER, DIMENSION(MAX_N) :: END_HOUR

      ! DOMAIN MENU
      INTEGER                   :: PROJ_CODE
      INTEGER, DIMENSION(MAX_N) :: DOM_ID
      REAL,  DIMENSION(MAX_N)   :: DX
      REAL,  DIMENSION(MAX_N)   :: DY
      INTEGER, DIMENSION(MAX_N) :: E_WE
      INTEGER, DIMENSION(MAX_N) :: E_SN
      INTEGER                   :: E_VERT
      REAL,  DIMENSION(MAX_N)   :: CORNER_LAT
      REAL,  DIMENSION(MAX_N)   :: CORNER_LON
      REAL                      :: STAND_LON
      REAL                      :: TRUELAT1
      REAL                      :: TRUELAT2

      ! EMISSION MENU
      LOGICAL, DIMENSION(MAX_N) :: LNEI
      LOGICAL, DIMENSION(MAX_N) :: LNEIPOINT_DAILY
      LOGICAL, DIMENSION(MAX_N) :: LNEIPOINT_HOURLY
      LOGICAL, DIMENSION(MAX_N) :: LNEIAREA_DAILY
      LOGICAL, DIMENSION(MAX_N) :: LNEIAREA_HOURLY
      LOGICAL, DIMENSION(MAX_N) :: LINTEX
      LOGICAL, DIMENSION(MAX_N) :: LFLAMBE
      LOGICAL, DIMENSION(MAX_N) :: LFINN
      LOGICAL, DIMENSION(MAX_N) :: LGBBEP
      LOGICAL, DIMENSION(MAX_N) :: LGFED
      LOGICAL, DIMENSION(MAX_N) :: LSEVIRI
      LOGICAL, DIMENSION(MAX_N) :: LGFAS
      LOGICAL, DIMENSION(MAX_N) :: LQFED

      ! INJECTION HEIGHT MENU
      REAL*8                    :: INJ_HEIGHT

      ! SMOKE SPECIES MENU
      INTEGER                   :: N_SPECIES
      CHARACTER (LEN = 120), DIMENSION(MAX_S) :: SPECIES
      REAL*8,  DIMENSION(MAX_S) :: P_SPECIES
      INTEGER                   :: NBASE_SPECIES
      CHARACTER (LEN = 120), DIMENSION(MAX_S) :: BASE_SPECIES

      ! DIRECTORY MENU
      CHARACTER (LEN = 255)     :: NEIDIR
      CHARACTER (LEN = 255)     :: INTEXDIR
      CHARACTER (LEN = 255)     :: FLAMBEDIR
      CHARACTER (LEN = 255)     :: FINNDIR
      CHARACTER (LEN = 255)     :: GBBEPDIR
      CHARACTER (LEN = 255)     :: GFEDDIR
      CHARACTER (LEN = 255)     :: SEVIRIDIR
      CHARACTER (LEN = 255)     :: GFASDIR
      CHARACTER (LEN = 255)     :: QFEDDIR
      CHARACTER (LEN = 255)     :: ERDIR
      CHARACTER (LEN = 255)     :: OUTPUTDIR

      ! DIAGNOSTIC MENU
      LOGICAL                   :: LDIAG

      ! DEBUG MENU
      LOGICAL                   :: LPRT
      LOGICAL                   :: LDEBUG_NEI
      LOGICAL                   :: LDEBUG_INTEX
      LOGICAL                   :: LDEBUG_FLAMBE
      LOGICAL                   :: LDEBUG_FINN
      LOGICAL                   :: LDEBUG_GBBEP
      LOGICAL                   :: LDEBUG_GFED
      LOGICAL                   :: LDEBUG_SEVIRI
      LOGICAL                   :: LDEBUG_GFAS
      LOGICAL                   :: LDEBUG_QFED

      ! END OF MODULE
      END MODULE NAMELIST_ARRAY_MOD
