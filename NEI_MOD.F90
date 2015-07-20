! $ID: NEI_MOD.F V01 06/09/2012 23:05 BRUCE EXP$
!
!******************************************************************************
!  MODULE NEI_MOD CONTAINS PROCESSING SUBROUTINES FOR NEI 2005, WHICH
!  ARE AS FOLLOWS
!
!  MODULE VARIABLES:
!  ============================================================================
!  (1 )
!
!  MODULE ROUTINES:
!  ============================================================================
!  (1 ) NEI            : READS AND COMPUTE NEI 2005 EMISSION
!  (2 ) VERTICAL_LEVEL : DETERMINE THE VERTICAL LEVEL OF EMISSIONS
!  (3 ) INIT_NEI_ARRAY : INITIALIZES NEI ARRAYS
!  (4 ) CLEANUP_NEI    : CLEAN UP NEI ARRAYS
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/09/2012)
!******************************************************************************

      MODULE NEI_MOD

      ! USE OTHER MODULES
      USE MAP_UTILS


      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO) :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES 
      ! AND ROUTINES FROM BEING SEEN OUTSIDE "NEI_MOD.F"
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC   :: NEI3RD, ZFA
      PUBLIC   :: ENAME, NRADM, NHR

      ! ... EXCEPT THESE ROUTINES
      PUBLIC   :: NEI
      PUBLIC   :: VERTICAL_LEVEL
      PUBLIC   :: INIT_NEI_ARRAY
      PUBLIC   :: CLEANUP_NEI

      !=================================================================
      ! PARAMETERS
      !=================================================================

      ! NEI 2005 AREA EMISSION
      INTEGER, PARAMETER :: IXA      = 1332
      INTEGER, PARAMETER :: JXA      = 1008

      ! NEI 2005 POINT EMISSION
      INTEGER, PARAMETER :: IPOINT   = 168516

      ! NUMBER OF PRIMARY(IPRIM),VOC(IVOC),AND PM2.5(IPM25) SPECIES
      INTEGER, PARAMETER :: IPRIM    = 7
      INTEGER, PARAMETER :: IVOC     = 41
      INTEGER, PARAMETER :: IPM25    = 5

      ! NUMBER OF SPECIES
      INTEGER, PARAMETER :: NAL2DO   = 48
      INTEGER, PARAMETER :: NRADM    = 30
      INTEGER, PARAMETER :: NAMFILE2 = 56

      ! NUMBER OF GRID LEVELS AND INTERVAL LEVELS FOR WIND CLIMATOLOGY
      ! DATA
      ! USED IN VERTICAL MOMENTUM LIFT CALCULATIONS
      INTEGER, PARAMETER :: KWIN     = 20
      INTEGER, PARAMETER :: KWINP    = KWIN + 1

      !=================================================================
      ! NHR     : TOTAL HOURS PER DAY
      ! HRTOSED : SECONDS IN ONE HOUR
      !=================================================================
      INTEGER, PARAMETER :: NHR      = 24
      REAL*8,  PARAMETER :: HRTOSED  = 3600.

      !=================================================================
      ! CHARACTER
      !=================================================================
      CHARACTER (LEN = 9 ), DIMENSION(NAMFILE2) :: NAM, NAMRAD
      CHARACTER (LEN = 9 ), DIMENSION(NAL2DO)   :: NAM2EM
      CHARACTER (LEN = 9 ), DIMENSION(NRADM)    :: ENAME
      CHARACTER (LEN = 10), DIMENSION(IPRIM)    :: NAMINF
      CHARACTER (LEN = 10), DIMENSION(IVOC)     :: NAMVOC
      CHARACTER (LEN = 10), DIMENSION(IPM25)    :: NAMPM2
      CHARACTER (LEN = 10)                      :: CHSPEC, CHSCRT

      !=================================================================
      ! INTEGER
      !=================================================================
      INTEGER                                   :: I, J, N, IHR
      INTEGER, DIMENSION(NAMFILE2)              :: MWT

      !=================================================================
      ! REAL
      !=================================================================
      REAL*8, DIMENSION(IXA, JXA)               :: EMT
      REAL*8, DIMENSION(IPOINT)                 :: EMP
      REAL*8, DIMENSION(24, KWIN)               :: WSPD
      REAL*8, DIMENSION(NAMFILE2)               :: FAC
      REAL*8, DIMENSION(KWINP)                  :: REFWZ
      REAL*8, DIMENSION(KWIN)                   :: ZFA
      REAL, ALLOCATABLE                         :: NEI3RD(:,:,:,:,:)
      REAL*8                                    :: ZTOP
      REAL*8                                    :: ZDIF, FRC
      ! VARIABLES FOR POINT STACK EMISSIONS
      REAL*8, DIMENSION(IPOINT)                 :: STKHGT, STKDIAM
      REAL*8, DIMENSION(IPOINT)                 :: STKTMP
      REAL*8, DIMENSION(IPOINT)                 :: STKVEL, STKFLOW
      REAL*8                                    :: XLTMX, XLNMX
      REAL*8                                    :: SHGTMX

      !=================================================================
      ! ASSIGN VALUES TO VARIABLES
      !=================================================================
      
      ! INPUT EMISSION SPECIES
      DATA NAM2EM  &
         /'CO  ', 'NOX ', 'SO2 ', 'NH3 ', 'HC02', 'HC03', 'HC04',      &
          'HC05', 'HC06', 'HC07', 'HC08', 'HC09', 'HC10', 'HC12',      & 
          'HC13', 'HC14', 'HC15', 'HC16', 'HC17', 'HC18', 'HC19',      & 
          'HC20', 'HC21', 'HC22', 'HC23', 'HC24', 'HC25', 'HC26',      & 
          'HC27', 'HC28', 'HC29', 'HC31', 'HC32', 'HC33', 'HC34',      & 
          'HC35', 'HC36', 'HC37', 'HC38', 'HC39', 'HC40', 'HC41',      & 
          'PM01', 'PM02', 'PM03', 'PM04', 'PM05', 'PM10-PRI'/

      DATA NAMVOC      & 
         /' METHANE  ', ' ALKANE1  ', ' ALKANE2  ', ' ALKANE3  ',      & 
          ' ALKANE4  ', ' ALKANE5  ', ' ETHYLENE ', ' OLEFIN1  ',      & 
          ' OLEFIN2  ', ' ISOPRENE ', ' TERPENES ', ' AROMATIC1',      & 
          ' AROMATIC2', ' CH2O     ', ' CH3CHO   ', ' HI_ALDEHY',      & 
          ' BENZALDHY', ' ACETONE  ', ' MEK      ', ' PRD2     ',      & 
          ' MEOH     ', ' GLYOXAL  ', ' METHGLYOX', ' BACL     ',      & 
          ' PHENOLS  ', ' CRESOLS  ', ' MACR     ', ' MVK      ',      & 
          ' IPRD     ', ' UNREACTIV', ' PROPYLENE', ' ACETYLENE',      & 
          ' BENZENE  ', ' BUTANES  ', ' PENTANES ', ' TOLUENE  ',      & 
          ' XYLENES  ', ' PROPANE  ', ' DIENES   ', ' STYRENES ',      & 
          ' ORGACIDS '/

      DATA NAMPM2 /'PMFINE', 'PSO4', 'PNO3', 'POA', 'PEC'/

      DATA NAMINF      & 
         /'VOC       ', 'NOX       ', 'CO        ', 'SO2       ',      & 
          'PM10-PRI  ', 'PM25-PRI  ', 'NH3       '/

      ! NAMES OF RADM2 EMISSION VARIABLES
      ! RADM2: the Regional Acid Deposition Model (version2).
      ! including 14 stable species, 4 reactive intermediates and
      !           3 abundant stable species (oxygen, nitrogen and water)
      ! total species: 30 (see ename)
      DATA ENAME      & 
         /'e_so2 ', 'e_no  ', 'e_ald ', 'e_hcho', 'e_ora2',            & 
          'e_nh3 ', 'e_hc3 ', 'e_hc5 ', 'e_hc8 ', 'e_eth ',            & 
          'e_co  ', 'e_ol2 ', 'e_olt ', 'e_oli ', 'e_tol ',            & 
          'e_xyl ', 'e_ket ', 'e_csl ', 'e_iso ', 'e_pm25i',           & 
          'e_pm25j','e_so4i', 'e_so4j', 'e_no3i', 'e_no3j',            & 
          'e_orgi', 'e_orgj', 'e_eci ', 'e_ecj ', 'e_pm10'/

      DATA NAM      & 
         /'CO       ', 'NOX      ', 'SO2      ', 'NH3      ',          & 
          'HC02     ', 'HC03     ', 'HC04     ', 'HC05     ',          & 
          'HC06     ', 'HC07     ', 'HC08     ', 'HC09     ',          & 
          'HC10     ', 'HC12     ', 'HC13     ', 'HC14     ',          & 
          'HC15     ', 'HC16     ', 'HC17     ', 'HC18     ',          & 
          'HC19     ', 'HC20     ', 'HC21     ', 'HC22     ',          & 
          'HC23     ', 'HC24     ', 'HC25     ', 'HC26     ',          & 
          'HC27     ', 'HC27     ', 'HC28     ', 'HC28     ',          & 
          'HC29     ', 'HC31     ', 'HC32     ', 'HC33     ',          & 
          'HC34     ', 'HC35     ', 'HC36     ', 'HC37     ',          & 
          'HC38     ', 'HC39     ', 'HC40     ', 'HC40     ',          & 
          'HC41     ', 'PM01     ', 'PM01     ', 'PM02     ',          & 
          'PM02     ', 'PM03     ', 'PM03     ', 'PM04     ',          & 
          'PM04     ', 'PM05     ', 'PM05     ', 'PM10-PRI '/

      DATA NAMRAD      & 
         /'e_co     ', 'e_no     ', 'e_so2    ', 'e_nh3    ',          & 
          'e_eth    ', 'e_hc3    ', 'e_hc3    ', 'e_hc5    ',          & 
          'e_hc8    ', 'e_ol2    ', 'e_olt    ', 'e_oli    ',          & 
          'e_iso    ', 'e_tol    ', 'e_xyl    ', 'e_hcho   ',          & 
          'e_ald    ', 'e_ald    ', 'e_ald    ', 'e_ket    ',          & 
          'e_ket    ', 'e_ket    ', 'e_hc3    ', 'e_ald    ',          & 
          'e_ald    ', 'e_ald    ', 'e_csl    ', 'e_csl    ',          & 
          'e_ald    ', 'e_olt    ', 'e_ket    ', 'e_olt    ',          & 
          'e_ket    ', 'e_olt    ', 'e_hc3    ', 'e_tol    ',          & 
          'e_hc3    ', 'e_hc5    ', 'e_tol    ', 'e_xyl    ',          & 
          'e_hc3    ', 'e_oli    ', 'e_olt    ', 'e_tol    ',          & 
          'e_ora2   ', 'e_pm25i  ', 'e_pm25j  ', 'e_so4i   ',          & 
          'e_so4j   ', 'e_no3i   ', 'e_no3j   ', 'e_orgi   ',          & 
          'e_orgj   ', 'e_eci    ', 'e_ecj    ', 'e_pm10   '/

      ! ELEVATION AT GRID CELL TOP (m)
      DATA ZFA      & 
         /0., 50., 100., 160., 230., 300., 380., 470., 590.,           & 
          720., 850., 1000., 1200., 1400., 1600., 1800., 2000.,        & 
          2250., 2500., 2750./

      DATA REFWZ      & 
         /0., 16.8, 50.5, 84.3, 127., 170., 256., 343., 476.,          & 
          610., 791., 883., 975., 1160., 1350., 1644., 1943.,          & 
          2252., 2677., 3010., 3350./
      ! WIND SPEED
      DATA WSPD      & 
         /4.27, 4.01, 4.16, 4.27, 4.30, 4.26, 4.20, 4.12, 4.11,        & 
          4.11, 4.15, 4.25, 4.44, 3.43, 3.61, 3.86, 4.12, 4.36,        & 
          4.62, 4.94, 5.12, 5.21, 5.24, 5.12, 5.30, 5.62, 5.79,        & 
          5.96, 6.02, 5.96, 5.87, 5.74, 5.73, 5.70, 5.75, 5.85,        & 
          5.96, 4.51, 4.67, 4.93, 5.26, 5.58, 5.95, 6.40, 6.70,        & 
          6.88, 6.99, 6.91, 6.14, 6.82, 7.07, 7.27, 7.36, 7.28,        & 
          7.19, 7.02, 6.99, 6.92, 6.93, 7.02, 7.08, 5.46, 5.43,        & 
          5.64, 5.97, 6.31, 6.72, 7.24, 7.60, 7.88, 8.09, 8.12,        & 
          6.87, 7.73, 8.08, 8.36, 8.45, 8.37, 8.28, 8.09, 8.03,        & 
          7.92, 7.87, 7.91, 7.91, 6.34, 6.15, 6.22, 6.50, 6.82,        & 
          7.25, 7.80, 8.22, 8.58, 8.88, 9.04, 7.53, 8.46, 8.87,        & 
          9.22, 9.34, 9.26, 9.16, 8.95, 8.86, 8.70, 8.60, 8.57,        & 
          8.51, 6.94, 6.75, 6.74, 6.94, 7.21, 7.61, 8.16, 8.62,        & 
          9.05, 9.43, 9.73, 8.17, 9.20, 9.72,10.14,10.30,10.26,        & 
          10.16,9.92, 9.77, 9.55, 9.36, 9.25, 9.11, 7.42, 7.28,        & 
          7.38, 7.52, 7.65, 7.95, 8.44, 8.92, 9.42, 9.92,10.38,        & 
          8.55, 9.65,10.28,10.81,11.03,11.02,10.92,10.65,10.44,        & 
          10.17, 9.93, 9.77, 9.59, 7.72, 7.72, 7.95,8.29, 8.37,        & 
          8.35, 8.64, 9.07, 9.57,10.12,10.71, 8.74, 9.82,10.49,        & 
          11.08,11.36,11.39,11.32,11.05,10.83,10.54,10.25,10.12,       & 
          10.02, 7.89, 7.96, 8.21, 8.69, 8.89, 8.78, 8.88, 9.24,       & 
          9.71,10.26,10.89, 8.82, 9.84,10.50,11.11,11.41,11.43,        & 
          11.33,11.02, 10.75,10.44,10.18,10.14,10.25,7.92,8.09,        & 
          8.39, 8.93, 9.32, 9.34, 9.32, 9.59, 9.96,10.45,11.06,        & 
          8.86, 9.81,10.42,10.98, 11.20,11.11,10.92,10.54,10.23,       & 
          9.91, 9.71, 9.77,10.02, 7.95, 8.20, 8.62, 9.23, 9.72,        & 
          9.93, 9.92,10.03,10.31,10.70,11.19, 8.89, 9.74,10.28,        & 
          10.76,10.88,10.72,10.46,10.07, 9.80, 9.55, 9.39, 9.50,       & 
          9.81, 8.09, 8.36, 8.83, 9.48,10.02,10.33,10.34, 10.40,       & 
          10.60,10.90,11.26, 8.90, 9.67,10.14,10.54,10.61,10.42,       & 
          10.16, 9.79, 9.56, 9.34, 9.24, 9.40, 9.73, 8.19, 8.45,       & 
          8.94, 9.63, 10.21,10.56,10.58,10.65,10.79,11.02,11.31,       & 
          8.96, 9.56, 9.91,10.21,10.18, 9.96, 9.71, 9.39, 9.25,        & 
          9.12, 9.12, 9.33, 9.69, 8.38, 8.63, 9.17, 9.91,10.54,        & 
          10.91,11.01,11.03,11.12,11.22,11.40, 9.10, 9.50, 9.66,       & 
          9.84, 9.77, 9.58, 9.42, 9.23, 9.27, 9.31, 9.36, 9.53,        & 
          9.83, 8.71, 8.98, 9.53,10.30,10.93,11.33,11.47,11.45,        & 
          11.48,11.51,11.52, 9.45, 9.64, 9.60, 9.62, 9.54, 9.44,       & 
          9.42, 9.39, 9.58, 9.74, 9.79, 9.89,10.08, 9.24, 9.49,        & 
          9.96,10.70,11.36,11.77,11.96,11.95,11.90,11.83,11.68,        & 
          10.05,10.08, 9.90, 9.80, 9.72, 9.74, 9.90,10.02,10.28,       & 
          10.45,10.47,10.52,10.63,10.05,10.23,10.54,11.15,11.76,       & 
          12.14,12.28,12.32,12.30,12.21,11.97,10.66,10.66,10.50,       & 
          10.37,10.32,10.42,10.66,10.81,11.05,11.22,11.25,11.29,       & 
          11.35,10.99,11.08,11.20,11.62,12.10,12.37,12.44,12.52,       & 
          12.62,12.62,12.43,11.40,11.53,11.45,11.42,11.43,11.50,       & 
          11.68,11.77,11.98,12.21,12.33,12.40,12.46,12.17,12.18,       & 
          12.08,12.23,12.44,12.54,12.61,12.84,13.13,13.33,13.34,       & 
          12.21,12.47,12.43,12.47,12.51,12.53,12.63,12.70,12.95,       & 
          13.26,13.44,13.56,13.64,13.37,13.36,13.07,12.91,12.82,       & 
          12.71,12.74,13.07,13.52,13.93,14.19,12.76,13.18,13.21,       & 
          13.32,13.40,13.43,13.54,13.64,13.97,14.32,14.53,14.68,       & 
          14.72,14.37,14.37,14.05,13.78,13.56,13.31,13.17,13.42,       & 
          13.87,14.33,14.73/

      DATA FAC /       & 
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.11, 0.97, 1.00, 1.00,    & 
        1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.33,    & 
        1.61, 1.61, 0.40, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.50,    & 
        0.50, 0.50, 1.00, 1.00, 0.40, 0.29, 1.11, 0.97, 1.00, 1.00,    & 
        0.57, 1.00, 1.00, 1.00, 1.00, 0.20, 0.80, 0.20, 0.80, 0.20,    & 
        0.80, 0.20, 0.80, 0.20, 0.80, 1.00   /

      DATA MWT /     & 
        28,  46,  64,  17,  00,  00,  00,  00,  00,  00,  00,  00,     &
        00,  00,  00,  00,  00,  00,  00,  00,  00,  00,  00,  00,     &
        00,  00,  00,  00,  00,  00,  00,  00,  00,  00,  00,  00,     &
        00,  00,  00,  00,  00,  00,  00,  00,  00,  01,  01,  01,     &
        01,  01,  01,  01,  01,  01,  01,  01  /


      ! ... EXCEPT THESE VARIABLES ...

      CONTAINS

!
!  $ID: NEI_POINT V01 06/10/2012 15:08 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE NEI_POINT READS NEI 2005 EMISSION POINT EMISSION INFO
!  (INCLUDING LATITUDE, LONGITUDE, HEIGHT, AND SO FORTH), POINT EMISION.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) MODIFIED FROM NEI PREPROCESSOR (emiss_v03.F) BY BRUCE. (06/10/2012)
!******************************************************************************
!

      SUBROUTINE NEI(IDOM)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD,   ONLY : NEIDIR
      USE NAMELIST_ARRAY_MOD,   ONLY : LNEIPOINT_DAILY
      USE NAMELIST_ARRAY_MOD,   ONLY : LNEIPOINT_HOURLY
      USE NAMELIST_ARRAY_MOD,   ONLY : LNEIAREA_DAILY
      USE NAMELIST_ARRAY_MOD,   ONLY : LNEIAREA_HOURLY
      USE PARAMETER_MOD,        ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,        ONLY : LATINC, LONINC

      !=================================================================
      ! CHARACTER
      !=================================================================
      CHARACTER (LEN = 128)                    :: POINTDIR, AREADIR
      CHARACTER (LEN = 128)                    :: CHR
      CHARACTER (LEN = 15 ), DIMENSION(IPOINT) :: SITEID
      CHARACTER (LEN = 6  ), DIMENSION(IPOINT) :: ERPTID, UNITID
      CHARACTER (LEN = 6  ), DIMENSION(IPOINT) :: PROCID

      !=================================================================
      ! INTEGER
      !=================================================================
      INTEGER                                  :: IDOM, IMAX
      INTEGER                                  :: ISP, IPN
      INTEGER                                  :: IIXA, IJXA
      INTEGER                                  :: STATUS, SYSTEM
      INTEGER                                  :: KK, K, I2, J2, X1, Y1
      INTEGER                                  :: IDHMX, KSTAK, KSTKDHM
      INTEGER                                  :: IBOT, ITOP
      INTEGER, DIMENSION(IPOINT)               :: ISTATE, ICOUN, IRTYP

      !=================================================================
      ! REAL
      !=================================================================
      REAL, DIMENSION(NAL2DO, NRADM)           :: TRANAL2R
      ! VARIABLES FOR POINT STACK EMISSIONS
      REAL, DIMENSION(IPOINT)                  :: XLNP, XLTP
      REAL, DIMENSION(IXA, JXA)                :: XLAT, XLON
      REAL                                     :: X, Y

      ! CALL MAP_SET TO SPECIFY THE WRF USED MAP, INCLUDING GRID BOX
      CALL MAP_SET                     &                               
            ( PROJ_LC, PROJ, CORNER_LAT(IDOM), CORNER_LON(IDOM),      &
              CORNER_LAT(IDOM), CORNER_LON(IDOM), KNOWNI, KNOWNJ,     &
              DX(IDOM), DY(IDOM), LATINC, LONINC, STAND_LON,          &
              TRUELAT1, TRUELAT2)

      ! SPECIFY POINT AND AREA EMISSION DIRECTORY
      POINTDIR = NEIDIR(1: LEN_TRIM(NEIDIR))//'point/'
      AREADIR  = NEIDIR(1: LEN_TRIM(NEIDIR))//'area4k/'

      !=================================================================
      ! READ POINT EMISSION INFO
      !=================================================================
      OPEN(8, FILE = POINTDIR(1:LEN_TRIM(POINTDIR))//'Relpnt_info.txt')
      SHGTMX = -1.E20
      DO I = 1, IPOINT
! STKHGT: [m], injection source height, e.g., chimney
! STKDIAM: [m}, diameter of injection source, e.g., that of chimney
! STKTMP: [k], temperature of injection source
! STKVEL: [m/s], vertical velocity of injection source
! STKFLOW: fluxes of injection source
       READ(8, 55) ISTATE(I),  ICOUN(I),  SITEID(I), ERPTID(I),      &
                   UNITID(I),  PROCID(I), IRTYP(I),  STKHGT(I),      &
                   STKDIAM(I), STKTMP(I), STKVEL(I),                 &
                   STKFLOW(I), XLNP(I),   XLTP(I)
       IF (STKHGT(I) .GT. SHGTMX) THEN
        IMAX   = I
        SHGTMX = STKHGT(I)
       ENDIF
      ENDDO
55    FORMAT(I2.2, I3.3, A15, 3A6, 1X, I2, 1X, F10.3, F10.3, F10.2,  &
             F10.3, F10.4, F11.5, F10.5)
      CLOSE(8)

      !=================================================================
      ! READ AREA EMISSION INFO
      !=================================================================
      OPEN(8,FILE = NEIDIR(1:LEN_TRIM(NEIDIR))//'grid_loc/LAT_xrs.txt')
      READ(8,56) XLAT
      CLOSE(8)
      OPEN(8,FILE = NEIDIR(1:LEN_TRIM(NEIDIR))//'grid_loc/LON_xrs.txt')
      READ(8,56) XLON
      CLOSE(8)
56    FORMAT(12F10.5)

      ! READ IN TABLE FOR MAPPING AL EMISSION INTO RADM2 EMISSIONS
      TRANAL2R = 0.0
      DO KK = 1, NAMFILE2
       DO K = 1, NAL2DO
        IF (NAM(KK) .EQ. NAM2EM(K)) THEN
         DO J = 1, NRADM
          IF(NAMRAD(KK) .EQ. ENAME(J))THEN
           IF(MWT(KK).EQ.0)THEN
      ! Units of VOC in mole/hr --> mole/km(2)/hr
            TRANAL2R( K, J ) = FAC(KK)/(((1.E-3)*DX(IDOM))*((1.E-3)*DY(IDOM)))
           ELSEIF(MWT(KK).EQ.1)THEN
      ! Units of aerosol in ton/hr --> microgram/m(2)/sec
            TRANAL2R( K, J ) = FAC(KK)*9.07184E11 /DX(IDOM)/     &
                               DY(IDOM)/HRTOSED
           ELSE
      ! Units of primary species ton/hr --> mole/km(2)/hr
            TRANAL2R( K, J ) = FAC(KK)*9.07184E5/FLOAT(MWT(KK))  &
                               /(((1.E-3)*DX(IDOM))*((1.E-3)*DY(IDOM)))
           ENDIF
          ELSE
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDDO

      !=================================================================
      ! READ NEI 2005 POINT DAILY EMISSION
      !=================================================================

      IF (LNEIPOINT_DAILY(IDOM)) THEN
       WRITE(6, '(A49I2.2)')    &
       'Starting to Process NEI POINT DAILY for Domain : ', IDOM
       IHR = 0
       EMP = 0.
       DO ISP = 1, IPRIM
        ! UNZIP POINT DAILY EMISSION FILE
        STATUS = SYSTEM('rm -f scratem*')
        CHSPEC = NAMINF(ISP)
        STATUS = SYSTEM('cp '//POINTDIR(1:LEN_TRIM(POINTDIR))//'dayav/' &
                        //CHSPEC(1:LEN_TRIM(CHSPEC))//'.gz scratem.gz')
        STATUS = SYSTEM('gunzip scratem')
        OPEN(8, FILE = 'scratem', FORM = 'FORMATTED')
        READ(8, '(12E9.2)')EMP
        CLOSE(8)
        ! INTERPOLATE POINT DAILY EMISSION INTO WRF GRID BOX
        DO IPN = 1, IPOINT
         CALL LATLON_TO_IJ(PROJ, XLTP(IPN), XLNP(IPN), X, Y)
         IF(X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.    &
            Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
         IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
         IF(X-INT(X) .LE. 0.5) X1 = INT(X)
         IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
         IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
         ! CALL SUBROUTINE VERTICAL_LEVEL TO FIND K LEVEL OF 
         ! THE NEI 2005 POINT DAILY EMISSION
         CALL VERTICAL_LEVEL(IPN, IHR, IBOT, ITOP, ZDIF)
         DO I = IBOT, ITOP
          FRC = (ZFA(I+1)-ZFA(I))/ZDIF
          DO N = 1, NRADM
           NEI3RD(X1, I, Y1, N, :) = NEI3RD(X1, I, Y1, N, :) +    &
                                    EMP(IPN)*TRANAL2R(ISP, N)*FRC
          ENDDO ! N LOOP
         ENDDO ! I LOOP
         ENDIF
        ENDDO ! IPN LOOP
       ENDDO ! ISP LOOP
      ! END OF NEI DAILY POINT
      ENDIF

      !=================================================================
      ! READ NEI 2005 AREA DAILY EMISSION
      !=================================================================

      IF(LNEIAREA_DAILY(IDOM)) THEN
       WRITE(6, '(A48I2.2)')'Starting to Process NEI AREA DAILY for Doamin : ', IDOM
       IHR = 0
       EMT = 0.
       DO ISP = 1, IPRIM
        ! UNZIP AREA DAILY EMISSION FILE
        STATUS = SYSTEM('rm -f scratem*')
        CHSPEC = NAMINF(ISP)
        STATUS = SYSTEM('cp '//AREADIR(1:LEN_TRIM(AREADIR))//'dayav/'  &
                        //CHSPEC(1:LEN_TRIM(CHSPEC))//'.gz scratem.gz')
        STATUS = SYSTEM('gunzip scratem')
        OPEN(8, FILE = 'scratem', FORM = 'FORMATTED')
        READ(8, '(12E9.2)')EMT
        CLOSE(8)
        ! INTERPOLATE AREA DAILY EMISSION INTO WRF GRID BOX
        DO IIXA = 1, IXA
         DO IJXA = 1, JXA
          CALL LATLON_TO_IJ(PROJ, XLAT(IIXA, IJXA), XLON(IIXA, IJXA), X, Y)
          IF(X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.     &
             Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
           IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
           IF(X-INT(X) .LE. 0.5) X1 = INT(X)
           IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
           IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
           DO N = 1, NRADM
            NEI3RD(X1, 1, Y1, N, :) = NEI3RD(X1, 1, Y1, N, :) +      &
                                     EMT(IIXA, IJXA)*TRANAL2R(ISP, N)
           ENDDO ! N LOOP
          ENDIF
         ENDDO ! IJXA LOOP
        ENDDO ! IIXA LOOP
       ENDDO ! ISP LOOP
      ! END OF NEI DAILY AREA
      ENDIF


      !=================================================================
      ! READ NEI 2005 POINT HOURLY EMISSION
      !=================================================================

      IF (LNEIPOINT_HOURLY(IDOM)) THEN
       WRITE(6,'(A50I2.2)') 'Starting to Process NEI POINT HOURLY for Domain : ', IDOM
       DO IHR = START_HOUR(IDOM), END_HOUR(IDOM)
       EMP = 0.
       WRITE(CHR, '(I2.2)') IHR
        DO ISP = 1, NAL2DO
         ! UNZIP POINT HOURLY EMISSION FILE
         STATUS = SYSTEM('rm -f scratem*')
         CHSPEC = NAM2EM(ISP)
         STATUS = SYSTEM('cp '//POINTDIR(1:LEN_TRIM(POINTDIR))//'HR'   &
                         //CHR(1:LEN_TRIM(CHR))//'/'//                 &
                         CHSPEC(1:LEN_TRIM(CHSPEC))//'.gz scratem.gz')
         STATUS = SYSTEM('gunzip scratem')
         OPEN(8, FILE = 'scratem', FORM = 'FORMATTED')
         READ(8, '(12E9.2)')EMP
         CLOSE(8)
         ! INTERPOLATE POINT HOURLY EMISSION INTO WRF GRID BOX
         DO IPN = 1, IPOINT
          CALL LATLON_TO_IJ(PROJ, XLTP(IPN), XLNP(IPN), X, Y)
          IF(X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.                  &
             Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
          IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
          IF(X-INT(X) .LE. 0.5) X1 = INT(X)
          IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
          IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
          ! CALL SUBROUTINE CERTICAL_LEVEL TO FIND K LEVEL OF 
          ! THE NEI 2005 POINT HOURLY EMISSION
          CALL VERTICAL_LEVEL(IPN, IHR, IBOT, ITOP, ZDIF)
          DO I = IBOT, ITOP
           FRC = (ZFA(I+1)-ZFA(I))/ZDIF
           DO N = 1, NRADM
            NEI3RD(X1, I, Y1, N, IHR) = NEI3RD(X1, I, Y1, N, IHR) +   &
                                      EMP(IPN)*TRANAL2R(ISP, N)*FRC
           ENDDO ! N LOOP
          ENDDO ! I LOOP
          ENDIF
         ENDDO ! IPN LOOP
        ENDDO ! ISP LOOP
       ENDDO ! IHR LOOP
      ! END OF NEI HOURLY POINT
      ENDIF

      !=================================================================
      ! READ NEI 2005 AREA HOURLY EMISSION
      !=================================================================

      IF(LNEIAREA_HOURLY(IDOM)) THEN
       WRITE(6, '(A49I2.2)')'Starting to Process NEI AREA HOURLY for Domain : ', IDOM
       DO IHR = START_HOUR(IDOM), END_HOUR(IDOM)
        EMT = 0.
        WRITE(CHR, '(I2.2)') IHR
        DO ISP = 1, NAL2DO
         ! UNZIP AREA HOURLY EMISSION FILE
         STATUS = SYSTEM('rm -f scratem*')
         CHSPEC = NAM2EM(ISP)
         STATUS = SYSTEM('cp '//AREADIR(1:LEN_TRIM(AREADIR))//'HR'    &
                         //CHR(1:LEN_TRIM(CHR))//'/'//                &
                         CHSPEC(1:LEN_TRIM(CHSPEC))//'.gz scratem.gz')
         STATUS = SYSTEM('gunzip scratem')
         OPEN(8, FILE = 'scratem', FORM = 'FORMATTED')
         READ(8, '(12E9.2)')EMT
         CLOSE(8)
         ! INTERPOLATE AREA HOURLY EMISSION INTO WRF GRID BOX
         DO IIXA = 1, IXA
          DO IJXA = 1, JXA
           CALL LATLON_TO_IJ(PROJ, XLAT(IIXA, IJXA), XLON(IIXA, IJXA), X, Y)
           IF(X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.                &
              Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
            IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
            IF(X-INT(X) .LE. 0.5) X1 = INT(X)
            IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
            IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
            DO N = 1, NRADM
             NEI3RD(X1, 1, Y1, N, IHR) = NEI3RD(X1, 1, Y1, N, IHR) + &
                                      EMT(IIXA, IJXA)*TRANAL2R(ISP, N)
            ENDDO ! N LOOP
           ENDIF
          ENDDO ! IJXA LOOP
         ENDDO ! IIXA LOOP
        ENDDO ! ISP LOOP
       ENDDO ! IHR LOOP
      ! END OF NEI HOURLY AREA
      ENDIF


      END SUBROUTINE NEI


!------------------------------------------------------------------------------
!  $ID: INIT_NEI_ARRAY V01 07/03/2012 09:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_NEI_ARRAY INITIALIZES NEI RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) NEI3RD    (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NRADM, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/03/2012)
!******************************************************************************
!
      SUBROUTINE INIT_NEI_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! NEI OUTPUT VARIABELS
      !========================================================================

      ! NEI3RD
      ALLOCATE(NEI3RD(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NRADM, NHR), STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('NEI3RD')
      NEI3RD = 0D0

      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_NEI_ARRAY


!------------------------------------------------------------------------------
!  $ID: VERTICAL_LEVEL V01 06/10/2012 17:11 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE VERTICAL_LEVEL DETERMINES HEIGHT LEVEL OF EMISSION,
!  DEPENDING ON EMISSION HEIGHT.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/10/2012)
!******************************************************************************
!

      SUBROUTINE VERTICAL_LEVEL (IPN, IHR, IBOT, ITOP, ZDIF)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : E_VERT

      ! ARGUMENTS
      INTEGER                   :: IPN, IHR, IBOT, ITOP
      REAL*8                    :: ZDIF

      ! LOCAL VARIABLES
      INTEGER                   :: I, K, KSTAK
      REAL*8                    :: ZTOP, BOT, TOP, DHM

      DHM = 0.
      K   = 1
      ZTOP = REFWZ(K)
      DO WHILE ( ZTOP .LE. STKHGT(IPN))
       K = K + 1
       IF (K .GT. KWINP) THEN
        STOP
       ENDIF
       ZTOP = REFWZ(K)
      ENDDO
      KSTAK = K - 1
      IF (IHR .EQ. 0) THEN
       DHM   = 3.0*STKDIAM(IPN)*STKVEL(IPN)/        &
               (SUM(WSPD(:, KSTAK))/MAX(1, SIZE(WSPD(:, KSTAK))))
      ELSE
       DHM   = 3.0*STKDIAM(IPN)*STKVEL(IPN)/WSPD(IHR, KSTAK)
      ENDIF
      TOP   = STKHGT(IPN) + 1.5*DHM
      BOT   = STKHGT(IPN) + 0.5*DHM
      DO I = E_VERT,1,-1
       IF(ZFA(I+1).GT.BOT)THEN
        IBOT = I
       ENDIF
      ENDDO
      DO I = E_VERT,1,-1
       IF(ZFA(I+1).GT.TOP)THEN
        ITOP = I
      ENDIF
      ENDDO
      IF(IBOT.GE.ITOP)THEN
       ITOP=IBOT+1
      ENDIF
      ZDIF = ZFA(ITOP+1)-ZFA(IBOT)

      END SUBROUTINE VERTICAL_LEVEL

!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_NEI V01 06/18/2012 08:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_NEI DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) NEI3RD    (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NRADM, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT.
!       (06/18/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_NEI

      !=================================================================
      ! CLEANUP_NEI BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( NEI3RD   ) ) DEALLOCATE( NEI3RD      )

      END SUBROUTINE CLEANUP_NEI



      END MODULE NEI_MOD      
