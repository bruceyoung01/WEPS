!  $ID: INTEX_MOD.F V01 07/04/2012 15:42 BRUCE EXP$
!
!******************************************************************************
!  MODULE INTEX_MOD READS AND COMPUTES INTEX-B EMISSION AND REWRITES IT
!  INTO THE FORMAT WHICH CAN BE READ BY convert_emiss.exe IN WRF-Chem.
!
!  MODULE VARIABLES:
!  ============================================================================
!  (1 )
!
!  MODULE ROUTINES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      MODULE INTEX_MOD


      ! USE OTHER MODULES
      USE MAP_UTILS
      USE NEI_MOD,     ONLY : ENAME, NRADM, NHR


      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO) :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES 
      ! AND ROUTINES FROM BEING SEEN OUTSIDE "INTEX_MOD.F"
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC   :: INTEX3RD

      ! ... EXCEPT THESE ROUTINES
      PUBLIC   :: INTEX
      PUBLIC   :: INIT_INTEX_ARRAY
      PUBLIC   :: CLEANUP_INTEX

      !=================================================================
      ! PARAMETERS
      !=================================================================



      REAL, ALLOCATABLE        :: INTEX3RD(:, :, :, :, :)

      ! VARIABLES FOR INTEX EMISSION
      INTEGER                  :: TMPROW, TMPCOL
      REAL                     :: TMPLAT, TMPLON
      REAL                     :: TMPBC, TMPCO, TMPNOX, TMPOC
      REAL                     :: TMPPM10, TMPPM25, TMPSO2,  TMPVOC
      REAL                     :: TMPACET, TMPALK1, TMPALK2, TMPALK3
      REAL                     :: TMPALK4, TMPALK5, TMPARO1, TMPARO2
      REAL                     :: TMPBACL, TMPBALD, TMPCCHO, TMPCRES
      REAL                     :: TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD
      REAL                     :: TMPISOP, TMPMACR, TMPMEK,  TMPMEOH
      REAL                     :: TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL
      REAL                     :: TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2
      REAL                     :: TMPRCHO, TMPTERP


      ! ... EXCEPT THESE VARIABLES ...

      CONTAINS


!------------------------------------------------------------------------------
!
!  $ID: INTEX V01 08/15/2012 15:18 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INTEX READS INTEX EMISSION AND REWRITES IT INTO WRF-CHEM
!  FORMAT.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (08/15/2012)
!******************************************************************************
!

      SUBROUTINE INTEX (IDOM)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD,   ONLY : INTEXDIR
      USE PARAMETER_MOD,        ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,        ONLY : LATINC, LONINC
      USE GEOMETRY_MOD,         ONLY : DEG_TO_DIS

      !=================================================================
      ! ARGUMENTS
      !=================================================================
      INTEGER                  :: IDOM

      !=================================================================
      ! PARAMETER
      !=================================================================
      INTEGER, PARAMETER       :: NLAT      = 134
      INTEGER, PARAMETER       :: NLON      = 196
      INTEGER, PARAMETER       :: NNAME     = 8
      INTEGER, PARAMETER       :: NCAT      = 4
      INTEGER, PARAMETER       :: NVOCNAME  = 30
      INTEGER, PARAMETER       :: NVOCCAT   = 6
      REAL,    PARAMETER       :: YEARTOHR  = 8760.
      REAL,    PARAMETER       :: HRTOSED   = 3600.
      REAL,    PARAMETER       :: YEARTOSEC = 31536000.
      REAL,    PARAMETER       :: TONTOUG   = 10.0E12
      REAL,    PARAMETER       :: DIS       = 0.5

      !=================================================================
      ! CHARACTER
      !=================================================================
      CHARACTER (LEN = 255)    :: DIR, VOCDIR, SUBDIR, VOCSUBDIR
      CHARACTER (LEN = 255)    :: HEADER_IND, HEADER_POW
      CHARACTER (LEN = 255)    :: HEADER_RES, HEADER_TRA
      CHARACTER (LEN = 255)    :: FIND, FPOW, FRES, FTRA
      CHARACTER (LEN = 255)    :: INPTF_IND, INPTF_POW
      CHARACTER (LEN = 255)    :: INPTF_RES, INPTF_TRA
      CHARACTER (LEN = 255)    :: HEADERVOC_DOB, HEADERVOC_DOF
      CHARACTER (LEN = 255)    :: HEADERVOC_DOP, HEADERVOC_IND
      CHARACTER (LEN = 255)    :: HEADERVOC_POW, HEADERVOC_TRA
      CHARACTER (LEN = 255)    :: FVOC_DOB, FVOC_DOF, FVOC_DOP
      CHARACTER (LEN = 255)    :: FVOC_IND, FVOC_POW, FVOC_TRA
      CHARACTER (LEN = 255)    :: INPTFVOC_DOB, INPTFVOC_DOF
      CHARACTER (LEN = 255)    :: INPTFVOC_DOP, INPTFVOC_IND
      CHARACTER (LEN = 255)    :: INPTFVOC_POW, INPTFVOC_TRA
      CHARACTER (LEN = 3  ), DIMENSION(NNAME)    :: EINTEXNAME
      CHARACTER (LEN = 4  ), DIMENSION(NVOCNAME) :: EVOCINTEXNAME

      !=================================================================
      ! INTEGER
      !=================================================================
      INTEGER                  :: I, J, K, L, M
      INTEGER                  :: X1, Y1
      INTEGER, DIMENSION(NNAME-1)    :: RADM2

      !=================================================================
      ! REAL
      !=================================================================
      REAL                     :: X,  Y, REALI, REALJ
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1) :: WLAT, WLON
      REAL, DIMENSION(NLAT, NLON) :: INTEXLAT, INTEXLON
      REAL, DIMENSION(NLAT, NLON) :: INTEXLAT_DIS, INTEXLON_DIS
      REAL, DIMENSION(NLAT, NLON, NNAME, NCAT)       :: EMINTEX
      REAL, DIMENSION(NLAT, NLON, NVOCNAME, NVOCCAT) :: EMVOCINTEX

      ! CALL MAP_SET TO SPECIFY THE WRF USED MAP, INCLUDING GRID BOX
      CALL MAP_SET   &
            ( PROJ_LC, PROJ, CORNER_LAT(IDOM), CORNER_LON(IDOM),   &
              CORNER_LAT(IDOM), CORNER_LON(IDOM), KNOWNI, KNOWNJ,  &
              DX(IDOM), DY(IDOM), LATINC, LONINC, STAND_LON,       &
              TRUELAT1, TRUELAT2)

      ! CALL SUBROUTINE IJ_TO_LATLON TO GET THE LATITUDE AND LONGITUDE OF 
      ! EACH GRID BOX IN WRF-CHEM
      DO I = 1, E_SN(IDOM)-1
       DO J = 1, E_WE(IDOM)-1
        REALI = REAL(I)
        REALJ = REAL(J)
        CALL IJ_TO_LATLON(PROJ, REALJ, REALI, WLAT(J, I), WLON(J, I))
       ENDDO
      ENDDO

      ! SPECIFY DIRECTORY FOR INTEX EMISSION, INCLUDING NORMAL EMISSION
      ! AND VOC
      SUBDIR   = '/EMISSIONS/STREETS.DAVID/2006_Asia_Inventory_updated/'
      VOCSUBDIR= '/EMISSIONS/STREETS.DAVID/2006_Asia_Inventory_NMVOC/'

      ! SPECIFY FILE NAME FOR EACH KIND OF EMISSION
      FIND     = '2006_Asia_industry.csv'
      FPOW     = '2006_Asia_power.csv'
      FRES     = '2006_Asia_residential.csv'
      FTRA     = '2006_Asia_transportation.csv'
      FVOC_DOB = '2006_Asia_VOC_dob.csv'
      FVOC_DOF = '2006_Asia_VOC_dof.csv'
      FVOC_DOP = '2006_Asia_VOC_dop.csv'
      FVOC_IND = '2006_Asia_VOC_ind.csv'
      FVOC_POW = '2006_Asia_VOC_pow.csv'
      FVOC_TRA = '2006_Asia_VOC_tra.csv'

      ! SPECIFY THE ORDER OF RADM SPECIES INTO RADM2
      RADM2 = (/29, 11, 2, 27, 30, 21, 1/)
      ! SPECIES NAME FOR NAMAL EMISSION
      EINTEXNAME = [CHARACTER (LEN = 3) ::       &
                    'BC',   'CO',   'NOX', 'OC', &
                    'PM10', 'PM25', 'SO2', 'VOC']

      ! SPECIES NAME FOR VOC EMISSION
      EVOCINTEXNAME = [CHARACTER (LEN = 4) ::                           &
                        'ACET', 'ALK1', 'ALK2', 'ALK3', 'ALK4', 'ALK5', &
                        'ARO1', 'ARO2', 'BACL', 'BALD', 'CCHO', 'CRES', &
                        'ETHE', 'GLY',  'HCHO', 'IPRD', 'ISOP', 'MACR', &
                        'MEK',  'MEOH', 'MGLY', 'MVK',  'NROG', 'NVOL', &
                        'OLE1', 'OLE2', 'PHEN', 'PRD2', 'RCHO', 'TERP']

      WRITE(6, '(/, A)') REPEAT('-', 79)
      WRITE(6, *)'PROCESSING INTEX EMISSION'
      WRITE(6, '(/, A)') REPEAT('-', 79)

      ! INITIALIZE INTEX3RD
      INTEX3RD = 0.0

      !=================================================================
      ! READ INTEX EMISSION
      !=================================================================

      ! OPEN NORMAL EMISSION
      DIR      = INTEXDIR(1: LEN_TRIM(INTEXDIR))//SUBDIR(1:LEN_TRIM(SUBDIR))
      ! INDUSTRY EMISSION
      INPTF_IND = DIR(1: LEN_TRIM(DIR))//FIND(1: LEN_TRIM(FIND))
      OPEN(1, FILE = INPTF_IND, STATUS = 'OLD')
      READ(1, *)HEADER_IND
      ! POWER EMISSION
      INPTF_POW = DIR(1: LEN_TRIM(DIR))//FPOW(1: LEN_TRIM(FPOW))
      OPEN(2, FILE = INPTF_POW, STATUS = 'OLD')
      READ(2, *)HEADER_POW
      ! RESIDENTIAL EMISSION
      INPTF_RES = DIR(1: LEN_TRIM(DIR))//FRES(1: LEN_TRIM(FRES))
      OPEN(3, FILE = INPTF_RES, STATUS = 'OLD')
      READ(3, *)HEADER_RES
      ! TRANSPORTATION EMISSION
      INPTF_TRA = DIR(1: LEN_TRIM(DIR))//FTRA(1: LEN_TRIM(FTRA))
      OPEN(4, FILE = INPTF_TRA, STATUS = 'OLD')
      READ(4, *)HEADER_TRA

      ! OPEN VOC EMISSION
      VOCDIR   = INTEXDIR(1: LEN_TRIM(INTEXDIR))//VOCSUBDIR(1:LEN_TRIM(VOCSUBDIR))
      ! RESIDENTIAL BIOFUEL (DOB) VOC EMISSION
      INPTFVOC_DOB = VOCDIR(1: LEN_TRIM(VOCDIR))//FVOC_DOB(1: LEN_TRIM(FVOC_DOB))
      OPEN(11, FILE = INPTFVOC_DOB, STATUS = 'OLD')
      READ(11, *)HEADERVOC_DOB
      ! RESIDENTIAL FOSSIL FUEL (DOF) VOC EMISSION
      INPTFVOC_DOF = VOCDIR(1: LEN_TRIM(VOCDIR))//FVOC_DOF(1: LEN_TRIM(FVOC_DOF))
      OPEN(12, FILE = INPTFVOC_DOF, STATUS = 'OLD')
      READ(12, *)HEADERVOC_DOF
      ! RESIDENTIAL NON-COMBUSTION (DOP) VOC EMISSION
      INPTFVOC_DOP = VOCDIR(1: LEN_TRIM(VOCDIR))//FVOC_DOP(1: LEN_TRIM(FVOC_DOP))
      OPEN(13, FILE = INPTFVOC_DOP, STATUS = 'OLD')
      READ(13, *)HEADERVOC_DOP
      ! INDUSTRY VOC EMISSION
      INPTFVOC_IND = VOCDIR(1: LEN_TRIM(VOCDIR))//FVOC_IND(1: LEN_TRIM(FVOC_IND))
      OPEN(14, FILE = INPTFVOC_IND, STATUS = 'OLD')
      READ(14, *)HEADERVOC_IND
      ! POWER VOC EMISSION
      INPTFVOC_POW = VOCDIR(1: LEN_TRIM(VOCDIR))//FVOC_POW(1: LEN_TRIM(FVOC_POW))
      OPEN(15, FILE = INPTFVOC_POW, STATUS = 'OLD')
      READ(15, *)HEADERVOC_POW
      ! TRANSPORTATION VOC EMISSION
      INPTFVOC_TRA = VOCDIR(1: LEN_TRIM(VOCDIR))//FVOC_TRA(1: LEN_TRIM(FVOC_TRA))
      OPEN(16, FILE = INPTFVOC_TRA, STATUS = 'OLD')
      READ(16, *)HEADERVOC_TRA


      DO I = 1, NLAT
       DO J = 1, NLON

        !===============================================================
        ! READ NORMAL EMISSION
        !===============================================================

        ! CALL SUBROUTINE INIT_INTEX_TMP TO INITIALIZE TEMPORARY
        ! VARIABLES
        CALL INIT_INTEX_TMP
        ! READ INDUSTRY EMISSION
        READ(1, *, END = 100)TMPROW, TMPCOL,  TMPLAT, TMPLON, &
                             TMPBC,  TMPCO,   TMPNOX, TMPOC,  &
                             TMPPM10,TMPPM25, TMPSO2, TMPVOC
        INTEXLAT(I, J) = TMPLAT
        INTEXLON(I, J) = TMPLON
        ! CALCULATE AREA IN EACH GRID BOX
        INTEXLAT_DIS(I, J) = DEG_TO_DIS(                               &
                             (INTEXLAT(I, J) - DIS/2), INTEXLON(I, J), &
                             (INTEXLAT(I, J) + DIS/2), INTEXLON(I, J))
        INTEXLON_DIS(I, J) = DEG_TO_DIS(                               &
                             INTEXLAT(I, J), INTEXLON(I, J) - DIS/2,   &
                             INTEXLAT(I, J), INTEXLON(I, J) + DIS/2)
        EMINTEX(I, J, :, 1) = (/TMPBC,  TMPCO,   TMPNOX, TMPOC,        &
                                TMPPM10,TMPPM25, TMPSO2, TMPVOC/)
        ! CONVERT UNIT TON/YEAR/GRID => UG/SEC/M2
        EMINTEX(I, J, :, 1) = EMINTEX(I, J, :, 1)*TONTOUG/(YEARTOSEC*  &
                              (INTEXLAT_DIS(I, J)*INTEXLON_DIS(I, J)))
        CALL INIT_INTEX_TMP
        ! READ POWER EMISSION
        READ(2, *, END = 200)TMPROW, TMPCOL,  TMPLAT, TMPLON,   &
                             TMPBC,  TMPCO,   TMPNOX, TMPOC,    &
                             TMPPM10,TMPPM25, TMPSO2, TMPVOC
        EMINTEX(I, J, :, 2) = (/TMPBC,  TMPCO,   TMPNOX, TMPOC, &
                                TMPPM10,TMPPM25, TMPSO2, TMPVOC/)
        ! CONVERT UNIT TON/YEAR/GRID => UG/SEC/M2
        EMINTEX(I, J, :, 2) = EMINTEX(I, J, :, 2)*TONTOUG/(YEARTOSEC*  &
                              (INTEXLAT_DIS(I, J)*INTEXLON_DIS(I, J)))

        CALL INIT_INTEX_TMP
        ! READ RESIDENTIAL EMISSION
        READ(3, *, END = 300)TMPROW, TMPCOL,  TMPLAT, TMPLON,   &
                             TMPBC,  TMPCO,   TMPNOX, TMPOC,    &
                             TMPPM10,TMPPM25, TMPSO2, TMPVOC
        EMINTEX(I, J, :, 3) = (/TMPBC,  TMPCO,   TMPNOX, TMPOC, & 
                                TMPPM10,TMPPM25, TMPSO2, TMPVOC/)
        ! CONVERT UNIT TON/YEAR/GRID => UG/SEC/M2
        EMINTEX(I, J, :, 3) = EMINTEX(I, J, :, 3)*TONTOUG/(YEARTOSEC*  & 
                              (INTEXLAT_DIS(I, J)*INTEXLON_DIS(I, J)))

        CALL INIT_INTEX_TMP
        ! READ TRANSPORTATION EMISSION
        READ(4, *, END = 400)TMPROW, TMPCOL,  TMPLAT, TMPLON,   &
                             TMPBC,  TMPCO,   TMPNOX, TMPOC,    &
                             TMPPM10,TMPPM25, TMPSO2, TMPVOC
        EMINTEX(I, J, :, 4) = (/TMPBC,  TMPCO,   TMPNOX, TMPOC, &
                                TMPPM10,TMPPM25, TMPSO2, TMPVOC/)
        ! CONVERT UNIT TON/YEAR/GRID => UG/SEC/M2
        EMINTEX(I, J, :, 4) = EMINTEX(I, J, :, 4)*TONTOUG/(YEARTOSEC*  &
                              (INTEXLAT_DIS(I, J)*INTEXLON_DIS(I, J)))

        !===============================================================
        ! READ VOC EMISSION
        !===============================================================

        CALL INIT_INTEX_VOC_TMP
        !READ RESIDENTIAL BIOFUEL (DOB) EMISSION
        READ(11, *, END = 1000) TMPACET, TMPALK1, TMPALK2, TMPALK3,   &
                                TMPALK4, TMPALK5, TMPARO1, TMPARO2,   &
                                TMPBACL, TMPBALD, TMPCCHO, TMPCRES,   &
                                TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,   &
                                TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,   &
                                TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,   &
                                TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,   &
                                TMPRCHO, TMPTERP
        EMVOCINTEX(I, J, :, 1) = (/TMPACET, TMPALK1, TMPALK2, TMPALK3,&
                                   TMPALK4, TMPALK5, TMPARO1, TMPARO2,&
                                   TMPBACL, TMPBALD, TMPCCHO, TMPCRES,&
                                   TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,&
                                   TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,&
                                   TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,&
                                   TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,&
                                   TMPRCHO, TMPTERP/)

        CALL INIT_INTEX_VOC_TMP
        !READ RESIDENTIAL FOSSIL FUEL (DOF) EMISSION
        READ(12, *, END = 2000) TMPACET, TMPALK1, TMPALK2, TMPALK3,   &
                                TMPALK4, TMPALK5, TMPARO1, TMPARO2,   &
                                TMPBACL, TMPBALD, TMPCCHO, TMPCRES,   &
                                TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,   &
                                TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,   &
                                TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,   &
                                TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,   &
                                TMPRCHO, TMPTERP
        EMVOCINTEX(I, J, :, 2) = (/TMPACET, TMPALK1, TMPALK2, TMPALK3,& 
                                   TMPALK4, TMPALK5, TMPARO1, TMPARO2,&
                                   TMPBACL, TMPBALD, TMPCCHO, TMPCRES,&
                                   TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,&
                                   TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,&
                                   TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,&
                                   TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,&
                                   TMPRCHO, TMPTERP/)
        CALL INIT_INTEX_VOC_TMP
        !READ RESIDENTIAL NON-COMBUSTION (DOP) EMISSION
        READ(13, *, END = 3000) TMPACET, TMPALK1, TMPALK2, TMPALK3,   &
                                TMPALK4, TMPALK5, TMPARO1, TMPARO2,   &
                                TMPBACL, TMPBALD, TMPCCHO, TMPCRES,   &
                                TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,   &
                                TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,   &
                                TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,   &
                                TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,   &
                                TMPRCHO, TMPTERP
        EMVOCINTEX(I, J, :, 3) = (/TMPACET, TMPALK1, TMPALK2, TMPALK3,&
                                   TMPALK4, TMPALK5, TMPARO1, TMPARO2,&
                                   TMPBACL, TMPBALD, TMPCCHO, TMPCRES,&
                                   TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,&
                                   TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,&
                                   TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,&
                                   TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,&
                                   TMPRCHO, TMPTERP/)
        CALL INIT_INTEX_VOC_TMP
        !READ RESIDENTIAL INDUSTRY (IND) EMISSION
        READ(14, *, END = 4000) TMPACET, TMPALK1, TMPALK2, TMPALK3,   &
                                TMPALK4, TMPALK5, TMPARO1, TMPARO2,   &
                                TMPBACL, TMPBALD, TMPCCHO, TMPCRES,   &
                                TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,   &
                                TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,   &
                                TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,   &
                                TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,   &
                                TMPRCHO, TMPTERP
        EMVOCINTEX(I, J, :, 4) = (/TMPACET, TMPALK1, TMPALK2, TMPALK3,&
                                   TMPALK4, TMPALK5, TMPARO1, TMPARO2,&
                                   TMPBACL, TMPBALD, TMPCCHO, TMPCRES,&
                                   TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,&
                                   TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,&
                                   TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,&
                                   TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,&
                                   TMPRCHO, TMPTERP/)
        CALL INIT_INTEX_VOC_TMP
        !READ RESIDENTIAL POWER (POW) EMISSION
        READ(15, *, END = 5000) TMPACET, TMPALK1, TMPALK2, TMPALK3,   &
                                TMPALK4, TMPALK5, TMPARO1, TMPARO2,   &
                                TMPBACL, TMPBALD, TMPCCHO, TMPCRES,   &
                                TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,   &
                                TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,   &
                                TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,   &
                                TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,   &
                                TMPRCHO, TMPTERP
        EMVOCINTEX(I, J, :, 5) = (/TMPACET, TMPALK1, TMPALK2, TMPALK3,&
                                   TMPALK4, TMPALK5, TMPARO1, TMPARO2,&
                                   TMPBACL, TMPBALD, TMPCCHO, TMPCRES,&
                                   TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,&
                                   TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,&
                                   TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,&
                                   TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,&
                                   TMPRCHO, TMPTERP/)
        CALL INIT_INTEX_VOC_TMP
        !READ RESIDENTIAL TRANSPORTATION (TRA) EMISSION
        READ(16, *, END = 6000) TMPACET, TMPALK1, TMPALK2, TMPALK3,   &
                                TMPALK4, TMPALK5, TMPARO1, TMPARO2,   &
                                TMPBACL, TMPBALD, TMPCCHO, TMPCRES,   &
                                TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,   &
                                TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,   &
                                TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,   &
                                TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,   &
                                TMPRCHO, TMPTERP
        EMVOCINTEX(I, J, :, 6) = (/TMPACET, TMPALK1, TMPALK2, TMPALK3,&
                                   TMPALK4, TMPALK5, TMPARO1, TMPARO2,&
                                   TMPBACL, TMPBALD, TMPCCHO, TMPCRES,&
                                   TMPETHE, TMPGLY,  TMPHCHO, TMPIPRD,&
                                   TMPISOP, TMPMACR, TMPMEK,  TMPMEOH,&
                                   TMPMGLY, TMPMVK,  TMPNROG, TMPNVOL,&
                                   TMPOLE1, TMPOLE2, TMPPHEN, TMPPRD2,&
                                   TMPRCHO, TMPTERP/)

        ! CALCULATE AREA IN EACH GRID BOX
        INTEXLAT_DIS(I, J) = DEG_TO_DIS(                              &
                             INTEXLAT(I, J) - DIS/2, INTEXLON(I, J),  &
                             INTEXLAT(I, J) + DIS/2, INTEXLON(I, J))
        INTEXLON_DIS(I, J) = DEG_TO_DIS(                              &
                             INTEXLAT(I, J), INTEXLON(I, J) - DIS/2,  &
                             INTEXLAT(I, J), INTEXLON(I, J) + DIS/2)

        !===============================================================
        ! CALCULATE THE TOTAL EMISSION MASS IN EACH WRF-CHEM GRID BOX
        !===============================================================

        DO K = 1, E_SN(IDOM)-1
         DO L = 1, E_WE(IDOM)-1
          IF (WLAT(L, K) .GE. INTEXLAT(I, J) - DIS/2 .AND.    &
              WLAT(L, K) .LE. INTEXLAT(I, J) + DIS/2 .AND.    &
              WLON(L, K) .GE. INTEXLON(I, J) - DIS/2 .AND.    &
              WLON(L, K) .LE. INTEXLON(I, J) + DIS/2) THEN
           ! DO INTEX SPECIES LOOP TO SELECT RADM2 EMISSION SPECIES
           DO M = 1, NNAME-1
            INTEX3RD(L, 1, K, RADM2(M), :) = EMINTEX(I, J, M, 1) +    &
                                             EMINTEX(I, J, M, 2) +    &
                                             EMINTEX(I, J, M, 3) +    &
                                             EMINTEX(I, J, M, 4)
           ENDDO
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
100   CONTINUE
200   CONTINUE
300   CONTINUE
400   CONTINUE
1000  CONTINUE
2000  CONTINUE
3000  CONTINUE
4000  CONTINUE
5000  CONTINUE
6000  CONTINUE
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)


      


      ! END OF SUBROUTINE
      END SUBROUTINE INTEX

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!  $ID: INIT_INTEX_ARRAY V01 07/03/2012 09:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_INTEX_ARRAY INITIALIZES INTEX RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) INTEX3RD    (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NRADM, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/10/2012)
!******************************************************************************
!
      SUBROUTINE INIT_INTEX_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! INTEX OUTPUT VARIABELS
      !========================================================================

      ! INTEX3RD
      ALLOCATE(INTEX3RD(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NRADM, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('INTEX3RD')
      INTEX3RD = 0D0

      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_INTEX_ARRAY


!------------------------------------------------------------------------------
!
!  $ID: INIT_INTEX_TMP V01 08/19/2012 11:42 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_INTEX_TMP INITIALIZE ALL THE TEMPORARY VARIABLES.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) TMPROW      (INTEGER) : TEMPORARY ROW                            [---]
!  (2 ) TMPCOL      (INTEGER) : TEMPORARY COLUMN                         [---]
!  (3 ) TMPLAT      (REAL)    : TEMPORARY LATITUDE                       [deg]
!  (4 ) TMPLON      (REAL)    : TEMPORARY LONGITUDE                      [deg]
!  (5 ) TMPBC       (REAL)    : TEMPORARY BC                   [ton/year/grid]
!  (6 ) TMPCO       (REAL)    : TEMPORARY CO                   [ton/year/grid]
!  (7 ) TMPNOX      (REAL)    : TEMPORARY NOX                  [ton/year/grid]
!  (8 ) TMPOC       (REAL)    : TEMPORARY OC                   [ton/year/grid]
!  (9 ) TMPPM10     (REAL)    : TEMPORARY PM10                 [ton/year/grid]
!  (10) TMPPM25     (REAL)    : TEMPORARY PM25                 [ton/year/grid]
!  (11) TMPSO2      (REAL)    : TEMPORARY SO2                  [ton/year/grid]
!  (12) TMPVOC      (REAL)    : TEMPORARY VOC                  [ton/year/grid]
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (08/19/2012)
!******************************************************************************
!
      SUBROUTINE INIT_INTEX_TMP

      TMPROW  = 0
      TMPCOL  = 0
      TMPLAT  = 0.0
      TMPLON  = 0.0
      TMPBC   = 0.0
      TMPCO   = 0.0
      TMPNOX  = 0.0
      TMPOC   = 0.0
      TMPPM10 = 0.0
      TMPPM25 = 0.0
      TMPSO2  = 0.0
      TMPVOC  = 0.0

      END SUBROUTINE INIT_INTEX_TMP


!------------------------------------------------------------------------------
!
!  $ID: INIT_INTEX_VOC_TMP V01 08/19/2012 21:34 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_INTEX_VOC_TMP INITIALIZE VOC TEMPORARY VARIABLES
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (08/19/2012)
!******************************************************************************
!
      SUBROUTINE INIT_INTEX_VOC_TMP
      TMPROW  = 0
      TMPCOL  = 0
      TMPLAT  = 0.0
      TMPLON  = 0.0
      TMPACET = 0.0
      TMPALK1 = 0.0
      TMPALK2 = 0.0
      TMPALK3 = 0.0
      TMPALK4 = 0.0
      TMPALK5 = 0.0
      TMPARO1 = 0.0
      TMPARO2 = 0.0
      TMPBACL = 0.0
      TMPBALD = 0.0
      TMPCCHO = 0.0
      TMPCRES = 0.0
      TMPETHE = 0.0
      TMPGLY  = 0.0
      TMPHCHO = 0.0
      TMPIPRD = 0.0
      TMPISOP = 0.0
      TMPMACR = 0.0
      TMPMEK  = 0.0
      TMPMEOH = 0.0
      TMPMGLY = 0.0
      TMPMVK  = 0.0
      TMPNROG = 0.0
      TMPNVOL = 0.0
      TMPOLE1 = 0.0
      TMPOLE2 = 0.0
      TMPPHEN = 0.0
      TMPPRD2 = 0.0
      TMPRCHO = 0.0
      TMPTERP = 0.0

      END SUBROUTINE INIT_INTEX_VOC_TMP


!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_INTEX V01 06/18/2012 08:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_INTEX DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) INTEX3RD    (REAL)  : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NRADM, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT. (07/10/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_INTEX

      !=================================================================
      ! CLEANUP_INTEX BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( INTEX3RD   ) ) DEALLOCATE( INTEX3RD      )

      END SUBROUTINE CLEANUP_INTEX


      END MODULE INTEX_MOD
