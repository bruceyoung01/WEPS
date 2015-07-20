! $ID: SEVIRI_MOD.F V01 07/16/2012 22:37 BRUCE EXP$
!
!******************************************************************************
!  MODULE SEVIRI_MOD READS SEVIRI FIRE EMISSION.
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
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/16/2012)
!******************************************************************************
!
      MODULE SEVIRI_MOD

      ! USE OTHER MODULES
      USE MAP_UTILS
      USE NAMELIST_MOD, ONLY : F_PM25
      USE GFED_MOD,     ONLY : FE, ER, SE
      USE GFED_MOD,     ONLY : VERTICAL_LEVEL_FIRE
      USE GFED_MOD,     ONLY : READ_NC_2D, READ_NC_3D

      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO) :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES
      ! AND ROUTINES FROM BEING SEEN OUTSIDE 'SEVIRI_MOD.F'
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC :: SEVIRI3RD_PM25,SEVIRI3RD_OC,SEVIRI3RD_BC

      ! ... AND THESE ROUTINES
      PUBLIC :: SEVIRI
      PUBLIC :: INIT_SEVIRI_ARRAY
      PUBLIC :: CLEANUP_SEVIRI

      ! MODULE VARIABLES
      REAL, ALLOCATABLE      :: SEVIRI3RD_PM25(:,:,:,:)
      REAL, ALLOCATABLE      :: SEVIRI3RD_OC(:,:,:,:)
      REAL, ALLOCATABLE      :: SEVIRI3RD_BC(:,:,:,:)

      INTEGER, PARAMETER     :: NHR = 24

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: SEVIRI V01 07/16/2012 22:40 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE SEVIRI READS AND REWRITES SEVIRI FIRE EMISSION.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/16/2012)
!******************************************************************************
!
      SUBROUTINE SEVIRI (IDOM, IYEAR, IMONTH, IDAY)

      !REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD,   ONLY : INJ_HEIGHT
      USE NAMELIST_ARRAY_MOD,   ONLY : SEVIRIDIR, ERDIR
      USE NAMELIST_ARRAY_MOD,   ONLY : OUTPUTDIR
      USE PARAMETER_MOD,        ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,        ONLY : LATINC, LONINC
      USE NEI_MOD,              ONLY : ZFA


      ! ARGUMENTS
      INTEGER                 :: IDOM, IYEAR, IMONTH, IDAY

      ! PARAMETER
      INTEGER, PARAMETER      :: AHLINE  = 3
      INTEGER, PARAMETER      :: ALINE   = 180
      INTEGER, PARAMETER      :: HLINE   = 3
      INTEGER, PARAMETER      :: MAXL    = 5000
      INTEGER, PARAMETER      :: NCOLUMN = 75
      REAL,    PARAMETER      :: DIS     = 1.0
      INTEGER, PARAMETER      :: MT      = 24
      REAL*8,  PARAMETER      :: HRTOSED = 3600.
      REAL,    PARAMETER      :: SPV     = -9999.

      ! LOCAL VARIABLES
      CHARACTER (LEN = 255)   :: CTIME, DIR, PREFIX, INPTF, INPTF_ER_BC
      CHARACTER (LEN = 255)   :: INPTF_ER_OC, INPTF_ER_PM25
      CHARACTER (LEN = 255)   :: CYEAR, CMONTH, CDAY
      CHARACTER (LEN = 2  )   :: CDOM
      CHARACTER (LEN = 255)   :: SYSCMD
      INTEGER                 :: I,J,ILINE,NL,IPF,ICOLUMN,IHR,JNP,K
      INTEGER                 :: IBOT, ITOP
      REAL                    :: X, Y
      REAL*8                  :: ZDIF, FRC

      ! VARIABLES OF GRID CELL AREA FOR SEVIRI EMISSION
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT) :: FFLUX
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1)     :: ER_BC
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1)     :: ER_OC
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1)     :: ER_PM25
      REAL                    :: FRACT
      real, dimension(12) :: ndays ! the total day numbers
      data ndays /31,28,31,30,31,30,31,31,30,31,30,31/

      WRITE(CTIME, 110) IYEAR, IMONTH, IDAY
      WRITE(CYEAR,  '(I4)') IYEAR
      WRITE(CDOM , '(I2.2)') IDOM
      IF (IMONTH .LT. 10) THEN
         WRITE(CMONTH, '(I1)') IMONTH
      ELSE
         WRITE(CMONTH, '(I2)') IMONTH
      ENDIF
      IF (IDAY .LT. 10) THEN
         WRITE(CDAY, '(I1)') IDAY
      ELSE
         WRITE(CDAY, '(I2)') IDAY
      ENDIF
110   FORMAT(I4.4I2.2I2.2)

      IF (MOD(IYEAR,4) == 0) THEN
         NDAYS(2) = 29
      ELSE
         NDAYS(2) = 28
      ENDIF

      ER_BC = 0.0; ER_OC = 0.0; ER_PM25 = 0.0

      ! DEFINE 
      PREFIX    = 'SEVIRI_Emissions_NSSA_'

      !==========================================================================
      ! CALL NCL FROM FORTRAN CODE
      ! NOTE: NCL CAN REALIZE THE ASSIGNING STATEMENT BEFORE NCL SCRIPT, SO HERE 
      ! WE USE THIS FEATURE TO SPECIFY THE VARIABLES IN THE NCL SCRIPT TO GET THE 
      ! THE CORRECT YEAR, MONTH, AND DAY FOR DIFFERENT LOOPS. THE FORMAT LIKE:
      ! $ncl iy=2010 im=2 id=1 REGRID_SEVIRI.ncl
      !==========================================================================
      SYSCMD = 'ncl '//'iy='//TRIM(CYEAR)//' im='//TRIM(CMONTH)// &
               ' id='//TRIM(CDAY)//' REGRID_SEVIRI.ncl'
      WRITE(6, *) "SYSCMD : ", SYSCMD
      CALL SYSTEM(SYSCMD)

      ! CALL SUBROUTINE VERTICAL_LEVEL TO FIND K LEVEL OF
      ! SEVIRI EMISSION
      CALL VERTICAL_LEVEL_FIRE(INJ_HEIGHT, IBOT, ITOP, ZDIF)

      ! the original SEVIRI FLUX UNIT IS kg/km^2/hr
      FRACT = 10.**3/HRTOSED


      !================================================================
      ! READS SEVIRI EMISSION
      !================================================================

      DIR = SEVIRIDIR(1: LEN_TRIM(SEVIRIDIR))//'HOURLY_WRF'//'/'//    &
            CTIME(1: 4)
      INPTF = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'regrid/SEVIRI/' &
              //trim(CYEAR)//'/'//PREFIX(1: LEN_TRIM(PREFIX)) &
              //CTIME(1: LEN_TRIM(CTIME))//'_d'//CDOM//'.nc'
      WRITE(6, *) INPTF
      ! reading in the pm2.5 emission data [kg/km2/hr]
      CALL READ_NC_3D(inptf, SE, e_we(idom)-1,e_sn(idom)-1,mt,FFLUX)
 
      IF (.not. F_PM25) THEN
         ! if f_PM25 = .false., then reading the emission ratios to
         ! PM25
         INPTF_ER_BC = ERDIR(1: LEN_TRIM(ERDIR))//'BC_2_TPM/ER'//&
                       CTIME(1:6)//'_d'//CDOM//'.nc'
         INPTF_ER_OC = ERDIR(1: LEN_TRIM(ERDIR))//'OC_2_TPM/ER'//&
                       CTIME(1:6)//'_d'//CDOM//'.nc'

         CALL READ_NC_2D(inptf_er_bc, ER, e_we(idom)-1, &
                         e_sn(idom)-1,ER_BC)
         CALL READ_NC_2D(inptf_er_oc, ER, e_we(idom)-1, &
                         e_sn(idom)-1,ER_OC)
      ELSE
         INPTF_ER_PM25 = ERDIR(1: LEN_TRIM(ERDIR))//'PM2p5_2_TPM/ER'//&
                         CTIME(1:6)//'_d'//CDOM//'.nc'

         CALL READ_NC_2D(inptf_er_pm25, ER, e_we(idom)-1, &
                         e_sn(idom)-1, ER_PM25)
      ENDIF


      DO IHR = START_HOUR(IDOM),END_HOUR(IDOM)
            DO J = 1, E_SN(IDOM)-1
            DO I = 1,E_WE(IDOM)-1
               DO K = IBOT,ITOP
                  FRC = (ZFA(K+1)-ZFA(K))/ZDIF

                  IF (FFLUX(I,J,IHR) /= SPV) THEN
                     IF (F_PM25) THEN
                        SEVIRI3RD_PM25(I,K,J,IHR) = &
                        FFLUX(I,J,IHR)*ER_PM25(I,J)*FRC*FRACT
                     ELSE
                        SEVIRI3RD_OC(I,K,J,IHR) = &
                        FFLUX(I,J,IHR)*ER_OC(I,J)*FRC*FRACT
                        SEVIRI3RD_BC(I,K,J,IHR) = &
                        FFLUX(I,J,IHR)*ER_BC(I,J)*FRC*FRACT
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            ENDDO
      ENDDO

      END SUBROUTINE SEVIRI


!------------------------------------------------------------------------------
!  $ID: INIT_SEVIRI_ARRAY V01 07/04/2012 17:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_SEVIRI_ARRAY INITIALIZES SEVIRI RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) SEVIRI3RD_*  (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1,
!  NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE INIT_SEVIRI_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! SEVIRI OUTPUT VARIABELS
      !========================================================================

      ! SEVIRI3RD_PM25
      ALLOCATE(SEVIRI3RD_PM25(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('SEVIRI3RD_PM25')
      SEVIRI3RD_PM25 = 0D0

      ! SEVIRI3RD_OC
      ALLOCATE(SEVIRI3RD_OC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('SEVIRI3RD_OC')
      SEVIRI3RD_OC = 0D0

      ! SEVIRI3RD_BC
      ALLOCATE(SEVIRI3RD_BC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('SEVIRI3RD_BC')
      SEVIRI3RD_BC = 0D0

      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_SEVIRI_ARRAY

!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_SEVIRI V01 07/04/2012 17:50 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_SEVIRI DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) SEVIRI3RD_*    (REAL*8) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1,
!  NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT.
!  (07/04/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_SEVIRI

      !=================================================================
      ! CLEANUP_SEVIRI BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( SEVIRI3RD_PM25   ) ) DEALLOCATE( SEVIRI3RD_PM25    )
      IF ( ALLOCATED( SEVIRI3RD_OC     ) ) DEALLOCATE( SEVIRI3RD_OC      )
      IF ( ALLOCATED( SEVIRI3RD_BC     ) ) DEALLOCATE( SEVIRI3RD_BC      )

      END SUBROUTINE CLEANUP_SEVIRI

      END MODULE SEVIRI_MOD
