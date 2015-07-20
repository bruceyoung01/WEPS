!  $ID: GFED_MOD.F V01 07/04/2012 15:42 BRUCE EXP$
!
!******************************************************************************
!  MODULE GFED_MOD READS AND COMPUTES GFED-B EMISSION AND REWRITES IT
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
!  (2 ) MODIFIED BY DR. FENG ZHANG (09/2013)
!  (3 ) ADD SUBROUTINE TO PREPROCESS GFED MONTHLY, DAILY, AND 3-HOURLY DATA.
!       BRUCE (11/24/2013)
!  (4 ) INTEGRATE NCL CODE WITH FORTRAN CODE, CONSIDERING FORTRAN CODE CAN 
!       CALL NCL CODE. BRUCE (11/18/2013)
!******************************************************************************
!
      MODULE GFED_MOD

      use netcdf

      ! USE OTHER MODULES
      USE NEI_MOD,      ONLY : NHR
      USE NAMELIST_MOD, ONLY : F_PM25
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES 
      ! AND ROUTINES FROM BEING SEEN OUTSIDE "GFED_MOD.F"
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC   :: FE
      PUBLIC   :: ER
      PUBLIC   :: SE
      PUBLIC   :: GFED3RD_OC
      PUBLIC   :: GFED3RD_BC
      PUBLIC   :: GFED3RD_PM25

      ! ... EXCEPT THESE ROUTINES
      PUBLIC   :: GFED
      PUBLIC   :: VERTICAL_LEVEL_FIRE
      PUBLIC   :: READ_NC_2D
      PUBLIC   :: READ_NC_3D
      PUBLIC   :: INIT_GFED_ARRAY
      PUBLIC   :: CLEANUP_GFED

      !=================================================================
      ! PARAMETERS
      !=================================================================
      LOGICAL, PUBLIC, PARAMETER :: F_BB = .FALSE. ! if .true., output biosphere fire

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      CHARACTER*255                             :: FE
      CHARACTER*255                             :: ER
      CHARACTER*255                             :: SE

      REAL, ALLOCATABLE                         :: GFED3RD_OC(:,:,:,:)
      REAL, ALLOCATABLE                         :: GFED3RD_BC(:,:,:,:)
      REAL, ALLOCATABLE                         :: GFED3RD_PM25(:,:,:,:)

      ! ... EXCEPT THESE VARIABLES ...

      !=================================================================
      ! MODULE ROUTINES --- FOLLOW BELOW THE "CONTAINS" STATEMENT
      !=================================================================

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: GFED_PREROCESS V01 11/24/2013 16:28 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE GFED_PREPROCESS READS GFEDV3 MONTHLY, DAILY, AND 3-HOURLY DATA 
!  TO CALCULATE DAILY 3-HOURLY DATA.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!ction_of_Emissions  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (11/24/2013)
!  (2 ) THE SCALING FACTOR FOR 3-HR, AND DAILY FROM .NC DATA FILES ARE FROM 
!       SOUTH TO NORTH, I.E., -89.75~89.75; HOWEVER, THE MONTHLY MEAN SMOKE 
!       EMISSION DATA IN DIRECTORY GFED3_MONTHLY/fields ARE FROM NORTH TO 
!       SOUTH, I.E., 89.75~-89.75.
!       OUTPUTS OF THIS PROCEDURE ARE FROM -89.75~89.75
!******************************************************************************
!
      SUBROUTINE GFED_PREPROCESS (ILOGIC, IYEAR, IMONTH, IDAY, NDAYS)

      !REFERENCES TO F90 MODULES
      USE NETCDF
      USE NAMELIST_ARRAY_MOD,   ONLY : GFEDDIR
      USE NAMELIST_ARRAY_MOD,   ONLY : N_SPECIES, SPECIES
      USE NAMELIST_ARRAY_MOD,   ONLY : NBASE_SPECIES, BASE_SPECIES

      !ARGUMENTS
      LOGICAL                  :: ILOGIC
      INTEGER                  :: IYEAR
      INTEGER                  :: IMONTH
      INTEGER                  :: IDAY
      REAL, DIMENSION(12)      :: NDAYS

      !DECARE PARAMETERS
      INTEGER, PARAMETER :: MX = 720, MY = 360, MT = 8

      !DECLARE VARIABLES
      INTEGER                  :: STRT(3),  COUT(3)
      INTEGER                  :: STRT1(2), COUT1(2)
      INTEGER                  :: NCID, AIRID, NX, NY, NT
      INTEGER                  :: NDIMS, NATTS, XTYPE, IS, LEN
      INTEGER                  :: DIMIDS(NF90_MAX_VAR_DIMS)
      INTEGER                  :: IY, IM, ID, IP, IPB, I, J, L
      INTEGER                  :: IT, N, NJ

      REAL, DIMENSION (MX, MY) ::  &
      SCALE_DAILY    ! DAILY SCALING FACTOR
      REAL, DIMENSION (MX, MY, N_SPECIES) ::  &
      FLUX           ! GFED MONTHLY EMISSION FOR SPECIES [g/m2/month]
      REAL, DIMENSION (MX, MY, NBASE_SPECIES) ::  &
      FLUXB          ! GFED MONTHLY EMISSION FOR BASE SPECIES [g/m2/month]
      REAL, DIMENSION (MX, MY, MT) ::  &
      SCALE_3HR      ! MONTHLY MEAN 3HR SCALING FACTOR
      REAL, DIMENSION (MX, MY, MT) :: &
      FLUX2          ! GFED 3HR EMISSION FOR SPECIES [g/m2/3hr]
      REAL, DIMENSION (MX, MY) :: &
      ER_B           ! MONTHLY EMISSION RATIO
      REAL                     :: SPV = -999.
      REAL, DIMENSION (MX, MY) :: TMP
      REAL                     :: FACT_DAILY, FACT_TOT
      REAL                     :: TMP1, TMP2
      CHARACTER*4              :: CYEAR
      CHARACTER*2              :: CMONTH, CDAY
      CHARACTER*255            :: FLNM_MONTH, FLNM_DAILY, FLNM_3HR
      CHARACTER*255            :: FLNM_3HR_OUT, FLNM_ER

      !=========================================================================
      ! SUBROUTINE STARTS HERE
      !=========================================================================
      WRITE(CYEAR,  '(I4.4)') IYEAR
      WRITE(CMONTH, '(I2.2)') IMONTH
      WRITE(CDAY,   '(I2.2)') IDAY

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! START TO SPECIES
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      !*************************************************************************
      ! NOTE : 
      ! START TO READ MONTHLY AVERAGE DATA AND 3-HOURLY FACTOR OF EMISSIONS DATA
      ! IF THIS IS THE FIRST TIME TO READ MONTHLY AND 3-HOURLY DATA IN THIS 
      ! MONTH, THEN START TO READ THE MONTHLY AND 3-HOURLY DATA.
      ! IF THIS IS NOT THE FIRST TIME TO READ MONTHLY AND 3-HOURLY DATA IN THIS 
      ! MONTH, THEN SKIP TO READ THE MONTHLY AND 3-HOURLY DATA.
      ! BECAUSE OF THE SAME MONTH, MONTHLY AND 3-HOURLY DATA THE SAME AS WELL.
      !*************************************************************************
      DO IP = 1, N_SPECIES
      IF (ILOGIC) THEN
         !READ GFEDV3 MONTHLY EMISSION DATA FOR SPECIES
         FLNM_MONTH = TRIM(GFEDDIR)//'GFED3_MONTHLY'//'/fields/'//&
                      TRIM(SPECIES(IP))//'/'//'GFED3.1_'//CYEAR//&
                      CMONTH//'_'//TRIM(SPECIES(IP))//'.txt'
         PRINT *, 'FLNM_MONTH = ', TRIM(FLNM_MONTH)
         OPEN (15, FILE = TRIM(FLNM_MONTH))
         DO J = 1, MY
            !FROM NORTH TO SOUTH
            READ(15, *) TMP(:, J)
            NJ = MY-J+1
            FLUX(:, NJ, IP) = TMP(:, J)
         ENDDO
         CLOSE(15)

         !================================================================
         ! READ GFEDV3 3HR SCALING FACTOR
         !================================================================
         FLNM_3HR = TRIM(GFEDDIR)//'GFED3_3HOURLY/'//CYEAR// &
                    '/fraction_emissions_'//CYEAR//CMONTH//'.nc'
         PRINT *, 'FLNM_3HR = ', TRIM(FLNM_3HR)
         ! CALL SUBROUTINE READ_NC_3D TO READ GFEDV3 3HR FRACTION OF 
         ! EMISSIONS
         CALL READ_NC_3D(FLNM_3HR, FE, MX, MY, MT, SCALE_3HR)
         WHERE (SCALE_3HR < -100.)
            SCALE_3HR = 1./MT
         ENDWHERE
      ENDIF ! ILOGIC

      FACT_DAILY = 1./NDAYS(IMONTH)
      FACT_TOT   = FACT_DAILY*(1./MT)

      !================================================================
      ! READ THE DAILY SCALING FACTOR
      !================================================================
      FLNM_DAILY = TRIM(GFEDDIR)//'GFED3_DAILY/'//CYEAR// &
                   '/fraction_emissions_'//CYEAR//CMONTH//CDAY//'.nc'
      PRINT *, 'FLNM_DAILY = ', TRIM(FLNM_DAILY)
      ! CALL SUBROUTINE READ_NC_2D TO READ GFEDV3 DAILY FRACTION OF 
      ! EMISSIONS
      CALL READ_NC_2D(FLNM_DAILY, FE, MX, MY, SCALE_DAILY)
      WHERE (SCALE_DAILY < -100.)
         SCALE_DAILY = FACT_DAILY
      ENDWHERE

      !================================================================
      ! OPEN FILE FOR WRITING CALCULATED 3-HOURLY EMISSION FLUX
      !================================================================
      FLNM_3HR_OUT = TRIM(GFEDDIR)//'emission_flux_3hourly/'//CYEAR// &
                     '/'//TRIM(SPECIES(IP))//'/GFED3.1_'//CYEAR//     &
                     CMONTH//CDAY//'_'//TRIM(SPECIES(IP))//'.txt'
      OPEN (20, FILE = TRIM(FLNM_3HR_OUT))

      DO N = 1, MT
         DO J = 1, MY
            DO I = 1, MX
               IF (SUM(SCALE_3HR(I, J, 1:MT)) == 0.0 .AND. &
                       SCALE_DAILY(I, J)      == 0.0) THEN
                  FLUX2(I, J, N) = &
                  FLUX (I, J, IP)*FACT_TOT
               ELSEIF (SUM(SCALE_3HR(I, J, 1:MT)) == 0.0) THEN
                  FLUX2(I, J, N) = &
                  FLUX(I, J, IP)*SCALE_DAILY(I, J)*(1./MT)
               ELSEIF (SCALE_DAILY(I, J) == 0.0) THEN
                  FLUX2(I, J, N) = &
                  FLUX(I, J, IP)*FACT_DAILY*SCALE_3HR(I, J, N)
               ELSE
                  FLUX2(I, J, N) = &
                  FLUX(I, J, IP)*SCALE_DAILY(I, J)*SCALE_3HR(I, J, N)
               ENDIF
            ENDDO
            WRITE(20, *) FLUX2(:, J, N)
         ENDDO
      ENDDO

      DO J = 1, MY
         DO I = 1, MX
            IF (FLUX(I, J, IP) /= 0.0 .AND. &
                SUM(FLUX2(I, J, 1:MT)) == 0.0) THEN
               PRINT *, 'ERROR HAPPENED!!!'
               PRINT *, 'FLUX  = ', FLUX(I, J, IP), &
                        'FLUX2 = ', FLUX2(I, J, 1:MT)
               STOP
            ENDIF
         ENDDO
      ENDDO
      ENDDO

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! START TO READ BASE SPECIES
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      DO IPB = 1, NBASE_SPECIES
      IF (ILOGIC) THEN
         !READ GFEDV3 MONTHLY EMISSION DATA FOR SPECIES
         FLNM_MONTH = TRIM(GFEDDIR)//'GFED3_MONTHLY'//'/fields/'//&
                      TRIM(BASE_SPECIES(IPB))//'/'//'GFED3.1_'//CYEAR//&
                      CMONTH//'_'//TRIM(BASE_SPECIES(IPB))//'.txt'
         PRINT *, 'FLNM_MONTH = ', TRIM(FLNM_MONTH)
         OPEN (15, FILE = TRIM(FLNM_MONTH))
         DO J = 1, MY
            !FROM NORTH TO SOUTH
            READ(15, *) TMP(:, J)
            NJ = MY-J+1
            FLUXB(:, NJ, IPB) = TMP(:, J)
         ENDDO
         CLOSE(15)
      ENDIF
      ENDDO

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! START TO CALCULATE EMISSION RATIO (SPECIES/BASE_SPECIES)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF (ILOGIC) THEN
      DO IPB = 1, NBASE_SPECIES
         DO IP = 1, N_SPECIES
            ER_B = 0.0
            ! OPEN A NEW FILE TO WRITE EMISSION RATIO
            FLNM_ER = TRIM(GFEDDIR)//'/emission_ratio/'// &
                      TRIM(SPECIES(IP))//'_2_'//          &
                      TRIM(BASE_SPECIES(IPB))//'/ER'//     &
                      CYEAR//CMONTH//'.txt'
            OPEN (23, FILE = TRIM(FLNM_ER))
            DO J = 1, MY
               DO I = 1, MX
                  IF (FLUXB(I, J, IPB) .GT. 0.0) THEN
                     ER_B(I, J) = FLUX(I, J, IP)/ &
                                  FLUXB(I, J, IPB)
                  ELSE
                     ER_B(I, J) = 0.0
                  ENDIF
               ENDDO ! MX
               ! WRITE EMISSION RATIO
               WRITE(23, *) ER_B(:, J)
            ENDDO !MY
            CLOSE(23)
         ENDDO ! IP
      ENDDO ! IPM
      ENDIF ! ILOGIC


      END SUBROUTINE GFED_PREPROCESS


!
!------------------------------------------------------------------------------
!
!  $ID: GFED V01 08/23/2012 11:31 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE GFED READS GFED GLOBAL FIRE EMISSION DATABASE AND REWRITES
!  IT INTO WRF-CHEM READABLE FORMAT.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (08/23/2012)
!******************************************************************************
!
      SUBROUTINE GFED(ILOGIC, IDOM, IYEAR, IMONTH, IDAY)

      !REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : INJ_HEIGHT
      USE NAMELIST_ARRAY_MOD,   ONLY : GFEDDIR
      USE NAMELIST_ARRAY_MOD,   ONLY : OUTPUTDIR
      USE NEI_MOD,              ONLY : ZFA

      ! ARGUMENTS
      LOGICAL                 :: ILOGIC
      INTEGER                 :: IDOM, IYEAR, IMONTH, IDAY

      ! PARAMETER
      INTEGER, PARAMETER      :: MT      = 8   
      REAL,    PARAMETER      :: HRTOSED = 10800.
      REAL,    PARAMETER      :: KGTOUG  = 10.0**6
      REAL,    PARAMETER      :: SPV     = -9999.

      CHARACTER (LEN = 255)   :: BBDIR, CTIME, DIR, CYEAR, CMONTH, CDAY
      CHARACTER (LEN = 255)   :: INPTF_BC, INPTF_OC, INPTF
      CHARACTER (LEN = 2  )   :: CDOM
      CHARACTER (LEN = 255)   :: SYSCMD           
      INTEGER                 :: I, J, K, L, M, ITOP, IBOT, N
      INTEGER                 :: IHR1, IHR2
      REAL*8                  :: FRC, ZDIF
      REAL                    :: FRACT, TEMP
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT)  :: FLUX_BC
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT)  :: FLUX_OC
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT)  :: FLUX_PM25

      REAL, DIMENSION(12)     :: NDAYS ! the total day numbers
      DATA NDAYS /31,28,31,30,31,30,31,31,30,31,30,31/


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

      !==========================================================================
      ! CALL SUBROUTINE GFED_PREPROCESS TO CALCULATE 3-HOURLY EMISSION DATA
      ! IN EVERY DAY BASED ON MONTHLY AVERAGED FLUX, DAILY FRACTION OF EMISSIONS, 
      ! AND 3-HOURLY FRACTION OF EMISSIONS.
      !==========================================================================
      CALL GFED_PREPROCESS(ILOGIC, IYEAR, IMONTH, IDAY, NDAYS)

      !==========================================================================
      ! CALL NCL FROM FORTRAN CODE
      ! NOTE: NCL CAN REALIZE THE ASSIGNING STATEMENT BEFORE NCL SCRIPT, SO HERE 
      ! WE USE THIS FEATURE TO SPECIFY THE VARIABLES IN THE NCL SCRIPT TO GET THE 
      ! THE CORRECT YEAR, MONTH, AND DAY FOR DIFFERENT LOOPS. THE FORMAT LIKE:
      ! $ncl iy=2010 im=2 id=1 REGRID_GFED.ncl
      !==========================================================================
      SYSCMD = 'ncl '//'iy='//TRIM(CYEAR)//' im='//TRIM(CMONTH)// &
               ' id='//TRIM(CDAY)//' REGRID_GFED.ncl'
      WRITE(6, *) "SYSCMD : ", SYSCMD
      CALL SYSTEM(SYSCMD)

      ! SPECIFY DIRECTORY FOR GFED FIRE EMISSION
      IF (F_BB) THEN
         BBDIR = '/GFED3_MONTHLY/GFED3.1_BB_wrf/'
      ELSE
         BBDIR = '/emission_flux_3hourly_wrf/'
      ENDIF

      ! READ GFED FIRE EMISSION DATA WITH DIFFERENT YEARS, MONTHS
      DIR   = GFEDDIR(1: LEN_TRIM(GFEDDIR))//BBDIR(1: LEN_TRIM(BBDIR))

      ! CALL SUBROUTINE CERTICAL_LEVEL TO FIND K LEVEL OF
      ! GFED EMISSION
      CALL VERTICAL_LEVEL_FIRE(INJ_HEIGHT, IBOT, ITOP, ZDIF)

      FRACT = KGTOUG/HRTOSED

!============================
      IF (F_PM25) THEN    ! for PM2.5
!============================
         IF (F_BB) THEN
            INPTF = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'regrid/GFED/' &
                    //'BB'//'/GFED3.1_'//CTIME(1:LEN_TRIM(CTIME))//   &
                    '_BB_d'//CDOM//'.nc'
         ELSE
            INPTF = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'regrid/GFED/' &
                    //trim(CYEAR)//'/PM2p5'// &
                    '/GFED3.1_'//CTIME(1:LEN_TRIM(CTIME))// &
                    '_PM2p5_d'//CDOM//'.nc'
         ENDIF

         ! reading in the pm2.5 emission data [g/m2/3hr]
         CALL READ_NC_3D(inptf, SE, e_we(idom)-1, e_sn(idom)-1, &
                         mt, FLUX_PM25)

         !=============================================================
         ! REDISTRIBUTE 3-HOURLY EMISSION DATA TO HOURLY.
         ! GFED 3-HOURLY DATA IS AVERAGE LIKE THIS :
         ! 1:00-3:00,   4:00-6:00,   7:00-9:00,   10:00-12:00, 
         ! 13:00-15:00, 16:00-18:00, 19:00-21:00, 22:00-24:00
         !=============================================================
         DO N = 1,MT
            IHR1 = (N-1)*3 + 1
            IHR2 = IHR1 + 2
            DO J = 1, E_SN(IDOM)-1
            DO I = 1, E_WE(IDOM)-1
               DO K = IBOT,ITOP
                  FRC = (ZFA(K+1)-ZFA(K))/ZDIF

                  ! UNIT OF FLUX_OC [G/M2/3HR]
                  ! UNIT OF gfed3rd_oc [uG/M2/s]
                  IF (FLUX_PM25(I,J,N) /= SPV) THEN
                     TEMP = FLUX_PM25(I,J,N)*FRC*FRACT
                        GFED3RD_PM25(I,K,J,IHR1:IHR2) = TEMP
                  ENDIF
               ENDDO
            ENDDO
            ENDDO
         ENDDO

      !============================
      ELSE                ! for OC, BC
      !============================
         INPTF_BC = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))// &
                    'regrid/GFED/'//trim(CYEAR)// &
                    '/GFED3.1_'//&
                    CTIME(1:LEN_TRIM(CTIME))//'_BC_d'//CDOM//'.nc'
         INPTF_OC = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))// &
                    'regrid/GFED/'//trim(CYEAR)// &
                    '/GFED3.1_'//&
                    CTIME(1:LEN_TRIM(CTIME))//'_OC_d'//CDOM//'.nc'
   
         ! reading in the BC emission data [g/m2/3hr]
         CALL READ_NC_3D(inptf_bc, SE, e_we(idom)-1, e_sn(idom)-1, &
                         mt, FLUX_BC)
         ! reading in the OC emission data [g/m2/3hr]
         CALL READ_NC_3D(inptf_oc, SE, e_we(idom)-1, e_sn(idom)-1, &
                         mt, FLUX_OC)

         !=============================================================
         ! REDISTRIBUTE 3-HOURLY EMISSION DATA TO HOURLY.
         ! GFED 3-HOURLY DATA IS AVERAGE LIKE THIS :
         ! 1:00-3:00,   4:00-6:00,   7:00-9:00,   10:00-12:00, 
         ! 13:00-15:00, 16:00-18:00, 19:00-21:00, 22:00-24:00
         !=============================================================

         DO N = 1, MT
            IHR1 = (N-1)*3 + 1
            IHR2 = IHR1 + 2
            DO J = 1, E_SN(IDOM)-1
            DO I = 1, E_WE(IDOM)-1
               DO K = IBOT,ITOP
                  FRC = (ZFA(K+1)-ZFA(K))/ZDIF

                  ! UNIT OF FLUX_OC [G/M2/3HR]
                  ! UNIT OF gfed3rd_oc [uG/M2/s]
                  IF (FLUX_OC(I,J,N) /= SPV) THEN
                     TEMP = FLUX_OC(I,J,N)*FRC*FRACT
                     GFED3RD_OC(I,K,J,IHR1:IHR2) = TEMP
                  ENDIF

                  ! UNIT OF FLUX_BC [G/M2/3HR]
                  ! UNIT OF gfed3rd_bc [uG/M2/s]
                  IF (FLUX_BC(I,J,N) /= SPV) THEN
                     TEMP = FLUX_BC(I,J,N)*FRC*FRACT
                     GFED3RD_BC(I,K,J,IHR1:IHR2) = TEMP
                  ENDIF
               ENDDO
            ENDDO
            ENDDO
         ENDDO
!============================
      ENDIF
!============================

      END SUBROUTINE GFED


!------------------------------------------------------------------------------
!
!  $ID: EMISSION_FACTOR V01 01/13/2014 16:34 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE EMISSION_FACTOR READS EMISSION FACTORS USED FOR DIFFERENT FIRE 
!  TYPES (g/kg).
!
!
!
!
!
!
!------------------------------------------------------------------------------
!
!  $ID: VERTICAL_LEVEL_FIRE V01 07/04/2012 17:28 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE VERTICAL_LEVEL_FIRE DETERMINES HEIGHT LEVEL OF EMISSION,
!  DEPENDING ON EMISSION HEIGHT.
!
!  VARIBALES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE VERTICAL_LEVEL_FIRE (HEIGHT, IBOT, ITOP, ZDIF)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : E_VERT
      USE NEI_MOD,            ONLY : ZFA

      ! ARGUMENTS
      REAL*8          :: HEIGHT, ZDIF
      INTEGER         :: IBOT, ITOP

      ! LOCAL VARIABLES
      REAL*8          :: BOT, TOP
      INTEGER         :: I

      TOP = HEIGHT
      BOT = 0.0
      DO I = E_VERT, 1, -1
       IF(ZFA(I+1).GT.BOT)THEN
        IBOT = I
       ENDIF
      ENDDO
      DO I = E_VERT, 1, -1
       IF(ZFA(I+1).GT.TOP)THEN
        ITOP = I
       ENDIF
      ENDDO
      IF(IBOT.GE.ITOP)THEN
       ITOP=IBOT+1
      ENDIF
      ZDIF = ZFA(ITOP)-ZFA(IBOT)

      END SUBROUTINE VERTICAL_LEVEL_FIRE

!------------------------------------------------------------------------------
!
!  $ID: READ_NC_2D V01 11/25/2013 15:10 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_NC_2D READS 2-D VARIABLES FROM NC FILE.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (11/25/2013)
!******************************************************************************

      SUBROUTINE READ_NC_2D(FLNM, VAR, MX, MY, FLUX)

      USE NETCDF
      INTEGER, INTENT(IN)       :: MX, MY
      CHARACTER*255, INTENT(IN) :: FLNM, VAR
      REAL, INTENT(OUT)         :: FLUX(MX, MY)

      CHARACTER*80 :: NAME
      INTEGER :: STRT(2), COUT(2)
      INTEGER :: NCID, AIRID, NX, NY, N
      INTEGER :: NDIMS, NATTS, XTYPE, IS, LEN
      INTEGER :: DIMIDS(NF90_MAX_VAR_DIMS)

      ! OPEN NC FILE
      IS = NF90_OPEN(TRIM(FLNM), NF90_NOWRITE, NCID)
      IF (IS /= NF90_NOERR) STOP

      IS = NF90_INQ_VARID(NCID, VAR, AIRID)
      IF (IS /= NF90_NOERR) STOP

      IS = NF90_INQUIRE_VARIABLE(NCID, AIRID, NAME, XTYPE, &
                                 NDIMS, DIMIDS, NATTS)
      IF (IS /= NF90_NOERR) STOP

      IS = NF90_INQUIRE_DIMENSION(NCID, DIMIDS(1), NAME, LEN = NX)
      IF (IS /= NF90_NOERR) STOP

      IS = NF90_INQUIRE_DIMENSION(NCID, DIMIDS(2), NAME, LEN = NY)
      IF (IS /= NF90_NOERR) STOP

      IF (MX /= NX .OR. MY /= NY) THEN
         PRINT *, 'ARRAY SIZES NOT RIGHT FOR DAILY SCALING DATA, STOP!'
         PRINT *, 'E_WE(IDOM)-1 = ', MX
         PRINT *, 'NX = ', NX
         PRINT *, 'E_SN(IDOM)-1 = ', MY
         PRINT *, 'NY = ', NY
         STOP
      ENDIF

      ! GET VARS
      STRT(1) = 1
      STRT(2) = 1
      COUT(1) = MX
      COUT(2) = MY
      IS = NF90_GET_VAR(NCID, AIRID, FLUX, STRT, COUT)
      IS = NF90_CLOSE(NCID)

      END SUBROUTINE READ_NC_2D


!------------------------------------------------------------------------------
!
!  $ID: READ_NC_3D V01 11/25/2013 15:10 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE READ_NC_3D READS 3-D VARIABLES FROM NC FILE.
!
!  VARIABLES:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (11/25/2013)
!******************************************************************************

     SUBROUTINE READ_NC_3D(FLNM, VAR, MX, MY, MT, FLUX)

      USE NETCDF
      INTEGER, INTENT(IN)       :: MX, MY, MT
      CHARACTER*255, INTENT(IN) :: FLNM, VAR
      REAL, INTENT(OUT)         :: FLUX(MX, MY, MT)

      CHARACTER*80 :: NAME
      INTEGER :: STRT(3), COUT(3)
      INTEGER :: NCID, AIRID, NX, NY, NT, N
      INTEGER :: NDIMS, NATTS, XTYPE, IS, LEN
      INTEGER :: DIMIDS(NF90_MAX_VAR_DIMS)

      IS = NF90_OPEN(TRIM(FLNM), NF90_NOWRITE, NCID)
      IF( IS /= NF90_NOERR ) STOP

      IS = NF90_INQ_VARID(NCID, VAR, AIRID)
      IF( IS /= NF90_NOERR ) STOP

      IS = NF90_INQUIRE_VARIABLE(NCID, AIRID, NAME, XTYPE, NDIMS, &
                                 DIMIDS, NATTS)
      IF( IS /= NF90_NOERR ) STOP

      IS = NF90_INQUIRE_DIMENSION(NCID, DIMIDS(1), NAME, LEN=NX)
      IF( IS /= NF90_NOERR ) STOP

      IS = NF90_INQUIRE_DIMENSION(NCID, DIMIDS(2), NAME, LEN=NY)
      IF( IS /= NF90_NOERR ) STOP

      IS = NF90_INQUIRE_DIMENSION(NCID, DIMIDS(3), NAME, LEN=NT)
      IF( IS /= NF90_NOERR ) STOP

      IF (NX /= MX .OR. NY /= MY) THEN
         PRINT *, 'xxxx, for BC, nx /= E_WE(IDOM)-1 .or.&
                   ny /= E_SN(IDOM)-1, stop'
         PRINT *, 'nx=', NX, 'E_WE(IDOM)-1=', MX,  &
                  'ny=', NY, 'E_SN(IDOM)-1=', MY
         STOP
      ENDIF

      DO N = 1, NT
         ! GET VARS
         STRT(1) = 1
         STRT(2) = 1
         STRT(3) = N

         COUT(1) = NX
         COUT(2) = NY
         COUT(3) = 1

         IS = NF90_GET_VAR(NCID, AIRID, FLUX(:,:,N), STRT, COUT)
      ENDDO
      IS = NF90_CLOSE(NCID)
 
      END SUBROUTINE READ_NC_3D


!------------------------------------------------------------------------------
!  $ID: INIT_GFED_ARRAY V01 07/03/2012 09:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_GFED_ARRAY INITIALIZES GFED RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) GFED3RD_BC    (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFED3RD_OC    (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFED3RD_PM25  (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/10/2012)
!******************************************************************************
!
      SUBROUTINE INIT_GFED_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! GFED OUTPUT VARIABELS
      !========================================================================

      ! GFED EMISSION TYPES
      FE = 'Fraction_of_Emissions'
      ER = 'emission_ratio'
      SE = 'smoke_emission'

      ! GFED3RD_BC
      ALLOCATE(GFED3RD_BC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GFED3RD_BC')
      GFED3RD_BC = 0D0

      ! GFED3RD_OC
      ALLOCATE(GFED3RD_OC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GFED3RD_OC')
      GFED3RD_OC = 0D0

      ! GFED3RD_PM25
      ALLOCATE(GFED3RD_PM25(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GFED3RD_PM25')
      GFED3RD_PM25 = 0D0


      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_GFED_ARRAY


!------------------------------------------------------------------------------
!
!  $ID: INTEGRATE_GFED V01 11/18/2013 15:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INTEGRATE_GFED INTEGRATES GFEDV3 MONTHLY, DAILY, AND 3-HOURLY 
!  DATA INTO 3-HOURLY DATA FOR EVERY DAY.
!  INSTRUCTION FOR GFEDV3 MONTHLY, DAILY, AND 3-HOURLY DATA:
!  (1 ) MONTHLY DATA : 
!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_GFED V01 06/18/2012 08:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_GFED DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) GFED3RD_BC (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFED3RD_OC (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFED3RD_PM25 (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT. (07/10/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_GFED

      !=================================================================
      ! CLEANUP_GFED BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( GFED3RD_OC) ) DEALLOCATE( GFED3RD_OC   )
      IF ( ALLOCATED( GFED3RD_BC) ) DEALLOCATE( GFED3RD_BC   )
      IF ( ALLOCATED( GFED3RD_PM25) ) DEALLOCATE( GFED3RD_PM25   )

      END SUBROUTINE CLEANUP_GFED

      END MODULE GFED_MOD
