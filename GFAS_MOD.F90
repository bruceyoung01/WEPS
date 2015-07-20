!  $ID: GFAS_MOD.F V01 07/04/2012 15:42 BRUCE EXP$
!
!******************************************************************************
!  MODULE GFAS_MOD READS AND COMPUTES GFAS-B EMISSION AND REWRITES IT
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
      MODULE GFAS_MOD

      use netcdf

      ! USE OTHER MODULES
      USE NEI_MOD,      ONLY : NHR
      USE NAMELIST_MOD, ONLY : F_PM25
      USE GFED_MOD,     ONLY : FE, ER, SE
      USE GFED_MOD,     ONLY : VERTICAL_LEVEL_FIRE
      USE GFED_MOD,     ONLY : READ_NC_3D
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES 
      ! AND ROUTINES FROM BEING SEEN OUTSIDE "GFAS_MOD.F"
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC   :: GFAS3RD_OC, GFAS3RD_BC, GFAS3RD_PM25

      ! ... EXCEPT THESE ROUTINES
      PUBLIC   :: GFAS
      PUBLIC   :: INIT_GFAS_ARRAY
      PUBLIC   :: CLEANUP_GFAS

      !=================================================================
      ! PARAMETERS
      !=================================================================
      REAL, ALLOCATABLE                         :: GFAS3RD_OC(:,:,:,:)
      REAL, ALLOCATABLE                         :: GFAS3RD_BC(:,:,:,:)
      REAL, ALLOCATABLE                         :: GFAS3RD_PM25(:,:,:,:)


      ! ... EXCEPT THESE VARIABLES ...

      CONTAINS

!------------------------------------------------------------------------------
!  $ID: GFAS V01 08/23/2012 11:31 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE GFAS READS GFAS GLOBAL FIRE EMISSION DATABASE AND REWRITES
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
      SUBROUTINE GFAS(IDOM, IYEAR, IMONTH, IDAY)

      !REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : INJ_HEIGHT
      USE NAMELIST_ARRAY_MOD,   ONLY : GFASDIR
      USE NAMELIST_ARRAY_MOD,   ONLY : OUTPUTDIR
      USE NEI_MOD,              ONLY : ZFA

      ! ARGUMENTS
      INTEGER                 :: IDOM, IYEAR, IMONTH, IDAY

      ! PARAMETER
      INTEGER, PARAMETER      :: MT      = 8   
      REAL,    PARAMETER      :: HRTOSED = 10800.
      REAL,    PARAMETER      :: KGTOUG  = 10.0**6
      REAL,    PARAMETER      :: SPV     = -9999.

      CHARACTER (LEN = 255)   :: BBDIR, CTIME, DIR, CYEAR, CMONTH, CDAY
      CHARACTER (LEN = 255)   :: INPTF_BC, INPTF_OC, INPTF
      CHARACTER*2             :: CDOM
      CHARACTER (LEN = 255)   :: SYSCMD
      INTEGER                 :: I, J, K, L, M, ITOP, IBOT, N
      INTEGER                 :: IHR1, IHR2
      REAL*8                  :: FRC, ZDIF
      REAL                    :: FRACT, TEMP, TEMP1
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT) :: FLUX_BC
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT) :: FLUX_OC
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1, MT) :: FLUX_PM25

      INTEGER                 :: IDAY1,IMONTH1,IYEAR1
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

      ! SPECIFY DIRECTORY FOR GFAS FIRE EMISSION
      BBDIR = '/gfas_wrf/'

      ! READ GFAS FIRE EMISSION DATA WITH DIFFERENT YEARS, MONTHS
      DIR   = GFASDIR(1: LEN_TRIM(GFASDIR))//BBDIR(1: LEN_TRIM(BBDIR))

      !==========================================================================
      ! CALL NCL FROM FORTRAN CODE
      ! NOTE: NCL CAN REALIZE THE ASSIGNING STATEMENT BEFORE NCL SCRIPT, SO HERE 
      ! WE USE THIS FEATURE TO SPECIFY THE VARIABLES IN THE NCL SCRIPT TO GET THE 
      ! THE CORRECT YEAR, MONTH, AND DAY FOR DIFFERENT LOOPS. THE FORMAT LIKE:
      ! $ncl iy=2010 im=2 id=1 REGRID_GFED.ncl
      !==========================================================================
      SYSCMD = 'ncl '//'iy='//TRIM(CYEAR)//' im='//TRIM(CMONTH)// &
               ' id='//TRIM(CDAY)//' REGRID_GFAS.ncl'
      WRITE(6, *) "SYSCMD : ", SYSCMD
      CALL SYSTEM(SYSCMD)

      ! CALL SUBROUTINE CERTICAL_LEVEL TO FIND K LEVEL OF
      ! GFAS EMISSION
      CALL VERTICAL_LEVEL_FIRE(INJ_HEIGHT, IBOT, ITOP, ZDIF)

      FRACT = KGTOUG/HRTOSED

!============================
      IF (F_PM25) THEN    ! for PM2.5
!============================
         INPTF = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))// &
                 'regrid/GFAS/'//trim(CYEAR)//       &
                 '/GFASv1.0_biomassBurning_'//       &
                 CTIME(1:LEN_TRIM(CTIME))//          &
                 '_PM2_5_d'//CDOM//'.nc'

         ! reading in the pm2.5 emission data [g/m2/3hr]
         CALL READ_NC_3D(inptf, SE, e_we(idom)-1,e_sn(idom)-1,mt,FLUX_PM25)

         DO N = 1,MT
            IHR1 = (N-1)*3 + 1
            IHR2 = IHR1 + 2
            DO J = 1, E_SN(IDOM)-1
            DO I = 1,E_WE(IDOM)-1
               DO K = IBOT,ITOP
                  FRC = (ZFA(K+1)-ZFA(K))/ZDIF

                  ! UNIT OF FLUX_OC [G/M2/3HR]
                  ! UNIT OF GFAS3rd_oc [uG/M2/s]
                  IF (FLUX_PM25(I,J,N) /= SPV) THEN
                     TEMP = FLUX_PM25(I,J,N)*FRC*FRACT
                     GFAS3RD_PM25(I,K,J,IHR1:IHR2) = TEMP
                  ENDIF
               ENDDO
            ENDDO
            ENDDO
         ENDDO

!============================
      ELSE                ! for OC, BC
!============================
         INPTF_BC = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//            &
                    'regrid/GFAS/'//trim(CYEAR)//                  &
                    '/GFASv1.0_biomassBurning_'//                  &
                    CTIME(1:LEN_TRIM(CTIME))//'_BC_d'//CDOM//'.nc'
         INPTF_OC = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//            &
                    'regrid/GFAS/'//trim(CYEAR)//                  &
                    '/GFASv1.0_biomassBurning_'//                  &
                    CTIME(1:LEN_TRIM(CTIME))//'_OC_d'//CDOM//'.nc'
   
         ! reading in the BC emission data [g/m2/3hr]
         CALL READ_NC_3D(inptf_bc, SE, e_we(idom)-1, e_sn(idom)-1, &
                         mt, FLUX_BC)
         ! reading in the OC emission data [g/m2/3hr]
         CALL READ_NC_3D(inptf_oc, SE, e_we(idom)-1, e_sn(idom)-1, &
                         mt, FLUX_OC)


         DO N = 1,MT
            IHR1 = (N-1)*3 + 1
            IHR2 = IHR1 + 2
            DO J = 1, E_SN(IDOM)-1
            DO I = 1,E_WE(IDOM)-1
               DO K = IBOT,ITOP
                  FRC = (ZFA(K+1)-ZFA(K))/ZDIF

                  ! UNIT OF FLUX_OC [G/M2/3HR]
                  ! UNIT OF GFAS3rd_oc [uG/M2/s]
                  IF (FLUX_OC(I,J,N) /= SPV) THEN
                     TEMP = FLUX_OC(I,J,N)*FRC*FRACT
                     GFAS3RD_OC(I,K,J,IHR1:IHR2) = TEMP
                  ENDIF

                  ! UNIT OF FLUX_BC [G/M2/3HR]
                  ! UNIT OF GFAS3rd_bc [uG/M2/s]
                  IF (FLUX_BC(I,J,N) /= SPV) THEN
                     TEMP = FLUX_BC(I,J,N)*FRC*FRACT
                     GFAS3RD_BC(I,K,J,IHR1:IHR2) = TEMP
                  ENDIF
               ENDDO
            ENDDO
            ENDDO
         ENDDO
!============================
      ENDIF
!============================

      END SUBROUTINE GFAS


!------------------------------------------------------------------------------
!  $ID: INIT_GFAS_ARRAY V01 07/03/2012 09:39 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_GFAS_ARRAY INITIALIZES GFAS RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) GFAS3RD_BC    (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFAS3RD_OC    (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFAS3RD_PM25  (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/10/2012)
!******************************************************************************
!
      SUBROUTINE INIT_GFAS_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! GFAS OUTPUT VARIABELS
      !========================================================================

      ! GFAS3RD_BC
      ALLOCATE(GFAS3RD_BC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GFAS3RD_BC')
      GFAS3RD_BC = 0D0

      ! GFAS3RD_OC
      ALLOCATE(GFAS3RD_OC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GFAS3RD_OC')
      GFAS3RD_OC = 0D0

      ! GFAS3RD_PM25
      ALLOCATE(GFAS3RD_PM25(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR), &
               STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GFAS3RD_PM25')
      GFAS3RD_PM25 = 0D0

      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_GFAS_ARRAY


!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_GFAS V01 06/18/2012 08:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_GFAS DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) GFAS3RD_BC (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFAS3RD_OC (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!  (2 ) GFAS3RD_PM25 (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT. (07/10/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_GFAS

      !=================================================================
      ! CLEANUP_GFAS BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( GFAS3RD_OC)   ) DEALLOCATE( GFAS3RD_OC   )
      IF ( ALLOCATED( GFAS3RD_BC)   ) DEALLOCATE( GFAS3RD_BC   )
      IF ( ALLOCATED( GFAS3RD_PM25) ) DEALLOCATE( GFAS3RD_PM25 )

      END SUBROUTINE CLEANUP_GFAS

!===============================================================================

END MODULE GFAS_MOD
