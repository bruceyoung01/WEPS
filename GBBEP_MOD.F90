! $ID: GBBEP_MOD.F V01 06/15/2012 15:46 BRUCE EXP$
!
!******************************************************************************
!  MODULE GBBEP_MOD READS GBBEP FIRE EMISSION AND REWRITES IT INTO THE 
!  FORMAT WHICH WRF-CHEM CAN READS.
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
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/15/2012)
!  (2 ) ADD GLOBAL GBBEP DATA SWITCH (L_GLOBAL) BY BRUCE. (09/20/2013)
!  (3 ) ADD GBBEP EMISSION RATIO SWITCH (L_GBBEP_ER) BY BRUCE.  (09/20/2013)
!       GBBEP DATA ARE TPM, THE EMISSION RATIO ARE AS FOLLOWING:
!       PM2.5 : 0.00804
!       BC    : 0.000481
!       OC    : 0.00497
!       WHICH ARE PROVIDED BY GBBEP DEVELOPER XIAOYANG ZHANG
!       (xiaoyang.zhang@noaa.gov)
!******************************************************************************

      MODULE GBBEP_MOD

      ! USE OTHER MODULES
      USE MAP_UTILS
      USE NAMELIST_MOD, ONLY : F_PM25
      USE GFED_MOD,     ONLY : FE, ER, SE
      USE GFED_MOD,     ONLY : VERTICAL_LEVEL_FIRE
      USE GFED_MOD,     ONLY : READ_NC_2D
      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE
      TYPE(PROJ_INFO) :: PROJ

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- KEEP CERTAIN INTERNAL VARIABLES
      ! AND ROUTINES FROM BEING SEEN OUTSIDE 'GBBEP_MOD.F'
      !=================================================================

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE VARIABLES
      PUBLIC :: GBBEP3RD_PM25, GBBEP3RD_OC, GBBEP3RD_BC

      ! ... AND THESE ROUTINES
      PUBLIC :: GBBEP
      PUBLIC :: INIT_GBBEP_ARRAY
      PUBLIC :: CLEANUP_GBBEP

      ! MODULE VARIABLES
      REAL, ALLOCATABLE      :: GBBEP3RD_PM25(:,:,:,:)
      REAL, ALLOCATABLE      :: GBBEP3RD_OC(:,:,:,:)
      REAL, ALLOCATABLE      :: GBBEP3RD_BC(:,:,:,:)

      INTEGER, PARAMETER     :: NHR = 24

      CONTAINS

!------------------------------------------------------------------------------
!
!  $ID: GBBEP V01 07/04/2012 23:37 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE GBBEP READS AND REWRITES GBBEP FIRE EMISSION.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE GBBEP (IDOM, IYEAR, JD)

      !REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD,   ONLY : START_HOUR, END_HOUR
      USE NAMELIST_ARRAY_MOD,   ONLY : PROJ_CODE, DX, DY
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD,   ONLY : CORNER_LAT, CORNER_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : STAND_LON
      USE NAMELIST_ARRAY_MOD,   ONLY : TRUELAT1, TRUELAT2
      USE NAMELIST_ARRAY_MOD,   ONLY : INJ_HEIGHT
      USE NAMELIST_ARRAY_MOD,   ONLY : GBBEPDIR, ERDIR
      USE PARAMETER_MOD,        ONLY : KNOWNI, KNOWNJ
      USE PARAMETER_MOD,        ONLY : LATINC, LONINC
      USE NEI_MOD,              ONLY : ZFA

      ! ARGUMENTS
      INTEGER                 :: IDOM, IYEAR, JD

      ! PARAMETER
      INTEGER, PARAMETER      :: MAXL       = 50000
      REAL*8,  PARAMETER      :: HRTOSED    = 3600.
      ! SWITCH TO SELECT GLOBAL OR USA SMOKE EMISSION DATA
      ! IF L_GLOBAL = .TRUE.,  THEN SELECT GLOBAL GBBEP DATA
      ! IF L_GLOBAL = .FALSE., THEN SELECT USA GBBEP DATA
      LOGICAL, PARAMETER      :: L_GLOBAL   = .TRUE.
      ! SWITCH TO SELECT EMISSION RATIO FROM GFED OR FROM GBBEP
      ! IF L_GBBEP_ER = .TRUE.,  THEN SELECT GBBEP EMISSION RATIO
      ! IF L_GBBEP_ER = .FALSE., THEN SELECT GFED EMISSION RATIO
      LOGICAL, PARAMETER      :: L_GBBEP_ER = .TRUE.

      ! LOCAL VARIABLES
      CHARACTER (LEN = 255)   :: CTIME, DIR, INPTF, CNPTS, HEADER, & 
                                 INPTF_ER_PM25, INPTF_ER_BC, INPTF_ER_OC
      CHARACTER (LEN = 4  )   :: CYEAR
      CHARACTER (LEN = 2  )   :: CDOM, CMONTH
      INTEGER                 :: I, ILINE, NL, IPF, IHR, IMONTH
      INTEGER                 :: X1, Y1
      INTEGER                 :: IBOT, ITOP
      INTEGER                 :: VIOSTAT
      REAL                    :: X, Y
      REAL*8                  :: ZDIF, FRC
      REAL                    :: FRACT

      ! GBBEP FILE VARIABLES
      INTEGER                 :: NPTS, TMPYEAR, TMPJD
      REAL*8                  :: TMPLON, TMPLAT
      REAL*8                  :: TMPPM25_01, TMPPM25_02, TMPPM25_03
      REAL*8                  :: TMPPM25_04, TMPPM25_05, TMPPM25_06
      REAL*8                  :: TMPPM25_07, TMPPM25_08, TMPPM25_09
      REAL*8                  :: TMPPM25_10, TMPPM25_11, TMPPM25_12
      REAL*8                  :: TMPPM25_13, TMPPM25_14, TMPPM25_15
      REAL*8                  :: TMPPM25_16, TMPPM25_17, TMPPM25_18
      REAL*8                  :: TMPPM25_19, TMPPM25_20, TMPPM25_21
      REAL*8                  :: TMPPM25_22, TMPPM25_23, TMPPM25_24
      REAL, DIMENSION(MAXL)   :: FLAT, FLON
      REAL*8, DIMENSION(MAXL, NHR) :: PM25
      REAL, DIMENSION(E_WE(IDOM)-1, E_SN(IDOM)-1)::ER_PM25, ER_BC, ER_OC
      real, dimension(12) :: ndays, sum_days ! the total day numbers
      data ndays /31,28,31,30,31,30,31,31,30,31,30,31/

      if (mod(iyear,4) == 0) then
         ndays(2) = 29
      else
         ndays(2) = 28
      endif

      sum_days(1) = ndays(1)
      do i = 2,12
         sum_days(i) = sum_days(i-1)+ndays(i)
      enddo

      do i = 1,12
         if (JD > sum_days(i)) cycle
         imonth = i
         exit
      enddo
            
      ! CALL MAP_SET TO SPECIFY THE WRF USED MAP, INCLUDING GRID BOX
      CALL MAP_SET                                                     &
           (proj_code = PROJ_LC, proj = PROJ, lat1 = CORNER_LAT(IDOM), &
            lon1 = CORNER_LON(IDOM), knowni = KNOWNI, knownj = KNOWNJ, &
            dx = DX(IDOM), stdlon = STAND_LON, truelat1 = TRUELAT1,    &
            truelat2 = TRUELAT2)

      ! INITIALIZE GBBEP3RD
      GBBEP3RD_PM25 = 0.0
      GBBEP3RD_OC = 0.0
      GBBEP3RD_BC = 0.0

      ER_PM25 = 0.0; ER_BC = 0.0; ER_OC = 0.0

      WRITE(CDOM  , '(I2.2)') IDOM
      WRITE(CYEAR , '(I4.4)') IYEAR
      WRITE(CMONTH, '(I2.2)') IMONTH

      IF (L_GLOBAL) THEN
         IF (F_PM25) THEN
            IF (L_GBBEP_ER) THEN
               ER_PM25 = 0.00804
            ELSE
               INPTF_ER_PM25 = ERDIR(1:LEN_TRIM(ERDIR))//  &
                               'PM2p5_2_TPM/ER'            &
                               //CYEAR//CMONTH//'_d'//CDOM//'.nc'
               CALL READ_NC_2D(INPTF_ER_PM25, ER, e_we(idom)-1, &
                               e_sn(idom)-1, ER_PM25)
            ENDIF
         ! if F_PM25 = .false., then reading the emission ratios to
         ! BC AND OC
         ELSE
            IF (L_GBBEP_ER) THEN
               ER_BC = 0.000481
               ER_OC = 0.00497
            ELSE
               INPTF_ER_BC = ERDIR(1: LEN_TRIM(ERDIR))//'BC_2_TPM/ER' &
                             //CYEAR//CMONTH//'_d'//CDOM//'.nc'
               INPTF_ER_OC = ERDIR(1: LEN_TRIM(ERDIR))//'OC_2_TPM/ER' &
                             //CYEAR//CMONTH//'_d'//CDOM//'.nc'
               CALL READ_NC_2D(inptf_er_bc, ER, e_we(idom)-1, &
                               e_sn(idom)-1, ER_BC)
               CALL READ_NC_2D(inptf_er_oc, ER, e_we(idom)-1, &
                               e_sn(idom)-1, ER_OC)
            ENDIF
         ENDIF
      ELSE
         IF (F_PM25) THEN
            IF (L_GBBEP_ER) THEN
               ER_PM25 = 1.0
            ENDIF
         ! if F_PM25 = .false., then reading the emission ratios to
         ! BC AND OC
         ELSE
            IF (L_GBBEP_ER) THEN
               ER_BC = 0.000481/0.00804
               ER_OC = 0.00497/0.00804
            ELSE
               INPTF_ER_BC = ERDIR(1: LEN_TRIM(ERDIR))//'BC_2_PM2p5/ER'&
                             //CYEAR//CMONTH//'_d'//CDOM//'.nc'
               INPTF_ER_OC = ERDIR(1: LEN_TRIM(ERDIR))//'OC_2_PM2p5/ER'&
                             //CYEAR//CMONTH//'_d'//CDOM//'.nc'
               CALL READ_NC_2D(inptf_er_bc, ER, e_we(idom)-1, &
                               e_sn(idom)-1, ER_BC)
               CALL READ_NC_2D(inptf_er_oc, ER, e_we(idom)-1, &
                               e_sn(idom)-1, ER_OC)
            ENDIF
         ENDIF
      ENDIF

      ! CALL SUBROUTINE CERTICAL_LEVEL TO FIND K LEVEL OF 
      ! GBBEP EMISSION
      CALL VERTICAL_LEVEL_FIRE(INJ_HEIGHT, IBOT, ITOP, ZDIF)

      ! DO HOUR LOOP
      FRACT = 10.**9/(DX(IDOM)*DY(IDOM))/HRTOSED

       !================================================================
       ! READS GBBEP EMISSION
       !================================================================

       ! DEFINE THE CHARACTER OF FILE NAME
       IF(JD .LT. 100)THEN
        WRITE(CTIME,'(I4.4A1I2.2)') IYEAR,'0',JD
       ELSE
        WRITE(CTIME,'(I4.4I3.3)')IYEAR,JD
       ENDIF

       ! GLOBAL GBBEP DATA DIRECTORY
       IF (L_GLOBAL) THEN
          DIR  = GBBEPDIR(1: LEN_TRIM(GBBEPDIR))//'GLOBAL/'//CTIME(1:4)
          INPTF= DIR(1: LEN_TRIM(DIR))//'/f'//CTIME(1: LEN_TRIM(CTIME))&
                  //'_HourlyMass_GOES11_13MET09MTS01.txt'
          OPEN(1, FILE=INPTF, IOSTAT = VIOSTAT, STATUS='OLD')
          IF (VIOSTAT == 0) THEN
             READ(1, *) CNPTS
             READ(1, *) HEADER

          ELSE
             WRITE(*, *) '--------------------------------------'
             WRITE(*, *) INPTF, ' DOES NOT EXIST'
             WRITE(*, *) 'ASSIGNING THE TOTAL MASS TO 0.0'
             WRITE(*, *) '--------------------------------------'
             ILINE               = 2
             FLON(ILINE-1)       = 0.0
             FLAT(ILINE-1)       = 0.0
             PM25(ILINE-1, 1:24) = 0.0
          ENDIF
       ELSE
          DIR  = GBBEPDIR(1: LEN_TRIM(GBBEPDIR))//'USA/'//CTIME(1:4)
          INPTF= DIR(1: LEN_TRIM(DIR))//'/f'//CTIME(1: LEN_TRIM(CTIME))&
                  //'_PM25_emis_hourly'
          OPEN(1, FILE=INPTF, IOSTAT = VIOSTAT, STATUS='OLD')
          IF (VIOSTAT == 0) THEN
          READ(1, *) HEADER
          ELSE
             WRITE(*, *) '--------------------------------------'
             WRITE(*, *) INPTF, ' DOES NOT EXIST'
             WRITE(*, *) 'ASSIGNING THE TOTAL MASS TO 0.0'
             WRITE(*, *) '--------------------------------------'
             ILINE               = 2
             FLON(ILINE-1)       = 0.0
             FLAT(ILINE-1)       = 0.0
             PM25(ILINE-1, 1:24) = 0.0
          ENDIF
       ENDIF
       IF (VIOSTAT == 0) THEN
        DO ILINE = 1, MAXL
           TMPYEAR    = 0
           TMPJD      = 0
           TMPLON     = 0.0
           TMPLAT     = 0.0
           TMPPM25_01 = 0.0
           TMPPM25_02 = 0.0
           TMPPM25_03 = 0.0
           TMPPM25_04 = 0.0
           TMPPM25_05 = 0.0
           TMPPM25_06 = 0.0
           TMPPM25_07 = 0.0
           TMPPM25_08 = 0.0
           TMPPM25_09 = 0.0
           TMPPM25_10 = 0.0
           TMPPM25_11 = 0.0
           TMPPM25_12 = 0.0
           TMPPM25_13 = 0.0
           TMPPM25_14 = 0.0
           TMPPM25_15 = 0.0
           TMPPM25_16 = 0.0
           TMPPM25_17 = 0.0
           TMPPM25_18 = 0.0
           TMPPM25_19 = 0.0
           TMPPM25_20 = 0.0
           TMPPM25_21 = 0.0
           TMPPM25_22 = 0.0
           TMPPM25_23 = 0.0
           TMPPM25_24 = 0.0
           READ(1, *, END=100)TMPYEAR, TMPJD, TMPLON, TMPLAT,      &
                              TMPPM25_01, TMPPM25_02, TMPPM25_03,  &
                              TMPPM25_04, TMPPM25_05, TMPPM25_06,  &
                              TMPPM25_07, TMPPM25_08, TMPPM25_09,  &
                              TMPPM25_10, TMPPM25_11, TMPPM25_12,  &
                              TMPPM25_13, TMPPM25_14, TMPPM25_15,  &
                              TMPPM25_16, TMPPM25_17, TMPPM25_18,  &
                              TMPPM25_19, TMPPM25_20, TMPPM25_21,  &
                              TMPPM25_22, TMPPM25_23, TMPPM25_24
           FLON(ILINE) = TMPLON
           FLAT(ILINE) = TMPLAT
           PM25(ILINE, 1:24) = (/TMPPM25_01, TMPPM25_02, TMPPM25_03, &
                                 TMPPM25_04, TMPPM25_05, TMPPM25_06, & 
                                 TMPPM25_07, TMPPM25_08, TMPPM25_09, &
                                 TMPPM25_10, TMPPM25_11, TMPPM25_12, & 
                                 TMPPM25_13, TMPPM25_14, TMPPM25_15, &
                                 TMPPM25_16, TMPPM25_17, TMPPM25_18, &
                                 TMPPM25_19, TMPPM25_20, TMPPM25_21, &
                                 TMPPM25_22, TMPPM25_23, TMPPM25_24/)

         ENDDO
100      CONTINUE
         CLOSE(1)
       ENDIF

       ! REAL # OF LINES IN THE FILES
       NL = ILINE - 1
       ! CALCULATE THE SMOKE EMISSION RATE OVER PER WRF/CHEM GRID BOX
       DO IPF = 1, NL
        CALL LATLON_TO_IJ(PROJ,FLAT(IPF),FLON(IPF),X,Y)
! for each value has been put at mid-point, e.g, at 1.5, 2.5 etc,
! instead of 1.0,2.0,...
! zf
!        IF(X .GE. 1 .AND. X .LE. E_WE(IDOM)-1 .AND.    &
!           Y .GE. 1 .AND. Y .LE. E_SN(IDOM)-1 )THEN
!        IF(X-INT(X) .GT. 0.5) X1 = INT(X+0.5)
!        IF(X-INT(X) .LE. 0.5) X1 = INT(X)
!        IF(Y-INT(Y) .GT. 0.5) Y1 = INT(Y+0.5)
!        IF(Y-INT(Y) .LE. 0.5) Y1 = INT(Y)
        IF (X .GE. 1 .AND. X .LT. E_WE(IDOM) .AND.     &
            Y .GE. 1 .AND. Y .LT. E_SN(IDOM) )THEN
           X1 = INT(X)
           Y1 = INT(Y)

        DO I = IBOT,ITOP
         FRC = (ZFA(I+1)-ZFA(I))/ZDIF
         IF (L_GLOBAL) THEN
            IF (F_PM25) THEN
               GBBEP3RD_PM25(X1,I,Y1,:)=GBBEP3RD_PM25(X1,I,Y1,:) + &
                                        PM25(IPF,:)*ER_PM25(X1,Y1) &
                                        *FRC*FRACT
            ELSE
               GBBEP3RD_OC(X1,I,Y1,:)  =GBBEP3RD_OC(X1,I,Y1,:)   + &
                                        PM25(IPF,:)*ER_OC(X1,Y1)   &
                                        *FRC*FRACT
               GBBEP3RD_BC(X1,I,Y1,:)  =GBBEP3RD_BC(X1,I,Y1,:)   + &
                                        PM25(IPF,:)*ER_BC(X1,Y1)   &
                                        *FRC*FRACT
            ENDIF
         ENDIF
        ENDDO ! END OF I LOOP
        ENDIF
       ENDDO ! END OF IPF LOOP

      ! END OF SUBROUTINE
      END SUBROUTINE GBBEP

!------------------------------------------------------------------------------
!  $ID: INIT_GBBEP_ARRAY V01 07/04/2012 17:48 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE INIT_GBBEP_ARRAY INITIALIZES GBBEP RELATED ARRAY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) GBBEP3RD_*  (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/04/2012)
!******************************************************************************
!
      SUBROUTINE INIT_GBBEP_ARRAY(IDOM)

      ! REFERENCES TO F90 MODULES
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE NAMELIST_ARRAY_MOD,   ONLY : E_WE, E_SN, E_VERT

      ! LOCAL VARIABLES
      INTEGER      :: AS, IDOM

      !========================================================================
      ! GBBEP OUTPUT VARIABELS
      !========================================================================

      ! GBBEP3RD_PM25
      ALLOCATE(GBBEP3RD_PM25(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GBBEP3RD_PM25')
      GBBEP3RD_PM25 = 0D0

      ! GBBEP3RD_OC
      ALLOCATE(GBBEP3RD_OC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GBBEP3RD_OC')
      GBBEP3RD_OC = 0D0

      ! GBBEP3RD_BC
      ALLOCATE(GBBEP3RD_BC(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR),STAT = AS)
      IF (AS /= 0) CALL ALLOC_ERR('GBBEP3RD_BC')
      GBBEP3RD_BC = 0D0

      ! RETURN TO THE CALLING ROUTINE
      END SUBROUTINE INIT_GBBEP_ARRAY

!------------------------------------------------------------------------------
!
!  $ID: CLEANUP_GBBEP V01 07/04/2012 17:50 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE CLEANUP_GBBEP DEALLOCATES ALL MODULE ARRAYS
!
!  VARIABLES:
!  ============================================================================
!  (1 ) GBBEP3RD_*  (REAL) : (E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1, NHR)
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE REFERING TO GEOS-CHEM FORMAT.(07/04/2012)
!******************************************************************************
!
      SUBROUTINE CLEANUP_GBBEP

      !=================================================================
      ! CLEANUP_GBBEP BEGINS HERE!
      !=================================================================
      IF ( ALLOCATED( GBBEP3RD_PM25   ) ) DEALLOCATE( GBBEP3RD_PM25      )
      IF ( ALLOCATED( GBBEP3RD_OC     ) ) DEALLOCATE( GBBEP3RD_OC        )
      IF ( ALLOCATED( GBBEP3RD_BC     ) ) DEALLOCATE( GBBEP3RD_BC        )

      END SUBROUTINE CLEANUP_GBBEP

      END MODULE GBBEP_MOD
