! $ID: WRITE_MOD.F V01 07/03/2012 20:58 BRUCE EXP$
!
!******************************************************************************
!  MODULE WRITE_MOD CONTAINS SUBROUTINES WHICH IS USED TO WRITE EMISSIONS TO 
!  THE FORMAT WHICH CAN BE READ BY convert_emiss.exe
!
!  VARIABLES:
!  ============================================================================
!  (1 ) 
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (07/03/2012)
!******************************************************************************
!
      MODULE WRITE_MOD

      ! REFERENCES TO F90 MODULES
      USE NEI_MOD,            ONLY : NEI3RD
      USE INTEX_MOD,          ONLY : INTEX3RD
      USE FLAMBE_MOD,         ONLY : FLAMBE3RD_BC, FLAMBE3RD_OC, FLAMBE3RD_PM25
      USE GBBEP_MOD,          ONLY : GBBEP3RD_BC, GBBEP3RD_OC, GBBEP3RD_PM25
      USE GFED_MOD,           ONLY : GFED3RD_BC, GFED3RD_OC, GFED3RD_PM25
      USE GFAS_MOD,           ONLY : GFAS3RD_BC, GFAS3RD_OC, GFAS3RD_PM25
      USE FINN_MOD,           ONLY : FINN3RD_BC, FINN3RD_OC, FINN3RD_PM25
      USE SEVIRI_MOD,         ONLY : SEVIRI3RD_BC, SEVIRI3RD_OC, SEVIRI3RD_PM25
      USE QFED_MOD,           ONLY : QFED3RD_BC, QFED3RD_OC, QFED3RD_PM25
      USE NAMELIST_MOD,       ONLY : F_PM25
      USE NAMELIST_ARRAY_MOD, ONLY : LNEI
      USE NAMELIST_ARRAY_MOD, ONLY : LINTEX
      USE NAMELIST_ARRAY_MOD, ONLY : LFLAMBE
      USE NAMELIST_ARRAY_MOD, ONLY : LFINN
      USE NAMELIST_ARRAY_MOD, ONLY : LGBBEP
      USE NAMELIST_ARRAY_MOD, ONLY : LGFED
      USE NAMELIST_ARRAY_MOD, ONLY : LGFAS
      USE NAMELIST_ARRAY_MOD, ONLY : LSEVIRI
      USE NAMELIST_ARRAY_MOD, ONLY : LQFED


      ! FORCE ALL VARIABLES TO BE DECLARED EXPLICITLY
      IMPLICIT NONE

      ! MAKE EVERYTHING PRIVATE ...
      PRIVATE

      ! ... EXCEPT THESE ROUTINES
      PUBLIC :: WRITE_RADM2

      ! ... AND THESE VARIABLES
      ! NONE

      !=================================================================
      ! MODULE ROUTINES -- FOLLOW BELOW THE 'CONTAINS' STATEMENT!
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
! $ID: WRITE_RADM2.F V01 06/16/2012 15:52 BRUCE EXP$
!
!******************************************************************************
!  SUBROUTINE WRITE_RADM2 WRITES THE EMISSION INTO THE FORMAT WHICH WRF-CHEM 
!  CAN READ IT DIRECTLY.
!
!  VARIABLES:
!  ============================================================================
!  (1 ) IDOM           (INTEGER) : # OF DOMAIN                         [---] 
!  (2 ) IYEAR          (INTEGER) : YEAR OF THE LOOP                    [---]
!  (3 ) IMONTH         (INTEGER) : MONTH OF THE LOOP                   [---]
!  (4 ) IDAY           (INTEGER) : DAY OF THE LOOP                     [---]
!  (5 ) STARTT         (INTEGER) : START TIME OF OUTPUT FILE           [---]
!  (6 ) ENDT           (INTEGER) : END TIEM OF OUTPUT FILE             [---]
!
!  NOTES:
!  ============================================================================
!  (1 ) ORIGINALLY WRITTEN BY BRUCE. (06/16/2012)
!  (2 ) MODIFED FROM write_sepa_v01.F. (06/16/2012)
!  (3 ) UNIT OF PM2.5 EMISSION IS ug/m3 m/s, FOR OTHER SPECIES, CAN REFER TO 
!       convert_emiss.exe OUTPUT FILE. (06/16/2012)
!  (4 ) OC : 90%, BC : 10%. (06/16/2012)
!  (5 ) REWRITTEN FROM SUBROUTINE WRITE_RADM2.F BY BRUCE. (07/03/2012)
!  (6 ) ADD DIFFERENT EMISSIONS WITH SWITCHES BY BRUCE. (07/10/2012)
!******************************************************************************

      SUBROUTINE WRITE_RADM2(IDOM, IYEAR, IMONTH, IDAY, STARTT, ENDT)

      ! REFERENCES TO F90 MODULES
      USE NAMELIST_ARRAY_MOD, ONLY : DOM_ID
      USE NAMELIST_ARRAY_MOD, ONLY : E_WE, E_SN, E_VERT
      USE NAMELIST_ARRAY_MOD, ONLY : SPECIES, N_SPECIES, P_SPECIES
      USE NAMELIST_ARRAY_MOD, ONLY : OUTPUTDIR
      USE NEI_MOD,            ONLY : NRADM, ENAME
      
      ! ARGUMENTS
      INTEGER                      :: IDOM, IYEAR, IMONTH, IDAY
      INTEGER                      :: STARTT, ENDT

      ! LOCAL VARIABLES
      INTEGER                      :: I, J, ILON, IHT, ILAT, ISPECIES
      INTEGER, DIMENSION(NRADM)    :: IMARK
      REAL*8,  DIMENSION(NRADM)    :: P_SPEC
      CHARACTER (LEN = 120)        :: CDAY, CDOMAIN, OUTPUTFILENAME
      CHARACTER (LEN = 120)        :: FPM25, BPM25
      CHARACTER (LEN = 120)        :: FOC,   BOC
      CHARACTER (LEN = 120)        :: FBC,   BBC
      REAL                         :: PERCENT
      REAL, DIMENSION(E_WE(IDOM)-1, E_VERT, E_SN(IDOM)-1) :: TMPEM

      BPM25   = 'e_pm25j'
      BOC     = 'e_orgj'
      BBC     = 'e_ecj'
      FPM25   = 'PM2p5'
      FOC     = 'OC'
      FBC     = 'BC'
      PERCENT = 0.01

      !=================================================================
      ! WRITE_RADM2 BEGINS HERE!
      !=================================================================

      ! SPECIFY DAY CHARACTER
      WRITE(CDAY,    100) IYEAR, '-', IMONTH, '-', IDAY
      WRITE(CDOMAIN, 120) 'd', DOM_ID(IDOM)
100   FORMAT(I4.4, A1, I2.2, A1, I2.2)
120   FORMAT(A1, I2.2)
130   FORMAT(A10, A100)

      !=================================================================
      ! NEI
      !=================================================================

      IF (LNEI(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI/'&
                             //CDAY(1: LEN_TRIM(CDAY))//              &
                             'wrfem_00to12z_'//                       &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSE IF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI/'&
                             //CDAY(1: LEN_TRIM(CDAY))//              &
                             'wrfem_12to24z_'//                       &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1,NRADM
               WRITE(19) NEI3RD(:, :, :, J, I)
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! INTEX
      !=================================================================

      IF (LINTEX(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'INTEX/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//   &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSE IF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'INTEX/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//   &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               WRITE(19)INTEX3RD(:, :, :, J, I)
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! to set index for species
      ! by zf
      !=================================================================
       IMARK = 0
       P_SPEC = 0.0
       DO ISPECIES = 1,N_SPECIES
          loop1: DO J = 1,NRADM
             IF ((SPECIES(ISPECIES) .EQ. FPM25 .AND. ENAME(J) .EQ. BPM25) .or. & ! PM2.5
                 (SPECIES(ISPECIES) .EQ. FOC   .AND. ENAME(J) .EQ. BOC)   .or. & ! OC
                 (SPECIES(ISPECIES) .EQ. FBC   .AND. ENAME(J) .EQ. BBC)) THEN    ! BC
                IMARK(J) = 1
                P_SPEC(J) = P_SPECIES(ISPECIES)
                exit loop1
             ENDIF
          ENDDO loop1
       ENDDO 

      !=================================================================
      ! NEI AND FLAMBE
      !=================================================================
       
      IF (LNEI(IDOM) .AND. LFLAMBE(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
          OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_FLAMBE/'// &
                           CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//       &
                           CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSE IF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
          OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_FLAMBE/'// &
                           CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//       & 
                           CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF
         
         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I 
            DO J = 1, NRADM 
               TMPEM = 0.0

               IF (IMARK(J) == 0) THEN
                  TMPEM(:,:,:) = NEI3RD(:,:,:,J,I)
               
               ELSE

               !==========================================================
               ! WRITE FLAMBE INTO NEI BACKGROUND SPECICES
               !==========================================================
                 
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 21, I) + FLAMBE3RD_PM25(:,:,:, I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 27, I) + FLAMBE3RD_OC(:,:,:, I) 

                  ! BC
                  ELSEIF(ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 29, I) + FLAMBE3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! NEI AND FINN
      !=================================================================

      IF (LNEI(IDOM) .AND. LFINN(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_FINN/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//      &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_FINN/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//      &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM   ! 30 aerosol species
               TMPEM = 0.0

               IF (IMARK(J) == 0) THEN
                  TMPEM(:,:,:) = NEI3RD(:,:,:,J,I)

               ELSE

               !==========================================================
               ! WRITE FINN INTO NEI BACKGROUND SPECICES
               !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 21, I) + FINN3RD_PM25(:,:,:, I) 

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 27, I) + FINN3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 29, I) + FINN3RD_BC(:,:,:, I) 

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF


      !=================================================================
      ! NEI AND GBBEP
      !=================================================================

      IF (LNEI(IDOM) .AND. LGBBEP(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_GBBEP/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//        &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_GBBEP/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//        &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) == 0) THEN
                  TMPEM(:,:,:) = NEI3RD(:,:,:,J,I)

               ELSE
                  !==========================================================
                  ! WRITE GBBEP INTO NEI BACKGROUND SPECICES
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 21, I) + GBBEP3RD_PM25(:,:,:, I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 27, I) + GBBEP3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 29, I) + GBBEP3RD_BC(:,:,:, I)
                  ENDIF
               ENDIF ! ILON
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF


      !=================================================================
      ! FLAMBE
      !=================================================================

      IF (LFLAMBE(IDOM)) THEN
       ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'FLAMBE/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//    &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'FLAMBE/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//    &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) then
               !==========================================================
               ! WRITE FLAMBE
               !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = FLAMBE3RD_PM25(:,:,:, I)

                  ! OC 
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = FLAMBE3RD_OC(:,:,:, I) 

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = FLAMBE3RD_BC(:,:,:, I) 

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF


      !=================================================================
      ! FINN
      !=================================================================

      IF (LFINN(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'FINN/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//  &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSE IF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'FINN/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//  &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) then
                  !==========================================================
                  ! WRITE FINN
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = FINN3RD_PM25(:,:,:, I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = FINN3RD_OC(:,:,:, I) 

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = FINN3RD_BC(:,:,:, I)
                     
                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF


      !=================================================================
      ! GBBEP
      !=================================================================

      IF (LGBBEP(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'GBBEP/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//   &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'GBBEP/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//   &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1,NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) then
                  !==========================================================
                  ! WRITE GBBEP
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = GBBEP3RD_PM25(:,:,:, I)

                  ! OC 
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = GBBEP3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = GBBEP3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF


      !=================================================================
      ! NEI AND GFED
      !=================================================================

      IF (LNEI(IDOM) .AND. LGFED(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_GFED/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//      &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_GFED/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//      &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) == 0) THEN
                  TMPEM(:,:,:) = NEI3RD(:,:,:,J,I)

               ELSE
                  !==========================================================
                  ! WRITE GFED INTO NEI BACKGROUND SPECICES
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 21, I) + GFED3RD_PM25(:,:,:,I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 27, I) + GFED3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 29, I) + GFED3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! NEI AND GFAS
      !=================================================================

      IF (LNEI(IDOM) .AND. LGFAS(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_GFAS/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'// &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'NEI_GFAS/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'// &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) == 0) THEN
                  TMPEM(:,:,:) = NEI3RD(:,:,:,J,I)

               ELSE
                  !==========================================================
                  ! WRITE GFAS INTO NEI BACKGROUND SPECICES
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 21, I) + GFAS3RD_PM25(:,:,:,I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 27, I) + GFAS3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = NEI3RD(:,:,:, 29, I) + GFAS3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! GFED
      !=================================================================

      IF (LGFED(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'GFED/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//   &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'GFED/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//   &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) THEN
                  !==========================================================
                  ! WRITE GFED
                  !==========================================================
                  ! PM2.5J
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = GFED3RD_PM25(:,:,:, I)
                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = GFED3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = GFED3RD_BC(:,:,:, I)
                     
                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! GFAS
      !=================================================================

      IF (LGFAS(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'GFAS/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'// &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'GFAS/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'// &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) THEN
                  !==========================================================
                  ! WRITE GFAS
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = GFAS3RD_PM25(:,:,:, I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = GFAS3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = GFAS3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! SEVIRI
      !=================================================================

      IF (LSEVIRI(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'SEVIRI/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'//    &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'SEVIRI/'// &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'//    &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) THEN
                  !==========================================================
                  ! WRITE SEVIRI
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = SEVIRI3RD_PM25(:,:,:, I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = SEVIRI3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = SEVIRI3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF

      !=================================================================
      ! QFED
      !=================================================================

      IF (LQFED(IDOM)) THEN
         ! DETERMINE FILE NAME DEPENDING ON STARTT AND ENDT
         IF (STARTT .EQ. 1 .AND. ENDT .EQ. 12) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'QFED/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_00to12z_'// &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ELSEIF (STARTT .EQ. 13 .AND. ENDT .EQ. 24) THEN
            OUTPUTFILENAME = OUTPUTDIR(1: LEN_TRIM(OUTPUTDIR))//'QFED/'//  &
                             CDAY(1: LEN_TRIM(CDAY))//'wrfem_12to24z_'// &
                             CDOMAIN(1: LEN_TRIM(CDOMAIN))
         ENDIF

         WRITE(6, 130) 'Writing : ', OUTPUTFILENAME
         OPEN(19, FILE = OUTPUTFILENAME, FORM = 'UNFORMATTED')
         WRITE(19) NRADM
         WRITE(19) ENAME

         DO I = STARTT, ENDT
            WRITE(19) I
            DO J = 1, NRADM
               TMPEM = 0.0

               IF (IMARK(J) /= 0) THEN
                  !==========================================================
                  ! WRITE QFED
                  !==========================================================
                  ! PM2.5I
                  IF (ENAME(J) .EQ. BPM25) THEN
                     TMPEM(:,:,:) = QFED3RD_PM25(:,:,:, I)

                  ! OC
                  ELSEIF (ENAME(J) .EQ. BOC) THEN
                     TMPEM(:,:,:) = QFED3RD_OC(:,:,:, I)

                  ! BC
                  ELSEIF (ENAME(J) .EQ. BBC) THEN
                     TMPEM(:,:,:) = QFED3RD_BC(:,:,:, I)

                  ENDIF
               ENDIF
               WRITE(19) TMPEM
            ENDDO ! J
         ENDDO ! I
         CLOSE(19)
      ENDIF
      ! RETURN TO THE CALLING PROGRAM
      END SUBROUTINE WRITE_RADM2


      END MODULE WRITE_MOD
