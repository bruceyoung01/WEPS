!  --------------------------------------------------------------------
!  $ID: tool_kit_mod.f
!  1/17/10
!  
!  INTERFACE:

      MODULE TOOL_KIT_MOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC   :: JD
      PUBLIC   :: ANNUAL_JD
      PUBLIC   :: JD_to_MMDD
      PUBLIC   :: DISTANCE
      PUBLIC   :: ZP_GRID
      PUBLIC   :: SCATTERING_ANGLE

      CONTAINS 

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      INTEGER FUNCTION JD (YEAR,MONTH,DAY)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).
!
      INTEGER YEAR,MONTH,DAY,I,J,K

      I= YEAR
      J= MONTH
      K= DAY

      JD= K-32075+1461*(I+4800+(J-14)/12)/4+367*(J-2-(J-14)/12*12) &
          /12-3*((I+4900+(J-14)/12)/100)/4

      RETURN
      END

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      INTEGER  FUNCTION ANNUAL_JD( YEAR,MONTH,DAY )

!  calculate julian day within one year

!  arguments

      INTEGER,INTENT( IN) :: YEAR,MONTH,DAY

!  local variables

      INTEGER             :: JD
      INTEGER             :: leap_day
      INTEGER             :: i

!  Check for leap year, and add extra day if necessary

      IF ( MOD(year,400) == 0 ) THEN
       leap_day = 1
      ELSE IF ( MOD(year,100) == 0 ) THEN
       leap_day = 0
      ELSE IF ( MOD(year,4) == 0 ) THEN
       leap_day = 1
      ELSE
       leap_day = 0
      ENDIF

! calculate julian day number

      IF ( month == 1 ) THEN
       JD = day
      ELSE  
       JD = day
       DO i = 1, month-1
        SELECT CASE(i)
        CASE(1,3,5,7,8,10,12)
          JD = JD + 31
        CASE(4,6,9,11)
          JD = JD + 30
        CASE(2)
          JD = JD + 28 + leap_day
        END SELECT
       ENDDO
      ENDIF
  
      annual_jd = JD
  
      RETURN

      END FUNCTION ANNUAL_JD

!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE GDATE (JD, YEAR,MONTH,DAY)
!
!---COMPUTES THE GREGORIAN CALENDAR DATE (YEAR,MONTH,DAY)
!   GIVEN THE JULIAN DATE (JD).
!
      INTEGER JD,YEAR,MONTH,DAY,I,J,K
!
      L= JD+68569
      N= 4*L/146097
      L= L-(146097*N+3)/4
      I= 4000*(L+1)/1461001
      L= L-1461*I/4+31
      J= 80*L/2447
      K= L-2447*J/80
      L= J/11
      J= J+2-12*L
      I= 100*(N-49)+I+L
!
      YEAR= I
      MONTH= J
      DAY= K
!
      RETURN
      END

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  retrun date given annule Julian day

      SUBROUTINE JD_to_MMDD (JD,YEAR,MONTH,DAY)

!  arguments

      INTEGER,INTENT(IN)  :: JD,YEAR
      INTEGER,INTENT(OUT) :: MONTH,DAY

      INTEGER :: JD_MON(13)

      INTEGER :: leap_day,I

!  Check for leap year, and add extra day if necessary

      IF ( MOD(year,400) == 0 ) THEN
       leap_day = 1
      ELSE IF ( MOD(year,100) == 0 ) THEN
       leap_day = 0
      ELSE IF ( MOD(year,4) == 0 ) THEN
       leap_day = 1
      ELSE
       leap_day = 0
      ENDIF

!  Start JD for each month

      JD_MON = (/0,31,     &
             59+leap_day,  &
             90+leap_day,  &
            120+leap_day,  &
            151+leap_day,  &
            181+leap_day,  &
            212+leap_day,  &
            243+leap_day,  &
            273+leap_day,  &
            304+leap_day,  & 
            334+leap_day,  & 
            365+leap_day /)

      DO I = 1, 12
       IF(JD > JD_MON(I) .AND. JD <= JD_MON(I+1)) THEN
        MONTH = I
        DAY   = JD - JD_MON(I)
       ENDIF
      ENDDO

      RETURN

      END SUBROUTINE JD_to_MMDD

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!   Given two points on the surface of a sphere, given by longitude
!   and latitude pairs, computes the linear distance along the surface 
!   of the sphere between the two points (i.e. the Great Circle distance). 
!   (xxu, 1/18/10)
! 

      REAL FUNCTION DISTANCE (in_pt1_lon, in_pt1_lat,  &
                              in_pt2_lon, in_pt2_lat,  &
                              in_sph_radius )

!  arguments

      REAL,INTENT(IN)          :: in_pt1_lon
      REAL,INTENT(IN)          :: in_pt1_lat
      REAL,INTENT(IN)          :: in_pt2_lon
      REAL,INTENT(IN)          :: in_pt2_lat
      REAL,INTENT(IN),OPTIONAL :: in_sph_radius

!  local arrays

      REAL*8                   :: pt1_lon
      REAL*8                   :: pt1_lat
      REAL*8                   :: pt2_lon
      REAL*8                   :: pt2_lat
      REAL*8                   :: sph_radius

      REAL*8,PARAMETER         :: DPI = 3.14159265358979d0
      REAL*8,PARAMETER         :: dradeg = 180.d0/DPI 
      REAL*8                   :: phi1, phi2, theta1, theta2
      REAL*8                   :: x1, x2, y1, y2, z1, z2
      REAL*8                   :: dot_prod, angle_center_arg
      REAL*8                   :: angle_center, dist_surf

!  using double precison 

      pt1_lon = in_pt1_lon
      pt1_lat = in_pt1_lat
      pt2_lon = in_pt2_lon
      pt2_lat = in_pt2_lat

      IF (PRESENT(in_sph_radius)) THEN 
       sph_radius = in_sph_radius
      ELSE
       sph_radius = 6.371d6
      ENDIF

!  calc. spherical coord.: 
!    + phi is pi-latitude
!    + theta is longitude

      phi1 = (DPI/2.d0) - (pt1_lat / dradeg)
      phi2 = (DPI/2.d0) - (pt2_lat / dradeg)

      theta1 = pt1_lon / dradeg
      theta2 = pt2_lon / dradeg

!  spherical to cartesian coord

      x1 = sph_radius * SIN(phi1) * COS(theta1)
      x2 = sph_radius * SIN(phi2) * COS(theta2)

      y1 = sph_radius * SIN(phi1) * SIN(theta1)
      y2 = sph_radius * SIN(phi2) * SIN(theta2)

      z1 = sph_radius * COS(phi1)
      z2 = sph_radius * COS(phi2)

      dot_prod         = (x1*x2) + (y1*y2) + (z1*z2)
      angle_center_arg = dot_prod / (sph_radius**2)

!  test to make sure argument for ACOS is within range;
!  this accts. for case pt1 & pt2 are the same

      IF ( ABS(angle_center_arg) > 1.00001) PRINT *, '*** error-bad arguments...'

      IF (angle_center_arg >  1.0) angle_center_arg =  1.d0
      IF (angle_center_arg < -1.0) angle_center_arg = -1.d0

!  angle at center of sphere

      angle_center = ACOS(angle_center_arg)

!  Great Circle distance

      dist_surf    = angle_center * sph_radius 

      DISTANCE = dist_surf

      RETURN
       
      END FUNCTION DISTANCE

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
      SUBROUTINE ZP_GRID(PS,ILAYERS,ZGRID,PGRID)

!  arguments

      REAL,   INTENT(IN)  :: PS
      INTEGER,INTENT(IN)  :: ILAYERS
      REAL,   INTENT(OUT) :: ZGRID(ILAYERS+1)
      REAL,   INTENT(OUT) :: PGRID(ILAYERS+1)

!  local variables

      INTEGER             :: K, IEDGES

      IEDGES = ILAYERS + 1

      ZGRID(1) = 7.3*LOG(1013.)  ! TOP OF THE ATMOSPHERE - 50.5km
      PGRID(1) =  .010001

      DO K = 2, IEDGES
         ZGRID(K) = 7.3*LOG(1013./((PS-1.)*(REAL(K-1)/REAL(ILAYERS))))

         PGRID(K) = PS * ( 1. / ILAYERS)*(K - 1)

         IF (ZGRID(K) .LT. 0.0 ) ZGRID(K) = 0.0
      ENDDO

      END SUBROUTINE ZP_GRID

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      REAL*8 FUNCTION SCATTERING_ANGLE(THETA_SUN, THETA_SAT, AZIMUTH) 

      ! Arguments
      REAL*8,  INTENT(IN) :: THETA_SUN
      REAL*8,  INTENT(IN) :: THETA_SAT
      REAL*8,  INTENT(IN) :: AZIMUTH

      ! Local Variables
      REAL*8,PARAMETER    :: DPIE = 3.1415927d0
      REAL*8              :: COS_THETA
      REAL*8              :: THETA

      ! Cosine of scattering angle
      COS_THETA = COS(THETA_SUN*DPIE/180d0) * COS(THETA_SAT*DPIE/180d0)  &
                + SIN(THETA_SUN*DPIE/180d0) * SIN(THETA_SAT*DPIE/180d0)  &
                * COS(AZIMUTH*DPIE/180.d0)

      THETA = 180.d0 - ACOS(COS_THETA)*180.d0/DPIE
      SCATTERING_ANGLE = THETA
       
      RETURN

      END FUNCTION SCATTERING_ANGLE

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      END MODULE TOOL_KIT_MOD


