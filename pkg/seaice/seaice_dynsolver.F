#include "SEAICE_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

CBOP
C     !ROUTINE: SEAICE_DYNSOLVER
C     !INTERFACE:
      SUBROUTINE SEAICE_DYNSOLVER( myTime, myIter, myThid )

C     *==========================================================*
C     | SUBROUTINE SEAICE_DYNSOLVER
C     | o Ice dynamics using LSR solver
C     |   Zhang and Hibler,   JGR, 102, 8691-8702, 1997
C     |   or EVP explicit solver by Hunke and Dukowicz, JPO 27,
C     |   1849-1867 (1997)
C     *==========================================================*
C     | written by Martin Losch, March 2006
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SEAICE_SIZE.h"
#include "SEAICE_PARAMS.h"
#include "SEAICE.h"

#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myTime     :: Simulation time
C     myIter     :: Simulation timestep number
C     myThid     :: my Thread Id. number
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
CEOP

#ifdef SEAICE_CGRID
C     !FUNCTIONS:
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE
#ifdef ALLOW_DIAGNOSTICS
      LOGICAL  DIAGNOSTICS_IS_ON
      EXTERNAL DIAGNOSTICS_IS_ON
#endif /* ALLOW_DIAGNOSTICS */

C     !LOCAL VARIABLES:
C     === In Local common block ===
C     TAUX   :: zonal      wind stress over seaice at U point
C     TAUY   :: meridional wind stress over seaice at V point
      COMMON/WIND_STRESS_ICE/TAUX,TAUY
      _RL TAUX   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL TAUY   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     === Local variables ===
C     i,j    :: Loop counters
C     bi,bj  :: tile counters
      INTEGER i, j, bi, bj
      _RL phiSurf(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL mask_uice, mask_vice

C     INTEGER l, m, n, lMin, lMax, mMin, mMax
      INTEGER iMin, iMax, jMin, jMax, n
      _RL tmparea (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpheff (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpmask (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpwght

#ifdef ALLOW_DIAGNOSTICS
# ifdef ALLOW_AUTODIFF
      _RL strDivX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL strDivY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
# endif /* ALLOW_AUTODIFF */
      _RL sig1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig2   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sig12  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigp, sigm, sigTmp, recip_prs
      _RL areaW, areaS, COSWAT
      _RS SINWAT
      INTEGER kSrf
      LOGICAL diag_SIsigma_isOn, diag_SIshear_isOn
      LOGICAL diag_SIenpi_isOn, diag_SIenpot_isOn
      LOGICAL diag_SIpRfric_isOn, diag_SIpSfric_isOn
#endif /* ALLOW_DIAGNOSTICS */


C Initialise local variables in common block
      IF ( myIter.EQ.nIter0 ) THEN
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           TAUX(i,j,bi,bj) = 0. _d 0
           TAUY(i,j,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF

# ifdef ALLOW_AUTODIFF
C Following re-initialisation breaks some "artificial" AD dependencies
C incured by IF (DIFFERENT_MULTIPLE ... statement
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          PRESS0(i,j,bi,bj) = SEAICE_strength*HEFF(i,j,bi,bj)
     &         *EXP(-SEAICE_cStar*(ONE-AREA(i,j,bi,bj)))
          ZMAX(i,j,bi,bj)   = SEAICE_zetaMaxFac*PRESS0(i,j,bi,bj)
          ZMIN(i,j,bi,bj)   = SEAICE_zetaMin
          PRESS0(i,j,bi,bj) = PRESS0(i,j,bi,bj)*HEFFM(i,j,bi,bj)
          TAUX(i,j,bi,bj)   = 0. _d 0
          TAUY(i,j,bi,bj)   = 0. _d 0
# ifdef SEAICE_ALLOW_FREEDRIFT
          uice_fd(i,j,bi,bj)= 0. _d 0
          vice_fd(i,j,bi,bj)= 0. _d 0
# endif
# ifdef ALLOW_DIAGNOSTICS
          strDivX(i,j,bi,bj)= 0. _d 0
          strDivY(i,j,bi,bj)= 0. _d 0
# endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO
# endif /* ALLOW_AUTODIFF */
# ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE uice    = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vice    = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE uicenm1 = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vicenm1 = comlev1, key=ikey_dynamics, kind=isbyte
# endif /* ALLOW_AUTODIFF_TAMC */

      IF (
     &  DIFFERENT_MULTIPLE(SEAICE_deltaTdyn,myTime,SEAICE_deltaTtherm)
     &   ) THEN

# ifdef ALLOW_AUTODIFF_TAMC
# ifdef SEAICE_ALLOW_EVP
CADJ STORE press0 = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE zmax   = comlev1, key=ikey_dynamics, kind=isbyte
# endif
# endif /* ALLOW_AUTODIFF_TAMC */

C--   NOW SET UP MASS PER UNIT AREA AND CORIOLIS TERM
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy+1,sNy+OLy
         DO i=1-OLx+1,sNx+OLx
          seaiceMassC(i,j,bi,bj)=SEAICE_rhoIce*HEFF(i,j,bi,bj)
          seaiceMassU(i,j,bi,bj)=SEAICE_rhoIce*HALF*(
     &          HEFF(i,j,bi,bj) + HEFF(i-1,j  ,bi,bj) )
          seaiceMassV(i,j,bi,bj)=SEAICE_rhoIce*HALF*(
     &          HEFF(i,j,bi,bj) + HEFF(i  ,j-1,bi,bj) )
         ENDDO
        ENDDO
        IF ( SEAICEaddSnowMass ) THEN
         DO j=1-OLy+1,sNy+OLy
          DO i=1-OLx+1,sNx+OLx
           seaiceMassC(i,j,bi,bj)=seaiceMassC(i,j,bi,bj)
     &          +                 SEAICE_rhoSnow*HSNOW(i,j,bi,bj)
           seaiceMassU(i,j,bi,bj)=seaiceMassU(i,j,bi,bj)
     &         +                  SEAICE_rhoSnow*HALF*(
     &          HSNOW(i,j,bi,bj) + HSNOW(i-1,j  ,bi,bj) )

           seaiceMassV(i,j,bi,bj)=seaiceMassV(i,j,bi,bj)
     &         +                  SEAICE_rhoSnow*HALF*(
     &          HSNOW(i,j,bi,bj) + HSNOW(i  ,j-1,bi,bj) )
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

#ifndef ALLOW_AUTODIFF
      IF ( SEAICE_maskRHS ) THEN
C     dynamic masking of areas with no ice, not recommended
C     and only kept for testing purposes
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy+1,sNy+OLy
          DO i=1-OLx+1,sNx+OLx
           seaiceMaskU(i,j,bi,bj)=AREA(i,j,bi,bj)+AREA(i-1,j,bi,bj)
           mask_uice=HEFFM(i,j,bi,bj)+HEFFM(i-1,j  ,bi,bj)
           IF ( (seaiceMaskU(i,j,bi,bj) .GT. 0. _d 0) .AND.
     &         (mask_uice .GT. 1.5 _d 0) ) THEN
            seaiceMaskU(i,j,bi,bj) = 1. _d 0
           ELSE
            seaiceMaskU(i,j,bi,bj) = 0. _d 0
           ENDIF
           seaiceMaskV(i,j,bi,bj)=AREA(i,j,bi,bj)+AREA(i,j-1,bi,bj)
           mask_vice=HEFFM(i,j,bi,bj)+HEFFM(i  ,j-1,bi,bj)
           IF ( (seaiceMaskV(i,j,bi,bj) .GT. 0. _d 0) .AND.
     &         (mask_vice .GT. 1.5 _d 0) ) THEN
            seaiceMaskV(i,j,bi,bj) = 1. _d 0
           ELSE
            seaiceMaskV(i,j,bi,bj) = 0. _d 0
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       CALL EXCH_UV_XY_RL( seaiceMaskU, seaiceMaskV, .FALSE., myThid )
      ENDIF
#endif /* ndef ALLOW_AUTODIFF */

C--   NOW SET UP FORCING FIELDS

C     initialise fields
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          TAUX (i,j,bi,bj)= 0. _d 0
          TAUY (i,j,bi,bj)= 0. _d 0
#ifdef ALLOW_AUTODIFF
# ifdef SEAICE_ALLOW_EVP
          stressDivergenceX(i,j,bi,bj) = 0. _d 0
          stressDivergenceY(i,j,bi,bj) = 0. _d 0
# endif
#endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO

# ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE uice = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vice = comlev1, key=ikey_dynamics, kind=isbyte
# endif /* ALLOW_AUTODIFF_TAMC */
C--   interface of dynamics with atmopheric forcing fields (wind/stress)
      CALL SEAICE_GET_DYNFORCING (
     I     uIce, vIce, AREA, SIMaskU, SIMaskV,
     O     TAUX, TAUY,
     I     myTime, myIter, myThid )

CNH   Smoothening of concentration and thickness fields for the ice strenght computation
C     loop boundaries
      iMin=1-OLx
      iMax=sNx+OLx
      jMin=1-OLy
      jMax=sNy+OLy

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO j=jMin,jMax
         DO i=iMin,iMax
          AREAMN(i,j,bi,bj)=AREA(i,j,bi,bj)
          HEFFMN(i,j,bi,bj)=HEFF(i,j,bi,bj)
          tmparea(i,j,bi,bj)= 0.0
          tmpheff(i,j,bi,bj)= 0.0
          tmpmask(i,j,bi,bj)= 0.0
          IF (AREA(i,j,bi,bj).GT.0) THEN
           tmpmask(i,j,bi,bj) = 1.0
          ENDIF
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      
      DO n=1,220
c$$$       DO bj=myByLo(myThid),myByHi(myThid)
c$$$        DO bi=myBxLo(myThid),myBxHi(myThid)
c$$$         DO j=jMin,jMax
c$$$          DO i=iMin,iMax
c$$$           tmpscal3(i,j,bi,bj)=0. _d 0
c$$$           tmpscal4(i,j,bi,bj)=0. _d 0
c$$$           tmpscal5=0. _d 0
c$$$           lMin=MAX(iMin-i,-4)
c$$$           lMax=MIN(iMax-i,4)
c$$$           mMin=MAX(jMin-j,-4)
c$$$           mMax=MIN(jMax-j,4)
c$$$C          print*,'loop parameters i',bi,bj,i,lMin,lMax,iMin,iMax
c$$$C          print*,'loop parameters j',bi,bj,j,mMin,mMax,jMin,jMax
c$$$           IF (AREA(i,j,bi,bj).GT.0) THEN
c$$$            DO l=lMin,lMax
c$$$             DO m=mMin,mMax
c$$$              tmpscal3(i,j,bi,bj)=tmpscal3(i,j,bi,bj)+
c$$$     &           HEFF(i+l,j+m,bi,bj)
c$$$              tmpscal4(i,j,bi,bj)=tmpscal4(i,j,bi,bj)+
c$$$     &           AREA(i+l,j+m,bi,bj)
c$$$              tmpscal5=tmpscal5 + 1. _d 0
c$$$C             IF (i.EQ.1) THEN
c$$$C              IF (j.EQ.1) THEN
c$$$C               print*,'seaice loop thickness',tmpscal3,
c$$$C     &          HEFF(i+l,j+m,bi,bj),l,m
c$$$C               print*,'seaice loop concentration',tmpscal4,
c$$$C     &          AREA(i+l,j+m,bi,bj),l,m
c$$$C              ENDIF
c$$$C             ENDIF
c$$$             ENDDO     
c$$$            ENDDO
c$$$            tmpscal3(i,j,bi,bj)=tmpscal3(i,j,bi,bj)/tmpscal5
c$$$            tmpscal4(i,j,bi,bj)=tmpscal4(i,j,bi,bj)/tmpscal5
c$$$            tmpscal3(i,j,bi,bj)=MAX(tmpscal3(i,j,bi,bj),HEFF(i,j,bi,bj))
c$$$            tmpscal4(i,j,bi,bj)=MAX(tmpscal4(i,j,bi,bj),AREA(i,j,bi,bj))
c$$$            IF (tmpscal4(i,j,bi,bj).GT.1) THEN
c$$$             print*,'WARNING Averaging lead to area lager 1',i,j,
c$$$     &            tmpscal4,tmpscal5
c$$$            ENDIF
c$$$           ENDIF
c$$$          ENDDO
c$$$         ENDDO
c$$$        ENDDO
c$$$       ENDDO
c$$$
c$$$       DO bj=myByLo(myThid),myByHi(myThid)
c$$$        DO bi=myBxLo(myThid),myBxHi(myThid)
c$$$         DO j=jMin,jMax
c$$$          DO i=iMin,iMax
c$$$           AREAMN(i,j,bi,bj)=tmpscal4(i,j,bi,bj)
c$$$           HEFFMN(i,j,bi,bj)=tmpscal3(i,j,bi,bj)
c$$$          ENDDO
c$$$         ENDDO
c$$$        ENDDO
c$$$       ENDDO

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=0,sNy
          DO i=0,sNx
           IF (AREA(i,j,bi,bj).GT.0) THEN
            tmpwght=0.16210282*tmpmask(i,j,bi,bj)+
     &            0.09832033*(tmpmask(i-1,j,bi,bj)+tmpmask(i,j-1,bi,bj)+
     &            tmpmask(i,j+1,bi,bj)+tmpmask(i+1,j,bi,bj))+
     &            0.0596343*(tmpmask(i-1,j-1,bi,bj)+
     &            tmpmask(i+1,j-1,bi,bj)+
     &            tmpmask(i+1,j+1,bi,bj)+tmpmask(i-1,j-1,bi,bj))+
     &            0.02193823*(tmpmask(i,j-2,bi,bj)+tmpmask(i+2,j,bi,bj)+
     &            tmpmask(i,j+2,bi,bj)+tmpmask(i-2,j,bi,bj))+
     &            0.00296902*(tmpmask(i-2,j-2,bi,bj)+
     &            tmpmask(i+2,j-2,bi,bj)+
     &            tmpmask(i-2,j+2,bi,bj)+tmpmask(i+2,j+2,bi,bj))+
     &            0.01330621*(tmpmask(i-2,j-1,bi,bj)+
     &            tmpmask(i-2,j+1,bi,bj)+
     &            tmpmask(i+2,j+1,bi,bj)+tmpmask(i+2,j-1,bi,bj)+
     &            tmpmask(i-1,j+2,bi,bj)+tmpmask(i+1,j+2,bi,bj)+
     &            tmpmask(i-1,j-2,bi,bj)+tmpmask(i+1,j-2,bi,bj))  
            tmparea(i,j,bi,bj)=0.16210282*AREAMN(i,j,bi,bj)+
     &             0.09832033*(AREAMN(i-1,j,bi,bj)+AREAMN(i,j-1,bi,bj)+
     &             AREAMN(i,j+1,bi,bj)+AREAMN(i+1,j,bi,bj))+
     &             0.0596343*(AREAMN(i-1,j-1,bi,bj)+
     &             AREAMN(i+1,j-1,bi,bj)+
     &             AREAMN(i+1,j+1,bi,bj)+AREAMN(i-1,j-1,bi,bj))+
     &             0.02193823*(AREAMN(i,j-2,bi,bj)+AREAMN(i+2,j,bi,bj)+
     &             AREAMN(i,j+2,bi,bj)+AREAMN(i-2,j,bi,bj))+
     &             0.00296902*(AREAMN(i-2,j-2,bi,bj)+
     &             AREAMN(i+2,j-2,bi,bj)+
     &             AREAMN(i-2,j+2,bi,bj)+AREAMN(i+2,j+2,bi,bj))+
     &             0.01330621*(AREAMN(i-2,j-1,bi,bj)+
     &             AREAMN(i-2,j+1,bi,bj)+
     &             AREAMN(i+2,j+1,bi,bj)+AREAMN(i+2,j-1,bi,bj)+
     &             AREAMN(i-1,j+2,bi,bj)+AREAMN(i+1,j+2,bi,bj)+
     &             AREAMN(i-1,j-2,bi,bj)+AREAMN(i+1,j-2,bi,bj))            
            tmpheff(i,j,bi,bj)=0.16210282*HEFFMN(i,j,bi,bj)+
     &             0.09832033*(HEFFMN(i-1,j,bi,bj)+HEFFMN(i,j-1,bi,bj)+
     &             HEFFMN(i,j+1,bi,bj)+HEFFMN(i+1,j,bi,bj))+
     &             0.0596343*(HEFFMN(i-1,j-1,bi,bj)+
     &             HEFFMN(i+1,j-1,bi,bj)+
     &             HEFFMN(i+1,j+1,bi,bj)+HEFFMN(i-1,j-1,bi,bj))+
     &             0.02193823*(HEFFMN(i,j-2,bi,bj)+HEFFMN(i+2,j,bi,bj)+
     &             HEFFMN(i,j+2,bi,bj)+HEFFMN(i-2,j,bi,bj))+
     &             0.00296902*(HEFFMN(i-2,j-2,bi,bj)+
     &             HEFFMN(i+2,j-2,bi,bj)+
     &             HEFFMN(i-2,j+2,bi,bj)+HEFFMN(i+2,j+2,bi,bj))+
     &             0.01330621*(HEFFMN(i-2,j-1,bi,bj)+
     &             HEFFMN(i-2,j+1,bi,bj)+
     &             HEFFMN(i+2,j+1,bi,bj)+HEFFMN(i+2,j-1,bi,bj)+
     &             HEFFMN(i-1,j+2,bi,bj)+HEFFMN(i+1,j+2,bi,bj)+
     &             HEFFMN(i-1,j-2,bi,bj)+HEFFMN(i+1,j-2,bi,bj))            
C            HEFFMN(i,j,bi,bj)=0.5*HEFFMN(i,j,bi,bj)+0.5*
C     &             MAX(HEFFMN(i,j,bi,bj),HEFF(i,j,bi,bj))
C            AREAMN(i,j,bi,bj)=0.5*AREAMN(i,j,bi,bj)+0.5*
C     &             MAX(AREAMN(i,j,bi,bj),AREA(i,j,bi,bj))
C     AREAMN(i,j,bi,bj)=MIN(AREAMN(i,j,bi,bj),1.0)
            tmpheff(i,j,bi,bj)=tmpheff(i,j,bi,bj)/tmpwght
            tmparea(i,j,bi,bj)=tmparea(i,j,bi,bj)/tmpwght
C            IF ( (j.EQ.1) .AND. (i.EQ.1) ) THEN 
CC            IF (tmparea(i,j,bi,bj).GT.1) THEN
C             print*,tmparea(i,j,bi,bj),AREAMN(i,j,bi,bj),AREA(i,j,bi,bj)
CC            ENDIF
C            ENDIF
            IF (tmpheff(i,j,bi,bj).GT.5.0) THEN
             tmpheff(i,j,bi,bj)=MAX(tmpheff(i,j,bi,bj),HEFF(i,j,bi,bj))
             tmparea(i,j,bi,bj)=MAX(tmparea(i,j,bi,bj),AREA(i,j,bi,bj))
            ENDIF
            tmparea(i,j,bi,bj)=MIN(tmparea(i,j,bi,bj),1.0)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=0,sNy
          DO i=0,sNx
           AREAMN(i,j,bi,bj)=tmparea(i,j,bi,bj)
           HEFFMN(i,j,bi,bj)=tmpheff(i,j,bi,bj)  
          ENDDO
         ENDDO
        ENDDO
       ENDDO

       _EXCH_XY_RL( HEFFMN,  myThid )
       _EXCH_XY_RL( AREAMN,  myThid )
      ENDDO
CNH   End of concentration and thickness smoothening      

      

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
C--   Compute surface pressure at z==0:
C-    use actual sea surface height for tilt computations
        DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            phiSurf(i,j) = Bo_surf(i,j,bi,bj)*etaN(i,j,bi,bj)
          ENDDO
        ENDDO
#ifdef ATMOSPHERIC_LOADING
C-    add atmospheric loading and Sea-Ice loading
        IF ( useRealFreshWaterFlux ) THEN
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            phiSurf(i,j) = phiSurf(i,j)
     &                   + ( pload(i,j,bi,bj)
     &                      +sIceLoad(i,j,bi,bj)*gravity
     &                     )*recip_rhoConst
           ENDDO
          ENDDO
        ELSE
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            phiSurf(i,j) = phiSurf(i,j)
     &                   + pload(i,j,bi,bj)*recip_rhoConst
           ENDDO
          ENDDO
        ENDIF
#endif /* ATMOSPHERIC_LOADING */
C--   basic forcing by wind stress
        IF ( SEAICEscaleSurfStress ) THEN
         DO j=1-OLy+1,sNy+OLy
          DO i=1-OLx+1,sNx+OLx
           FORCEX0(i,j,bi,bj)=TAUX(i,j,bi,bj)
     &          * 0.5 _d 0*(AREA(i,j,bi,bj)+AREA(i-1,j,bi,bj))
           FORCEY0(i,j,bi,bj)=TAUY(i,j,bi,bj)
     &          * 0.5 _d 0*(AREA(i,j,bi,bj)+AREA(i,j-1,bi,bj))
          ENDDO
         ENDDO
        ELSE
         DO j=1-OLy+1,sNy+OLy
          DO i=1-OLx+1,sNx+OLx
           FORCEX0(i,j,bi,bj)=TAUX(i,j,bi,bj)
           FORCEY0(i,j,bi,bj)=TAUY(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDIF

        IF ( SEAICEuseTILT ) then
        DO j=1-OLy+1,sNy+OLy
         DO i=1-OLx+1,sNx+OLx
C--   now add in tilt
          FORCEX0(i,j,bi,bj)=FORCEX0(i,j,bi,bj)
     &         -seaiceMassU(i,j,bi,bj)*_recip_dxC(i,j,bi,bj)
     &         *( phiSurf(i,j)-phiSurf(i-1,j) )
          FORCEY0(i,j,bi,bj)=FORCEY0(i,j,bi,bj)
     &         -seaiceMassV(i,j,bi,bj)* _recip_dyC(i,j,bi,bj)
     &         *( phiSurf(i,j)-phiSurf(i,j-1) )
         ENDDO
        ENDDO
        ENDIF

        CALL SEAICE_CALC_ICE_STRENGTH( bi, bj, myTime, myIter, myThid )

       ENDDO
      ENDDO

#ifdef SEAICE_ALLOW_DYNAMICS
      IF ( SEAICEuseDYNAMICS ) THEN
#ifdef SEAICE_ALLOW_FREEDRIFT
       IF ( SEAICEuseFREEDRIFT .OR. SEAICEuseEVP
     &                         .OR. LSR_mixIniGuess.EQ.0 ) THEN
        CALL SEAICE_FREEDRIFT( myTime, myIter, myThid )
       ENDIF
       IF ( SEAICEuseFREEDRIFT ) THEN
        DO bj=myByLo(myThid),myByHi(myThid)
         DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            uIce(i,j,bi,bj) = uIce_fd(i,j,bi,bj)
            vIce(i,j,bi,bj) = vIce_fd(i,j,bi,bj)
            stressDivergenceX(i,j,bi,bj) = 0. _d 0
            stressDivergenceY(i,j,bi,bj) = 0. _d 0
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDIF
#endif /* SEAICE_ALLOW_FREEDRIFT */

#ifdef ALLOW_OBCS
       IF ( useOBCS ) THEN
         CALL OBCS_APPLY_UVICE( uIce, vIce, myThid )
       ENDIF
#endif /* ALLOW_OBCS */

#ifdef SEAICE_ALLOW_EVP
# ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE uice    = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vice    = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE uicenm1 = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vicenm1 = comlev1, key=ikey_dynamics, kind=isbyte
# endif /* ALLOW_AUTODIFF_TAMC */
       IF ( SEAICEuseEVP ) THEN
C     Elastic-Viscous-Plastic solver, following Hunke (2001)
        CALL SEAICE_EVP( myTime, myIter, myThid )
       ENDIF
#endif /* SEAICE_ALLOW_EVP */

       IF ( SEAICEuseLSR ) THEN
C     Picard solver with LSR scheme (Zhang-J/Hibler 1997), ported to a C-grid
        CALL SEAICE_LSR( myTime, myIter, myThid )
       ENDIF

#ifdef SEAICE_ALLOW_KRYLOV
# ifdef ALLOW_AUTODIFF
       STOP 'Adjoint does not work with Picard-Krylov solver.'
# else
       IF ( SEAICEuseKrylov ) THEN
C     Picard solver with Matrix-free Krylov solver (Lemieux et al. 2008)
        CALL SEAICE_KRYLOV( myTime, myIter, myThid )
       ENDIF
# endif /*  ALLOW_AUTODIFF */
#endif /* SEAICE_ALLOW_KRYLOV */

#ifdef SEAICE_ALLOW_JFNK
# ifdef ALLOW_AUTODIFF
       STOP 'Adjoint does not work with JFNK solver.'
# else
       IF ( SEAICEuseJFNK ) THEN
C     Jacobian-free Newton Krylov solver (Lemieux et al. 2010, 2012)
        CALL SEAICE_JFNK( myTime, myIter, myThid )
       ENDIF
# endif /*  ALLOW_AUTODIFF */
#endif /* SEAICE_ALLOW_JFNK */

C End of IF (SEAICEuseDYNAMICS ...
      ENDIF
#endif /* SEAICE_ALLOW_DYNAMICS */

C End of IF (DIFFERENT_MULTIPLE ...
      ENDIF

#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE uice  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vice  = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE stressDivergenceX  = comlev1,
CADJ &     key=ikey_dynamics, kind=isbyte
CADJ STORE stressDivergenceY  = comlev1,
CADJ &     key=ikey_dynamics, kind=isbyte
CADJ STORE DWATN  = comlev1, key=ikey_dynamics, kind=isbyte
#ifdef SEAICE_ALLOW_BOTTOMDRAG
CADJ STORE CbotC  = comlev1, key=ikey_dynamics, kind=isbyte
#endif /* SEAICE_ALLOW_BOTTOMDRAG */
#endif /* ALLOW_AUTODIFF_TAMC */

C Update ocean surface stress
      IF ( SEAICEupdateOceanStress ) THEN
        CALL SEAICE_OCEAN_STRESS (
     I              TAUX, TAUY, myTime, myIter, myThid )
      ENDIF

#ifdef SEAICE_ALLOW_DYNAMICS
#ifdef SEAICE_ALLOW_CLIPVELS
      IF ( SEAICEuseDYNAMICS .AND. SEAICE_clipVelocities) THEN
#ifdef ALLOW_AUTODIFF_TAMC
CADJ STORE uice = comlev1, key=ikey_dynamics, kind=isbyte
CADJ STORE vice = comlev1, key=ikey_dynamics, kind=isbyte
#endif /* ALLOW_AUTODIFF_TAMC */
c Put a cap on ice velocity
c limit velocity to 0.40 m s-1 to avoid potential CFL violations
c in open water areas (drift of zero thickness ice)
       DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           uIce(i,j,bi,bj)=
     &          MAX(MIN(uIce(i,j,bi,bj),0.40 _d +00),-0.40 _d +00)
           vIce(i,j,bi,bj)=
     &          MAX(MIN(vIce(i,j,bi,bj),0.40 _d +00),-0.40 _d +00)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif /* SEAICE_ALLOW_CLIPVELS */
#endif /* SEAICE_ALLOW_DYNAMICS */

#ifdef ALLOW_DIAGNOSTICS
C     diagnostics related to mechanics/dynamics/momentum equations
      IF ( useDiagnostics .AND. SEAICEuseDYNAMICS ) THEN
       CALL DIAGNOSTICS_FILL(zeta   ,'SIzeta  ',0,1,0,1,1,myThid)
       CALL DIAGNOSTICS_FILL(eta    ,'SIeta   ',0,1,0,1,1,myThid)
       CALL DIAGNOSTICS_FILL(press  ,'SIpress ',0,1,0,1,1,myThid)
       CALL DIAGNOSTICS_FILL(deltaC ,'SIdelta ',0,1,0,1,1,myThid)
       IF ( DIAGNOSTICS_IS_ON('SItensil',myThid) ) THEN
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
C     use sig1 as a temporary field
          DO j=1,sNy
           DO i=1,sNx
            IF ( tensileStrFac(i,j,bi,bj) .EQ. oneRL ) THEN
C     This special case of tensile strength equal to compressive strength
C     is not very physical and should actually not happen but you never know;
C     in this case, press = P-T = P*(1-k) = 0. and we have to use press0 to
C     get something
             sig1(i,j) = press0(i,j,bi,bj)
            ELSE
C     This is more complicated than you think because press = P-T = P*(1-k),
C     but we are looking for T = k*P = k*press/(1-k)
             sig1(i,j) = tensileStrFac(i,j,bi,bj)
     &            *press(i,j,bi,bj)/(1. _d 0 - tensileStrFac(i,j,bi,bj))
            ENDIF
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig1,'SItensil',0,1,2,bi,bj,myThid)
         ENDDO
        ENDDO
       ENDIF
C     If any of the stress or energy diagnostics are required,
C     first recompute strainrates from up-to-date velocities
       diag_SIsigma_isOn  = DIAGNOSTICS_IS_ON('SIsig1  ',myThid)
     &                 .OR. DIAGNOSTICS_IS_ON('SIsig2  ',myThid)
       diag_SIshear_isOn  = DIAGNOSTICS_IS_ON('SIshear ',myThid)
       diag_SIenpi_isOn   = DIAGNOSTICS_IS_ON('SIenpi  ',myThid)
       diag_SIenpot_isOn  = DIAGNOSTICS_IS_ON('SIenpot ',myThid)
       diag_SIpRfric_isOn = DIAGNOSTICS_IS_ON('SIpRfric',myThid)
       diag_SIpSfric_isOn = DIAGNOSTICS_IS_ON('SIpSfric',myThid)
       IF ( diag_SIsigma_isOn  .OR. diag_SIshear_isOn .OR.
     &      diag_SIenpi_isOn   .OR. diag_SIenpot_isOn .OR.
     &      diag_SIpRfric_isOn .OR. diag_SIpSfric_isOn ) THEN
        CALL SEAICE_CALC_STRAINRATES(
     I       uIce, vIce,
     O       e11, e22, e12,
     I       0, myTime, myIter, myThid )
C     but use old viscosities and pressure for the
C     principle stress components
CML     CALL SEAICE_CALC_VISCOSITIES(
CML  I     e11, e22, e12, zMin, zMax, HEFFM, press0, tensileStrFac,
CML  O     eta, etaZ, zeta, zetaZ, press, deltaC,
CML  I     0, myTime, myIter, myThid )
       ENDIF
C
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
C
C     stress diagnostics
C
         IF ( diag_SIsigma_isOn ) THEN
#ifdef SEAICE_ALLOW_EVP
C     This could go directly into EVP, but to keep a better eye on it,
C     I would like to keep it here.
          IF ( SEAICEuseEVP ) THEN
C     for EVP compute principle stress components from recent
C     stress state and normalize with latest
C     PRESS = PRESS(n-1), n = number of sub-cycling steps
           DO j=1,sNy
            DO i=1,sNx
             sigp = seaice_sigma1(i,j,bi,bj)
             sigm = seaice_sigma2(i,j,bi,bj)
             sig12(i,j) = 0.25 _d 0 *
     &            ( seaice_sigma12(i  ,j  ,bi,bj)
     &            + seaice_sigma12(i+1,j  ,bi,bj)
     &            + seaice_sigma12(i+1,j+1,bi,bj)
     &            + seaice_sigma12(i  ,j+1,bi,bj) )
             sigTmp = SQRT( sigm*sigm + 4. _d 0*sig12(i,j)*sig12(i,j) )
             recip_prs = 0. _d 0
             IF ( press0(i,j,bi,bj) .GT. 1. _d -13 )
     &            recip_prs = 1. _d 0 / press0(i,j,bi,bj)
             sig1(i,j) = 0.5 _d 0*(sigp + sigTmp)*recip_prs
             sig2(i,j) = 0.5 _d 0*(sigp - sigTmp)*recip_prs
            ENDDO
           ENDDO
          ELSE
#else
          IF ( .TRUE. ) THEN
#endif /* SEAICE_ALLOW_EVP */
           CALL SEAICE_CALC_STRESS(
     I          e11, e22, e12, press, zeta, eta, etaZ,
     O          sig1, sig2, sig12,
     I          bi, bj, myTime, myIter, myThid )
           DO j=1,sNy
            DO i=1,sNx
             sigp   = sig1(i,j) + sig2(i,j)
             sigm   = sig1(i,j) - sig2(i,j)
C     This should be the way of computing sig12 at C-points,
C            sigTmp = 0.25 _d 0 *
C    &            ( sig12(i  ,j  ) + sig12(i+1,j  )
C    &            + sig12(i  ,j+1) + sig12(i+1,j+1) )
C     but sig12 = 2*etaZ*e12, and because of strong gradients in eta,
C     etaZ can be very large for a cell with small eta and the straightforward
C     way of averaging mixes large etaZ with small press0, so we have to do it
C     in different way to get meaningfull sig12C (=sigTmp):
             sigTmp = 2. _d 0 * eta(i,j,bi,bj) * 0.25 _d 0 *
     &            (e12(i,j,bi,bj) + e12(i+1,j,bi,bj)
     &            +e12(i,j+1,bi,bj)+e12(i+1,j+1,bi,bj))
             sigTmp = SQRT( sigm*sigm + 4. _d 0*sigTmp*sigTmp )
             recip_prs = 0. _d 0
             IF ( press0(i,j,bi,bj) .GT. 1. _d -13 )
     &            recip_prs = 1. _d 0 / press0(i,j,bi,bj)
             sig1(i,j) = 0.5 _d 0*(sigp + sigTmp)*recip_prs
             sig2(i,j) = 0.5 _d 0*(sigp - sigTmp)*recip_prs
            ENDDO
           ENDDO
          ENDIF
          CALL DIAGNOSTICS_FILL(sig1,'SIsig1  ',0,1,2,bi,bj,myThid)
          CALL DIAGNOSTICS_FILL(sig2,'SIsig2  ',0,1,2,bi,bj,myThid)
         ENDIF
C
         IF ( diag_SIshear_isOn  ) THEN
          DO j=1,sNy
           DO i=1,sNx
            sigm = e11(i,j,bi,bj) - e22(i,j,bi,bj)
            sigTmp =
     &           ( e12(i  ,j  ,bi,bj)**2 + e12(i+1,j  ,bi,bj)**2
     &           + e12(i+1,j+1,bi,bj)**2 + e12(i  ,j+1,bi,bj)**2 )
C     shear deformation as sqrt((e11-e22)**2 + 4*e12**2); the 4 pulled into
C     the average
            sig1(i,j) = sqrt(sigm*sigm + sigTmp)
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig1,'SIshear ',0,1,2,bi,bj,myThid)
         ENDIF
C
C     most of the energy diagnostics, re-use sig1, sig2, sig12 as tmp-arrays
C
         IF ( diag_SIenpi_isOn ) THEN
C     compute internal stresses with updated ice velocities
C     TAF gets confused when we use stressDivergenceX/Y as temporary arrays
C     therefore we need to use new arrays
          IF ( .NOT. SEAICEuseFREEDRIFT )
     &         CALL SEAICE_CALC_STRESSDIV(
     I         e11, e22, e12, press, zeta, eta, etaZ,
# ifdef ALLOW_AUTODIFF
     O         strDivX, strDivY,
# else /* not ALLOW_AUTODIFF */
     O         stressDivergenceX, stressDivergenceY,
# endif /* ALLOW_AUTODIFF */
     I         bi, bj, myTime, myIter, myThid )
          DO j=1,sNy+1
           DO i=1,sNx+1
# ifdef ALLOW_AUTODIFF
            sig1(i,j) = uIce(i,j,bi,bj)*strDivX(i,j,bi,bj)
            sig2(i,j) = vIce(i,j,bi,bj)*strDivY(i,j,bi,bj)
# else /* not ALLOW_AUTODIFF */
            sig1(i,j) = uIce(i,j,bi,bj)*stressDivergenceX(i,j,bi,bj)
            sig2(i,j) = vIce(i,j,bi,bj)*stressDivergenceY(i,j,bi,bj)
# endif /* ALLOW_AUTODIFF */
           ENDDO
          ENDDO
C     average from velocity points to pressure points
          DO j=1,sNy
           DO i=1,sNx
            sig12(i,j)= 0.5 _d 0 * ( sig1(i,j) + sig1(i+1,j)
     &                             + sig2(i,j) + sig2(i,j+1) )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig12,'SIenpi  ',0,1,2,bi,bj,myThid)
         ENDIF
         IF ( diag_SIenpot_isOn ) THEN
          DO j=1,sNy
           DO i=1,sNx
            sig1(i,j) = 0.5 _d 0 * press0(i,j,bi,bj) *
     &           ( e11(i,j,bi,bj) + e22(i,j,bi,bj) )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig1,'SIenpot ',0,1,2,bi,bj,myThid)
         ENDIF
         IF ( diag_SIpRfric_isOn ) THEN
          DO j=1,sNy
           DO i=1,sNx
            sig1(i,j) = - zeta(i,j,bi,bj)
     &           * ( e11(i,j,bi,bj) + e22(i,j,bi,bj) )
     &           * ( e11(i,j,bi,bj) + e22(i,j,bi,bj) )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig1,'SIpRfric',0,1,2,bi,bj,myThid)
         ENDIF
         IF ( diag_SIpSfric_isOn ) THEN
          DO j=1,sNy
           DO i=1,sNx
            sig1(i,j) = - eta(i,j,bi,bj)
     &        * ( (e11(i,j,bi,bj) - e22(i,j,bi,bj))**2 +
     &            ( e12(i  ,j  ,bi,bj)**2 + e12(i+1,j  ,bi,bj)**2
     &            + e12(i+1,j+1,bi,bj)**2 + e12(i  ,j+1,bi,bj)**2 )
     &          )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig1,'SIpSfric',0,1,2,bi,bj,myThid)
         ENDIF
         IF ( DIAGNOSTICS_IS_ON('SIenpa  ',myThid) ) THEN
          DO j=1,sNy+1
           DO i=1,sNx+1
            areaW = 0.5 _d 0 * (AREA(i,j,bi,bj) + AREA(i-1,j,bi,bj))
     &           * SEAICEstressFactor
            areaS = 0.5 _d 0 * (AREA(i,j,bi,bj) + AREA(i,j-1,bi,bj))
     &           * SEAICEstressFactor
            sig1(i,j) = TAUX(i,j,bi,bj)*areaW * uIce(i,j,bi,bj)
            sig2(i,j) = TAUY(i,j,bi,bj)*areaS * vIce(i,j,bi,bj)
           ENDDO
          ENDDO
C     average from velocity points to pressure points
          DO j=1,sNy
           DO i=1,sNx
            sig12(i,j)= 0.5 _d 0 * ( sig1(i,j) + sig1(i+1,j)
     &                             + sig2(i,j) + sig2(i,j+1) )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig12,'SIenpa  ',0,1,2,bi,bj,myThid)
         ENDIF
         IF ( DIAGNOSTICS_IS_ON('SIenpw  ',myThid) ) THEN
C     surface level
          kSrf = 1
C     introduce turning angle (default is zero)
          SINWAT=SIN(SEAICE_waterTurnAngle*deg2rad)
          COSWAT=COS(SEAICE_waterTurnAngle*deg2rad)
          DO j=1,sNy+1
           DO i=1,sNx+1
            areaW = 0.5 _d 0 * (AREA(i,j,bi,bj) + AREA(i-1,j,bi,bj))
     &         * SEAICEstressFactor
            areaS = 0.5 _d 0 * (AREA(i,j,bi,bj) + AREA(i,j-1,bi,bj))
     &         * SEAICEstressFactor
            sig1(i,j) = areaW *
     &           HALF*( DWATN(i,j,bi,bj)+DWATN(i-1,j,bi,bj) )*
     &         COSWAT *
     &         ( uIce(i,j,bi,bj)-uVel(i,j,kSrf,bi,bj) )
     &         - SIGN(SINWAT, _fCori(i,j,bi,bj)) * 0.5 _d 0 *
     &         ( DWATN(i  ,j,bi,bj) *
     &         0.5 _d 0*(vIce(i  ,j  ,bi,bj)-vVel(i  ,j  ,kSrf,bi,bj)
     &                  +vIce(i  ,j+1,bi,bj)-vVel(i  ,j+1,kSrf,bi,bj))
     &         + DWATN(i-1,j,bi,bj) *
     &         0.5 _d 0*(vIce(i-1,j  ,bi,bj)-vVel(i-1,j  ,kSrf,bi,bj)
     &                  +vIce(i-1,j+1,bi,bj)-vVel(i-1,j+1,kSrf,bi,bj))
     &         )
            sig1(i,j) = sig1(i,j) * uIce(i,j,bi,bj)
            sig2(i,j) = areaS *
     &           HALF*( DWATN(i,j,bi,bj)+DWATN(i,j-1,bi,bj) )*
     &         COSWAT *
     &         ( vIce(i,j,bi,bj)-vVel(i,j,kSrf,bi,bj) )
     &         + SIGN(SINWAT,  _fCori(i,j,bi,bj)) * 0.5 _d 0 *
     &         ( DWATN(i,j  ,bi,bj) *
     &         0.5 _d 0*(uIce(i  ,j  ,bi,bj)-uVel(i  ,j  ,kSrf,bi,bj)
     &                  +uIce(i+1,j  ,bi,bj)-uVel(i+1,j  ,kSrf,bi,bj))
     &         + DWATN(i,j-1,bi,bj) *
     &         0.5 _d 0*(uIce(i  ,j-1,bi,bj)-uVel(i  ,j-1,kSrf,bi,bj)
     &                  +uIce(i+1,j-1,bi,bj)-uVel(i+1,j-1,kSrf,bi,bj))
     &         )
            sig2(i,j) = sig2(i,j) * vIce(i,j,bi,bj)
           ENDDO
          ENDDO
C     average from velocity points to pressure points
          DO j=1,sNy
           DO i=1,sNx
            sig12(i,j)= - 0.5 _d 0 * ( sig1(i,j) + sig1(i+1,j)
     &                               + sig2(i,j) + sig2(i,j+1) )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig12,'SIenpw  ',0,1,2,bi,bj,myThid)
         ENDIF
         IF ( SEAICEuseTilt
     &        .AND. DIAGNOSTICS_IS_ON('SIenpg  ',myThid) ) THEN
          DO j=1-OLy+1,sNy+OLy
           DO i=1-OLx+1,sNx+OLx
            sig1(i,j)=
     &           - seaiceMassU(i,j,bi,bj)*_recip_dxC(i,j,bi,bj)
     &           *( phiSurf(i,j)-phiSurf(i-1,j) ) * uIce(i,j,bi,bj)
            sig2(i,j)=
     &           - seaiceMassV(i,j,bi,bj)* _recip_dyC(i,j,bi,bj)
     &           *( phiSurf(i,j)-phiSurf(i,j-1) ) * vIce(i,j,bi,bj)
           ENDDO
          ENDDO
C     average from velocity points to pressure points
          DO j=1,sNy
           DO i=1,sNx
            sig12(i,j) = 0.5 _d 0 * ( sig1(i,j) + sig1(i+1,j)
     &                              + sig2(i,j) + sig2(i,j+1) )
           ENDDO
          ENDDO
          CALL DIAGNOSTICS_FILL(sig12,'SIenpg  ',0,1,2,bi,bj,myThid)
         ENDIF
C     bi/bj-loop
        ENDDO
       ENDDO
C     useDiagnostics
      ENDIF

#endif   /* ALLOW_DIAGNOSTICS */
#endif  /* SEAICE_CGRID */

      RETURN
      END
