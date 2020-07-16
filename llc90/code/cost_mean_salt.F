#include "COST_OPTIONS.h"
c#ifdef ALLOW_CTRL
c# include "CTRL_OPTIONS.h"
c#endif

      SUBROUTINE COST_MEAN_SALT( myThid )
C     *==========================================================*
C     | SUBROUTINE COST_TEMP
C     | o the subroutine computes the sum of the squared errors
C     |   relatively to the Levitus climatology
C     *==========================================================*

       IMPLICIT NONE

C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "cost.h"
#include "ctrl_weights.h"

C     ======== Routine arguments ======================
C     myThid - Thread number for this instance of the routine.
      INTEGER myThid

#ifdef ALLOW_COST_MEAN_SALT
C     ========= Local variables =========================
      INTEGER i, j, k
      INTEGER bi, bj
      INTEGER Nk
      _RL locfc,tmp
      _RL mask_arc_oc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      Nk = 2
C Read annual mean Levitus temperature

      CALL READ_FLD_XY_RL('mask_Arctic_Ocean_depth_mask_500.data',
     &                    ' ',mask_arc_oc,0,myThid)

C  Total number of wet temperature point
      tmp  = 0. _d 0
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
          DO j=1,sNy
           DO i=1,sNx
              IF ( mask_arc_oc(i,j,bi,bj).GT.0. ) THEN
                 DO k=1, Nk
                    IF ( cMeanSalt(i,j,k,bi,bj).LT.35. ) THEN
                       tmp = tmp + maskC(i,j,k,bi,bj)*
     &                      hFacC(i,j,k,bi,bj)*rA(i,j,bi,bj)
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
          ENDDO
       ENDDO
      ENDDO
      _GLOBAL_SUM_RL( tmp , myThid )
      IF ( tmp.GT.0. ) tmp = 1. _d 0 / tmp

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

          locfc = 0. _d 0
          DO j=1,sNy
           DO i=1,sNx
              IF ( mask_arc_oc(i,j,bi,bj).GT.0. ) THEN
                 DO k=1, Nk
                    print*, 'mean salinity: ',cMeanSalt(i,j,k,bi,bj)
                    IF ( cMeanSalt(i,j,k,bi,bj).LT.35. ) THEN
                       locfc = locfc + tmp*maskC(i,j,k,bi,bj)*
     &                      hFacC(i,j,k,bi,bj)*rA(i,j,bi,bj)*
     &                      cMeanSalt(i,j,k,bi,bj)
                       print*, 'salinity smaller 35'
                       print*, 'mean of tile:',locfc
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
          ENDDO

          objf_mean_salt(bi,bj) = locfc
c         print*,'objf_mean_salt =',locfc,startTime,endTime,tmp

       ENDDO
      ENDDO

#endif /* ALLOW_COST_MEAN_SALT */
      RETURN
      END
