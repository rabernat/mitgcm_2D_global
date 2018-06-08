C $Header: 1.0  2017/6/1 net $
C $Name: flt_lavd $

#include "FLT_OPTIONS.h"

CBOP 0
C !ROUTINE: FLT_LADIV

C !INTERFACE:
      SUBROUTINE FLT_LADIV (
     I                      myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE FLT_LADIV
C     | o This routine computes the Lagrangian-averaged horiztontal divergence of
C     | the velocity field for each particle 
C     *==========================================================*

C     !USES:
      IMPLICIT NONE
      
C     == global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "FLT_SIZE.h"
#include "FLT.h"
#include "FLT_BUFF.h"
#include "GRID.h"
#ifdef ALLOW_EXCH2
#include "W2_EXCH2_SIZE.h"
#include "W2_EXCH2_TOPOLOGY.h"
#endif

C     === Routie Arguments ===
C     myTime :: current time in simulation
C     myIter :: current iteration number
C     myThid :: my Thread Id number
      _RL myTime
      INTEGER myIter, myThid

C     === Local Variables === 
      INTEGER bi, bj
      INTEGER ip, kp
      INTEGER i, j
      _RL ix, jy 
      
C     divergence on float position
      _RL float_div
C     horizontal divergence (Eulerian field) 
      _RL divUV(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr, nSx, nSy)
      _RL hFacZ(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RL r_hFacZ(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

      CHARACTER*(MAX_LEN_MBUF) msgBuf

CEOP



C--   loop over all the subgrids and compute the horizontal divergence field

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
         IF ( flt_selectTrajOutp.GE.3 ) THEN
C --     Calculate horizontal divergence
          DO kp = 1,Nr
             CALL MOM_CALC_HFACZ( bi,bj,kp,hFacZ,r_hFacZ,myThid )
             CALL MOM_CALC_HDIV( bi, bj, kp, 1,
     &                           uVel, vVel,
     &                           divUV(1-OLx,1-OLy,kp, bi, bj),
     &                           myThid)
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
                IF ( hFacZ(i,j).EQ.0. ) THEN
                  divUV(i,j,kp,bi,bj)  = 0.D0
                ENDIF
              ENDDO
             ENDDO
           ENDDO
         ENDIF
       ENDDO
      ENDDO      

C--   loop over all subgrids and interpolate div field to find ladiv for each flt
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO ip=1,npart_tile(bi,bj)
            ix = ipart(ip,bi,bj)
            jy = jpart(ip,bi,bj)
            kp = NINT(kpart(ip,bi,bj))

            CALL FLT_BILINEAR(ix,jy,float_div,divUV, kp, 4,bi,bj,myThid)

            ladiv(ip,bi,bj) = (float_div*flt_deltaT +
     &      (myTime - flt_deltaT)*ladiv(ip,bi,bj) )/myTime
        ENDDO
       ENDDO
      ENDDO

      RETURN

      END
