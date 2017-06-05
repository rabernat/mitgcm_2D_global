C $Header: 1.0  2017/6/1 net $
C $Name: flt_lavd $

#include "FLT_OPTIONS.h"

CBOP 0
C !ROUTINE: FLT_LAVD

C !INTERFACE:
      SUBROUTINE FLT_LAVD (
     I                      myTime, myIter, myThid )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE FLT_LAVD
C     | o This routine computes the LAVD field for each particle 
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
      _RL ix, jy, i0x, j0y
      
C     relative vorticity on float position
      _RL vo
C     relative vorticity (Eulerian field)
      _RL vort3(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, Nr)
      _RS hFacZ(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
      _RS r_hFacZ(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      
C     the global mean vort and inst. vorticity deviation 
      _RL mean_vort 
      _RL ivd

      CHARACTER*(MAX_LEN_MBUF) msgBuf

#ifdef ALLOW_EXCH2
      INTEGER nT
#endif
CEOP

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)

#ifdef ALLOW_EXCH2
         nT = W2_myTileList(bi,bj)
         i0x = DFLOAT( exch2_txGlobalo(nT) - 1 )
         j0y = DFLOAT( exch2_tyGlobalo(nT) - 1 )
#else
         i0x = DFLOAT( myXGlobalLo-1 + (bi-1)*sNx )
         j0y = DFLOAT( myYGlobalLo-1 + (bj-1)*sNy )
#endif
         IF ( flt_selectTrajOutp.GE.3 ) THEN
C --     Calculate relative vorticity
           DO kp = 1,Nr
             CALL MOM_CALC_HFACZ( bi,bj,kp,hFacZ,r_hFacZ,myThid )
             CALL MOM_CALC_RELVORT3( bi, bj, kp,
     &                                 uVel, vVel, hFacZ,
     &                                 vort3(1-OLx,1-OLy,kp),
     &                                 myThid)
             DO j=1-OLy,sNy+OLy
              DO i=1-OLx,sNx+OLx
                IF ( hFacZ(i,j).EQ.0. ) THEN
                  vort3(i,j,kp)  = 0. _d 0
                ENDIF
              ENDDO
             ENDDO
           ENDDO
         ENDIF

C use new subroutine to compute area-weighted average
C     SUBROUTINE COMPUTE_GLOBAL_MEAN_RL(
C     I                myNr, arr, arrMask, arrArea,
C     I                mean, bi, bj,
C     I                ,myTime, myThid)
         
         CALL COMPUTE_GLOBAL_MEAN_RL(kp, vort3, maskInC, rA,
     &        mean_vort, 1, 1, myTime, myThid)

C output global mean to STDOUT for reference          
         WRITE(msgBuf,'(A,F)')
     &        'global mean vorticity =', mean_vort
         CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
        
         DO ip=1,npart_tile(bi,bj)
            ix = ipart(ip,bi,bj)
            jy = jpart(ip,bi,bj)
            kp = NINT(kpart(ip,bi,bj))
              
C-- RPA+NET hack: calculate relative vorticity and lavd
            CALL FLT_BILINEAR(ix,jy,vo,vort3, kp, 4,1,1,myThid)
            
C-- mean_vort is currently buggy, so we set to zero as tmp fix 
            mean_vort = 0
            ivd = ABS(vo-mean_vort) 

            lavd(ip,bi,bj) = ( ivd*flt_deltaT +
     &      (myTime - flt_deltaT)*lavd(ip,bi,bj) )/myTime
         ENDDO

       ENDDO
      ENDDO

      RETURN

      END
