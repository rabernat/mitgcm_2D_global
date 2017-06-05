C $Header: 
C $Name: compute_global_mean

#include "CPP_OPTIONS.h"

CBOP 0
C !ROUTINE: FLT_LAVD

C !INTERFACE:

      SUBROUTINE COMPUTE_GLOBAL_MEAN_RL (
     I                myNr, arr, arrMask, arrArea,
     I                mean, bi, bj, 
     I                myTime, myThid )
C     !DESCRIPTION: \bv
C     /==========================================================\
C     | SUBROUTINE COMPUTE_GLOBAL_MEAN_RL                        |
C     | o Calculate the area wieghted global  mean of array      |
C     | "_RL arr" and returns it as mean                         |
C     |                                                          |
C     | o based on the code in remove_mean.F and mon_vort3.F     | 
C      \==========================================================/
C     \ev

      IMPLICIT NONE

C     === Global data ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     === Functions ====
      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     === Routine arguments ===
      INTEGER myNr
      _RL arr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS arrMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS arrArea(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      _RL mean
      _RL myTime
      INTEGER myThid, bi, bj 

C#ifdef ALLOW_BALANCE_FLUXES
C     === Local variables ====
      INTEGER I,J,K
      _RL tmpArea
      _RL theMean
      _RL theArea
      _RL tmpVal

      CHARACTER*(MAX_LEN_MBUF) msgBuf
      
      theMean=0.
      theArea=0.
      K = myNr



C     DO bj=myByLo(myThid),myByHi(myThid)
C       DO bi=myBxLo(myThid),myBxHi(myThid)
C        DO K=1,myNr
         DO J=1,sNy
          DO I=1,sNx
           tmpVal=arr(I,J,K,bi,bj)
           IF (arrMask(I,J,K,bi,bj).NE.0.) THEN
            tmpArea = arrArea(I,J,bi,bj)
            theArea = theArea   + tmpArea
            theMean = theMean + tmpArea*tmpVal
           ENDIF
          ENDDO
         ENDDO
C       ENDDO
C       ENDDO
C     ENDDO

C     output to help in debugging

        WRITE(msgBuf,'(A)')
     &    '----- DEBUG COMPUTE_GLOBAL_MEAN_RL -------- '
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT, myThid )

        WRITE(msgBuf,'(A,I,A,F)')
     &    'myThid= ', myThid, 'last temp area  =', tmpArea
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         
        WRITE(msgBuf,'(A,I,A,F)')
     &    'myThid= ', myThid, 'last temp vorticity  =', tmpVal
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )
         
        WRITE(msgBuf,'(A,I,A,F)')
     &    'myThid= ', myThid, 'thread total vorticity*area =', theMean
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                       SQUEEZE_RIGHT, myThid )

        WRITE(msgBuf,'(A,I,A,F)')
     &       'myThid= ', myThid, 'thread total area =', theArea
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT, myThid )

        WRITE(msgBuf,'(A)')
     &    '------ END DEBUG COMPUTE_GLOBAL_MEAN_RL -------- '
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &               SQUEEZE_RIGHT, myThid )
        
      _GLOBAL_SUM_RL(theArea,myThid)
      _GLOBAL_SUM_RL(theMean,myThid)

      IF (theArea.GT.0.) THEN
       theMean=theMean/theArea

      ENDIF

      mean = theMean 

      RETURN
      END

