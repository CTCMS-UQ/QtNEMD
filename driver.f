C     3-DIMENSIONAL WCA FOR A BOUNDARY DRIVEN SHEARING SYSTEM OF
C     PARITICLES OR MOLECULES  
C     
C     Last changes: OCTOBER 2010
C     Based on program kawascii
C
C     NEMD transient trajectory simulation 
C     with a constant colour field applied to the waslls after equilibration
C     Colour field equations - 4th order Runge-Kutta
C     Peculiar internal energy or kinetic energy
c     thermostatted with Gaussian thermostat
C     Periodic boundary conditions in the x & y direction
C     WCA potential
C     FENE potential for intramolecular forces
C     Included NCELL AND ABILITY TO ALTER CUBE LENGTH IN Z DIRECTION
C
C    (i) Equilibration is carried out
C    (ii) A constant colour field to the walls is then applied and 
C         various properties
C        and tcf calculated. Some tau-averaged properties are
C        stored to disk for tau=MAXTAU
C    (ii) The time-correlation function is integrated (if NTYPE.gt.1).
C
C     The accumulated properties from are stored in MEAN 
C   and the number of data accumulated in NMEAN 
C     The accumulated autocorrelation functions are stored in TCF(I,0..NTCF)
C   and the number of data accumulated in KOUNT 
C     The integrated autocorrelation function is stored in INTCF(I,0..NTCF) 
C     The integral of properties along a trajectory is stored in TAUAV(0..NTCF)
C     The average of properties at points on a trajectory is stored in
C   AVER(I,0..NTCF)
C
C-----------------------------------------------------------------------
C This file contains a "driver" subroutine to be called from python 
C which handles user input, sets up calculation parameters and 
C launches the MD subroutines. 
C-----------------------------------------------------------------------
C
      SUBROUTINE SETUP()

      USE OMP_LIB
      INCLUDE "SSGK.inc"

C
C***** LOCAL VARIABLES
C
      REAL (KIND=sp) TEST,R6,R12,R4,RSI,TOTPX,TOTPY,TOTPZ
      REAL (KIND=sp) AR1PX,AR1PY,AR1PZ,XX,XX2,RX,RXEQ
      REAL (KIND=sp) AR2PX,AR2PY,AR2PZ,TESTL,TESTH
      REAL (KIND=sp) KENER
      INTEGER I,J,K,M,N,IXX
      INTEGER RANK, NUM_PROCESSORS, IERR
      CHARACTER*9 DATNOW
      REAL (KIND=sp) TIMNOW,SECNDS,RZERO

C
C***** DISPLAY TITLE
C
C      WRITE(6,'(//''NEMD WITH COLOUR: RK4,LJ'')')
C      CALL DATE(DATNOW)
      RZERO=0.0E0
C      TIMNOW=SECNDS(RZERO)
C      WRITE(6,*) DATNOW,TIMNOW
      CONSER=0.0D0
      E00=0.0D0
      KB=0
      TOTPX=0.0D0
      TOTPY=0.0D0
      TOTPZ=0.0D0
      AR1PX=0.0D0
      AR1PY=0.0D0
      AR1PZ=0.0D0
      AR2PX=0.0D0
      AR2PY=0.0D0
      AR2PZ=0.0D0

      FIELD=0.0D0

      DO_PRESSURE = .FALSE.
C
C***** ANALYSE INPUT
C
       IF (MAXTAU.GT.NK) THEN
          WRITE(6,*)'MAXTAU TOO HIGH FOR DIMENSION:',NK
          STOP
       END IF
       IF (NTYPE.LT.1.OR.NTYPE.GT.3) THEN
          WRITE(*,*) NTYPE
          WRITE(6,*)'NO SUCH NTYPE IN THIS PROGRAM, CHOOSE:1,2,3'
          STOP
       END IF
C
       TEST=DBLE(NPART)
       IF(NTYPE.EQ.1) THEN
        IF(LATT.EQ.1)THEN
            IF (MOD(TEST,4.0D0).NE.0.0D0) THEN
               WRITE(6,*)'NPART,MUST BE DIVISIBLE BY 4'
               STOP
            ELSE
               TEST=(NPART/4.0D0*YZDIVX)**(1.0D0/3.0D0)
               TESTL=TEST-0.0000001
               TESTH=TEST+0.0000001
              IF (ANINT(TEST).LT.TESTL.OR.ANINT(TEST).GT.TESTH) THEN
                 WRITE(6,*)'NPART MUST BE 4*N^3/YZDIVX; N INTEGER'
                 WRITE(6,*)TESTL,TESTH
                 STOP
              ENDIF
            ENDIF
        ENDIF 
C C
       TEST=DBLE(NPART*YZDIVX)
         IF(LATT.EQ.2)THEN
           IF (MOD(TEST,2.0D0).NE.0.0D0) THEN
              WRITE(6,*)'NPART,MUST BE AN EVEN NO'
              STOP
           ELSE
              TEST=(NPART*YZDIVX/2.0D0)**(1.0D0/3.0D0)
              TESTL=TEST-0.0000001
              TESTH=TEST+0.0000001

             IF (ANINT(TEST).LT.TESTL.OR.ANINT(TEST).GT.TESTH) THEN
                WRITE(6,*)'NPART MUST BE 2*N^3, WHERE N IS AN INTEGER'
                WRITE(6,*)TESTL,TESTH
                STOP
             END IF
           END IF
          END IF
       END IF

C*****DEFINE NPART2
       NPART2 = NPART/2
C***** ASSIGN MOLECULE NO TO FLUID PARTICLES
    
      DO J=1,NPART2
         IMOL(J)=1
      ENDDO

      DO J=NPART2+1, NPART 
         IMOL(J)=2
      ENDDO
C       
C C***** CHECK ALL FLUID PARTICLES ARE PART OF A MOLECULE
C        IF (MOD(NFLUID,LIMOL).NE.0.0D0) THEN
C           WRITE(6,*)'NFLUID,MUST BE DIVISIBLE BY LIMOL NO REMAINDER'
C           STOP
C        END IF
C
C***** ADD CHARGE TO WALLS PARTICLES ONLY
C
      DO I =1,NPART
         IKOL(I)= (-1)**(I)
      ENDDO
      DO I =1,NPART2
         IMOL(I)= 1
      ENDDO
       DO I =NPART2+1,NPART
         IMOL(I)= -1
      ENDDO

      IF (NGAUS.EQ.6) THEN
         DF=3.0D0*DBLE(NPART)-3.0D0
      ELSE 
         DF=3.0D0*DBLE(NPART)-4.0D0 
      ENDIF
C
C***** ZERO ARRAYS
C
      DO I = 0, NK
          DO J= 0, BINSX
            DENSAV1(I,J) = 0.0D0
            DENSAV2(I,J) = 0.0D0
              DO K = 1, NTCF
                TCFL(I,J,K) = 0.0D0
                TCFLL(I,J,K) = 0.0D0
                INTCFL(I,J,K) = 0.0D0
                INTCFLL(I,J,K) = 0.0D0
              ENDDO    
          ENDDO       
      ENDDO

      DO 1 I = 0, NK
         KOUNT(I) = 0
         DO 2 J = 1, NTCF
            AVER(I,J)=0.0D0
            TCF(I,J) = 0.0D0
            INTCF(I,J) = 0.0D0
            TCF1(I,J) = 0.0D0
            INTCF1(I,J) = 0.0D0
 2       CONTINUE
 1    CONTINUE
      NMEAN = 0
      DO 3 I=1,NPROP
         MEAN(I) = 0.0D0
 3    CONTINUE
      DO 44 I=1,NTCF
         M0(I)=0.0D0
         M1(I)=0.0D0
         M2(I)=0.0D0
         M3(I)=0.0D0
 44   CONTINUE
C
C***** SET INITIAL VALUES
C      RCUT SET TO 2**(1/6) FOR WCA POTENTIAL
C      SEED set as unix time 
        ISEED = 10101
C       ISEED = TIME()
#ifdef DEBUG
       WRITE(6,*) 'SEED', ISEED  
#endif
C
C***** CALCULATE RUN PARAMETERS
C
      IF (NTYPE.EQ.1) THEN
         VOL = NPART/DRF
      ELSE 
         VOL = NPART/DRF
      END IF
C      write(6,*)'VOLF,VOL,DRW,DRF',VOLF,VOL,DRW,DRF
C***** CALCULATE LENGTHS BASED ON YZDIVX RATIO 
      CUBEZ = (YZDIVX*VOL)**(1.0D0/3.0D0)
      CUBEY = CUBEZ
      CUBEX = CUBEZ/YZDIVX
      CUBEX2 = CUBEX/2.0D0
      CUBEY2 = CUBEY/2.0D0                                     
      CUBEZ2 = CUBEZ/2.0D0
      BINS = FLOOR(CUBEX)
      TAUBINSX = FLOOR(CUBEX)
C
C****** IF RCUT TOO LARGE RESET TO HALF BOXLENGTH
C      MIX applies mixing rules  based on EPS1 and EPS2 
      IF (RCUT.GT.CUBEY2) THEN
         RCUT = CUBEX2
         WRITE(6,*)'RCUT RESET TO HALF BOXLENGTH',RCUT
      END IF
      RMAX  = RCUT**2
      RSI    = 1.D0/RMAX
      R4     = RSI*RSI
      R6     = R4*RSI
      R12    = R6*R6
      SHIFT  = 4.0D0*(R12-R6)
C
C***** WRITE SIMULATION PARAMETERS
C
#ifdef DEBUG
      WRITE(6,150) NPART
      WRITE(6,151) TR,E0,DRW,DRF,RCUT,FE0,NTYPE,
     &             DELTA,NON,RCUT,NGAUS,CUBE,MAXTAU,
     &             EQTIM,NCYC,NLAYER,NLP,KH,LATT,NPRINT,
     &             KF,R0,LIMOL,YZDIVX,CUBEX,CUBEY,CUBEZ,
     &             DXXDIV
  150 FORMAT(/,10X,'SIMULATION PARAMETERS',I6,' PARTICLES')
  151 FORMAT(/,5X,'TR      =',F14.7,10X,'E/N     =',F14.7,10X,
     &       /,5X,'WALL DENSITY =',F14.7,10X,'FL DEN=',F14.7,10X,
     &       /,5X,'RCUT    =',F14.7,
     &       /,5X,'FE0   =',F14.7,10X,'NTYPE  =',I6,
     &       /,5X,'DELTA  =',F14.7,10X,'NON    =',I6,
     &       /,5X,'RCUT   =',F14.7,10X,'NGAUS  =',I6,
     &       /,5X,'CUBE   =',F14.7,10X,'MAXTAU =',I6,
     &       /,5X,'EQTIM  =',I6,10X,'NCYC    =',I6, 
     &       /,5X,'NLAYER =',F14.7,10X,'NLP =',I6,
     &       /,5X,'KH     =',F14.7,10X,'LATT =',I6,
     &       /,5X,'NPRINT =',I6
     &       /,5X,'KF     =',F14.7,10X,'R0  =',F14.7,
     &       /,5X,'LIMOL =',I6
     &       /,5X,'YZDIVX =',F14.7,10X,'CUBEX =',F14.7,
     &       /,5X,'CUBEY =',F14.7,10X,'CUBEZ =',F14.7,
     &       /,5X,'DXXDIV =',F14.7)
#endif

C
C****************************************************************
C
C***** INITIAL STARTUP
C
C
C***** NTYPE=1  INITIAL FROM FCC
C
         CALL FCC

C
C***** CHECK ALL FLUID PARTICLES ARE PART OF A MOLECULE
C        IF (MOD(NFLUID,LIMOL).NE.0.0D0) THEN
C           WRITE(6,*)'NPART DIFFERENT FROM INITIAL INPUT, NFLUID',
C      &               ' MUST BE DIVISIBLE BY LIMOL NO REMAINDER'
C           STOP
C        END IF
C
C***** UPDATE INPUT UNIT
C
C
#ifdef DEBUG
      WRITE(6,*) TR,DRW,DRF,DELTA,LATT,NPART,NLP
      WRITE(6,*) FE0,RCUT,NLAYER,KH,NPRINT
      WRITE(6,*) MIX, EPS1, EPS2, QVOL
      WRITE(6,*) KF,R0,LIMOL,YZDIVX
      WRITE(6,*) DXXDIV
      WRITE(6,*) NTYPE,NON,NGAUS,E0
      WRITE(6,*) NPLOT,MAXTAU,EQTIM,NCYC
      WRITE(6,*) 'TR,DRW,DRF,DELTA,LATT,NPART,NLP'
      WRITE(6,*) 'FE0,RCUT,NLAYER,KH,NPRINT'
      WRITE(6,*) 'KF,R0,LIMOL,YZDIVX'
      WRITE(6,*) 'MIX, EPS1, EPS2, QVOL'
      WRITE(6,*) 'DXXDIV'
      WRITE(6,*) 'NTYPE,NON,NGAUS,E0'
      WRITE(6,*) 'NPLOT,MAXTAU,EQTIM,NCYC'

C***** CHECK IMOL IKOL AND SWITCHES
        WRITE(6,*)'I,IMOL(I),IKOL(I),S1(I),S2(I),S3(I),S4(I)'
        DO I=1,NPART
         WRITE(6,*)I,IMOL(I),IKOL(I),S1(I),S2(I),S3(I),S4(I)
        ENDDO
#endif
C
C CALCULATE PROPERTIES AT TIME ZERO
C
         KE2=0.0D0
         DO 565 I=1, NPART
            KE2 = KE2+(PX(I)**2+PY(I)**2+PZ(I)**2)
 565     CONTINUE
c         write(6,*)'force 1'
         CALL FORCECELL 
         E00=UPOT+KE2/2.0D0
         FIELD=0.0D0

C***** UPDATE INITIAL TEMPERATURE
         TEMP = (KE2/2.0D0)/(DF/2.0D0)
       
 45   FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 75   FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 55   FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 85    FORMAT(/,2X,'EQUILIBRATION - INSTANTANEOUS VALUES',/,8X,
     &  'STEP    TEMP    UPOT/NP    TOTE/NP  ALPHA')
 65   FORMAT(6E14.6)
      END SUBROUTINE

