      IMPLICIT NONE
C Parameters controlling floating point precision
      integer, parameter :: sp = selected_real_kind(6, 37)
      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: qp = selected_real_kind(33, 4931)
      integer, parameter :: prec = dp
      INTEGER NP,NK,NPROP,NTCF,NL,NBINS,BINS,DIM
      PARAMETER (NP=3000,NK=10000,NTCF=22,NPROP=11,NBINS=150)
      PARAMETER (DIM=3)
      INTEGER NTYPE,IKOL,IMOL,LIMOL,NPART2,NLP,LATT
      REAL(kind=prec) DPXY,DPXZ,DPYZ,DAL,DPB,DBETA,DF
      REAL(kind=prec) X,Y,Z,X0,Y0,Z0,XEQ,YEQ,ZEQ,XR,YR,ZR
      REAL(kind=prec) PX,PY,PZ,FX,FY,FZ,CONSER
      REAL(kind=prec) FXF,FYF,FZF,FXW,FYW,FZW
      REAL(kind=prec) TR,E0,UPOT,UWALL,PT,PTF,DRW,DRF,RCUT,EPS12
      REAL(kind=prec) EPS1,EPS2,MIX
      REAL(kind=prec) S1,S2,S3,S4,SL,S0L,KH,NLAYER,KF,R0
      REAL(kind=prec) DXX,DXXDIV
      REAL(kind=prec) YZDIVX,CUBE,CUBE2,CUBEX,CUBEY,CUBEZ
      REAL(kind=prec) CUBEX2,CUBEY2,CUBEZ2
      REAL(kind=prec) DELTA,FIELD,FE0,RMAX,SHIFT
      REAL(kind=prec) SHEAR,DXD,DXDSAV
      INTEGER ISEED,ISEED1,NGAUS
      INTEGER KOUNT,NMEAN,NTAUAV,NPLOT,KB,JFLAG,NPRINT
      INTEGER NPART,MAXTAU,NON,PTAU,EQTIM,NCYC,TAUNBINS
      REAL(kind=prec) ALPH,ALPHB,VOL,VOLF,MEAN,TCF,INTCF
      REAL(kind=prec) MSDX,MSDY,MSDZ, MSD
      REAL(kind=prec) TAUAV,TAUAVL,T0,T0L,T0V
      REAL(kind=prec) M0,M1,M2,M3,M4,M5,KE2,E00,NHALPH
      REAL(kind=prec) XSAV,YSAV,ZSAV,PXSAV,PYSAV,PZSAV,UPOTSAV
      REAL(kind=prec) UWALLSAV,KE2SAV,NHSAV,DENSAV1,DENSAV2 
      REAL(kind=prec) FXSAV,FYSAV,FZSAV,ALSAV,ALBSAV,TCF1,INTCF1
      REAL(kind=prec) TCFL,TCFLL,INTCFL,DISSL, INTCFLL 
      REAL(kind=prec) XEQSAV,YEQSAV,ZEQSAV
      REAL(kind=prec) FXFSAV,FYFSAV,FZFSAV,FXWSAV,FYWSAV,FZWSAV
      REAL(kind=prec) AVER,AVERL,MAXV,MINV,DEL,ENAV,DENS1,DENS2,QVOL
      REAL(kind=prec) TAUDENS1,TAUDENS2
      REAL(kind=prec) TEMP,TOTE
      INTEGER RIJ_HIST
      REAL(kind=prec) RBIN_INV
      LOGICAL DO_PRESSURE

      INTEGER NENAV
      COMMON/FLUX/  UPOT,UWALL,PT(3,3),PTF(3,3)
      COMMON/NOPART/ NPART,IKOL(NP),IMOL(NP),LIMOL,NPART2
      COMMON/IPARM/ MAXTAU,NON,NPLOT,ISEED,NGAUS,KB,PTAU,BINS,
     &              EQTIM,NCYC,NTYPE,JFLAG,NPRINT,LATT,TAUNBINS
      COMMON/COORD/ X(NP),Y(NP),Z(NP),XEQ(NP),YEQ(NP),ZEQ(NP),
     &              QVOL,DF,
     &              X0(NP),Y0(NP),Z0(NP),XR(NP),YR(NP),ZR(NP)
      COMMON/MOMENT/ PX(NP),PY(NP),PZ(NP)
      COMMON/FORSE/ FX(NP),FY(NP),FZ(NP),FXF(NP),FYF(NP),FZF(NP)
      COMMON/FORSE/ FXW(NP),FYW(NP),FZW(NP)
      COMMON/SWITCH/ S1(NP),S2(NP),S3(NP),S4(NP),SL(NP),S0L(NP)
      COMMON/XSAV/ XSAV(NP),YSAV(NP),ZSAV(NP),PXSAV(NP),PYSAV(NP),
     &             PZSAV(NP),XEQSAV(NP),YEQSAV(NP),ZEQSAV(NP),
     &             FXSAV(NP),FYSAV(NP),FZSAV(NP),
     &             FXFSAV(NP),FYFSAV(NP),FZFSAV(NP),FXWSAV(NP),
     &             FYWSAV(NP),FZWSAV(NP),NHSAV,
     &             UPOTSAV,UWALLSAV,KE2SAV,ALSAV,ALBSAV 
      COMMON/INENER/  TR,E0
      COMMON/PARM/  FIELD,FE0,YZDIVX,CUBE,CUBE2,CUBEX,CUBEY,CUBEZ,
     &              CUBEX2,CUBEY2,CUBEZ2,RMAX,SHIFT,CONSER,
     &              KH,KF,R0,DRW,DRF
      COMMON/PARM/  RCUT,EPS12,EPS1,EPS2,MIX
      COMMON/PARM/  DXX,DXXDIV
      COMMON/SIMUL/ DELTA
      COMMON/IRANDI/ ISEED1
      COMMON/NUMAV/ KOUNT(0:NK),NMEAN,NTAUAV
      COMMON/COUNT/ T0(DIM),TCF(0:NK,NTCF),T0V(NP,DIM),
     &              MSDX,MSDY,MSDZ, MSD, 
     &              TCF1(0:NK,NTCF),INTCF1(0:NK,NTCF),T0L(DIM),
     &              TCFL(0:NK,0:NBINS,NTCF),INTCFL(0:NK,0:NBINS,NTCF), 
     &              INTCF(0:NK,NTCF),DISSL(0:NBINS), 
     &              TCFLL(0:NK,0:NBINS,NTCF),INTCFLL(0:NK,0:NBINS,NTCF)  
      COMMON/AVERG/ MEAN(NPROP),AVER(0:NK,NTCF),DENS1(0:NBINS),
     &              DENS2(0:NBINS), AVERL(0:NK,0:NBINS,NTCF),
     &              TAUDENS1(0:NBINS), TAUDENS2(0:NBINS),
     &              DENSAV1(0:NK,0:NBINS),DENSAV2(0:NK,0:NBINS),
     &              TEMP,TOTE,
     &              RIJ_HIST(NBINS), RBIN_INV
      COMMON/TAUAV/ TAUAV(NTCF), TAUAVL(0:NBINS,NTCF) 
      COMMON/PROPS/ ALPH,ALPHB,VOL,VOLF,KE2,E00,NHALPH
      COMMON/MOM/   M0(NTCF),M1(NTCF),M2(NTCF),M3(NTCF),M4(NTCF),
     &              M5(NTCF)
      COMMON/IOFLAGS/ DO_PRESSURE
      COMMON/SLLOD/ SHEAR,DXD,DXDSAV
C

