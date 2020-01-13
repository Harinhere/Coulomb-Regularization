
      implicit double precision(a-h,o-z)
      parameter(nmax=20000)
      dimension tt(nmax),xx(nmax,3),vv(nmax,3)
      dimension y(6),yout(6),fd(4)
      common/path/tt,xx,vv
      open(11,file='clb3.out')
      pi=2.*asin(1.d0)

c   initialize
      rho1=10.
	  z1=-500.
	  phi1=0.
	  vrho1=0.
	  vz1=0.11
c let's put in some lz
	  vphi1=0.
c zr is the residual charge
	  zr=1.
c f1 is the amplitude of the static field which is in -z direction
	  f1=0.
c f2 is the amplitude of the time rf field initially oriented in -z direction
	  f2=0.05338
c angular frequency of the rf field
	  omega=0.0147
c   coulomb initialize
	  y(1)=rho1
	  y(2)=vrho1
	  y(3)=z1
	  y(4)=vz1
	  y(5)=phi1
	  y(6)=vphi1
c  field initialize
      fd(1)=zr
      fd(2)=f1
	  fd(3)=f2
	  fd(4)=omega

	  t1=0.
	  t2=-1.
c if t2 is unknown at the begining, set it to a negative value and ajust stp.
c if t2 is known(>0), subroutine ignores rcl

	  rcoul=10.
	  zfar=2.*abs(z1)
c nt is the number of time steps
c iflag returns 1 for a successfull evaluation,iflag=0, otherwise.
      call coulomb(y,fd,t1,t2,rcoul,zfar,nt,iflag)
      if(iflag.eq.0) pause 'try a new trajectory'
      print*,nt,iflag
      do 3 i=1,nt
c      print*,tt(i),xx(i,1),xx(i,2),xx(i,3)
c check the conservation of energy
	  et=.5*(vv(i,1)**2+vv(i,2)**2+xx(i,1)*xx(i,1)*vv(i,3)**2)-
	 1 zr/sqrt(xx(i,1)**2+xx(i,2)**2)-f1*xx(i,2)-f2*xx(i,2)*
	 2 cos(omega*tt(i))
c convert to x,y,z coordinates
      xt=xx(i,1)*cos(xx(i,3))
	  yt=xx(i,1)*sin(xx(i,3))
	  zt=xx(i,2)
	  write(11,*),tt(i),xx(i,1),xx(i,2),et
 3    continue


      end
C****************************************************************************
C REGULARIZED COULOMB TRAJECTORIES IN 3D.
C VERSION: 04 (2018 DEC 20)
C AUTHORS: H.B AMBLAMPITIYA AND I.I FABRIKANT
C****************************************************************************
C THIS PROGRAM COMPUTES THE TRAJECTORIES OF AN ELECTRON IN COULOMB PLUS EXTERNAL
C ELECTRIC FIELDS. REGULARIZATION OF TRAJECTORIES 
C IS ACHIVED VIA A TIME TRANSFORMATION IN THE EXTENDED PHASE SPACE.
C
C FOR GIVEN COORDINATES AND VELOCITES AT T1, THIS SUBROUTINE PROPAGATES THE
C SOLUTION  UPTO T2. TRAJECTORY IS STORED IN A COMMON BLOCK
C
C DESCRIPTION OF THE METHOD, SEE: "REGULARIZED COULOMB TRAJECTORIES.PDF"
C
C INPUTS:
C       Y:  INITIAL POSITION AND VELOCITIES IN CYLYNDRICAL COORDINATES
C           Y(1)-RHO, Y(2)-RHO_DOT, Y(3)-Z, Y(4)-Z_DOT, Y(5)-PHI
C           Y(6)-PHI_DOT
C
C       FD: CONFIGURATION OF THE EXTERNAL ELECTRIC FIELDS + RESIDUE CHARGE
C           ALL FIELDS ARE DIRECTED ALONG THE NEGATIVE Z AXIS
C           
C           FD(1)-CHARGE OF THE ATOMIC RESIDUE
C           FD(2)-AMPLITUDE OF THE STATIC FIELD, FD(3)-AMPLITUDE OF THE
C           OSCILLATING FIELD, FD(4)-ANGULAR FREQUENCY (OMEGA) OF THE
C           OSCILLATING FIELD
C
C       T1: START-TIME
C       TF: END-TIME
C    RCOUL: IF TF IS NEGATIVE, "RCOUL" MATTERS. RCOUL IS THE RADIUS OF THE
C           COULOMB SPHERE
C     ZFAR: IF TF IS NEGATIVE, "ZFAR" MATTERS. ZFAR IS THE LIMITING DISTANCE
C           IN Z COORDINATE 
C OUTPUTS:
C
C      NT : NUMBER OF TIME STEPS TAKEN IN BETWEEN T2-T1
C
C      TT,XX,VV  : ARRAYS CONTAINING THE DETAILS OF THE TRAJECTORY
C                  TT- TIME STEPS, XX-COORDINATES, VV-VELOCITIES
C      USER MUST INCLUDE COMMON/PATH/TT,XX,VV IN THE MAIN PROGRAM
C      IFLAG : TYPES OF POSSIBLE OUTPUT TRAJECTORIES
C              IF IFLAG EQUALS
C
C             1. TRAJECTORY ENTERS AND LEAVES THE SPHERE (COULD BE MANY TIMES)
C                AND REACH THE END (IN Z COORDINATE) BEFORE THE TIME ITERATION
C                ENDS
C
C             2. TRAJECTORY DOES NOT ENTER OR LEAVE THE SPHERE BUT STILL
C                REACHES THE END (IN Z COORDINATE) BEFORE THE TIME ITERATION
C                ENDS
C
C             3. TRAJECTORY IS TRAPPED INSIDE THE SPHERE. IT MAY NOT LEAVE
C                THE SPHERE EVEN WHEN THE TIME ITERATION STOPS.
C
C             4. TRAJECTORY ENTERS AND LEAVES THE SPHERE (COULD BE MANY TIMES)
C                BUT DOES NOT REACH THE END (IN Z COORDINATE)BEFORE THE TIME ITERATION
C                ENDS. 
C
C             5. TRAJECTORY DOES NOT ENTER OR LEAVE THE SPHERE BUT DOES NOT
C                REACHE THE END (IN Z COORDINATE) BEFORE THE TIME ITERATION
C                ENDS
C 
C
C NOTE: IF IFLAG 3,4,5 HAPPENS, TRY A DIFFERENT INITIAL POSITION OR A DIFFERNT
C       FIELD CONFIGURATION OR INCREASE TMAX IN ACCORDANCE WITH NMAX
C****************************************************************************
C ADJUSTABLE PARAMETERS :
C      TSTEP : TSTEP IS THE TIME STEP IN TRANSFROMED TIME VARIABLE. 
C              TSTEP CAN BE CHANGED TO SMALLER VALUES IF ONE NEEDS MORE
C              SMOOTHER TRAJECTORIES NEAR THE CENTER. 
C 			   DEFAULT VALUE: TSTEP=0.05
C
C****************************************************************************
C****************************************************************************
C CODE BEGINS!

      SUBROUTINE COULOMB(Y,FD,T1,T2,RCOUL,ZFAR,NT,IFLAG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NMAX=20000)
      DIMENSION TT(NMAX),XX(NMAX,3),VV(NMAX,3),APPR(NMAX)
      DIMENSION Y(6),FD(4)
      DIMENSION Y1(7),DYDT(7),YOUT1(7)
      COMMON/SYS/F1,F2,OMEGA,ALZ,ZR
      COMMON/PATH/TT,XX,VV
      EXTERNAL DERIVS,BSSTEP
	  
	  HMIN=1.D-14
      HTRY=1.D-2
      EPSB=1.D-8
	  IFLAG=0
      
	  
      RHO1=Y(1)
	  VRHO1=Y(2)
	  Z1=Y(3)
      VZ1=Y(4)
	  PHI1=Y(5)
	  ALZ=RHO1*RHO1*Y(6)
        

      
C FIELD CONFIGURATION AND RESIDUAL CHARGE
      ZR=FD(1)
      F1=FD(2)
	  F2=FD(3)
	  OMEGA=FD(4)
C INITIAL CONJUGATED MOMENTUM FOR TIME
      PT1=-.5*(VRHO1*VRHO1+VZ1*VZ1)+ZR/SQRT(RHO1*RHO1+Z1*Z1)+
	 1 F1*Z1+F2*Z1*COS(OMEGA*T1)-.5*ALZ*ALZ/RHO1**2
C TRANSFORM TO THE SEMIPARABOLIC COORDINATES
      R1=SQRT(RHO1*RHO1+Z1*Z1)
      CHI1=SQRT(R1+Z1)
	  ETA1=SQRT(R1-Z1)
	  VCHI1=ETA1*VRHO1+CHI1*VZ1
	  VETA1=CHI1*VRHO1-ETA1*VZ1
C INITIALIZE RK4 STEPPER FOR REPS(REGULARIZED EXTENDED PHASE SPACE)
      Y1(1)=CHI1
	  Y1(2)=VCHI1
	  Y1(3)=ETA1
	  Y1(4)=VETA1
	  Y1(5)=T1
	  Y1(6)=PT1
	  Y1(7)=PHI1


	  T0=0.
      TSTEP=0.01
	  TMAX=200.
	  ITM=TMAX/TSTEP
	  IENT=0
	  ILEV=0
      DO 1 IT=1,ITM
	   TF=IT*TSTEP
       CALL ODEINT(Y1,7,T0,TF,EPSB,HTRY,HMIN,NOK,NBAD,DERIVS,
	 1 BSSTEP)

	   CHIT=Y1(1)
	   VCHIT=Y1(2)
	   ETAT=Y1(3)
	   VETAT=Y1(4)
	   TCY=Y1(5)
       IF(T2.GE.0. .AND. TCY.GE.T2)GOTO 99
	   
	   PT=Y1(6)
	   GT=CHIT*CHIT+ETAT*ETAT
C CONVERT TO THE ORIGINAL COORDINATES AND VELOCITIES
	   RHOT=CHIT*ETAT
	   VRHOT=(CHIT*VETAT+ETAT*VCHIT)/GT

	   ZT=.5*(CHIT*CHIT-ETAT*ETAT)
	   VZT=(CHIT*VCHIT-ETAT*VETAT)/GT

	   PHIT=Y1(7)
	   VPHIT=ALZ*(1./CHIT**2+1./ETAT**2)
C CHECK WHETHER THE TRAJECTORY ENTERS THE PREDEFINED COULOMB SPHERE
       APPR(IT)=SQRT(RHOT*RHOT+ZT*ZT)
       IF(IT.GT.1)THEN
	   RSP1=APPR(IT-1)
	   RSP2=APPR(IT)
C COUNTS THE ENTRANCES TO THE COULOMB SPHERE
        IF(RSP2.LE.RCOUL .AND. RSP1.GT.RCOUL)THEN
         IENT=IENT+1
C COUNT THE DEPARTURES FROM THE COULOMB SPHERE
        ELSE IF(RSP2.GT.RCOUL .AND. RSP1.LE.RCOUL)THEN
		 ILEV=ILEV+1
        ELSE
         CONTINUE
        ENDIF
       ENDIF
	   
       IF(IENT.GE.1 .AND. ILEV.GE.1 .AND. ABS(ZT).GE.ZFAR)THEN
	    IFLAG=1
        GOTO 99
       ENDIF
       IF(IENT.EQ.0 .AND. ILEV.EQ.0 .AND. ABS(ZT).GE.ZFAR)THEN
        IFLAG=2
        GOTO 99
       ENDIF
	   TORG=TCY
	   TT(IT)=TORG
	   XX(IT,1)=RHOT
	   XX(IT,2)=ZT
	   XX(IT,3)=PHIT
	   VV(IT,1)=VRHOT
	   VV(IT,2)=VZT
	   VV(IT,3)=VPHIT
C REINITIALIZE
       T0=TF
       
 1    CONTINUE
      IF(IENT.EQ.1 .AND. ILEV.EQ.0)THEN
       IFLAG=3
       WRITE(*,*)'TRAJECTOYR IS TRAPPED INSIDE THE COULOMB SPHERE'
      ELSE IF(IENT.GE.1 .AND. ILEV.GE.1 .AND. ABS(ZT).LT.ZFAR)THEN      
	   IFLAG=4
       WRITE(*,*)'COULOMB SPHERE MET BUT NOT REACHED TO ZFAR'
      ELSE IF(IENT.EQ.0 .AND. ILEV.EQ.0 .AND. ABS(ZT).LT.ZFAR)THEN      
	   IFLAG=5
       WRITE(*,*)'COULOMB SPHERE MISSED BUT NOT REACHED TO ZFAR'
      ELSE 
       CONTINUE
      ENDIF
 99   NT=IT-1
      RETURN
      END
C PROVIDES DERIVATIVES FOR RK4 STEPPER
      SUBROUTINE DERIVS(T,Y,DYDT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(7),DYDT(7)
      COMMON/SYS/F1,F2,OMEGA,ALZ,ZR
      TCY=Y(5)
      PT=Y(6)
	  DYDT(1)=Y(2)
	  DYDT(2)=2.*F2*COS(OMEGA*TCY)*Y(1)**3-2.*Y(1)*PT+
	 1 ALZ*ALZ/Y(1)**3+2.*F1*Y(1)**3
	  DYDT(3)=Y(4)
	  DYDT(4)=-2.*F2*COS(OMEGA*TCY)*Y(3)**3-2.*Y(3)*PT+
	 1 ALZ*ALZ/Y(3)**3-2.*F1*Y(3)**3
	  DYDT(5)=Y(1)*Y(1)+Y(3)*Y(3)
	  DYDT(6)=-.5*F2*OMEGA*SIN(OMEGA*TCY)*(Y(1)**4-Y(3)**4)
	  DYDT(7)=ALZ*(1./Y(1)**2+1./Y(3)**2)
      RETURN
      END

C REAL BULISCH-STOER 
      SUBROUTINE ODEINT(YSTART,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,
     1 RKQS)
C FOR BACKWARD INTEGRATION INITIAL DERIVATIVE IS IN POSITIVE DIRECTION!
      IMPLICIT NONE
      INTEGER NBAD,NOK,NVAR,KMAXX,MAXSTP,NMAX
      DOUBLE PRECISION EPS,H1,HMIN,X1,X2,TINY
      EXTERNAL DERIVS,RKQS
      PARAMETER (MAXSTP=100000,NMAX=50,KMAXX=200,TINY=1.D-30)
      INTEGER I,KMAX,KOUNT,NSTP
      DOUBLE PRECISION DXSAV,H,HDID,HNEXT,X,XSAV,XP(KMAXX),
     1 YSCAL(NMAX)
C      COMPLEX*16 DYDX(NMAX),Y(NMAX),YP(NMAX,KMAXX),YSTART(NVAR)  
      REAL*8 DYDX(NMAX),Y(NMAX),YP(NMAX,KMAXX),YSTART(NVAR)  
C      COMMON /PATH/ KMAX,KOUNT,DXSAV,XP,YP
C !USER STORAGE FOR INTERMEDIATE RESULTS. PRESET DXSAV AND KMAX.
C        PRINT 100,YSTART(2)
  100 FORMAT(E20.5)  
      X=X1
      H=SIGN(H1,X2-X1)
      NOK=0.D0
      NBAD=0.D0
      KOUNT=0.D0
      DO I=1,NVAR
        Y(I)=YSTART(I)
      ENDDO
      IF (KMAX.GT.0) XSAV=X-2.*DXSAV
      DO NSTP=1,MAXSTP
        CALL DERIVS(X,Y,DYDX)
        DO I=1,NVAR
C !SCALING USED TO MONITOR ACCURACY. THIS GENERAL-PURPOSE CHOICE CAN 
C BE MODIFIED IF NEED BE.
            YSCAL(I)=ABS(Y(I))+ABS(H*DYDX(I))+TINY
        ENDDO 
        IF(KMAX.GT.0)THEN
            IF(ABS(X-XSAV).GT.ABS(DXSAV)) THEN
                IF(KOUNT.LT.KMAX-1)THEN
                    KOUNT=KOUNT+1
                    XP(KOUNT)=X
                    DO I=1,NVAR
                        YP(I,KOUNT)=Y(I)
                    ENDDO 
                    XSAV=X
                ENDIF
            ENDIF
        ENDIF
        IF((X+H-X2)*(X+H-X1).GT.0.) H=X2-X
            CALL RKQS(Y,DYDX,NVAR,X,H,EPS,YSCAL,HDID,HNEXT,DERIVS)
            IF(HDID.EQ.H)THEN
                NOK=NOK+1
            ELSE
                NBAD=NBAD+1
            ENDIF
            IF((X-X2)*(X2-X1).GE.0.)THEN
            DO I=1,NVAR
                YSTART(I)=Y(I)
            ENDDO 
            IF(KMAX.NE.0)THEN
                KOUNT=KOUNT+1 
                XP(KOUNT)=X
                DO I=1,NVAR
                    YP(I,KOUNT)=Y(I)
                ENDDO 
            ENDIF
            RETURN !NORMAL EXIT.
        ENDIF
        IF(ABS(HNEXT).LT.HMIN) PAUSE 
     1 'STEPSIZE SMALLER THAN MINIMUM IN ODEINT'
            H=HNEXT
      ENDDO 
      PAUSE 'TOO MANY STEPS IN ODEINT'
      RETURN
      END
C END SUBROUTINE ODEINT


      SUBROUTINE BSSTEP(Y,DYDX,NV,X,HTRY,EPS,YSCAL,HDID,HNEXT,DERIVS)
      IMPLICIT NONE
      INTEGER NV,NMAX,KMAXX,IMAX
      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,YSCAL(NV),
     1 SAFE1,SAFE2,REDMAX,REDMIN,TINY,SCALMX
C      COMPLEX*16 DYDX(NV),Y(NV)
      REAL*8 DYDX(NV),Y(NV)
C	  DIMENSION DYDX(NV),Y(NV)
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7,
     1 REDMAX=1.E-5,REDMIN=.7,TINY=1.E-30,SCALMX=.1)
C !C USES DERIVS,MMID,PZEXTR
C !BULIRSCH-STOER STEP WITH MONITORING OF LOCAL TRUNCATION ERROR TO 
C ENSURE ACCURACY AND ADJUST
      INTEGER I,IQ,K,KK,KM,KMAX,KOPT,NSEQ(IMAX)
      DOUBLE PRECISION EPS1,EPSOLD,ERRMAX,FACT,H,RED,SCALE,WORK,
     1 WRKMIN,XEST,XNEW
      DOUBLE PRECISION A(IMAX),ALF(KMAXX,KMAXX),ERR(KMAXX)
C      COMPLEX*16 YSAV(NMAX),YSEQ(NMAX),YERR(NMAX)
      REAL*8 YSAV(NMAX),YSEQ(NMAX),YERR(NMAX)
      LOGICAL FIRST,REDUCT
      SAVE A,ALF,EPSOLD,FIRST,KMAX,KOPT,NSEQ,XNEW
      EXTERNAL DERIVS
      DATA FIRST/.TRUE./,EPSOLD/-1./
      DATA NSEQ /2,4,6,8,10,12,14,16,18/
      IF(EPS.NE.EPSOLD)THEN 
        HNEXT=-1.D290       
        XNEW=-1.D290
        EPS1=SAFE1*EPS
        A(1)=NSEQ(1)+1    
        DO K=1,KMAXX
             A(K+1)=A(K)+NSEQ(K+1)
        ENDDO
        DO IQ=2,KMAXX   
            DO K=1,IQ-1
         ALF(K,IQ)=EPS1**((A(K+1)-A(IQ+1))/((A(IQ+1)-A(1)+1.)*(2*K+1)))
            ENDDO 
        ENDDO 
        EPSOLD=EPS
        DO KOPT=2,KMAXX-1 
           IF(A(KOPT+1).GT.A(KOPT)*ALF(KOPT-1,KOPT))GOTO 1
        ENDDO
1       KMAX=KOPT
      ENDIF
      H=HTRY
      DO I=1,NV !SAVE THE STARTING VALUES.
        YSAV(I)=Y(I)
      ENDDO 
      IF(H.NE.HNEXT.OR.X.NE.XNEW)THEN
        FIRST=.TRUE.
        KOPT=KMAX
      ENDIF
      REDUCT=.FALSE.
2     DO K=1,KMAX
       XNEW=X+H
       IF(XNEW.EQ.X) THEN
	    PAUSE 'STEP SIZE UNDERFLOW IN BSSTEP'
       END IF
       CALL MMID(YSAV,DYDX,NV,X,H,NSEQ(K),YSEQ,DERIVS)
       XEST=(H/NSEQ(K))**2 
       CALL PZEXTR(K,XEST,YSEQ,Y,YERR,NV) 
       IF(K.NE.1)THEN 
          ERRMAX=TINY
          DO I=1,NV
             ERRMAX=MAX(ERRMAX,ABS(YERR(I)/YSCAL(I)))
          ENDDO 
          ERRMAX=ERRMAX/EPS 
          KM=K-1
          ERR(KM)=(ERRMAX/SAFE1)**(1./(2*KM+1))
       ENDIF
       IF(K.NE.1.AND.(K.GE.KOPT-1.OR.FIRST))THEN 
          IF(ERRMAX.LT.1.)GOTO 4 !CONVERGED.
          IF(K.EQ.KMAX.OR.K.EQ.KOPT+1)THEN 
             RED=SAFE2/ERR(KM)
             GOTO 3
          ELSE IF(K.EQ.KOPT)THEN
             IF(ALF(KOPT-1,KOPT).LT.ERR(KM))THEN
                RED=1./ERR(KM)
                GOTO 3
             ENDIF
          ELSE IF(KOPT.EQ.KMAX)THEN
             IF(ALF(KM,KMAX-1).LT.ERR(KM))THEN
                RED=ALF(KM,KMAX-1)*SAFE2/ERR(KM)
                GOTO 3
             ENDIF
          ELSE IF(ALF(KM,KOPT).LT.ERR(KM))THEN
             RED=ALF(KM,KOPT-1)/ERR(KM)
             GOTO 3
          ENDIF
       ENDIF
      ENDDO 
3     RED=MIN(RED,REDMIN)
      RED=MAX(RED,REDMAX)
      H=H*RED
      REDUCT=.TRUE.
      GOTO 2 
4     X=XNEW
      HDID=H
      FIRST=.FALSE.
      WRKMIN=1.D300
      DO KK=1,KM
       FACT=MAX(ERR(KK),SCALMX)
       WORK=FACT*A(KK+1)
       IF(WORK.LT.WRKMIN)THEN
          SCALE=FACT
          WRKMIN=WORK
          KOPT=KK+1
       ENDIF
      ENDDO
      HNEXT=H/SCALE
      IF(KOPT.GE.K.AND.KOPT.NE.KMAX.AND..NOT.REDUCT)THEN 
       FACT=MAX(SCALE/ALF(KOPT-1,KOPT),SCALMX)
       IF(A(KOPT+1)*FACT.LE.WRKMIN)THEN
          HNEXT=H/FACT
          KOPT=KOPT+1 
       ENDIF
      ENDIF
      RETURN
      END
C END SUBROUTINE BSSTEP


      SUBROUTINE PZEXTR(IEST,XEST,YEST,YZ,DY,NV)
      IMPLICIT NONE
      INTEGER IEST,NV,IMAX,NMAX
      DOUBLE PRECISION XEST
C      COMPLEX*16 YZ(NV),YEST(NV),DY(NV)
      REAL*8 YZ(NV),YEST(NV),DY(NV)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER J,K1
      DOUBLE PRECISION X(IMAX)
C      COMPLEX*16 D(NMAX),QCOL(NMAX,IMAX),DELTA,F1,F2,Q
      REAL*8 D(NMAX),QCOL(NMAX,IMAX),DELTA,F1,F2,Q
      SAVE QCOL,X
      X(IEST)=XEST
      DO J=1,NV
       DY(J)=YEST(J)
       YZ(J)=YEST(J)
      ENDDO 
      IF(IEST.EQ.1) THEN 
       DO J=1,NV
         QCOL(J,1)=YEST(J)
       ENDDO 
      ELSE
       DO J=1,NV
         D(J)=YEST(J)
       ENDDO 
       DO K1=1,IEST-1
         DELTA=1./(X(IEST-K1)-XEST)
         F1=XEST*DELTA
         F2=X(IEST-K1)*DELTA
         DO J=1,NV 
            Q=QCOL(J,K1)
            QCOL(J,K1)=DY(J)
            DELTA=D(J)-Q
            DY(J)=F1*DELTA
            D(J)=F2*DELTA
            YZ(J)=YZ(J)+DY(J)
         ENDDO 
       ENDDO 
       DO J=1,NV
         QCOL(J,IEST)=DY(J)
       ENDDO 
      ENDIF
      RETURN
      END
C END SUBROUTINE PZEXTR
C
      SUBROUTINE MMID(Y,DYDX,NVAR,XS,HTOT,NSTEP,YOUT,DERIVS)
      IMPLICIT NONE
      INTEGER NSTEP,NVAR,NMAX
      DOUBLE PRECISION HTOT,XS
C      COMPLEX*16 DYDX(1:NVAR),Y(1:NVAR),YOUT(1:NVAR)
      REAL*8 DYDX(1:NVAR),Y(1:NVAR),YOUT(1:NVAR)
      EXTERNAL DERIVS
      PARAMETER (NMAX=50)
!MODIFIED MIDPOINT STEP. DEPENDENT VARIABLE VECTOR Y(1:NVAR) AND ITS DERIVATIVE VECTOR
!DYDX(1:NVAR) ARE INPUT AT XS. ALSO INPUT IS HTOT, THE TOTAL STEP TO BE MADE, AND NSTEP,
!THE NUMBER OF SUBSTEPS TO BE USED. THE OUTPUT IS RETURNED AS YOUT(1:NVAR), WHICH NEED
!NOT BE A DISTINCT ARRAY FROM Y; IF IT IS DISTINCT, HOWEVER, THEN Y AND DYDX ARE RETURNED
!UNDAMAGED.
      INTEGER I,N
      DOUBLE PRECISION H,H2,X
C      COMPLEX*16 YM(NMAX),YN(NMAX),SWAP
      REAL*8 YM(NMAX),YN(NMAX),SWAP
      H=HTOT/NSTEP !STEPSIZE THIS TRIP.
      DO I=1,NVAR
      YM(I)=Y(I)
      YN(I)=Y(I)+H*DYDX(I) !FIRST STEP.
      ENDDO 
      X=XS+H
      CALL DERIVS(X,YN,YOUT) !WILL USE YOUT FOR TEMPORARY STORAGE OF DERIVATIVES.
      H2=2.*H
      DO N=2,NSTEP !GENERAL STEP.
      DO I=1,NVAR
         SWAP=YM(I)+H2*YOUT(I)
         YM(I)=YN(I)
         YN(I)=SWAP
      ENDDO 
      X=X+H
      CALL DERIVS(X,YN,YOUT)
      ENDDO 
      DO I=1,NVAR !LAST STEP.
      YOUT(I)=0.5D0*(YM(I)+YN(I)+H*YOUT(I))
      ENDDO 
      RETURN
      END



