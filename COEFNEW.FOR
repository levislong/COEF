C  PROGRAM COEF                   Version am 30.10.1989 von ZOU ZAOJIAN
C                               geaendert am 21.10.1993 von ZOU ZAOJIAN
C  Dieses Hauptprogramm dient zur Berechnung der dimensionslosen
C  hydrodynamischen Koeffizienten des manoevrierenden Schiffes 
C  auf tiefem und/oder flachem Wasser bei endlicher Froudezahl. 
C
* -------------------------------------------------------------------- *
*  Aufgerufene Unterprogramme:                                         *
*  TUCK: Berechnet die Absenkung und die Vertrimmung infolge der       *
*        Fahrgeschwindigkeit nach dem Verfahren von E.O.TUCK           *
*  CSS:  Berechnet die Versperrungskoeffizienten der symmetrischen     *
*        Schiffsspanten                                                *
*  QGVS: Berechnet die Queranstroemungsgeschwindigkeit im Nahfeld      *
*        infolge der reinen Quer- bzw. Drehbewegung des Schiffes       *
* -------------------------------------------------------------------- *
C
      PARAMETER (NSP=21,NV=10)
      DIMENSION XX(NSP),T(NSP),B(NSP),A(NSP)
      DIMENSION AA(NSP),TT(NSP),HMY(NSP),FX(NSP),AX(NSP),BX(NSP)
      DIMENSION CC(NSP),BB(NSP,NSP+1),B0(NSP-1,NSP),
     &          VSTAR(NSP),ALFA(NSP)
      DIMENSION FN(NV),VER(NV)
      CHARACTER*32 SHIPNAME 
      REAL KG
C-----------------------------------------------------------------------
C  !! Reading the data of the ship:                                        
C     SL     - Length
C     DX     - Distance between the frames 
C     N      - Number of frames (N < 22) 
C     XG     - X-coordinate of the center of gravity
C     KG     - Distance between the center of gravity and the keel.
C              If KG is unknown, KG = 0 should be used.
C     XX     - X-coordinate of the frames
C     T,B,A  - Draft, width and area of the frames
C     Dieses Programm gilt fuer den Fall 
C     -SL/2 < X(1) < X(2) ... < X(N) < SL/2  und
C     X(N)-X(N-1)=X(N-1)-X(N-2)=...=X(2)-X(1)=DX, 
C     worin X(0)=-SL/2, X(N+1)=SL/2.
C-----------------------------------------------------------------------
      OPEN (unit = 5,file = 'COEFEIN.DAT',STATUS='OLD',ERR=1)
      READ (5,*,ERR=2) SHIPNAME
      READ (5,*,ERR=2) SL,DX,N,XG,KG
      READ (5,*,ERR=2) (XX(I),I=1,N)
      READ (5,*,ERR=2) (T(I),I=1,N)
      READ (5,*,ERR=2) (B(I),I=1,N)
      READ (5,*,ERR=2) (A(I),I=1,N)
      GOTO 3
1     STOP '=> Eingabedatei COEFEIN.DAT nicht gefunden !'
2     STOP '=> Fehlerhafte Eingabedatei COEFEIN.DAT !'
3     CONTINUE
C
C  => Calculation of the water-line area Aw, the X-coordinate of the WL center of gravity XW 
C     , The displacement VD, the length-of-metacentric height GML:
      SUM=0.5*(XX(1)+SL/2.)*B(1)+0.5*(SL/2.-XX(N))*B(N)
      AW=SUM+SIMP(B,DX,N)
      SUM=0.5*(XX(1)+SL/2.)*A(1)+0.5*(SL/2.-XX(N))*A(N)
      VD=SUM+SIMP(A,DX,N)
      DO 4 I=1,N
      FX(I)=XX(I)*B(I)
4     CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      XW=SUM/AW
C  !! GML=BML+KB-KG, naeherungsweise: GML=BML+0.5*(T(N)+T(1))/2.-KG
      DO 5 I=1,N
      FX(I)=(XX(I)-XW)*(XX(I)-XW)*B(I)
5     CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      IF (KG.EQ.0.) THEN
       GML=SUM/VD
      ELSE 
       GML=SUM/VD+0.5*(T(N)+T(1))/2.-KG
      END IF
C
      OPEN (6,file='COEFRES.DAT')
      WRITE (6,'(A/)') SHIPNAME
      WRITE (6,*) ' Aw=',AW,'   Xw=',XW,'   Verd.=',VD,'   GML=',GML
C  => Berechnung der Ableitung DB/DX und der Ableitung DA/DX
      DO 6 I=1,N-1
      BX(I)=(B(I+1)-B(I))/DX
      AX(I)=(A(I+1)-A(I))/DX
6     CONTINUE
      BX(N)=BX(N-1)
      AX(N)=AX(N-1)          
C
C  !! Eingaben
C     FN:  Froude Number, FN=U/SQRT(g*SL)      
C     VER: Water depth-draft-Relations, VER=H/T
      READ (5,*,ERR=2) NFN,NVER
      IF (NFN > NV) STOP '=> Too many FN input !'
      IF (NVER > NV) STOP '=> To many VER input !'
      READ (5,*,ERR=2) (FN(I),I=1,NFN)
      READ (5,*,ERR=2) (VER(I),I=1,NVER)
* -------------------------------------------------------------------- *
C
      DO 100 IFN=1,NFN
      DO 100 JVER=1,NVER
      H=VER(JVER)*(T(N)+T(1))/2.
      FNH=FN(IFN)*SQRT(SL/H)
      WRITE(6,*)
      WRITE(6,*) 
      WRITE (6,*) 'Founde Number Fn=', FN(IFN)
      WRITE (6,*) 'Depth Founde Number Fnh=',FNH
      WRITE (6,*) 'Depth-Draft ratio VER=',VER(JVER)
!10    FORMAT (//1X,'Fn=',F6.4,',    Fnh=',F6.4,',    H/T=',F4.2)
C										
C  => Berechnung der Absenkung DD und der Vertrimmung SITA(grad) des 
C     Schiffes infolge der Geradeausfahrt nach dem Verfahren von Tuck:
      CALL TUCK(XX,B,BX,AX,DX,N,AW,XW,VD,GML,FNH,DD,SITA)
      WRITE (6,*) 'Absenkung DD=',DD
      WRITE (6,*) 'Vertrimmung SITA(grad)=',SITA
      WRITE (6,*)
C  => Berechnung der Tiefgaenge TT(I) und der Spantflaechen AA(I), 
C     worin TT=T+DT und AA=A+DA (DT, DA: Veraenderungen der Tiefgaenge 
C     und der Spantflaechen infolge der Absenkung und der Vertrimmung):
      DO 15 I=1,N
      TT(I)=T(I)+DD-(XX(I)-XW)*SITA*3.1415927/180.
      AA(I)=A(I)+B(I)*(DD-(XX(I)-XW)*SITA*3.1415927/180.)
15    CONTINUE
C
C  => Berechnung der Versperrungskoeffizienten CC  
      DO 20 I=1,N
      CALL CSS(B(I),TT(I),AA(I),H,CC(I))
20    CONTINUE
*      WRITE (6,*) 'Versperrungskoeffizienten CC:'
*      WRITE (6,*) CC
C
C  => Berechnung ALFA und VSTAR infolge der reinen Querbewegung:
      CALL QGVS(1.,0.,FNH,CC,BB,B0,XX,DX,SL,N-1,N,N+1)
      DO 30 I=1,N
      VSTAR(I)=BB(I,N+1)
      ALFA(I)=BB(I,N+1)*CC(I)
30    CONTINUE
*      WRITE (6,*) 'ALFA infolge der reinen Querbewegung:'
*      WRITE (6,*) ALFA
*      WRITE (6,*) 'VSTAR infolge der reinen Querbewegung:'
*      WRITE (6,*) VSTAR
C  => Berechnung der Verteilung der hydrodynamischen Masse HMY(X) 
C     infolge der reinen Querbewegung und der Stelle XX(IMAX), wo 
C     der hydrodynamische Impuls infolge der reinen Querbewegung 
C     am groessten ist.
      DO 31 I=1,N
      HMY(I)=-AA(I)+2.*H*ALFA(I)
31    CONTINUE
*      WRITE (6,*) 'My(X) infolge der reinen Querbewegung'
*      WRITE (6,*) HMY
      FMAX=ABS(HMY(1))
      IMAX=1
      DO 32 I=2,N
      IF (ABS(HMY(I)).GT.FMAX) THEN
      FMAX=ABS(HMY(I))
      IMAX=I
      END IF
32    CONTINUE
      WRITE (6,33) IMAX
33    FORMAT(1X,'IMAXV=',I2)
C  => Berechnung der dimensionslosen hydrodynamischen Masse 
C     HM22' (-Ydv') und der dimensionslosen hydrodynamischen 
C     Koeffizienten Yv',Nv',Ndv':
      DO 40 I=1,N
      FX(I)=-(AA(I)-A(I))+2.*H*ALFA(I)
40    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      HSM2=2.*SUM/(SL**3)
      SUM=0.5*(XX(1)+SL/2.)*HMY(1)+0.5*(SL/2.-XX(N))*HMY(N)
      SUM=SUM+SIMP(HMY,DX,N)
      HM22=2.*SUM/(SL**3)
      SM=HSM2-HM22
      WRITE (6,*) 'M22''+M''=',HSM2,'   M22''=',HM22,'   M''=',SM
      YV=-2.*HMY(IMAX)/(SL*SL)
      DO 41 I=1,N
      FX(I)=XX(I)*HMY(I)
41    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      HNDV=-2.*SUM/(SL**4)
      NIMAX=(N-IMAX)/2
      IF (2*NIMAX.EQ.N-IMAX) THEN
C     N-IMAX ist eine gerade Zahl.
      FX(1)=HMY(IMAX)
      DO 42 I=1,N-IMAX
      FX(I+1)=HMY(IMAX+I)
42    CONTINUE
      SUM=0.5*(SL/2.-XX(N))*HMY(N)
      SUM=SUM+SIMP(FX,DX,N-IMAX+1)
      ELSE
      DO 43 I=1,N-IMAX
      FX(I)=HMY(IMAX+I)
43    CONTINUE
      SUM=0.5*DX*(HMY(IMAX)+HMY(IMAX+1))
      SUM=SUM+0.5*(SL/2.-XX(N))*HMY(N)
      SUM=SUM+SIMP(FX,DX,N-IMAX)
      END IF
      HNV=2.*(-XX(IMAX)*HMY(IMAX)-SUM)/(SL**3)
      WRITE (6,*) 'Yv''=',YV,'   Ndv''=',HNDV,'   Nv''=',HNV 
C
C  => Berechnung ALFA und VSTAR infolge der reinen Drehbewegung
      CALL QGVS(0.,1.,FNH,CC,BB,B0,XX,DX,SL,N-1,N,N+1)
      DO 50 I=1,N
      VSTAR(I)=BB(I,N+1)
      ALFA(I)=BB(I,N+1)*CC(I)
50    CONTINUE
*      WRITE (6,*) 'ALFA infolge der reinen Drehbewegung:'
*      WRITE (6,*) ALFA
*      WRITE (6,*) 'VSTAR infolge der reinen Drehbewegung:'
*      WRITE (6,*) VSTAR
C  => Berechnung der Verteilung der hydrodynamischen Masse HMY(X) 
C     infolge der reinen Drehbewegung und der Stelle XX(IMAX), wo 
C     der hydrodynamische Impuls infolge der reinen Drehbewegung 
C     am groessten ist.
      DO 51 I=1,N
      IF (XX(I).EQ.0.) THEN
      HMY(I)=0.
      ELSE  
      HMY(I)=-AA(I)+2.*H*ABS(ALFA(I)/XX(I))
      IF (HMY(I).LT.0.) HMY(I)=0.
      END IF
51    CONTINUE
      SUM=0.
      NSUM=N
      DO 52 I=1,N
      SUM=SUM+HMY(I)
      IF (HMY(I).EQ.0.) NSUM=NSUM-1
52    CONTINUE
      HMY0=SUM/NSUM
      DO 53 I=1,N
      IF (HMY(I).EQ.0..OR.HMY(I).GT.SUM-HMY(I)) HMY(I)=HMY0
53    CONTINUE
      DO 54 I=3,N-2
      IF (HMY(I).GT.HMY(I-1)+HMY(I+1)) HMY(I)=(HMY(I-1)+HMY(I+1))/2.
54    CONTINUE
*      DO 55 I=1,N
*      FX(I)=XX(I)*HMY(I)
*55    CONTINUE  
*      WRITE (6,*) 'My(X) infolge der reinen Drehbewegung'
*      WRITE (6,*) HMY
*      WRITE (6,*) 'IMPULSy(X) infolge der reinen Drehbewegung'
*      WRITE (6,*) FX
      FMAX=ABS(HMY(1))
      IMAX=1
      DO 56 I=2,N
      IF (ABS(HMY(I)).GT.FMAX) THEN
      FMAX=ABS(HMY(I))
      IMAX=I
      END IF
56    CONTINUE
      WRITE (6,57) IMAX
57    FORMAT(1X,'IMAXR=',I2)  
C  => Berechnung des dimensionslosen hydrodynamischen Massentraegheits-
C     momentes HM66'(-Ndr') und der dimensionslosen hydrodynamischen 
C     Koeffizienten Yr', Ydr', Nr':
      DO 60 I=1,N
      FX(I)=XX(I)*XX(I)*(A(I)+HMY(I))
60    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      HSM6=2.*SUM/(SL**5)
      DO 61 I=1,N
      FX(I)=XX(I)*XX(I)*HMY(I)
61    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      HM66=2.*SUM/(SL**5)
      SIZ=HSM6-HM66
      WRITE (6,*) 'M66''+Izz''=',HSM6,'   M66''=',HM66,'   Izz''=',SIZ
      DO 62 I=1,N
      FX(I)=XX(I)*HMY(I)
62    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+0.5*(SL/2.-XX(N))*FX(N)
      SUM=SUM+SIMP(FX,DX,N)
      YDR=-2.*SUM/(SL**4)
      I1=(IMAX-1)/2
      IF (2*I1.EQ.IMAX-1) THEN
C     IMAX-1 ist eine gerade Zahl.
      DO 63 I=1,IMAX
      FX(I)=HMY(I)
63    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)+SIMP(FX,DX,IMAX)
      DO 64 I=1,IMAX 
      FX(I)=XX(I)*HMY(I)
64    CONTINUE
      SUM1=0.5*(XX(1)+SL/2.)*FX(1)+SIMP(FX,DX,IMAX)
      ELSE
      DO 65 I=1,IMAX-1
      FX(I)=HMY(I)
65    CONTINUE
      SUM=0.5*(XX(1)+SL/2.)*FX(1)
      SUM=SUM+0.5*DX*(FX(IMAX-1)+HMY(IMAX))  
      SUM=SUM+SIMP(FX,DX,IMAX-1)
      DO 66 I=1,IMAX-1
      FX(I)=XX(I)*HMY(I)
66    CONTINUE
      SUM1=0.5*(XX(1)+SL/2.)*FX(1)
      SUM1=SUM1+0.5*DX*(FX(IMAX-1)+XX(IMAX)*HMY(IMAX))
      SUM1=SUM1+SIMP(FX,DX,IMAX-1)
      END IF
      YR=2.*(-XX(IMAX)*HMY(IMAX)+SUM)/(SL**3)
      NIMAX=(N-IMAX)/2
      IF (2*NIMAX.EQ.N-IMAX) THEN
C     N-IMAX ist eine gerade Zahl.
      FX(1)=XX(IMAX)*HMY(IMAX)
      DO 67 I=1,N-IMAX
      FX(I+1)=XX(IMAX+I)*HMY(IMAX+I)
67    CONTINUE
      SUM2=0.5*(SL/2.-XX(N))*XX(N)*HMY(N)
      SUM2=SUM2+SIMP(FX,DX,N-IMAX+1)
      ELSE
      DO 68 I=1,N-IMAX
      FX(I)=XX(IMAX+I)*HMY(IMAX+I)
68    CONTINUE
      SUM2=0.5*DX*(XX(IMAX)*HMY(IMAX)+FX(1))
      SUM2=SUM2+0.5*(SL/2.-XX(N))*XX(N)*HMY(N)
      SUM2=SUM2+SIMP(FX,DX,N-IMAX)
      END IF
      HNR=2.*(-XX(IMAX)*XX(IMAX)*HMY(IMAX)+SUM1-SUM2)/(SL**4)
      WRITE (6,*) 'Yr''=',YR,'   Ydr''=',YDR,'   Nr''=',HNR
C
C  => Berechnung der Stabilitaet !
      VL=HNV/YV
      RL=(HNR-SM*XG/SL)/(YR-SM)
      CL=RL-VL
      WRITE (6,*)
      WRITE (6,*) 'VL=',VL,'   RL=',RL,'   CL=RL-VL=',CL
      IF (CL) 70,80,90
70    WRITE (6,*) 'Das Schiff ist gierinstabil !'
      GOTO 100
80    WRITE (6,*) 'Das Schiff ist grenzstabil !'
      GOTO 100
90    WRITE (6,*) 'Das Schiff ist gierstabil.'
C
100   CONTINUE
      END
  
C*********************************************************************
      FUNCTION SIMP(FX,DX,N)
C  Diese Funktion dient zur numerischen Integralrechnung der Funktion FX
C  zwischen den Stuetzstellen X(1) und X(N) mit der Simpsonschen Regel.
C  N: Anzahl der Stuetzstellen (N muss eine ungerade Zahl sein).
C  DX: Abstand zwischen Stuetzstellen.
      DIMENSION FX(N)
      N0=(N-1)/2
      SIMP=0.
      IF (N.EQ.1) GOTO 30
      DO 10 I=1,N0
      SIMP=SIMP+4.*FX(2*I)
10    CONTINUE
      DO 20 I=1,N0-1
      SIMP=SIMP+2.*FX(2*I+1)
20    CONTINUE
      SIMP=DX*(FX(1)+FX(N)+SIMP)/3.
30    END  
C*********************************************************************
      SUBROUTINE TUCK(X,B,BX,AX,DX,N0,AW,XW,VD,GML,FNH,DD,SITA)
C  Diese Subroutine dient zur Berechnung der Absenkung DD und der 
C  Vertrimmung SITA (grad) infolge der Geradeausfahrt nach dem 
C  Verfahren von E.O.TUCK.
      DIMENSION X(N0),B(N0),BX(N0),AX(N0)
      PI=3.1415927
      SUM=0.
      DO 10 I=1,N0
      DO 10 J=1,N0
      IF (I.NE.J) THEN
       SUM=SUM+BX(I)*AX(J)*ALOG(ABS(X(I)-X(J)))*DX*DX
      ELSE
       IF (I.NE.1.AND.I.NE.N0) THEN
        C1=BX(I)*AX(J-1)*ALOG(ABS(X(I)-X(J-1)))*DX*DX
        C2=BX(I)*AX(J+1)*ALOG(ABS(X(I)-X(J+1)))*DX*DX
        SUM=SUM+(C1+C2)/2.
       ELSE
        SUM=SUM+0.
       END IF
      END IF
10    CONTINUE
      DD=-FNH*FNH*SUM/(2.*PI*AW*SQRT(1.-FNH*FNH))
      SUM=0.
      DO 20 I=1,N0
      DO 20 J=1,N0
      IF (I.NE.J) THEN
       SUM=SUM+(B(I)+(X(I)-XW)*BX(I))*AX(J)*ALOG(ABS(X(I)-X(J)))*DX*DX
      ELSE 
       IF (I.NE.1.AND.I.NE.N0) THEN
        C3=(B(I)+(X(I)-XW)*BX(I))*AX(J-1)*ALOG(ABS(X(I)-X(J-1)))*DX*DX
        C4=(B(I)+(X(I)-XW)*BX(I))*AX(J+1)*ALOG(ABS(X(I)-X(J+1)))*DX*DX
        SUM=SUM+(C3+C4)/2.
       ELSE
        SUM=SUM+0.
       END IF
      END IF
20    CONTINUE
      SITA=FNH*FNH*SUM/(2.*PI*VD*GML*SQRT(1.-FNH*FNH))
      SITA=SITA*180./PI
      RETURN
      END
C*********************************************************************
      SUBROUTINE LEWIS(B,T,A,Y,Z,FAI,NN)
C  Diese Subroutine dient zur Berechnung der Koordinaten (Y,Z) der 
C  Punkte auf der Kontur des symmetrischen Lewis-Spantes.
C  B: Breite, T: Tiefgang, A: Spantflaeche, NN: Anzahl der Punkte
      DIMENSION Y(NN),Z(NN),FAI(NN)  
      ALFA=A/(B*T)
      BETA=T+0.5*B
      GAMA=T-0.5*B
      COE=0.75*BETA-0.25*SQRT(BETA**2+4.*B*T*(1.-4.*ALFA/3.14159))
      A1=0.5*GAMA/COE   
      A3=0.5*BETA/COE-1.
      DO 10 I=1,NN
      FAI(I)=(FLOAT(NN-I)/FLOAT(NN-1))*0.5*3.14159
      Y(I)=COE*((1-A1)*SIN(FAI(I))-A3*SIN(3*FAI(I)))
      Z(I)=COE*((1+A1)*COS(FAI(I))+A3*COS(3*FAI(I)))
10    CONTINUE
      RETURN
      END
C*********************************************************************
      SUBROUTINE GAUSS(N,N1,A,EPS,*)
C  Diese Subroutine dient zur Loesung des linearen Gleichungssystems.
C  Inputs sind Anzahl der Gleichungen N(N1=N+1), Koeffizientenmatrix 
C  A(I,J) (I=1,2,...,N, J=1,2,...,N1) und EPS (eine kleine Zahl, z.B. 
C  EPS=1.E-6). Output ist die Loesung A(I,N1) (I=1,2,...,N). Wenn das 
C  Gleichungssystem singulaer ist, "Return to" der '*' entsprechenden 
C  Reihe.
      DIMENSION A(N,N1)
      DO 10 K=1,N
      DO 20 IO=K,N
      IF (ABS(A(IO,K))-EPS) 20,20,100
20    CONTINUE
      RETURN 1
100   IF (IO.NE.K) THEN
      DO 30 J=K,N1
      T=A(K,J)
      A(K,J)=A(IO,J)
      A(IO,J)=T
30    CONTINUE
      END IF
200   T=1/A(K,K)
      DO 40 J=K,N
40    A(K,J+1)=T*A(K,J+1)
      IN=N-1
      IF (K.EQ.N) GOTO 300
      DO 10 I=IO,IN
      DO 10 J=K,N
10    A(I+1,J+1)=A(I+1,J+1)-A(I+1,K)*A(K,J+1)
300   DO 50 IK=2,N
      I=N1-IK
      DO 50 J=I,IN
50    A(I,N1)=A(I,N1)-A(I,J+1)*A(J+1,N1)
      RETURN
      END
C*********************************************************************
      SUBROUTINE CSS(B,T,A,H,C)
C   Diese Subroutine dient zur Berechnung der Versperrungskoeffizienten
C   C symmetrischer Spanten.
C   B: Breite, T: Tiefgang, A: Spantflaeche, H: Wassertiefe.
      DIMENSION Y(91),Z(91),FAI(91)
      DIMENSION YK(45),ZK(45),YQ(45),ZQ(45),YNK(45),ZNK(45),AA(45,46)
      DATA NN/91/,N/45/
C   Berechnung der Koordinaten der Punkte auf Kontur der "Lewis-Spanten"
C   NN: Anzahl der gesamten Punkte auf der Kontur
C   N: Anzahl der Kollokations- bzw. Quellpunkte  
      CALL LEWIS(B,T,A,Y,Z,FAI,NN)
C   Berechnung der Kollokationspunktkoordinaten und der Normalenvektor
      DO 10 I=1,N
      YK(I)=Y(2*I) 
      ZK(I)=Z(2*I)
      YNK(I)=-Z(2*I+1)+Z(2*I-1)
      ZNK(I)=Y(2*I+1)-Y(2*I-1)
10    CONTINUE
C   Berechnung der Quellpunktkoordinaten   
      II=5
      DO 20 I=1,N-1
      YQ(I)=(YK(I)+YK(I+1))/2.+II*(ZK(I)-ZK(I+1))/2.
      ZQ(I)=(ZK(I)+ZK(I+1))/2.+II*(YK(I+1)-YK(I))/2. 
20    CONTINUE
      YQ(N)=(YK(N)+Y(NN))/2.+II*(ZK(N)-Z(NN))/2.      
      ZQ(N)=Z(NN)
C   Berechnung der Koeffizientenmatrix
      PI=3.1415927
      DATA EPS/1.E-6/
      DO 30 I=1,N
      AA(I,N+1)=YNK(I)
      DO 30 J=1,N
      SH1=SINH(PI*(YK(I)-YQ(J))/H)
      SH2=SINH(PI*(YK(I)+YQ(J))/H)
      S1=SIN(PI*(ZK(I)-ZQ(J))/H)
      S2=SIN(PI*(ZK(I)+ZQ(J))/H)
      CH1=COSH(PI*(YK(I)-YQ(J))/H)
      CH2=COSH(PI*(YK(I)+YQ(J))/H)
      C1=COS(PI*(ZK(I)-ZQ(J))/H) 
      C2=COS(PI*(ZK(I)+ZQ(J))/H)
      IF (ABS(CH1-C1).LT.EPS.AND.ABS(CH2-C1).LT.EPS) THEN
      A1=SH1/(CH1-C2)-SH2/(CH2-C2) 
      A2=S2/(CH1-C2)-S2/(CH2-C2)
      ELSE
      A1=SH1/(CH1-C1)+SH1/(CH1-C2)-SH2/(CH2-C1)-SH2/(CH2-C2)
      A2=S1/(CH1-C1)+S2/(CH1-C2)-S1/(CH2-C1)-S2/(CH2-C2)
      END IF
      AA(I,J)=A1*YNK(I)+A2*ZNK(I)
30    CONTINUE
C   Loesung des linearen Gleichungssystems mit der Subroutine "GAUSS"
      CALL GAUSS(N,N+1,AA,EPS,*50)
      SUM=0.
      DO 40 I=1,N        
      SUM=SUM+AA(I,N+1)*YQ(I)
40    CONTINUE
      C=4.*SUM
      GOTO 60
50    WRITE (6,*) 'Das lineare Gleichungssystem ist singulaer!'
60    RETURN
      END
C*********************************************************************
      SUBROUTINE QGVS(V,R,FNH,C,B,B0,X,DX,L,N0,N,N1)
C  Diese Subroutine dient zur Berechnung der Queranstroemungsgeschwin-
C  digkeit VSTAR unter Einfluss des Versperrungseffektes.
C  V: Quergeschwindigkeit des Schiffes, R: Drehgeschwindigkeit des 
C  Schiffes. 
C  FNH: Mit Wassertiefe gebildete Froudezahl.
C  C: Versperrungskoeffizienten.
C  B,B0: Koeffizientenmatrix des linearen Gleichungssystems fuer die 
C  Queranstroemungsgeschwindigkeit VSTAR.
C  N: Anzahl der Spanten.
C  X: X-Koordinaten der Spanten. L: Schiffslaenge.
C  DX: Abstand zwischen den Spanten. 
C  Diese Subroutine gilt fuer den Fall -L/2<X(1)<X(2)...<X(N)<L/2 und
C  X(N)-X(N-1)=X(N-1)-X(N-2)=...=X(2)-X(1)=DX. X(0)=-L/2, X(N+1)=L/2. 
C  Eingaben sind V,R,FNH,C,X,DX,L,N0(N0=N-1),N,N1(N1=N+1).
C  (R=0 fuer reine Querbewegung, V=0 fuer reine Drehbewegung)
C  Ausgaben sind B(I,N1)=VSTAR(I) (I=1,2,...,N).
      DIMENSION C(N),X(N),B(N,N1),B0(N0,N)
      REAL L
      IF (R.EQ.0.) THEN
      KK=1
      ELSE IF (V.EQ.0.) THEN
      KK=2
      ELSE
      KK=3
      END IF
C  KK=1 fuer reine Querbewegung, KK=2 fuer reine Drehbewegung, KK=3 fuer
C  kombinierte Quer- und Drehbewegung.
* --------------------------------------------------------------------- *
      DO 10 K=1,N
      DO 10 I=1,N1
      B(K,I)=0.
10    CONTINUE
C  Berechnung der rechten Seite des linearen Gleichungssystems.   
      IF (KK-2) 11,12,13
11    DO 1 K=2,N-1
      B(K,N1)=V*DX
1     CONTINUE
      B(1,N1)=V*(DX+X(1)+L/2.)/2.
      B(N,N1)=V*(L/2.-X(N)+DX)/2.
      GOTO 20
12    DO 2 K=2,N-1
      B(K,N1)=R*X(K)*DX
2     CONTINUE
      B(1,N1)=R*(DX+X(1)+L/2.)*(DX+2.*X(1)-L/2.)/6.
      B(N,N1)=R*(L/2.-X(N)+DX)*(L/2.+2.*X(N)-DX)/6.
      GOTO 20
13    DO 3 K=2,N-1
      B(K,N1)=V*DX+R*X(K)*DX
3     CONTINUE
      B(1,N1)=V*(DX+X(1)+L/2.)/2.+R*(DX+X(1)+L/2.)*(DX+2.*X(1)-L/2.)/6.
      B(N,N1)=V*(L/2.-X(N)+DX)/2.+R*(L/2.-X(N)+DX)*(L/2.+2.*X(N)-DX)/6.
C
C  Berechnung der Koeffizienten B1 des linearen Gleichungssystems 
20    EP=1.E-2
      DO 21 K=2,N
      IF (ABS(C(K)-C(K-1))/C(K).LT.EP) THEN
      B(K,K-1)=DX/(6.*C(K))
      ELSE
      C1=(C(K)*C(K)-C(K-1)*C(K-1))/2.-C(K)*C(K-1)*ALOG(C(K)/C(K-1))
      C2=(C(K)-C(K-1))**3
      B(K,K-1)=DX*C1/C2
      END IF
21    CONTINUE
      DO 22 K=2,N-1
      C11=3.*C(K-1)*C(K-1)/2.-2.*C(K)*C(K-1)+C(K)*C(K)/2.
      C12=(C(K)-C(K-1))**3
      C13=C(K-1)*C(K-1)*ALOG(C(K)/C(K-1))
      C21=2.*C(K)*C(K+1)-C(K)*C(K)/2.-3.*C(K+1)*C(K+1)/2.
      C22=(C(K+1)-C(K))**3
      C23=C(K+1)*C(K+1)*ALOG(C(K+1)/C(K))
      IF (ABS(C(K)-C(K-1))/C(K).LT.EP.AND.ABS(C(K)-C(K+1))/C(K).LT.EP)
     &THEN 
      B(K,K)=2.*DX/(3.*C(K))
      ELSE IF (ABS(C(K)-C(K-1))/C(K).LT.EP) THEN
      B(K,K)=DX/(3.*C(K))+DX*(C21+C23)/C22
      ELSE IF (ABS(C(K)-C(K+1))/C(K).LT.EP) THEN
      B(K,K)=DX*(C11+C13)/C12+DX/(3.*C(K))
      ELSE
      B(K,K)=DX*(C11+C13)/C12+DX*(C21+C23)/C22
      END IF
22    CONTINUE
      C11=2.*C(1)*C(2)-C(1)*C(1)/2.-3.*C(2)*C(2)/2.
      C12=(C(2)-C(1))**3
      C13=C(2)*C(2)*ALOG(C(2)/C(1))
      B(1,1)=(X(1)+L/2.)/(2*C(1))+DX*(C11+C13)/C12
      C21=3.*C(N-1)*C(N-1)/2.-2.*C(N)*C(N-1)+C(N)*C(N)/2.
      C22=(C(N)-C(N-1))**3
      C23=C(N-1)*C(N-1)*ALOG(C(N)/C(N-1))
      B(N,N)=DX*(C21+C23)/C22+(L/2.-X(N))/(2.*C(N))
      DO 23 K=1,N-1
      IF (ABS(C(K)-C(K+1))/C(K).LT.EP) THEN
      B(K,K+1)=DX/(6.*C(K))
      ELSE
      C1=(C(K+1)*C(K+1)-C(K)*C(K))/2.-C(K)*C(K+1)*ALOG(C(K+1)/C(K))
      C2=(C(K+1)-C(K))**3
      B(K,K+1)=DX*C1/C2
      END IF
23    CONTINUE
C
C  Berechnung der Koeffizienten B2 und B=B1+B2 des linearen 
C  Gleichungssystems
      COE=SQRT(1.-FNH*FNH)/3.1415927   
      X1=X(1)+L/2.
      X2=X(2)+L/2.
      X3=X(3)+L/2.
      X4=X(4)+L/2.
      XN=L/2.-X(N)
      XN1=L/2.-X(N-1)
      XN2=L/2.-X(N-2)
      XN3=L/2.-X(N-3)
      DO 30 K=2,N-1
      DO 30 I=2,N-1
      IF (I.EQ.K) THEN
      B(K,I)=B(K,I)+COE*(4.*ALOG(2.))
      ELSE IF (I.EQ.K-1.OR.I.EQ.K+1) THEN
      B(K,I)=B(K,I)+COE*(4.5*ALOG(3.)-8.*ALOG(2.))
      ELSE IF (I.EQ.K-2.OR.I.EQ.K+2) THEN
      B(K,I)=B(K,I)+COE*(28.*ALOG(2.)-18.*ALOG(3.))
      ELSE
      S1=3.*(I-K)*(I-K)*ALOG(ABS(I-K)*DX)
      S2=2.*(I-K-1)*(I-K-1)*ALOG(ABS(I-K-1)*DX)
      S3=2.*(I-K+1)*(I-K+1)*ALOG(ABS(I-K+1)*DX)
      S4=0.5*(I-K+2)*(I-K+2)*ALOG(ABS(I-K+2)*DX)
      S5=0.5*(I-K-2)*(I-K-2)*ALOG(ABS(I-K-2)*DX)
      B(K,I)=B(K,I)+COE*(S1-S2-S3+S4+S5)
      END IF
30    CONTINUE
      S=X2*X2*ALOG(X2)/(X1*DX)-X2*ALOG(X1)/DX-X2*ALOG(DX)/X1
      B(1,1)=B(1,1)+COE*S
      S1=(1.+X2/X1)*ALOG(DX)
      S2=2.*X2*ALOG(2.*DX)/X1
      S11=0.5*X1*ALOG(X1)/DX
      S12=X2*X2*ALOG(X2)/(DX*X1)
      S13=0.5*X3*X3*ALOG(X3)/(DX*X1)
      B(1,2)=B(1,2)+COE*(S1-S2+S11-S12+S13)
      B(2,1)=B(2,1)+COE*(S1-S2+S11-S12+S13)
      S1=(1.+0.5*X2/X1)*ALOG(DX)
      S2=(2.+4.*X2/X1)*ALOG(2.*DX)
      S3=4.5*X2*ALOG(3.*DX)/X1
      S12=0.5*X2*X2*ALOG(X2)/(DX*X1)
      S13=X3*X3*ALOG(X3)/(DX*X1)
      S14=0.5*X4*X4*ALOG(X4)/(DX*X1)
      B(1,3)=B(1,3)+COE*(S2-S1-S3+S12-S13+S14)
      B(3,1)=B(3,1)+COE*(S2-S1-S3+S12-S13+S14)
      DO 31 I=4,N-1
      S=0.5*X2*I*I*ALOG(I*DX)/X1
      S1=(0.5+X2/X1)*(I-1)*(I-1)*ALOG((I-1)*DX)
      S2=(1.+0.5*X2/X1)*(I-2)*(I-2)*ALOG((I-2)*DX)
      S3=0.5*(I-3)*(I-3)*ALOG((I-3)*DX)
      S11=0.5*(X(I-1)+L/2.)*(X(I-1)+L/2.)*ALOG(X(I-1)+L/2.)/(DX*X1)
      S12=(X(I)+L/2.)*(X(I)+L/2.)*ALOG(X(I)+L/2.)/(DX*X1)
      S13=0.5*(X(I+1)+L/2.)*(X(I+1)+L/2.)*ALOG(X(I+1)+L/2.)/(DX*X1)
      B(1,I)=B(1,I)+COE*(S1-S2+S3-S+S11-S12+S13)
      B(I,1)=B(I,1)+COE*(S1-S2+S3-S+S11-S12+S13)
31    CONTINUE
      S=0.5*L*L*ALOG(L)/(X1*XN)
      S1=0.5*X2*(X(1)-L/2.)*(X(1)-L/2.)*ALOG(L/2.-X(1))/(DX*X1*XN)
      S2=0.5*(X(2)-L/2.)*(X(2)-L/2.)*ALOG(L/2.-X(2))/(DX*XN)
      S3=0.5*(L/2.+X(N-1))*(L/2.+X(N-1))*ALOG(L/2.+X(N-1))/(DX*X1)
      S4=0.5*XN1*(X(N)+L/2.)*(X(N)+L/2.)*ALOG(X(N)+L/2.)/(DX*X1*XN)
      S11=0.5*XN1*X2*(N-1)*(N-1)*ALOG((N-1)*DX)/(X1*XN)
      S12=0.5*(N-2)*(N-2)*(X2/X1+XN1/XN)*ALOG((N-2)*DX)
      S13=0.5*(N-3)*(N-3)*ALOG((N-3)*DX)
      B(1,N)=B(1,N)+COE*(S-S1+S2+S3-S4+S11-S12+S13)
      B(N,1)=B(N,1)+COE*(S-S1+S2+S3-S4+S11-S12+S13)
      DO 32 K=2,N-3
      S1=(L/2.-X(K))*(L/2.-X(K))*ALOG(L/2.-X(K))/(DX*XN)
      S2=0.5*(X(K-1)-L/2.)*(X(K-1)-L/2.)*ALOG(L/2.-X(K-1))/(DX*XN)
      S3=0.5*(X(K+1)-L/2.)*(X(K+1)-L/2.)*ALOG(L/2.-X(K+1))/(DX*XN)
      S11=(N-K)*(N-K)*(0.5+XN1/XN)*ALOG((N-K)*DX)
      S12=(N-K-1)*(N-K-1)*(1.+0.5*XN1/XN)*ALOG((N-K-1)*DX)
      S13=0.5*(N-K-2)*(N-K-2)*ALOG((N-K-2)*DX)
      S14=0.5*(N-K+1)*(N-K+1)*XN1*ALOG((N-K+1)*DX)/XN
      B(K,N)=B(K,N)+COE*(S3+S2-S1+S11-S12+S13-S14)
      B(N,K)=B(N,K)+COE*(S3+S2-S1+S11-S12+S13-S14)
32    CONTINUE
      S1=(1.+0.5*XN1/XN)*ALOG(DX)
      S2=(2.+4.*XN1/XN)*ALOG(2.*DX)
      S3=4.5*XN1*ALOG(3.*DX)/XN
      S11=0.5*XN1*XN1*ALOG(XN1)/(DX*XN)
      S12=XN2*XN2*ALOG(XN2)/(DX*XN)
      S13=0.5*XN3*XN3*ALOG(XN3)/(DX*XN)
      B(N-2,N)=B(N-2,N)+COE*(S2-S1-S3+S11-S12+S13)
      B(N,N-2)=B(N,N-2)+COE*(S2-S1-S3+S11-S12+S13)
      S1=(1.+XN1/XN)*ALOG(DX)
      S2=2.*XN1*ALOG(2.*DX)/XN
      S11=XN1*XN1*ALOG(XN1)/(DX*XN)
      S12=0.5*XN2*XN2*ALOG(XN2)/(DX*XN)
      S13=0.5*XN*ALOG(XN)/DX
      B(N-1,N)=B(N-1,N)+COE*(S1-S2-S11+S12+S13)
      B(N,N-1)=B(N,N-1)+COE*(S1-S2-S11+S12+S13)
      S=XN1*XN1*ALOG(XN1)/(DX*XN)-XN1*ALOG(DX)/XN-XN1*ALOG(XN)/DX
      B(N,N)=B(N,N)+COE*S
C
C  Loesung des linearen Gleichungssystems mit der Subroutine "Gauss"
C     
C  Wenn die Fahrgeschwindigkeit nicht gleich 0 ist, muss man am 
C  Hintersteven eine Kutta-Bedingung (glatte Abstoemung-Bedingung) 
C  ansetzen. Diese Bedingung wird erfuellt durch die Forderung, dass
C  am Hinterende die Wirbelstaerke gleich 0 ist. d.h. (C*VSTAR)'=0 an 
C  Stelle X=-L/2. naeherungsweise: C(1)*VSTAR(1)=C(2)*VSTAR(2) oder 
C  ALFA(1)=ALFA(2).     
C
      DO 40 I=1,N0
      B0(I,N)=B(I+1,N1)
      B0(I,1)=B(I+1,1)+B(I+1,2)
      DO 40 J=2,N0
      B0(I,J)=B(I+1,J+1)
40    CONTINUE
      CALL GAUSS(N0,N,B0,1.E-6,*90)
      DO 41 I=1,N0
      B(I+1,N1)=B0(I,N)
41    CONTINUE
      B(1,N1)=B0(1,N)
      DO 50 I=1,N
      B(I,N1)=B(I,N1)/C(I)
50    CONTINUE
      GO TO 100
90    WRITE (6,*) 'Das lineare Gleichungssystem ist singulaer!'
100   RETURN
      END
C*********************************************************************

