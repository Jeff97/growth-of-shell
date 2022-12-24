C234567890123456789012345678901234567890123456789012345678901234567890
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
C LOCAL ARRAYS
C ----------------------------------------------------------------
C A1 - Elastic tensor at the end of the increment
C BA1 - Left Cauchy-green elastic tensor at the end of the increment
C BA1B - Deviatoric left Cauchy-green elastic tensor at the end of the increment
C TRBA1B - Trace of BA1B
C ----------------------------------------------------------------
C
        DIMENSION A1(3,3), BA1(6), BA1B(6)
C
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, 
     &    FOUR=4.D0, SIX=6.D0)
        REAL  TotalT,MMOD,G11,G12,G21,G22,G33,
     &  Theta,SQMMOD,Pi,
     &  Lambda1St1,Lambda2St1,Lambda3St1,
     &  Lambda1St2,Lambda2St2,Lambda3St2,
     &  CarteG11St1,CarteG12St1,CarteG22St1,
     &  CarteG11St2,CarteG12St2,CarteG22St2
C ----------------------------------------------------------------
C UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C CANNOT BE USED FOR PLANE STRESS
C CANNOT BE USED FOR REDUCED INTEGRATE ELEMENTS
C HYBRID ELEMENTS ARE RECOMMENDED
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C STATEV(1) - Increment of G11 or Lambda1
C STATEV(2) - Increment of G23 or Lambda2
C STATEV(3) - Initial coordinate COORDS(1)
C STATEV(4) - Initial coordinate COORDS(2)
C ----------------------------------------------------------------
C
C ELASTIC PROPERTIES
        EMOD=PROPS(1)
        ENU=PROPS(2)
        C10=EMOD/(FOUR*(ONE+ENU))
        D1=SIX*(ONE-TWO*ENU)/EMOD
        Pi=3.14159265359
C
        IF (STATEV(3) .EQ. 0) THEN
                STATEV(3)=COORDS(1)
        END IF
        IF (STATEV(4) .EQ. 0) THEN
                STATEV(4)=COORDS(2)
        END IF
        IF (STATEV(5) .EQ. 0) THEN
                STATEV(5)=COORDS(3)
        END IF
C G is growth tensor
        MMOD=STATEV(3)*STATEV(3)+STATEV(4)*STATEV(4)
        SQMMOD=Sqrt(MMOD)
        ! IF (STATEV(3) .GE. 0 .AND. STATEV(4) .GE. 0) THEN
        !     Theta=ASin((STATEV(4)/SQMMOD))
        ! END IF
        ! IF (STATEV(3) .LE. 0 .AND. STATEV(4) .GE. 0) THEN
        !     Theta=ASin((-STATEV(3)/SQMMOD))+0.5*Pi
        ! END IF
        ! IF (STATEV(3) .LE. 0 .AND. STATEV(4) .LE. 0) THEN
        !     Theta=ASin((-STATEV(4)/SQMMOD))+Pi
        ! END IF
        ! IF (STATEV(3) .GE. 0 .AND. STATEV(4) .LE. 0) THEN
        !     Theta=ASin((STATEV(3)/SQMMOD))+1.5*Pi
        ! END IF
        Theta=ASin((STATEV(4)/SQMMOD))
        STATEV(9)=Theta
        TotalT=10.0
C
        Lambda1St1=1.0
        Lambda2St1=2.0*Cos(9.0/40.0*Pi*(-2.0+STATEV(5)))
        Lambda3St1=9.0/40.0*Pi*SQMMOD
C
        STATEV(6)=Lambda1St1
        STATEV(7)=Lambda2St1
        STATEV(8)=Lambda3St1
C
C
        CarteG11St1=(Lambda1St1*STATEV(3)*STATEV(3)
     &      +Lambda2St1*STATEV(4)*STATEV(4))/MMOD
        CarteG22St1=(Lambda1St1*STATEV(4)*STATEV(4)
     &      +Lambda2St1*STATEV(3)*STATEV(3))/MMOD
        CarteG12St1=(Lambda1St1-Lambda2St1)*STATEV(3)*STATEV(4)/MMOD
C
        G11=1.0+(CarteG11St1-1.0)*(TIME(1)+DTIME)/TotalT
        G22=1.0+(CarteG22St1-1.0)*(TIME(1)+DTIME)/TotalT
        G12=CarteG12St1*(TIME(1)+DTIME)/TotalT
        G21=G12
        G33=1.0+(Lambda3St1-1.0)*(TIME(1)+DTIME)/TotalT
C
C Elastic deformation tensor A=F.G^(-1)
C
        A1(1, 1)=(DFGRD1(1, 2)*G21-DFGRD1(1, 1)*G22)/
     &           (G12*G21-G11*G22)
        A1(1, 2)=(DFGRD1(1, 2)*G11-DFGRD1(1, 1)*G12)/
     &           (G11*G22-G12*G21)
        A1(2, 1)=(DFGRD1(2, 2)*G21-DFGRD1(2, 1)*G22)/
     &           (G12*G21-G11*G22)
        A1(2, 2)=(DFGRD1(2, 2)*G11-DFGRD1(2, 1)*G12)/
     &           (G11*G22-G12*G21)
        A1(3, 3)=DFGRD1(3, 3)/G33
!       IF (NTENS .EQ. 3) THEN
        A1(1, 3)=DFGRD1(1, 3)/G33
        A1(2, 3)=DFGRD1(2, 3)/G33
        A1(3, 1)=(DFGRD1(3, 2)*G21-DFGRD1(3, 1)*G22)/
     &           (G12*G21-G11*G22)
        A1(3, 2)=(DFGRD1(3, 2)*G11-DFGRD1(3, 1)*G12)/
     &           (G11*G22-G12*G21)
        ! END IF
C
C Determinent of the elastic deformation tensor DETA=det(A)
C
        DETA1=A1(1, 1)*A1(2, 2)*A1(3, 3)
     &     -A1(1, 2)*A1(2, 1)*A1(3, 3)
!       IF (NTENS .EQ. 3) THEN
        DETA1=DETA1+A1(1, 2)*A1(2, 3)*A1(3, 1)
     &     +A1(1, 3)*A1(3, 2)*A1(2, 1)
     &     -A1(1, 3)*A1(3,1)*A1(2, 2)
     &     -A1(2, 3)*A1(3, 2)*A1(1, 1)
        ! END IF
C
C Left Cauchy-Green strain tensor BA=A*A^(T)
C
        BA1(1)=A1(1, 1)**2+A1(1, 2)**2
        BA1(2)=A1(2, 1)**2+A1(2, 2)**2
        BA1(3)=A1(3, 3)**2
        BA1(4)=A1(1, 1)*A1(2, 1)+A1(1, 2)*A1(2, 2)
!       IF (NTENS .EQ. 3) THEN
          BA1(1)=BA1(1)+A1(1, 3)**2
          BA1(2)=BA1(2)+A1(2, 3)**2
          BA1(3)=BA1(3)+A1(3, 1)**2+A1(3, 2)**2
          BA1(4)=BA1(4)+A1(1, 3)*A1(2, 3)
          BA1(5)=A1(1, 1)*A1(3, 1)+A1(1, 2)*A1(3, 2)+A1(1, 3)*A1(3, 3)
          BA1(6)=A1(2, 1)*A1(3, 1)+A1(2, 2)*A1(3, 2)+A1(2, 3)*A1(3, 3)
        ! END IF
C
C Deviatoric left Cauchy-Green strain tensor BAB=DetBA^(-1/3)*BA
        DETBA1=BA1(1)*BA1(2)*BA1(3)
     &     -BA1(4)*BA1(4)*BA1(3)
!       IF (NTENS .EQ. 3) THEN
        DETBA1=DETBA1+BA1(4)*BA1(6)*BA1(5)
     &     +BA1(5)*BA1(6)*BA1(4)
     &     -BA1(5)*BA1(5)*BA1(2)
     &     -BA1(6)*BA1(6)*BA1(1)
        ! END IF
        BA1B(1)=DETBA1**(-ONE/THREE)*BA1(1)
        BA1B(2)=DETBA1**(-ONE/THREE)*BA1(2)
        BA1B(3)=DETBA1**(-ONE/THREE)*BA1(3)
        BA1B(4)=DETBA1**(-ONE/THREE)*BA1(4)
!       IF (NTENS .EQ. 3) THEN
        BA1B(5)=DETBA1**(-ONE/THREE)*BA1(5)
        BA1B(6)=DETBA1**(-ONE/THREE)*BA1(6)
        ! END IF
C
C Calculate Cauchy stress
C
        TRBA1B=BA1B(1)+BA1B(2)+BA1B(3)
        EG=TWO*C10/DETA1
        PR=TWO/D1*(DETA1-ONE)
        DO K1=1,NDI
        STRESS(K1)=EG*(BA1B(K1)-TRBA1B/THREE)+PR
        END DO
        DO K1=NDI+1,NDI+NSHR
        STRESS(K1)=EG*BA1B(K1)
        END DO
C
C Calculate the consistent Jacobian
C
        EG23=EG*TWO/THREE
        EK=TWO/D1*(TWO*DETA1-ONE)
        DDSDDE(1, 1)= EG23*(BA1B(1)+TRBA1B/THREE)+EK
        DDSDDE(1, 2)=-EG23*(BA1B(1)+BA1B(2)-TRBA1B/THREE)+EK
        DDSDDE(1, 3)=-EG23*(BA1B(1)+BA1B(3)-TRBA1B/THREE)+EK
        DDSDDE(1, 4)= EG23*BA1B(4)/TWO
        DDSDDE(2, 2)= EG23*(BA1B(2)+TRBA1B/THREE)+EK
        DDSDDE(2, 3)=-EG23*(BA1B(2)+BA1B(3)-TRBA1B/THREE)+EK
        DDSDDE(2, 4)= EG23*BA1B(4)/TWO
        DDSDDE(3, 3)= EG23*(BA1B(3)+TRBA1B/THREE)+EK
        DDSDDE(3, 4)=-EG23*BA1B(4)
        DDSDDE(4, 4)= EG*(BA1B(1)+BA1B(2))/TWO
!       IF (NTENS .EQ. 3) THEN
        DDSDDE(1, 5)= EG23*BA1B(5)/TWO
        DDSDDE(1, 6)=-EG23*BA1B(6)
        DDSDDE(2, 5)=-EG23*BA1B(5)
        DDSDDE(2, 6)= EG23*BA1B(6)/TWO
        DDSDDE(3, 5)= EG23*BA1B(5)/TWO
        DDSDDE(3, 6)= EG23*BA1B(6)/TWO
        DDSDDE(4, 5)= EG*BA1B(6)/TWO
        DDSDDE(4, 6)= EG*BA1B(5)/TWO
        DDSDDE(5, 5)= EG*(BA1B(1)+BA1B(3))/TWO
        DDSDDE(5, 6)= EG*BA1B(4)/TWO
        DDSDDE(6, 6)= EG*(BA1B(2)+BA1B(3))/TWO
        ! END IF
        DO K1=1, NTENS
          DO K2=1, K1-1
            DDSDDE(K1, K2)=DDSDDE(K2, K1)
          END DO
        END DO
C
        RETURN
       END