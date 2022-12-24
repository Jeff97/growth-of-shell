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
        REAL  TotalT,MMOD,Theta,SQMMOD,Pi,
     &  FinalG11St1,FinalG12St1,FinalG13St1,
     &  FinalG21St1,FinalG22St1,FinalG23St1,
     &  FinalG31St1,FinalG32St1,FinalG33St1,
     &  LamrRSt1, LamrTSt1, LamrZSt1, 
     &  LamtRSt1, LamtTSt1, LamtZSt1,
     &  LamzRSt1, LamzTSt1, LamzZSt1,
     &  FinalG11St2,FinalG12St2,FinalG13St2,
     &  FinalG21St2,FinalG22St2,FinalG23St2,
     &  FinalG31St2,FinalG32St2,FinalG33St2,
     &  LamrRSt2, LamrTSt2, LamrZSt2, 
     &  LamtRSt2, LamtTSt2, LamtZSt2,
     &  LamzRSt2, LamzTSt2, LamzZSt2,
     &  G11,G12,G13,G21,G22,G23,G31,G32,G33
C
C ----------------------------------------------------------------
C UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C CANNOT BE USED FOR PLANE STRESS
C CANNOT BE USED FOR REDUCED INTEGRATE ELEMENTS
C HYBRID ELEMENTS ARE RECOMMENDED
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C STATEV(1) - Increment of G11 or LambdaR
C STATEV(2) - Increment of G23 or LambdaT
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
        Theta=ASin((STATEV(4)/SQMMOD))

        ! STATEV(9)=Theta
        LamrRSt1=1.0
        LamtTSt1=0.8
        LamtZSt1=0.0
        LamzTSt1=0.0
        LamzZSt1=1.0

        LamrRSt2=1.0
        LamtTSt2=(((8.0-SQMMOD)/Pi-(4.0*(-4.0+SQMMOD))/
     &   (1.0+(Pi+2.0*Theta)**2.0)))/(8.0*Sqrt(2.0))
        LamtZSt2=-((Pi**2.0*(-8.0+SQMMOD)+4.0*Pi
     &   *(-4.0+SQMMOD+(-8.0+SQMMOD)*Theta)
     &   +(-8.0+SQMMOD)*(1.0+4.0*Theta**2.0))/
     &   (4.0*Sqrt(2.0)*Sqrt(1.0+(Pi+2.0*Theta)**2.0))) + 
     &   ((-4+SQMMOD)*Sqrt(1+(Pi+2.0*Theta)**2))/(4.0*Sqrt(2.0))
        LamzTSt2=-(((1.0+Pi*(-4.0+SQMMOD)+(Pi+2.0*Theta)**2.0))
     &   /(2.0*Sqrt(2.0)*Pi*(1.0+(Pi+2.0*Theta)**2.0))) + 
     &   (-4.0+SQMMOD)/(8.0*Sqrt(2.0)*Pi)
        LamzZSt2=(1.0+Pi*(-4.0+SQMMOD)+(Pi+2.0*Theta)**2.0)
     &   /(Sqrt(2.0)*Sqrt(1.0+(Pi+2.0*Theta)**2.0))

!         LamrRSt2=1.0
!         LamtTSt2=(8.0-4*Pi*(-4+SQMMOD)-SQMMOD-(-8+SQMMOD)*
!      &  Cosh(2.0*ASinh(Pi+2.0*Theta)))/(8.0*Sqrt(2.0)*Pi
!      &  *(1.0+(Pi+2*Theta)**2))
!         LamtZSt2=(8.0-4*Pi*(-4+SQMMOD)-SQMMOD-(-8+SQMMOD)
!      &  *Cosh(2.0*ASinh(Pi+2.0*Theta)))/(4*Sqrt(2.0)
!      &  *Sqrt(1.0+(Pi+2*Theta)**2))
!         LamzTSt2=-((1.0+Pi*(-4+SQMMOD)+Cosh(2.0*ASinh(Pi+2*Theta)))/
!      &  (2.0*Sqrt(2.0)*Pi*(1+(Pi+2*Theta)**2)))
!         LamzZSt2=(1.0+Pi*(-4+SQMMOD)+Cosh(2.0*ASinh(Pi+2.0*Theta)))/
!      &  (Sqrt(2.0)*Sqrt(1.0+(Pi+2*Theta)**2))

        TotalT=10.0

        STATEV(6)=LamrRSt2
        STATEV(7)=LamtTSt2
        STATEV(8)=LamtZSt2
        STATEV(9)=LamzTSt2
        STATEV(10)=LamzZSt2
C
C Stage 1
        FinalG11St1=(LamrRSt1*STATEV(3)*STATEV(3)
     &      +LamtTSt1*STATEV(4)*STATEV(4))/MMOD
        FinalG22St1=(LamrRSt1*STATEV(4)*STATEV(4)
     &      +LamtTSt1*STATEV(3)*STATEV(3))/MMOD
        FinalG12St1=(LamrRSt1-LamtTSt1)*STATEV(3)*STATEV(4)/MMOD
        FinalG21St1=FinalG12St1

        FinalG31St1=-(LamzTSt1*STATEV(4))/SQMMOD
        FinalG32St1=(LamzTSt1*STATEV(3))/SQMMOD

        FinalG13St1=-(LamtZSt1*STATEV(4))/SQMMOD
        FinalG23St1=(LamtZSt1*STATEV(3))/SQMMOD

        FinalG33St1=LamzZSt1
C Stage 2
        FinalG11St2=(LamrRSt2*STATEV(3)*STATEV(3)
     &      +LamtTSt2*STATEV(4)*STATEV(4))/MMOD
        FinalG22St2=(LamrRSt2*STATEV(4)*STATEV(4)
     &      +LamtTSt2*STATEV(3)*STATEV(3))/MMOD
        FinalG12St2=(LamrRSt2-LamtTSt2)*STATEV(3)*STATEV(4)/MMOD
        FinalG21St2=FinalG12St2

        FinalG31St2=-(LamzTSt2*STATEV(4))/SQMMOD
        FinalG32St2=(LamzTSt2*STATEV(3))/SQMMOD

        FinalG13St2=-(LamtZSt2*STATEV(4))/SQMMOD
        FinalG23St2=(LamtZSt2*STATEV(3))/SQMMOD

        FinalG33St2=LamzZSt2
! During the calculation the G is related to time TIME(1)
      IF ( (TIME(1)+DTIME) .LE. 5.0) THEN
        G11=1.0+(FinalG11St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
        G12=FinalG12St1*2.0*(TIME(1)+DTIME)/TotalT
        G13=FinalG13St1*2.0*(TIME(1)+DTIME)/TotalT
        G21=FinalG21St1*2.0*(TIME(1)+DTIME)/TotalT
        G22=1.0+(FinalG22St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
        G23=FinalG23St1*2.0*(TIME(1)+DTIME)/TotalT
        G31=FinalG31St1*2.0*(TIME(1)+DTIME)/TotalT
        G32=FinalG32St1*2.0*(TIME(1)+DTIME)/TotalT
        G33=1.0+(FinalG33St1-1.0)*2.0*(TIME(1)+DTIME)/TotalT
      ELSE
        G11=2.0*FinalG11St1-FinalG11St2+(FinalG11St2-FinalG11St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G12=2.0*FinalG12St1-FinalG12St2+(FinalG12St2-FinalG12St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G13=2.0*FinalG13St1-FinalG13St2+(FinalG13St2-FinalG13St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G21=2.0*FinalG21St1-FinalG21St2+(FinalG21St2-FinalG21St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G22=2.0*FinalG22St1-FinalG22St2+(FinalG22St2-FinalG22St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G23=2.0*FinalG23St1-FinalG23St2+(FinalG23St2-FinalG23St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G31=2.0*FinalG31St1-FinalG31St2+(FinalG31St2-FinalG31St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G32=2.0*FinalG32St1-FinalG32St2+(FinalG32St2-FinalG32St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
        G33=2.0*FinalG33St1-FinalG33St2+(FinalG33St2-FinalG33St1)
     &  *2.0*(TIME(1)+DTIME)/TotalT
      END IF
C
C Elastic deformation tensor A=F.G^(-1)
C
!         A1(1, 1)=(DFGRD1(1, 2)*G21-DFGRD1(1, 1)*G22)/
!      &           (G12*G21-G11*G22)
        A1(1,1)=(DFGRD1(1,3)*G22*G31-DFGRD1(1,2)*G23*G31
     &  -DFGRD1(1,3)*G21*G32+DFGRD1(1,1)*G23*G32+DFGRD1(1,2)*G21*G33
     &  -DFGRD1(1,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(1,2)=(DFGRD1(1,3)*G12*G31-DFGRD1(1,2)*G13*G31
     &  -DFGRD1(1,3)*G11*G32+DFGRD1(1,1)*G13*G32+DFGRD1(1,2)*G11*G33
     &  -DFGRD1(1,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        A1(2,1)=(DFGRD1(2,3)*G22*G31-DFGRD1(2,2)*G23*G31
     &  -DFGRD1(2,3)*G21*G32+DFGRD1(2,1)*G23*G32+DFGRD1(2,2)*G21*G33
     &  -DFGRD1(2,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(2,2)=(DFGRD1(2,3)*G12*G31-DFGRD1(2,2)*G13*G31
     &  -DFGRD1(2,3)*G11*G32+DFGRD1(2,1)*G13*G32+DFGRD1(2,2)*G11*G33
     &  -DFGRD1(2,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
        A1(3,3)=(DFGRD1(3,3)*G12*G21-DFGRD1(3,2)*G13*G21
     &  -DFGRD1(3,3)*G11*G22+DFGRD1(3,1)*G13*G22+DFGRD1(3,2)*G11*G23
     &  -DFGRD1(3,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
!       IF*(NTENS*.EQ.*3)*THEN
        A1(1,3)=(DFGRD1(1,3)*G12*G21-DFGRD1(1,2)*G13*G21
     &  -DFGRD1(1,3)*G11*G22+DFGRD1(1,1)*G13*G22+DFGRD1(1,2)*G11*G23
     &  -DFGRD1(1,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(2,3)=(DFGRD1(2,3)*G12*G21-DFGRD1(2,2)*G13*G21
     &  -DFGRD1(2,3)*G11*G22+DFGRD1(2,1)*G13*G22+DFGRD1(2,2)*G11*G23
     &  -DFGRD1(2,1)*G12*G23)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(3,1)=(DFGRD1(3,3)*G22*G31-DFGRD1(3,2)*G23*G31
     &  -DFGRD1(3,3)*G21*G32+DFGRD1(3,1)*G23*G32+DFGRD1(3,2)*G21*G33
     &  -DFGRD1(3,1)*G22*G33)/(G13*G22*G31-G12*G23*G31
     &  -G13*G21*G32+G11*G23*G32+G12*G21*G33-G11*G22*G33)
        A1(3,2)=(DFGRD1(3,3)*G12*G31-DFGRD1(3,2)*G13*G31
     &  -DFGRD1(3,3)*G11*G32+DFGRD1(3,1)*G13*G32+DFGRD1(3,2)*G11*G33
     &  -DFGRD1(3,1)*G12*G33)/(-G13*G22*G31+G12*G23*G3
     &  1+G13*G21*G32-G11*G23*G32-G12*G21*G33+G11*G22*G33)
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