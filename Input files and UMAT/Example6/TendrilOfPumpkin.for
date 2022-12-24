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
        DIMENSION A1(3,3), BA1(6), BA1B(6),NG11Z0(10201),
     &  NG22Z0(10201),NG11Z1(10201),NG22Z1(10201),
     &  NodeCoordTh(10201),NodeCoordZ(10201),NG12Z0(10201),
     &  NG21Z0(10201),NG12Z1(10201),NG21Z1(10201)
C
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, 
     &    FOUR=4.D0, SIX=6.D0)
        REAL  TotalT,MMOD,Theta,SQMMOD,Pi,
     &  FinalG11,FinalG12,FinalG13,
     &  FinalG21,FinalG22,FinalG23,
     &  FinalG31,FinalG32,FinalG33,
     &  LamrR, LamrT, LamrZ, 
     &  LamtR, LamtT, LamtZ,
     &  LamzR, LamzT, LamzZ,
     &  G11,G12,G13,G21,G22,G23,G31,G32,G33

        double precision NodeCoordTh,NodeCoordZ,NG11Z0,NG22Z0,
     &  NG22Z1,NG11Z1, NG12Z0,NG21Z0, NG21Z1,NG12Z1

        SAVE NodeCoordTh
        SAVE NodeCoordZ
        SAVE NG11Z0
        SAVE NG22Z0
        SAVE NG22Z1
        SAVE NG11Z1
        SAVE NG12Z0
        SAVE NG21Z0
        SAVE NG21Z1
        SAVE NG12Z1
        data iread /1/
        SAVE iread
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
C Calculate R and Theta
        MMOD=STATEV(3)*STATEV(3)+STATEV(4)*STATEV(4)
        SQMMOD=Sqrt(MMOD)
        IF (STATEV(3) .GE. 0 .AND. STATEV(4) .GE. 0) THEN
            Theta=ASin((STATEV(4)/SQMMOD))
        END IF
        IF (STATEV(3) .LE. 0 .AND. STATEV(4) .GE. 0) THEN
            Theta=ASin((-STATEV(3)/SQMMOD))+0.5*Pi
        END IF
        IF (STATEV(3) .LE. 0 .AND. STATEV(4) .LE. 0) THEN
            Theta=ASin((-STATEV(4)/SQMMOD))+Pi
        END IF
        IF (STATEV(3) .GE. 0 .AND. STATEV(4) .LE. 0) THEN
            Theta=ASin((STATEV(3)/SQMMOD))+1.5*Pi
        END IF
        STATEV(6)=Theta
        TotalT=10.0

        LamrR=1.0

C Read data file, call Mutex to ensure thread safety
        call MutexInit( 1 )
        call MutexLock( 1 )

        IF (iread .EQ. 1) THEN

        open(201,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NodeCoordTh.CSV',status="old")
        read(201,*) NodeCoordTh
        close(201)

        open(202,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NodeCoordZ.CSV',status="old")
        read(202,*) NodeCoordZ
        close(202)

        open(301,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG11Z0.CSV',status="old")
        read(301,*) NG11Z0
        close(301)

        open(302,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG11Z1.CSV',status="old")
        read(302,*) NG11Z1
        close(302)

        open(303,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG12Z0.CSV',status="old")
        read(303,*) NG12Z0
        close(303)

        open(304,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG12Z1.CSV',status="old")
        read(304,*) NG12Z1
        close(304)

        open(301,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG21Z0.CSV',status="old")
        read(301,*) NG21Z0
        close(301)

        open(302,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG21Z1.CSV',status="old")
        read(302,*) NG21Z1
        close(302)

        open(303,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG22Z0.CSV',status="old")
        read(303,*) NG22Z0
        close(303)

        open(304,FILE='F:\OneDrive - mail.scut.edu.cn\'//
     &  'AbaqusRemote\Shell\HelicalRod\NG22Z1.CSV',status="old")
        read(304,*) NG22Z1
        close(304)

        iread=2
        END IF
        call MutexUnLock( 1 )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C calculate growth functions in Gauss points
        du=2*Pi/100.0
        dv=8.0/100.0
        iNx=100
          ! last time is 8181
          do K1 = 1, 10201

        CoordLeft=NodeCoordTh(K1)-du/2.0
        CoordRight=NodeCoordTh(K1)+du/2.0
        CoordBottom=NodeCoordZ(K1)-dv/2.0
        CoordTop=NodeCoordZ(K1)+dv/2.0
C Select 4 integration points near a node
            if ( STATEV(6) .GE. CoordLeft .and.
     &           STATEV(6) .LE. CoordRight .and.
     &           STATEV(5) .GE. CoordBottom .and.
     &           STATEV(5) .LE. CoordTop ) THEN

C point NE
              if (STATEV(6)>NodeCoordTh(K1) .and. 
     &               STATEV(5)>NodeCoordZ(K1) ) THEN
              xi =STATEV(6)-NodeCoordTh(K1)
              eta=STATEV(5)-NodeCoordZ(K1)
              STATEV(10) = NG11Z0(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z0(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG11Z0(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG11Z0(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(11) = NG11Z1(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z1(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG11Z1(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG11Z1(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(12) = NG12Z0(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z0(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG12Z0(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG12Z0(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(13) = NG12Z1(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z1(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG12Z1(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG12Z1(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(14) = NG21Z0(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z0(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG21Z0(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG21Z0(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(15) = NG21Z1(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z1(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG21Z1(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG21Z1(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(16) = NG22Z0(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z0(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG22Z0(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG22Z0(K1+iNx+2)*(xi/du)*(eta/dv)
              STATEV(17) = NG22Z1(K1)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z1(K1+1)*(xi/du)*(1-eta/dv)
     &        + NG22Z1(K1+iNx+1)*(1-xi/du)*(eta/dv)
     &        + NG22Z1(K1+iNx+2)*(xi/du)*(eta/dv)
              end if 
C point SE
              if (STATEV(6)>NodeCoordTh(K1) .and. 
     &               STATEV(5)<NodeCoordZ(K1) ) THEN
              xi =STATEV(6)-NodeCoordTh(K1-iNx-1)
              eta=STATEV(5)-NodeCoordZ(K1-iNx-1)
              STATEV(10) = NG11Z0(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z0(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG11Z0(K1)*(1-xi/du)*(eta/dv)
     &        + NG11Z0(K1+1)*(xi/du)*(eta/dv)
              STATEV(11) = NG11Z1(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z1(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG11Z1(K1)*(1-xi/du)*(eta/dv)
     &        + NG11Z1(K1+1)*(xi/du)*(eta/dv)
              STATEV(12) = NG12Z0(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z0(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG12Z0(K1)*(1-xi/du)*(eta/dv)
     &        + NG12Z0(K1+1)*(xi/du)*(eta/dv)
              STATEV(13) = NG12Z1(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z1(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG12Z1(K1)*(1-xi/du)*(eta/dv)
     &        + NG12Z1(K1+1)*(xi/du)*(eta/dv)
              STATEV(14) = NG21Z0(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z0(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG21Z0(K1)*(1-xi/du)*(eta/dv)
     &        + NG21Z0(K1+1)*(xi/du)*(eta/dv)
              STATEV(15) = NG21Z1(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z1(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG21Z1(K1)*(1-xi/du)*(eta/dv)
     &        + NG21Z1(K1+1)*(xi/du)*(eta/dv)
              STATEV(16) = NG22Z0(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z0(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG22Z0(K1)*(1-xi/du)*(eta/dv)
     &        + NG22Z0(K1+1)*(xi/du)*(eta/dv)
              STATEV(17) = NG22Z1(K1-iNx-1)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z1(K1-iNx)*(xi/du)*(1-eta/dv)
     &        + NG22Z1(K1)*(1-xi/du)*(eta/dv)
     &        + NG22Z1(K1+1)*(xi/du)*(eta/dv)
              end if 
C point NW
              if (STATEV(6)<NodeCoordTh(K1) .and. 
     &               STATEV(5)>NodeCoordZ(K1) ) THEN
              xi =STATEV(6)-NodeCoordTh(K1-1)
              eta=STATEV(5)-NodeCoordZ(K1-1)
              STATEV(10) = NG11Z0(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z0(K1)*(xi/du)*(1-eta/dv)
     &        + NG11Z0(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG11Z0(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(11) = NG11Z1(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z1(K1)*(xi/du)*(1-eta/dv)
     &        + NG11Z1(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG11Z1(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(12) = NG12Z0(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z0(K1)*(xi/du)*(1-eta/dv)
     &        + NG12Z0(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG12Z0(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(13) = NG12Z1(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z1(K1)*(xi/du)*(1-eta/dv)
     &        + NG12Z1(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG12Z1(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(14) = NG21Z0(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z0(K1)*(xi/du)*(1-eta/dv)
     &        + NG21Z0(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG21Z0(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(15) = NG21Z1(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z1(K1)*(xi/du)*(1-eta/dv)
     &        + NG21Z1(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG21Z1(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(16) = NG22Z0(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z0(K1)*(xi/du)*(1-eta/dv)
     &        + NG22Z0(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG22Z0(K1+iNx+1)*(xi/du)*(eta/dv)
              STATEV(17) = NG22Z1(K1-1)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z1(K1)*(xi/du)*(1-eta/dv)
     &        + NG22Z1(K1+iNx)*(1-xi/du)*(eta/dv)
     &        + NG22Z1(K1+iNx+1)*(xi/du)*(eta/dv)
              end if 
C point SW
              if (STATEV(6)<NodeCoordTh(K1) .and. 
     &               STATEV(5)<NodeCoordZ(K1) ) THEN
              xi =STATEV(6)-NodeCoordTh(K1-iNx-2)
              eta=STATEV(5)-NodeCoordZ(K1-iNx-2)
              STATEV(10) = NG11Z0(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z0(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG11Z0(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG11Z0(K1)*(xi/du)*(eta/dv)
              STATEV(11) = NG11Z1(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG11Z1(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG11Z1(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG11Z1(K1)*(xi/du)*(eta/dv)
              STATEV(12) = NG12Z0(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z0(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG12Z0(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG12Z0(K1)*(xi/du)*(eta/dv)
              STATEV(13) = NG12Z1(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG12Z1(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG12Z1(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG12Z1(K1)*(xi/du)*(eta/dv)
              STATEV(14) = NG21Z0(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z0(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG21Z0(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG21Z0(K1)*(xi/du)*(eta/dv)
              STATEV(15) = NG21Z1(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG21Z1(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG21Z1(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG21Z1(K1)*(xi/du)*(eta/dv)
              STATEV(16) = NG22Z0(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z0(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG22Z0(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG22Z0(K1)*(xi/du)*(eta/dv)
              STATEV(17) = NG22Z1(K1-iNx-2)*(1-xi/du)*(1-eta/dv)
     &        + NG22Z1(K1-iNx-1)*(xi/du)*(1-eta/dv)
     &        + NG22Z1(K1-1)*(1-xi/du)*(eta/dv)
     &        + NG22Z1(K1)*(xi/du)*(eta/dv)
              end if 

            end if
          end do

        LamtT=STATEV(10)+(SQMMOD-4.0)*STATEV(11)
        LamtZ=STATEV(12)+(SQMMOD-4.0)*STATEV(13)
        LamzT=STATEV(14)+(SQMMOD-4.0)*STATEV(15)
        LamzZ=STATEV(16)+(SQMMOD-4.0)*STATEV(17)

        STATEV(21) = LamtT
        STATEV(22) = LamtZ
        STATEV(23) = LamzT
        STATEV(24) = LamzZ
C
C
        FinalG11=(LamrR*STATEV(3)*STATEV(3)
     &      +LamtT*STATEV(4)*STATEV(4))/MMOD
        FinalG22=(LamrR*STATEV(4)*STATEV(4)
     &      +LamtT*STATEV(3)*STATEV(3))/MMOD
        FinalG12=(LamrR-LamtT)*STATEV(3)*STATEV(4)/MMOD
        FinalG21=FinalG12

        FinalG31=-(LamzT*STATEV(4))/SQMMOD
        FinalG32=(LamzT*STATEV(3))/SQMMOD

        FinalG13=-(LamtZ*STATEV(4))/SQMMOD
        FinalG23=(LamtZ*STATEV(3))/SQMMOD

        FinalG33=LamzZ
! During the calculation the G is related to time TIME(1)
        G11=1.0+(FinalG11-1.0)*(TIME(1)+DTIME)/TotalT
        G12=FinalG12*(TIME(1)+DTIME)/TotalT
        G13=FinalG13*(TIME(1)+DTIME)/TotalT
        G21=FinalG21*(TIME(1)+DTIME)/TotalT
        G22=1.0+(FinalG22-1.0)*(TIME(1)+DTIME)/TotalT
        G23=FinalG23*(TIME(1)+DTIME)/TotalT
        G31=FinalG31*(TIME(1)+DTIME)/TotalT
        G32=FinalG32*(TIME(1)+DTIME)/TotalT
        G33=1.0+(FinalG33-1.0)*(TIME(1)+DTIME)/TotalT
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