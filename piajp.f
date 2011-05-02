! program using Bosted's routines to calculate pion production cross sections
! based on fits to SLAC data to calculate rates of 
! a) pi+ and pi- into the BigCal
! b) gammas from pi0 decays entering BigCal
! copy necessary subroutines from pi0.f, write my own program to calculate 
! rates into BigCal

      program piajp
      implicit none 
      real e0
      real pmin,pmax,p,dp
      integer ip
      real cald
      integer npbin
      real xcal,xcalmin,xcalmax,xcal0,dxcal,zcal,zcalmin,zcalmax,zcal0
      real dzcal
      real ycal,ycalmin,ycalmax,dycal,rcal
      real theta,phi
      integer nxbin/32/,nybin/56/
      integer ix,iy
      real calangledeg/68.0/
      real calwidth,calheight
      real rl,radlen
      real sigpip,sigpim,sigpi0,ratepip,ratepim,ratepi0
      real tgtinversebarns
      real elrate
      real PI/3.14159265359/
      OPEN(unit=8,FILE='pical.out')
      OPEN(unit=9,FILE='pixy.out');
      WRITE(8,FMT='(A6,3A15)') 'p(GeV)','pi+, Hz','pi-, Hz','pi0, Hz'
      radlen=.03
      rl = (0.015 + radlen/2.) * 100.
      e0 = 6.0
      cald=440.
      pmin=.05
      pmax=1.35
      npbin=13
      dp = (pmax - pmin)/float(npbin)
      calwidth=121.984
      calheight=218.272
      xcal0 = cald*sin(calangledeg*PI/180.)
      zcal0 = cald*cos(calangledeg*PI/180.)
      xcalmin = xcal0 - .5*calwidth*cos(calangledeg*PI/180.)
      xcalmax = xcal0 + .5*calwidth*cos(calangledeg*PI/180.)
      ycalmin = -.5*calheight
      ycalmax =  .5*calheight
      zcalmin = zcal0 - .5*calwidth*sin(calangledeg*PI/180.)
      zcalmax = zcal0 + .5*calwidth*sin(calangledeg*PI/180.)
      dxcal = (xcalmax - xcalmin)/float(nxbin)
      dycal = (ycalmax - ycalmin)/float(nybin)
      dzcal = (zcalmax - zcalmin)/float(nxbin)
      tgtinversebarns = .63631E-9
      elrate = 6.25E14
      do ip=1,npbin
         p = pmin + (float(ip)-.5)*dp
         ratepip=0.
         ratepim=0.
         ratepi0=0.
         do ix=1,nxbin
            do iy=1,nybin
               xcal = xcalmin + (float(ix) - .5)*dxcal
               ycal = ycalmin + (float(iy) - .5)*dycal
               zcal = zcalmax - (float(ix) - .5)*dzcal
               rcal = sqrt(xcal**2 + ycal**2 + zcal**2)
               theta = acos(zcal/rcal)
               phi   = atan2(ycal,xcal)
               
               call WISER_ALL_SIG(e0,p,theta*180./PI,rl,1,sigpip)
               call WISER_ALL_SIG(e0,p,theta*180./PI,rl,2,sigpim)
               sigpi0 = .5*(sigpip + sigpim)
! WISER_ALL_SIG returns d^2 sigma/ dp domega for a given p, theta
!     in nanobarns/GeV/sr
               sigpim = sigpim*dp*dxcal*dycal/rcal**2
               sigpip = sigpip*dp*dxcal*dycal/rcal**2
               sigpi0 = sigpi0*dp*dxcal*dycal/rcal**2
               ratepip = ratepip + sigpip*tgtinversebarns*elrate
               ratepim = ratepim + sigpim*tgtinversebarns*elrate
               ratepi0 = ratepi0 + sigpi0*tgtinversebarns*elrate
               
               write(9,fmt='(2I4,E15.3E2)') ix,iy,(sigpip+sigpim)*
     &              tgtinversebarns*elrate
            enddo
         enddo
         write(*,*) p,ratepip,ratepim,ratepi0
         write(8,FMT='(F6.3,3E15.3E2)') p,ratepip,ratepim,ratepi0    
      enddo
      end

C SUBROUTINE TO DECAY A PARTICLE INTO TWO OTHER PARTICLES
C AT angle COSTH in THE C.M. FRAME WHERE THE DIRECTION OF THE
C INITIAL PARTICLE NEED NOT BE ALONG THE Z-AXIS. PHI SHOULD BE 0 TO 2 PI
C vectors are (E, px, py, pz)
      SUBROUTINE DECAY (V1,V2,V3,M1,M2,M3,COSTH,PHI)
      IMPLICIT NONE
      REAL V1(4),V2(4),V3(4),KV(4),M1,M2,M3,V1P,G,K,COSTH,SINTH,PHI
      REAL CHK,CHK1,CHK2,KCHK,BETA,VT(4),VTT(4),VDM
! Find magnitude of momentum of initial particle
      V1P=SQRT(V1(2)**2+V1(3)**2+V1(4)**2)
      G=V1(1)/M1 ! gamma = E/M
      BETA=V1P/V1(1) ! v/c = beta = p/E

! Find c.m. momentum and pick decay isotropically
      K=SQRT(( ((M1**2-M2**2-M3**2)/2.)**2 - (M2*M3)**2 )/M1**2)
      KCHK=SQRT((M1**2-(M2+M3)**2)*(M1**2-(M2-M3)**2))/2./M1
      IF(ABS(K-KCHK).GT.0.001) WRITE(6,'(1X,''ERROR,K,KCHK='',2F10.4)')
     >  K,KCHK
      IF(ABS(SQRT(K*K+M2*M2)+SQRT(K*K+M3*M3)-M1).GT.0.001)
     >  WRITE(6,'(1X,''ERROR IN CM ENERGY'')')
      SINTH=SQRT(1.-COSTH**2)
! Put k into vector for first decay particle
      KV(1)=SQRT(K*K+M2*M2)
      KV(2)=K*SINTH*COS(PHI)
      KV(3)=K*SINTH*SIN(PHI)
      KV(4)=K*COSTH
! Transform k to lab frame for first decay particle
! First find vector as though decay were along Z axis
      V2(1)=G*(KV(1)+BETA*KV(4))
      VT(2)=KV(2)
      VT(3)=KV(3)
      VT(4)=G*(BETA*KV(1)+KV(4))
! Now rotate vector in xz plane
      VDM=SQRT(V1(2)**2+V1(4)**2)
      VTT(2)= VT(2)*V1(4)/VDM+VT(4)*V1(2)/VDM
      VTT(3)= VT(3)
      VTT(4)=-VT(2)*V1(2)/VDM+VT(4)*V1(4)/VDM
! Now rotate vector in yz plane
      VDM=SQRT(V1(3)**2+V1(4)**2)
      V2(2)= VTT(2)
      V2(3)= VTT(3)*V1(4)/VDM + VTT(4)*V1(3)/VDM
      V2(4)=-VTT(3)*V1(3)/VDM + VTT(4)*V1(4)/VDM
! Transform k to lab frame for first  decay particle
      CHK=V2(1)**2-V2(2)**2-V2(3)**2-V2(4)**2-M2**2
      IF(ABS(CHK).GT.0.005) WRITE(6,'(1X,''ERROR, CHK='',E9.2)') CHK
! Put k into vector for second decay particle
      KV(1)=SQRT(K*K+M3*M3)
      KV(2)=-K*SINTH*COS(PHI)
      KV(3)=-K*SINTH*SIN(PHI)
      KV(4)=-K*COSTH
! Transform k to lab frame for second decay particle
! First find vector as though decay were along Z axis
      V3(1)=G*(KV(1)+BETA*KV(4))
      VT(2)=KV(2)
      VT(3)=KV(3)
      VT(4)=G*(BETA*KV(1)+KV(4))
! Now rotate vector in xz plane
      VDM=SQRT(V1(2)**2+V1(4)**2)
      VTT(2)= VT(2)*V1(4)/VDM+VT(4)*V1(2)/VDM
      VTT(3)= VT(3)
      VTT(4)=-VT(2)*V1(2)/VDM+VT(4)*V1(4)/VDM
! Now rotate vector in yz plane
      VDM=SQRT(V1(3)**2+V1(4)**2)
      V3(2)= VTT(2)
      V3(3)= VTT(3)*V1(4)/VDM + VTT(4)*V1(3)/VDM
      V3(4)=-VTT(3)*V1(3)/VDM + VTT(4)*V1(4)/VDM
      CHK=V3(1)**2-V3(2)**2-V3(3)**2-V3(4)**2-M3**2
      IF(ABS(CHK).GT.0.005) WRITE(6,'(1X,''ERROR, CHK='',E9.2)') CHK
      CHK1=V1(1)-V2(1)-V3(1)
      CHK2=V1(2)-V2(2)-V3(2)+V1(3)-V2(3)-V3(3)+V1(4)-V2(4)-V3(4)
!      IF(ABS(CHK1).GT.0.005.OR.ABS(CHK2).GT.0.005) WRITE(6,
!     >  '(1X,''E12='',14F5.2)') CHK1,CHK2,V1,V2,V3
      RETURN
      END

      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL)   1717.   
      REAL*8 FUNCT,LOWER,UPPER,EPSLON                                   1718.   
         INTEGER NLVL                                                   1719.   
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  1720.   
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            1721.   
     1   FMRX(50), ESTRX(50), EPSX(50)                                  1722.   
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       1723.   
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS                         1724.   
         LEVEL = 0                                                      1725.   
         NLVL = 0                                                       1726.   
         ABAREA = 0.0                                                   1727.   
         L = LOWER                                                      1728.   
         R = UPPER                                                      1729.   
         FL = FUNCT(L)                                                  1730.   
         FM = FUNCT(0.5*(L+R))                                          1731.   
         FR = FUNCT(R)                                                  1732.   
         EST = 0.0                                                      1733.   
         EPS = EPSLON                                                   1734.   
  100 LEVEL = LEVEL+1                                                   1735.   
      M = 0.5*(L+R)                                                     1736.   
      COEF = R-L                                                        1737.   
      IF(COEF.NE.0) GO TO 150                                           1738.   
         ROMBRG = EST                                                   1739.   
         GO TO 300                                                      1740.   
  150 FML = FUNCT(0.5*(L+M))                                            1741.   
      FMR = FUNCT(0.5*(M+R))                                            1742.   
      ESTL = (FL+4.0*FML+FM)*COEF                                       1743.   
      ESTR = (FM+4.0*FMR+FR)*COEF                                       1744.   
      ESTINT = ESTL+ESTR                                                1745.   
      AREA=DABS(ESTL)+DABS(ESTR)                                        1746.   
      ABAREA=AREA+ABAREA-DABS(EST)                                      1747.   
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     1748.   
         NLVL = NLVL+1                                                  1749.   
         ROMBRG = ESTINT                                                1750.   
         GO TO 300                                                      1751.   
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                         1752.   
     1         (LEVEL.LT.MINLVL))  GO TO 400                            1753.   
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             1754.   
  300    LEVEL = LEVEL-1                                                1755.   
         I = RETRN(LEVEL)                                               1756.   
         VALINT(LEVEL, I) = ROMBRG                                      1757.   
         GO TO (500, 600), I                                            1758.   
  400    RETRN(LEVEL) = 1                                               1759.   
         MX(LEVEL) = M                                                  1760.   
         RX(LEVEL) = R                                                  1761.   
         FMX(LEVEL) = FM                                                1762.   
         FMRX(LEVEL) = FMR                                              1763.   
         FRX(LEVEL) = FR                                                1764.   
         ESTRX(LEVEL) = ESTR                                            1765.   
         EPSX(LEVEL) = EPS                                              1766.   
         EPS = EPS/1.4                                                  1767.   
         R = M                                                          1768.   
         FR = FM                                                        1769.   
         FM = FML                                                       1770.   
         EST = ESTL                                                     1771.   
         GO TO 100                                                      1772.   
  500    RETRN(LEVEL) = 2                                               1773.   
         L = MX(LEVEL)                                                  1774.   
         R = RX(LEVEL)                                                  1775.   
         FL = FMX(LEVEL)                                                1776.   
         FM = FMRX(LEVEL)                                               1777.   
         FR = FRX(LEVEL)                                                1778.   
         EST = ESTRX(LEVEL)                                             1779.   
         EPS = EPSX(LEVEL)                                              1780.   
         GO TO 100                                                      1781.   
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          1782.   
      IF(LEVEL.GT.1) GO TO 300                                          1783.   
      QUADMO = ROMBRG /12.0D0                                           1784.   
      RETURN                                                            1785.   
      END                                                               1786.   
      Subroutine WISER_ALL_SIG(E0,P,THETA_DEG,RAD_LEN,TYPE,SIGMA)

!------------------------------------------------------------------------------
! Calculate pi,K,p  cross section for electron beam on a proton target
! IntegrateQs over function WISER_FIT using integration routine QUADMO
! E0         is electron beam energy, OR max of Brem spectra
! P,E       is scattered particle  momentum,energy
! THETA_DEG  is kaon angle in degrees
! RAD_LEN (%)is the radiation length of target, including internal
!                (typically 5%)
!               = .5 *(target radiation length in %) +5.
!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA ***
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! SIGMA      is output cross section in nanobars/GeV-str
!------------------------------------------------------------------------------

      IMPLICIT NONE       
      REAL E0,P,THETA_DEG,RAD_LEN,SIGMA
      INTEGER TYPE
      COMMON/WISER_ALL/ E,P_COM,COST,P_T,TYPE_COM,PARTICLE ,M_X,U_MAN
      REAL E,P_COM,COST,P_T,M_X,U_MAN
      INTEGER TYPE_COM,PARTICLE
!  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
      REAL A5(6)/-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53/
      REAL A6(6)/-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45/
      REAL MASS2(3)/.019488, .2437, .8804/
      REAL MASS(3)/.1396, .4973, .9383/ 
      REAL MP/.9383/,  MP2/.8804/, RADDEG/.0174533/
      REAL  M_L,SIG_E
      REAL*8 E_GAMMA_MIN,WISER_ALL_FIT,QUADMO,E08,EPSILON/.003/
      EXTERNAL WISER_ALL_FIT                        
      INTEGER N,CHARGE
                                            
      P_COM = P
      TYPE_COM = TYPE
      PARTICLE = (TYPE+1)/2       ! 1= pi, 2= K, 3 =P
      CHARGE = TYPE -2*PARTICLE +2  ! 1 for + charge, 2 for - charge
      E08 =E0
                
      E =SQRT(MASS2(PARTICLE) + P**2)

      COST = COS(RADDEG * THETA_DEG)
      P_T = P * SIN(RADDEG * THETA_DEG)
      IF(TYPE.LE.4) THEN  !mesons
       IF(CHARGE.EQ.1) THEN   ! K+ n final state
        M_X = MP
       ELSE   ! K- K+ P final state
        M_X = MP+ MASS(PARTICLE)
       ENDIF
      ELSE  ! baryons 
       IF(CHARGE.EQ.1) THEN   ! pi p  final state
        M_X = MASS(1)  ! pion mass
       ELSE   ! P P-bar  P final state
        M_X = 2.*MP
       ENDIF
      ENDIF
      E_GAMMA_MIN = (M_X**2 -MASS2(PARTICLE ) -MP2+2.*MP*E)/
     >  (2.*(MP -E +P*COST))
!      WRITE(10,'(''E_GAMMA_MIN='',F10.2,''  p_t='',F8.2)')
!     >     E_GAMMA_MIN,P_T
!      E_GAMMA_MIN = MP *(E + MASS(PARTILCE))/(MP -P*(1.-COST))
      
      IF(E_GAMMA_MIN.GT..1) THEN !Kinematically allowed?
       M_L = SQRT(P_T**2 + MASS2(PARTICLE))    

       IF(TYPE.NE.5) THEN  ! everything but proton
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) *EXP(A6(TYPE) *P_T**2/E)
       ELSE ! proton

        U_MAN = ABS(MP2 + MASS2(PARTICLE) -2.*MP*E)
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) 
       ENDIF
       SIGMA = P**2/E * 1000. * RAD_LEN/100. *SIG_E 
      ELSE   ! Kinematically forbidden
       SIGMA = 0.
      ENDIF

      RETURN
      END

      REAL*8 FUNCTION WISER_ALL_FIT(E_GAMMA)

!---------------------------------------------------------
! Calculates  pi, k, p  cross section for gamma + p -> k
!  It is already divided by E_GAMMA, the bremstrulung spectra
! David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
! Can be called from WISER_SIG using integration routine QUADMO
! E,P are KAON energy and momentum
! P_t is KAON transverse momentum
! P_CM is KAON center of mass momentum
! P_CM_L is KAON center of mass longitudinal momentum
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! E_GAMMA is photon energy.
!             Steve Rock 2/21/96
!---------------------------------------------------------
                           
      IMPLICIT NONE       
      COMMON/WISER_ALL/ E,P,COST,P_T,TYPE,PARTICLE,M_X,U_MAN

      REAL  E,P,COST,P_T,M_X,U_MAN
      INTEGER  TYPE  !  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
      INTEGER PARTICLE   ! 1= pi, 2= K, 3 =P
!  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
      REAL A1(6)/566.,  486.,   368., 18.2,  1.33E5,  1.63E3 / 
      REAL A2(6)/829.,  115.,   1.91, 307.,  5.69E4, -4.30E3 / 
      REAL A3(6)/1.79,  1.77,   1.91, 0.98,  1.41,    1.79 / 
      REAL A4(6)/2.10,  2.18,   1.15, 1.83,   .72,    2.24 /
      REAL A6/1.90/,A7/-.0117/ !proton only
      REAL MASS2(3)/.019488, .2437, .8804/
      REAL MASS(3)/.1396, .4973, .9383/ 
      REAL MP2/.8804/,MP/.9383/, RADDEG/.0174533/
      REAL X_R, S,B_CM, GAM_CM,  P_CM
      REAL P_CM_MAX, P_CM_L
      REAL*8 E_GAMMA
                                            

!Mandlestam variables                                                
      S = MP2 + 2.* E_GAMMA * MP    

!Go to Center of Mass to get X_R
      B_CM = E_GAMMA/(E_GAMMA+MP)
      GAM_CM = 1./SQRT(1.-B_CM**2)
      P_CM_L = -GAM_CM *B_CM *E + 
     >          GAM_CM * P * COST
      P_CM = SQRT(P_CM_L**2 + P_T**2)  


      P_CM_MAX =SQRT (S +(M_X**2-MASS2(PARTICLE))**2/S 
     >    -2.*(M_X**2 +MASS2(PARTICLE)) )/2.
      X_R =  P_CM/P_CM_MAX   
       IF(X_R.GT.1.) THEN  ! Out of kinematic range
        WISER_ALL_FIT = 0.
       ELSEIF(TYPE.NE.5) THEN  ! not the proton
        WISER_ALL_FIT = (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)/E_GAMMA  
       ELSE ! special formula for proton
        WISER_ALL_FIT = ( (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)          /
     >   (1.+U_MAN)**(A6+A7*S) )/E_GAMMA  
       ENDIF
      
      RETURN
      END


       SUBROUTINE INEFT(QQ,W,W1,W2,amuM)                                      
                                                                        
C Modified 6feb87 by lww to accept target information passed through    
C common block /targt/.                                                 
                                                                        
C This program takes the old slac structure function model (Atwood,     
C Bodek, et.al.) and outputs values for W1 and W2 at given kinematics. 
! As of 11/3/95 this version is per NEUCLEON   ! Steve Rock

! amuM is atomic number, ie. 1. 2.xxx etc.

      Implicit None 
!      COMMON       /TARGT/ iZ, iA, avgN, avgA, avgM, amuM               
      REAL*8 QQ,W,W1,W2,amuM,WW,V,VV,OMEGAP,SP,UNIV,BRES,SLACF2,B
      REAL*8 VW2,X,EMCFAC
      REAL*8    C(24),CF(11),CD(24),CFD(11)                          
      REAL*8    EF(7) 
      REAL FITEMC
                                                  
      DATA         EF / -0.00136693,-.00510425,-.0375986,-.0946004,     
     +                  -.122435,-.0112751,0.406435/                    
      REAL*8         PM / .938256/,PMPM/.880324/,TPM/1.876512/            
      REAL*8         R /  .18/,ALPHAX/137.0388/,THCONST/0.0174533/
      LOGICAL GOODFIT
                                                                        
C FINAL HYDROGEN COEFFS. FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=3995,NO. OF POINTS=2533,FREE PARAMS=28,CHISQ/D.F.=1.59.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C C(24)=HYDROGEN COEFFICIENTS FOR B (BACKGROUND AND RESONANCE TERMS)    
                                                                        
      DATA   C(1) / 0.10741163E 01/,  C(2) / 0.75531124E 00/,           
     *       C(3) / 0.33506491E 01/,  C(4) / 0.17447015E 01/,           
     *       C(5) / 0.35102405E 01/,  C(6) / 0.10400040E 01/,           
     *       C(7) / 0.12299128E 01/,  C(8) / 0.10625394E 00/,           
     *       C(9) / 0.48132786E 00/,  C(10)/ 0.15101467E 01/,           
     *       C(11)/ 0.81661975E-01/,  C(12)/ 0.65587179E 00/,           
     *       C(13)/ 0.17176216E 01/,  C(14)/ 0.12551987E 00/,           
     *       C(15)/ 0.74733793E 00/,  C(16)/ 0.19538129E 01/,           
     *       C(17)/ 0.19891522E 00/,  C(18)/-0.17498537E 00/,           
     *       C(19)/ 0.96701919E-02/,  C(20)/-0.35256748E-01/,           
     *       C(21)/ 0.35185207E 01/,  C(22)/-0.59993696E 00/,           
     *       C(23)/ 0.47615828E 01/,  C(24)/ 0.41167589E 00/            
                                                                        
C CF(11)=HYDROGEN COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION) OMEGAW FIT   
                                                                        
      DATA CF(1) / 0.25615498E 00/,  CF(2) / 0.21784826E 01/,           
     *     CF(3) / 0.89783738E 00/,  CF(4) /-0.67162450E 01/,           
     *     CF(5) / 0.37557472E 01/,  CF(6) / 0.16421119E 01/,           
     *     CF(7) / 0.37635747E 00/,  CF(8) / 0.93825625E 00/,           
     *     CF(9) / 0.10000000E 01/,  CF(10)/ 0.0           /,           
     *     CF(11)/ 0.50000000E 00/                                      
                                                                        
C FINAL DEUTERIUM COEFFS FROM FITTING WITH NOTUSE=TRUE,ISLIDE=FALSE     
C CHISQ=4456,NO. OF POINTS 2303,FREE PERAMS=26,CHISQ/D.F.=1.96.         
C THIS DATA PROVIDED BY ARIE BODEK FOR E139 (10/3/83)                   
                                                                        
C CD(24)=DEUTERIUM COEFFICIENTS FOR B (BACKGROUND AND RESONANT TERMS)   
                                                                        
      DATA  CD(1) / 0.10521935E 01/, CD(2) / 0.76111537E 00/,           
     *      CD(3) / 0.41469897E 01/, CD(4) / 0.14218146E 01/,           
     *      CD(5) / 0.37119053E 01/, CD(6) / 0.74847487E 00/,           
     *      CD(7) / 0.12399742E 01/, CD(8) / 0.12114898E 00/,           
     *      CD(9) / 0.11497852E-01/, CD(10)/ 0.14772317E 01/,           
     *      CD(11)/ 0.69579815E-02/, CD(12)/ 0.12662466E 00/,           
     *      CD(13)/ 0.15233427E 01/, CD(14)/ 0.84094736E-01/,           
     *      CD(15)/ 0.74733793E 00/, CD(16)/ 0.19538129E 01/,           
     *      CD(17)/ 0.19891522E 00/, CD(18)/-0.24480414E 00/,           
     *      CD(19)/ 0.14502846E-01/, CD(20)/-0.35256748E-01/,           
     *      CD(21)/ 0.35185207E 01/, CD(22)/-0.21261862E 00/,           
     *      CD(23)/ 0.69690531E 01/, CD(24)/ 0.40314293E 00/            
                                                                        
C CFD(11) ARE DEUTERIUM COEFFICIENTS FOR F2 (UNIVERSAL FUNCTION)        
C OMEGAW FIT                                                            
                                                                        
      DATA CFD(1) / 0.47708776E 00/, CFD(2) / 0.21601918E 01/,          
     *     CFD(3) / 0.36273894E 01/, CFD(4) /-0.10470367E 02/,          
     *     CFD(5) / 0.49271691E 01/, CFD(6) / 0.15120763E 01/,          
     *     CFD(7) / 0.35114723E 00/, CFD(8) / 0.93825625E 00/,          
     *     CFD(9) / 0.10000000E 01/, CFD(10)/ 0.0           /,          
     *     CFD(11)/ 0.50000000E 00/                                     
                                                                        
C COMPUTE SOME KINEMATIC QUANTITIES                                     
                                                                        
      WW     = W**2                                                     
      V      = (WW+QQ-PMPM)/2.D0/PM                                     
      VV     = V*V                                                      
      OMEGAP = TPM*V/QQ+PMPM/QQ                                         
                                                                        
C OVERCOME RISK OF UNDERFLOW IN THE EXPONENTIATION                      
      OMEGAP = DMIN1(20.0D0,OMEGAP)                                       
                                                                        
      SP = 1.0-EXP(-7.7*(OMEGAP-1.0))                                   
      IF (amuM.LE.1.5) THEN !hydrogen
C          UNIVERSAL AND RESONANCE FIT FOR HYDROGEN                     
           UNIV = SLACF2(W,QQ,CF)                                       
           BRES = B(W,QQ,C)                                             
      ELSE                                                              
C          UNIVERSAL AND RESONANCE FIT FOR DEUTERIUM                    
           UNIV = SLACF2(W,QQ,CFD)/SP
           BRES = B(W,QQ,CD)          
      ENDIF                                                             
                                                                        
C COMPUTE VW2,W2,W1                                                     
                                                                        
      VW2    = UNIV*BRES 
      IF (amuM.GE.1.5) VW2=VW2/2.  !*****  per nucleon 11/3/95   ***********
      W2     = VW2/V                                                    
      W1     = (1.0D0+VV/QQ)/(V*(1.0D0+R))*VW2                          
      IF (amuM.LE.2.5) RETURN                                               
      X      = QQ/2./PM/V
      EMCFAC= FITEMC(REAL(X),REAL(amuM),GOODFIT)
C$$      SUMEF  = EF(1)                                                    
C$$      DO 11 J=2,7                                                       
C$$      ZZ     = J-1.                                                     
C$$11    SUMEF  = SUMEF+EF(J)*X**ZZ                                        
C$$      EMCFAC = 1+SUMEF*DLOG(amuM)                                       
                                                                        
      W2     = W2*EMCFAC                                                
      W1     = W1*EMCFAC                                                
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION SLACF2(WM,QSQ,CF)                                        
                                                                        
C UNIVERSAL FUNCTION FOR ATWOOD'S FIT                                   

      Implicit none                                      
      REAL*8    WM,QSQ,CF(11)                                               
      REAL*8    PM2/1.876512/, PMSQ/.880324/, PHTR/.61993/
      REAL*8    V,OMEGA,XX,XPX,OMEGAW,ARG

                                                                        
C OMEGAW FIT...NO PHOTO-PRODUCTION COUPLING                             
                                                                        
      V      = (WM**2+QSQ-PMSQ)/PM2                                     
      OMEGA  = 2.*CF(8)*V/QSQ                                           
      XX     = 1./OMEGA                                                 
      XPX    = CF(9)+CF(10)*(XX-CF(11))**2                              
      OMEGAW = (2.D0*CF(8)*V+CF(6))/(QSQ+CF(7))                         
      ARG    = 1.-1./OMEGAW                                             
                                                                        
      SLACF2 = OMEGAW/OMEGA*ARG**3*(CF(1)+CF(2)*ARG+                    
     >         CF(3)*ARG**2+CF(4)*ARG**3+CF(5)*ARG**4)                  
      SLACF2 = SLACF2*XPX                                               
                                                                        
      RETURN                                                            
      END                                                               
                                                                        
C-----------------------------------------------------------------------
                                                                        
      REAL*8 FUNCTION B(WM,QSQ,C)                                              
                                                                        
C BACKGROUND AND RESONANCE CONTRIBUTION FOR ATWOOD'S FIT                

      Implicit none
      REAL*8  WM,QSQ,C(24),WSQ,OMEGA,X,XPX,PIEMSQ,B1,EB1,B2,BBKG,BRES
      REAL*8  RAM,RMA,RWD,QSTARN,QSTARO,TERM,TERMO,GAMRES,BRWIG,RES
      REAL*8  RESSUM,EB2
      INTEGER   LSPIN(4),INDEX,J,K                                             
      DATA      LSPIN/1,2,3,2/                                       
      REAL*8    PMSQ/.880324/, PM2/1.876512/, PM/.938256/            
      INTEGER   NRES/4/, NBKG/5/,I                                     
                                                                        
C KINEMATICS                                                            
                                                                        
      WSQ    = WM**2                                                    
      OMEGA  = 1.+(WSQ-PMSQ)/QSQ                                        
      X      = 1./OMEGA                                                 
      XPX    = C(22)+C(23)*(X-C(24))**2                                 
      PIEMSQ = (C(1)-PM)**2                                             
                                                                        
C COLLECT BACKGROUND TERMS AND CHECK FOR EXPONENTIAL UNDERFLOWS BEFORE  
C THEY HAPPEN                                                           
                                                                        
      B1 = 0.                                                           
      IF (WM.GT.C(1)) B1 = C(2)                                         
      EB1 = C(3)*(WM-C(1))                                              
      IF (EB1.LE.25.) B1 = B1*(1.-EXP(-EB1))                            
      B2 = 0.                                                           
      IF (WM.GT.C(4)) B2 = (1.-C(2))                                    
      EB2 = C(5)*(WSQ-C(4)**2)                                          
      IF (EB2.LE.25.0) B2 = B2*(1.-EXP(-EB2))                           
      BBKG = B1+B2                                                      
      BRES = C(2)+B2                                                    
                                                                        
C COLLECT RES. CONTRIBUTION                                             
                                                                        
      RESSUM = 0.                                                       
      DO 30 I=1,NRES                                                    
           INDEX  = (I-1)*3+1+NBKG                                      
           RAM    = C(INDEX)                                            
           IF (I.EQ.1) RAM=C(INDEX)+C(18)*QSQ+C(19)*QSQ**2              
           RMA    = C(INDEX+1)                                          
           IF (I.EQ.3) RMA=RMA*(1.D0+C(20)/(1.D0+C(21)*QSQ))            
           RWD    = C(INDEX+2)                                          
           QSTARN =SQRT(DMAX1(0.D0,((WSQ+PMSQ-PIEMSQ)/(2.*WM))**2-PMSQ)) 
           QSTARO = SQRT(DMAX1(0.D0,((RMA**2-PMSQ+PIEMSQ)/                
     >              (2.*RMA))**2-PIEMSQ))                               
                                                                        
           RES = 0.                                                     
           IF (QSTARO.NE.0.) THEN                                       
                TERM   = 6.08974*QSTARN                                 
                TERMO  = 6.08974*QSTARO                                 
                J      = 2*LSPIN(I)                                     
                K      = J+1                                            
                GAMRES = RWD*(TERM/TERMO)**K*(1.+TERMO**J)/(1.+TERM**J) 
                GAMRES = GAMRES/2.                                      
                BRWIG  = GAMRES/((WM-RMA)**2+GAMRES**2)/3.1415926       
                RES    = RAM*BRWIG/PM2                                  
           ENDIF                                                        
           RESSUM = RESSUM+RES                                          
30    CONTINUE                                                          
                                                                        
C FORM VW2/F2                                                           
                                                                        
      B = BBKG*(1.+(1.-BBKG)*XPX)+RESSUM*(1.-BRES)                      
                                                                        
      RETURN                                                            
      END                                                               
      REAL FUNCTION FITEMC_N(X,A,Z,GOODFIT)                                         
!---------------------------------------------------------------------  
! Modified FITEMC.F with Neutron excess correction and proton=1 added 8/19/98
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! Z= number of protons
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
! 8/19/98 **  If proton, return 1.
! 8/19/98 **  Add neutron excess correction.
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL ALPHA, C,LN_C,X,A,Z ,X_U,SIG_N_P,F_IS
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points
!Term    Coeficient     Error
      REAL*8 ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965E-03,                                   
     >  2.18888887D+00,    3.792E-01,                                   
     > -2.46673765D+01,    6.302E+00,                                   
     >  1.45290967D+02,    4.763E+01,                                 
     > -4.97236711D+02,    1.920E+02,                                   
     >  1.01312929D+03,    4.401E+02,                                   
     > -1.20839250D+03,    5.753E+02,                                   
     >  7.75766802D+02,    3.991E+02,                                   
     > -2.05872410D+02,    1.140E+02 /                                  

             !
!Chisq=         22.    for 30 points
!Term    Coeficient     Error 
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        

      IF(A.LT.1.5) THEN    ! Added 8/19/98
       FITEMC_N=1.
       GOODFIT=.TRUE.
       RETURN
      ENDIF                                                                                
      IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
       IF(X.GT. 0.88) X_U =.88
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC_N  =  C *A**ALPHA    !isoscaler
      SIG_N_P = 1.-0.8*X_U
      F_IS = .5*(1.+ SIG_N_P)/(Z/A +(1.-Z/A)*SIG_N_P)
      FITEMC_N = FITEMC_N/F_IS
      RETURN                                                            
      END                                                               
      REAL FUNCTION FITEMC(X,A,GOODFIT)                                         
!---------------------------------------------------------------------  
! Fit to EMC effect.  Steve Rock 8/3/94                                 
! Funciton returns value of sigma(A)/sigma(d) for isoscaler nucleus     
!  with A/2 protons=neutrons.                                           
! A= atomic number                                                      
! x = Bjorken x.                                                        
!                                                                       
! Fit of sigma(A)/sigma(d) to form C*A**alpha where A is atomic number  
! First data at each x was fit to form C*A**alpha.  The point A=2(d)    
!  was included with a value of 1 and an error of about 2%.             
! For x>=.125 Javier Gomez fit of 7/93 to E139 data was used.           
! For .09 >=x>=.0085 NMC data from Amaudruz et al Z. Phys C. 51,387(91) 
!  Steve did the fit for alpha and C to the He/d. C/d and Ca/d NMC data.
! Alpha(x) was fit to a 9 term polynomial a0 +a1*x +a2*x**2 + a3*x**3 ..
! C(x) was fit to a 3 term polynomial in natural logs as                
!  Ln(C) = c0 + c1*Ln(x) + c2*[Ln(x)]**2.                               

! 6/2/98 *****  Bug (which set x= .00885 if x was out of range) fixed
!                    also gave value at x=.0085 if x>.88
!-----------------------------------------------------------------------
                                                                        
                                                                        
      IMPLICIT NONE                                                     
      INTEGER I                                                         
      REAL*4 ALPHA, C,LN_C,X,A ,X_U
      LOGICAL GOODFIT                                  
                                                                        
!Chisq=         19.   for 30 points                                     
!Term    Coeficient     Error                                           

      REAL*8  ALPHA_COEF(2,0:8)    /                                     
     > -6.98871401D-02,    6.965D-03,                                   
     >  2.18888887D+00,    3.792D-01,                                   
     > -2.46673765D+01,    6.302D+00,                                   
     >  1.45290967D+02,    4.763D+01,                                   
     > -4.97236711D+02,    1.920D+02,                                   
     >  1.01312929D+03,    4.401D+02,                                   
     > -1.20839250D+03,    5.753D+02,                                   
     >  7.75766802D+02,    3.991D+02,                                   
     > -2.05872410D+02,    1.140D+02 /                                  
                                                     
                              
!Chisq=         22.    for 30 points                                   
!Term    Coeficient     Error                                          
      REAL*8 C_COEF(2,0:2) /        ! Value and error for 6 term fit to 
     >  1.69029097D-02,    4.137D-03,                                   
     >  1.80889367D-02,    5.808D-03,                                   
     >  5.04268396D-03,    1.406D-03   /                                
                                                                        
                                                                        
      IF( (X.GT.0.88).OR.(X.LT. 0.0085) ) THEN   !Out of range of fit   
       IF(X.LT. 0.0085) X_U =.0085
       IF(X.GT. 0.88) X_U =.88
       GOODFIT=.FALSE.
      ELSE
       X_U=X
       GOODFIT=.TRUE.
      ENDIF                                                            
                                                                        
      LN_C = C_COEF(1,0)                                                
      DO I =1,2                                                         
       LN_C = LN_C + C_COEF(1,I) * (ALOG(X_U))**I                         
      ENDDO                                                             
                                                                        
      C = EXP(LN_C)                                                     
                                                                        
      ALPHA = ALPHA_COEF(1,0)                                           
      DO I=1,8                                                          
       ALPHA = ALPHA + ALPHA_COEF(1,I) * X_U**I                           
      ENDDO                                                             
                                                                        
      FITEMC  =  C *A**ALPHA                                            
      RETURN                                                            
      END

