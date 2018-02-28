C=======================================================================
C
C        INTEGRALJENJE MATRICA 3/D ELEMENATA
C
C   SUBROUTINE K21DMP
C              READM3
C              SIST3M
C              ELTM3
C
C=======================================================================
      SUBROUTINE K21DMP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C    GLAVNI PROGRAM ZA POZIVANJE PROGRAMA ZA RACUNANJE MATRICA ELEMENATA
C
      include 'paka.inc'
      
C
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' K21DMP'
C
      LAU=LMAX
      CALL READD3(A(LAU))
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      LAE=LMAX
      NCVE3=NCVE*3
      MXAE = NCVE3+(13*NCVE+NCVE3*(NCVE3+1)/2)*IDVA
      LMAX = LAE + MXAE
      IF(LMAX.LT.MTOT) GO TO 70
      WRITE(IZLAZ,2009) LMAX,MTOT
      STOP
C
C     FORMIRANJE MATRICE KRUTOSTI ELEMENATA I PAKOVANJE U SISTEM
C
   70 CALL SIST3D(A(LAE),A(LAU))
C
      RETURN
C-----------------------------------------------------------------------
 2009 FORMAT(///' NEDOVOLJNA DIMENZIJA U VEKTORU A ZA MATRICE ELEMENATA'
     1/' POTREBNA DIMENZIJA, LMAX=',I10/
     2' RASPOLOZIVA DIMENZIJA, MTOT=',I10)
C-----------------------------------------------------------------------
      END
C======================================================================
      SUBROUTINE READD3(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
C
      include 'paka.inc'
      
C
      COMMON /IZOL4B/ NGS12,ND,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /SIMO2D/ LALFE,LHAEM,LHINV,LGEEK,IALFA
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /INCONF/ LINDBEL,LBIRTHC
      COMMON /INTGRA/ INCOTX,INCOTY,INCOTZ
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /GAUSVR/ LTEMGT,LCORGT,ICORGT
      COMMON /CVSILE/ NSILA,LESILA
      COMMON /COEFSM/ COEF(3),ICOEF
      COMMON /IKOVAR/ INDKOV
      COMMON /CEPMAT/ INDCEP
      COMMON /LEVDES/ ILEDE,NLD,ICPM1
      COMMON /CDEBUG/ IDEBUG
C
      DIMENSION AU(*)
      REAL AU
C
C     POZIVANJE PROGRAMA ZA ULAZNE PODATKE .
C
      IF(IDEBUG.GT.0) PRINT *, ' READD3'
C
      LSTAZA(1)=LMAX8
      READ(IELEM,REC=LMAX8)
     1NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,ICORGT,LCEL,LELC,NMA,NMI,
     1MXAU,LNEL,LNMAT,LIPGC,LIPRC,LISNA,LTHID,LLMEL,IPODT,LNNOD,
     1ND,NGS12,MSLOJ,MXS,MSET,LNSLOJ,LMATSL,LDSLOJ,LBBET,NSILA,LESILA,
     1INCOTX,INCOTY,INCOTZ,LBET0,((CPP(J,I),J=1,3),I=1,3),
     1((TSG(J,I),J=1,6),I=1,6),BETA,IBB0,LALFE,LHAEM,LHINV,LGEEK,IALFA,
     1INDBTH,INDDTH,LTBTH,LTDTH,INDKOV,INDCEP,ILEDE,NLD,ICPM1,
     1COEF(3),ICOEF,LINDBEL,LBIRTHC
      LSTAZA(2)=LMAX8+1
      CALL READDD(AU,MXAU/IDVA,IELEM,LMAX8,LDUZI)
      LMAX=LAU+MXAU
      CALL DELJIV(LMAX,2,INDL)
      IF(INDL.EQ.0) LMAX=LMAX+1
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE SIST3D(AE,AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     GLAVNI UPRAVLJACKI PROGRAM  ZA MATRICE ELEMENATA I SISTEMA
C     RACUNANJE NAPONA
C
      include 'paka.inc'
      
C
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BTHDTH/ INDBTH,INDDTH,LTBTH,LTDTH
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /CDEBUG/ IDEBUG
      DIMENSION AE(*),AU(*)
      REAL AE,AU
C
C     REPERI U VEKTORU ELEMENATA AE
C
      IF(IDEBUG.GT.0) PRINT *, ' SIST3D'
C
      NCVE3=NCVE*3
      NWE=NCVE3*(NCVE3+1)/2
      LHE=1
      LBET=LHE+4*NCVE*IDVA
      LSKE=LBET+3*NCVE3*IDVA
      LLM=LSKE+NWE*IDVA
      MXAE1=LLM+NCVE3-1
      IF(MXAE1.NE.MXAE) THEN
        WRITE(IZLAZ,2000)MXAE1,MXAE
 2000   FORMAT(//' NE POKLAPAJU SE MXAE1 I MXAE'/' PROGRAM STOP'/2I10//)
      STOP ' PROGRAM STOP - PAK35 - SIST3D'
      ENDIF
      IZBR=13*NCVE+NWE
C
C     OSNOVNA PETLJA PO ELEMENTIMA
C
      DO 100 NLM=1,NE
C
CS       NASTAJANJE I NESTAJANJE ELEMENATA
CE       ELEMENT BIRTH AND DEATH OPTION
C
         IBD=0
         CALL DTHBTH(AU(LTBTH),AU(LTDTH),VREME,NLM,IBD)
         IF(IBD.EQ.1) GO TO 100
C
      CALL CLEAR(AE,IZBR)
C
C     RACUNANJE MATRICE ELEMENATA I/ILI NAPONA
C
      KORD=LCORD
      IF(IATYP.EQ.3) KORD=LCORUL
      CALL ELTD3(AE(LSKE),AE(LLM),AU(LNEL),AU(LNMAT),AU(LTHID),
     1AE(LHE),AE(LBET),A(KORD),AU(LIPGC),AU(LLMEL),A(LGUSM),NCVE3,
     1AU(LCEL))
  100 CONTINUE
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE ELTD3(SKE,LM,NEL,NMAT,TSGE,
     1                 HE,BET,CORD,IPGC,LMEL,GUSM,NCVE3,MCVEL)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     INTEGRACIJA MATRICE MASA 3D ELEMENATA
C
      include 'paka.inc'
      
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /REPERM/ MREPER(4)
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /CDEBUG/ IDEBUG
      COMMON /NEWMARK/ ALFAM,BETAK,DAMPC,NEWACC
      COMMON /ELEM3D/ DC(3,3)
      COMMON /INDABC/ INDVB
C
      DIMENSION SKE(*),LM(*),NEL(NE,*),NMAT(*),TSGE(6,6,*),
     1CORD(NP,*),HE(NCVE,*),BET(3,*),IPGC(*),
     1LMEL(NCVE3,*),GUSM(50,*),CMC(21),MCVEL(*)
C
      DIMENSION XG(55),WGT(55),NREF(11),XGG(15)
C
      DATA NREF/0,1,3,6,10,15,21,28,36,45,55/
      DATA WGT/            2.D0,               1.D0,               1.D0,
     1       .555555555555556D0, .888888888888889D0, .555555555555556D0,
     2       .347854845137454D0, .652145154862546D0, .652145154862546D0,
     3       .347854845137454D0, .236926885056189D0, .478628670499366D0,
     4       .568888888888889D0, .478628670499366D0, .236926885056189D0,
     5       .171324492379170D0, .360761573048139D0, .467913934572691D0,
     6       .467913934572691D0, .360761573048139D0, .171324492379170D0,
     7       .129484966168870D0, .279705391489277D0, .381830050505119D0,
     8       .417959183673469D0, .381830050505119D0, .279705391489277D0,
     9       .129484966168870D0, .101228536290376D0, .222381034453374D0,
     9       .313706645877887D0, .362683783378362D0, .362683783378362D0,
     1       .313706645877887D0, .222381034453374D0, .101228536290376D0,
     2       .081274388361574D0, .180648160694857D0, .260610696402935D0,
     3       .312347077040003D0, .330239355001260D0, .312347077040003D0,
     4       .260610696402935D0, .180648160694857D0, .081274388361574D0,
     5       .066671344308688D0, .149451349150581D0, .219086362515982D0,
     6       .269266719309996D0, .295524224714753D0, .295524224714753D0,
     7       .269266719309996D0, .219086362515982D0, .149451349150581D0,
     8       .066671344308688D0/
      DATA XG /            0.D0,-.577350269189626D0, .577350269189626D0,
     1      -.774596669241483D0,               0.D0, .774596669241483D0,
     2      -.861136311594053D0,-.339981043584856D0, .339981043584856D0,
     3       .861136311594053D0,-.906179845938664D0,-.538469310105683D0,
     4                     0.D0, .538469310105683D0, .906179845938664D0,
     5      -.932469514203152D0,-.661209386466265D0,-.238619186083197D0,
     6       .238619186083197D0, .661209386466265D0, .932469514203152D0,
     7      -.949107912342759D0,-.741531185599394D0,-.405845151377397D0,
     8                     0.D0, .405845151377397D0, .741531185599394D0,
     9       .949107912342759D0,-.960289856497536D0,-.796666477413627D0,
     9      -.525532409916329D0,-.183434642495650D0, .183434642495650D0,
     1       .525532409916329D0, .796666477413627D0, .960289856497536D0,
     2      -.968160239507626D0,-.836031107326636D0,-.613371432700590D0,
     3      -.324253423403809D0,               0.D0, .324253423403809D0,
     4       .613371432700590D0, .836031107326636D0, .968160239507626D0,
     5      -.973906528517172D0,-.865063366688985D0,-.679409568299024D0,
     6      -.433395394129247D0,-.148874338981631D0, .148874338981631D0,
     7       .433395394129247D0, .679409568299024D0, .865063366688985D0,
     8       .973906528517172D0/
C
      DATA XGG/ 0.000000000000000,-1.000000000000000, 1.000000000000000,
     1         -1.000000000000000, 0.000000000000000, 1.000000000000000,
     2         -1.000000000000000,-0.333333333333333, 0.333333333333333,
     3          1.000000000000000,-1.000000000000000,-0.500000000000000,
     4          0.000000000000000, 0.500000000000000, 1.000000000000000/
C
C     FORMIRANJE VEKTORA LM
C
      IF(IDEBUG.GT.0) PRINT *, ' ELTD3'
C     prebaceno na pocetak celine
c      IF(IDAMP.EQ.3) RETURN
C
      DO 10 NC=1,NCVE3
      LM(NC)=LMEL(NC,NLM)
   10 CONTINUE
C
C     PETLJA PO GAUSOVIM TACKAMA
C
C      CALL WRR(RTDT,JEDN,'RTDT')
      IPGCT=IPGC(NLM)
      MAT=NMAT(NLM)
      GUST=GUSM(NMODM,MAT)
      AMASA=0.D0
c      GO TO 999
C
C     Potrebni su podaci o materijalu za izracunavanje a i b i E, ni i ro
C     broj elementa i broj povrsine na kojoj je granica
C
C     konstatnno R

C
         GO TO (  1,  2,999, 40,999,999,999,999,999, 40,
     1          999,999,999,999,999,999, 40, 40, 40,999,
     2           40,999,999,999,999,999,999,999,999,999,
     3          999,999,999,999,999,999,999,999,999,999,
     4          999,999,999,999,999,999,999,999,999,999,
     5          999,999,999,999,999,999,999,999,999,999,
     6          999,999,999,999,999,999,999,999,999,999,
     7          999,999,999,999,999,999,999,999,999,999,
     8          999,999,999,999,999,999,999,999,999,999,
     9          999,999,999,999,999,999,999,999,999,999),NMODM
C
CE       ELAST(6,6): ELASTIC COEFFICIENT MATRIX IN CONSTITUTIVE EQUATION
CE                   FOR ISOTROPIC CASE
CE       LFUN: POINTER FOR MATERIAL DATA
    1    LFUN=MREPER(1)
         CALL MDV31(A(LFUN),GUST)
         GO TO 999
C
CE       ELAST(6,6): ELASTIC COEFFICIENT MATRIX IN CONSTITUTIVE EQUATION
CE                   FOR ANISOTROPIC CASE
CE       LFUN: POINTER FOR MATERIAL DATA
    2    LFUN=MREPER(1)
CE       COMPUTE ELASTIC COEFFICIENT MATIX USING ELASTIC MODULUS AND 
CE       POISSON'S RATIO IN X,Y,Z DIRECTIONS
         CALL MEL32(A(LFUN))
C
CE       TAKE THE ORIENTATION OF MATERIAL INTO ACCOUNT
C
CE       IBB0: INDICATOR FOR FIBER DIRECTION
CE       (=0-IN LOCAL R,S PLANE (BETA); =1-IN SPACE [CPP(3,3)],
CE        SEE CARDS /13-3/A AND /13-3/A1 IN USER MANUAL)
CE       TBETA(6,6): CONSTITUTIVE TRANSFORMATION MATRIX FOR ANGLE (BETA)
CE       BETA: ANGLE BETWEEN FIBER AND LOCAL AXIS R IN R,S PLANE
CE       TSGE(6,6,NE): CONSTITUTIVE TRANSFORMATION MATRIX 
CE                 [LOCAL CARTESIAN - MATERIAL (FIBER) DIRECTIONS]
CE       TSS(6,6): CONSTITUTIVE TRANSFORMATION MATRIX 
CE                 [GLOBAL CARTESIAN - LOCAL CARTESIAN]
CE       TSG(6,6): CONSTITUTIVE TRANSFORMATION MATRIX FOR (NMODM)=2
CE                 [GLOBAL CARTESIAN - MATERIAL (FIBER) DIRECTIONS]
CE       XJJ(3,3): JACOBIAN MATRIX
C        (Cg=QeT*Cl*Qe)
         IF(IBB0.EQ.1) THEN
            IF(IATYP.GT.1) THEN
               CALL JACTE3(NOP,CORD,HE,R,S,T,0)
               CALL TRANAL(XJJ,TSS,0)
               CALL MNOZM1(TSG,TSGE(1,1,NLM),TSS,6,6,6)
            ENDIF
         ELSE
            IF(IATYP.LE.1) THEN
               CALL JEDNA1(TSG,TSGE(1,1,NLM),36)
            ELSE
               CALL JACTE3(NOP,CORD,HE,R,S,T,0)
               CALL TRANAL(XJJ,TSS,0)
               CALL MNOZM1(TSG,TBETA,TSS,6,6,6)
            ENDIF
         ENDIF
CE       TRANSFORM  CONSTITUTIVE MATRIX (Cl-LOCAL) TO (Cg-GLOBAL)
C        (Cg=QeT*Cl*Qe)
         CALL TRAETP(ELAST,ELAST,TSG)
         GO TO 999
C
CE       TSG(6,6): CONSTITUTIVE TRANSFORMATION MATRIX FOR (NMODM)=4
   40    IF(IBB0.EQ.1) THEN
            IF(IATYP.GT.1) THEN
               CALL JACTE3(NOP,CORD,HE,R,S,T,0)
               CALL TRANAL(XJJ,TSS,0)
               CALL MNOZM1(TSG,TSGE(1,1,NLM),TSS,6,6,6)
            ENDIF
         ELSE
            IF(IATYP.LE.1) THEN
               CALL JEDNA1(TSG,TSGE(1,1,NLM),36)
            ELSE
               CALL JACTE3(NOP,CORD,HE,R,S,T,0)
               CALL TRANAL(XJJ,TSS,0)
               CALL MNOZM1(TSG,TBETA,TSS,6,6,6)
            ENDIF
         ENDIF
C INTEGRACIJA MATRICE D ZA UCITANO b
  999 DO 500 NGR=1,NGAUSX
      JGR=NREF(NGAUSX)+NGR
      R=XG(JGR)
      IF(IPGCT.EQ.1) R=XGG(JGR)
      WR=WGT(JGR)
C
      DO 50 NGS=1,NGAUSY
      JGS=NREF(NGAUSY)+NGS
      S=XG(JGS)
      IF(IPGCT.EQ.1) S=XGG(JGS)
      WS=WGT(JGS)
C
      DO 55 NGT=1,NGAUSZ
      JGT=NREF(NGAUSZ)+NGT
      T=XG(JGT)
      IF(IPGCT.EQ.1) T=XGG(JGT)
      WT=WGT(JGT)
CS     JAKOBIJAN
      CALL JACTE3(NEL,CORD,HE,R,S,T,0)
      WD=WR*WS*WT*DET*DAMPC
      AMASA= AMASA+WD
C
C     INTEGRACIJA KONZISTENTNE MATRICE PRIGUSENJA
C
      IF(IDAMP.EQ.1)THEN
        CALL BETVDY(BET,HE,NCVE,3)
        CALL INTEGV(SKE,BET,LM,WD,NCVE3,3)
      ENDIF
C
   55 CONTINUE
   50 CONTINUE
  500 CONTINUE
C
C     VISCOUS BOUNDARY CONDITIONS
C
      NMM=NLM
      IF(ICVEL.EQ.1) NMM=MCVEL(NLM)
C     element na granici u kome se prigusuje talas
      IF(NMM.EQ.INDVB) THEN
      DO 530 NGR=1,NGAUSX
      JGR=NREF(NGAUSX)+NGR
      R=XG(JGR)
      IF(IPGCT.EQ.1) R=XGG(JGR)
      WR=WGT(JGR)
C
      DO 530 NGS=1,NGAUSY
      JGS=NREF(NGAUSY)+NGS
      S=XG(JGS)
      IF(IPGCT.EQ.1) S=XGG(JGS)
      WS=WGT(JGS)
C
      DO 530 NGT=1,NGAUSZ
      JGT=NREF(NGAUSZ)+NGT
      T=XG(JGT)
      IF(IPGCT.EQ.1) T=XGG(JGT)
      WT=WGT(JGT)
CS     JAKOBIJAN
       CALL JACTE3(NEL,CORD,HE,R,S,T,0)
       WD=WR*WS*WT*DET
C
C     INTEGRACIJA KONZISTENTNE MATRICE VISKOZNOG PRIGUSENJA (6.10)
C
        CALL BETVDY(BET,HE,NCVE,3)
        CALL INTEGK(SKE,BET,DC,LM,WD,NCVE3,3)
C
  530 CONTINUE
      ENDIF
C
C     KONCETRISANA MATRICA PRIGUSENJA I PAKOVANJE
C
      IF(IDAMP.EQ.2)THEN
        CALL PODMA3(NEL,CMC)
        CALL LUMMAS(ALSK,A(LMAXA),NEL,LM,CMC,AMASA,NCVE,NE,NLM,3)
      ENDIF
C
C     PAKOVANJE KONZISTENTNE MATRICE ELEMENTA U MATRICU SISTEMA
C
      IF(IDAMP.EQ.1) CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,NCVE3)
C      
      RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE MDV31(FUN,RO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       TO FORM ELASTIC MATRIX FOR MATERIAL MODEL 1, SEE /11/
C .
C ......................................................................
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /MATIZO/ E,V
      COMMON /ELEM3D/ DC(3,3)
C
      DIMENSION FUN(2,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' MEL31 '
C
      E=FUN(1,MAT)
      V=FUN(2,MAT)
      G=E/2./(1.+V)
      S=DSQRT((1.-2.*V)/2./(1.-V))
      VS=DSQRT(G/RO)
      VP=VS/S
      CALL CLEAR(DC,9)
      A=1.
      B=1.
      DC(1,1)=RO*A*VS
      DC(2,2)=RO*A*VS
      DC(3,3)=RO*B*VP
      RETURN
      END
C======================================================================
      SUBROUTINE JACTE3V(NOP,CORD,H,R,S,T,KFIX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.   P R O G R A M
CE.       FOR INTERPOLATION FUNCTIONS AND JACOBIAN MATRIX IN CURRENT
CE.       INTEGRATION POINT (R,S,T - ARE NATURAL COORDINATES)
C .
C ......................................................................
C
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NOP(NE,*),CORD(NP,*),H(NCVE,*),IPERM(8),NCVOR(8),
     1          HV(21,4)
      DATA IPERM/2,3,4,1,6,7,8,5/
C
      IF(IDEBUG.GT.0) PRINT *, ' JACTE3'
      RP=1.0+R
      SP=1.0+S
      TP=1.0+T
      RM=1.0-R
      TM=1.0-T
      SM=1.0-S
      RR=1.0-R*R
      SS=1.0-S*S
      TT=1.0-T*T
C
C     INTERPOLACIJSKE FINKCIJE I NJIHOVI IZVODI
C
C
C     PRVIH 8 CVOROVA
C
      H(1,1)=0.125*RP*SP*TP
      H(2,1)=0.125*RM*SP*TP
      H(3,1)=0.125*RM*SM*TP
      H(4,1)=0.125*RP*SM*TP
      H(5,1)=0.125*RP*SP*TM
      H(6,1)=0.125*RM*SP*TM
      H(7,1)=0.125*RM*SM*TM
      H(8,1)=0.125*RP*SM*TM
C
      H(1,2)=0.125*SP*TP
      H(2,2)=-H(1,2)
      H(3,2)=-0.125*SM*TP
      H(4,2)=-H(3,2)
      H(5,2)=0.125*SP*TM
      H(6,2)=-H(5,2)
      H(7,2)=-0.125*SM*TM
      H(8,2)=-H(7,2)
C
      H(1,3)=0.125*RP*TP
      H(2,3)=0.125*RM*TP
      H(3,3)=-H(2,3)
      H(4,3)=-H(1,3)
      H(5,3)=0.125*RP*TM
      H(6,3)=0.125*RM*TM
      H(7,3)=-H(6,3)
      H(8,3)=-H(5,3)
C
      H(1,4)=0.125*RP*SP
      H(2,4)=0.125*RM*SP
      H(3,4)=0.125*RM*SM
      H(4,4)=0.125*RP*SM
      H(5,4)=-H(1,4)
      H(6,4)=-H(2,4)
      H(7,4)=-H(3,4)
      H(8,4)=-H(4,4)
C
      IF(NCVE.EQ.8) GO TO 50
      NND9=NCVE-8
      I=0
    5 I=I+1
      IF(I.GT.NND9) GO TO 80
      NN=I+8
      IF(NOP(NLM,NN).EQ.0) GO TO 5
      GO TO (9,10,11,12,13,14,15,16,17,18,19,20,21),I
C
C     STEPENE SLOBODE ZA CVOROVE PREKO 8
C
C     DEVETI CVOR
    9 H(9,1)=0.25*RR*SP*TP
      H(9,2)=-0.50*R*SP*TP
      H(9,3)=0.25*RR*TP
      H(9,4)=0.25*RR*SP
      GO TO 5
C     DESETI CVOR
   10 H(10,1)=0.25*RM*SS*TP
      H(10,2)=-0.25*SS*TP
      H(10,3)=-0.50*RM*S*TP
      H(10,4)=0.25*RM*SS
      GO TO 5
C     JEDANAESTI CVOR
   11 H(11,1)=0.25*RR*SM*TP
      H(11,2)=-0.50*R*SM*TP
      H(11,3)=-0.25*RR*TP
      H(11,4)=0.25*RR*SM
      GO TO 5
C     DVANAESTI CVOR
   12 H(12,1)=0.25*RP*SS*TP
      H(12,2)=0.25*SS*TP
      H(12,3)=-0.50*RP*S*TP
      H(12,4)=0.25*RP*SS
      GO TO 5
C     TRINAESTI CVOR
   13 H(13,1)=0.25*RR*SP*TM
      H(13,2)=-0.50*R*SP*TM
      H(13,3)=0.25*RR*TM
      H(13,4)=-0.25*RR*SP
      GO TO 5
C     CETRNAESTI CVOR
   14 H(14,1)=0.25*RM*SS*TM
      H(14,2)=-0.25*SS*TM
      H(14,3)=-0.50*RM*S*TM
      H(14,4)=-0.25*RM*SS
      GO TO 5
C     PETNAESTI CVOR
   15 H(15,1)=0.25*RR*SM*TM
      H(15,2)=-0.50*R*SM*TM
      H(15,3)=-0.25*RR*TM
      H(15,4)=-0.25*RR*SM
      GO TO 5
C     SESNAESTI CVOR
   16 H(16,1)=0.25*RP*SS*TM
      H(16,2)=0.25*SS*TM
      H(16,3)=-0.50*RP*S*TM
      H(16,4)=-0.25*RP*SS
      GO TO 5
C     SEDAMNAESTI CVOR
   17 H(17,1)=0.25*RP*SP*TT
      H(17,2)=0.25*SP*TT
      H(17,3)=0.25*RP*TT
      H(17,4)=-0.50*RP*SP*T
      GO TO 5
C     OSAMNAESTI CVOR
   18 H(18,1)=0.25*RM*SP*TT
      H(18,2)=-0.25*SP*TT
      H(18,3)=0.25*RM*TT
      H(18,4)=-0.50*RM*SP*T
      GO TO 5
C     DEVETNAESTI CVOR
   19 H(19,1)=0.25*RM*SM*TT
      H(19,2)=-0.25*SM*TT
      H(19,3)=-0.25*RM*TT
      H(19,4)=-0.50*RM*SM*T
      GO TO 5
C     DVADESETI CVOR
   20 H(20,1)=0.25*RP*SM*TT
      H(20,2)=0.25*SM*TT
      H(20,3)=-0.25*RP*TT
      H(20,4)=-0.50*RP*SM*T
      GO TO 5
C     DVADESETPRVI CVOR
   21 H(21,1)=RR*SS*TT
      H(21,2)=-2.0*R*SS*TT
      H(21,3)=-2.0*S*RR*TT
      H(21,4)=-2.0*T*RR*SS
C
C     KOREKCIJE PRVIH 20 FINKCIJA AKO JE UPOTREBLJEN CVOR 21
C
      DO 150 I=1,8
      DO 150 J=1,4
  150 H(I,J)=H(I,J)-0.125*H(21,J)
      DO 160 I=9,20
      IF(NOP(NLM,I).EQ.0) GO TO 160
      DO 155 J=1,4
  155 H(I,J)=H(I,J)-0.25*H(21,J)
  160 CONTINUE
C
C     KOREKCIJE PRVIH 8 FUNKCIJA AKO SU UPOTREBLJENI CVOROVI PREKO 8
C
C  MEDJUCVOROVI OD 9 DO 16
   80 I16=16
      IF(NND9.LT.8) I16=NCVE
      DO 200 I=9,I16
      IF(NOP(NLM,I).EQ.0) GO TO 200
      DO 210 J=1,4
      I1=I-8
      I2=IPERM(I1)
      H(I1,J)=H(I1,J)-0.5*H(I,J)
  210 H(I2,J)=H(I2,J)-0.5*H(I,J)
  200 CONTINUE
C  MEDJUCVOROVI OD 17 DO 20
      IF(NCVE.LE.16) GO TO 50
      I20=20
      IF(NND9.LT.20) I20=NCVE
      DO 250 I=17,I20
      IF(NOP(NLM,I).EQ.0) GO TO 250
      DO 260 J=1,4
      I1=I-16
      I2=I1+4
      H(I1,J)=H(I1,J)-0.5*H(I,J)
  260 H(I2,J)=H(I2,J)-0.5*H(I,J)
  250 CONTINUE
C     HV FUNKCIJE
      CALL CLEAR(HV,84)
      HV(5,2)=H(5,2)*4./TM/TM
      HV(6,2)=H(6,2)*4./TM/TM
      HV(7,2)=H(7,2)*4./TM/TM
      HV(8,2)=H(8,2)*4./TM/TM
      HV(5,3)=H(5,3)*4./TM/TM
      HV(6,3)=H(6,3)*4./TM/TM
      HV(7,3)=H(7,3)*4./TM/TM
      HV(8,3)=H(8,3)*4./TM/TM
      HV(5,4)=H(5,4)*4./TM/TM+H(5,1)*8./TM/TM/TM
      HV(6,4)=H(6,4)*4./TM/TM+H(6,1)*8./TM/TM/TM
      HV(7,4)=H(7,4)*4./TM/TM+H(7,1)*8./TM/TM/TM
      HV(8,4)=H(8,4)*4./TM/TM+H(8,1)*8./TM/TM/TM
      IF(NCVE.GT.8) THEN
        DO J=13,NCVE
         HV(J,2)=H(J,2)/TM/TM
         HV(J,3)=H(J,3)/TM/TM
         HV(J,4)=H(J,4)/TM/TM+H(J,1)*2./TM/TM/TM
        ENDDO
      ENDIF
C
C     JAKOBIJAN U TACKI R,S,T
C
   50 DO 60 I=1,3
      DO 60 J=1,3
      XJ(I,J)=0.
      DO 60 KM=1,NCVE
      K=NOP(NLM,KM)
      IF(K.EQ.0) GO TO 60
      XJ(I,J)=XJ(I,J)+HV(KM,I+1)*CORD(K,J)
   60 CONTINUE
      CALL JEDNA1(XJJ,XJ,9) 
      IF(KFIX.GT.0) GO TO 70
C
C     DERERMINANTA JAKOBIJANA U TACKI R,S,T
C
      CALL WRR3(XJ,9,'XJV ')
      CALL MINV3(XJ,DET)
CZ    DET=DABS(DET)
      IF(DET.GT.1.D-13) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NLM,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NLM,KFIX,R,S,T,DET
	write(izlaz,*) 'ncve',ncve
C      IF(NCVE.GT.8) NCVE=8
      I0=0
      I1=1
      I7=7
      I8=8
      I19=19
      I115=115
      DO 61 KM=1,NCVE
      K=NOP(NLM,KM)
      IF(K.EQ.0) GO TO 61
      NCVOR(KM)=K
      WRITE(IZLAZ,1000) K,I0,I0,I8,(CORD(K,J),J=1,3)
C      write(izlaz,1200) (h(km,j),j=1,4)
   61 CONTINUE
      WRITE(IZLAZ,1100) NLM,I19,I115,I1,I1,I7,NCVE
      WRITE(IZLAZ,1100) (NCVOR(J),J=1,NCVE)

      STOP 'PROGRAM STOP - PAK32 - JACTE3'
C
C     DETERMINATA POVRSINSKOG JAKOBIJANA
   70 GO TO (71,72,73),KFIX
C     KONSTANTNO KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      GO TO 74
C     KONSTANTNO ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      GO TO 74
C     KONSTANTNO ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
   74 DET=DSQRT(DET)
      IF(DET.GT.1.D-13) RETURN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000) NLM,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000) NLM,KFIX,R,S,T,DET
      STOP 'PROGRAM STOP - PAK32 - JACTE3'
C
 1000 FORMAT(4I10,3(1PE13.5))
 1100 FORMAT(8I10)
 1200 FORMAT(' h',4(1PE13.5))
C-----------------------------------------------------------------------
 2000 FORMAT(' ** GRESKA **: NEGATIVNA ILI NULA DETERMINATA JAKOBIJANA',
     1       ' ZA ELEMENT BR.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F12.5)
C-----------------------------------------------------------------------
 6000 FORMAT(/' ','ZERO OR NEGATIVE JACOBIAN DETERMINANTE'/
     1' ','ELEMENT NUM. =',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',D12.5)
C-----------------------------------------------------------------------
      END
