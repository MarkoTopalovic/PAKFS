C* SILE93 ----> RETURN
C=======================================================================
C
C   PENALTI METOD
C
C   SUBROUTINE LMMH93P
C              SIS93P
C              ELTE93P
C              SILE93P
C			 ALFC3P
C=======================================================================
      SUBROUTINE LMMH93P(ID,NEL,NCVSF,ITSRF,NELSF,ISNA)
C
C	kod Lagranza ovaj podprogram je u pak931
CS     FORMIRANJE VEKTORA LM I VISINA STUBOVA (MHT)
CE     OBTAIN VECTOR  LM   AND COLUMN HEIGHTS
C

      include 'paka.inc'
       
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
C
      DIMENSION ID(NP,*),NEL(NE,*),NCVSF(*),ITSRF(*),NELSF(*),
     &          LMEL(30),ISNA(*)
      IF(ISRPS.EQ.0.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,2025)
      IF(ISRPS.EQ.1.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,6025)
C  FORMIRANJE DOPUNSKIH JEDNACINA ZA LAGRANZEV MNOZIOC
C  PETLJA PO KONTAKTNIM ELEMENTIMA (CVOROVIMA KONTAKTORA)
      CALL ICLEAR(LMEL,30)
C      ND=NCVE*3
      DO 100 NLM=1,NE
      ISRF=NEL(NLM,2)
      MAXCE=ISNA(ISRF)
C
        IF(MAXCE.EQ.2) THEN
           K2=2
         ELSE
           K2=3
         ENDIF
C
       JJ=NEL(NLM,1)
       DO 5 I=1,K2
C       IF(ID(JJ,I).GT.0)THEN
C         JEDN=JEDN+1
C         IDUM=JEDN
C       ELSE
C         IDUM=0
C       ENDIF
C       IDC(NLM,I)=IDUM
c    5  LMEL(I)=IDC(NLM,I)
    5   LMEL(I)=ID(JJ,I)
C    5  WRITE(3,*)'LMEL(KK),KK',LMEL(I),I
C  VEZA SA KONTAKTOR TELOM
C      DO 7 I=1,3
C   7  LMEL(I+3)=ID(JJ,I)
C SAMO KONTROLA U FAZI TESTIRANJA
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1   WRITE(IZLAZ,5030)JJ,(ID(JJ,I),I=1,3)
C     IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
C    1   WRITE(IZLAZ,5030)JJ,(IDC(NLM,I),I=1,3)
C     FORMIRANJE VISINA STUBOVA
C       CALL VISINE(A(LMHT),ND,LMEL)
C
C  PETLJA PO CVROVIMA CILJNOG SEGMENTA
C       ISRF=NEL(NLM,2)
       NNC =NCVSF(ISRF)
       LL  =(IABS(ITSRF(ISRF))-1)*NCVE
        NC3=0
        DO 10 NC=1,NNC
        DO 9  J3=1,NCVE
         NC3=NC3+1
         JJ=NELSF(LL+NC3)
         K1=K2*J3
         IF(JJ.EQ.0)GO TO 9
         DO 20 I=1,K2
           KK=K1+I
           LMEL(KK)=ID(JJ,I)
c          WRITE(3,*)'LMEL(KK),KK',LMEL(KK),KK
   20    CONTINUE
    9   CONTINUE
      ND=(MAXCE+1)*K2
C     FORMIRANJE VISINA STUBOVA
        CALL VISINE(A(LMHT),ND,LMEL)
   10   CONTINUE
C      ND=(MAXCE+1)*K2
CC     FORMIRANJE VISINA STUBOVA
C        CALL VISINE(A(LMHT),ND,LMEL)
  100 CONTINUE
C
C     NEQC=JEDN-NEQ
C
C      WRITE(3,*)'***NEQ,NEQC,JEDN',NEQ,NEQC,JEDN
 5030 FORMAT(11X,I5,3I5)
 2025 FORMAT(//11X,' CVOR    JEDNACINE')
 6025 FORMAT(//11X,' NODE    EQUATIONS')
      RETURN
      END
C=======================================================================
      SUBROUTINE SIS93P(AE,AU,IND)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     GLAVNI UPRAVLJACKI PROGRAM  ZA MATRICE ELEMENATA I SISTEMA
CE     MAIN MANAGEMENT  PROGRAM  FOR ELEMENT MATRIX
C

      include 'paka.inc'
       
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ELEMAU/ MXAU,LAU,LLMEL,LNEL,LNMAT,LTHID,LIPGC,LIPRC,LISNA,
     1 LMXAU,LAPRS
      COMMON /ELEMAE/ MXAE,LAE,LMXAE,LHE,LBET,LBED,LRTHE,LSKE,LLM
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /ZAPISI/ LSTAZA(5)
      COMMON /DUPLAP/ IDVA
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /KONTK3/ FMSTAT,FMDIN,EPSIL,NTSF,NTTOT,LNCVSF,LITSRF,
     &                LNELSF,LIDC,LIK,LIK1,LFSFD,LMASE,LNELAB,NWKCDY
      COMMON /EPUREP/ LPUU,LPUV,LPUA,IPUU,IPUV,IPUA,ISUU,ISUV,ISUA
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /UPRIRI/ LUPRI
      COMMON /PROBAS/ IILS
      COMMON /ITERBR/ ITER
      COMMON /CDEBUG/ IDEBUG
      COMMON /PENALTY/ AKN,AKS,IPENALTY
C
      DIMENSION AE(*),AU(*)
      REAL AE,AU
      IF(IDEBUG.GT.0) PRINT *, ' SIS93'
C  DOPUNSKI REPERI U VEKTORU A
      LIK2 =LIK1 +NE
      LALFK=LIK2 +NE
      KORD=LCORD
      LNGLOB=LMASE +(NE+NTTOT)*IDVA
      LJDIAG=LNGLOB+NTTOT
      LXYZ=1
      LHE =LXYZ+3*(NCVE+1)*IDVA
      LTRA=LHE+27*IDVA
      LAE0=LTRA+9*IDVA
C... PROVERITI ZASTO ZEZA KOD IMPACT I ZA STA JE NEOPHODNO !!!!!
CC      IF(KOR.EQ.1.AND.ITER.EQ.0) CALL CLEAR(A(LRCTDT),JEDN)
CC      CALL CLEAR(A(LRCTDT),JEDN)
      CALL CLEAR(AE,LAE0/IDVA)
C
CS            P E T L J A    P O    E L E M N T I M A
CE            L O O P    O V E R    E L E M E N T S
C
      NCVE36=NCVE*3+6
      IF(IND.EQ.2.OR.IND.EQ.3)THEN
C       REPERI U VEKTORU ELEMENATA AE
       NWE =120
       LSKE=LAE0
       LLM =LSKE+NWE*IDVA
       MXAE=LLM +NCVE36-1
C
C       call WRR(A(LRTDT),JEDN,'U   ')
C       call WRR(A(LUPRI),JEDN,'UPRI')
       ICCMOV=0
       DO 100 NLM=1,NE
       CALL CLEAR(AE,MXAE/IDVA)
       CALL ELTE93P(AE(LSKE),AU(LNCVSF),AU(LITSRF),AU(LNELSF),AU(LISNA),
     &             AE(LLM),AU(LIDC),AU(LNEL),A(KORD),A(LID),A(LRTDT),
     &             A(LRCTDT),A(LUPRI),A(LIK),A(LIK1),A(LALFK),AU(LFSFD),
     &             A(LSIGMA),AU(LNELAB),AE(LXYZ),AE(LHE),AE(LTRA),EPSIL,
     &             ICCMOV)
  100  CONTINUE
C       call WRR(A(LFTDT),JEDN,'FTDT')
C       call WRR(A(LRTDT),JEDN,'RTDT')
C       call WRR(A(LRCTDT),JEDN,'RC  ')
C       WRITE(3,*)'NWK,NWKC,NEQC,NEQ',NWK,NWKC,NEQC,NEQ
C       CALL IWRR(A(LMAXA),JEDN+1,'MAXA')
C       call WRR(A(LSK),NWK,'SKC ')
CE    STORE DISPLACEMENTS ON DISK
CS    ZAPISIVANJE POMERANJA NA DISK
C070794
      IF(IILS.NE.-1)CALL WSTAZ(A(LIPODS),LRTDT,52)
      ELSEIF(IND.EQ.5)THEN
C
CS     PRIPADAJUCE MASE CVOROVA
CE     NODES MASSES
C
C       call WRR(A(LSK),NWK,'MASE')
C PROMENITI AKO TREBA ZA PENALTY  !!!
       CALL MASE93(ALSK,A(LMAXA),A(LID),AU(LMASE),AU(LNEL),AU(LNGLOB),
     &             NP,NE,NTTOT)
C
       LMA8=LSTAZA(2)-1
       CALL WRITDD(AU(LFSFD),MXAU/IDVA,IELEM,LMA8,LDUZI)
       RETURN
      ELSEIF(IND.EQ.7.OR.IND.EQ.8)THEN
C       REPERI U VEKTORU ELEMENATA AE
       LS  =LAE0
       LA  =LS+NCVE*NCVE*IDVA
       LC  =LA+NWKCDY*IDVA
       LB  =LC+NWKCDY*IDVA
       LLM =LB+NTTOT*3*IDVA
       LLD =LLM+NCVE36
       MXAE=LLD+NCVE-1
       MCLR=(LLM-LA)/IDVA
       CALL CLEAR(AE(LA),MCLR)
C
CS     KOREKCIJA UBRZANJA I BRZINA PRI UDARU
CE     UPDATE ACCELERATION AND VELOCITY AFFTER IMPACT
C
C   LPUU --> UBRZANJE ITERACIJA (I-1), BRZINA TRENUTAK (T)
C   LRR  --> ( LPUU ) + KOREKCIJE
       LRR =LPUU
       LSIG=LSIGMA
       LII =LIK2
       IF(IND.EQ.8)THEN
C*
         LPUU =LPUV
C010795         LPUU =LPUV
C*
         LRR =LPUV
         LSIG=LSIGMA+(NE+NTTOT)*3*IDVA
         LII =LIK
       ENDIF
C       call WRR(A(LPUU),JEDN,'PUU-')
C       call WRR(A(LRCTDT),JEDN,'RC- ')
       CALL VAUPD3P(A(LRR),A(LPUU),A(LID),A(LRCTDT),AU(LMASE),AE(LA),
     &             AE(LC),AE(LB),AU(LISNA),A(LII),A(LIK1),A(LALFK),
     &             AU(LNELSF),AU(LITSRF),AU(LIDC),AU(LNEL),AU(LNCVSF),
     &             A(LSIG),AE(LLM),A(LRTDT),NP,NE,NTTOT,NTSF,DT,IND,
     &             A(LIK),AU(LJDIAG),AU(LNGLOB),AU(LNELAB),
     &             AE(LXYZ),AE(LHE),AE(LTRA),AE(LS),AE(LLD),NCVE)
C       call WRR(A(LPUU),JEDN,'PUU+')
C       call WRR(A(LRCTDT),JEDN,'RC+ ')
       CALL IJEDN1(A(LIK2),A(LIK1),NE)
C
CS     KOREKCIJA PRIRASTAJA SILA
CE     UPDATE FORCE INCREMENT
C
      ELSEIF(IND.EQ.9)THEN
        LLL=LUPRI
        IF(METOD.GT.5) LLL=N9
        CALL SILE93(AU(LIDC),A(LLL),A(LIK1),A(LALFK),NE)
      ENDIF
C
C      IF(AKS.GT.0.D0) THEN
C	  NPROS=(NE*3+1)/2*2/IDVA+8*NE
C	ELSE
        NPROS=(NE*3+1)/2*2/IDVA+5*NE
C      ENDIF
	LMA8=LSTAZA(3)-1
      IF(IILS.EQ.-1)RETURN
      CALL WRITDD(A(LIK),NPROS,IELEM,LMA8,LDUZI)
C
      NPROS =(NE+NTTOT)*6
      LMA8=LSTAZA(5)-1
C      call WRR(A(LSIGMA),NPROS,'SIG2')
      CALL WRITDD(A(LSIGMA),NPROS,IELEM,LMA8,LDUZI)
C
      IF(IND.EQ.2)
     & WRITE(IZLAZ,*)'*** CONTACT NODES CHANCHED STATUS:',ICCMOV,'  ***'
C
      RETURN
      END
C=======================================================================
      SUBROUTINE ELTE93P(SKE,NCVSF,ITSRF,NELSF,ISNA,LM,IDC,
     &                  NEL,CORD,ID,U,RC,UPRI,IK,IK1,ALFK,FSFD,SILE,
     &                  NELAB,XYZ,HE,TRA,EPSIL,ICCMOV)
      USE MATRICA
      USE DRAKCE8
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     FORMIRANJE MATRICA ELEMENATA I SISTEMA
CE     FORM ELEMENT MATRIX
C
C  IK     INDIKATOR KONTAKTA U TRENUTKU (T)
C  IK1    INDIKATOR KONTAKTA U TRENUTKU (T+DT) (U ITERACIJI)
C  NCAA   REDNI BROJ POLIGONA NA CILJNOJ POVRSINI

      include 'paka.inc'
       
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /ELEMEN/ ELAST(6,6),XJ(3,3),ALFA(6),TEMP0,DET,NLM,KK
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /PROBAS/ IILS
      COMMON /ITERBR/ ITER
      COMMON /PENALTY/ AKN,AKS,IPENALTY
C
      DIMENSION SKE(*),NCVSF(*),ITSRF(*),NELSF(NCVE,*),LM(*),NEL(NE,*),
     1          CORD(NP,*),U(*),RC(*),ID(NP,*),FSFD(NE,*),IK(*),IK1(*),
     2          ALFK(5,*),SILE(3,*),NELAB(NCVE,*),ISNA(*),XYZ(3,*),
     3         HE(9,*),TRA(3,*),UPRI(*),TSG(5,15),ETP(5,5),ETPG(15,15)
	DIMENSION ETPV(6,6),TSGV(6,6),EN(15),AMETRIC(2,2),DT1(15),DT2(15)
      DIMENSION DC(3),D(3),DL(3),LM3(15),ILM(3),UPRI0(3),TTRRSN(3)
     1,EKTSL(15,15),RNTR(3),RNTR1(3),EKT(15),DT(15,2),TNR1(3),TNR(3)
      DIMENSION TAURSN(3,5),TAU123(3,5),AMETRIC1(2,2),TTR123(3),DG(2)
      DATA IT1/1000000/,SMALL/1.D-10/,BIG/1.D+10/,TCMAX/0.D0/
C      call IWRR(IK,NE,'IK  ')
C      call IWRR(IK1,NE,'IK1 ')
C
C  CVOR KONTAKTOR
C
      NPOL=1
      NCK=NEL(NLM,1)
C        WRITE(3,*)'*********NCK',NCK
C
C  KOEFICIJENTI TRENJA
C
      FS=FSFD(NLM,1)
      FD=FSFD(NLM,2)
C indikator ulaska u kontakt 
C  INDS=1 cvor prvi put ulazi u kontakt
C  INDS=0 cvor je bio u kontakt u prethodnoj iteraciji/koraku
C
	INDS=1
C
C   JEDINICNI VEKTOR PRIRASTAJA POMERANJA
      DO 2 K=1,3
      KK=ID(NCK,K)
      UPRI0(K)=0.D0
2     IF(KK.GT.0) UPRI0(K)=UPRI(KK)
      IF(DABS(UPRI0(1)).GT.SMALL.OR.DABS(UPRI0(2)).GT.SMALL.OR.
     &   DABS(UPRI0(3)).GT.SMALL) CALL JEDV(UPRI0(1),UPRI0(2),UPRI0(3))
C      WRITE(3,*)' UPRI0 ',UPRI0
C
C
CS...... PETLJA PO CILJNIM POLIGONIMA
CE...... LOOP OVER TARGET SEGMENTS
C
C  ISRF - CILJNI SEGMENT, NNC - BROJ CVOROVA SEGMENTA
      ISRF=NEL(NLM,2)
      MAXCE=ISNA(ISRF)
      NODMAX=MAXCE+1
      NNC =NCVSF(ISRF)
      LL =IABS(ITSRF(ISRF))-1
C  PROVERA STANJA U PRETHODNOJ ITERACIJI
      IND =IK1(NLM)/IT1
      IND1=IK(NLM)/IT1
C  KOREKCIJA KOORDINATA ZA CVOR KONTAKTOR
      ICAL=2
      CALL UPXYZ(XYZ,NCK,NELSF,CORD,U,ID,NP,MAXCE,ICAL)
      ICAL=1
      ICHEL=0
      IF(IND.GT.0) THEN
	  INDS=0
        NCAA=IK1(NLM)-IND*IT1
        NCA =NCAA+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
        CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
      NNOD=0
      XL2M=1.D15
      XL2=0.D0
       DO 4 NN=1,MAXCE
         DO 3 L=1,3
           D(L)=XYZ(L,NODMAX)-XYZ(L,NN)
    3    CONTINUE
         XL2=D(1)*D(1)+D(2)*D(2)+D(3)*D(3)
         IF(XL2.LT.XL2M)THEN
          XL2M=XL2
          NNOD=NN
         ENDIF
    4  CONTINUE
        IFLAG=1
        CALL ALFC3P (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &        TOLPD,NNOD,TCMAX,NPOL,NPOL,NCTC,NPOLMX,INDD,IFLAG,AMETRIC)
        IF(INDD.EQ.0)THEN
          IND=0
          ICHEL=1
 	    INDS=1
        ENDIF
        IOVR=0
      ENDIF
C      WRITE(3,*)'IND,IND1,INDD,IK1(NLM)/IT1',IND,IND1,INDD,IK1(NLM)/IT1
C
C   NA KRAJU PRETHODNE ITERACIJE NIJE BILO KONTAKTA
C
      IF(IND.EQ.0) THEN
C
C  1) NALAZENJE NAJBLIZEG CVORA   (NNOD)
C
      NNOD=0
      XL2M=1.D15
      DO 10 NC=1,NNC
       NC2 =NC+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
       CALL UPXYZ(XYZ,NCK,NELSF(1,NC2),CORD,U,ID,NP,MAXCE,ICAL)
       XL2=0.D0
C rastojanje slave cvora do master cvora 
       DO 7 NN=1,MAXCE
         DO 5 L=1,3
           D(L)=XYZ(L,NODMAX)-XYZ(L,NN)
    5    CONTINUE
         XL2=D(1)*D(1)+D(2)*D(2)+D(3)*D(3)
         IF(XL2.LT.XL2M)THEN
          XL2M=XL2
          NNOD=NELSF(NN,NC2)
         ENDIF
    7  CONTINUE
   10 CONTINUE
C
C  2)     NALAZENJE POLIGONA   (NCAMIN)
C
      IFLAG=1
C BROJ ELEMENATA KOJE SADRZE NAJBLIZI CVOR
	NBROJNN=0
      NCAMIN=1
      NNMIN=1
      XL2M=1.D15
	TRA1=0.D0
	TRA2=0.D0
	DC1=0.D0
	DC2=0.D0
	DC3=0.D0
C
      DO 15 NCAA=1,NNC
       NCA =NCAA+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
       CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
       DO 12 NN=1,MAXCE
         IF(NELSF(NN,NCA).EQ.NNOD)THEN
C         CALL AMIN (TRA,XYZ,MAXCE,DC,PROJEK,
C     &              NN,NCAA,NPOL,NCTC,NPOLMX)
           CALL ALFC3P (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &           TOLPD,NN,TCMAX,NCAA,NPOL,NCTC,NPOLMX,IND,IFLAG,AMETRIC)
c           IF(IND.NE.0)GO TO 16
      IF(NNC.EQ.1.AND.IND.NE.0)GO TO 16
C SNEZA   
c           IF(IND.NE.0)THEN
           IF(IND.EQ.1)THEN
C      WRITE(3,*)'ro,so',R0,S0
              NBROJNN=NBROJNN+1
              IF(PRODOR.LT.XL2M)THEN
               XL2M=PRODOR
               NCAMIN=NCAA
               NNMIN=NN
	         r0min=R0
	         s0min=S0
c      WRITE(3,*)'pre NBROJNN,NNMIN,NCAMIN',NBROJNN,NNMIN,NCAMIN
c      WRITE(3,*)'romin,somin',R0min,S0min
              ENDIF
           ENDIF
C SNEZA
	   ENDIF
   12  CONTINUE
   15 CONTINUE
c      IF(NNC.EQ.1.AND.IND.NE.0)GO TO 16
c      WRITE(3,*)'pre NBROJNN,NNMIN,NCAMIN',NBROJNN,NNMIN,NCAMIN
c      WRITE(3,*)'ro,so',R0min,S0min
       IF(IFLAG.EQ.2)THEN
         NCA =NPOL+LL
         CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
         CALL ALFC3P (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &           TOLPD,NN,TCMAX,NCAA,NPOL,NCTC,NPOLMX,IND,IFLAG,AMETRIC)
       ENDIF
C SNEZA
C 
C ODREDJIVANJE TACNIH KOORDINATA PRODORA NA "AKTIVNOM" SEGMENTU (NCAMIN)
C
	IF (NBROJNN.GT.0) THEN
c      WRITE(3,*)'NBROJNN,NNMIN,NNOD,NCAMIN',NBROJNN,NNMIN,NNOD,NCAMIN
	 NCAA=NCAMIN
       NCA =NCAA+LL
C   KOREKCIJA KOORDINATA ZA CILJNE CVOROVE
       CALL UPXYZ(XYZ,NCK,NELSF(1,NCA),CORD,U,ID,NP,MAXCE,ICAL)
       CALL ALFC3P (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &        TOLPD,NNMIN,TCMAX,NCAA,NPOL,NCTC,NPOLMX,IND,IFLAG,AMETRIC)
C       WRITE(3,*)'PRODOR,XL2M',PRODOR, XL2M
C       WRITE(3,*)'posle NBROJNN,NNMIN,NCAMIN',NBROJNN,NNMIN,NCAMIN
C       WRITE(3,*)'ro,so',R0,S0
      ENDIF
C SNEZA
   16 IOVR=IND
C       WRITE(3,*)'IND, NCAA',IND,NCAA
C.. ICHEL > 0   AKO JE PROMENJEN PODSEGMENT
      ICHEL=IND*ICHEL
C      WRITE(3,*)'nlm,NCK,DC,IND',NLM,NCK,DC,IND
      ENDIF
c       LM3(1)=ID(NLM,1)
c       LM3(2)=ID(NLM,2)
c       LM3(3)=ID(NLM,3)
C   20  IF((ICHEL.NE.0.OR.IND.EQ.0).AND.IILS.NE.-1)THEN
C       ENDIF
C
C  3)     RAZLICITI KONTAKTNI USLOVI
C
C
C-----    (A)   NEMA KONTAKTA       -  za PENALTY netreba nista
C Sneza
C      MDIM=MAXCE*3+6
C      IF(IND.EQ.0)THEN
C       IF(IILS.EQ.-1)RETURN
C       MDIM=3
C       SKE(1)=-1.D0
C       SKE(4)=-1.D0
C       SKE(6)=-1.D0
C
C       CALL SPAKUJ(A(LSK),A(LMAXA),SKE,LM3,MDIM)
C
C       DO 30 J=1,3
C        IJ=LM3(J)
C        IF(IJ.GT.0)U(IJ)=0.D0
C        SILE(J,NLM)=0.D0
C   30  CONTINUE
C       IK (NLM)=0
C       IK1(NLM)=0
C       RETURN
C      ENDIF
C Sneza
C  PROVERA ZA OSLOBADJANJE IZ KONTAKTA
C
C      IF(TAURSN(3,NODMAX).LT.-SMALL.AND.PRODOR.GT.-TOLPD)THEN
c      IF(PRODOR.GT.-TOLPD)THEN
C 22.01.06 promenjeno kontakt je uvek kad ima prodora
C       PSMALL=9.99999998D-12
      IF(PRODOR.GT.PSMALL)THEN
c      IF(PRODOR.GT.0.D0)THEN
        IND=0
C        GO TO 20
C       WRITE(3,*)'izlaz IND, prodor',IND,PRODOR
       IK (NLM)=0
       IK1(NLM)=0
       RETURN
      ENDIF
C       WRITE(3,*)'IND, prodor',IND,PRODOR
C
C-----    (B)   KONTAKT SA ILI BEZ KLIZANJA
C
c      KK=-MAXCE-1
c      KK=-MAXCE
c      LMILM=0
c      DO 35 K=1,3
Cc        KK=KK+MAXCE+2
c        KK=KK+MAXCE+1
c       LM(KK) =ID(NCK,K)
Cc Sneza   ne treba za Penalty
Cc        ILM(K)  =IDC(NLM,K)
Cc        LM(KK+1)=ILM(K)
Cc        IF(LMILM.EQ.0) LMILM=ILM(K)
Cc        K1=KK+1
Cc Sneza
c        K1=KK
c         DO 32 L=1,MAXCE
c          NCV=NELSF(L,NCA)
cC      IF(K.EQ.1)WRITE(3,*)'NCV',NCV
c   32     LM(K1+L)=ID(NCV,K)
c         DO 34 J=1,NODMAX
c   34     TAU123(K,J)=0.D0
c   35 CONTINUE
C
C       
C           PROVERA TIPA KONTAKTA
C       
C  	za Penalty ako je unet koeficijent trenja kontakt je sa klizanjem
C	proverava se da li je stick ili slip tipa kontakta 
C
c      IND=1
CC... RRS = REZULTANTA U  RS RAVNI
C      RRS=DSQRT(TAURSN(1,NODMAX)*TAURSN(1,NODMAX)+
C     &          TAURSN(2,NODMAX)*TAURSN(2,NODMAX))
CCC      IF(RRS.GE.FS*DABS(TAURSN(3,NODMAX))) IND=2
CC070794      IF(RRS.GE.FS*DABS(TAURSN(3,NODMAX)).AND.DABS(FS).LT.BIG) IND=2
C      RN=DABS(TAURSN(3,NODMAX))
CC100794      IF(RN.GT.SMALL.AND.RRS.GE.FS*RN.AND.DABS(FS).LT.BIG) IND=2
C      IF(ICHEL.EQ.0.AND.
C     &   RN.GT.SMALL.AND.RRS.GE.FS*RN.AND.DABS(FS).LT.BIG) IND=2
C      IF(DABS(FS).LT.SMALL) IND=2
CC      CALL WRR(TAURSN(1,NODMAX),3,'*RSN')
CC      WRITE(3,*)' IND',IND
CC
C
C          B1)      KONTAKT BEZ KLIZANJA	 za Penalty za sada samo ovo (30.03.04)
C
C
C      IF(IND.EQ.1) THEN
C                        za Penalty za sada radi sve kao kontakt bez klizanja (12.04.04)
      IF (IND.EQ.2) IND=1
      IF (IND.EQ.1) THEN
C
C RACUNANJE TRENUTNE DUZINE SEGMENTA (Sneza)
       IF(MAXCE.EQ.2) THEN
         DO L=1,2
           DL(L)=XYZ(L,2)-XYZ(L,1)
         ENDDO
         XLT=DL(1)*DL(1)+DL(2)*DL(2)
         XLT=DSQRT(XLT)
c	WRITE(3,*)'XLT',XLT
       ENDIF
C
C
         IF(MAXCE.EQ.2) THEN
           K2=2
         ELSE
           K2=3
         ENDIF
        ii=MAXCE+1
        jj=ii*K2
C	WRITE(3,*)'II,JJ,K2',ii,jj,K2
      DO 35 K=1,K2        
C     DO 35 K=1,3        
         LM(K) =ID(NCK,K)  
C       WRITE(3,*)'LM(KK),KK',LM(K),K
         DO 34 J=1,NODMAX
   34     TAU123(K,J)=0.D0
   35 CONTINUE
        DO 32 L=1,MAXCE
          NCV=NELSF(L,NCA)
          K1=K2*L
          DO 33 K=1,K2
           KK=K1+K
   33        LM(KK)=ID(NCV,K)
C   33   WRITE(3,*)'LM(KK),KK',LM(KK),KK
   32  CONTINUE

C          KONTAKTNA KRUTOST 
CC       MDIM=MAXCE+2
C       MDIM=MAXCE+1
      IF(IILS.NE.-1)THEN
C       KK=-MAXCE
C       DO 56 K=1,3
C        DUM=TRA(K,3)
C      WRITE(3,*)'DUM, NLM, K',DUM,NLM,K
C        SKE(1) =DUM*AKN
C      WRITE(3,*)'SKE(1), 1',SKE(1)
CC       SKE(1)=-AKN
CC       LLL=kk+MDIM
C       LLL=1
C       DO 55 L=1,MAXCE
C        LLL=LLL+1
C       SKE(LLL)=-AKN*DUM*HE(L,1)
C      WRITE(3,*)'SKE(LLL), LLL',SKE(LLL), LLL
C   55 CONTINUE
CC
C       DO L=1,MAXCE
C        KL=L
C        DO M=KL,MAXCE
C          LLL=LLL+1
C          SKE(LLL)=AKN*DUM*HE(L,1)*HE(M,1)
C      WRITE(3,*)'SKE(LLL), LLL, L, M',SKE(LLL), LLL, L, M
C        ENDDO
C       ENDDO
   
C      CALL WRR(HE(1,1),MAXCE,'HE  ')
C
C       KK=-MAXCE-1
C       KK=-MAXCE
C       DO 56 K=1,3
C         KK=KK+MAXCE+2
C         KK=KK+MAXCE+1
C         CALL SPAKUJ(A(LSK),A(LMAXA),SKE,LM(KK),MDIM)
C   56  CONTINUE
C stara varijanta sa matricom transformacije
c	DO I=1,5
c	  DO J=1,5
c	    ETP(I,J)=0
c	  ENDDO
c	ENDDO
C     KONTAKTNA KRUTOST U LOKALNOM SISTEMU (NNtransponovano)      
C
c       KL=0
cc       ETP(1,1)=AKN
c       ETP(1,1)=1.0
c       DO L=1,MAXCE
c	    KL=L+1
c          ETP(1,KL)=-HE(L,1)
cc          ETP(1,KL)=-AKN*HE(L,1)
Cc          WRITE(3,*)'ETP(1,KL),KL',ETP(1,KL),KL
c          DO LL=L,MAXCE
c	       KLL=LL+1
c	       KL=L+1
c             ETP(KL,KLL)=HE(L,1)*HE(LL,1)
cc             ETP(KL,KLL)=AKN*HE(L,1)*HE(LL,1)
Cc           WRITE(3,*)'ETP(KL,KLL),KL,KLL',ETP(KL,KLL),KL,KLL
c          ENDDO
c       ENDDO
c      DO L=2,MAXCE+1
c          DO LL=1,L-1
c             ETP(L,LL)=ETP(LL,L)
Cc       WRITE(3,*)'ETP(L,LL),L,LL',ETP(L,LL),L,LL
c          ENDDO
c       ENDDO
C    MATRICA TRANSFORMACIJE TSG
c      DO I=1,ii
c         DO J=1,jj
c         TSG(I,J)=0.0
c         ENDDO
c      ENDDO
c      K1=0
c      DO L=1,MAXCE+1
Cc            K1=K1+K2
c         DO I=1,K2
c	      K1=K1+1
c            TSG(L,K1)=TRA(I,3)
CC          WRITE(3,*)'TSG(L,K1),L,K1',TSG(L,K1),L,K1
c         ENDDO
c      ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
c      CALL TRAETPP(ETP,ETPG,TSG,ii,jj)      
C
	 DO J=1,jj
	    EN(J)=0.D0
	 ENDDO
       DO I=1,K2
         EN(I)=TRA(I,3)
c          WRITE(3,*)'EN,K',EN(I),I
       ENDDO
c        K1=K2
       DO L=1,MAXCE
          K1=K2*L
        DO K=1,K2
           KK=K1+K
           EN(KK)=-HE(L,1)*TRA(K,3)
c          WRITE(3,*)'EN,K',EN(KK),KK
        ENDDO
	 ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
      CALL MNOZI(EN,ETPG,jj)      
C   FORMIRANJE MATRICE SKE
      KL=0
      DO 600 L=1,jj
        DO 600 K=L,jj
           KL=KL+1
c           SKE(KL)=ETPG(L,K)
C        WRITE(3,*)'ETPG(I,J),I,J',ETPG(L,K),L,K
           SKE(KL)=AKN*ETPG(L,K)
C       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  600  CONTINUE
      MDIM=jj
C   PAKOVANJE KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
      IF (TIPTACKANJA.EQ.1) THEN
         CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
         ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C
C VELIKE DEFORMACIJE SE UZIMAJU U OBZIR (2D kontaktna krutost u pravcu noramle) 16.01.06
        IF((IATYP.GE.3).AND.(MAXCE.EQ.2)) THEN
	    DO I=1,6
	      DO J=1,6
	        ETPV(I,J)=0.0
	        TSGV(I,J)=0.0
	      ENDDO
	   ENDDO
C AKN*gn/L*(N0Tt+TN0t)
         TNL=PRODOR*AKN/XLT
         ETPV(1,5)=TNL
         ETPV(5,1)=ETPV(1,5)
         ETPV(1,6)=-TNL
         ETPV(6,1)=ETPV(1,6)
         ETPV(2,5)=-HE(1,1)*TNL
         ETPV(5,2)=ETPV(2,5)
         ETPV(2,6)=HE(1,1)*TNL
         ETPV(6,2)=ETPV(2,6)
         ETPV(3,5)=-HE(2,1)*TNL
         ETPV(5,3)=ETPV(3,5)
         ETPV(3,6)=HE(2,1)*TNL
         ETPV(6,3)=ETPV(3,6)
C AKN*gn/L*gn/L*N0N0t
         ETPV(5,5)=-TNL*(PRODOR/XLT)
         ETPV(6,6)=ETPV(5,5)
         ETPV(5,6)=TNL*(PRODOR/XLT)
         ETPV(6,5)=ETPV(5,6)
C
      K1=0
      DO L=1,MAXCE+1
         DO I=1,K2
	      K1=K1+1
            TSGV(L,K1)=TRA(I,2)
            TSGV(L+ii,K1)=TRA(I,3)
C          WRITE(3,*)'TSGV(L,K1),L,K1',TSGV(L,K1),L,K1
         ENDDO
      ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
         CALL TRAETPP(ETPV,ETPG,TSGV,jj,jj)      
C   FORMIRANJE MATRICE SKE
         KL=0
         DO 700 L=1,jj
           DO 700 K=L,jj
             KL=KL+1
             SKE(KL)=-ETPG(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  700  CONTINUE
      MDIM=jj
C   PAKOVANJE KONTAKTNE KRUTOSTI U GLOBALNI SISTEM 
      IF (TIPTACKANJA.EQ.1) THEN
         CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
         ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C
	  ENDIF
C KRAJ DELA ZA VELIKE DEFORMACIJE (2D)
C
      ENDIF
C
C klizanje/slepljivanje za penalty (AKS>0)
C
      IF (AKS.GT.0.D0)THEN
	  ISLIP=0
c
       write(3,*) 'inds', INDS
C EKT TANGENTNA "MATRICA KRUTOSTI"
 	  DO J=1,jj
	    DT1(J)=0.D0
	    DT1(J)=0.D0
c	    EKT(J)=0.D0
C          DO K=1,2
c	      DT(J,K)=0.D0
c	    ENDDO
	    DO I=1,jj
c            EKTSL(I,J)=0.D0
	      ETPG(I,J)=0.D0
	      EKTSL(I,J)=0.D0
	    ENDDO
	  ENDDO
        DO J=1,3
c         TTRRSN(I)=0.D0
          RNTR(J)=0.D0
		RNTR1(J)=0.D0
		TNR1(J)=0.D0
		TNR(J)=0.D0
        ENDDO
C
C 2D KONTAKT
C
	  AMETRIC1(1,1)=0.0
	  AMETRIC1(1,2)=0.0
	  AMETRIC1(2,1)=0.0
	  AMETRIC1(2,2)=0.0
      IF (MAXCE.EQ.2) THEN
	    DETM=1.0
          AMETRIC1(1,1)=1.0/AMETRIC(1,1)
C provera da li je prva iteracija ili nije bilo u prethodnoj iteraciji kontakta
	  IF ((ITER.EQ.0).OR.(INDS.EQ.1)) THEN
C STICK
C IZRACUNAVANJE Tbeta=EKTJ
          DO I=1,K2
            EKT(I)=TRA(I,2)
          ENDDO
          DO L=1,MAXCE
            K1=K2*L
           DO K=1,K2
            KK=K1+K
            EKT(KK)=-HE(L,1)*TRA(K,2)
C          WRITE(3,*)'ETP(1,KL),KL',ETP(1,KL),KL
           ENDDO 
	    ENDDO
C IZRACUNAVANJE Dalfa=DT
          DO K=1,jj
	        DT1(K)=DT1(K)+AMETRIC1(1,1)*EKT(K)
	    ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
c        CALL MNOZIT(DT,ETPG,AMETRIC,jj)      
        DO 88 I=1,jj
        DO 88 J=1,jj
           ETPG(I,J)=ETPG(I,J)+AMETRIC(1,1)*DT1(I)*DT1(J)
   88   CONTINUE
C   FORMIRANJE MATRICE SKE
        KL=0
       WRITE(3,*)'Stick1 2d /n'           
        DO 620 L=1,jj
          DO 620 K=L,jj
            KL=KL+1
            SKE(KL)=AKS*ETPG(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  620   CONTINUE
         MDIM=jj
         IF (TIPTACKANJA.EQ.1) THEN
         CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
         ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C IZRACUNAVANJE SILE U TANGENCIJALNOJ RAVNI
	   DG(1)=R0-ALFK(1,NLM)
C         DG(2)=S0-ALFK(2,NLM)
             TTRRSN(1)= 0.0
             TTRRSN(2)= 0.0
c             TTRRSN(2)= AKS*AMETRIC(1,1)*DG(1)
c             TTRRSN(2)= TTRRSN(2)+AKS*AMETRIC(1,1)*DG(1)
c       WRITE(3,*)'R0, S0',R0,S0           
c       WRITE(3,*)'DG(1), DG(2)',DG(1),DG(2)           
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
c	  ENDIF
C
       ELSE
C IZRACUNAVANJE "TRIAL STATE" - probno stanje
C
         DTR=0
	   DG(1)=R0-ALFK(1,NLM)
C         DG(2)=S0-ALFK(2,NLM)
         TTR123(1)=ALFK(3,NLM)
         TTR123(2)=ALFK(4,NLM)
         TTR123(3)=ALFK(5,NLM)
	   CALL RSN123P(TTRRSN(1),TTR123(1),TRA,2)
C jedinicni vektor normale na povrs je prvi vektor
         TTRRSN(2)= TTRRSN(2)+AKS*AMETRIC(1,1)*DG(1)
	   DTR=0.D0
c       WRITE(3,*)'R0, S0',R0,S0           
c       WRITE(3,*)'DG(1), DG(2)',DG(1),DG(2)           
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
         DTR=DTR+TTRRSN(2)*AMETRIC1(1,1)*TTRRSN(2)
         DTR=DSQRT(DTR)
	   FTRIAL=DTR-FS*AKN*PRODOR
	   IF(FTRIAL.LE.0.D0)THEN
C STICK
C IZRACUNAVANJE Tbeta=EKT
             DO I=1,K2
               EKT(I)=TRA(I,2)
             ENDDO
             DO L=1,MAXCE
               K1=K2*L
              DO K=1,K2
               KK=K1+K
               EKT(KK)=-HE(L,1)*TRA(K,2)
              ENDDO
	       ENDDO
C IZRACUNAVANJE Dalfa=DT
          DO K=1,jj
	        DT1(K)=DT1(K)+AMETRIC1(1,1)*EKT(K)
	    ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
c        CALL MNOZIT(DT,ETPG,AMETRIC,jj)      
        DO 89 I=1,jj
        DO 89 J=1,jj
           ETPG(I,J)=ETPG(I,J)+AMETRIC(1,1)*DT1(I)*DT1(J)
   89   CONTINUE
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
C          CALL MNOZIT(DT,ETPG,AMETRIC,jj)      
C   FORMIRANJE MATRICE SKE
          KL=0
       WRITE(3,*)'Stick 2d /n'           
          DO 622 L=1,jj
            DO 622 K=L,jj
              KL=KL+1
              SKE(KL)=AKS*ETPG(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  622     CONTINUE
            MDIM=jj
            IF (TIPTACKANJA.EQ.1) THEN
            CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
            ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C sila
             TTRRSN(2)= TTRRSN(2)*AMETRIC1(1,1)
C
	   ELSE
C
C SLIP
C SILE U TANGENTNOJ RAVNI
C vektor napona za slucaj klizanja
         RNTR(1)=0.0
         RNTR(2)=TTRRSN(2)/DTR
C TRANSFORMACIJA VEKTORA NAPONA U GLOBALNI SISTEM
	   CALL RSN123P(RNTR(1),TNR(1),TRA,1)
	   TTRRSN(1)=FS*AKN*PRODOR*RNTR(1)*AMETRIC1(1,1)
	   TTRRSN(2)=FS*AKN*PRODOR*RNTR(2)*AMETRIC1(1,1)
C	   TTRRSN(1)=FS*AKN*PRODOR*RNTR(1)
C	   TTRRSN(2)=FS*AKN*PRODOR*RNTR(2)
         RNTR1(1)=0.D0
         RNTR1(2)=0.D0
C "inverzni" vektor napona
	       RNTR1(2)=RNTR1(2)+RNTR(2)*AMETRIC1(1,1)
C TRANSFORMACIJA inverznog VEKTORA NAPONA U GLOBALNI SISTEM
	   CALL RSN123P(RNTR1(1),TNR1(1),TRA,1)
C MATRICA "KRUTOSTI"
C IZRACUNAVANJE Tbeta=EKT
             DO I=1,K2
               EKT(I)=TRA(I,2)
             ENDDO
             DO L=1,MAXCE
               K1=K2*L
              DO K=1,K2
               KK=K1+K
               EKT(KK)=-HE(L,1)*TRA(K,2)
              ENDDO
	       ENDDO
C IZRACUNAVANJE Dalfa=DT
          DO K=1,jj
	        DT1(K)=DT1(K)+AMETRIC1(1,1)*EKT(K)
	    ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
         DO 82 I=1,jj
         DO 82 J=1,jj
             EKTSL(I,J)=EKTSL(I,J)+FS*AKN*RNTR(2)*DT1(I)*EN(J)
     1+FS*AKN*PRODOR*AKS/DTR*AMETRIC(1,1)*(DELTA(1,1)
     1        -RNTR(2)*RNTR1(2))*DT1(I)*DT1(J)
c             EKTSL(I,J)=EKTSL(I,J)+FS*AKN*TNR(2)*DT1(I)*EN(J)
c     1+FS*AKN*PRODOR*AKS/DTR*AMETRIC(1,1)*(DELTA(1,1)
c     1        -TNR(2)*TNR1(2))*DT1(I)*DT1(J)
C     1        -RNTR(2)*RNTR1(2))*DT1(I)*DT1(J)
   82    CONTINUE
C   FORMIRANJE MATRICE SKE
          KL=0
       WRITE(3,*)'Slip 2d /n'           
          DO 633 L=1,jj
            DO 633 K=L,jj
              KL=KL+1
              SKE(KL)=EKTSL(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  633     CONTINUE
            MDIM=jj
            IF (TIPTACKANJA.EQ.1) THEN
            CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
            ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
	   ENDIF
	  ENDIF
c	 ENDIF
C
C 3D KONTAKT STICK/SLIP
C
	 ELSEIF(MAXCE.GT.2)THEN
	    DETM=AMETRIC(1,1)*AMETRIC(2,2)-AMETRIC(2,1)*AMETRIC(1,2)
	    AMETRIC1(1,1)=AMETRIC(2,2)/DETM
	    AMETRIC1(1,2)=-AMETRIC(2,1)/DETM
	    AMETRIC1(2,1)=-AMETRIC(1,2)/DETM
	    AMETRIC1(2,2)=AMETRIC(1,1)/DETM
c         WRITE(3,*)'DETM',DETM
C provera da li je prva iteracija ili nije bilo u prethodnoj iteraciji kontakta
	  IF ((ITER.EQ.0).OR.(INDS.EQ.1))THEN
C STICK
        DO J=1,K2-1
C IZRACUNAVANJE Tbeta=EKTJ
          DO I=1,K2
            EKT(I)=TRA(I,J)
          ENDDO
          DO L=1,MAXCE
            K1=K2*L
            DO K=1,K2
             KK=K1+K
             EKT(KK)=-HE(L,1)*TRA(K,J)
C            WRITE(3,*)'ETP(1,KL),KL',ETP(1,KL),KL
            ENDDO
	    ENDDO
C IZRACUNAVANJE Dalfa=DT
          DO K=1,jj
            DO M=1,2
	        DT(K,M)=DT(K,M)+AMETRIC1(J,M)*EKT(K)
C	        DT(K,M)=DT(K,M)+AMETRIC1(M,J)*EKT(K)
	      ENDDO
	    ENDDO
	  ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
        CALL MNOZIT(DT,ETPG,AMETRIC,jj)      
C   FORMIRANJE MATRICE SKE
        KL=0
        WRITE(3,*)'Stick1 3d /n'           
        DO 621 L=1,jj
          DO 621 K=L,jj
            KL=KL+1
            SKE(KL)=AKS*ETPG(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  621   CONTINUE
         MDIM=jj
         IF (TIPTACKANJA.EQ.1) THEN
         CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
         ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C	  ENDDO
C IZRACUNAVANJE SILE U TANGENCIJALNOJ RAVNI
C R0 i S0 = KSI(n+1)
             TTRRSN(1)= 0.0
             TTRRSN(2)= 0.0
	   DG(1)=R0-ALFK(1,NLM)
         DG(2)=S0-ALFK(2,NLM)
C         DO I=1,2
C	     DO J=1,2
C             TTRRSN(I)= TTRRSN(I)+AKS*AMETRIC(I,J)*DG(J)
C	     ENDDO
C         ENDDO
c	  ENDIF
C
       ELSE
C IZRACUNAVANJE "TRIAL STATE" 
C
         DTR=0
	   DG(1)=R0-ALFK(1,NLM)
         DG(2)=S0-ALFK(2,NLM)
         TTR123(1)=ALFK(3,NLM)
         TTR123(2)=ALFK(4,NLM)
         TTR123(3)=ALFK(5,NLM)
	   CALL RSN123P(TTRRSN(1),TTR123(1),TRA,2)
c       WRITE(3,*)'R0, S0',R0,S0           
c       WRITE(3,*)'DG(1), DG(2)',DG(1),DG(2)           
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
         DO I=1,2
	     DO J=1,2
             TTRRSN(I)= TTRRSN(I)+AKS*AMETRIC(I,J)*DG(J)
	     ENDDO
         ENDDO
c       WRITE(3,*)'3d slip'
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
	   DTR=0.D0
         DO I=1,2
	     DO J=1,2
             DTR=DTR+TTRRSN(I)*AMETRIC1(I,J)*TTRRSN(J)
	     ENDDO
         ENDDO
         DTR=DSQRT(DTR)
	   FTRIAL=DTR-FS*AKN*PRODOR
c       WRITE(3,*)'3d slip, ftrial',FTRIAL           
	   IF(FTRIAL.LE.0.D0)THEN
C STICK
           DO J=1,K2-1
C IZRACUNAVANJE Tbeta=EKT
             DO I=1,K2
               EKT(I)=TRA(I,J)
             ENDDO
             DO L=1,MAXCE
               K1=K2*L
              DO K=1,K2
               KK=K1+K
               EKT(KK)=-HE(L,1)*TRA(K,J)
              ENDDO
	       ENDDO
C IZRACUNAVANJE Dalfa=DT
            DO K=1,jj
              DO M=1,2
	         DT(K,M)=DT(K,M)+AMETRIC1(J,M)*EKT(K)
	        ENDDO
	      ENDDO
	    ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
          CALL MNOZIT(DT,ETPG,AMETRIC,jj)      
C   FORMIRANJE MATRICE SKE
          KL=0
       WRITE(3,*)'Stick 3d /n'           
          DO 623 L=1,jj
            DO 623 K=L,jj
              KL=KL+1
              SKE(KL)=AKS*ETPG(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  623     CONTINUE
            MDIM=jj
            IF (TIPTACKANJA.EQ.1) THEN
            CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
            ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C
         DO I=1,2
	     DO J=1,2
             TTRRSN(I)=TTRRSN(I)+TTRRSN(I)*AMETRIC1(I,J)
	     ENDDO
         ENDDO
c       WRITE(3,*)'R0, S0',R0,S0           
c       WRITE(3,*)'DG(1), DG(2)',DG(1),DG(2)           
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
	   ELSE
C
C SLIP
C
C vektor normale
         RNTR(1)=TTRRSN(1)/DTR
         RNTR(2)=TTRRSN(2)/DTR
C SILE U TANGENCIJALNOJ RAVNI
	   TTRRSN(1)=FS*AKN*PRODOR*RNTR(1)
	   TTRRSN(2)=FS*AKN*PRODOR*RNTR(2)
         RNTR1(1)=0.D0
         RNTR1(2)=0.D0
c       WRITE(3,*)' slip 3d'
c       WRITE(3,*)'RNTR(1), RNTR(2)',RNTR(1),RNTR(2)           
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
C
         DO I=1,2
	     DO J=1,2
             TTRRSN(I)=TTRRSN(I)+TTRRSN(I)*AMETRIC1(I,J)
	     ENDDO
         ENDDO
C TRANSFORMACIJA VEKTORA NAPONA U GLOBALNI SISTEM
	   CALL RSN123P(RNTR(1),TNR(1),TRA,1)
C inverzni vektor napona
	   DO I=1,2
	     DO J=1,2
	       RNTR1(I)=RNTR1(I)+RNTR(J)*AMETRIC1(J,I)
	     ENDDO
	   ENDDO
C TRANSFORMACIJA inverznog VEKTORA NAPONA U GLOBALNI SISTEM
	   CALL RSN123P(RNTR1(1),TNR1(1),TRA,1)
C MATRICA "KRUTOSTI"
           DO J=1,K2-1
C IZRACUNAVANJE Tbeta=EKT
             DO I=1,K2
               EKT(I)=TRA(I,J)
             ENDDO
             DO L=1,MAXCE
               K1=K2*L
              DO K=1,K2
               KK=K1+K
               EKT(KK)=-HE(L,1)*TRA(K,J)
              ENDDO
	       ENDDO
C IZRACUNAVANJE Dalfa=DT
            DO K=1,jj
              DO M=1,2
	        DT(K,M)=DT(K,M)+AMETRIC1(J,M)*EKT(K)
C	          DT(K,M)=DT(K,M)+AMETRIC1(M,J)*EKT(K)
	        ENDDO
	      ENDDO
	    ENDDO
C   TRANSFORMACIJA KONTAKTNE KRUTOSTI U GLOBALNI SISTEM
         DO 81 I=1,jj
         DO 81 J=1,jj
	     DO K=1,2
             EKTSL(I,J)=EKTSL(I,J)+FS*AKN*RNTR(K)*DT(I,K)*EN(J)
c	      DO MB=1,2
	      DO MG=1,2
	      EKTSL(I,J)=EKTSL(I,J)+FS*AKN*PRODOR/DTR*AKS*AMETRIC(K,MG)*
     1        (DELTA(K,MB)-RNTR(K)*RNTR1(MB))*DT(I,K)*DT(J,MB)
c	      EKTSL(I,J)=EKTSL(I,J)+FS*AKN*PRODOR/DTR*AKS*AMETRIC(K,MG)*
c     1        (DELTA(K,MB)-TNR(K)*TNR1(MB))*DT(I,K)*DT(J,MB)
  	      ENDDO
c  	      ENDDO
           ENDDO
   81    CONTINUE
C   FORMIRANJE MATRICE SKE
          KL=0
       WRITE(3,*)'Slip 3d /n'           
          DO 603 L=1,jj
            DO 603 K=L,jj
              KL=KL+1
              SKE(KL)=EKTSL(L,K)
c       WRITE(3,*)'SKE(KL), KL',SKE(KL),KL           
  603     CONTINUE
            MDIM=jj
            IF (TIPTACKANJA.EQ.1) THEN
            CALL SPAKUJ(ALSK,A(LMAXA),SKE,LM,MDIM)
            ELSE
         CALL sparseassembler_addelemmatrix(MDIM,LM,SKE)
            ENDIF
C
         DO I=1,2
	     DO J=1,2
             TTRRSN(I)=TTRRSN(I)+TTRRSN(I)*AMETRIC1(I,J)
	     ENDDO
         ENDDO
c       WRITE(3,*)'TTRRSN(1), TTRRSN(2)',TTRRSN(1),TTRRSN(2)           
	   ENDIF
	 ENDIF
	 ENDIF
	ENDIF
C           KONTAKTNE SILE  - PRAVAC X,Y,Z U KONTAKTORU
C  Gruja
C      KK=-MAXCE
C      DO 45 K=1,3
C       KK=KK+MAXCE+2
C       LMKK=LM(KK)
C       TAU123(K,NODMAX)=0.D0
CC       WRITE(3,*)'IOVR',IOVR
CC       IF(ICHEL.EQ.0)THEN
C         IF(LMKK.NE.0) TAU123(K,NODMAX)=U(LMKK)
CC       ENDIF
C   45 CONTINUE
CC      CALL WRR(TAU123(1,NODMAX),3,'123 ')
CC           TRANSFORMACIJA U LOKALNI SISTEM I UVODJENJE USLOVA KLIZANJA
C***
C*7
CC   AKO JE PRETHODNI STATUS BIO KLIZANJE LOKALNO = GLOBALNO
CC       WRITE(3,*)'************ IND',IND
C*290694       IF(IND.NE.2)THEN
CC       IF(IND1.EQ.1)THEN
C          CALL RSN123(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,2)
CC       ELSEIF(IND1.EQ.2)THEN
CC*9
CC        TAURSN(3,NODMAX)=TAU123(1,NODMAX)
CC*** OSTALA DVA CLANA TREBA DA BUDU IZ PRETHODNOG PROLAZA !!!
CC*9
CC       ENDIF
C
C
C*7 Kraj Gruja
C Sneza
      IF(AKS.GT.0.0) THEN
          TAURSN(1,NODMAX)=TTRRSN(1)
          TAURSN(2,NODMAX)=TTRRSN(2)
C        IF(ITER.EQ.0)THEN
C          TAURSN(1,NODMAX)=0.0
C          TAURSN(2,NODMAX)=0.0
C	  ENDIF
	ELSE
          TAURSN(1,NODMAX)=0.0
          TAURSN(2,NODMAX)=0.0
      ENDIF
	    TAURSN(3,NODMAX)=-AKN*PRODOR
c      IF((-TAURSN(3,NODMAX).LT.SMALL).AND.(PRODOR.LT.SMALL)
c     1.and.(iter.gt.0))THEN
c      IF(PRODOR.GT.0.D0)THEN
c        IND=0
C        GO TO 20
c       WRITE(3,*)'izlaz',TAURSN(3,NODMAX)
c       RETURN
c      ENDIF
c      IF(AKS.GT.0.0) THEN
c          CALL RSN123P(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,1)
c	ELSE
C       do i=1,K2
c        WRITE(3,*)'TAURSN(i,NODMAX),i',TAURSN(i,NODMAX),i
C	 enddo
C TRANSFORMACIJA SILA U GLOBALNI SISTEM
          CALL RSN123(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,1)
c      ENDIF
c         WRITE(3,*)'TAURSN(3,NODMAX)',TAURSN(3,NODMAX)          
C       do i=1,K2
c        WRITE(3,*)'TAU123(i,NODMAX),i',TAU123(i,NODMAX),i
C	 enddo
C           SILE NA CILJNIM CVOROVIMA
c      KK=-MAXCE
      DO 52 K=1,K2
C       KK=KK+MAXCE+2
c       KK=KK+MAXCE+1
c       IF(LM(KK).NE.0) THEN
       IF(LM(K).NE.0) THEN
        DUM=TAU123(K,NODMAX)
c	DUM=DC(K)
        SILE(K,NLM)=DUM
c        WRITE(3,*)'K,NLM,SILE(K,NLM)',K,NLM,SILE(K,NLM)
        DO 51 L=1,MAXCE
          NELL=NE+NELAB(L,NCA)
          DUM1=SILE(K,NELL)-HE(L,1)*DUM
c        WRITE(3,*)'HE',HE(L,1)
          TAU123(K,L)=DUM1
          SILE(K,NELL)=DUM1
C kontrolna stampa (sledeca dva reda)
c          WRITE(3,*)'K,NELL,SILE(K,NELL)',K,NELL,SILE(K,NELL)
C          WRITE(3,*)'K,L,TAU123(K,L)',K,L,TAU123(K,L)
   51   CONTINUE
       ENDIF
   52 CONTINUE
C
      ENDIF
C
C
C
C             B2)      KONTAKT SA KLIZANJEM
C
C
C      IF(IND.EQ.2) THEN
CC          KONTAKTNA KRUTOST 
C       MDIM=MAXCE+2
C      IF(IILS.NE.-1)THEN
C       KK=-MAXCE-1
C       DO 258 K=1,3
C         DUM=TRA(K,3)
C         SKE(2) =-DUM
C         LLL=MDIM+1
C           DO 255 L=1,MAXCE
C           LLL=LLL+1
C  255      SKE(LLL)= HE(L,1)*DUM
C         KK=KK+MAXCE+2
C         KK1=KK+1
C         LMDUM=LM(KK1)
C         IF(LMDUM.NE.0) LM(KK1)=LMILM
CC         CALL IWRR(LM(KK),MDIM,'LM**')
C         CALL SPAKUJ(A(LSK),A(LMAXA),SKE,LM(KK),MDIM)
C         LM(KK1)=LMDUM
C  258  CONTINUE
C       MDIM=2
C       KK=0
C       DO 260 M=1,3
C       LM3(M)=0
C         IF(IDC(NLM,M).NE.0.AND.IDC(NLM,M).NE.LMILM.AND.KK.LT.2)THEN
C           KK=KK+1
C           LM3(KK)=IDC(NLM,M)
C         ENDIF
C  260  CONTINUE
C       SKE(1) =-1.D0
C       SKE(2) = 0.D0
C       SKE(3) =-1.D0
CC         CALL IWRR(LM3,MDIM,'LM3*')
C       CALL SPAKUJ(A(LSK),A(LMAXA),SKE,LM3,MDIM)
C      ENDIF
CC
CC          KONTAKTNE SILE PRI KLIZANJU
C       IF(RRS.GT.SMALL)THEN
C         DUM=FD*DABS(TAURSN(3,NODMAX))/RRS
C         DO 57 I=1,2
C   57    TAURSN(I,NODMAX)=DUM*TAURSN(I,NODMAX)
C       ENDIF
C       CC=0.D0
C       DO 96 K=1,3
C   96  CC=CC+TRA(K,3)*DC(K)
C       IF(DABS(CC).LT.SMALL)CC=0.D0
C       DC(1)=CC
CC           TRANSFORMACIJA SILA U GLOBALNI SISTEM
C       CALL RSN123(TAURSN(1,NODMAX),TAU123(1,NODMAX),TRA,1)
CC      CALL WRR(TAURSN(1,NODMAX),3,'@RSN ')
CC      CALL WRR(TAU123(1,NODMAX),3,'@123 ')
CCC      CALL WRR(TRA,9,'TRA ')
CC           SILE NA CILJNIM CVOROVIMA
C      KK=-MAXCE
C      DO 65 K=1,3
C       KK=KK+MAXCE+2
C       IF(LM(KK).GT.0) THEN
C        DUM=TAU123(K,NODMAX)
C        SILE(K,NLM)=DUM
C        DO 64 L=1,MAXCE
C        NELL=NE+NELAB(L,NCA)
C        DUM1=SILE(K,NELL)-HE(L,1)*DUM
C        TAU123(K,L)=DUM1
C        SILE(K,NELL)=DUM1
CC        WRITE(3,*)'K,NELL,SILE(K,NELL)',K,NELL,SILE(K,NELL)
CC        WRITE(3,*)'K,L,TAU123(K,L)',K,L,TAU123(K,L)
C   64   CONTINUE
C      ENDIF
C   65 CONTINUE
C      ENDIF
C........... END  IND=2
C
C  3)  RASPOREDJIVANJE KONTAKTNIH SILA
C
      ISLD=1
C      WRITE(3,*)'RC'
C      KK=-MAXCE
      DO 80 K=1,K2
C       KK=KK+MAXCE+2
C       KK=KK+MAXCE+1
       KKLL=LM(K)
C	 RC(KKLL)=0.0
       IF(KKLL.NE.0) THEN
  	  RC(KKLL)=TAU123(K,NODMAX)
C  	  RC(KKLL)=DC(K)
c        WRITE(3,*)'KKLL,RC(KKLL)',KKLL,RC(KKLL)
        DO 70 L=1,MAXCE
          K2L=K+L*K2
          KKL=LM(K2L)
C          RC(KKL)=0.0
          IF(KKL.NE.0) RC(KKL)=TAU123(K,L)
c        WRITE(3,*)'KKL,RC(KKL)',KKL,RC(KKL)
   70   CONTINUE
C           KONTAKTNO PRODIRANJE ne treba za Penalty
C       IF(IND.EQ.2.AND.KKLL.NE.LMILM) GOTO 80
C        DUM=DC(K)
C        IF(IND.EQ.2)DUM=DC(1)
C        RC(KKLL)=DUM
CC        WRITE(3,*)'KKLL,DC',KKLL,DUM
       ENDIF
   80 CONTINUE
C
C..  KONTAKTOR KOD GOVORI O TIPU KONTAKTA I PRVOM CVORU CILJNOG PODSEG.
C
      IK (NLM)=IND*IT1+NCAA
c        WRITE(3,*)'NA KRAJU NCAA',NCAA
C      IK (NLM)=IND*IT1+NCAMIN
	IF(IK(NLM).NE.IK1(NLM)) ICCMOV=ICCMOV+1
      IK1(NLM)=IK(NLM)
C KOORDINATE TACKE PRODORA
         ALFK(1,NLM)=R0
         ALFK(2,NLM)=S0
C PODACI ZA KLIZANJE/SLEPLJIVANJE
      IF (AKS.GT.0) THEN
c         ALFK(3,NLM)=-DC(3)+XYZ(3,NODMAX)
C SILA U TANGENCIJALNOJ RAVNI
         ALFK(3,NLM)=TAU123(1,NODMAX)
         ALFK(4,NLM)=TAU123(2,NODMAX)
         ALFK(5,NLM)=TAU123(3,NODMAX)
C NORMALA U TACKI PRODORA
C         ALFK(6,NLM)=TRA(1,3)
C         ALFK(7,NLM)=TRA(2,3)
C         ALFK(8,NLM)=TRA(3,3)
      ELSE
         ALFK(3,NLM)=TRA(1,3)
         ALFK(4,NLM)=TRA(2,3)
         ALFK(5,NLM)=TRA(3,3)
	ENDIF
C      write(3,*)'status kontakta'
C      write(3,*)'NLM,IK1,R0,S0,ICCMOV',NLM,IK1(NLM),R0,S0,ICCMOV
C      write(3,*)'koord. sistem u tacki prodora'
C      write(3,*)'TRA(1,3),TRA(2,3),TRA(3,3)',TRA(1,3),TRA(2,3),TRA(3,3)
      RETURN
      END
C=======================================================================
      SUBROUTINE SILE93P(IDC,U,IK,ALFK,NE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     KOREKCIJA PRIRASTAJA SILA
C  AKO JE BILO KLIZANJE U PRETHODNOJ ITERACIJI 
C  SILU TRANSFORMISE U GLOBALNI SISTEM KAO DA NEMA KLIZANJA
CE     UPDATE FORCE INCREMENT
C
      DIMENSION U(*),IDC(NE,*),IK(*),ALFK(5,*)
      DATA IT1/1000000/
      DO 100 NLM=1,NE
C      CALL WRR(ALFK(1,NLM),5,'ALFK')
C... SAMO ZA CVOROVE KOJI KLIZAJU
        IND =IK(NLM)/IT1
        IF(IND.NE.2)GO TO 100
        LM=0
        DO 10 K=1,3
   10   IF(LM.EQ.0) LM=IDC(NLM,K)
        DF=U(LM)
        DO 20 K=1,3
          LM=IDC(NLM,K)
C          WRITE(3,*)'*LM',LM
          IF(LM.EQ.0)GO TO 20
          K2=2+K
          U(LM)=DF*ALFK(K2,NLM)
   20 CONTINUE
  100 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE TRAETPP(ETP,ELAST,TSG,ii,jj)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS      TRANSFORMACIJA MATRICE ETP(LOKALNO)  U  ELAST(GLOBALNO)
CS      TRANSPONOVANO TSG  *  ETP  *  TSG  =  ELAST
CE      TRANSFORM  MATRIX  ETP(LOCAL)  TO  ELAST(GLOBAL)
CE      TRANSPOSED TSG  *  ETP  *  TSG  =  ELAST
C
C      DIMENSION ETP(ii,ii),ELAST(jj,jj),TSG(ii,jj),P(15,5)
      DIMENSION ETP(5,5),ELAST(15,15),TSG(5,15),P(15,5)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRAETP'
      DO 80 I=1,jj
      DO 80 J=1,ii
        P(I,J)=0.D0
          DO 75 K=1,ii
   75     P(I,J)=P(I,J)+TSG(K,I)*ETP(K,J)
c      WRITE(3,*)'P(I,J),I,J',P(I,J),I,J
   80 CONTINUE
      DO 84 I=1,jj
      DO 84 J=I,jj
        ELAST(I,J)=0.D0
          DO 82 K=1,ii
   82     ELAST(I,J)=ELAST(I,J)+P(I,K)*TSG(K,J)
        ELAST(J,I)=ELAST(I,J)
C      WRITE(3,*)'ETPG(I,J),I,J',ELAST(I,J),I,J
   84 CONTINUE
      RETURN
      END
C=======================================================================
C=======================================================================
      SUBROUTINE MNOZI(EN,ELAST,ii)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS      TRANSFORMACIJA MATRICE ETP(LOKALNO)  U  ELAST(GLOBALNO)
CS      TRANSPONOVANO TSG  *  ETP  *  TSG  =  ELAST
CE      TRANSFORM  MATRIX  ETP(LOCAL)  TO  ELAST(GLOBAL)
CE      TRANSPOSED TSG  *  ETP  *  TSG  =  ELAST
C
      DIMENSION EN(ii),ELAST(15,15)
C      DIMENSION ETP(5,5),ELAST(15,15),TSG(5,15),P(15,5)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' TRAETP'
      DO 81 I=1,ii
      DO 81 J=1,ii
        ELAST(I,J)=0.0
   81 CONTINUE
      DO 80 I=1,ii
      DO 80 J=1,ii
        ELAST(I,J)=ELAST(I,J)+EN(I)*EN(J)
C      WRITE(3,*)'ELAST(I,J),I,J',ELAST(I,J),I,J
   80 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE MNOZIT(DT,ELAST,AMETRIC,ii)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS      Kt=a(alfa,beta)*Dalfa*DbetaT
CE      Kt=a(alfa,beta)*Dalfa*DbetaT
C
      DIMENSION DT(ii,2),ELAST(ii,ii),AMETRIC(2,2)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' STICK KT'
c        DO 80 I=1,ii
c        DO 80 J=1,ii
c          ELAST(I,J)=0.0
c   80   CONTINUE
        DO 81 I=1,ii
        DO 81 J=1,ii
          ELAST(I,J)=0.0
         DO K=1,2
	    DO M=1,2
           ELAST(I,J)=ELAST(I,J)+AMETRIC(K,M)*DT(I,K)*DT(J,M)
c      WRITE(3,*)'P(I,J),I,J',P(I,J),I,J
	    ENDDO
         ENDDO
   81   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE ALFC3P (HE,TRA,XYZ,UPRI0,MAXCE,DC,R0,S0,EPSIL,PRODOR,
     &            TOLPD,NN,TCMAX,NPL,NPOL,NCTC,NPOLMX,IND,IFLAG,AMETRIC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     NALAZENJE POLOZAJA KONTAKTA
CE     FIND CONTACT LOCATION
C
C  IFLAG = 1   PRETRAZIVANJE PO POVRSINAMA (TEST A,B)
C  IFLAG = 2   PRETRAZIVANJE PO LINIJAMA (TEST C)
C  NN          CVOR POLIGONA NAJBLIZI KONTAKTORU
C  IND         INDIKATOR KONTAKTA (0-NEMA, 1-IMA)
C  DC()        KOMPONENTE PRODORA
C  TRA()       MATRICA TRANSFORMACIJE
C.........  ZA IFLAG=2   ....................................
C  NPOL        REDNI BROJ POLIGONA SA NAJVECIM (TESTC)
C  NPOLMX      UKUPAN BROJ POLIGONA NA SEGMENTU
C  NCTC        CVOR POLIGONA (NPOL) LINIJE NA KOJOJ JE PRODOR
C  TCMAX       MAKSIMALNA PROJEKCIJA
C
      DIMENSION HE(9,*),TRA(3,*),XYZ(3,*),DC(*),UPRI0(*)
      DIMENSION D(3),C1(3),C2(3),P(3),AN(3),C1P(3),PC2(3)
C      DATA RR/1.D0,-1.D0,-1.D0,1.D0/,SS/1.D0,1.D0,-1.D0,-1.D0/
C      DATA W0/0.D0/,W1/1.D0/
C
C      WRITE(3,*)'STIGO 1, IFLAG',IFLAG
      R0=0.D0
      S0=0.D0
      NODMAX=MAXCE+1
      DO 10 I=1,3
       D(I) =XYZ(I,NODMAX)-XYZ(I,NN)
c       WRITE(3,*)'I,NODMAX,NN',I,NODMAX,NN
c       WRITE(3,*)'XYZ(I,NODMAX),XYZ(I,NN)',XYZ(I,NODMAX),XYZ(I,NN)
       DC(I)=0.D0
   10 CONTINUE
C
C--- 2D CONTACT SURFACE
C
      IF(MAXCE.EQ.2)THEN
        MD=2
        IND=0
        DUM=0.D0
        DO 520 I=1,2
          C1(I)=XYZ(I,2)-XYZ(I,1)
          DUM=DUM+C1(I)*C1(I)
  520   CONTINUE        
        DUM=DSQRT(DUM)
C... PROBLEM DEPENDENT TOLERANCE
      TOLPD=EPSIL*DUM
        DO 525 I=1,2
  525   C1(I)=C1(I)/DUM
        R0=(C1(1)*D(1)+C1(2)*D(2))/DUM
C100794
c        WRITE(3,*)'C(1),C(2),D(1),D(2)',C1(1),C1(2),D(1),D(2)
c        WRITE(3,*)'****DUM,NN,R0',DUM,NN,R0
C... POKUSAJ NA DRUGI NACIN
        IF(DABS(R0).GT.1.D08)
     &    R0=(UPRI0(2)*D(1)-UPRI0(1)*D(2))/
     &       (UPRI0(2)*C1(1)-UPRI0(1)*C1(2))/DUM
C100794
        IF(NN.EQ.2) R0=1.+R0
C        WRITE(3,*)'>>>>NN,R0',NN,R0
C  PEGLANJE
        IF(DABS(R0).LT.EPSIL)    R0=0.D0
        IF(DABS(R0-1.).LT.EPSIL) R0=1.D0
C*1        IF(R0.LT.-EPSIL.OR.R0.GT.1.+EPSIL) RETURN
        IF(R0.LT.0.D0.OR.R0.GT.1.D0) RETURN
        IND=1
        CALL TRAHEP(XYZ,HE,TRA,R0,S0,MAXCE,AMETRIC)
C        HE(1,1)=1.-R0
C        HE(2,1)=R0
C        TRA(1,1)= W0
C        TRA(2,1)= W0
C        TRA(3,1)= W1
C        TRA(1,2)= C1(1)
C        TRA(2,2)= C1(2)
C        TRA(3,2)= W0
C        TRA(1,3)=-C1(2)
C        TRA(2,3)= C1(1)
C        TRA(3,3)= W0
CC      CALL WRR(TRA,9,'TRA ')
      ELSE
C
C--- 3D CONTACT SURFACE
C
      MD=3
      IND=1
C*      IF(IFLAG.EQ.2)THEN
C*        R0=RR(NN)+TCMAX*(RR(NCTC)-RR(NN))
C*        S0=SS(NN)+TCMAX*(SS(NCTC)-SS(NN))
C*        IFLAG=1
C*      ELSE
C
C  0) ... PROBLEM DEPENDENT TOLERANCE   (TOLPD)
C
      N1=NN+1
      N2=NN-1
      IF(NN.EQ.MAXCE) N1=1
      IF(NN.EQ.1) N2=MAXCE
      DO 20 I=1,3
       C1(I)=XYZ(I,N1)-XYZ(I,NN)
       C2(I)=XYZ(I,N2)-XYZ(I,NN)
   20 CONTINUE
      CINT2=AOBS(C1,C1)
C
      TOLPD=EPSIL*DSQRT(CINT2)
C
C  1) CVOR KONTAKTOR I CILJNI CVOR SE POKLAPAJU
C
      R0=0.D0
      S0=0.D0
      DI=DSQRT(AOBS(D,D))
C      WRITE(3,*)'DI,EPSIL,IND',DI,EPSIL,IND
      IF(DI.LT.EPSIL)THEN
        IND=-1
C        DO 15 I=1,3
C   15   XYZ(I,NODMAX)=XYZ(I,NN)
        GO TO 50
      ENDIF
C
C  2) CVOR KONTAKTOR I CILJNI CVOR SE NE POKLAPAJU
C
C
C   TEST  (C) :  D * C0    ,    C0 = CK / I CK I
C
      TESTC = AOBS( D, C1 )/CINT2
      IF(TESTC.GT.TCMAX)THEN
        TCMAX=TESTC
        NPOL =NPL
        NCTC =N1
      ENDIF
      CINT2=AOBS(C2,C2)
      TESTC = AOBS( D, C2 )/CINT2
      IF(TESTC.GT.TCMAX)THEN
        TCMAX=TESTC
        NPOL =NPL
        NCTC =N2
      ENDIF
C      WRITE(3,*)'TESTC',TESTC
C
C   NORMALA (MALO N) AN = (C1 X C2) / I C1 X C2 I
C 
      CALL AXBV( C1, C2, AN )
      CALL JEDV( AN(1), AN(2), AN(3) )
C
C   VEKTOR (MALO A)     =  - (AN * D) AN   SMESTA SE U  P
C			     
      AND = -AOBS(AN,D)      
      CALL JEDNAK(P,AN,AND,3)
C
C   VEKTOR (MALO P)  P = D + MALO A
C
      CALL ZBIRM1(P,D,3)
C
C   TEST  (A) :  (C1 X P)*AN > 0   ,   C1P = C1 X P 
C
      CALL AXBV( C1, P, C1P )
      TESTA = AOBS( AN, C1P )
C      WRITE(3,*)'TESTA',TESTA
      IF(TESTA.LT.-TOLPD)THEN
        IND=0
        RETURN
      ENDIF
C
C   TEST  (B) :  (C1 X P)*(P X C2) > 0   ,   PC2 = P X C2
C
      CALL AXBV( P, C2, PC2 )
      TESTB = AOBS( C1P, PC2 )
C      WRITE(3,*)'TESTB',TESTB
      IF(TESTB.LT.-TOLPD)THEN
        IND=0
        RETURN
      ENDIF
C
C   PROVERA DA LI PRECI NA IFLAG=2
C
      IF(DABS(TESTA).LE.TOLPD.OR.DABS(TESTB).LE.TOLPD) IND=0
      IF(IND.EQ.0)THEN
        IND=1
C        R0=RR(NN)+TCMAX*(RR(NCTC)-RR(NN))
C        S0=SS(NN)+TCMAX*(SS(NCTC)-SS(NN))
C      WRITE(3,*)'**R0,S0,IND',R0,S0,IND
C        GO TO 55
        GO TO 50
C*        IF(NPOL.EQ.NPOLMX)IFLAG=2
C*        RETURN
      ENDIF
C
C  3) NALAZENJE TACNIH KOORDINATA (R,S) PRODORA
C
   50 CALL PRODR3(XYZ,R0,S0,MAXCE,IND)
      IF(IND.EQ.-1)THEN
        IND=IABS(IND)
        R0=R0/DABS(R0)
        S0=S0/DABS(S0)
      ENDIF
C      WRITE(3,*)'R0,S0,IND',R0,S0,IND
      IF(IND.EQ.0)RETURN
C   55 CONTINUE
C
C  ... KRAJ ZA IFLAG=1
C
C*      ENDIF
C
C  4) MATRICA TRANSFORMACIJE  TRA() (LOKALNI <-> GLOBALNI SISTEM)
C
      CALL TRAHEP(XYZ,HE,TRA,R0,S0,MAXCE,AMETRIC)
C
C  ... KRAJ ZA MAXCE=2
C
      ENDIF
C
C  5) VEKTOR     (KONTAKTOR_CVOR    CILJNA_PRODORNA_TACKA_P)
C
      IF (IND.EQ.1)THEN
       DO 210 J=1,MD
        DO 200 I=1,MAXCE
  200   DC(J)=DC(J)+HE(I,1)*XYZ(J,I)
        DC(J)=-DC(J)+XYZ(J,NODMAX)
C290694       IF(DABS(DC(J)).LT.1.D-10) DC(J)=0.D0
  210  CONTINUE
c      WRITE(3,*)'DC',(DC(K),K=1,3)
C
C  6) PROVERA DA LI SE RADI O PRODORU ILI ZAZORU
C
      PRODOR = AOBS( DC, TRA(1,3) )
C SNEZA
C       TOLPD=9.9999999998D-12
C SNEZA
      IF(PRODOR.GT.TOLPD) IND=0
C      IF(PRODOR.GT.0.D0) IND=0
C      WRITE(3,*)'93 - IND,PRODOR',IND,PRODOR
      ENDIF
      RETURN
      END      
C=======================================================================
      SUBROUTINE TRAHEP(XYZ,HE,TRA,R,S,MAXCE,AMETRIC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  INTERPOLACIONE FUNKCIJE I MATRICA TRANSFORMACIJE
C
      DIMENSION XYZ(3,*),HE(9,*),TRA(3,*),AMETRIC(2,2)
      DATA C0/0.D0/,C1/1.D0/,C25/.25D0/
C
C 1)  INTERPOLACIJE
C
      IF(MAXCE.EQ.4)THEN
      RP=C1+R
      RM=C1-R
      RK=C1-R*R
      SP=C1+S
      SM=C1-S
      SK=C1-S*S
      HE(1,1)=C25*RP*SP
      HE(2,1)=C25*RM*SP
      HE(3,1)=C25*RM*SM
      HE(4,1)=C25*RP*SM
      HE(1,2)=C25*SP
      HE(2,2)=-C25*SP
      HE(3,2)=-C25*SM
      HE(4,2)=C25*SM
      HE(1,3)=C25*RP
      HE(2,3)=C25*RM
      HE(3,3)=-C25*RM
      HE(4,3)=-C25*RP
C
      ELSEIF(MAXCE.EQ.3)THEN
      HE(1,1)=C1-R-S
      HE(2,1)=R
      HE(3,1)=S
      HE(1,2)=-C1
      HE(2,2)= C1
      HE(3,2)= C0
      HE(1,3)=-C1
      HE(2,3)= C0
      HE(3,3)= C1
      ELSEIF(MAXCE.EQ.2)THEN
        HE(1,1)=C1-R
        HE(2,1)=R
      ENDIF
C
C 2)  BAZNI VEKTORI
C
      IF(MAXCE.GT.2)THEN
      DO 20 I=1,3
      DO 20 J=1,2
      TRA(I,J)=C0
       DO 10 K=1,MAXCE
        TRA(I,J)=TRA(I,J)+HE(K,J+1)*XYZ(I,K)
   10  CONTINUE
   20 CONTINUE
C... NORMALA  V3
        CALL AXBV( TRA(1,1), TRA(1,2), TRA(1,3) )
C        CALL JEDV( TRA(1,1), TRA(2,1), TRA(3,1) )
        CALL JEDV( TRA(1,3), TRA(2,3), TRA(3,3) )
C JEDINICNI VEKTOR TANGENTE, NIJE ORTOGONALNI SISTEM (Sneza)
C        CALL JEDV( TRA(1,2), TRA(2,2), TRA(3,2) )
C... DRUGI JEDINICNI VEKTOR PRAVOUGLOG SISTEMA  V2
C        CALL AXBV( TRA(1,3), TRA(1,1), TRA(1,2) )
C SNEZA
C IZRACUNAVANJE METRICKOG TENZORA
	  DO I=1,2
	    DO J=1,2
	    AMETRIC(I,J)=C0
            DO K=1,3
	        AMETRIC(I,J)=AMETRIC(I,J)+TRA(K,I)*TRA(K,J)
	      ENDDO
	    ENDDO
	  ENDDO
      ELSE
        TRA(1,1)= C0
        TRA(2,1)= C0
        TRA(3,1)= C1
        DO 50 I=1,2
   50   TRA(I,2)=XYZ(I,2)-XYZ(I,1)
        TRA(3,2)= C0
      WRITE(3,*)'tra',TRA(1,2), TRA(2,2)
C SNEZA
        TRA(1,3)=-TRA(2,2)
        TRA(2,3)= TRA(1,2)
        TRA(3,3)= C0
        CALL JEDV( TRA(1,3), TRA(2,3), TRA(3,3) )
C SNEZA
C IZRACUNAVANJE METRICKOG TENZORA
C	  DO I=1,2
C	    DO J=1,2
	    AMETRIC(1,1)=C0
            DO K=1,3
	        AMETRIC(1,1)=AMETRIC(1,1)+TRA(K,2)*TRA(K,2)
	      ENDDO
C	    ENDDO
C	  ENDDO
C        CALL JEDV( TRA(1,2), TRA(2,2), TRA(3,2) )
C
C        TRA(1,3)=-TRA(2,2)
C        TRA(2,3)= TRA(1,2)
C        TRA(3,3)= C0
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE RSN123P(TAURSN,TAU123,TRA,ICAL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C   TRANSFORMACIJA SILA   RSN   NA  123
C   ICAL = 1   -   RSN ---> 123
C   ICAL = 2   -   123 ---> RSN 
C
      DIMENSION TAURSN(*),TAU123(*),TRA(3,*)        
      GO TO (1,2), ICAL
1     DO 20 J=1,3
        CC=0.D0
        DO 15 K=1,3
   15   CC=CC+TRA(J,K)*TAURSN(K)
        TAU123(J)=CC
   20 CONTINUE
      RETURN
C
2     DO 30 J=1,3
        CC=0.D0
        DO 25 K=1,3
   25   CC=CC+TRA(K,J)*TAU123(K)
        TAURSN(J)=CC
   30 CONTINUE
      RETURN
      END
C=======================================================================
