C=======================================================================
C
C   SUBROUTINE UL92GL
C              UL92EK
C              LMMH92
C              TGRF92
C
C=======================================================================
      SUBROUTINE UL92GL
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*250 ACOZ
C
CS     GLAVNI UPRAVLJACKI PROGRAM ZA ULAZNIE PODATAKE 2/D KONTAKTA
CE     MAIN PROGRAM FOR INPUT DATA ABOUT 2/D CONTACT ELEMENTS
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEM99/ MXAU,LAU,LLMEL,LNEL,LNMAT,LAPRS,LIPGC,LIPRC,LISNA,
     1                LBETA
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /DUPLAP/ IDVA
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SRPSKI/ ISRPS
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /RSN   / DETER,IPIVOT,IDETER
      EQUIVALENCE (IETYP,NTSF),(NGAUSX,NTTOT),(LNMAT,LNCVSF),
     &            (LIPRC,LITSRF),(LIPGC,LNELSF),(LLMEL,LIDC),
     &            (LBETA,LFSFD),(LAPRS,LMASE),
     &            (FMSTAT,CPP(1,1)),(FMDIN,CPP(1,2))

C
C... PRIVREMENO RESENJE ZA NEGATIVAN PIVOT, POSLE UBACITI U RESEN
C    PROVERU ZA ICONT=1 I BROJ JEDNACINE VECI OD NEQ
      IPIVOT=1
      ICONT=1
C
CE    BASIC DATA ABOUT ELEMENTS
CS    OSNOVNI PODACI O ELEMENTIMA
C
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) NTSF,NTTOT,FMSTAT,FMDIN
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1010) NTSF,NTTOT,FMSTAT,FMDIN
      IF(NTSF.LE.0) NTSF=1
      IF(NTTOT.LT.2)THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2005)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6005)
      STOP
      ENDIF
      IF(ISRPS.EQ.0.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,2000)NTSF,NTTOT,FMSTAT,FMDIN
      IF(ISRPS.EQ.1.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,6000)NTSF,NTTOT,FMSTAT,FMDIN
C
CS     REPERI U VEKTORU ZA ULAZNE PODATKE
CE     POINTERS IN INPUT VECTOR
C
      NCVE=2
      LMASE=1
      LFSFD=LMASE  +(NE+NTTOT)*IDVA
      IF(NDIN.EQ.0) LFSFD=LMASE
      LNCVSF=LFSFD +NE*2*IDVA
      LITSRF=LNCVSF+NTSF
      LISNA =LITSRF+NTSF
      LNELSF=LISNA +NTSF
      LNEL  =LNELSF+NTTOT+NTSF
      LIDC  =LNEL  +NE*NCVE
      LMXAU =LIDC  +NE*2
      LAU=LMAX
C
      MXAU = LMXAU - 1
      IF(MOD(MXAU,2).NE.0) MXAU=MXAU+1
      LMAX=LAU+MXAU
      IF (LMAX.LE.MTOT) GO TO 5
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2010) LMAX,MTOT
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6010) LMAX,MTOT
      STOP
C
CS     POZIVANJE PROGRAMA ZA ULAZNE PODATKE U VEKTOR AU
CE     CALL ROUTINES FOR INPUT DATA IN VECTOR  AU
C
    5 CALL UL92EK(A(LAU))
      MXAE = MAX0(NTTOT*5*IDVA+8,36*IDVA+8)
      MXAE3= NE+1+(4*NE+1)*IDVA
      IF(IATYP.EQ.3) MXAE3=MXAE3+NE*IDVA
      IF(MXAE3.GT.MXAE) MXAE=MXAE3
      LMAX = LMAX + MXAE
      RETURN
C
 1010 FORMAT(2I5,2F10.0)
C-----------------------------------------------------------------------
 2000 FORMAT(///
     611X,'UKUPAN BROJ CILJNIH 2/D POVRSINA ..............   NTSF =',I5/
     616X,'EQ.0; NTSF = 1'///
     611X,'UKUPAN BROJ CVOROVA NA CILJNIM POVRSINAMA......  NTTOT =',I5/
     7//
     711X,'REFERENTNI STATICKI KOEFICIJENT TRENJA ... FMSTAT =',1PD10.3/
     711X,'REFERENTNI DINAMICKI KOEFICIJENT TRENJA .. FMDIN  =',1PD10.3)
 2005 FORMAT(///' BROJ CVOROVA CILJNE POVRSI MORA BITI    .GE. 2 ')
 2010 FORMAT(///' NEDOVOLJNA DIMENZIJA U VEKTORU A ZA PODATKE O STAPOVIM
     1A'/' POTREBNA DIMENZIJA , LMAX=',I10/' RASPOLOZIVA DIMENZIJA,
     2NTOT =',I10)
C-----------------------------------------------------------------------
 6000 FORMAT(///
     611X,'TOTAL NUMBER OF 2/D TARGET SURFACES ...........   NTSF =',I5/
     616X,'EQ.0; NTSF = 1'///
     611X,'TOTAL NUMBER OF NODES ON TARGET SURFACES ......  NTTOT =',I5/
     7//
     711X,'REFERENT STATIC FRICTION COEFFICIENT ..... FMSTAT =',1PD10.3/
     711X,'REFERENT DYNAMIC FRICTION COEFFICIENT .... FMDIN  =',1PD10.3)
 6005 FORMAT(///' NUMBER OF TARGET SURFACE NODES MUST BE  .GE. 2 ')
 6010 FORMAT(///' NOT ENOUGH SPACE IN WORKING VECTOR  A'
     1/' REQUESTED DIMENSION , LMAX=',I10
     2/' AVAILABLE DIMENSION , MTOT=',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE UL92EK(AU)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     GLAVNI UPRAVLJACKI PROGRAM ZA UCITAVANJE ULAZNIH PODATAKA U AU
CE     MENAGEMENT PROGRAM FOR INPUT DATA IN  AU  VECTOR
C
      include 'paka.inc'
      
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /ELEM99/ MXAU,LAU,LLMEL,LNEL,LNMAT,LAPRS,LIPGC,LIPRC,LISNA,
     1                LBETA
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUPLAP/ IDVA
      COMMON /PLASTI/ LPLAST,LPLAS1,LSIGMA
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /RESTAR/ TSTART,IREST
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      EQUIVALENCE (IETYP,NTSF),(NGAUSX,NTTOT),(LNMAT,LNCVSF),
     &            (LIPRC,LITSRF),(LIPGC,LNELSF),(LLMEL,LIDC),
     &            (LIK,LPLAST),(LIK1,LPLAS1),(LBETA,LFSFD),
     &            (LAPRS,LMASE)
C
      DIMENSION AU(*)
      REAL AU
C
CS     POZIVANJE PROGRAMA ZA ULAZNE PODATKE O 2/D KONTAKT ELEMENTIMA
CE     CALL ROUTINE FOR INPUT DATA ABOUT 2/D CONTACT ELEMENTS
C
      CALL ICLEAR(AU(LNEL),NE*NCVE)
      CALL UL921(AU(LNCVSF),AU(LITSRF),AU(LISNA),AU(LNELSF),AU(LNEL),
     &           AU(LFSFD))
      IF(NBLGR.GE.0) CALL TGRF99(AU(LNCVSF),AU(LITSRF),AU(LNELSF))
C
CS     FORMIRANJE VEKTORA LM
CE     FORM  LM  VECTOR
C
      CALL ICLEAR(AU(LLMEL),NE*NCVE*3)
      CALL LMMH92(A(LID),AU(LNEL),AU(LIDC),AU(LNCVSF),AU(LITSRF),
     &            AU(LNELSF))
C
      LMAX8=LMAX8+1
      WRITE(IELEM,REC=LMAX8)
     1 NTSF,NTTOT,NCVE,MXAU,LNCVSF,LITSRF,
     1 LISNA,LNELSF,LNEL,LLMEL,LFSFD,LMASE
      CALL WRITDD(AU(LMASE),MXAU/IDVA,IELEM,LMAX8,LDUZI)
      IF(MOD(LMAX,2).EQ.0) LMAX=LMAX+1
CS....  PROSTOR ZA VELICINE PRETHODNOG KORAKA
CE....  SPACE FOR PREVIOUS STEP VALUES
      NPROS=(NE*3+1)/2*2/IDVA+NE
      LIK  =LMAX
      LIK1 =LIK +NE
      IF(MOD(LIK1,2).EQ.0) LIK1=LIK1+1
      LMAX =LIK1+2*NE+NE*IDVA
      IF(LMAX.GT.MTOT) CALL ERROR(1)
CS..... PROSTOR ZA VELICINE TEKUCEG KORAKA
      IF(IREST.NE.1) THEN
        CALL WRITDD(A(LIK),NPROS,IELEM,LMAX8,LDUZI)
      ELSE
        CALL READDD(A(LIK),NPROS,IELEM,LMAX8,LDUZI)
      ENDIF
CS      KONTAKTNE SILE
CE      CONTACT FORCES
      LSIGMA=LMAX
      NPROS =(NE+NTTOT)*2
      NPRO2 =2*NPROS
      LMAX  =LSIGMA+NPRO2*IDVA
      IF(IREST.NE.1) THEN
        CALL CLEAR (A(LSIGMA),NPRO2)
        CALL WRITDD(A(LSIGMA),NPRO2,IELEM,LMAX8,LDUZI)
      ELSE
        CALL READDD(A(LSIGMA),NPRO2,IELEM,LMAX8,LDUZI)
      ENDIF
      RETURN
      END
C======================================================================
      SUBROUTINE UL921(NCVSF,ITSRF,ISNA,NELSF,NEL,FSFD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     PODPROGRAM ZA UCITAVANJE PODATAKA O 2-D KONTAKTNIM ELEMENTIMA
CE     SUBROUTINE FOR READING DATA ABOUT 2-D CONTACT ELEMENTS
C
      CHARACTER*250 ACOZ
      COMMON /BROJUK/ KARTIC,INDFOR,NULAZ
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ORIENT/ CPP(3,3),XJJ(3,3),TSG(6,6),BETA,LBET0,IBB0
      COMMON /SRPSKI/ ISRPS
      EQUIVALENCE (IETYP,NTSF),(NGAUSX,NTTOT),(LNMAT,LNCVSF),
     &            (LIPRC,LITSRF),(LIPGC,LNELSF),(LLMEL,LIDC),
     &            (FMSTAT,CPP(1,1)),(FMDIN,CPP(1,2))
C
      DIMENSION NCVSF(*),ITSRF(*),ISNA(*),NELSF(*),NEL(NE,*),FSFD(NE,*)
C
CS     P O D A C I   O   C I L J N I M    2-D    P O V R S I M A
CE     D A T A   A B O U T   T A R G E T   2-D   S U R F A C E S
C
      IF(ISRPS.EQ.0.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,2000) NGE
      IF(ISRPS.EQ.1.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,6000) NGE
C
      LL=1
      DO 10 ITS=1,NTSF
      CALL ISPITA(ACOZ)
      IF(ITS.EQ.1) KARTI=KARTIC
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*)   NN,NCVSF(NN),ITSRF(NN),ISNA(NN),KORC
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1000) NN,NCVSF(NN),ITSRF(NN),ISNA(NN),KORC
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,5001) NN,NCVSF(NN),ITSRF(NN),ISNA(NN),KORC
      IF(ITSRF(NN).EQ.0)THEN
         ITSRF(NN)= LL
      ELSE
         ITSRF(NN)=-LL
      ENDIF
      NC=NCVSF(NN)
      K1=0
    5 K1=K1+14
      IF(NC.LE.K1.OR.KORC.NE.0)THEN
        IF(KORC.EQ.0)THEN
          L1=LL+13-K1+NC
          CALL ISPITA(ACOZ)
          IF(INDFOR.EQ.1)
     1    READ(IULAZ,*)   (NELSF(L),L=LL,L1)
          IF(INDFOR.EQ.2)
     1    READ(ACOZ,1000) (NELSF(L),L=LL,L1)
          LL=L1+1
        ELSE
          READ(IULAZ,*) NELSF(LL)
          LL=LL+1
          DO 7 I=2,NC
          NELSF(LL)=NELSF(LL-1)+KORC
    7     LL=LL+1
        ENDIF
      ELSE
        L1=LL+13
        CALL ISPITA(ACOZ)
        IF(INDFOR.EQ.1)
     1  READ(IULAZ,*)   (NELSF(L),L=LL,L1)
        IF(INDFOR.EQ.2)
     1  READ(ACOZ,1000) (NELSF(L),L=LL,L1)
        LL=L1+1
        GO TO 5
      ENDIF
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,1000) (NELSF(L),L=IABS(ITSRF(NN)),LL-1)
   10 CONTINUE
C
CS     P O D A C I   O   K O N T A K T N I M   E L E M E N T I M A
CE     D A T A   A B O U T   C O N T A C T   E L E M E N T S
C
      IF(ISRPS.EQ.0.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,2010)
      IF(ISRPS.EQ.1.AND.(NULAZ.EQ.1.OR.NULAZ.EQ.3))
     1WRITE(IZLAZ,6010)
C
      DO 20 NN=1,NE
      CALL ISPITA(ACOZ)
      IF(INDFOR.EQ.1)
     1READ(IULAZ,*) (NEL(NN,J),J=1,2),(FSFD(NN,J),J=1,2)
      IF(INDFOR.EQ.2)
     1READ(ACOZ,1020) (NEL(NN,J),J=1,2),(FSFD(NN,J),J=1,2)
      FSFD(NN,1)=FMSTAT*FSFD(NN,1)
      FSFD(NN,2)=FMDIN *FSFD(NN,2)
      IF(NULAZ.EQ.1.OR.NULAZ.EQ.3)
     1WRITE(IZLAZ,1030) (NEL(NN,J),J=1,2),(FSFD(NN,J),J=1,2)
   20 CONTINUE
      RETURN
C
 1000 FORMAT(14I5)
 1020 FORMAT(2I5,2F10.3)
 1030 FORMAT(2I10,8X,D12.5,3X,D12.5)
 5001 FORMAT(11X,I5,6X,I6,3X,I6,4X,I6,4X,I6)
C-----------------------------------------------------------------------
 2000 FORMAT(///,' UCITANI PODACI O 2/D KONTAKT ELEMENTIMA GRUPE ELEMENA
     1TA',I6/1X,61('-')///12X,
     3'SEGMENT',6X,'BROJ',6X,'TIP',5X,'STAMPA',5X,
     4'KORAK'/14X,'BROJ',5X,'CVOROVA',3X,'POVRSI')
 2010 FORMAT(//6X,'KONTAKTNI PAROVI            KOEFICIJENTI TRENJA'
     1       //6X,'CVOR   CILJNI SEGMENT     STATICKI      DINAMICKI')
C-----------------------------------------------------------------------
 6000 FORMAT(///,' UCITANI PODACI O 2/D KONTAKT ELEMENTIMA GRUPE ELEMENA
     1TA',I6/1X,61('-')///12X,
     3'SEGMENT',6X,'BROJ',6X,'TIP',5X,'STAMPA',5X,
     4'KORAK'/14X,'BROJ',5X,'CVOROVA',3X,'POVRSI')
 6010 FORMAT(//6X,'KONTAKTNI PAROVI           KOEFICIJENTI TRENJA'
     1       //6X,'CVOR   CILJNI SEGMENT     STATICKI      DINAMICKI')
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE LMMH92(ID,NEL,IDC,NCVSF,ITSRF,NELSF)
C
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
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
C
      DIMENSION ID(NP,*),NEL(NE,*),IDC(NE,*),NCVSF(*),ITSRF(*),NELSF(*),
     &          LMEL(4)
C  FORMIRANJE DOPUNSKIH JEDNACINA ZA LAGRANZEV MNOZIOC
C  PETLJA PO KONTAKTNIM ELEMENTIMA (CVOROVIMA KONTAKTORA)
      DO 100 NLM=1,NE
       JJ=NEL(NLM,1)
       DO 5 I=1,2
       IF(ID(JJ,I).GT.0)THEN
         JEDN=JEDN+1
         IDUM=JEDN
       ELSE
         IDUM=0
       ENDIF
       IDC(NLM,I)=IDUM
    5  LMEL(I)=IDC(NLM,I)
C  VEZA SA KONTAKTOR TELOM
       DO 7 I=1,2
    7  LMEL(I+2)=ID(JJ,I)
C     FORMIRANJE VISINA STUBOVA
C       WRITE(3,*)'LMEL',LMEL
        ND=4
        CALL VISINE(A(LMHT),ND,LMEL)
      IF (IABS(ICCGG).EQ.1) THEN
         WRITE(IDRAKCE) ND,(LMEL(I),I=1,ND)
         NELUK=NELUK+1
      ENDIF
C
C  PETLJA PO CVROVIMA CILJNOG SEGMENTA
       ISRF=NEL(NLM,2)
       NNC =NCVSF(ISRF)
       LL  =IABS(ITSRF(ISRF))-1
C       WRITE(3,*)'* NNC,LL',NNC,LL
        DO 10 NC=1,NNC
         JJ=NELSF(LL+NC)
         DO 20 I=1,2
   20     LMEL(I+2)=ID(JJ,I)
C     FORMIRANJE VISINA STUBOVA
C       WRITE(3,*)'* LMEL',LMEL
        CALL VISINE(A(LMHT),ND,LMEL)
      IF (IABS(ICCGG).EQ.1) THEN
         WRITE(IDRAKCE) ND,(LMEL(I),I=1,ND)
         NELUK=NELUK+1
      ENDIF
   10   CONTINUE
  100 CONTINUE
C
      NEQC=JEDN-NEQ
C
CC      WRITE(3,*)'***NEQ,NEQC,JEDN',NEQ,NEQC,JEDN
      RETURN
      END
C=======================================================================
      SUBROUTINE TGRF99(NCVSF,ITSRF,NELSF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS     PROGRAM ZA STAMPANJE STAPA  U UNIVERZALNI FILE
CE     PROGRAM FOR PRINTOUT DATA IN UNIVERSAL GRAPHICS FILE
C
      COMMON /ELEIND/ NGAUSX,NGAUSY,NGAUSZ,NCVE,ITERME,MAT,IETYP
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /SUMELE/ ISUMEL,ISUMGR
      EQUIVALENCE (IETYP,NTSF)
      DIMENSION NCVSF(*),ITSRF(*),NELSF(*)
C
C     GRAFICKI OPIS STAPA: SA 2 CVORA = 1
      IFGD=1
C     
      IFDI=14
C     TABELA FIZICKIH OSOBINA
      IPTN=NGE
C     TABELA MATERIJALA
      MPTN=NMODM
C     BOJA  
      ICOL=8
C     BROJ CVOROVA NA ELEMENTU
      NNODS=2
      IND=-1
      ITYP=71
      WRITE(IGRAF,1100) IND
      WRITE(IGRAF,1100) ITYP
      DO 10 ITS=1,NTSF
C     REDNI BROJ ELEMENTA
      IEL=ITS+ISUMEL
C
C  ITS - CILJNI SEGMENT, NNC - BROJ CVOROVA SEGMENTA
      NNC =NCVSF(ITS)
      LL  =IABS(ITSRF(ITS))-1
      DO 8 NC=1,NNC-1
      IC1=NELSF(LL+NC)
      IC2=NELSF(LL+NC+1)
      WRITE(IGRAF,1000) IEL,IFGD,IFDI,IPTN,MPTN,ICOL,NNODS
    8 WRITE(IGRAF,1000) IC1,IC2
   10 CONTINUE
      WRITE(IGRAF,1100) IND
      ISUMEL=ISUMEL+NTSF
      RETURN
C
 1000 FORMAT(8I10)
 1100 FORMAT(I6)
      END
