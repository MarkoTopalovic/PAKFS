C=======================================================================
C
CE       LINEAR ANALYSIS
CS       LINEARNA ANALIZA
C
C   SUBROUTINE FORMGR
C              INTKMK
C              INTKK
C              PERIOD
C              LEVSTR
C              DESNAL
C              UKUPNA
C              NULDIS
C              OPTLIN
C              FORMSS
C              FTDTRC
C              FORMCT
C              DISKCT
C              FTDTCT
C              ZADDES
C              INTNMK
C              DESSTR
C              INTNAP
C              UCITAM
C
C=======================================================================
      SUBROUTINE PODDAT(NPODS,IND)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .                                                                     
CE.    P R O G R A M                                                    
CE.        TO BACK READ BASIC DATA ABOUT SUBSTRUCTURES
CS.    P R O G R A M                                                    
CS.        ZA VRACANJE UCITANIH OPSTIH PODATAKA O PODSTRUKTURAMA        
C .                                                                     
CE.    V E C T O R
CE.      I=1,JPS1   (JPS - TOTAL SUBSTRUCTURE NUMBER)                  
CE.        NPODS(I,*) - BASIC DATA ABOUT SUBSTRUCTURE
CS.    V E K T O R                                                     
CS.      I=1,JPS1   (JPS - BROJ PODSTRUKTURA)                          
CS.        NPODS(I,*) - OPSTI PODACI O PODSTRUKTURAMA I BROJEVI STAZA   
C .                                                                     
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /DIM   / N9,N10,N11,N12,MAXUP
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /OPTERE/ NCF,NPP2,NPP3,NPGR,NPGRI,NPLJ,NTEMP
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /UPRIRI/ LUPRI
      COMMON /GRUPER/ LIGRUP
      COMMON /NELINE/ NGENN
      COMMON /DSTAZE/ NSTAZ
      COMMON /DUPLAP/ IDVA
      DIMENSION NPODS(JPS1,*)
      COMMON /CDEBUG/ IDEBUG
C
      IF(IDEBUG.GT.0) PRINT *, ' PODDAT'
      GO TO (1,2,3,2),IND
C
    1 IF(JPBR.EQ.JPS1) THEN
        NPG   =NPODS(JPBR,3)
        NGA   =NPODS(JPBR,4)
        NGI   =NPODS(JPBR,5)
        JEDNG =NPODS(JPBR,6)
        NCF   =NPODS(JPBR,15)
        NTEMP =NPODS(JPBR,21)
        NZADP =NPODS(JPBR,22)
        NWG   =NPODS(JPBR,24)
        LMAXM =NPODS(JPBR,26)
        NRAD  =NPODS(JPBR,33)
        IDIREK=NPODS(JPBR,69)
        NSTAZ =NPODS(JPBR,83)
        NP=NPG
        NPA=NGA
        NPI=NGI
        JEDN=JEDNG
        NWK=NWG
      ELSE
        NP    =NPODS(JPBR,3)
        NPA   =NPODS(JPBR,4)
        NPI   =NPODS(JPBR,5)
        JEDN  =NPODS(JPBR,6)
        NPK   =NPODS(JPBR,7)
        NGA   =NPODS(JPBR,8)
        NGI   =NPODS(JPBR,9)
        NGELEM=NPODS(JPBR,11)
        ITERG =NPODS(JPBR,13)
        NCF   =NPODS(JPBR,15)
        NPP2  =NPODS(JPBR,16)
        NPP3  =NPODS(JPBR,17)
        NPGR  =NPODS(JPBR,18)
        NPGRI =NPODS(JPBR,19)
        NPLJ  =NPODS(JPBR,20)
        NTEMP =NPODS(JPBR,21)
        NZADP =NPODS(JPBR,22)
        JEDNP =NPODS(JPBR,23)
        NWP   =NPODS(JPBR,24)
        NWK   =NPODS(JPBR,25)
        LMAXM =NPODS(JPBR,26)
        NGEL  =NPODS(JPBR,28)
        NGENL =NPODS(JPBR,29)
        LGEOM =NPODS(JPBR,30)
        NGEOM =NPODS(JPBR,31)
        ITERM =NPODS(JPBR,32)
        NRAD  =NPODS(JPBR,33)
        IDIREK=NPODS(JPBR,69)
        NSTAZ =NPODS(JPBR,83)
      ENDIF
      RETURN
C
    2 NPP=NP
      IF(ITERM.NE.0) THEN
         IF(JPS.GT.1.AND.JPBR.EQ.JPS1) THEN
         ELSE
            LTECV=NPODS(JPBR,40)
            LMAX13=NPODS(JPBR,41)-1
            CALL READDD(A(LTECV),NPP,IPODS,LMAX13,LDUZI)
         ENDIF
      ENDIF
C
      IF(LUL.NE.0) THEN
         IF(JPS.GT.1.AND.JPBR.EQ.JPS1) THEN
         ELSE
            LCORUL=NPODS(JPBR,42)
            LMAX13=NPODS(JPBR,43)-1
            CALL READDD(A(LCORUL),NPP*3,IPODS,LMAX13,LDUZI)
         ENDIF
      ENDIF
C
C      IF(NGENN.NE.0) THEN
         LFTDT=NPODS(JPBR,44)
         LMAX13=NPODS(JPBR,45)-1
         CALL READDD(A(LFTDT),JEDN,IPODS,LMAX13,LDUZI)
C      ENDIF
C
      IF(NGENN.NE.0) THEN
         LUPRI=NPODS(JPBR,46)
         LMAX13=NPODS(JPBR,47)-1
         CALL READDD(A(LUPRI),JEDN,IPODS,LMAX13,LDUZI)
      ENDIF
C
      LRTDT=NPODS(JPBR,48)
      LMAX13=NPODS(JPBR,38)-1
      IF(IND.EQ.4) LMAX13=NPODS(JPBR,39)-1
      CALL READDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
      IF(IND.EQ.4) THEN
         LMAX13=NPODS(JPBR,52)-1
         IF(NGENN.GT.0) LMAX13=NPODS(JPBR,47)-1
         LRTDB=LRTDT+JEDNP*IDVA
         LSKG=LRTDT+JEDN*IDVA
         LSKB=LSKG+JEDNP*IDVA
C         CALL READDD(A(LSKG),JEDN,IPODS,LMAX13,LDUZI)
         NWKP=JEDN-JEDNP
         CALL JEDNA1(A(LRTDB),A(LSKB),NWKP)
      ENDIF
C
      IF(IDIREK.LE.0) THEN
         LDRV0=NPODS(JPBR,66)
         LDRV1=NPODS(JPBR,67)
         LDRV =NPODS(JPBR,68)
         DO 10 I=1,9
            LSTAZZ(I)=NPODS(JPBR,69+I)
   10    CONTINUE
C
C        U C I T A V A N J E   S A   D I S K A
C
         NPROS=3*NP
C        VEKTOR  VN  U TRENUTKU  T=TAU  NA DISKU
         LMAX13=LSTAZZ(1)
         CALL READDD(A(LDRV), NPROS,IPODS,LMAX13,LDUZI)
C        VEKTOR  V1  U TRENUTKU  T=TAU  NA DISKU
         LMAX13=LSTAZZ(2)
         CALL READDD(A(LDRV1),NPROS,IPODS,LMAX13,LDUZI)
         IF(NGENN.NE.0) THEN
C           VEKTOR  V0  U TRENUTKU  T=TAU  NA DISKU
            LMAX13=LSTAZZ(5)
            CALL READDD(A(LDRV0),2*NPROS,IPODS,LMAX13,LDUZI)
         ENDIF
      ENDIF
C
CE    ALLOCATE SPACES FOR ITERATIONS
CS    REZERVISANJE PROSTORA ZA ITERACIJE
C
      IF(JPS.GT.1.AND.JPBR.EQ.JPS1) THEN
         IF(ITERGL.GT.0) CALL MEMITR(A(LIPODS))
      ENDIF
C
      NBLOCK=NPODS(JPBR,49)
      LSK=NPODS(JPBR,50)
      IF(IND.EQ.4) THEN
         LMAX13=NPODS(JPBR,60)-1
         IF(NBLOCK.EQ.1) THEN
            IF(JPBR.EQ.JPS1) THEN
C               CALL READDD(A(LSK),NWG,IPODS,LMAX13,LDUZI)
            ELSE
C               CALL READDD(A(LSK),NWP,IPODS,LMAX13,LDUZI)
               LMAX13=NPODS(JPBR,61)-1
               LSKG=LSK+NWP*IDVA
               NWKP=NWK-NWP
C               CALL READDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
            ENDIF
            LRAD=NPODS(JPBR,51)
         ENDIF
      ELSE
         LMAX13=NPODS(JPBR,35)-1
         IF(NBLOCK.EQ.1) THEN
C            CALL READDD(A(LSK),NWK,IPODS,LMAX13,LDUZI)
            LRAD=NPODS(JPBR,51)
         ENDIF
      ENDIF
      RETURN
C
    3 IF(JPBR.EQ.JPS1) THEN
        LID=LRAD
        JMAXA=LRAD
        LMAX=JMAXA+JEDNG+1
        IF(LMAX.GT.MTOT) CALL ERROR(1)
        LMAX13=NPODS(JPBR,12)-1
        CALL IREADD(A(JMAXA),JEDNG+1,IPODS,LMAX13,LDUZI)
        JEDN=JEDNG
        NWK=NWG
        NP=NPG
        LRAD=LMAX
      ELSE
        LMAX13=NPODS(JPBR,1)-1
        NMM=NGA-NGI+1
        NPP=NP+NPK
        NP=NPP
        LID=LRAD
        LMAX=LID+NPP*6
        NPM=NPA-NPI+1
        LCVEL=LMAX
        LELCV=LCVEL+NPP
        LMAX=LELCV+NPM+NMM
        NP6=NPP*7+NPM+NMM
        CALL DELJIV(LMAX,2,INDL)
        IF(INDL.EQ.0) LMAX=LMAX+1
        LCORD=LMAX
        LMAX=LCORD+NPP*3*IDVA
        IF(LMAX.GT.MTOT) CALL ERROR(1)
        CALL READDD(A(LCORD),NPP*3,IPODS,LMAX13,LDUZI)
        CALL IREADD(A(LID),NP6,IPODS,LMAX13,LDUZI)
        LIGRUP=LMAX
        LMAXA=LIGRUP+NGELEM*5
        LMAX=LMAXA+JEDN+1
        IF(LMAX.GT.MTOT) CALL ERROR(1)
        NP6=NGELEM*5+JEDN+1 
        LMAX13=NPODS(JPBR,12)-1
        CALL IREADD(A(LIGRUP),NP6,IPODS,LMAX13,LDUZI)
        LRAD=LMAX
      ENDIF
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTKMK(NPODS)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        TO CALL PROGRAM TO FORM LINEAR PART OF THE STIFFNESS MATRIX
CS.    P R O G R A M
CS.        ZA POZIVANJE PROGRAMA
CS.        ZA FORMIRANJE KONSTANTNE MATRICE KRUTOSTI I
CS.        ZA ZAPISIVANJE MATRICA NA DISK 4 I 9
C .
C ......................................................................
C
      CHARACTER*250 ACOZ
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GRUPER/ LIGRUP
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /CDEBUG/ IDEBUG
      DIMENSION NPODS(JPS1,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' INTKMK'
      IF(ITEST.EQ.0) GO TO 10
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      CALL ISPITA(ACOZ)
C      CALL READTE(A(LSK),NWK)
C     
CS    WRITE LINEAR MATRIX K ON DISK
CS    ZAPISIVANJE LINEARNE MATRICE K NA DISK
C     MATRICA JE SADA U MODULU
C      CALL WSTAZK(NPODS,LSK,35)
      GO TO  20
C
C
   10 NUL=NWK
      IF(NBLOCK.GT.1) NUL=KC
!      RAD SA BLOKOVIMA JE ZASTAREO     
      CALL CLEARB(ALSK,A(LMAXA),A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,NUL)
CE    FORM LINEAR PART OF THE STIFFNESS MATRIX
      CALL INTKK(A(LIGRUP),NPODS)
20    RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(///' M A T R I C A   S K')
C-----------------------------------------------------------------------
 6000 FORMAT(///' M A T R I X     S K')
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE INTKK(IGRUP,NPODS)
      USE MATRICA
      USE STIFFNESS
      USE DRAKCE8
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        WITH LOOP OVER LINEAR GROUP OF ELEMENTS
CE.        TO CALL ROUTINE TO FORM LINEAR STIFFNESS MATRIX
CS.    P R O G R A M
CS.        SA PETLJOM PO GRUPAMA ELEMENATA
CS.        ZA POZIVANJE PROGRAMA
CS.        ZA FORMIRANJE KONSTANTNE MATRICE KRUTOSTI
CS.
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /ELEALL/ NETIP,NE,IATYP,NMODM,NGE,ISKNP,LMAX8
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /MATERM/ LMODEL,LGUSM
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      COMMON /SCRATC/ ISCRC
      COMMON /SKDISK/ ISKDSK
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IGRUP(NGELEM,*)
C
CE    FORM CONSTANT STIFFNESS MATRIX
CS    RACUNANJE KONSTANTNE MATRICE KRUTOSTI
C
      IF(IDEBUG.GT.0) PRINT *, ' INTKK'
      ISKNP=1
      OPEN (ISCRC,FILE='ZSKLIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND ISCRC
C
      ISKDSK=1
CE    LOOP OVER LINEAR ELEMENT GROUPS
      DO 100 NGE = 1,NGELEM
         IATYP = IGRUP(NGE,3)
         IF(IATYP.GT.0) GO TO 100
         NETIP = IGRUP(NGE,1)
         NE    = IGRUP(NGE,2)
         NMODM = IGRUP(NGE,4)
         LMAX8 = IGRUP(NGE,5)
CE       READING POINTERS FOR MATERIAL MODEL
CS       UCITAVANJE REPERA ZA MODEL MATERIJALA
         CALL UCITAM(A(LMODEL),NMODM)
         LMAX=LRAD
C
CE       CALL ROUTINE FOR INTEGRATION OF THE ELEMENT STIFFNESS MATIX 
CE       AND ADD TO GLOBAL STIFFNESS MATRIX
C
         CALL ELEME(NETIP,2)
C
  100 CONTINUE
      
      IF (TIPTACKANJA.NE.1) THEN
           
      if(.not.allocated(rows)) then
            call sparseassembler_getnz(nonzeros)
            !if (nonzeros.ge.((2**32)-1)) stop 'izmeni imaxa u 64bit'
            allocate(rows(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating rows'
            allocate(iirows(nonzeros),STAT=istat)
            allocate(columns(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating columns'
            allocate(iicolumns(nonzeros),STAT=istat)
            allocate(IMAXA(JEDN+1),STAT=istat)
            if(istat.ne.0) stop 'error allocating IMAXA'
            allocate(ALSK(nonzeros),STAT=istat)
            if(istat.ne.0) stop 'error allocating ALSK'
             ALLOCATE (AIROWS(nonzeros*2), STAT = istat)
         if(istat.ne.0) stop 'error allocating AIROWS'
              IF(NDIN.GT.0.OR.ISOPS.GT.0) THEN 
                IF(IMASS.GE.1) THEN
              ALLOCATE (ALSM(nonzeros), STAT = iAllocateStatus)
          IF (iAllocateStatus /= 0) write(3,*)'ALSM Not enough memory *'
          IF (iAllocateStatus /= 0) STOP '*** ALSM Not enough memory *'
                ENDIF
             ENDIF
             IF(NDIN.GT.0) THEN
                IF(IDAMP.GT.0) THEN
                   ALLOCATE (ALSC(nonzeros), STAT = iAllocateStatus)
          IF (iAllocateStatus /= 0) write(3,*)'ALSC Not enough memory *'
          IF (iAllocateStatus /= 0) STOP '*** ALSC Not enough memory *'
                ENDIF
             ENDIF
          endif

      CALL sparseassembler_getsparse(nonzeros,rows,columns,ALSK,IMAXA)
      IMAXA(JEDN+1) = nonzeros+1
      do i=1,nonzeros
            iirows(i)= rows(i)
            iicolumns(i)=columns(i)
            AIROWS(i)= rows(i)
            AIROWS(nonzeros+i)=columns(i)
      enddo
   
          CALL sparseassembler_kill()
      ENDIF
      
      
      
      ISKDSK=0
      IF(NBLOCK.GT.1)THEN
        LLM =LRAD
        LSKE=LLM+100
        CALL SPAKUA(ALSK,A(LMAXA),A(LSKE),A(LLM),ND,0,
     &              A(LMNQ),A(LLREC),NBLOCK,LR,IBLK,A(LCMPC),A(LMPC))
      ENDIF
      CLOSE (ISCRC,STATUS='KEEP')
C
CS    WRITE LINEAR MATRIX K ON DISK
CS    ZAPISIVANJE LINEARNE MATRICE K NA DISK
C     SADA SE ZAPISUJE U MODUL
C      CALL WSTAZK(NPODS,LSK,35)
      if(jedn.le.30) CALL WRR6(ALSK,NWK,'K07W')
      RETURN
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE PPAKS(NPODS,KOCID,IPDT,DT0,VREM0,KORB0,KOJPAK)
      USE MATRICA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        FOR INTEGRATION OF SYSTEM VECTORS AND MATRICES, SOLUTION OF
CE.        THE SYSTEM, EQUILIBRIUM ITERATIONS AND RESULTS PRINTING
CS.    P R O G R A M
CS.        SA PETLJOM PO VREMENSKIM PERIODIMA I KORACIMA
C .
CE.    V A R I A B L E S
CE.        METOD - ITERATION METHOD EMPLOYED, SEE CARD /9/
CE.        KOR   - TIME STEP NUMBER
CE.        ITER  - NUMBER OF ITERATION
CE.        IREST - INDICATOR FOR RESTART, /7/
CE.        NGEL  - NUMBER OF LINEAR GROUP OF ELEMENTS
CE.        NGENN - NUMBER OF NONLINEAR GROUP OF ELEMENTS
CE.        NDIN  - TYPE OF ANALYSIS (=0, STATICS; =1, DYNAMICS), /4/
C .
C ......................................................................
C
      include 'paka.inc'
      
      CHARACTER*6 IME,IM1,IMD,IM2
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /ITERBR/ ITER
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /OPTERE/ NCF,NPP2,NPP3,NPGR,NPGRI,NPLJ,NTEMP
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /IMPERF/ NMODS,LIDIM,LSCIM,MODES
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /UPRIRI/ LUPRI
      COMMON /RESTAR/ TSTART,IREST
      COMMON /NEBACK/ INDBACK
      COMMON /NELINE/ NGENN
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /GRUPER/ LIGRUP
      COMMON /LINSBR/ LIN,KORBR
      COMMON /INDKON/ IKONV
      COMMON /INDNAP/ NAPON
      COMMON /RADIZA/ INDBG
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /ZAPSIL/ UBXYZ(3,10),NVFUB(10),INDZS,LZAPS,LNPRZ
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /PRIT2D/ LFAKP,LTHICV,LITIPE,LNFUN,LIPRAV,LNODPR,ISRBA
      COMMON /JCERNI/ LKOLKO,ICERNI
      COMMON /STOSZP/ STOSZ
      COMMON /GEORGE/ TOLG,ALFAG,ICCGG
      COMMON /DRAKCE/ IDRAKCE,NELUK,NZERO,NEED1,NEED2,NEED3,NNZERO
     1                ,IROWS,LAILU,LUCG,LVCG,LWCG,LPCG,LRCG
      COMMON /FORSPA/ Lr1,Lz1,Lp1,LM0,LMaxM1,LColM
      COMMON /CDEBUG/ IDEBUG
      INCLUDE 'mpif.h'
      integer myid, ierr
      DIMENSION NPODS(JPS1,*)
C
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      IF (myid.ne.0) goto 43
      IF(IDEBUG.GT.0) PRINT *, ' PPAKS'
C
      DT=DT0
      VREME=VREM0
      KOR=KORB0
      KORBR=KOR
C
      ITER=0
      SMKOR=0.
      NSMKOR=0
      NGRKOR=0
      IF(ITERGL.GT.0.OR.KOR.EQ.NDT) KOCID=1
C
      IF(ISTKO.GT.0.AND.NBLPR.GE.0) THEN
         CLOSE (IZLAZ,STATUS='KEEP')
         IME(1:2)='LS'
         WRITE(IME(3:6),1900) KOR
         DO 6 J=1,6
            IF(IME(J:J).EQ.' ') IME(J:J)='0' 
    6    CONTINUE
         OPEN (IZLAZ,FILE=IME,STATUS='UNKNOWN',FORM='FORMATTED',
     1                          ACCESS='SEQUENTIAL')
         CLOSE (IZLAZ,STATUS='DELETE')
         OPEN (IZLAZ,FILE=IME,STATUS='UNKNOWN',FORM='FORMATTED',
     1                          ACCESS='SEQUENTIAL')
         REWIND IZLAZ
      ENDIF
      IF(ISTKO.GT.0.AND.NBLGR.GE.0) THEN
         CLOSE (IGRAF,STATUS='KEEP')
         IM1(1:2)='UN'
         WRITE(IM1(3:6),1900) KOR
         DO 7 J=1,6
            IF(IM1(J:J).EQ.' ') IM1(J:J)='0' 
    7    CONTINUE
         OPEN (IGRAF,FILE=IM1,STATUS='UNKNOWN',FORM='FORMATTED',
     1                          ACCESS='SEQUENTIAL')
         CLOSE (IGRAF,STATUS='DELETE')
         OPEN (IGRAF,FILE=IM1,STATUS='UNKNOWN',FORM='FORMATTED',
     1                          ACCESS='SEQUENTIAL')
         REWIND IGRAF
      ENDIF
C     rakic 29.05.2012.
      IF(ISTKO.GT.0.AND.NBLGR.GE.0) THEN
         CLOSE (49,STATUS='KEEP')
         IM2(1:2)='NE'
         WRITE(IM2(3:6),1900) KOR
         DO 8 J=1,6
            IF(IM2(J:J).EQ.' ') IM2(J:J)='0' 
    8    CONTINUE
         OPEN (49,FILE=IM2,STATUS='UNKNOWN',FORM='FORMATTED',
     1                          ACCESS='SEQUENTIAL')
         CLOSE (49,STATUS='DELETE')
         OPEN (49,FILE=IM2,STATUS='UNKNOWN',FORM='FORMATTED',
     1                          ACCESS='SEQUENTIAL')
         REWIND 49
      ENDIF
 1900 FORMAT(I4)
C
      IF(JPS.EQ.1) JPBR=JPS1
      LIN=0
      IKONV=0
      NAPON=0
      IF(INDBG.LE.0) THEN
         IF(ISRPS.EQ.0)
     1WRITE(*,2000) KOR,NDT,VREME
         IF(ISRPS.EQ.1)
     1WRITE(*,6000) KOR,NDT,VREME
      ENDIF
C
CZILESK
C     POMERANJA IZ PREDHODNOG KORAKA
      CALL WSTAZ(NPODS,LRTDT,87)
C      CALL WRR6(A(LRTDT),JEDN,'POMP')
      IF(KOR.GT.1.AND.NDIN.EQ.0.AND.ISOPS.GT.0) THEN
C     POMERANJA IZ PREDHODNOG KORAKA
C        CALL WSTAZ(NPODS,LRTDT,87)
C     NORMALE IZ PREDHODNOG KORAKA
         IF(IDIREK.LE.0) THEN
            LMAX13=NPODS(JPBR,88)-1
            NPROS=3*NP
            CALL WRITDD(A(LDRV),2*NPROS,IPODS,LMAX13,LDUZI)
         ENDIF
      ENDIF
CZILESK
C
      IF(INDZS.GT.0) THEN
         CALL CLEAR(A(LZAPS),JEDN)
         CALL ICLEAR(A(LNPRZ),JEDN)
      ENDIF
C
      IF(ICCGG.EQ.1) CALL ICLEAR(A(Lz1),LAILU-Lz1)
      IF(ICCGG.EQ.-1) CALL ICLEAR(A(LAILU),Lz1-LAILU)
C
CE    C A S E   O F   A U T O M A T I C   L O A D I N G
CS    S L U C A J   A U T O M A T S K O G   O P T E R E C I V A N J A
C
C
43    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(METOD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(KOR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF(METOD.GT.5) THEN
         IF (myid.ne.0) then
            IF(ISOPS.GT.0) THEN
               IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2100) METOD
               IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6100) METOD
              STOP 'PROGRAM STOP - PAK07 - PERIOD'
            ENDIF
         ENDIF
         IF(KOR.EQ.1) THEN

            IF (myid.ne.0) then
C
CE             INTEGRATION OF LINEAR PART OF THE STIFFNESS MATRIX K
CS             INTEGRACIJE KONSTANTNE MATRICE K PO GRUPAMA ELEMENATA
C
               IF(NGEL.NE.0) CALL INTKMK(NPODS)
C
CE             FORM MASS AND DAMP MATRIX (M AND C)
CS             PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE MATRICE M I C
C
               IF(NDIN.NE.0) CALL INTKMM
            endif

            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(NDIN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(ISOPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF(NDIN.GT.0.AND.ISOPS.GT.0) RETURN
            IF (myid.ne.0) goto 44
C
CE          INITIALIZE VECTOR OF RIGHT SIDE
CS          INICIJALIZACIJA VEKTORA DESNE STRANE
C
            CALL CLEAR(A(LRTDT),JEDN)
            IF(ICONT.NE.0) CALL CLEAR(A(LRCTDT),JEDN)
C           RESTART
            IF(IREST.EQ.1) CALL BACKUP(2)
            IF(NDIN.EQ.0) THEN
CE             IMPOSE INITIAL IMPERFECTION
CS             INICIJALNA IMPERFEKCIJA - UPOTREBA SOPSTVENIH VEKTORA
               IF(NMODS.GT.0.AND.IREST.EQ.0) THEN
                  CALL IPRF(A(LRTDT),A(LFTDT),A(LID),A(LIDIM),
     1                      A(LSCIM),A(LELCV),NPI,ICVEL)
CE                UPDATE OF COORDINATES FOR IMPERFECTION
CS                KOREKCIJA KOORDINATA ZA IMPERFEKCIJU
                  CALL ULCOR(A(LCORD),A(LCORD),A(LRTDT),A(LID),NP)
C                 CALL WRR3(A(LCORD),NP*3,'CORD')
C                 CALL WRR6(A(LRTDT),JEDN,'RTDT')
C                 CALL WRR6(A(LFTDT),JEDN,'FTDT')
               ENDIF
            ENDIF
C
CE             EQUATION NUMBER FOR AUTOMATIC LOADING 
CS             BROJ JEDNACINE SA KONTROLISANIM POMERANJEM
C   
            IF(IAUTO.NE.0) CALL RBJED(A(LID),A(LELCV),NPI,ICVEL)
C
         ENDIF
         IF (myid.ne.0) goto 44
C
         IF(IREST.EQ.1.AND.METOD.EQ.-1.AND.
     +                  KOR.EQ.1.AND.ITER.EQ.0.AND.NBLOCK.EQ.1) THEN
C            CALL READD(A(LSK),NWK,ILDLT)
C           CALL WRR(A(LSK),NWK,'RLDL')
         ENDIF
         IF(NDIN.EQ.0.AND.NGENL.GT.0) GO TO 20
C
C         IF(KOR.EQ.IPDT) THEN
44       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(KOR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(IPDT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(NGENL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF(KOR.EQ.IPDT.OR.NGENL.GT.0) THEN
            IF (myid.ne.0) goto 10
            IF(NDIN.EQ.0) GO TO 10
            IF(KOR.GT.1.AND.NGENL.GT.0) GO TO 10
C
CE          FORM DYNAMIC CONSTANTS FOR CONSTANT TIME STEP
CS          FORMIRANJE DINAM. KONSTANTI ZA KONSTANTNI VREMENSKI KORAK
C
            CALL CONSTD
C
CE          FORM INITIAL VECTOR OF DISPLACEMENTS INCREMENT AT CENTRAL
CE          DIFFERENCE
CS          FORMIRANJE POCETNOG VEKTORA PRIRASTAJA POMERANJA KOD 
CS          CENTALNIH RAZLIKA
C
            IF(MDVI.EQ.3) CALL CENPOC
C
CE          FORM MATRIX OF LEFT SIDE: KE = K + A0*M +A1*C
CS          FORMIRANJE MATRICE LEVE STRANE: KE = K + A0*M +A1*C
C
   10       CALL LEVSTR(NPODS)
C
   20    ENDIF
C
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(IREST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(VREME,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(TSTART,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF(IREST.EQ.1.AND.VREME.LE.TSTART) RETURN
C
CE       LOOP OVER GROUPS OF ELEMENT TO INTEGRATION OF NONLIN. MATRIX K
CS       PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE NELIN. MATRICE K
C
CSKDISK...
C         IF(NGENL.NE.0.AND.NGEL.NE.0) CALL INTKK(A(LIGRUP),NPODS)
CSKDISK...
         CALL INTNMK(A(LIGRUP),NPODS)
         IF (myid.ne.0) goto 46
C
CE       FORM CONSTANT PART OF LOAD VECTOR
CS       FORMIRANJE KONSTANTNOG DELA VEKTORA OPTERECENJA
C
         IF(KOR.EQ.1) CALL DESNAL(NPODS,KOJPAK)
C
CE       ITERATION METODS
CS       ITERATIVNI METODI
C
46       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(ITERGL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF(ITERGL.GT.0) THEN
            IF (myid.ne.0) goto 47
C CERNI
            IF(ICERNI.EQ.1) THEN
C              ZAPISIVANJE POMERANJA U ITERACIJAMA
               IDUM=27
               IMD='ZDIVER'
               III=KOR
               KOR=ITER+1
               VRE=VREME
               VREME=KOR
               STOSZ=VREME
               OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',FORM='FORMATTED',
     1               ACCESS='SEQUENTIAL')
               CALL NEUTRC(MAXIT,IDUM,1)
C
               IF(ICERNI.EQ.1) CALL ICLEAR(A(LKOLKO),NP)
C               CALL STAGP1(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,IDUM,0,
C     +                     A(LNCVP),NCVPR)
               CALL STAU09(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,IDUM,1,
     +                     A(LNCVP),NCVPR)
               KOR=III
               VREME=VRE
            ENDIF 
C CERNI
47          CALL IMETOD
            IF (myid.ne.0) goto 48
C CERNI
            IF(ICERNI.EQ.1) THEN
               IDUM=27
               CALL NEUTRC(MAXIT,IDUM,2)
               CLOSE(IDUM,STATUS='KEEP')
            ENDIF   
C CERNI
            CALL RSTAZ(NPODS,LRTDT,52)
            NAPON=1
            CALL INTNAP(A(LIGRUP))
48       ENDIF
C
CE       R E S U L T S
CS       R E Z U L T A T I
C
CE       PRINT DISPLACEMENTS AND STRESSES
CS       STAMPANJE POMERANJA I NAPONA
C
         if (myid.eq.0) CALL STAMPA(A(LIPODS))
C
      ELSE
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(JPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
C
  888    DO 200 JPBB=1,JPS
C
            IF (myid.eq.0) then
               IF(JPS.GT.1) THEN
C                 LRAD=JMAXA
                  JPBR=JPBB 
                  CALL PODDAT(NPODS,1)
                  CALL PODDAT(NPODS,3)
                  CALL PODDAT(NPODS,2)
               ELSE
                  JPBR=JPS1
               ENDIF
            endif
C
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(KOR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(ITER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF(KOR.EQ.1.AND.ITER.EQ.0) THEN
               if (myid.ne.0) goto 55
C
CE             INTEGRATION OF LINEAR PART OF THE STIFFNESS MATRIX K
CS             INTEGRACIJE KONSTANTNE MATRICE K PO GRUPAMA ELEMENATA
C
               IF((IREST.NE.1.OR.METOD.NE.-1).AND.NGEL.NE.0) 
     +         CALL INTKMK(NPODS)
C              CALL WRR(A(LSK),NWK,'SK  ')
C
CE             FORM MASS MATRIX M 
CS             PETLJA PO GRUPAMA ELEM. RADI  INTEGRACIJE MATRICE M I C
C
               IF(NDIN.NE.0) CALL INTKMM
55             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(NDIN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
               CALL MPI_BCAST(ISOPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
               IF(NDIN.GT.0.AND.ISOPS.GT.0) THEN
                  IF(JPS.GT.1.and.myid.eq.0) LRAD=LID
                  GO TO 200
               ENDIF
               if (myid.ne.0) goto 56
C
CE             INITIALIZE VECTOR OF RIGHT SIDE
CS             INICIJALIZACIJA VEKTORA DESNE STRANE
C
               CALL CLEAR(A(LRTDT),JEDN)
               IF(ICONT.NE.0) CALL CLEAR(A(LRCTDT),JEDN)
C              RESTART
               IF(IREST.EQ.1) CALL BACKUP(2)
               IF(NDIN.EQ.0) THEN
CE                IMPOSE INITIAL IMPERFECTION
CS                INICIJALNA IMPERFEKCIJA - UPOTREBA SOPSTVENIH VEKTORA
                  IF(NMODS.GT.0.AND.IREST.EQ.0) THEN
                     CALL IPRF(A(LRTDT),A(LFTDT),A(LID),A(LIDIM),
     1                         A(LSCIM),A(LELCV),NPI,ICVEL)
CE                   UPDATE OF COORDINATES FOR IMPERFECTION
CS                   KOREKCIJA KOORDINATA ZA IMPERFEKCIJU
                     CALL ULCOR(A(LCORD),A(LCORD),A(LRTDT),A(LID),NP)
                  ENDIF
               ENDIF
56          ENDIF
            if (myid.ne.0) goto 49
C           CALL WRR(A(LRTDT),JEDN,'RKO1')
C
            IF(IREST.EQ.1.AND.METOD.EQ.-1.AND.
     1          KOR.EQ.1.AND.ITER.EQ.0.AND.NBLOCK.EQ.1) THEN
C               CALL READD(A(LSK),NWK,ILDLT)
C               CALL WRR(A(LSK),NWK,'RLDL')
            ENDIF
C            IF(IREST.EQ.1.AND.VREME.LE.TSTART) RETURN
C            CALL WRR6(A(LSK),NWK,'SK00')
49          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(NDIN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(NGENL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(KOR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(IPDT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF(NDIN.EQ.0.AND.NGENL.GT.0) GO TO 40
C
C            IF(KOR.EQ.IPDT) THEN
            IF(KOR.EQ.IPDT.OR.NGENL.GT.0) THEN
               IF (myid.ne.0) goto 30
               IF(ITER.GT.0) STOP 'ITER > 0'
               IF(NDIN.EQ.0) GO TO 30
C
CE             FORM DYNAMIC CONSTANTS FOR CONSTANT TIME STEP
CS             FORMIRANJE DINAM. KONSTANTI ZA KONSTANTNI VREM. KORAK
C
               IF(KOR.EQ.IPDT) CALL CONSTD
               IF(KOR.GT.1.AND.NGENL.GT.0) GO TO 30
C
CE             FORM INITIAL VECTOR OF DISPLACEMENTS INCREMENT AT CENTRAL
CE             DIFFERENCE
CS             FORMIRANJE POCETNOG VEKTORA PRIRASTAJA POMERANJA KOD 
CS             CENTALNIH RAZLIKA
C
               IF(MDVI.EQ.3) CALL CENPOC
C
CE             FORM EFFECTIVE LINEAR LEFT-HAND-SIDE STIFFNESS MATRIX
CE             KE = K + A0*M +A1*C
CE             FACTORIZATION OF THE SYSTEM MATRIX IN CASE OF LINEAR
CE             ANALYSIS
CS             FORMIRANJE MATRICE LEVE STRANE: KE = K + A0*M +A1*C
C
 30            CALL LEVSTR(NPODS)
C
   40       ENDIF
C
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(IREST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(VREME,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(TSTART,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF(IREST.EQ.1.AND.VREME.LE.TSTART) RETURN
            IF (myid.ne.0) goto 50
C
CE          FORM CONSTANT PART OF LOAD (RIGHT-HAND-SIDE) VECTOR
CS          FORMIRANJE KONSTANTNOG DELA VEKTORA OPTERECENJA
C
            CALL DESNAL(NPODS,KOJPAK)
C           CALL WRR(A(LRTDT),JEDN,'RDEL')
C
CE          LOOP OVER NONLINEAR GROUP OF ELEMENTS TO INTEGRATE NONLINEAR
CE          STIFFNESS MATRIX KNL
CS          PETLJA PO GRUPAMA ELEM. RADI INTEGRACIJE NELIN. MATRICE K
C
CSKDISK...
C         IF(NGENL.NE.0.AND.NGEL.NE.0) CALL INTKK(A(LIGRUP),NPODS)
CSKDISK...
C            CALL WRR6(A(LSK),NWK,'SK01')
50          CALL INTNMK(A(LIGRUP),NPODS)
            IF (myid.ne.0) goto 51
C
CE          FORM VECTOR OF RIGHT-HAND-SIDE
CE          SUBTRACT INTERNAL FROM EXTERNAL FORCES VECTOR,RTDT=RTDT-FTDT
CS          FORMIRANJE KONACNOG VEKTORA DESNE STRANE:
CS          NA UCITANE KONSTANTNE SILE DODAJU SE GEOMETRIJSKI
CS          NELINEARNE I ODUZMU SE UNUTRASNJE SILE, RTDT=RTDT-FTDT
C
            CALL DESSTR(NPODS)
C           CALL WRR(A(LRTDT),JEDN,'RDSS')
C
CE          PUT ON RIGHT SIDE PRODUCT OF LARGE NUMBER AND PRESCRIBED
CE          DISPLACEMENT
CS          POSTAVLJANJE NA DESNOJ STRANI PROIZVODA VELIKOG BROJA I
CS          ZADATOG POMERANJA
C
            IF(NZADP.GT.0) CALL ZADATD(NPODS,1)
C
51          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(JPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF(JPS.GT.1) THEN
               CALL RESEN(ALSK,A(LRTDT),A(LMAXA),JEDN,2)
               IF (myid.eq.0) then
                  LMAX13=NPODS(JPBR,39)-1
                  CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
                  LRAD=LID
               endif
            ENDIF
C
  200    CONTINUE
C         CALL WRR(A(LRTDT),JEDN,'R200')
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(NDIN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(ISOPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(JSP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF(NDIN.GT.0.AND.ISOPS.GT.0) RETURN
C
         IF(JPS.GT.1) THEN
            IF (myid.ne.0) goto 52
            DO 300 JPBR=1,JPS
               CALL PODDAT(NPODS,1)
               LSK=LRAD
               LSKG=LSK+(NWK-NWP)*IDVA
               JED=JEDN-JEDNP
               NWKP=JED*(JED+1)/2
               LIGRUP=LSKG+NWKP*IDVA
               LMAXA=LIGRUP+NGELEM*5
               LMAX=LMAXA+JEDN+1
               IF(LMAX.GT.MTOT) CALL ERROR(1)
               NP6=NGELEM*5+JEDN+1 
               LMAX13=NPODS(JPBR,12)-1
               CALL IREADD(A(LIGRUP),NP6,IPODS,LMAX13,LDUZI)
               LMAX13=NPODS(JPBR,61)-1
!                 SADA JE SACUVANO U MODULU               
C               CALL READDD(A(LSK),NWK-NWP,IPODS,LMAX13,LDUZI)
C               CALL TROUGO(A(LSK),A(LSKG),A(LMAXA+JEDNP),JED,NWP)
              LMAX13=NPODS(JPBR,37)-1
C              CALL WRITDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
  300       CONTINUE
C
            JPBR=JPS1
            CALL PODDAT(NPODS,1)
            CALL PODDAT(NPODS,3)
            CALL PODDAT(NPODS,2)
C
CE          SET DISK ON ZERO
CS          POSTAVLJANJE DISKA  NA NULU
C
            CALL NULDIS(NPODS)
            CALL CLEAR(ALSK,NWG)
C
CE          FORM VECTOR OF CONCENTRATED FORCES AT TIME T+DT
CS          FORMIRANJE VEKTORA KONCENTRISANIH SILA U TRENUTKU T+DT
C
            IF(NCF.NE.0) THEN 
               LMAX13=NPODS(JPBR,34)-1
               CALL FORMSS
               CALL WSTAZ(NPODS,LRTDT,38)
               IF(NBLGR.GE.0) THEN
                  LID=LRAD
                  LCVEL=LID+NP*6
                  LMAX=LCVEL+NP
                  NP6=NP*7
                  LMAX13=NPODS(JPBR,2)-1
                  CALL IREADD(A(LID),NP6,IPODS,LMAX13,LDUZI)
                  CALL STAGKS(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP
     1                           ,IGRAF,A(LNCVP),NCVPR)
                  CALL STAU09(A(LRTDT),A(LID),A(LCVEL),ICVEL,
     1                           NP,49,41,A(LNCVP),NCVPR)
               ENDIF
            ENDIF
C
            DO 400 JPBR=1,JPS
               JEDN=NPODS(JPBR,6)
               JEDNP=NPODS(JPBR,23)
               JED=JEDN-JEDNP
               NWKP=JED*(JED+1)/2
               LSKG=LRAD
               LRTG=LSKG+NWKP*IDVA
               LLMG=LRTG+JEDN*IDVA
               LMAX=LLMG+JED
               IF(LMAX.GT.MTOT) CALL ERROR(1)
               LMAX13=NPODS(JPBR,37)-1
               CALL READDD(A(LSKG),NWKP,IPODS,LMAX13,LDUZI)
               LMAX13=NPODS(JPBR,39)-1
               CALL READDD(A(LRTG),JEDN,IPODS,LMAX13,LDUZI)
               LMAX13=NPODS(JPBR,27)-1
               CALL IREADD(A(LLMG),JED,IPODS,LMAX13,LDUZI)
               CALL SPAKUJ(ALSK,A(JMAXA),A(LSKG),A(LLMG),JED)
               LRTG=LRTG+JEDNP*IDVA
               CALL SPAKUD(A(LRTDT),A(LRTG),A(LLMG),JED)
  400       CONTINUE
            JPBR=JPS1
            JEDN=JEDNG
            IF(NGENN.GT.0) CALL JEDNA1(A(LFTDT),A(LRTDT),JEDN)
C            CALL WSTAZK(NPODS,LSK,35)
52          CALL RESEN(ALSK,A(LRTDT),A(JMAXA),JEDN,1)
            IF (myid.ne.0) goto 53
            CALL WSTAZK(NPODS,LSK,60)
            LMAXA=JMAXA
53       ENDIF
         IF (myid.ne.0) goto 54
C
CS       RACUNANJE UKUPNIH OPTERECENJA
C
         IF(KOR.EQ.1.AND.ITER.EQ.0)
     1   CALL SILMOM(A(LID),A(LRTDT),NP,IZLAZ,ISRPS)
C
CE       BACKSUBSTITUTION - SOLUTION OF SYSTEM EQUATIONS
CS       ZAMENA UNAZAD - RESAVANJE :SISTEMA JEDNACINA
C
C
C        CALL WRR6(A(LSK),NWK,'SR-1')
C        CALL WRR6(A(LRTDT),JEDN,'RTDT')
54       CALL RESEN(ALSK,A(LRTDT),A(LMAXA),JEDN,2)
C        CALL WRR6(A(LRTDT),JEDN,'RESE')
C
C
CE       N O N L I N E A R      A N A L Y S I S
CS       N E L I N E A R N A      A N A L I Z A
C
C
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(NGENN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         IF(NGENN.GT.0) THEN

C
CE          COMPUTE OF ENERGY ON THE BEGINING OF STEP: E1=U(0)*(RTDT-FT)
CS          RACUNANJE ENERGIJE NA POCETKU KORAKA : ENE1=U(0)*(RTDT-FT)
C
C CERNI
            IF (myid.eq.0) THEN
               IF(ICERNI.EQ.1) THEN
C                 ZAPISIVANJE POMERANJA U ITERACIJAMA
                  IDUM=27
                  IMD='ZDIVER'
                  OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',FORM='FORMATTED',
     1               ACCESS='SEQUENTIAL')
                  CALL NEUTRC(MAXIT,IDUM,1)
C
                  IF(ICERNI.EQ.1) CALL ICLEAR(A(LKOLKO),NP)
               ENDIF   
C CERNI
               IF(ITERGL.GT.0.AND.KONVE.GT.0) CALL KONVER(IKONV)
C
CE             CALCULATE TOTAL DISPLACEMENT AT TIME T+DT:  UTDT=UT+U(0)
CE             STORE DISPLACEMENTS ON DISK
CS             RACUNANJE UKUPNOG POMERANJA U TRENUTKU T+DT:  UTDT=UT+U(0)
CS             ZAPISIVANJE POMERANJA NA DISK
C
               CALL UKUPNA
C
               IF(JPS.GT.1) THEN
                  IF(ITERGL.EQ.0) THEN
                     LID=LRAD
                     LCVEL=LID+NP*6
                     NPM=NPA-NPI+1
                     LELCV=LCVEL+NP
                     LMAX=LELCV+NPM
                     NP6=NP*7+NPM
                     LMAX13=NPODS(JPBR,2)-1
                     CALL IREADD(A(LID),NP6,IPODS,LMAX13,LDUZI)
C
CE                   PRINT CONTOUR DISPLACEMENTS 
CS                   STAMPANJE KONTURNIH POMERANJA
C
                     CALL STAMPA(A(LIPODS))
                  ENDIF
                  DO 650 JPBR=1,JPS
                     JEDN=NPODS(JPBR,6)
                     JEDNP=NPODS(JPBR,23)
                     JED=JEDN-JEDNP
                     LRTG=LRAD
                     LLMG=LRTG+JEDN*IDVA
                     LMAX=LLMG+JED
                     IF(LMAX.GT.MTOT) CALL ERROR(1)
                     LMAX13=NPODS(JPBR,27)-1
                     CALL IREADD(A(LLMG),JED,IPODS,LMAX13,LDUZI)
                     CALL RSTAZ(NPODS,LRTG,47)
                     LRTGG=LRTG+JEDNP*IDVA
                     CALL SPAKUU(A(LUPRI),A(LRTGG),A(LLMG),JED)
                     CALL WSTAZ(NPODS,LRTG,47)
  650             CONTINUE
                  LRAD=JMAXA
               ENDIF
            ENDIF
C
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(JPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            DO 850 JPBB=1,JPS
               IF(JPS.GT.1) THEN
                  IF (myid.eq.0) THEN
                     JPBR=JPBB
                     CALL PODDAT(NPODS,1)
                     CALL PODDAT(NPODS,3)
                     CALL PODDAT(NPODS,4)
                  ENDIF
                  write(*,*) 'pre drugog poziva resen 
     1                          iz pak07 sabrutine',myid
                  CALL RESEN(ALSK,A(LRTDT),A(LMAXA),JEDN,3)
                  IF (myid.eq.0) CALL UKUPNA
               ELSE
                  IF (myid.eq.0) JPBR=JPS1
               ENDIF
               IF (myid.eq.0) THEN
                  IF(ITERGL.EQ.0) THEN
C
CE                   CALCULATE STRESSES
CS                   RACUNANJE NAPONA 
C
                     NAPON=1
                     CALL INTNAP(A(LIGRUP))
C
CE                   R E S U L T S
CS                   R E Z U L T A T I
C
CE                   PRINT DISPLACEMENTS AND STRESSES
CS                   STAMPANJE POMERANJA I NAPONA
C
                     CALL STAMPA(A(LIPODS))
                  ENDIF
C
                  IF(JPS.GT.1) LRAD=LID
               ENDIF
  850       CONTINUE
C
CE          ITERATION METHODS
CS          ITERATIVNI METODI
C
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(ITERGL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            IF(ITERGL.GT.0) THEN
C CERNI
               IF (myid.eq.0) THEN
                  IF(ICERNI.EQ.1) THEN
C                    ZAPISIVANJE POMERANJA U ITERACIJAMA
                     IDUM=27
                     III=KOR
                     KOR=ITER+1
                     VRE=VREME
                     VREME=KOR
                     STOSZ=VREME
C                     CALL STAGP1(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,IDUM,0,
C     +                        A(LNCVP),NCVPR)
                     CALL STAU09(A(LRTDT),A(LID),A(LCVEL),ICVEL,NP,IDUM,
     +                        1,A(LNCVP),NCVPR)
                     KOR=III
                     VREME=VRE
                  ENDIF
               ENDIF
C CERNI
C
               CALL IMETOD
C CERNI
               IF (myid.eq.0) THEN
                  IF(ICERNI.EQ.1) THEN
                     IDUM=27
                     CALL NEUTRC(MAXIT,IDUM,2)
                     CLOSE(IDUM,STATUS='KEEP')
                  ENDIF  
               ENDIF 
C CERNI
            ENDIF
C
         ELSE
C
C
CE          L I N E A R     A N A L Y S I S
CS          L I N E A R N A      A N A L I Z A
C
            IF (myid.eq.0) THEN
               IF(JPS.GT.1) THEN
                  LID=LRAD
                  LCVEL=LID+NP*6
                  NPM=NPA-NPI+1
                  LELCV=LCVEL+NP
                  LMAX=LELCV+NPM
                  NP6=NP*7+NPM
                  LMAX13=NPODS(JPBR,2)-1
                  CALL IREADD(A(LID),NP6,IPODS,LMAX13,LDUZI)
C
CE                PRINT CONTOUR DISPLACEMENTS 
CS                STAMPANJE KONTURNIH POMERANJA
C
                  CALL STAMPA(A(LIPODS))
C
                  DO 600 JPBR=1,JPS
                     JEDN=NPODS(JPBR,6)
                     JEDNP=NPODS(JPBR,23)
                     JED=JEDN-JEDNP
                     LRTG=LRAD
                     LLMG=LRTG+JEDN*IDVA
                     LMAX=LLMG+JED
                     IF(LMAX.GT.MTOT) CALL ERROR(1)
                     LMAX13=NPODS(JPBR,27)-1
                     CALL IREADD(A(LLMG),JED,IPODS,LMAX13,LDUZI)
                     LMAX13=NPODS(JPBR,52)-1
                     CALL READDD(A(LRTG),JEDN,IPODS,LMAX13,LDUZI)
                     LRTGG=LRTG+JEDNP*IDVA
                     CALL SPAKUU(A(LRTDT),A(LRTGG),A(LLMG),JED)
                     LMAX13=NPODS(JPBR,52)-1
                     CALL WRITDD(A(LRTG),JEDN,IPODS,LMAX13,LDUZI)
  600             CONTINUE
                  LRAD=JMAXA
               ENDIF
            ENDIF
c
            CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(JPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
C
            DO 800 JPBB=1,JPS
C           OVO SU PODSTRUKTURE KOJE SE NE KORISTE
               IF(JPS.GT.1) THEN
                  IF (myid.eq.0) THEN
                     JPBR=JPBB
                     CALL PODDAT(NPODS,1)
                     IF(NPUP.EQ.1) GO TO 800
                     CALL PODDAT(NPODS,3)
                     CALL PODDAT(NPODS,4)
                  ENDIF
C                  CALL RESEN(A(LSK),A(LRTDT),A(LMAXA),JEDN,3)
                  IF (myid.eq.0) THEN
                     LMAX13=NPODS(JPBR,52)-1
                     CALL WRITDD(A(LRTDT),JEDN,IPODS,LMAX13,LDUZI)
                  ENDIF
               ELSE
                  IF (myid.eq.0) JPBR=JPS1
               ENDIF
C
CE             CALCULATE OF STRESSES FOR LINEAR ANALYSIS
CS             RACUNANJE NAPONA ZA LINEARNU ANALIZU
C
               IF (myid.eq.0) THEN
                  NAPON=1
                  IF(ITEST.NE.1) CALL INTNAP(A(LIGRUP))
CSRBA
                  IF(ISRBA.EQ.1) THEN
C UCITAVANJE POMERANJA IZ PRETHODNOG KORAKA
C LFTDT - REPER ZA POMERANJA IZ PREDHODNOG KORAKA
                     CALL RSTAZ(NPODS,LFTDT,87)
C         CALL WRR6(A(LFTDT),JEDN,'POMR')
C POZIVANJE PROGRAMA ZA RACUNANJE BRZINE U CVOROVIMA I RACUNANJE NOVIH
C VREDNOSTI (P) I (TAU) I ISPITIVANJE KONVERGENCIJE
C LRTDT - REPER ZA TEKUCA POMERANJA 
C LID - REPER ZA MATRICU VEZE IZMEDJU CVOROVA I BROJEVA JEDNACINA 
C DT - VREMENSKI KORAK
C LFAKP - REPER ZA VREDNOSTI (P) ILI (TAU)
C LIPRAV - REPER ZA INDIKATORE PRAVCA (0-P,-1-TAU)
C LNODPR - REPER ZA BROJEVE CVOROVA NA GRANICNOJ LINIJI
C NPP2 - BROJ UCITANIH GRANICNIH LINIJA SA (P) ILI (TAU)
C KSRBA - INDIKATOR KONVERGENCIJE (0-NE,1-KONVERGIRAO)
                     KSRBA=0
                     CALL SRBAP(A(LRTDT),A(LFTDT),A(LID),DT,NP,JEDN,
     1                       A(LFAKP),A(LIPRAV),A(LNODPR),NPP2,KSRBA)
C MORA DA OBRISE OVAJ PROSTOR JER SE KORISTI ZA UNUTRASNJE SILE
                     CALL CLEAR(A(LFTDT),JEDN)
C POVECAVA BROJ ITERACIJA U KORAKU
                     ITER=ITER+1
                     IF(KSRBA.EQ.0) GO TO 888
                  ENDIF
CSRBA
C
CE             R E S U L T S
CS             R E Z U L T A T I
C
CE             PRINT DISPLACEMENTS AND STRESSES
CS             STAMPANJE POMERANJA I NAPONA
C
                  CALL STAMPA(A(LIPODS))
C
                  IF(JPS.GT.1) LRAD=LID
               ENDIF
  800       CONTINUE
C
         ENDIF
C
      ENDIF
C
      IF (myid.eq.0) THEN
         IF(KOCID.EQ.1.AND.INDBACK.EQ.0) THEN
C
CE          PREPEAR FOR RESTART
CS          PRIPREMA ZA RESTART
C
C CERNI
C           PREPISIVANJE DISKA ZA ELEMENTE POSLE KONVERGENCIJE ZBOG ALFI
            IDUM=28
            IMD='ZBELEM'
            OPEN (IDUM,FILE=IMD,STATUS='UNKNOWN',
     1         FORM='UNFORMATTED',ACCESS='DIRECT',RECL=LDUZI)
            CALL PREPEL(A(LRAD),IELEM,IDUM,LDUZI)
            CLOSE(IDUM,STATUS='KEEP')
C CERNI
            CALL BACKUP(1)
         ENDIF
C
c        IF(SMKOR.GT.1.D-10) THEN
         IF(NBLPR.GE.0) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2200) KOR,NGRKOR,NSMKOR,SMKOR
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6200) KOR,NGRKOR,NSMKOR,SMKOR
         ENDIF
      ENDIF
c        ENDIF
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(' *** K O R A K    B R O J ',I5,' /',I4,
     +       ',  VREME ',1PD12.5,' ***')
 2100 FORMAT(//' PROGRAM STOP'/' ZA METOD =',I5/
     1' NE MOGU SE RACUNATI SOPSTVENE VREDNOSTI'/)
 2200 FORMAT(//' I Z V E S T A J   N A   K R A J U   K O R A K A ',I5//
     +5X,' MAKSIMALNI EFEKTIVNI NAPON:'/
     +5X,' GRUPA =',I5,',  ELEMENT =',I6,',  MAX.EFE.NAPON =',1PD11.3//)
C-----------------------------------------------------------------------
 6000 FORMAT(' *** S T E P    N U M B E R ',I5,' /',I4,
     +       ',  TIME ',1PD12.5,' ***')
 6100 FORMAT(//' PROGRAM STOP'/' FOR METHOD =',I5/
     1' NO SOLVE EIGEN VALUE'/)
 6200 FORMAT(//' R E P O R T   O N   E N D   O F   S T E P ',I5//
     +5X,' MAXIMUM EFFECTIVE STRESS:'/
     +5X,' GROUP =',I5,', ELEMENT =',I6,',  MAX.EFF.STRESS =',1PD11.3//)
C-----------------------------------------------------------------------
      END
C=======================================================================
C
C=======================================================================
      SUBROUTINE EPAKS(NPODS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.        WITH LOOP OVER TIME PERIODS AND STEPS
CS.    P R O G R A M
CS.        SA PETLJOM PO VREMENSKIM PERIODIMA I KORACIMA
C .
C ......................................................................
C
      include 'paka.inc'
      
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /PERKOR/ LNKDT,LDTDT,LVDT,NDT,DT,VREME,KOR
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /SOPSVR/ ISOPS,ISTYP,NSOPV,ISTSV,IPROV,IPROL
      COMMON /MAXREZ/ PMALL,BMALL,AMALL,SMKOR,SMALL,
     +                NPMALL,NBMALL,NAMALL,KPMALL,KBMALL,KAMALL,
     +                NSMKOR,NSMALL,NGRKOR,NGRALL,KSMALL
      COMMON /BEAMAX/ FZATM,FSMICM,FTORM,FSAVM,NEB1KOR,NPB1KOR,NEB2KOR,
     1                NPB2KOR,NEB4KOR,NPB4KOR,NEB5KOR,NPB5KOR,
     1                NMB1KOR,NMB2KOR,NMB4KOR,NMB5KOR,IMAGREDA 
      COMMON /ITERBR/ ITER
      COMMON /DINAMI/ IMASS,IDAMP,PIP,DIP,MDVI
      COMMON /SISTEM/ LSK,LRTDT,NWK,JEDN,LFTDT
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /OPTERE/ NCF,NPP2,NPP3,NPGR,NPGRI,NPLJ,NTEMP
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /ITERAC/ METOD,MAXIT,TOLE,TOLS,TOLM,KONVE,KONVS,KONVM
      COMMON /AUTSTP/ DTUK,ALFG,DELS,IAUTO,ITEOPT,KPNOD,KPDIR,KEQ
      COMMON /CVOREL/ ICVEL,LCVEL,LELCV,NPA,NPI,LCEL,LELC,NMA,NMI
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /DUZINA/ LMAX,MTOT,LMAXM,LRAD,NRAD
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /ZADATA/ LNZADJ,LNZADF,LZADFM,NZADP
      COMMON /ECLANM/ AMAXK,AMINK,AMAXF,AMINF
      COMMON /ANALIZ/ LINEAR,ITERGL,INDDIN
      COMMON /IMPERF/ NMODS,LIDIM,LSCIM,MODES
      COMMON /UPDLAG/ LUL,LCORUL
      COMMON /UPRIRI/ LUPRI
      COMMON /RESTAR/ TSTART,IREST
      COMMON /NELINE/ NGENN
      COMMON /DUPLAP/ IDVA
      COMMON /SRPSKI/ ISRPS
      COMMON /GRUPER/ LIGRUP
      COMMON /LINSBR/ LIN,KORBR
      COMMON /INDKON/ IKONV
      COMMON /INDNAP/ NAPON
      COMMON /RADIZA/ INDBG
      COMMON /STAMKO/ ISTKO,NCVPR,LNCVP,LNCVZ,
     +                ISTEM,ISTVN,ISTSI,ISTDE,ISTNA
      COMMON /ZAPSIL/ UBXYZ(3,10),NVFUB(10),INDZS,LZAPS,LNPRZ
      COMMON /KONTKT/ ICONT,NEQC,NEQ,NWKC,LMAXAC,LRCTDT
      COMMON /BLOCKS/ NBMAX,IBLK,NBLOCK,LMNQ,LICPL,LLREC,KC,LR
      COMMON /PRIT2D/ LFAKP,LTHICV,LITIPE,LNFUN,LIPRAV,LNODPR,ISRBA
      COMMON /JCERNI/ LKOLKO,ICERNI
      COMMON /STOSZP/ STOSZ
      COMMON /CDEBUG/ IDEBUG
      INCLUDE 'mpif.h'
      integer myid,ierr
      DIMENSION NPODS(JPS1,*)
C
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      IF(IDEBUG.GT.0) PRINT *, ' EPAKS'
C
CE    E I G E N    V A L U E S  -  S T A B I L I T Y    A N A L Y S I S
CS    S O P S T V E N E    V R E D N O S T I  -  S T A B I L N O S T
C
CZILESK
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(ISOPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      IF(ISOPS.GT.0) THEN
            if (myid.ne.0) goto 10
C
            IPROL=1
C
CE          READING OF RIGHT SIDE VECTOR FROM DISK 2 MECHANICAL LOADS, RTDT
CS          UCITAVANJE VEKTORA DESNE STRANE SA DISKA 2 MEHANICKE SILE, RTDT
C
            CALL RSTAZ(NPODS,LRTDT,38)
C
CE          FORM VECTOR OF GEOMETRY NONLINEAR PRESSURE LOADING IN (T+DT)
CS          FORMIR. VEKTORA PRITISAKA GEOMETRIJSKI NELINEARNIH U (T+DT)
C
            IF(NGEOM.NE.0) CALL OPTLIN(NPODS)
C
            CALL WSTAZ(NPODS,LRTDT,39)
C
CE          LOOP OVER GROUPS OF ELEMENT TO INTEGRATION OF NONLIN. MATRIX K
CS          PETLJA PO GRUPAMA ELEMENATA RADI INTEGRACIJE NELIN. MATRICE K
C
10       CALL INTNMK(A(LIGRUP),NPODS)
         if (myid.ne.0) return
         IF(NBLOCK.EQ.1) THEN
            IPROL=2
            LSK2=LSK
            LSK=LSK+NWK*IDVA
            CALL INTNMK(A(LIGRUP),NPODS)
            LSK=LSK2
         ENDIF
C
      ENDIF
      if (myid.ne.0) return
CZILESK
      IF(NBLPR.GE.0) THEN
         IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2300)
         IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6300)
         IF(PMALL.GT.1.D-10) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2400) KPMALL,NPMALL,PMALL
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6400) KPMALL,NPMALL,PMALL
         ENDIF
         IF(NDIN.EQ.1) THEN
            IF(BMALL.GT.1.D-10) THEN
               IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2500) KBMALL,NBMALL,BMALL
               IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6500) KBMALL,NBMALL,BMALL
            ENDIF
            IF(AMALL.GT.1.D-10) THEN
               IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2600) KAMALL,NAMALL,AMALL
               IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6600) KAMALL,NAMALL,AMALL
            ENDIF
         ENDIF
         IF(SMALL.GT.1.D-10) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2700) KSMALL,NGRALL,NSMALL,SMALL
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6700) KSMALL,NGRALL,NSMALL,SMALL
         ENDIF
      ENDIF
C
      IF(IMAGREDA.EQ.1) THEN
         IF(DABS(FZATM).GT.1.D-10) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2800) NMB1KOR,NEB1KOR,NPB1KOR,FZATM
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6800) NMB1KOR,NEB1KOR,NPB1KOR,FZATM
         ENDIF
         IF(DABS(FSMICM).GT.1.D-10) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2810) NMB2KOR,NEB2KOR,NPB2KOR,FSMICM
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6810) NMB2KOR,NEB2KOR,NPB2KOR,FSMICM
         ENDIF
         IF(DABS(FTORM).GT.1.D-10) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2820) NMB4KOR,NEB4KOR,NPB4KOR,FTORM
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6820) NMB4KOR,NEB4KOR,NPB4KOR,FTORM
         ENDIF
         IF(DABS(FSAVM).GT.1.D-10) THEN
            IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2830) NMB5KOR,NEB5KOR,NPB5KOR,FSAVM
            IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6830) NMB5KOR,NEB5KOR,NPB5KOR,FSAVM
         ENDIF
      ENDIF
      RETURN
C-----------------------------------------------------------------------
 2300 FORMAT(//' I Z V E S T A J   N A   K R A J U   P R O G R A M A')
 2400 FORMAT(/5X,' MAKSIMALNO POMERANJE:'/
     +5X,' KORAK =',I5,',  CVOR =',I6,',   MAX.POM. =',1PD12.4)
 2500 FORMAT(/' MAKSIMALNA BRZINA:'/
     +5X,' KORAK =',I5,',  CVOR =',I6,',   MAX.BRZ. =',1PD12.4)
 2600 FORMAT(/' MAKSIMALNO UBRZANJE:'/
     +5X,' KORAK =',I5,',  CVOR =',I6,',   MAX.UBR. =',1PD12.4)
 2700 FORMAT(/5X,' MAKSIMALNI EFEKTIVNI NAPON:'/5X,' KORAK =',I5,
     +', GRUPA =',I5,', ELEMENT =',I6,',  MAX.EFE.NAPON =',1PD11.3)
 2800 FORMAT(/5X,' MAKSIMALNA NORMALNA SILA:'/5X,' KORAK =',I5,
     +', ELEMENT =',I5,', CVOR =',I6,',  MAX. NOR. SILA =',1PD11.3)
 2810 FORMAT(/5X,' MAKSIMALNA SMICUCA SILA:'/5X,' KORAK =',I5,
     +', ELEMENT =',I5,', CVOR =',I6,',  MAX. SMIC. SILA =',1PD11.3)
 2820 FORMAT(/5X,' MAKSIMALNI TORZIONI MOMENT:'/5X,' KORAK =',I5,
     +', ELEMENT =',I5,', CVOR =',I6,',  MAX. TOR. MOMENT =',1PD11.3)
 2830 FORMAT(/5X,' MAKSIMALNI SAVOJNI MOMENT:'/5X,' KORAK =',I5,
     +', ELEMENT =',I5,', CVOR =',I6,',  MAX. SAV. MOMENT =',1PD11.3)
C-----------------------------------------------------------------------
 6300 FORMAT(//' R E P O R T   O N   E N D   O F   P R O G R A M')
 6400 FORMAT(/5X,' MAXIMUM DISPLACEMENT:'/
     +5X,' STEP =',I5,',  NODE =',I6,',   MAX.DISPL. =',1PD12.4)
 6500 FORMAT(/5X,' MAXIMUM VELOCITY:'/
     +5X,' STEP =',I5,',  NODE =',I6,',   MAX.VELOC. =',1PD12.4)
 6600 FORMAT(/5X,' MAXIMUM ACCELERATION:'/
     +5X,' STEP =',I5,',  NODE =',I6,',   MAX.ACCEL. =',1PD12.4)
 6700 FORMAT(/5X,' MAXIMUM EFFECTIVE STRESS:'/5X,' STEP =',I5,
     +', GROUP =',I5,', ELEMENT =',I6,',  MAX.EFF.STRESS =',1PD11.3)
 6800 FORMAT(/5X,' MAXIMUM NORMAL FORCE:'/5X,' STEP =',I5,
     +', ELEMENT =',I5,', NODE =',I6,',  MAX. NOR. FORCE =',1PD11.3)
 6810 FORMAT(/5X,' MAXIMUM SHEAR FORCE:'/5X,' STEP =',I5,
     +', ELEMENT =',I5,', NODE =',I6,',  MAX. SHEAR FORCE =',1PD11.3)
 6820 FORMAT(/5X,' MAXIMUM TORSION MOMENT:'/5X,' STEP =',I5,
     +', ELEMENT =',I5,', NODE =',I6,',  MAX. TOR. MOMENT =',1PD11.3)
 6830 FORMAT(/5X,' MAXIMUM BENDING MOMENT:'/5X,' STEP =',I5,
     +', ELEMENT =',I5,', NODE =',I6,',  MAX. BEND. MOMENT =',1PD11.3)
C-----------------------------------------------------------------------
      END
C=======================================================================
      SUBROUTINE REWDRV
C
C  NOVE VREDNOSTI ZA REFERENTNE JEDINICNE VEKTORE
C  TAU = TRENUTAK FORMIRANJA JEDINICNE MATRICE
C
      include 'paka.inc'
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' REWDRV'
C
      IF(IDIREK.GT.0.OR.NGENL.EQ.0) RETURN
      NPROS=NP*3
C
C....   UCITAVANJE  VN(I)  ,   V1(I)   
C       U           A(LDRV),   A(LDRV1)
      LMA8=LSTAZZ(3)
      CALL READDD(A(LDRV), NPROS,IPODS,LMA8,LDUZI)
      LMA8=LSTAZZ(4)
      CALL READDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
C
C
C.....  ZAPIS  VN(TAU) I  V1(TAU)
C       TAU = REFERENTNI TRENUTAK KADA SE FORMIRA MATR. KRUTOSTI
C
      LMA8=LSTAZZ(1)
      CALL WRITDD(A(LDRV), NPROS,IPODS,LMA8,LDUZI)
      LMA8=LSTAZZ(2)
      CALL WRITDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
      RETURN
      END
C=======================================================================
      SUBROUTINE REWDRA (IOPT)
C
C  NOVE VREDNOSTI ZA REFERENTNE JEDINICNE VEKTORE
C     IOPT = 1    REFERENCA POSTAJE  V(T)
C     IOPT = 2    REFERENCA POSTAJE  V(T+DT)
C  RUTINA SE KORISTI ZA METODE AUTOMATSKOG INKREMENTIRANJA
C
      include 'paka.inc'
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' REWDRA'
C
      IF(IDIREK.GT.0.OR.NGENL.EQ.0) RETURN
      NPROS=NP*3
C
C....   UCITAVANJE  VN     ,   V1
C       U           A(LDRV),   A(LDRV1)
      LMA8=LSTAZZ(7)
      IF(IOPT.EQ.2) LMA8=LSTAZZ(3)
      CALL READDD(A(LDRV), NPROS,IPODS,LMA8,LDUZI)
      LMA8=LSTAZZ(8)
      IF(IOPT.EQ.2) LMA8=LSTAZZ(4)
      CALL READDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
C
C
C.....  ZAPIS  VN(TAU) I  V1(TAU)
C       TAU = REFERENTNI TRENUTAK KADA SE FORMIRA MATR. KRUTOSTI
C
      LMA8=LSTAZZ(1)
      IF(IOPT.EQ.2) LMA8=LSTAZZ(7)
      CALL WRITDD(A(LDRV), NPROS,IPODS,LMA8,LDUZI)
      LMA8=LSTAZZ(2)
      IF(IOPT.EQ.2) LMA8=LSTAZZ(8)
      CALL WRITDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
C
      IF(IOPT.EQ.1) RETURN
      LMA8=LSTAZZ(3)
      CALL WRITDD(A(LDRV), NPROS,IPODS,LMA8,LDUZI)
      LMA8=LSTAZZ(4)
      CALL WRITDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE ROTDRV (IOPT)
C
C  ROTIRANJE JEDINICNIH VEKTORA U ITERACIJI
C
C
C   FORMIRANJE VEKTORA NORMALE U TRENUTKU T+DT(I) NA OSNOVU 
C   VEKTORA U TRENUTKU T+DT(I-1) - PRETHODNE ITERACIJE  I ROTACIJA 
C   REFERENTNOG VEKTORA U TRENUTKU (TAU) U ITERACIJI
C   RAZLIKA VN(T+DT)-VN(0) SMESTA SE U VN(0)
C
      include 'paka.inc'
      COMMON /GLAVNI/ NP,NGELEM,NMATM,NPER,
     1                IOPGL(6),KOSI,NDIN,ITEST
      COMMON /GRUPEE/ NGEL,NGENL,LGEOM,NGEOM,ITERM
      COMMON /REPERI/ LCORD,LID,LMAXA,LMHT
      COMMON /UPRIRI/ LUPRI
      COMMON /TRAKEJ/ IULAZ,IZLAZ,IELEM,ISILE,IRTDT,IFTDT,ILISK,ILISE,
     1                ILIMC,ILDLT,IGRAF,IDINA,IPOME,IPRIT,LDUZI
      COMMON /OPSTIP/ JPS,JPBR,NPG,JIDG,JCORG,JCVEL,JELCV,NGA,NGI,NPK,
     1                NPUP,LIPODS,IPODS,LMAX13,MAX13,JEDNG,JMAXA,JEDNP,
     1                NWP,NWG,IDF,JPS1
      COMMON /DIREKT/ LSTAZZ(9),LDRV0,LDRV1,LDRV,IDIREK
      COMMON /DUPLAP/ IDVA
      COMMON /CDEBUG/ IDEBUG
      IF(IDEBUG.GT.0) PRINT *, ' ROTDRV'
C
      IF(IDIREK.GT.0.OR.NGENL.EQ.0) RETURN
      NPROS=NP*3
      LDRV2=LDRV0+NPROS*IDVA
C
C
C.... UCITATI REFERENTNE VEKTORE  VN(I-1) I V1(I-1)
      LMA8=LSTAZZ(1)
      CALL READDD(A(LDRV), NPROS,IPODS,LMA8,LDUZI)
      LMA8=LSTAZZ(2)
      CALL READDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
C
C
C.... ROTIRANJE  VN(I-1) I V1(I-1)  OKO REFERENTNIH VEKTORA
      CALL DRGTDT(A(LDRV),A(LDRV1),A(LUPRI),A(LID),NP)
C
C.... ZAPIS TEKUCE VREDNOSTI VN(I) I V1(I)
      LMA8=LSTAZZ(3)
      IF(IOPT.EQ.2) LMA8=LSTAZZ(6)
      CALL WRITDD(A(LDRV),NPROS,IPODS,LMA8,LDUZI)
      IF(IOPT.EQ.1) THEN
         LMA8=LSTAZZ(4)
         CALL WRITDD(A(LDRV1),NPROS,IPODS,LMA8,LDUZI)
      ENDIF
C
C
C.... UCITAVANJE     VN(0) 
C     U              A(LDRV0)
      LMA8=LSTAZZ(5)
      CALL READDD(A(LDRV0),2*NPROS,IPODS,LMA8,LDUZI)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE DRGTDT(DRV,DRV1,UP,ID,NP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     FORMIRANJE VEKTORA NORMALE I VEKTORA V1 U TRENUTKU T+DT NA OSNOVU 
C     VEKTORA U TRENUTKU T  I ROTACIJA U KORAKU
C
      include 'paka.inc'
      
      COMMON /MPOINC/ MMP,NMPC,NEZAV,LCMPC,LMPC,NEZA1
      DIMENSION DRV(NP,*),DRV1(NP,*),UP(*),ID(NP,*)
      DIMENSION EF123(3,3),QMTX(3,3),S(3,3)
C
      DO 47 IG=1,NP
C
C......   UGLOVI ROTACIJE U KORAKU
      IDF=ID(IG,4)
      IF(IDF)1,2,3
    1 ALFA=CONDOF(UP,A(LCMPC),A(LMPC),IDF)
      GO TO 4
    2 ALFA=0.D0
      GO TO 4
    3 ALFA=UP(IDF)
C
    4 IDF=ID(IG,5)
      IF(IDF)5,6,7
    5 BETA=CONDOF(UP,A(LCMPC),A(LMPC),IDF)
      GO TO 8
    6 BETA=0.D0 
      GO TO 8
    7 BETA=UP(IDF)
C
    8 IDF=ID(IG,6)
      IF(IDF)9,10,11
    9 GAMA=CONDOF(UP,A(LCMPC),A(LMPC),IDF)
      GO TO 12
   10 GAMA=0.D0 
      GO TO 12
   11 GAMA=UP(IDF)
C
   12 DO 48 I=1,3
      EF123(I,1)=DRV1(IG,I)
   48 EF123(I,3)=DRV(IG,I)
      CALL V1V2(EF123(1,1),EF123(1,2),EF123(1,3),2)
C
      IF(DABS(ALFA).LT.1.D-10.AND.DABS(BETA).LT.1.D-10.AND.
     1   DABS(GAMA).LT.1.D-10)GO TO 47
C
         CALL QROTUP(ALFA,BETA,GAMA,QMTX,S)
C
C        ROTIRANJE  S = EF123 * QMTX
         CALL MNOZM1(S,EF123,QMTX,3,3,3)
C
C...     V3=VN I V1
         DO 51 K=1,3
            DRV(IG,K)=S(K,3)
            DRV1(IG,K)=S(K,1)
   51    CONTINUE
C
   47 CONTINUE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE QROTUP(ALFA,BETA,GAMA,QMTX,S)
C
C   HUGHES/WINGET MATRICA ROTACIJE
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION QMTX(3,*),S(3,*)
C   INTEZITET PSEUDO VEKTORA ROTACIJE
      THETA=ALFA*ALFA+BETA*BETA+GAMA*GAMA
      CC=1./(2.+0.5*THETA)
C   TENZOR ROTACIJE
      S(1,1)=0.D0
      S(2,2)=0.D0
      S(3,3)=0.D0
      S(1,2)=-GAMA
      S(1,3)= BETA
      S(2,3)=-ALFA
      S(2,1)= GAMA
      S(3,1)=-BETA
      S(3,2)= ALFA
C   QMTX  =  1 + CC*(S*S + 2*S)
      DO 20 I=1,3
      DO 10 J=1,3
       QMTX(I,J)=0.D0
       DO 5 K=1,3
        QMTX(I,J)=QMTX(I,J)+S(I,K)*S(K,J)
    5  CONTINUE
       QMTX(I,J)=CC*(QMTX(I,J)+2.*S(I,J))
   10 CONTINUE
   20 CONTINUE
      DO 30 I=1,3
       QMTX(I,I)=1.+QMTX(I,I)
   30 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SILMOM(ID,RTDT,NP,IZLAZ,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CS       RACUNANJE UKUPNIH OPTERECENJA
CE       TOTAL LOADS
C
      COMMON /POSTPR/ LNDTPR,LNDTGR,NBLPR,NBLGR,INDPR,INDGR
      DIMENSION ID(NP,*),RTDT(*),SIL(6)
C
      CALL CLEAR(SIL,6)
      DO 10 I=1,NP
      DO 10 J=1,6
         NJ=ID(I,J)
         IF(NJ.EQ.0) GO TO 10
         SIL(J)=SIL(J)+RTDT(NJ)
   10 CONTINUE
      IF(NBLPR.GE.0) THEN
      IF(ISRPS.EQ.0)
     1WRITE(IZLAZ,2000)
      IF(ISRPS.EQ.1)
     1WRITE(IZLAZ,6000)
      WRITE(IZLAZ,5000) (SIL(J),J=1,6)
      ENDIF
      RETURN
 5000 FORMAT(11X,'FX',11X,'FY',11X,'FZ',11X,'MX',11X,'MY',11X,'MZ'/
     1' ',6(1PD12.5,1X))
C-----------------------------------------------------------------------
 2000 FORMAT(///15X,' U K U P N E   S I L E   I   M O M E N T I')
C-----------------------------------------------------------------------
 6000 FORMAT(///12X,' T O T A L   F O R C E S   A N D   M O M E N T S')
C-----------------------------------------------------------------------
      END
C======================================================================
C
C======================================================================
      SUBROUTINE PREPEL(A,IDR,IDW,LD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CE.    P R O G R A M
CE.       TO COPY DIRECT ACCESS DISK
CS.    P R O G R A M
CS        ZA KOPIRANJE DISKA SA DIREKTNIM PRISTUPOM
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION A(*)
C
      IF(IDEBUG.GT.0) PRINT *, ' PREPEL'
      NN=LD/8
      LS=1
C NDP  
C   10 READ(IDR,END=100,ERR=100,REC=LS) (A(I),I=1,NN)
C FTN77
   10 READ(IDR,ERR=100,REC=LS) (A(I),I=1,NN)
      WRITE(IDW,REC=LS) (A(I),I=1,NN)
      LS=LS+1
      GO TO 10
  100 RETURN
      END
C======================================================================
C
C======================================================================
      SUBROUTINE SRBAP(RTDT,FTDT,ID,DT,NP,JEDN,
     1                 FAKP,IPRAV,NODPR,NPP2,KSRBA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C ......................................................................
C .
CS.    P R O G R A M
CS.       ZA RACUNANJE BRZINE U CVOROVIMA I RACUNANJE NOVIH
CS.       VREDNOSTI (P) I (TAU) I ISPITIVANJE KONVERGENCIJE
C .
C . RTDT(JEDN) - VEKTOR TEKUCA POMERANJA 
C . FTDT(JEDN) - VEKTOR POMERANJA IZ PRETHODNOG KORAKA
C . ID(NP,6) - MATRICA VEZE IZMEDJU CVOROVA I BROJEVA JEDNACINA 
C . DT - VREMENSKI KORAK
C . NP - UKUPAN BROJ CVOROVA
C . JEDN - UKUPAN BROJ JEDNCINA
C . FAKP(NPP2,2) - VREDNOSTI (P) ILI (TAU)
C . IPRAV(NPP2) - REPER ZA INDIKATORE PRAVCA (0-P,-1-TAU)
C . NODPR(NPP2,3) -  BROJEVI CVOROVA NA GRANICNOJ LINIJI
C . NPP2 - BROJ UCITANIH GRANICNIH LINIJA SA (P) ILI (TAU)
C . KSRBA - INDIKATOR ZA KONVERGENCIJU (0-NE,1-DA)
C .
C ......................................................................
C
      COMMON /CDEBUG/ IDEBUG
      DIMENSION IPRAV(*),FAKP(NPP2,*),NODPR(NPP2,*),
     1          RTDT(*),FTDT(*),ID(NP,*)
C
      IF(IDEBUG.GT.0) PRINT *, ' SRBAP'
C
C                  CALL WRR(RTDT,JEDN,'RTDU')
C                  CALL WRR(FTDT,JEDN,'FTDU')
      DO 10 I=1,NPP2
C CVOROVI NA LINIJI 
         N1=NODPR(I,1)
         N2=NODPR(I,2)
C         WRITE(3,*) 'N1,N2,NP,I,NPP2',N1,N2,NP,I,NPP2
C JEDNACINE ZA CVOROVE U PRAVCU X I Y OSE
         J1X=ID(N1,1)
         J1Y=ID(N1,2)
         J2X=ID(N2,1)
         J2Y=ID(N2,2)
C         WRITE(3,*) 'J1X,J1Y,J2X,J2Y',J1X,J1Y,J2X,J2Y
C BRZINE U CVOROVIMA U (X) I (Y) PRAVCU
         V1X=0.
         IF(J1X.GT.0) V1X=(RTDT(J1X)-FTDT(J1X))/DT
         V1Y=0.
         IF(J1Y.GT.0) V1Y=(RTDT(J1Y)-FTDT(J1Y))/DT
         V2X=0.
         IF(J2X.GT.0) V2X=(RTDT(J2X)-FTDT(J2X))/DT
         V2Y=0.
         IF(J2Y.GT.0) V2Y=(RTDT(J2Y)-FTDT(J2Y))/DT
         WRITE(3,*) 'V1X,V1Y,V2X,V2Y',V1X,V1Y,V2X,V2Y
C INDIKATOR PRAVCA (IPR=0 (P),IPR=-1 (TAU))
         IPR=IPRAV(I)
C STARA VREDNOST P ILI TAU U CVOROVIMA
         P1=FAKP(I,1)
         P2=FAKP(I,2)
         WRITE(3,*) 'IPR,P1,P2',IPR,P1,P2
C OVDE TREBA POZVATI PROGRAM ZA RACUNANJE NOVIH VREDNOSTI ZA P ILI TAU
C        CALL ......(V1X,V1Y,V2X,V2Y,P1,P2,IPR....)
C
C NOVA VREDNOST P ILI TAU U CVOROVIMA
         FAKP(I,1)=P1
         FAKP(I,2)=P2
   10 CONTINUE
C INDIKATOR ZA KONVERGENCIJU (0-NE,1-DA)
      KSRBA=1
C
      RETURN
      END

