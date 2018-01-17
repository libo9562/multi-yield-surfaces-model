C     MAIN PROGRAME
C     VARIABLE CHART (IN THE MAIN PROGRAM)
C
C     AS      THE CENTER OF YIELD SURFACES
C     B       BULK MODULUS
C     C_SS    CURRENT STRESS
C     C_ST    CURRENT STRAIN
C     H_0     THE ELASTIC MODULUS
C     H       THE PLASTIC MODULUS
C     K       SIZE OF YIELD SURFACE
C     NLS     NUMBER OF LOAD STEPS
C     NYS     NUMBER OF YIELD SURFACES
C     N_S2S   STRESS OF STRAIN INFORMATION AS A TEMP VARIABLE
C     
      REAL *8 B,H_0,C_SS(7)/7*0.D0/,C_ST(7)/7*0.D0/,N_S2S(7)/7*0.D0/
      COMMON /S2S/ N_S2S
C     ---------------------------------------------------------------------
C     IMPORT DATA
      CALL IMPORT(H_0,B,NYS)
      OPEN(86,FILE='STRAIN.TXT')
      OPEN(87,FILE='STRESS.TXT')
      NLS=0
    3 READ(86,*,END=2,ERR=1),(N_S2S(J),J=1,6)
      CALL DECOMPOSE(C_ST)
      CALL MODEL(H_0,B,NYS,C_SS)
      WRITE(87,901),(N_S2S(J),J=1,6)
  901 FORMAT(6(1X,F13.5))
      NLS=NLS+1
      GOTO 3
    1 PRINT *,'ERR AT LINE',NLS
    2 CLOSE(86)
      CLOSE(87)
      PRINT *,NLS,'LINES STRESS HISTORY HAVE BEEN CALCUALTED'
      STOP
      END

C     THIS SUBROUTINE CALCULATE AND OUTPUT THE RESULT FORM M3
      SUBROUTINE MODEL(H_0,B,NYS,C_SS)
      COMMON /A/ AS /PARA/ K,H /S2S/ N_S2S
      REAL *8 AS(6,50),K(50),H(2,50)
      REAL *8 N_S2S(7),C_SS(7),KTR,KPL
      REAL *8 B,H_0,LA,V(6),QM(6),QMP(6)
C     --------------------------------------------------------------
C     STEP 1 INITIALIZE
      M=0
C     STEP 2 ELASTIC PREDICTOR
      CALL ELASTIC_P(H_0,B,C_SS)
C     STEP 3 CHECK FOR OVERSHOOTING OF NEXT YIELD SURFACE
      DO 62 IC=1,NYS+1
        IF(M.EQ.NYS) GOTO 63
        CALL GET_Q(M+1,QMP,KTR)
        IF(KTR.LE.K(M+1)) GO TO 63
        QM=QMP
        M=M+1
C       STEP 4 COMPUTE STRESS CORRECTION
        LA=(KTR-K(M))/H(2,M)
        IF(M.EQ.1) THEN
          KPL=H_0*LA
        ELSE
          LA=LA*H(2,M-1)/H(1,M-1)
          KPL=LA*H_0*(1-H(2,M)/H(2,M-1))
        ENDIF
C       STEP 5 UPDATES 
        CALL UPDATEDS((KTR-KPL),QM,M)
   62 CONTINUE
C     STEP 6 FINAL UPDATE
   63 CALL UPDATESS(C_SS)
      IF(M.GT.0) THEN 
        CALL GET_NA(NYS,M,H(1,M)*LA,K(M+1)-K(M),QM)
        IF(M.GT.1) THEN
          DO 65 I=1,M-1
            DO 66 J=1,6
              AS(J,I)=C_SS(J)-K(I)/K(M)*(C_SS(J)-AS(J,M))
   66       CONTINUE
   65     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END

C     THIS SUBROUTINE IMPORT YIELD SURFACE PRETREAMENT 
      SUBROUTINE IMPORT(H_0,B,NYS)
      COMMON /PARA/ K,H /A/ AS
      REAL *8 GMAX,RR,V,TEMP,GR,B,H_0,K(50),H(2,50),AS(6,50)
C     IMPORT SYSTEM INFORMATION AND MODEL PARAMENTERS
      OPEN(98, FILE='SYS.TXT')
      READ(98,*),GMAX,RR,V,TEMP
      CLOSE(98)
C     ELASTIC PARAMETERS
      NYS=IFIX(SNGL(TEMP))
      GR=GMAX*RR
      DELTA=GR/(TEMP+1.D0)
      H_0=2.D0*(GR-DELTA)/RR
      B=H_0*(1.D0+V)/3.D0/(1-2.D0*V)
      DO 108 I=1,NYS
        DO 109 J=1,6
          AS(J,I)=0.D0
  109   CONTINUE
  108 CONTINUE
C     FIND OTHER PARAMETERS H
      DO 110 I=1,NYS
        K(I)=DELTA*I*SQRT(2.D0)
        IF(I.LT.NYS)THEN
          H(1,I)=2.D0*(GR-(I+1)*DELTA)*(GR-I*DELTA)/GR/RR
        ELSE
          H(1,I)=0.D0
        ENDIF
        H(1,I)=H(1,I)*H_0/(H_0-H(1,I))
        H(2,I)=H_0+H(1,I)
  110 CONTINUE
      RETURN
      END
C
C     GET THE DECOMPOSED STRAIN INCREAMENT AND UPDATE THE TOTAL TRAIN
      SUBROUTINE DECOMPOSE(C_ST)
      COMMON /S2S/ IN_ST
      REAL *8 O_ST(7),C_ST(7),IN_ST(7),EV
      IN_ST(7)=IN_ST(1)+IN_ST(2)+IN_ST(3)
      DO 11 I=1,3
        IN_ST(I)=IN_ST(I)-IN_ST(7)/3.D0
   11 CONTINUE
      DO 12 I=1,7
        O_ST(I)=C_ST(I)
        C_ST(I)=IN_ST(I)
        IN_ST(I)=C_ST(I)-O_ST(I)
   12 CONTINUE  
      RETURN
      END

C     THIS SUBROUTINE GIVES THE ELASTIC PREDICTOR
      SUBROUTINE ELASTIC_P(H_0,B,C_SS)
      COMMON /S2S/ N_S2S
      REAL *8 H_0,B,C_SS(7),N_S2S(7)
      N_S2S(7)=N_S2S(7)*B/3.D0+C_SS(7)
      DO 21 I=1,6
        N_S2S(I)=N_S2S(I)*H_0+C_SS(I)
   21 CONTINUE
      RETURN
      END

C     THIS SUBROUTINE GIVES Q THE NORMAL OF YIELD SURFACE AND KTR
      SUBROUTINE GET_Q(MP,QM,KTR)
      COMMON /S2S/ N_S2S /A/ AS
      REAL *8 QM(6),KTR,N_S2S(7),AS(6,50)
      KTR=0.D0
      DO 41 I=1,3
        QM(I)=N_S2S(I)-AS(I,MP)
        QM(I+3)=N_S2S(I+3)-AS(I+3,MP)
        KTR=KTR+QM(I)**2+2*QM(I+3)**2
   41 CONTINUE
      IF(KTR.GT.0.D0) THEN
        KTR=SQRT(KTR)
        DO 42 I=1,6
          QM(I)=QM(I)/KTR
   42   CONTINUE
      ELSE
        KTR=0.D0
      ENDIF
      RETURN
      END

C     THIS SUBROUTINE GIVES A ASSOCIATE WITH M
      SUBROUTINE GET_NA(NYS,M,C,D,QM)
      COMMON /A/AS
      REAL *8 QM(6),AS(6,50),C,D,E,V(6)
      E=0.D0
      IF(M.LT.NYS) THEN
        DO 53 I=1,3
          V(I)=D*QM(I)+AS(I,M+1)-AS(I,M)
          V(I+3)=D*QM(I+3)+AS(I+3,M+1)-AS(I+3,M)
          E=E+V(I)*QM(I)+2*V(I+3)*QM(I+3)
   53   CONTINUE
        DO 54 I=1,6
          AS(I,M)=AS(I,M)+C/E*V(I)
   54   CONTINUE
      ELSE
        DO 55 I=1,6
          AS(I,M)=AS(I,M)+C*QM(I)
   55   CONTINUE
      ENDIF
      RETURN
      END

C     THIS SUBROUTINE UPDATE THE CALCULATED STRESS RESULTS
      SUBROUTINE UPDATESS(C_SS)
      COMMON /S2S/ N_S2S
      REAL *8 N_S2S(7),C_SS(7)
      C_SS(7)=N_S2S(7)
      DO 52 I=1,3
        C_SS(I)=N_S2S(I)
        N_S2S(I)=C_SS(I)+C_SS(7)
        C_SS(I+3)=N_S2S(I+3)
   52 CONTINUE
      RETURN
      END

      SUBROUTINE UPDATEDS(KTR,QM,M)
      COMMON /S2S/ N_S2S /A/AS
      REAL *8 KTR,QM(6),AS(6,50),N_S2S(7)
      DO 51 I=1,6
        N_S2S(I)=KTR*QM(I)+AS(I,M)
   51 CONTINUE
      RETURN
      END