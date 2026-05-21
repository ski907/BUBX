
C  LATEST PROGRAM MODIFICATIONS BY JCT- OCTOBER 1988
C
C  PORTED TO 300 SERIES AUG 94 ,ZABILANSKY
C  PORTED TO FORTRAN 77 MARCH 1997, R.B. HAEHNEL AND T. KERR
C
C Bug fixes to Program 17 Nov 2015, R. Haehnel
C
C  DIRECT QUESTIONS REGARDING THIS CODE TO:
C	ROBERT B. HAEHNEL
C	RESEARCH MECHANICAL ENGINEER
C	US ARMY COLD REGIONS RESEARCH AND ENGINEERING LABORATORY
C	72 LYME ROAD 
C	HANOVER, NH 03755
C	(603)646-4325
C	rhaehnel@crrel.usace.army.mil
C
C  PROGRAM GIVEN COMPRESSOR DISCHARGE VOLUME AND PRESSURE 'GVQCPC'
C
C  PROGRAM IS INTENDED FOR ADJUSTING PARAMETERS
C
C  PROGRAM IS AN IMPLIMANTATION OF G.ASHTON ITERATIVE SCHEME
C  FOR DETERMINING AIR DISCHARGE RATE AT EACH ORFICE. IN THE PROGRAM
C  THE USER IS AT LIBERITY TO USES VARIOUS PIPE SIZES AND
C  LENGTHS IN THE SUPPLY MANUFOLD. DIFFRENT PIPE SIZES AND
C  ORFICE SPACING CAN ALSO BE USED IN THE DIFFUSSER LINE.
C  USING AN INTERACTIVE APPROCH THE USER CAN VARY A SINGLE ELEMENT
C  OR THE ENTIRE PIPING LAYOUT WHILE EVALUATING VARIOUS PIPING
C  LAYOUTS.
C
C  ROUTINE STARTS WITH A GIVEN COMPRESSOR DISCHARGE VOLUME
C  AND PRESSURE.  WORKING DOWN THE PIPING THE PRESSURE DROPS
C  ARE DEDUCTED AT EACH PIPE CONNECTION AND THE AIR DISCHARGE
C  OF THE ORFICE SUBTRACTED FOR THE AIR AVAILABLE.
C
C  WRITTEN BY:  L.ZABILANSKY
C  23,MARCH  84
C  UPDATE:      4 APRIL
C
C  VARIABLE       USAGE
C    Aorf         AREA OF ORIFICE
C    Cdis         DISCHARGE COEFFICIENT FOR THE ORIFICE
C                 DEFAULT VALUE IS 0.63 FOR A SQUARE EDGE
C    Cgrv         GRAVITATIONAL CONSTANT (9.81 M^2/S)
C    Denair0      DENSITY OF AIR AT 0oC, 1 ATM (1.294 Kg/M^3)
C    Denair       DENSITY OF AIR AT COMPRESSOR PRESSURE
C    Den          DENSITY OF AIR AT LOCAL PIPE PRESSURE
C    Deno         DENSITY OF AIR AT SUBMERGENCE PRESSURE
C    Denwtr       DENSITY OF WATER (1000 Kg/M^3)
C    Diff(I,J)    DIFFUSER LAYOUT
C                 I SECTION ID
C                   I = Imax IS USED FOR TOTAL LENGTH
C                 J 1 = LENGTH OF PIPE
C                   2 = DIAMETER OF PIPE
C                   3 = AREA OF PIPE
C    Difftyp      DEFINE TYPE OF DIFFUSER
C                   1= CONSTANT DIFFUSER DIAMETER WITH EQUALLY SPACED
C                        ORIFICES
C                   2= VARYING PIPE DIA AND/OR ORIFICE SPACING
C    Dorf         DIAMETER OF ORIFICES
C    Exp          1: ADIABATIC, 2:ISOTHERMAL AIR EXPANISION AT NOZZLES
C    F            FRICTION FACTOR (ROUSE 1946)
C    Hl           TOTAL HEAD LOSS IN SUPPLY PIPE
C    Hwtr         WATER DEPTH
C    I            FOR- NEXT COUNTER
C    Imax         MAX NUMBER OF ORIFICES IN DIFFUSER PIPE AND NUMBER OF
C                   PIPE SECTIONS OF DIFFERENT DIAMETER IN SUPPLY HEADER
C    K            PIPE SECTION No.
C    Kr           PIPE ROUGHNESS (MATERIAL)
C    Kage         0=NEW SYTEM; 1=OLD SYSTEM
C    Lockstr      LOCK NAME
C    Matestr      PIPE MATERIAL (SAME FOR SUPPLY AND DIFFUSER LINES)
C    Nodia        NUMBER OF PIPE SIZES IN SUPPLY HEADER
C    Noorf        NUMBER OF ORIFICES IN DIFFUSER PIPE
C    P(K)         PRESSURE DIFFERENCE BETWEEN INSIDE AND OUTSIDE OF
C                   DIFFUSER
C    Pc           NOMINAL COMPRESSOR PRESSURE
C    Pcalc        CALCULATE PRESSURE
C    Pdrpd(K)     PRESSURE DROP ASSOCIATED W/ LENGTH K IN DIFFUSER PIPE
C    Pdrps(K)     PRESSURE DROP ASSOCIATED W/ LENGTH K IN SUPPLY PIPE
C    Pexc         DEAD END PRESSURE
C    Pwtr         HYDROSTATIC WATER PRESSURE
C    Qcin         RATED DISCHARGE OF COMPRESSOR
C    Qc           UPDATED COMPRESSOR DISCHARGE
C    Qtot         CALCULATED DISCHARGE
C    Supp(I,J)    ARRAY FOR SUPPLY PIPE HEADER
C                 I= PIPE ID
C                   I= Imax IS USED FOR TOTAL LENGTH
C                 J 1= LENGTH
C                   2= DIAMETER
C                   3= AREA
C    Tpdrp        TOTAL PRESSURE DROP IN SUPPLY PIPE
C    Typ          TYPE OF BUBBLER SYSTEM BEING DESIGNED
C    Typesstr(*)  ARRAY OF TYPES OF BUBBLER SYSTEM
C    Unit         UNITS USED FOR INPUT
C                   1= ENGLISH
C                   2= S-I (CALCULATIONS ARE DONE IN S-I)
C    Vel          VELOCITY IN PIPE
C    Visair       DYNAMIC VISCOSITY OF AIR AT 0oC
C    Q(K)         VOL OF AIR DISCHARGE IN ORIFICE K
C    Qc           COMP VOL
C    Qperl        VOL OF AIR DISCHARGE PER UNIT LENGTH OF PIPE
C    Qsum         TOTAL VOL OF AIR DISCHARGED

      integer endpt
      real FNArea
      external endpt, FNArea

      real P(200), Pdrpd(200), Pdrps(200), Q(200), Qhyd(200)
      real Diff(200,3), Supp(200,3)
      real Kr, Mc2, Mout, Msum, Qsum
      integer Choice, Conti, Difftyp, Exp, Opt, Run, Typ, Unit, Noorf
      character Lockstr*80,Matestr(6)*50, Pstr, Pausestr,Typesstr(4)*30
      character Ystr
      character filename*12, fni*8, fn*8, vdate*13
C  INITIALIZE CONSTANTS
      PARAMETER (Imax=200, Cgrv=9.81, Visair=1.71e-5, Denair0=1.294,
     *     Denwtr=1000, Patm=98100)
      Noorf=0                   ! initialize number of orifaces
10    Cdis=1
      Run=0
      fni='BubblOut'
C Program Version date
      vdate = '9 Sept 2005'
C TYPES OF BUBBLERS
      Typesstr(1)='Deflector screen'
      Typesstr(2)='High velocity screen'
      Typesstr(3)='Gate recess'
      Typesstr(4)='Deicing bubbler'
C PIPE MATERIAL
      Matestr(1)='New Stainless Steel'
      Matestr(2)='Old Stainless Steel'
      Matestr(3)='New Galvanized Steel'
      Matestr(4)='Old Galvanized Steel'
      Matestr(5)='PVC Pipe'

C Program start-up header
      write(*,*)' '
      write(*,*)' '
      write(*,*)' '
      write(*,*)' '
      write(*,*)' '
      write(*,*)' *****************************************************'
      write(*,*)' *                                                   *'
      write(*,*)' *            BUBBLER DESIGN PROGRAM                 *'
      write(*,*)' *                                                   *'
      write(*,*)' *****************************************************'
      write(*,*)' '
      write(*,*)' '
      write(*,*)' '
      write(*,*)' Developed by:'
      write(*,*)'  US Army Corps of Engineers'
      write(*,*)'  Cold Regions Research and Engineering Laboratory'
      write(*,*)'  Ice Engineering Research Division'
      write(*,*)'  Hanover, NH'
      write(*,101)vdate
      write(*,*)' '
      write(*,*)' '
101   format('   Version = ',a)
      
C  SET-UP OUTPUT FILE TO WRITE INPUT PARAMETERS AND ANALYSIS RESULTS
      PRINT *, 'Enter name of output file for Bubbler Design Data'
      PRINT '(A37, A8, A1)', 'or press <ENTER> to accept default (', fni
     *, ')'
      READ '(A8)', fn
      IF (fn.ne.'        ') then
        fni=fn
      END IF
      filename=fni(:endpt(fni))//'.txt'
      OPEN (1, FILE=filename)
C  MAIN PROGRAM
      CALL Menu(Typesstr, Matestr, Typ, Lockstr, Km, Kr,Run,Unit,Opt)
      Run=1
      Opt=0
C  INPUT INITIAL CONDITIONS
      call Diffcon(Difftyp)
c      write(*,*)'***Difftype', Difftype
      IF (Difftyp.EQ.1) THEN
        CALL Lendiff(Unit, Imax, Noorf, Diff)
        CALL Diadiff(Unit, Imax, Diff)
        CALL Orfspc(Unit, Imax, Noorf, Diff)
      ELSE
        Diff(1,1)=0
        CALL Diffpipe(Unit, Imax, Noorf, Diff)
      END IF
      CALL Orfdia(Unit, Dorf, Aorf)
      PRINT *, 'ENTER NUMBER OF PIPE SIZES IN SUPPLY LINE'
      READ *, Nodia
      IF (Nodia.EQ.1) THEN
        CALL Lensupp(Unit, Imax, Supp)
        CALL Diasupp(Unit, Imax, Supp)
      ELSE
        CALL Supppipe(Unit, Imax, Nodia, Supp)
      END IF
      CALL Comp(Unit, Pc, Qcin, Pstr)
      Denair=Denair0*(Pc+Patm)/Patm
      Qc=Qcin
      IF (Pstr.EQ.'A') Qc=Qcin*Denair0/Denair
      CALL H20dpth(Unit, Hwtr)
      CALL Expansion(Exp)
20    PRINT *
      PRINT *, ' 1  SYSTEM TYPE, MATERIAL, OR CONFIGURATION'
      PRINT *, ' 2  DIFFUSER LENGTH:   EQUALLY SPACED NOZZLES, CONSTANT 
     *DIAMETER'
      PRINT *, ' 3  DIFFUSER DIAMETER: EQUALLY SPACED NOZZLES, CONSTANT 
     *DIAMETER'
      PRINT *, ' 4  DIFFUSER WITH UNEQUAL NOZZLE SPACINGS OR PIPE DIAMET
     *ER'
      PRINT *
      PRINT *, ' 5  NOZZLE DIAMETER'
      PRINT *, ' 6  NOZZLE SPACING'
      PRINT *, ' 7  NOZZLE DISCHARGE COEFFICIENT'
      PRINT *
      PRINT *, ' 8  NUMBER OF SECTIONS IN SUPPLY LINE WITH DIFFERENT DIA
     *METER'
      PRINT *, ' 9  LENGTH OF SUPPLY LINE (CONSTANT DIAMETER)'
      PRINT *, '10  DIAMETER OF SUPPLY LINE (CONSTANT DIAMETER)'
      PRINT *, '11  CONFIGURATION OF SUPPLY LINE'
      PRINT *
      PRINT *, '12  RATED COMPRESSOR PRESSURE AND DISCHARGE'
      PRINT *, '13  DEPTH OF SUBMERGENCE OF DIFFUSER'
      PRINT *, '14  CHOICE OF EXPANSION PROCESS AT NOZZLE'
      PRINT *, '15  PRINT INPUT PARAMETERS TO OUTPUT FILE'
      PRINT *, '16  CALCULATION OF DIFFUSER PERFORMANCE (Results output
     * to file)'
      PRINT *, '17  NEW RUN'
      PRINT *, '18  EXIT PROGRAM'
      PRINT *
      PRINT *, 'SELECT OPTION FROM ABOVE MENU'
      READ *, Choice
      IF ((Choice.LT.1).OR.(Choice.GT.18)) GOTO 20
C  DEFINE TYPE OF DIFFUSER
      IF (Choice.EQ.1) THEN
        CALL Menu(Typesstr, Matestr, Typ, Lockstr,Km,Kr,Run,Unit,Opt)
        IF ((Opt.EQ.1).OR.(Opt.EQ.3)) call Diffcon(Difftyp)
        GOTO 20
C  DEFINE LENGTH OF DIFFUSER
      ELSEIF (Choice.EQ.2) THEN
        IF (Difftyp.EQ.1) CALL Lendiff(Unit, Imax, Noorf, Diff)
C  NEED TO DEFINE LAYOUT
        IF (Difftyp.NE.1) CALL Diffpipe(Unit, Imax, Noorf, Diff)
        GOTO 20
C  DEFINE DIAMETER OF DIFFUSER PIPE
      ELSEIF (Choice.EQ.3) THEN
        CALL Diadiff(Unit, Imax, Diff)
        GOTO 20
C  DIFFUSER GEOMETRY IS VARIABLE
      ELSEIF (Choice.EQ.4) THEN
        CALL Diffpipe(Unit, Imax, Noorf, Diff)
        GOTO 20
C  DEFINE DIAMETER OF ORFICES
      ELSEIF (Choice.EQ.5) THEN
        CALL Orfdia(Unit, Dorf, Aorf)
        GOTO 20
C  DEFINE ORFICE SPACING
      ELSEIF (Choice.EQ.6) THEN
        IF (Difftyp.EQ.1) CALL Orfspc(Unit, Imax, Noorf, Diff)
        IF (Difftyp.NE.1) CALL Diffpipe(Unit, Imax, Noorf, Diff)
        GOTO 20
C  CHANGE DISCHARGE COEFFICIENT
      ELSEIF (Choice.EQ.7) THEN
        PRINT *, 'ENTER DISCHARGE COEFFICIENT (DEFAULT ', Cdis, ')'
        READ *, Cdis
        GOTO 20
C  DEFINE SUPPLY HEADER
      ELSEIF (Choice.EQ.8) THEN
        PRINT *, 'ENTER NUMBER OF PIPE SIZES IN SUPPLY HEADER'
        READ *, Nodia
        GOTO 20
C  DEFINE LENGTH OF SUPPLY HEADER
      ELSEIF (Choice.EQ.9) THEN
        IF (Nodia.EQ.1) CALL Lensupp(Unit, Imax, Supp)
        IF (Nodia.NE.1) CALL Supppipe(Unit, Imax, Nodia, Supp)
        GOTO 20
C  DEFINE DIAMETER OF SUPPLY PIPE
      ELSEIF (Choice.EQ.10) THEN
        CALL Diasupp(Unit, Imax, Supp)
        GOTO 20
C  DEFINE SUPPLY HEADER PIPING
      ELSEIF (Choice.EQ.11) THEN
        CALL Supppipe(Unit, Imax, Nodia, Supp)
        GOTO 20
C  INPUT NOMINAL COMPRESS PRESSURE
      ELSEIF (Choice.EQ.12) THEN
        CALL Comp(Unit, Pc, Qcin, Pstr)
        Denair=Denair0*(Pc+Patm)/Patm
        Qc=Qcin
        IF (Pstr.EQ.'A') Qc=Qcin*Denair0/Denair
        GOTO 20
C  DEFINE DEPTH OF SUBMERGENCE
      ELSEIF (Choice.EQ.13) THEN
        CALL H20dpth(Unit, Hwtr)
        GOTO 20
C  DEFINE THE EXPANSION PROCESS AT NOZZLE
      ELSEIF (Choice.EQ.14) THEN
        CALL Expansion(Exp)
        GOTO 20
C  PRINT INPUT DATA
      ELSEIF (Choice.EQ.15) THEN
        WRITE (1, '(32X, A8, A30)') 'PROJECT ', Lockstr
        WRITE (1, '(22X, A16, A30)') 'TYPE OF SYSTEM: ', Typesstr(Typ)
        WRITE (1, *)
        WRITE (1, '(34X, A16)') 'INPUT PARAMETERS'
        WRITE (1, *)
        WRITE (1, *) '** SUPPLY PIPE PARAMETERS'
        WRITE (1, *)
        IF (Nodia.EQ.1) THEN
          IF (Unit.EQ.1) THEN
            WRITE (1, 30) 'TOTAL LENGTH: ', Supp(Imax,1)/0.3048, ' ft'
30          FORMAT (4X, A14, 10X, F5.1, A3)
            WRITE (1, 40) 'DIAMETER: ', Supp(Imax,2)/0.0254,' in'
40          FORMAT (4X, A10, 14X, F5.1, A3)
          ELSE
            WRITE (1, 30) 'TOTAL LENGTH: ', Supp(Imax,1),' M'
            WRITE (1, 40) 'DIAMETER: ', Supp(Imax,2)*100.,' cm'
          END IF
        ELSE
          IF (Unit.EQ.1) WRITE (1, 50) 'SECT. #','LENGTH (ft)', 
     +'DIAMETER (in)'
          IF (Unit.EQ.2) WRITE (1, 50) 'SECT. #','LENGTH (M)', 
     +'DIAMETER (cm)'
50        FORMAT (2X, A10, A20, 2X, A20)
          DO 70, I=1, Nodia
            IF (Unit.EQ.1) WRITE (1, 60) I, Supp(I,1)/0.3048, Supp(I,2)/
     *0.0254
            IF (Unit.EQ.2) WRITE (1, 60) I, Supp(I,1), Supp(I,2)*100.
60          FORMAT (7X, I3, 14X, F5.1, 16X, F5.2)
70        CONTINUE
          IF (Unit.EQ.1) WRITE (1, 80) 'TOTAL LENGTH:', Supp(Imax,1)/
     +0.3048,' ft'
          IF (Unit.EQ.2) WRITE (1, 80) 'TOTAL LENGTH:', Supp(Imax,1), 
     +'M'
80        FORMAT (4X, A13, 7X, F5.1, A3)
        END IF
        WRITE (1, *)
        WRITE (1, *) '** DIFFUSER LINE PARAMETERS'
        WRITE (1, *)
        WRITE (1, '(4X, A20, 7X, I3)') 'NUMBER OF ORIFICES: ', Noorf
90      FORMAT (4X, A17, 10X, F5.3, A3)
        IF (Unit.EQ.1) WRITE (1, 90) 'ORIFICE DIAMETER:', Dorf/0.0254, '
     * in'
        IF (Unit.EQ.2) WRITE (1, 90) 'ORIFICE DIAMETER:', Dorf*100, ' cm
     *'
        WRITE (1, '(4X, A22, 5X, F5.3)') 'DISCHARGE COEFFICIENT:',
     * Cdis
        WRITE (1, *)
        IF (Unit.EQ.1) WRITE (1, 100) 'SECT. #', 'SPACING (ft)', 
     +'DIAMETER (in)'
        IF (Unit.EQ.2) WRITE (1, 100) 'SECT. #', 'SPACING (M)', 
     +'DIAMETER (cm)'
100     FORMAT (2X, A10, A20, 2X, A20)
110     FORMAT (7X, I3, 15X, F4.1, 16X, F5.2)
        DO 120, I=1, Noorf-1
          IF (Unit.EQ.1) WRITE (1, 110) I, Diff(I,1)/0.3048, Diff(I,2)/
     + 0.0254
          IF (Unit.EQ.2) WRITE (1, 110) I, Diff(I,1), Diff(I,2)*100
120     CONTINUE
        IF (Unit.EQ.1) WRITE (1, 130) 'TOTAL LENGTH:', Diff(Imax,1)/
     + 0.3048,' ft'
        IF (Unit.EQ.2) WRITE (1, 130) 'TOTAL LENGTH:', Diff(Imax,1),' M'
130     FORMAT (4X, A13, 7X, F5.1, A)
        WRITE (1, *)
        IF (Km.EQ.0) THEN
          WRITE (1, *) 'ALL PIPES ARE "HYDRAULICALLY SMOOTH"'
        ELSE
          IF (Unit.EQ.1) WRITE (1, 140) 'PIPE MATERIAL: ', Matestr(Km),
     *' of ROUGHNESS K= ', Kr/0.0254, ' in'
          IF (Unit.EQ.2) WRITE (1, 140) 'PIPE MATERIAL: ', Matestr(Km),
     *' of ROUGHNESS K= ', Kr*1000, ' mm'
140       FORMAT (A16, A20, A17, F5.3, A3)
        END IF
        WRITE (1, *)
        WRITE (1, *) '** COMPRESSOR RATINGS'
        IF (Unit.EQ.1) THEN
          WRITE (1, 150) 'RATED COMPRESSOR PRESSURE:', INT(Pc/6894.757+
     + 0.5), ' psi'
          IF (Pstr.EQ.'A') WRITE (1, 160) 'RATED COMPRESSOR DISCHARGE:', 
     *INT(Qcin*2118.9+0.5), ' CFM at ATMOSPHERIC PRESSURE'
          IF (Pstr.EQ.'P') WRITE (1, 170) 'RATED COMPRESSOR DISCHARGE:', 
     *INT(Qcin*2118.9+0.5), ' CFM at RATED PRESSURE'
          WRITE (1, *)
          WRITE (1, '(4X, A21, 9X, I3, A3)') 'DEPTH OF SUBMERGENCE:',
     * INT(Hwtr/0.3048+0.5), ' ft'
          WRITE (1, '(A32, F5.3, A9)') 'AIR DENSITY AT RATED PRESSURE: '
     *, Denair/16.02, ' lbm/ft^3'
        ELSE
          WRITE (1, 150) 'RATED COMPRESSOR PRESSURE:', INT(Pc/1000+0.5),
     * ' kPa'
          IF (Pstr.EQ.'A') WRITE (1, 160) 'RATED COMPRESSOR DISCHARGE:', 
     *INT(Qcin*1000+0.5), ' L/S at ATMOSPHERIC PRESSURE'
          IF (Pstr.EQ.'P') WRITE (1, 170) 'RATED COMPRESSOR DISCHARGE:', 
     *INT(Qcin*1000+0.5), ' L/S at RATED PRESSURE'
          WRITE (1, *)
          WRITE (1, '(4X, A21, 8X, F4.1, A2)') 'DEPTH OF SUBMERGENCE:',
     * Hwtr, ' M'
          WRITE (1, '(A32, F7.3, A7)') 'AIR DENSITY AT RATED PRESSURE: '
     * , Denair, ' kg/M^3'
        END IF
150     FORMAT (4X, A26, 2X, I5, A4)
160     FORMAT (4X, A27, 1X, I5, A28)
170     FORMAT (4X, A27, 1X, I5, A22)
        WRITE (1, *)
        IF (Exp.EQ.0) WRITE (1, *) 'AIR IS CONSIDERED INCOMPRESSIBLE THR
     *OUGHOUT'
        IF (Exp.GT.0) WRITE (1, *) 'AIR FLOW THROUGH PIPES IS ISOTHERMAL
     * @ 0oC'
        IF(Exp.EQ.1) WRITE(1, *)'AIR EXPANSION AT NOZZLES IS ADIABATIC'
        IF(Exp.EQ.2) WRITE(1, *)'AIR EXPANSION AT NOZZLES IS ISOTHERMAL'
        WRITE (1, *)
        GOTO 20

C     PROCESS DATA
      ELSEIF (Choice.EQ.16) THEN
        Pwtr=Denwtr*Cgrv*Hwtr
        IF (Pwtr.GT.Pc) THEN
          PRINT *, 'COMPRESSOR PRESSURE IS INADEQUATE TO OVERCOME HYDROS
     *TATIC PRESSURE'
          PRINT *, 'Press ENTER when ready'
          READ '(A)', Pausestr
          GOTO 20
        END IF
C  CHECK ORIFICE AREA VS PIPE AREA
        Areq=0
        DO 190, I=1, Noorf-1
C     IS AREA OF DIFFUSER LESS THAN DOWN STREAM ORIFICE AREA?
          Dnsta=(Noorf-I)*Aorf
          IF (Diff(I,3).LT.Dnsta) THEN
C  DIFFUSER PIPING IS UNDERSIZED
            PRINT *, 'TOTAL ORIFICE AREA DOWNSTREAM FROM ORIFICE ', I,
     * ' EXCEEDS AREA OF DIFFUSER PIPE'
            PRINT *, 'Joint ', I,'  Diffuser X sect  ', Diff(I,3),
     * '  Req. Orfice area ', Dnsta
            PRINT *, 'Press ENTER when ready'
            READ '(A)', Pausestr
            GOTO 20
          ELSEIF (Diff(I,3).LT.(4*Dnsta)) THEN
            PRINT *
            PRINT '(A8, I3, A44)', 'SECTION', I,' VIOLATES RULE FOR UNIF
     *ORM FLOW DISTRIBUTION'
            write(*,*)
            write(*,*)' Either orfice diameter needs to be reduced or' 
            write(*,*)'  diffuser diameter needs to be increased'
            write(*,*)
            PRINT *,' ENTER  1  IF YOU WANT TO CONTINUE CALCULATIONS'
            PRINT *,'   OR   2  IF YOU WANT TO MODIFY SYSTEM PARAMETERS'
180         PRINT *, 'YOUR CHOICE? (1 or 2)'
            READ *, Conti
            IF ((Conti.NE.1).AND.(Conti.NE.2)) GOTO 180
            IF (Conti.EQ.2) GOTO 20
            write(*,*)
          END IF
190     CONTINUE
C   ASSUME SIZE OF SUPPLY LINE IS THE SAME AS THE FIRST LENGTH OF
C     DIFFUSER PIPE
        Dnsta=Dnsta+Aorf
C  DETERMINE PRESSURE DROP IN SUPPLY PIPE
200     PRINT *, 'DO YOU WISH TO SEE INTERMEDIATE CALCULATIONS? (Y/N)'
        READ '(A)', Ystr
        IF (Ystr.EQ.'y') Ystr='Y'
        IF (Ystr.EQ.'n') Ystr='N'
        IF ((Ystr.NE.'Y').AND.(Ystr.NE.'N')) GOTO 200
210    Hl=0
        Mc2=Qc*Denair
        Qtot=0
        Den=Denair
        Msum=Mc2
        Pcalc=Pc
C  DETERMINE PRESSURE DROPS IN SUPPLY LINE STARTING AT THE SOURCE
        DO 220, I=1, Nodia
          P1=Pcalc
          Vel=Msum/Den/Supp(I,3)
          Rey=Den*Vel*Supp(I,2)/Visair
          CALL Fric(Rey, Supp(I,2), Kr, F)
          Delpi=F*Supp(I,1)*Den*Vel**2/(2*Supp(I,2))
          P2=P1-Delpi
          Delpc=Delpi*2.*(1+2*Supp(I,2)/Supp(I,1)/F*LOG(P1/P2)/
     +          (1+P2/P1))
C  PRESSURE DROP FOR COMPRESSIBLE FLOW
          Pdrps(I)=Delpc
          Hl=Hl+Pdrps(I)
          Pcalc=Pcalc-Pdrps(I)
          IF (Exp.NE.0) Den=Denair*(Pcalc+Patm)/(Pc+Patm)
220     CONTINUE
C  DETERMINE PRESSURE DROPS IN DIFFUSER
C  STARTING AT THE SOURCE END REDUCING AIR VOL BY THE ORIFICE
C  DISCHARGE AND PRESSURE BY THE PRESSURE DROP
        K=1
C  DETERMINE DELTA PRESSURE AT EACH ORIFICE
230     P(K)=Pcalc-Pwtr
        IF (P(K).LE.0) THEN
          PRINT '(A36, I3)', 'RAN OUT OF AIR PRESSURE AT ORIFICE:', K
          Qc=Qc/2.
          PRINT '(A23, F6.4)', 'REDUCED AIR TO SUPPLY:', Qc
          GOTO 210
        END IF
        IF (Exp.EQ.0) Deno=Den
        IF (Exp.EQ.1) Deno=Den*((Pwtr+Patm)/(Pcalc+Patm))**(1./1.4)
        IF (Exp.EQ.2) Deno=Denair0*(Pwtr+Patm)/Patm
C  Deno IS AIR DENSITY AT DEPTH OF SUBMERGENCE
        Vout=SQRT(P(K)*2/Deno)
        Qhyd(K)=Cdis*Aorf*Vout
        Mout=Deno*Qhyd(K)
        Q(K)=Qhyd(K)*Deno/Denair
        Msum=Msum-Mout
        Qtot=Qtot+Q(K)
        IF (Ystr.EQ.'Y') THEN
          PRINT '(A11, I3)', 'NOZZLE #: ', K
          PRINT '(A16, F6.1, A5)', '     EXIT VEL: ', Vout/0.3048, 
     +' ft/s'
          PRINT '(A31, F5.3, A9)', '     AIR DENSITY IN DIFFUSER: ',
     * Den*0.3048**3/0.453, ' lbm/ft^3'
          PRINT '(A34, F4.1, A4)', '     DISCHARGE at COMP.PRESSURE: ', 
     *Q(K)*2118.9, ' CFM'
          PRINT '(A16, F5.3, A4)', '     MASS OUT: ', Mout/0.453, ' lbm'
          PRINT '(A22, F5.3, A4)', '     REMAINING MASS: ', Msum/0.453, 
     *' lbm'
          PRINT '(A26, F5.1, A4)', '     SUBTOTAL DISCHARGE: ', 
     *Qtot*2118.9, ' CFM'
          PRINT *
          PRINT *, 'Press ENTER when ready'
          READ '(A)', Pausestr
        END IF
        IF (Msum.LE.0) THEN
          PRINT '(A27, I2, A24, F5.3, A16, F6.4)', 
     +'RAN OUT OF AIR AT ORIFICE ', K, ' CURRENT COMPRESSOR VOL ', Qc, 
     +' TOTAL MASS AIR ', Msum
C  Qc=1.1*Qtot*Noorf/K
          Qc=1.1*Qc
          IF (Unit.EQ.2) PRINT '(A26, F5.1, A6)', 
     +'NEW COMPRESSOR DISCHARGE ', Qc, ' M^3/S'
          IF (Unit.EQ.1) PRINT '(A26, F5.1, A4)', 
     +'NEW COMPRESSOR DISCHARGE ', Qc*2118.9, ' CFM'
          GOTO 210
        END IF
        IF (K.NE.Noorf) THEN
C  DETERMINE PRESSURE DROP IN DIFFUSER PIPE BETWEEN ORIFICES
C  STARTING AT THE SOURCE END REDUCING AIR VOL BY THE ORIFICE
C  DISCHARGE AND PRESSURE BY THE PRESSURE DROP
          P1=Pcalc
          Vel=Msum/Den/Diff(K,3)
          Rey=Den*Vel*Diff(K,2)/Visair
          CALL Fric(Rey, Diff(K,2), Kr, F)
          Delpi=F*Diff(K,1)*Den*Vel**2./(2*Diff(K,2))
          P2=P1-Delpi
          IF (Exp.GT.0) Delpc=Delpi*2.*(1+2*Diff(K,2)/Diff(K,1)/F*
     +         LOG(P1/P2)/(1+P2/P1))
          IF (Exp.EQ.0) Delpc=Delpi
C  PRESSURE DROP FOR COMPRESSIBLE FLOW
          Pdrpd(K)=Delpc
          Pcalc=Pcalc-Pdrpd(K)
          IF (Exp.GT.0) Den=Denair*(Pcalc+Patm)/(Pc+Patm)
          K=K+1
          GOTO 230
        END IF
C  CHECK CALCULATED DISCHARGE AGAINST ASSUMED DISCHARGE
        Qsum=Qc-Qtot
        IF ((Qsum*2118.9).GT.15) THEN
C          Qc=0.95*Qc
          Qc=Qc-0.3*Qsum
          IF (Unit.EQ.1) PRINT '(A28, F5.1, A23, F5.1, A4)', 
     +'EXCESS COMP DISCHARGE (CFM)', Qsum*2118.9, 
     +'. DISCHARGE REDUCED TO ', Qc*2118.9, ' CFM'
          IF (Unit.EQ.2) PRINT '(A30, F6.4, A23, F5.4, A6)', 
     +'EXCESS COMP DISCHARGE (M^3/S)', Qsum,'. DISCHARGE REDUCED TO ', 
     +Qc, ' M^3/S'
          GOTO 210
        END IF
C  PRINT OUT RESULTS
        WRITE (1, '(24X, A31)') ' CALCULATED SYSTEM PERFORMANCE '
        WRITE (1, *)
        IF (Nodia.GT.1) THEN
          WRITE (1, *) '** PRESSURE LOSSES IN SUPPLY LINE'
          WRITE (1, '(11X, A10, 10X, A13)') 'SECTION #','PRESS. DROP'
          IF (Unit.EQ.1) WRITE (1, 250) 'psi'
          IF (Unit.EQ.2) WRITE (1, 250) 'kPa'
          DO 240, I=1, Nodia
            IF (Unit.EQ.1) WRITE (1, 260) I, Pdrps(I)/6894.757
            IF (Unit.EQ.2) WRITE (1, 260) I, Pdrps(I)/1000
240       CONTINUE
          WRITE (1, *)
        END IF
250     FORMAT (36X, A4)
260     FORMAT (14X, I2, 20X, F4.1)
        WRITE (1, *) '** DIFFUSER: PRESSURE AND FLOW DISTRIBUTION'
        WRITE (1, *)
        IF (Unit.EQ.1) THEN
          WRITE (1, 270) 'LOCATION', 'PRESS. DROP', 'Pin-Pout', 'ORIFICE
     * DISCHARGE (CFM)'
          WRITE (1, 280) 'ft', 'psi', 'psi', 'at', INT(Pc/6894.757+0.5),
     * 'psi', 'at Pdepth', 'at Patm'
        ELSE
          WRITE (1, 270) 'LOCATION', 'PRESS. DROP', 'Pin-Pout',
     * 'ORIFICE DISCHARGE (L/S)'
          WRITE (1, 280) 'M', 'kPa', 'kPa', 'at', INT(Pc/1000+0.5), 'kPa
     *', 'at Pdepth', 'at Patm'
        END IF
270     FORMAT (A9, A14, A10, 7X, A23)
280     FORMAT (4X, A2, 10X, A3, 9X, A3, 5X, A2, 1X, I3, 1X, A3, 5X, A9,
     * 5X, A7)
        Shl=0
        WRITE (1, *)
        Qsum=0
        DO 290, I=1, Noorf
          IF (I.EQ.1) THEN
            X=0
            IF (Unit.EQ.1) WRITE (1, 300) X/0.3048, P(I)/6894.757,
     * Q(I)/4.719E-4, Qhyd(I)/4.719E-4, Q(I)*(1+Pc/Patm)/4.719E-4
            IF (Unit.EQ.2) WRITE (1, 300) X, P(I)/1000, Q(I)*1000,
     * Qhyd(I)*1000, Q(I)*(1+Pc/Patm)*1000
            Dhl=0
          ELSE
            X=Diff(I-1,1)+X
            IF (Unit.EQ.1) WRITE (1, 310) X/0.3048, Pdrpd(I-1)/6894.757,
     * P(I)/6894.757, Q(I)/4.719E-4, Qhyd(I)/4.719E-4,
     * Q(I)*(1+Pc/Patm)/4.719E-4
            IF (Unit.EQ.2) WRITE (1, 310) X, Pdrpd(I-1)/1000, P(I)/1000,
     * Q(I)*1000, Qhyd(I)*1000, Q(I)*(1+Pc/Patm)*1000
            Dhl=Dhl+Pdrpd(I-1)
         END IF
         Qsum=Qsum+Qhyd(I)
290     CONTINUE
300     FORMAT (1X, F5.1, 20X, F6.1, 6X, F7.1, 6X, F7.1, 8X, F7.1)
310     FORMAT (1X, F5.1, 9X, F6.1, 5X, F6.1, 6X, F7.1, 6X, F7.1, 8X, F7
     *.1)
        WRITE (1, *)
        IF (Unit.EQ.1) THEN
          WRITE (1, 320) 'TOTAL PRESSURE DROP IN SUPPLY LINE:', Hl/6894.
     *757, ' psi'
          WRITE (1, 330) 'TOTAL PRESSURE DROP IN DIFFUSER:', Dhl/6894.75
     *7, ' psi'
          WRITE (1, 340) 'TOTAL PRESSURE DROP IN SYSTEM:', (Dhl+Hl)/6894
     *.757, ' psi'
        ELSE
          WRITE (1, 320) 'TOTAL PRESSURE DROP IN SUPPLY LINE:', Hl/1000 
     *, ' kPa'
          WRITE (1, 330) 'TOTAL PRESSURE DROP IN DIFFUSER:', Dhl/1000, '
     * kPa'
          WRITE (1, 340) 'TOTAL PRESSURE DROP IN SYSTEM:', (Dhl+Hl)/1000
     *, ' kPa'
        END IF
320     FORMAT (A36, 5X, F6.1, A4)
330     FORMAT (A33, 8X, F6.1, A4)
340     FORMAT (A31, 10X, F6.1, A4)
        WRITE (1, *)
C  DETERMINE EXCESS END OF DIFFUSER PRESSURE
        Pexc=Pcalc-Pwtr
C  DETERMINE EXCESS DISCHARGE
        IF (Unit.EQ.1) THEN
          WRITE (1, 350) 'PRESSURE AT DIFFUSER END:', Pcalc/6894.757, ' 
     +psi'
          WRITE (1, 360) 'HYDROSTATIC PRESSURE:', Pwtr/6894.757, ' psi'
          IF (Pstr.EQ.'P') WRITE (1, 370) 'CALCULATED COMPRESSOR DISCHAR
     + GE AT ', INT(Pc/6894.757+0.5), ' psi: ', 
     + INT(Qtot*2118.9+0.5),' CFM'
          IF (Pstr.EQ.'A') WRITE (1, 380) 'CALCULATED COMPRESSOR DISCHAR
     +GE AT ATMOSPHERIC PRESSURE:', INT(Qtot*2118.9*Denair/Denair0+0.5),
     + ' CFM'
          IF (Pstr.EQ.'P') WRITE (1, 380) 'RATED minus CALCULATED DISCHA
     +RGE (at RATED Pressure) ', INT((Qcin-Qtot)*2118.9+0.5), ' CFM'
          IF (Pstr.EQ.'A') WRITE (1, 380) 'RATED minus CALCULATED DISCHA
     +RGE (at ATM Pressure) ', 
     +         INT((Qcin-Qtot*Denair/Denair0)*2118.9+0.5),' CFM'
          write(1,*)
          write(1,*)'Total Discharge from nozzles at depth ',
     +     Qsum/4.719E-4,' CFM'
          write(1,*)'Average nozzle discharge ',Qsum/Noorf/4.719E-4,
     +     ' CFM'
          write(1,*)'Coefficient of uniformity ',
     +        (Qhyd(1)-Qhyd(Noorf))/(Qsum/Noorf)
        ELSE
          WRITE (1, 350) 'PRESSURE AT DIFFUSER END:', Pcalc/1000, ' kPa'
          WRITE (1, 360) 'HYDROSTATIC PRESSURE:', Pwtr/1000, ' kPa'
          IF (Pstr.EQ.'P') WRITE (1, 370) 'CALCULATED COMPRESSOR DISCHARGE
     * AT ', INT(Pc/1000+0.5), ' kPa: ', INT(Qtot*1000+0.5), ' L/S'
          IF (Pstr.EQ.'A') WRITE (1, 380) 'CALCULATED COMPRESSOR DISCHARGE
     * AT ATMOSPHERIC PRESSURE:', INT(Qtot*Denair/Denair0*1000+0.5),
     * ' L/S'
          IF (Pstr.EQ.'P') WRITE (1, 380) 'RATED minus CALCULATED DISCHA
     +RGE (at RATED Pressure) ', INT((Qcin-Qtot)*1000+0.5), ' L/S'
          IF (Pstr.EQ.'A') WRITE (1, 380) 'RATED minus CALCULATED DISCHA
     +RGE (at ATM Pressure) ', INT((Qcin-Qtot*Denair/Denair0)*1000+0.5),
     +         ' L/S'
          write(1,*)
          write(1,*)'Total Discharge from nozzles at depth ',
     +     Qsum*1000,' L/S'
          write(1,*)'Average nozzle discharge ',Qsum/Noorf*1000,
     +     ' CFM'
          write(1,*)'Coefficient of uniformity ',
     +        (Qhyd(1)-Qhyd(Noorf))/(Qsum/Noorf)
        END IF
350     FORMAT (A26, 15X, F6.1, A4)
360     FORMAT (A22, 19X, F6.1, A4)
370     FORMAT (A36, I5, A6, I5, A4)
380     FORMAT (1X, A, I5, A4)
        GOTO 20
      ELSEIF (Choice.EQ.17) THEN
        CLOSE (1)
        Diff(1,1)=0
        Supp(1,1)=0
        GOTO 10
      END IF
      CLOSE (1)
      STOP
      END

      SUBROUTINE Menu(Typesstr,Matestr,Typ,Lockstr,Km,Kr,Run,Unit,
     * Opt)
C  SUBROUTINE MENU
C
C  ROUTINE IS USED TO INPUT VARIABLES THAT ARE CONSTANT WHILE RUNNING
C  THE PROGRAM.
C
      real Kr
      integer Opt, Run, Typ, Unit
      character Typesstr(4)*30, Lockstr*80, Matestr(6)*50
      IF (Run.NE.0) THEN
        PRINT *, 'OPTIONS'
        PRINT *, '1. CHANGE SYSTEM CONFIGURATION'
        PRINT *, '2. CHANGE PIPE MATERIAL'
        PRINT *, '3. CHANGE BOTH'
390     PRINT *, 'SELECT ONE OF THE ABOVE OPTIONS (1,2,or 3)'
        READ *, Opt
        IF ((Opt.LT.1).OR.(Opt.GT.3)) GOTO 390
        IF (Opt.EQ.1) RETURN
      ELSE
400     PRINT *, 'SELECT UNITS FOR INPUT (1) ENGLISH  (2) METRIC'
        READ *, Unit
        IF ((Unit.NE.1).AND.(Unit.NE.2)) GOTO 400
        PRINT *, 'ENTER NAME OF PROJECT'
        READ '(A80)', Lockstr
C  DEFINE TYPE OF BUBBLER SYSTEM
        PRINT *, 'TYPES OF BUBBLER SYSTEMS'
        DO 410, I=1, 4
          write(*,450) I,'= ', Typesstr(I)
410     CONTINUE
420     PRINT *, 'SELECT TYPE OF BUBBLER SYSTEM FROM ABOVE'
        READ *, Typ
        IF ((Typ.LT.1).OR.(Typ.GT.4)) GOTO 420
      END IF
430   PRINT *, ' PIPE MATERIAL '
      PRINT 460, '0. "HYDRAULICALLY SMOOTH" PIPE'
      DO 440, I=1, 5
        PRINT 450, I,'. ',  Matestr(I)
440   CONTINUE
450   FORMAT (4X, I1, 2A, 50A)
460   FORMAT (4X, 15A)
      PRINT 460, '6. OTHER MATERIAL'
      PRINT *, 'PIPE MATERIAL? (0 TO 6)'
      READ *, Km
      IF (Km.EQ.0) Kr=0
      IF ((Km.EQ.1).OR.(Km.EQ.3)) Kr=0.000025
      IF ((Km.EQ.2).OR.(Km.EQ.4)) Kr=0.00025
      IF (Km.EQ.5) Kr=0.00000025
      IF (Km.EQ.6) THEN
        PRINT *, 'WHAT IS THE PIPE MATERIAL?'
        READ '(A50)', Matestr(6)
        IF (Unit.EQ.1) THEN
          PRINT *, 'ENTER PIPE ROUGHNESS K (in)'
          READ *, Kr
          Kr=Kr*0.0254
        ELSE
          PRINT *, 'ENTER PIPE ROUGHNESS K (mm)'
          READ *, Kr
          Kr=Kr/1000
        END IF
      END IF
      IF ((Km.LT.0).OR.(Km.GT.6)) THEN
        PRINT *, 'WHAT KIND OF PIPE ??, LETS TRY AGAIN'
        GOTO 430
      END IF
      RETURN
      END

      SUBROUTINE Lendiff(Unit, Imax, Noorf, Diff)
C  SUBROUTINE LENGTH OF DIFFUSER
C
C  ROUTINE IS USED FOR ENTERING TOTAL LENGTH OF DIFFUSER EITHER IN
C  ENGLISH OR S-I UNITS.  TOTAL LENGTH IS STORE IN S-I UNITS AS THE
C  LAST ENTRY IN THE DIFFUSER ARRAY.  THIS VARIABLE IS NOT ALTERED
C  DURING THE PROGRAM.
C
      implicit none
      real Diff(200,3)
      integer Unit, Imax, Noorf

      real L
      integer I
      IF (Unit.EQ.1) THEN
        PRINT *, 'ENTER LENGTH OF DIFFUSER PIPE (ft)'
        READ *, L
        Diff(Imax,1)=0.3048*L
      ELSE
        PRINT *, 'ENTER LENGTH OF DIFFUSER PIPE (M)'
        READ *, Diff(Imax,1)
      END IF
C  SPACE OUT ORIFICES ,ONE EACH END REST EQUALLY SPACED IN THE MIDDLE
      IF (Noorf.EQ.0) RETURN
      L=Diff(Imax,1)/(Noorf-1)
      DO 470, I=1, Noorf-1
        Diff(I,1)=L
 470  CONTINUE
      RETURN
      END

      SUBROUTINE Diadiff(Unit, Imax, Diff)
C
C  SUBROUTINE DIFFUSER DIAMETER
C
C  ROUTINE IS USED TO ENTER THE DIAMETER OF THE DIFFUSER PIPE IN EITHER
C  ENGLISH OR  S-I UNITS.  THE INPUT IS CONVERTED TO S-I FOR
C  STORAGE.  THE ROUTINE ALSO OBTAINS THE AREA OF THE DIFFUSER PIPE
C  IN S-I UNITS.
C
      real Diff(200,3)
      integer Unit
      IF (Unit.EQ.1) THEN
        PRINT *, 'ENTER DIAMETER OF DIFFUSER PIPE (in)'
        READ *, Dia
        Dia=0.0254*Dia
      ELSE
        PRINT *, 'ENTER DIAMETER OF DIFFUSER PIPE (cm)'
        READ *, Dia
        Dia=Dia/100
      END IF
C  OBTAIN THE AREA
      AREA=FNArea(Dia)
C  STORE THE DIA AND AREA IN THE DIFFUSER ARRAY
      DO 480, I=1, Imax
        Diff(I,2)=Dia
        Diff(I,3)=AREA
480   CONTINUE
      RETURN
      END

      SUBROUTINE Diffpipe(Unit, Imax, Noorf, Diff)
C  SUBROUTINE DIFFUSER PIPE
C
C  ROUTINE IS USED FOR ENTERING DATA FOR A DIFFUSER HAVING
C  UNEQUAL ORIFICE SPACING AND/OR DIFFERENT PIPE DIAMETERS.
C  AGAIN THE DATA IS STORED IN S-I UNITS.
C
      real Diff(200,3), L
      integer Unit
C  PRINT CURRENT VALUES IF THEY EXIST
      IF (Diff(1,1).EQ.0) THEN
        PRINT *, 'ENTER NUMBER OF ORIFICES'
        READ *, Noorf
C  DATA ENTRY
        PRINT *, 'NOTE'
        PRINT *
        PRINT *, 'INPUT DATA STARTING FROM END OF SUPPLY LINE!!'
        PRINT *, '    SL ^-----^-----^---.........--^-----^-----^| (PLUG
     *)'
        PRINT *, 'NOZZLE 1     2     3                          N'
        PRINT *
        DO 490, I=1, Noorf-1
          PRINT '(A54, I2, A5, I2)', 'ENTER LENGTH AND DIA OF PIPE SECTI
     *ON BETWEEN NOZZLES ', I,' AND ', I+1
          IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
            PRINT *, 'ENTER LEN(ft), DIA(in)'
            READ *, L, D
            L=0.3048*L
            D=0.0254*D
          ELSE
C  S-I INPUT
            PRINT *, 'ENTER LEN(M), DIA(cm)'
            READ *, L, D
            D=D/100
          END IF
          Diff(I,1)=L
          Diff(I,2)=D
          Diff(I,3)=FNArea(D)
490     CONTINUE
      END IF
C  EDITING DATA
C  PRINT CURRENT VALUES
      PRINT *, 'DIFFUSER LINE WITH UNEQUAL ORIFICE SPACING OR PIPE DIA'
      PRINT '(A15, 5X, A22, 8X, A12)', 'PIPE SECTION #', 
     *'ORIFICE SPACING ft (M)', 'DIA  in (cm)'

      DO 500, I=1, Noorf-1
        PRINT 510, I, Diff(I,1)/0.3048, ' (', Diff(I,1), ')',
     * Diff(I,2)/0.0254, ' (', Diff(I,2)*100, ')'
500   CONTINUE
510   FORMAT (7X, I2, 11X, F7.2, A2, F9.3, A1, 5X, F9.3, A2, F9.3, A1)
520   PRINT *, 'ENTER PIPE SECTION # (Neg No. to quit)'
      READ *, K
      IF (K.GT.0) THEN
        IF (Unit.EQ.1) THEN
          PRINT '(A21, F5.2, 5X, F5.3)', 'OLD VALUES len, dia:', Diff(K,
     *1)/0.3048, Diff(K,2)/0.0254
          PRINT *, 'ENTER NEW VALUES len(ft), dia(in)'
          READ *, Diff(K,1), Diff(K,2)
          PRINT '(A44, I2, 3X, F5.2, 3X, F5.3)', 'NEW VALUES OF SECTION 
     *K, len(ft), dia(in):', K, Diff(K,1), Diff(K,2)
          Diff(K,1)=Diff(K,1)*0.3048
          Diff(K,2)=Diff(K,2)*0.0254
        ELSE
          PRINT '(A21, F5.3, 5X, F5.3)', 'OLD VALUES len, dia:', Diff(K,
     *1), Diff(K,2)*100
          PRINT *, 'ENTER NEW VALUES len(M), dia(cm)'
          READ *, Diff(K,1), Diff(K,2)
          PRINT '(A42, I2, 3X, F5.3, 3X, F5.3)', 'NEW VALUES OF SECTION 
     *K, len(M), dia(cm):', K, Diff(K,1), Diff(K,2)
          Diff(K,2)=Diff(K,2)/100
        END IF
      GOTO 520
      END IF
C  DETERMINE SIZE AND TOTAL LENGTH
      Diff(Imax,1)=0
      DO 530, I=1, Noorf-1
        Diff(Imax,1)=Diff(Imax,1)+Diff(I,1)
        Diff(I,3)=FNArea(Diff(I,2))
530   CONTINUE
      RETURN
      END

      SUBROUTINE Orfdia(Unit, Dorf, Aorf)
C  SUBROUTINE ORIFICE DIAMETER
C
C  ROUTINE IS USED FOR ENTERING ORIFICE DIAMETER. CONVERTS TO S-I.
C  THE AREA OF THE PORT IS ALSO OBTAINED IN S-I UNITS.
C
      integer Unit
      IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
        PRINT *, 'ENTER DIAMETER OF ORIFICES (in)'
        READ *, Dorf
        Dorf=0.0254*Dorf
      ELSE
C  S-I INPUT
        PRINT *, 'ENTER DIAMETER OF ORIFICES (cm)'
        READ *, Dorf
        Dorf=Dorf/100
      END IF
C  AREA CALCULATION
      Aorf=FNArea(Dorf)
      RETURN
      END

      SUBROUTINE Orfspc(Unit, Imax, Noorf, Diff)
C  SUBROUTINE ORIFICE SPACING
C
C  GIVEN EITHER THE # OF ORIFICES OR AN ESTIMATED ORIFICE SPACING
C  THE ROUTINE CALCULATES THE OTHER VARIABLE.  IN THE CASE OF A GIVEN
C  SPACING THE NUMBER OF ORIFICES REQUIRED ARE ROUNDED OFF AND THEN
C  EVENLY DISTRIBUTED ALONG THE LENGTH OF THE DIFFUSER.
C
      real Diff(200,3), L
      integer Cmd, Unit
540   PRINT *, 'SELECT EITHER (1) ORIFICE SPACING OR (2) NUMBER OF ORIFI
     *CES'
      READ *, Cmd
      IF ((Cmd.NE.1).AND.(Cmd.NE.2)) GOTO 540
      IF (Cmd.EQ.1) THEN
        IF (Unit.EQ.1) THEN
          PRINT *, 'ENTER SPACING OF ORIFICES (ft)'
          READ *, L
          L=0.3048*L
        ELSE
          PRINT *, 'ENTER SPACING OF ORIFICES (M)'
          READ *, L
        END IF
C  ADJUST ORIFICE SPACING TO FIT EQUALLY INTO DIFFUSER LENGTH
        Noorf=INT(Diff(Imax,1)/L+0.5)+1
        L=Diff(Imax,1)/(Noorf-1)
      ELSE
        PRINT *, 'ENTER NUMBER OF ORIFICES'
        READ *, Noorf
        L=Diff(Imax,1)/(Noorf-1)
      END IF
      write(*,*),'Noorf ',Noorf
C  SAVE DATA
      DO 550, I=1, Noorf-1
        Diff(I,1)=L
550   CONTINUE
      RETURN
      END

      SUBROUTINE Lensupp(Unit, Imax, Supp)
C  SUBROUTINE LENGTH OF SUPPLY LINE
C
C  ROUTINE IS USED FOR ENTERING TOTAL LENGTH OF SUPPLY EITHER IN
C  ENGLISH OR S-I UNITS.  TOTAL LENGTH IS STORE IN S-I UNITS AS THE
C  LAST ENTRY IN THE SUPPLY ARRAY.  THIS VARIABLE IS NOT ALTERED
C  DURING THE PROGRAM.
C
      real Supp(200,3)
      integer Unit
      IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
        PRINT *, 'ENTER LENGTH OF SUPPLY PIPE (ft)'
        READ *, Supp(1,1)
        Supp(1,1)=0.3048*Supp(1,1)
        Supp(Imax,1)=Supp(1,1)
      ELSE
C  S-I INPUT
        PRINT *, 'ENTER LENGTH OF SUPPLY PIPE (M)'
        READ *, Supp(1,1)
        Supp(Imax,1)=Supp(1,1)
      END IF
      RETURN
      END

      SUBROUTINE Diasupp(Unit, Imax, Supp)
C  SUBROUTINE SUPPLY DIAMETER
C
C  ROUTINE IS USED TO ENTER THE DIAMETER OF THE SUPPLY PIPE IN
C  EITHER ENGLISH OR  S-I UNITS.  THE INPUT IS CONVERTED TO S-I FOR
C  STORAGE.  THE ROUTINE ALSO OBTAINS THE AREA OF THE SUPPLY PIPE
C  IN S-I UNITS.
C
      real Supp(200,3)
      integer Unit
      IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
        PRINT *, 'ENTER DIAMETER OF SUPPLY PIPE (in)'
        READ *, Dia
        Dia=0.0254*Dia
      ELSE
C  S-I INPUT
        PRINT *, 'ENTER DIAMETER OF SUPPLY PIPE (cm)'
        READ *, Dia
        Dia=Dia/100
      END IF
C  OBTAIN THE AREA
      Area=FNArea(Dia)
C  STORE THE DIA AND AREA IN THE SUPPLY ARRAY
      DO 560, I=1, Imax
        Supp(I,2)=Dia
        Supp(I,3)=Area
560   CONTINUE
      RETURN
      END

      SUBROUTINE Supppipe(Unit, Imax, Nodia, Supp)
C  SUBROUTINE SUPPLY PIPELINE
C
C  ROUTINE IS USED FOR ENTERING DATA FOR A SUPPLY LINE
C  HAVING DIFFERENT SIZE PIPE DIAMETERS. AGAIN THE DATA IS
C  STORED IN S-I UNITS.
C
      REAL Supp(200,3), L
      integer Unit
C  PRINT CURRENT VALUES IF THEY EXIST
      IF (Supp(1,1).EQ.0) THEN
C  DATA ENTRY
        DO 570, I=1, Nodia
          PRINT *, 'ENTER LENGTH AND DIA IN PIPE SECTION ', I
          IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
            PRINT *, 'ENTER LEN (ft), DIA (in)'
            READ *, L, D
            L=0.3048*L
            D=0.0254*D
          ELSE
C  S-I INPUT
            PRINT *, 'ENTER LEN (M), DIA (cm)'
            READ *, L, D
            D=D/100
          END IF
          Supp(I,1)=L
          Supp(I,2)=D
          Supp(I,3)=FNArea(D)
570     CONTINUE
      END IF
C  EDITING DATA
C  PRINT CURRENT VALUES
      PRINT *, 'SUPPLY MANIFOLD PIPING'
      PRINT 580, 'PIPE # ', 'PIPE LENGTH ft (M)', 'DIA  in (cm)'
580   FORMAT (A8, 5X, A18, 10X, A12)
      DO 590, I=1, Nodia
        PRINT 600, I, Supp(I,1)/0.3048, ' (', Supp(I,1), ')',
     * Supp(I,2)/0.0254, ' (', Supp(I,2)*100, ')'
590   CONTINUE
600   FORMAT (2X, I2, 8X, F7.2, A2, F9.3, A1, 5X, F8.3, A2, F8.3, A1)
610   PRINT *, 'ENTER PIPE # (Neg no. to quit)'
      READ *, K
      IF (K.GT.0) THEN
        IF (Unit.EQ.1) THEN
          PRINT '(A21, 2X, F5.2, 5X, F5.3)', 'OLD VALUES len, dia:', Sup
     *p(K,1)/0.3048, Supp(K,2)/0.0254
          PRINT *, 'ENTER NEW VALUES len(ft), dia(in)'
          READ *, Supp(K,1), Supp(K,2)
          PRINT '(A42, 2X, I3, 5X, F5.2, 5X, F5.3)', 'NEW VALUES OF SECT
     *ION K, len(ft), dia(in)', K, Supp(K,1), Supp(K,2)
          Supp(K,1)=Supp(K,1)*0.3048
          Supp(K,2)=Supp(K,2)*0.0254
        ELSE
          PRINT '(A21, 2X, F5.3, 5X, F5.3)', 'OLD VALUES len, dia:', Sup
     *p(K,1), Supp(K,2)*100
          PRINT *, 'ENTER NEW VALUES len(M), dia(cm)'
          READ *, Supp(K,1), Supp(K,2)
          PRINT '(A41, 2X, I3, 5X, F5.3, 5X, F5.3)', 'NEW VALUES OF SECT
     *ION K, len(M), dia(cm)', K, Supp(K,1), Supp(K,2)
          Supp(K,2)=Supp(K,2)/100
        END IF
        GOTO 610
      END IF
C  DETERMINE SIZE AND TOTAL LENGTH
      Supp(Imax,1)=0
      DO 620, I=1, Nodia
        Supp(Imax,1)=Supp(Imax,1)+Supp(I,1)
        Supp(I,3)=FNArea(Supp(I,2))
620   CONTINUE
      RETURN
      END
      SUBROUTINE Comp(Unit, Pc, Qc, Pstr)
C  SUBPROGRAM COMPRESSOR
C
C  ROUTINE IS USED TO ENTER THE NOMINAL COMPRESSOR PRESSURE
C  IN EITHER ENGLISH OR S-I UNITS.  AGAIN THE ENGLISH INPUT IS
C  CONVERTED TO S-I FOR USE BY THE PROGRAM.
C
      integer Unit
      character Pstr
      IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
        PRINT *, 'ENTER RATED COMPRESSOR PRESSURE AND DISCHARGE (psi, 
     *CFM)'
        READ *, Pc, Qc
        PRINT 650, 'ENTER a P if the Rated Discharge is at Pressure of',
     * Pc, ' psi or'
        PRINT *, 'ENTER A if the Rated Discharge is at ATMOSPHERIC Press
     *ure.'
630     PRINT *, 'ENTER CHARACTER'
        READ '(A)', Pstr
        IF (Pstr.EQ.'p') Pstr='P'
        IF (Pstr.EQ.'a') Pstr='A'
        IF ((Pstr.NE.'P').AND.(Pstr.NE.'A')) GOTO 630
        Pc=6894.757*Pc
        Qc=Qc/2118.9
      ELSE
C  S-I INPUT
        PRINT *, 'ENTER INITIAL COMPRESSOR PESSURE AND DISCHARGE (kPa, L
     */S)'
        READ *, Pc, Qc
        PRINT 650, 'ENTER a P if the Rated Discharge is at Pressure of',
     * Pc, ' kPa or'
        PRINT *, 'ENTER A if the Rated Discharge is at ATMOSPHERIC 
     *Pressure.'
640     PRINT *, 'ENTER CHARACTER '
        READ '(A)', Pstr
        IF (Pstr.EQ.'p') Pstr='P'
        IF (Pstr.EQ.'a') Pstr='A'
        IF ((Pstr.NE.'P').AND.(Pstr.NE.'A')) GOTO 640
        Pc=Pc*1000
        Qc=Qc/1000
      END IF
650   FORMAT (A51, F6.1, A7)
      RETURN
      END
      SUBROUTINE H20dpth(Unit, Hwtr)
C  SUBPROGRAM WATER DEPTH
C
C  ROUTINE IS USED TO ENTER THE DEPTH OF SUBMERGENCE IN
C  EITHER ENGLISH OR S-I UNITS. THE INPUT IS STORED IN S-I
C  FOR USE BY THE PROGRAM.
C
      integer Unit
      IF (Unit.EQ.1) THEN
C  ENGLISH INPUT
        PRINT *, 'ENTER DEPTH OF SUBMERGENCE  (ft)'
        READ *, Hwtr
        Hwtr=0.3048*Hwtr
      ELSE
C  S-I INPUT
        PRINT *, 'ENTER DEPTH OF SUBMERGENCE  (M)'
        READ *, Hwtr
      END IF
      RETURN
      END

      real FUNCTION FNArea(Dia)
C  FUNCTION AREA
C
C  FUNCTION DETERMINES THE CROSS SECTIONAL AREA OF THE PIPE
C  GIVEN THE DIAMETER.
C
      FNArea=3.14159*Dia**2/4
      RETURN
      END
      
      subroutine Diffcon(difftyp)
      integer difftyp
      PRINT *, ' CONFIGURATION OF DIFFUSER PIPE'
      PRINT *, '  1. EQUALLY SPACED ORIFICES AND CONSTANT PIPE SIZE'
      PRINT *, '  2. UNEQUALLY SPACED ORIFICES OR DIFFERENT PIPE SECTION 
     +DIAMETERS'
660   PRINT *, 'SELECT TYPE OF DIFFUSER FROM ABOVE'
      READ *, difftyp
      IF ((difftyp.LT.1).OR.(difftyp.GT.2)) GOTO 660
      RETURN
      END
      
      SUBROUTINE Expansion(Exp)
      integer Exp
      PRINT *, ' HOW SHOULD I TREAT THE AIR EXPANSION ?'
      PRINT *, '      0= AIR IS TO BE TREATED AS INCOMPRESSIBLE THROUGHO
     *UT'
      PRINT *, '      1= AIR EXPANSION AT NOZZLE IS ADIABATIC'
      PRINT *, '      2= AIR EXPANSION AT NOZZLE IS ISOTHERMAL'
670   PRINT *, 'WHICH ASSUMPTION SHOULD I USE ?'
      READ *, Exp
      IF ((Exp.LT.0).OR.(Exp.GT.2)) GOTO 670
      RETURN
      END
      
      SUBROUTINE Fric(Rey, Dia, K, F)
C     CALCULATION OF FRICTION FACTOR
      implicit none
      REAL K,Rey, Dia, F, F1, Rok
      Integer count
      count = 1
      !write(*,*)'Rey ',Rey
      F=0.3164/Rey**0.25
      !write(*,*)'F',F
680   IF (K.EQ.0) THEN
        K=0
        F1=2*LOG10(Rey*SQRT(F)-0.8)
      ELSE
        Rok=Dia/K/2
        F1=1.74+2*LOG10(Rok)-2*LOG10(1+18.7*Rok/Rey/SQRT(F))
      END IF
      if (F1.lt.1e-20) then
         write(*,*)' F1 ',F1
         stop
      endif
      F1=1/F1**2
      IF (ABS(1-F1/F).LT.0.05) RETURN
      if (count.gt.50)then
         write(*,*)
         write(*,*)'***************************************************'
         write(*,*)'WARNING: Friction factor calculation not converging'
         write(*,*)' Increase compressor pressure.'
         write(*,*)'***************************************************'
         write(*,*)
         stop
      endif
      count = count + 1
      F=F1
      GOTO 680
      END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ENDPT
c
c  Determines the end point of a string that has trailing blanks. If
c the entire string is blanks then it returns zero. This is taken from 
c Fortran With Engineering Applications, Koffman and Friedman (1993).
 
      integer function endpt(string)

      implicit none

c Passed variables
      character string*(*)

c Local variables
      character*1 blank
      parameter (blank=' ')
      integer i

c Start at the las charcater and find the first nonblank character
      do i=len(string),1,-1
         if(string(i:i).ne.blank)then
            endpt=i
            return
         endif
      enddo

c If all characters are blank
      endpt=0

      return
      end
