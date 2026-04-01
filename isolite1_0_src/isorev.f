C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF 
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP1R (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE COMMON BLOCK VARIABLES *********************************
C
      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      IF (RH.GE.DRNH42S4) THEN         ! WET AEROSOL, NEED NH4 AT SRATIO=2.0
         SULRATW = GETASR(WAER(2), RHI)     ! AEROSOL SULFATE RATIO
      ELSE
         SULRATW = 2.0D0                    ! DRY AEROSOL SULFATE RATIO
      ENDIF
      SULRAT  = WAER(3)/WAER(2)         ! SULFATE RATIO
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (SULRATW.LE.SULRAT) THEN
C

         SCASE = 'S2'
         CALL CALCS2                 ! Only liquid (metastable)

C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN
      W(2) = WAER(2)
      W(3) = WAER(3)
C
         SCASE = 'B4'
         CALL CALCB4                 ! Only liquid (metastable)
         SCASE = 'B4'
C
      CALL CALCNH3P          ! Compute NH3(g)
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
      W(2) = WAER(2)
      W(3) = WAER(3)
C
         SCASE = 'C2'
         CALL CALCC2                 ! Only liquid (metastable)
         SCASE = 'C2'

C 
      CALL CALCNH3P
C
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE ISRP1R *****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP2R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF 
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP2R (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      TRYLIQ = .TRUE.             ! Assume liquid phase, sulfate poor limit 
C
10    CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN ! *** WET AEROSOL
         SULRATW = GETASR(WAER(2), RHI)     ! LIMITING SULFATE RATIO
      ELSE
         SULRATW = 2.0D0                    ! *** DRY AEROSOL
      ENDIF
      SULRAT = WAER(3)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (SULRATW.LE.SULRAT) THEN                
C
         SCASE = 'N3'
         CALL CALCN3                 ! Only liquid (metastable)

C *** SULFATE RICH (NO ACID)
C
C     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
C     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE 
C     AEROSOL EQUILIBRIUM.
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)
C
         SCASE = 'B4'
         CALL CALCB4                 ! Only liquid (metastable)
         SCASE = 'B4'
C
C *** Add the NO3 to the solution now and calculate partitioning.
C
      MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
      MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
      CALL CALCNAP            ! HNO3, NH3 dissolved
      CALL CALCNH3P
C
C *** SULFATE RICH (FREE ACID)
C
C     FOR SOLVING THIS CASE, NITRIC ACID AND AMMONIA IN THE GAS PHASE ARE
C     ASSUMED A MINOR SPECIES, THAT DO NOT SIGNIFICANTLY AFFECT THE 
C     AEROSOL EQUILIBRIUM.
C
      ELSEIF (SULRAT.LT.1.0) THEN             
      W(2) = WAER(2)
      W(3) = WAER(3)
      W(4) = WAER(4)
C
         SCASE = 'C2'
         CALL CALCC2                 ! Only liquid (metastable)
         SCASE = 'C2'

C
C *** Add the NO3 to the solution now and calculate partitioning.
C
      MOLAL(7) = WAER(4)             ! There is always water, so NO3(aer) is NO3-
      MOLAL(1) = MOLAL(1) + WAER(4)  ! Add H+ to balance out
C
      CALL CALCNAP                   ! HNO3, NH3 dissolved
      CALL CALCNH3P
      ENDIF
C
C *** IF SULRATW < SULRAT < 2.0 and WATER = 0 => SULFATE RICH CASE.
C
      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  
     &                                    .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE ISRP2R *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP3R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP3R (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ
ccC
ccC *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
ccC
cc      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
cc      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** INITIALIZE ALL VARIABLES ******************************************
C
      TRYLIQ = .TRUE.             ! Use liquid phase sulfate poor limit 
C
10    CALL ISOINIT3 (WI, RHI, TEMPI) ! COMMON block variables
ccC
ccC *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
ccC
cc      REST = 2.D0*WAER(2) + WAER(4) + WAER(5) 
cc      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
cc         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
cc         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
cc      ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
      IF (TRYLIQ .AND. RH.GE.DRNH4NO3) THEN  ! ** WET AEROSOL
         FRSO4   = WAER(2) - WAER(1)/2.0D0     ! SULFATE UNBOUND BY SODIUM
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR(FRSO4, RHI)          ! SULFATE RATIO FOR NH4+
         SULRATW = (WAER(1)+FRSO4*SRI)/WAER(2) ! LIMITING SULFATE RATIO
         SULRATW = MIN (SULRATW, 2.0D0)
      ELSE
         SULRATW = 2.0D0                     ! ** DRY AEROSOL
      ENDIF
      SULRAT = (WAER(1)+WAER(3))/WAER(2)
      SODRAT = WAER(1)/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM POOR
C
      IF (SULRATW.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
C
         SCASE = 'Q5'
         CALL CALCQ5                 ! Only liquid (metastable)
         SCASE = 'Q5'
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.SULRATW .AND. SODRAT.GE.2.0) THEN                
C
         SCASE = 'R6'
         CALL CALCR6                 ! Only liquid (metastable)
         SCASE = 'R6'
C
C *** SULFATE RICH (NO ACID) 
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.SULRATW) THEN 
      DO 100 I=1,NCOMP
         W(I) = WAER(I)
100   CONTINUE
C
         SCASE = 'I6'
         CALL CALCI6                 ! Only liquid (metastable)
         SCASE = 'I6'
C
      CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
      CALL CALCNH3P
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
      DO 200 I=1,NCOMP
         W(I) = WAER(I)
200   CONTINUE
C
         SCASE = 'J3'
         CALL CALCJ3                 ! Only liquid (metastable)
         SCASE = 'J3'
C
      CALL CALCNHP                ! HNO3, NH3, HCL in gas phase
      CALL CALCNH3P
C
      ENDIF
C
C *** IF AFTER CALCULATIONS, SULRATW < SULRAT < 2.0  
C                            and WATER = 0          => SULFATE RICH CASE.
C
      IF (SULRATW.LE.SULRAT .AND. SULRAT.LT.2.0  
     &                      .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE ISRP3R *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE ISRP4R
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE REVERSE PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTTASIUM-MAGNESIUM AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP4R (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
      LOGICAL   TRYLIQ
ccC
ccC *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
ccC
cc      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
cc      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** INITIALIZE ALL VARIABLES ******************************************
C
      TRYLIQ  = .TRUE.             ! Use liquid phase sulfate poor limit
      IPROB   = 1            ! SOLVE REVERSE PROBLEM
C
10    CALL INIT4 (WI, RHI, TEMPI) ! COMMON block variables
ccC
ccC *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
ccC
cc      REST = 2.D0*WAER(2) + WAER(4) + WAER(5)
cc      IF (WAER(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
cc         WAER(1) = (ONE-1D-6)*REST         ! Adjust Na amount
cc         CALL PUSHERR (0050, 'ISRP3R')     ! Warning error: Na adjusted
cc      ENDIF
C
C *** CALCULATE SULFATE, CRUSTAL & SODIUM RATIOS ***********************
C
      IF (TRYLIQ) THEN                               ! ** WET AEROSOL
         FRSO4   = WAER(2) - WAER(1)/2.0D0
     &           - WAER(6) - WAER(7)/2.0D0 - WAER(8) ! SULFATE UNBOUND BY SODIUM,CALCIUM,POTTASIUM,MAGNESIUM
         FRSO4   = MAX(FRSO4, TINY)
         SRI     = GETASR(FRSO4, RHI)                ! SULFATE RATIO FOR NH4+
         SULRATW = (WAER(1)+FRSO4*SRI+WAER(6)
     &              +WAER(7)+WAER(8))/WAER(2)       ! LIMITING SULFATE RATIO
         SULRATW = MIN (SULRATW, 2.0D0)
      ELSE
         SULRATW = 2.0D0                     ! ** DRY AEROSOL
      ENDIF
      SO4RAT = (WAER(1)+WAER(3)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
      CRNARAT = (WAER(1)+WAER(6)+WAER(7)+WAER(8))/WAER(2)
      CRRAT  = (WAER(6)+WAER(7)+WAER(8))/WAER(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR ; SODIUM+CRUSTALS POOR
C
      IF (SULRATW.LE.SO4RAT .AND. CRNARAT.LT.2.0) THEN
C
         SCASE = 'V7'
         CALL CALCV7                 ! Only liquid (metastable)
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
      ELSEIF (SO4RAT.GE.SULRATW .AND. CRNARAT.GE.2.0) THEN
C
       IF (CRRAT.LE.2.0) THEN
C
         SCASE = 'U8'
         CALL CALCU8                 ! Only liquid (metastable)
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
       ELSEIF (CRRAT.GT.2.0) THEN
C
         SCASE = 'W13'
         CALL CALCW13                 ! Only liquid (metastable)

C        CALL CALCNH3
       ENDIF
C
C *** SULFATE RICH (NO ACID): 1<Rso4<2;
C
      ELSEIF (1.0.LE.SO4RAT .AND. SO4RAT.LT.SULRATW) THEN
      DO 800 I=1,NCOMP
         W(I) = WAER(I)
 800  CONTINUE
C
         SCASE = 'L9'
         CALL CALCL9                 ! Only liquid (metastable)

C
      CALL CALCNHP                ! MINOR SPECIES: HNO3, HCl
      CALL CALCNH3P               !                NH3
C
C *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
C
      ELSEIF (SO4RAT.LT.1.0) THEN
      DO 900 I=1,NCOMP
         W(I) = WAER(I)
 900  CONTINUE
C
         SCASE = 'K4'
         CALL CALCK4                 ! Only liquid (metastable)

C
      CALL CALCNHP                  ! MINOR SPECIES: HNO3, HCl
      CALL CALCNH3P                 !                NH3
C
      ENDIF
C
C *** IF AFTER CALCULATIONS, SULRATW < SO4RAT < 2.0
C                            and WATER = 0          => SULFATE RICH CASE.
C
      IF (SULRATW.LE.SO4RAT .AND. SO4RAT.LT.2.0
     &                      .AND. WATER.LE.TINY) THEN
          TRYLIQ = .FALSE.
          GOTO 10
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE ISRP4R *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCS2
C *** CASE S2
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCS2
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION NH4I, NH3GI, NH3AQ
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU   =.TRUE.     ! Outer loop activity calculation flag
      FRST     =.TRUE.
      CALAIN   =.TRUE.
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLALR(4)= MIN(WAER(2), 0.5d0*WAER(3))
      WATER    = MOLALR(4)/M0(4)  ! ZSR correlation
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
CC         A21  = XK21*WATER*R*TEMP
         A2   = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
         AKW  = XKW *RH*WATER*WATER
C
         NH4I = WAER(3)
         SO4I = WAER(2)
         HSO4I= ZERO
C
         CALL CALCPH (2.D0*SO4I - NH4I, HI, OHI)    ! Get pH
C
         NH3AQ = ZERO                               ! AMMONIA EQUILIBRIUM
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ (NH4I, OHI, DEL)
            NH4I  = MAX (NH4I-DEL, ZERO) 
            OHI   = MAX (OHI -DEL, TINY)
            NH3AQ = DEL
            HI    = AKW/OHI
         ENDIF
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)         ! SULFATE EQUILIBRIUM
         SO4I  = SO4I - DEL
         HI    = HI   - DEL
         HSO4I = DEL
C
         NH3GI = NH4I/HI/A2   !    NH3AQ/A21
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL(1) = HI
         MOLAL(3) = NH4I
         MOLAL(5) = SO4I
         MOLAL(6) = HSO4I
         COH      = OHI
         GASAQ(1) = NH3AQ
         GNH3     = NH3GI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE
C
20    RETURN
C
C *** END OF SUBROUTINE CALCS2 ****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCN3
C *** CASE N3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCN3
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION NH4I, NO3I, NH3AQ, NO3AQ
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU =.TRUE.              ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** AEROSOL WATER CONTENT
C
      MOLALR(4) = MIN(WAER(2),0.5d0*WAER(3))       ! (NH4)2SO4
      AML5      = MAX(WAER(3)-2.D0*MOLALR(4),ZERO) ! "free" NH4
      MOLALR(5) = MAX(MIN(AML5,WAER(4)), ZERO)     ! NH4NO3=MIN("free",NO3)
      WATER     = MOLALR(4)/M0(4) + MOLALR(5)/M0(5)
      WATER     = MAX(WATER, TINY)
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
         A2    = XK2 *R*TEMP/XKW/RH*(GAMA(8)/GAMA(9))**2.
CC         A21   = XK21*WATER*R*TEMP
         A3    = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4    = XK7*(WATER/GAMA(4))**3.0
         AKW   = XKW *RH*WATER*WATER
C
C ION CONCENTRATIONS
C
         NH4I  = WAER(3)
         NO3I  = WAER(4)
         SO4I  = WAER(2)
         HSO4I = ZERO
C
         CALL CALCPH (2.D0*SO4I + NO3I - NH4I, HI, OHI)
C
C AMMONIA ASSOCIATION EQUILIBRIUM
C
         NH3AQ = ZERO
         NO3AQ = ZERO
         GG    = 2.D0*SO4I + NO3I - NH4I
         IF (HI.LT.OHI) THEN
            CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
            HI    = AKW/OHI
         ELSE
            HI    = ZERO
            CALL CALCNIAQ2 (GG, NO3I, HI, NO3AQ) ! HNO3
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
            CALL CALCHS4 (HI, SO4I, ZERO, DEL)
            SO4I  = SO4I  - DEL
            HI    = HI    - DEL
            HSO4I = DEL
            OHI   = AKW/HI
         ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
         MOLAL (1) = HI
         MOLAL (3) = NH4I
         MOLAL (5) = SO4I
         MOLAL (6) = HSO4I
         MOLAL (7) = NO3I
         COH       = OHI
C
         CNH42S4   = ZERO
         CNH4NO3   = ZERO
C
         GASAQ(1)  = NH3AQ
         GASAQ(3)  = NO3AQ
C
         GHNO3     = HI*NO3I/A3
         GNH3      = NH4I/HI/A2   !   NH3AQ/A21 
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP ******************
C
         IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
            CALL CALCACT     
         ELSE
            GOTO 20
         ENDIF
10    CONTINUE
C
C *** RETURN ***********************************************************
C
20    RETURN
C
C *** END OF SUBROUTINE CALCN3 *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCQ5
C *** CASE Q5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); SODIUM POOR (SODRAT < 2.0)
C     2. LIQUID AND SOLID PHASES ARE POSSIBLE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCQ5
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      FRST    =.TRUE.
      CALAIN  =.TRUE. 
      CALAOU  =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
      CALL CALCQ1A
C
      PSI1   = CNA2SO4      ! SALTS DISSOLVED
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
C
      CALL CALCMR           ! WATER
C
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
      AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
C
C ION CONCENTRATIONS
C
      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
C
C SOLUTION ACIDIC OR BASIC?
C
      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
         HSO4I = ZERO
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
ccc      CALL PUSHERR (0002, 'CALCQ5')    ! WARNING ERROR: NO CONVERGENCE
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4
C
      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ
C
      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CNA2SO4 = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCQ5 ******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCQ1A
C *** CASE Q1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCQ1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE SOLIDS **************************************************
C
      CNA2SO4 = 0.5d0*WAER(1)
      FRSO4   = MAX (WAER(2)-CNA2SO4, ZERO)
C
      CNH42S4 = MAX (MIN(FRSO4,0.5d0*WAER(3)), TINY)
      FRNH3   = MAX (WAER(3)-2.D0*CNH42S4, ZERO)
C
      CNH4NO3 = MIN (FRNH3, WAER(4))
CCC      FRNO3   = MAX (WAER(4)-CNH4NO3, ZERO)
      FRNH3   = MAX (FRNH3-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRNH3, WAER(5))
CCC      FRCL    = MAX (WAER(5)-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCQ1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCR6
C *** CASE R6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCR6
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      CALL CALCR1A
C
      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
C
      FRST   = .TRUE.
      CALAIN = .TRUE. 
      CALAOU = .TRUE. 
C
C *** CALCULATE WATER **************************************************
C
      CALL CALCMR
C
C *** SETUP LIQUID CONCENTRATIONS **************************************
C
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
      AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+      
C
      NAI    = WAER(1)
      SO4I   = WAER(2)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
C
C SOLUTION ACIDIC OR BASIC?
C
      GG  = 2.D0*WAER(2) + NO3I + CLI - NAI - NH4I
      IF (GG.GT.TINY) THEN                        ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.LT.OHI) THEN
         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
         HI    = AKW/OHI
      ELSE
         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I, ZERO)
         GGCL  = MAX(GG-GGNO3, ZERO)
         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
         IF (GGNO3.GT.TINY) THEN
            IF (GGCL.LE.TINY) HI = ZERO
            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
         SO4I  = SO4I  - DEL
         HI    = HI    - DEL
         HSO4I = DEL
         OHI   = AKW/HI
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
ccc      CALL PUSHERR (0002, 'CALCR6')    ! WARNING ERROR: NO CONVERGENCE
C 
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
20    A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3     = NH4I/HI/A2
      GHNO3    = HI*NO3I/A3
      GHCL     = HI*CLI /A4
C
      GASAQ(1) = NH3AQ
      GASAQ(2) = CLAQ
      GASAQ(3) = NO3AQ
C
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = ZERO 
C
      RETURN
C
C *** END OF SUBROUTINE CALCR6 ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCR1A
C *** CASE R1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4NO3, NH4CL, NANO3, NA2SO4, NANO3, NACL
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCR1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE SOLIDS **************************************************
C
      CNA2SO4 = WAER(2)
      FRNA    = MAX (WAER(1)-2*CNA2SO4, ZERO)
C
      CNH42S4 = ZERO
C
      CNANO3  = MIN (FRNA, WAER(4))
      FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
      FRNA    = MAX (FRNA-CNANO3, ZERO)
C
      CNACL   = MIN (FRNA, WAER(5))
      FRCL    = MAX (WAER(5)-CNACL, ZERO)
      FRNA    = MAX (FRNA-CNACL, ZERO)
C
      CNH4NO3 = MIN (FRNO3, WAER(3))
      FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
      FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRCL, FRNH3)
      FRCL    = MAX (FRCL-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCR1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCV7
C *** CASE V7
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCV7
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
      CALL CALCV1A
C
      CHI9   = CCASO4
C
      PSI1   = CNA2SO4      ! SALTS DISSOLVED
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI6   = CNH42S4
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4
C
      CALL CALCMR           ! WATER
C
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
C
C ION CONCENTRATIONS
C
      NAI    = WAER(1)
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      KI     = WAER(7)
      MGI    = WAER(8)
C
C SOLUTION ACIDIC OR BASIC?
C
      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
     &       - 2.D0*CAI - KI - 2.D0*MGI
      IF (GG.GT.TINY) THEN                    ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                     ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF

C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.GT.OHI) THEN
C         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
C         HI    = AKW/OHI
C         HSO4I = ZERO
C      ELSE
C         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
C     &           - KI - 2.D0*MGI, ZERO)
C         GGCL  = MAX(GG-GGNO3, ZERO)
C         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
C         IF (GGNO3.GT.TINY) THEN
C            IF (GGCL.LE.TINY) HI = ZERO
C            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
C         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL
C         IF (HI.LE.TINY) HI = SQRT(AKW)
      OHI   = AKW/HI
C
      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
ccc      CALL PUSHERR (0002, 'CALCV7')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4
C
      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ
C
      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNA2SO4 = ZERO
      CMGSO4  = ZERO
      CK2SO4  = ZERO
      CCASO4  = MIN (WAER(6), WAER(2))
C
      RETURN
C
C *** END OF SUBROUTINE CALCV7 ******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCV1A
C *** CASE V1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SO4RAT > 2.0), Cr+NA poor (CRNARAT < 2)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3, NH4Cl, NA2SO4, K2SO4, MGSO4, CASO4
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCV1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE SOLIDS **************************************************
C
      CCASO4  = MIN (WAER(6), WAER(2))                     ! CCASO4
      SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
      CAFR    = MAX (WAER(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)                 ! CK2SO4
      FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*WAER(1), SO4FR)                 ! CNA2SO4
      NAFR    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CNA2SO4, ZERO)
      CMGSO4  = MIN (WAER(8), SO4FR)                       ! CMGSO4
      FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
      CNH42S4 = MAX (MIN (SO4FR , 0.5d0*WAER(3)) , TINY)
      FRNH3   = MAX (WAER(3) - 2.D0*CNH42S4, ZERO)
C
      CNH4NO3 = MIN (FRNH3, WAER(4))
CCC      FRNO3   = MAX (WAER(4) - CNH4NO3, ZERO)
      FRNH3   = MAX (FRNH3 - CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRNH3, WAER(5))
CCC      FRCL    = MAX (WAER(5) - CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3 - CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCV1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCU8
C *** CASE U8
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCU8
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      CALL CALCU1A
C
      CHI9   = CCASO4        ! SALTS
C
      PSI1   = CNA2SO4
      PSI2   = CNANO3
      PSI3   = CNACL
      PSI4   = CNH4CL
      PSI5   = CNH4NO3
      PSI7   = CK2SO4
      PSI8   = CMGSO4
      PSI9   = CCASO4
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      CALAOU = .TRUE.
C
C *** CALCULATE WATER **************************************************
C
      CALL CALCMR
C
C *** SETUP LIQUID CONCENTRATIONS **************************************
C
      HSO4I  = ZERO
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      AKW = XKW*RH*WATER*WATER                        ! H2O    <==> H+
C
      NAI    = WAER(1)
      SO4I   = MAX(WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      KI     = WAER(7)
      MGI    = WAER(8)

C
C SOLUTION ACIDIC OR BASIC?
C
      GG  = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
     &       - 2.D0*CAI - KI - 2.D0*MGI
      IF (GG.GT.TINY) THEN                        ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF
      IF (HI.LE.TINY) HI = SQRT(AKW)
C      OHI   = AKW/HI
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.GT.OHI) THEN
C         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
C         HI    = AKW/OHI
C         HSO4I = ZERO
C      ELSE
C         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
C     &           - KI - 2.D0*MGI, ZERO)
C         GGCL  = MAX(GG-GGNO3, ZERO)
C         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
C         IF (GGNO3.GT.TINY) THEN
C            IF (GGCL.LE.TINY) HI = ZERO
C            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
C         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL
C         IF (HI.LE.TINY) HI = SQRT(AKW)
      OHI   = AKW/HI
C
      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
ccc      CALL PUSHERR (0002, 'CALCU8')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
20    A2       = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3       = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4       = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3     = NH4I/HI/A2
      GHNO3    = HI*NO3I/A3
      GHCL     = HI*CLI /A4
C
      GASAQ(1) = NH3AQ
      GASAQ(2) = CLAQ
      GASAQ(3) = NO3AQ
C
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4CL   = ZERO
      CNACL    = ZERO
      CNANO3   = ZERO
      CNA2SO4  = ZERO
      CMGSO4   = ZERO
      CK2SO4   = ZERO
      CCASO4   = MIN (WAER(6), WAER(2))
C
      RETURN
C
C *** END OF SUBROUTINE CALCU8 ******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCU1A
C *** CASE U1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0); CRUSTAL+SODIUM RICH (CRNARAT >= 2.0); CRUSTAL POOR (CRRAT<2)
C     2. THERE IS ONLY A SOLID PHASE
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCU1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE SOLIDS *************************************************
C
      CCASO4  = MIN (WAER(6), WAER(2))                 ! CCASO4
      SO4FR   = MAX(WAER(2) - CCASO4, ZERO)
      CAFR    = MAX(WAER(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*WAER(7), SO4FR)             ! CK2SO4
      FRK     = MAX(WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX(SO4FR - CK2SO4, ZERO)
      CMGSO4  = MIN (WAER(8), SO4FR)                   ! CMGSO4
      FRMG    = MAX(WAER(8) - CMGSO4, ZERO)
      SO4FR   = MAX(SO4FR - CMGSO4, ZERO)
      CNA2SO4 = MAX (SO4FR, ZERO)                      ! CNA2SO4
      FRNA    = MAX (WAER(1) - 2.D0*CNA2SO4, ZERO)
C
      CNH42S4 = ZERO
C
      CNANO3  = MIN (FRNA, WAER(4))
      FRNO3   = MAX (WAER(4)-CNANO3, ZERO)
      FRNA    = MAX (FRNA-CNANO3, ZERO)
C
      CNACL   = MIN (FRNA, WAER(5))
      FRCL    = MAX (WAER(5)-CNACL, ZERO)
      FRNA    = MAX (FRNA-CNACL, ZERO)
C
      CNH4NO3 = MIN (FRNO3, WAER(3))
      FRNO3   = MAX (FRNO3-CNH4NO3, ZERO)
      FRNH3   = MAX (WAER(3)-CNH4NO3, ZERO)
C
      CNH4CL  = MIN (FRCL, FRNH3)
      FRCL    = MAX (FRCL-CNH4CL, ZERO)
      FRNH3   = MAX (FRNH3-CNH4CL, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCU1A *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCW13
C *** CASE W13
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
C                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCW13
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION NH4I, NAI, NO3I, NH3AQ, NO3AQ, CLAQ, CAI, KI, MGI
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** SETUP PARAMETERS ************************************************
C
      FRST    =.TRUE.
      CALAIN  =.TRUE.
      CALAOU  =.TRUE.
C
C *** CALCULATE INITIAL SOLUTION ***************************************
C
      CALL CALCW1A
C
      CHI11   = CCASO4
C
      PSI1   = CNA2SO4      ! SALTS DISSOLVED
      PSI5   = CNH4CL
      PSI6   = CNH4NO3
      PSI7   = CNACL
      PSI8   = CNANO3
      PSI9   = CK2SO4
      PSI10  = CMGSO4
      PSI11  = CCASO4
      PSI12  = CCANO32
      PSI13  = CKNO3
      PSI14  = CKCL
      PSI15  = CMGNO32
      PSI16  = CMGCL2
      PSI17  = CCACL2
C
      CALL CALCMR           ! WATER
C
      NH3AQ  = ZERO
      NO3AQ  = ZERO
      CLAQ   = ZERO
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP

      AKW = XKW*RH*WATER*WATER               ! H2O       <==> H+
C
C ION CONCENTRATIONS
C
      NAI    = WAER(1)
      SO4I   = MAX (WAER(2) - WAER(6), ZERO)
      NH4I   = WAER(3)
      NO3I   = WAER(4)
      CLI    = WAER(5)
      CAI    = ZERO
      KI     = WAER(7)
      MGI    = WAER(8)
C
C SOLUTION ACIDIC OR BASIC?
C
      GG   = 2.D0*SO4I + NO3I + CLI - NAI - NH4I
     &       - 2.D0*CAI - KI - 2.D0*MGI
      IF (GG.GT.TINY) THEN                        ! H+ in excess
         BB =-GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         HI = 0.5D0*(-BB + SQRT(DD))
         OHI= AKW/HI
      ELSE                                        ! OH- in excess
         BB = GG
         CC =-AKW
         DD = BB*BB - 4.D0*CC
         OHI= 0.5D0*(-BB + SQRT(DD))
         HI = AKW/OHI
      ENDIF
C
C UNDISSOCIATED SPECIES EQUILIBRIA
C
      IF (HI.GT.OHI) THEN
C         CALL CALCAMAQ2 (-GG, NH4I, OHI, NH3AQ)
C         HI    = AKW/OHI
C         HSO4I = ZERO
C      ELSE
C         GGNO3 = MAX(2.D0*SO4I + NO3I - NAI - NH4I - 2.D0*CAI
C     &           - KI - 2.D0*MGI, ZERO)
C         GGCL  = MAX(GG-GGNO3, ZERO)
C         IF (GGCL .GT.TINY) CALL CALCCLAQ2 (GGCL, CLI, HI, CLAQ) ! HCl
C         IF (GGNO3.GT.TINY) THEN
C            IF (GGCL.LE.TINY) HI = ZERO
C            CALL CALCNIAQ2 (GGNO3, NO3I, HI, NO3AQ)              ! HNO3
C         ENDIF
C
C CONCENTRATION ADJUSTMENTS ; HSO4 minor species.
C
         CALL CALCHS4 (HI, SO4I, ZERO, DEL)
      else
        del= zero
      ENDIF
      SO4I  = SO4I  - DEL
      HI    = HI    - DEL
      HSO4I = DEL
C         IF (HI.LE.TINY) HI = SQRT(AKW)
      OHI   = AKW/HI
C
      IF (HI.LE.TINY) THEN
      HI = SQRT(AKW)
      OHI   = AKW/HI
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL(1) = HI
      MOLAL(2) = NAI
      MOLAL(3) = NH4I
      MOLAL(4) = CLI
      MOLAL(5) = SO4I
      MOLAL(6) = HSO4I
      MOLAL(7) = NO3I
      MOLAL(8) = CAI
      MOLAL(9) = KI
      MOLAL(10)= MGI
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT
      ELSE
         GOTO 20
      ENDIF
10    CONTINUE
ccc      CALL PUSHERR (0002, 'CALCW13')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
20    A2      = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2. ! NH3  <==> NH4+
      A3      = XK4 *R*TEMP*(WATER/GAMA(10))**2.        ! HNO3 <==> NO3-
      A4      = XK3 *R*TEMP*(WATER/GAMA(11))**2.        ! HCL  <==> CL-
C
      GNH3    = NH4I/HI/A2
      GHNO3   = HI*NO3I/A3
      GHCL    = HI*CLI /A4
C
      GASAQ(1)= NH3AQ
      GASAQ(2)= CLAQ
      GASAQ(3)= NO3AQ
C
      CNH42S4 = ZERO
      CNH4NO3 = ZERO
      CNH4CL  = ZERO
      CNACL   = ZERO
      CNANO3  = ZERO
      CMGSO4  = ZERO
      CK2SO4  = ZERO
      CCASO4  = MIN (WAER(6), WAER(2))
      CCANO32 = ZERO
      CKNO3   = ZERO
      KCL     = ZERO
      CMGNO32 = ZERO
      CMGCL2  = ZERO
      CCACL2  = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCW13 ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCW1A
C *** CASE W1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : CaSO4, CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
C                          MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCW1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE SOLIDS **************************************************
C
      CCASO4  = MIN (WAER(2), WAER(6))              !SOLID CASO4
      CAFR    = MAX (WAER(6) - CCASO4, ZERO)
      SO4FR   = MAX (WAER(2) - CCASO4, ZERO)
      CK2SO4  = MIN (SO4FR, 0.5D0*WAER(7))          !SOLID K2SO4
      FRK     = MAX (WAER(7) - 2.D0*CK2SO4, ZERO)
      SO4FR   = MAX (SO4FR - CK2SO4, ZERO)
      CMGSO4  = SO4FR                               !SOLID MGSO4
      FRMG    = MAX (WAER(8) - CMGSO4, ZERO)
      CNACL   = MIN (WAER(1), WAER(5))              !SOLID NACL
      FRNA    = MAX (WAER(1) - CNACL, ZERO)
      CLFR    = MAX (WAER(5) - CNACL, ZERO)
      CCACL2  = MIN (CAFR, 0.5D0*CLFR)              !SOLID CACL2
      CAFR    = MAX (CAFR - CCACL2, ZERO)
      CLFR    = MAX (WAER(5) - 2.D0*CCACL2, ZERO)
      CCANO32 = MIN (CAFR, 0.5D0*WAER(4))           !SOLID CA(NO3)2
      CAFR    = MAX (CAFR - CCANO32, ZERO)
      FRNO3   = MAX (WAER(4) - 2.D0*CCANO32, ZERO)
      CMGCL2  = MIN (FRMG, 0.5D0*CLFR)              !SOLID MGCL2
      FRMG    = MAX (FRMG - CMGCL2, ZERO)
      CLFR    = MAX (CLFR - 2.D0*CMGCL2, ZERO)
      CMGNO32 = MIN (FRMG, 0.5D0*FRNO3)             !SOLID MG(NO3)2
      FRMG    = MAX (FRMG - CMGNO32, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CMGNO32, ZERO)
      CNANO3  = MIN (FRNA, FRNO3)                   !SOLID NANO3
      FRNA    = MAX (FRNA - CNANO3, ZERO)
      FRNO3   = MAX (FRNO3 - CNANO3, ZERO)
      CKCL    = MIN (FRK, CLFR)                     !SOLID KCL
      FRK     = MAX (FRK - CKCL, ZERO)
      CLFR    = MAX (CLFR - CKCL, ZERO)
      CKNO3   = MIN (FRK, FRNO3)                    !SOLID KNO3
      FRK     = MAX (FRK - CKNO3, ZERO)
      FRNO3   = MAX (FRNO3 - CKNO3, ZERO)
C
C *** OTHER PHASES ******************************************************
C
      WATER   = ZERO
C
      GNH3    = ZERO
      GHNO3   = ZERO
      GHCL    = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCW1A *****************************************
C
      END
