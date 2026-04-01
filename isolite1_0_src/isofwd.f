C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP1F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY 
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP1F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT1 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (2.0.LE.SULRAT) THEN 
      DC   = W(3) - 2.001D0*W(2)  ! For numerical stability
      W(3) = W(3) + MAX(-DC, ZERO)
C
         SCASE = 'A2'
         CALL CALCA2                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID)
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
         SCASE = 'B4'
         CALL CALCB4                 ! Only liquid (metastable)
  
      CALL CALCNH3
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
         SCASE = 'C2'
         CALL CALCC2                 ! Only liquid (metastable)

      CALL CALCNH3
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP1F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP2F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FOREWARD PROBLEM OF 
C     AN AMMONIUM-SULFATE-NITRATE AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE RATIO AND BY
C     THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP2F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT2 (WI, RHI, TEMPI)
C
C *** CALCULATE SULFATE RATIO *******************************************
C
      SULRAT = W(3)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************
C
C *** SULFATE POOR 
C
      IF (2.0.LE.SULRAT) THEN                
C
         SCASE = 'D3'
         CALL CALCD3                 ! Only liquid (metastable)
   
C
C *** SULFATE RICH (NO ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM.
C     SUBROUTINES CALCB? ARE CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
         SCASE = 'B4'
         CALL CALCB4                 ! Only liquid (metastable)
         SCASE = 'E4'
C
      CALL CALCNA                 ! HNO3(g) DISSOLUTION
C
C *** SULFATE RICH (FREE ACID)
C     FOR SOLVING THIS CASE, NITRIC ACID IS ASSUMED A MINOR SPECIES, 
C     THAT DOES NOT SIGNIFICANTLY PERTURB THE HSO4-SO4 EQUILIBRIUM
C     SUBROUTINE CALCC? IS CALLED, AND THEN THE NITRIC ACID IS DISSOLVED
C     FROM THE HNO3(G) -> (H+) + (NO3-) EQUILIBRIUM.
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
         SCASE = 'C2'
         CALL CALCC2                 ! Only liquid (metastable)
         SCASE = 'F2'
C
      CALL CALCNA                 ! HNO3(g) DISSOLUTION
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP2F *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE ISRP3F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM AEROSOL SYSTEM. 
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM 
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE ISRP3F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
C
      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
      ENDIF
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL ISOINIT3 (WI, RHI, TEMPI)
C
C *** CHECK IF TOO MUCH SODIUM ; ADJUST AND ISSUE ERROR MESSAGE *********
C
      REST = 2.D0*W(2) + W(4) + W(5) 
      IF (W(1).GT.REST) THEN            ! NA > 2*SO4+CL+NO3 ?
         W(1) = (ONE-1D-6)*REST         ! Adjust Na amount
         CALL PUSHERR (0050, 'ISRP3F')  ! Warning error: Na adjusted
      ENDIF
C
C *** CALCULATE SULFATE & SODIUM RATIOS *********************************
C
      SULRAT = (W(1)+W(3))/W(2)
      SODRAT = W(1)/W(2)
C
C *** FIND CALCULATION REGIME FROM (SULRAT,RH) **************************

C *** SULFATE POOR ; SODIUM POOR
C
      IF (2.0.LE.SULRAT .AND. SODRAT.LT.2.0) THEN                
C
         SCASE = 'G5'
         CALL CALCG5                 ! Only liquid (metastable)
   
C
C *** SULFATE POOR ; SODIUM RICH
C
      ELSE IF (SULRAT.GE.2.0 .AND. SODRAT.GE.2.0) THEN                
C
         SCASE = 'H6'
         CALL CALCH6                 ! Only liquid (metastable)
C
C *** SULFATE RICH (NO ACID) 
C
      ELSEIF (1.0.LE.SULRAT .AND. SULRAT.LT.2.0) THEN 
C
         SCASE = 'I6'
         CALL CALCI6                 ! Only liquid (metastable)
C                                    
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL CALCNH3                !                NH3 
C
C *** SULFATE RICH (FREE ACID)
C
      ELSEIF (SULRAT.LT.1.0) THEN             
C
         SCASE = 'J3'
         CALL CALCJ3                 ! Only liquid (metastable)

C                                    
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl       
      CALL CALCNH3                !                NH3 
      ENDIF
C
C *** RETURN POINT
C
      RETURN
C
C *** END OF SUBROUTINE ISRP3F *****************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE ISRP4F
C *** THIS SUBROUTINE IS THE DRIVER ROUTINE FOR THE FORWARD PROBLEM OF
C     AN AMMONIUM-SULFATE-NITRATE-CHLORIDE-SODIUM-CALCIUM-POTASSIUM-MAGNESIUM
C     AEROSOL SYSTEM.
C     THE COMPOSITION REGIME IS DETERMINED BY THE SULFATE & SODIUM
C     RATIOS AND BY THE AMBIENT RELATIVE HUMIDITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ISRP4F (WI, RHI, TEMPI)
      INCLUDE 'isrpia.inc'
      DIMENSION WI(NCOMP)
      DOUBLE PRECISION NAFRI, NO3FRI
C
C *** ADJUST FOR TOO LITTLE AMMONIUM AND CHLORIDE ***********************
C
C      WI(3) = MAX (WI(3), 1.D-10)  ! NH4+ : 1e-4 umoles/m3
C      WI(5) = MAX (WI(5), 1.D-10)  ! Cl-  : 1e-4 umoles/m3
C
C *** ADJUST FOR TOO LITTLE SODIUM, SULFATE AND NITRATE COMBINED ********
C
C      IF (WI(1)+WI(2)+WI(4) .LE. 1d-10) THEN
C         WI(1) = 1.D-10  ! Na+  : 1e-4 umoles/m3
C         WI(2) = 1.D-10  ! SO4- : 1e-4 umoles/m3
C      ENDIF
C
C *** INITIALIZE ALL VARIABLES IN COMMON BLOCK **************************
C
      CALL INIT4 (WI, RHI, TEMPI)
C
C *** CHECK IF TOO MUCH SODIUM+CRUSTALS ; ADJUST AND ISSUE ERROR MESSAGE
C
      REST = 2.D0*W(2) + W(4) + W(5)
C
      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
C
      CCASO4I  = MIN (W(2),W(6))
      FRSO4I   = MAX (W(2) - CCASO4I, ZERO)
      CAFRI    = MAX (W(6) - CCASO4I, ZERO)
      CCANO32I = MIN (CAFRI, 0.5D0*W(4))
      CAFRI    = MAX (CAFRI - CCANO32I, ZERO)
      NO3FRI   = MAX (W(4) - 2.D0*CCANO32I, ZERO)
      CCACL2I  = MIN (CAFRI, 0.5D0*W(5))
      CLFRI    = MAX (W(5) - 2.D0*CCACL2I, ZERO)
      REST1    = 2.D0*FRSO4I + NO3FRI + CLFRI
C
      CNA2SO4I = MIN (FRSO4I, 0.5D0*W(1))
      FRSO4I   = MAX (FRSO4I - CNA2SO4I, ZERO)
      NAFRI    = MAX (W(1) - 2.D0*CNA2SO4I, ZERO)
      CNACLI   = MIN (NAFRI, CLFRI)
      NAFRI    = MAX (NAFRI - CNACLI, ZERO)
      CLFRI    = MAX (CLFRI - CNACLI, ZERO)
      CNANO3I  = MIN (NAFRI, NO3FRI)
      NO3FR    = MAX (NO3FRI - CNANO3I, ZERO)
      REST2    = 2.D0*FRSO4I + NO3FRI + CLFRI
C
      CMGSO4I  = MIN (FRSO4I, W(8))
      FRMGI    = MAX (W(8) - CMGSO4I, ZERO)
      FRSO4I   = MAX (FRSO4I - CMGSO4I, ZERO)
      CMGNO32I = MIN (FRMGI, 0.5D0*NO3FRI)
      FRMGI    = MAX (FRMGI - CMGNO32I, ZERO)
      NO3FRI   = MAX (NO3FRI - 2.D0*CMGNO32I, ZERO)
      CMGCL2I  = MIN (FRMGI, 0.5D0*CLFRI)
      CLFRI    = MAX (CLFRI - 2.D0*CMGCL2I, ZERO)
      REST3    = 2.D0*FRSO4I + NO3FRI + CLFRI
C
         IF (W(6).GT.REST) THEN                       ! Ca > 2*SO4+CL+NO3 ?
             W(6) = (ONE-1D-6)*REST              ! Adjust Ca amount
             W(1)= ZERO                          ! Adjust Na amount
             W(7)= ZERO                          ! Adjust K amount
             W(8)= ZERO                          ! Adjust Mg amount
             CALL PUSHERR (0051, 'ISRP4F')       ! Warning error: Ca, Na, K, Mg in excess
C
         ELSE IF (W(1).GT.REST1) THEN                 ! Na > 2*FRSO4+FRCL+FRNO3 ?
             W(1) = (ONE-1D-6)*REST1             ! Adjust Na amount
             W(7)= ZERO                          ! Adjust K amount
             W(8)= ZERO                          ! Adjust Mg amount
             CALL PUSHERR (0052, 'ISRP4F')       ! Warning error: Na, K, Mg in excess
C
         ELSE IF (W(8).GT.REST2) THEN                 ! Mg > 2*FRSO4+FRCL+FRNO3 ?
             W(8) = (ONE-1D-6)*REST2             ! Adjust Mg amount
             W(7)= ZERO                          ! Adjust K amount
             CALL PUSHERR (0053, 'ISRP4F')       ! Warning error: K, Mg in excess
C
         ELSE IF (W(7).GT.REST3) THEN                 ! K > 2*FRSO4+FRCL+FRNO3 ?
             W(7) = (ONE-1D-6)*REST3             ! Adjust K amount
             CALL PUSHERR (0054, 'ISRP4F')       ! Warning error: K in excess
         ENDIF
      ENDIF
C
C *** CALCULATE RATIOS *************************************************
C
      SO4RAT  = (W(1)+W(3)+W(6)+W(7)+W(8))/W(2)
      CRNARAT = (W(1)+W(6)+W(7)+W(8))/W(2)
      CRRAT   = (W(6)+W(7)+W(8))/W(2)
C
C *** FIND CALCULATION REGIME FROM (SO4RAT, CRNARAT, CRRAT, RRH) ********
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) POOR: R(Cr+Na)<2
C
      IF (2.0.LE.SO4RAT .AND. CRNARAT.LT.2.0) THEN
C
         SCASE = 'O7'
         CALL CALCO7                 ! Only liquid (metastable)

C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
      ELSEIF (SO4RAT.GE.2.0 .AND. CRNARAT.GE.2.0) THEN
C
       IF (CRRAT.LE.2.0) THEN
C
         SCASE = 'M8'
         CALL CALCM8                 ! Only liquid (metastable)

C        CALL CALCHCO3
C
C *** SULFATE POOR: Rso4>2; (DUST + SODIUM) RICH: R(Cr+Na)>2; DUST POOR: Rcr<2.
C
        ELSEIF (CRRAT.GT.2.0) THEN
C
         SCASE = 'P13'
         CALL CALCP13                 ! Only liquid (metastable)

C        CALL CALCHCO3
       ENDIF
C
C *** SULFATE RICH (NO ACID): 1<Rso4<2;
C
      ELSEIF (1.0.LE.SO4RAT .AND. SO4RAT.LT.2.0) THEN
C
         SCASE = 'L9'
         CALL CALCL9                ! Only liquid (metastable)

C
      CALL CALCNHA                ! MINOR SPECIES: HNO3, HCl
      CALL CALCNH3                !                NH3
C
C *** SULFATE SUPER RICH (FREE ACID): Rso4<1;
C
      ELSEIF (SO4RAT.LT.1.0) THEN
C
         SCASE = 'K4'
         CALL CALCK4                 ! Only liquid (metastable)

C
      CALL CALCNHA                  ! MINOR SPECIES: HNO3, HCl
      CALL CALCNH3                  !                NH3
C
      ENDIF
C
      RETURN
      END
C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCA2
C *** CASE A2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT >= 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED TOWARDS X, THE
C     AMOUNT OF HYDROGEN IONS (H+) FOUND IN THE LIQUID PHASE.
C     FOR EACH ESTIMATION OF H+, FUNCTION FUNCB2A CALCULATES THE
C     CONCENTRATION OF IONS FROM THE NH3(GAS) - NH4+(LIQ) EQUILIBRIUM.
C     ELECTRONEUTRALITY IS USED AS THE OBJECTIVE FUNCTION.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCA2
      INCLUDE 'isrpia.inc'
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU    =.TRUE.       ! Outer loop activity calculation flag
      OMELO     = TINY        ! Low  limit: SOLUTION IS VERY BASIC
      OMEHI     = 2.0D0*W(2)  ! High limit: FROM NH4+ -> NH3(g) + H+(aq)
C
C *** CALCULATE WATER CONTENT *****************************************
C
      MOLAL(5) = W(2)
      MOLAL(6) = ZERO
      CALL CALCMR
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = OMEHI
      Y1 = FUNCA2 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (OMEHI-OMELO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = MAX(X1-DX, OMELO)
         Y2 = FUNCA2 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
      IF (ABS(Y2).LE.EPS) THEN
         RETURN
      ELSE
         CALL PUSHERR (0001, 'CALCA2')    ! WARNING ERROR: NO SOLUTION
         RETURN
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCA2 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCA2')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCA2 (X3)
      RETURN
C
C *** END OF SUBROUTINE CALCA2 ****************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCA2
C *** CASE A2 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE A2 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCA2.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCA2 (OMEGI)
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA
C
C *** SETUP PARAMETERS ************************************************
C
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI    = W(2)         ! INITIAL AMOUNT OF (NH4)2SO4 IN SOLUTION
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
         A1    = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         A2    = XK2*R*TEMP/XKW*(GAMA(8)/GAMA(9))**2.
         A3    = XKW*RH*WATER*WATER
C
         LAMDA = PSI/(A1/OMEGI+ONE)
         ZETA  = A3/OMEGI
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = OMEGI                                        ! HI
C    ! slc.question
         MOLAL (5) = MAX(PSI-LAMDA,TINY)                          ! SO4I
         MOLAL (3) = MAX(W(3)/(ONE/A2/OMEGI + ONE), 2.*MOLAL(5))  ! NH4I
         MOLAL (6) = LAMDA                                        ! HSO4I
         GNH3      = MAX (W(3)-MOLAL(3), TINY)                    ! NH3GI
         COH       = ZETA                                         ! OHI
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
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    DENOM = (2.0*MOLAL(5)+MOLAL(6))
      FUNCA2= (MOLAL(3)/DENOM - ONE) + MOLAL(1)/DENOM
      RETURN
C
C *** END OF FUNCTION FUNCA2 ********************************************
C
      END
C=======================================================================

C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB4
C *** CASE B4 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. LIQUID AEROSOL PHASE ONLY POSSIBLE
C
C     FOR CALCULATIONS, A BISECTION IS PERFORMED WITH RESPECT TO H+.
C     THE OBJECTIVE FUNCTION IS THE DIFFERENCE BETWEEN THE ESTIMATED H+
C     AND THAT CALCULATED FROM ELECTRONEUTRALITY.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCB4
      INCLUDE 'isrpia.inc'
C
C *** SOLVE EQUATIONS **************************************************
C
      FRST       = .TRUE.
      CALAIN     = .TRUE.
      CALAOU     = .TRUE.
C
C *** CALCULATE WATER CONTENT ******************************************
C
      CALL CALCB1A         ! GET DRY SALT CONTENT, AND USE FOR WATER.
      MOLALR(13) = CLC       
      MOLALR(9)  = CNH4HS4   
      MOLALR(4)  = CNH42S4   
      CLC        = ZERO
      CNH4HS4    = ZERO
      CNH42S4    = ZERO
      WATER      = MOLALR(13)/M0(13)+MOLALR(9)/M0(9)+MOLALR(4)/M0(4)
C
      MOLAL(3)   = W(3)   ! NH4I
C
      DO 20 I=1,NSWEEP
         AK1   = XK1*((GAMA(8)/GAMA(7))**2.)*(WATER/GAMA(7))
         BET   = W(2)
         GAM   = MOLAL(3)
C
         BB    = BET + AK1 - GAM
         CC    =-AK1*BET
         DD    = BB*BB - 4.D0*CC
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (5) = MAX(TINY,MIN(0.5*(-BB + SQRT(DD)), W(2))) ! SO4I
         MOLAL (6) = MAX(TINY,MIN(W(2)-MOLAL(5),W(2)))         ! HSO4I
         MOLAL (1) = MAX(TINY,MIN(AK1*MOLAL(6)/MOLAL(5),W(2))) ! HI
         CALL CALCMR                                           ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT
20    CONTINUE
C
30    RETURN
C
C *** END OF SUBROUTINE CALCB4 ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCB1A
C *** CASE B1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH
C     2. THERE IS NO LIQUID PHASE
C     3. SOLIDS POSSIBLE: LC, { (NH4)2SO4  XOR  NH4HSO4 } (ONE OF TWO
C                         BUT NOT BOTH)
C
C     A SIMPLE MATERIAL BALANCE IS PERFORMED, AND THE AMOUNT OF LC
C     IS CALCULATED FROM THE (NH4)2SO4 AND NH4HSO4 WHICH IS LEAST
C     ABUNDANT (STOICHIMETRICALLY). THE REMAINING EXCESS OF SALT
C     IS MIXED WITH THE LC.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON
C UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
      SUBROUTINE CALCB1A
      INCLUDE 'isrpia.inc'
C
C *** SETUP PARAMETERS ************************************************
C
      X = 2*W(2)-W(3)       ! Equivalent NH4HSO4
      Y = W(3)-W(2)         ! Equivalent (NH4)2SO4
C
C *** CALCULATE COMPOSITION *******************************************
C
      IF (X.LE.Y) THEN      ! LC is the MIN (x,y)
         CLC     = X        ! NH4HSO4 >= (NH4)2S04
         CNH4HS4 = ZERO
         CNH42S4 = Y-X
      ELSE
         CLC     = Y        ! NH4HSO4 <  (NH4)2S04
         CNH4HS4 = X-Y
         CNH42S4 = ZERO
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE CALCB1
C ******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCC2
C *** CASE C2 
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCC2
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA, KAPA
C
      CALAOU =.TRUE.         ! Outer loop activity calculation flag
      FRST   =.TRUE.
      CALAIN =.TRUE.
C
C *** SOLVE EQUATIONS **************************************************
C
      LAMDA  = W(3)           ! NH4HSO4 INITIALLY IN SOLUTION
      PSI    = W(2)-W(3)      ! H2SO4 IN SOLUTION
      DO 20 I=1,NSWEEP
         PARM  = WATER*XK1/GAMA(7)*(GAMA(8)/GAMA(7))**2.
         BB    = PSI+PARM
         CC    =-PARM*(LAMDA+PSI)
         KAPA  = 0.5*(-BB+SQRT(BB*BB-4.0*CC))
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL(1) = PSI+KAPA                               ! HI
         MOLAL(3) = LAMDA                                  ! NH4I
         MOLAL(5) = KAPA                                   ! SO4I
         MOLAL(6) = MAX(LAMDA+PSI-KAPA, TINY)              ! HSO4I
         CH2SO4   = MAX(MOLAL(5)+MOLAL(6)-MOLAL(3), ZERO)  ! Free H2SO4
         CALL CALCMR                                       ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
         IF (.NOT.CALAIN) GOTO 30
         CALL CALCACT     
20    CONTINUE
C 
30    RETURN
C    
C *** END OF SUBROUTINE CALCC2 *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD3
C *** CASE D3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. THERE IS OLNY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCD3
      INCLUDE 'isrpia.inc'
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL CALCD1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4NO3               ! Save from CALCD1 run
      CHI2 = CNH42S4
      CHI3 = GHNO3
      CHI4 = GNH3
C
      PSI1 = CNH4NO3               ! ASSIGN INITIAL PSI's
      PSI2 = CHI2
      PSI3 = ZERO   
      PSI4 = ZERO  
C
      MOLAL(5) = PSI2              ! Include initial amount in water calc
      MOLAL(6) = ZERO
      MOLAL(3) = PSI1
      MOLAL(7) = PSI1
      CALL CALCMR                  ! Initial water
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      PSI4LO = TINY                ! Low  limit
      PSI4HI = CHI4                ! High limit
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
60    X1 = PSI4LO
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y1 = FUNCD3 (X1)
      IF (ABS(Y1).LE.EPS) RETURN
      YLO= Y1                 ! Save Y-value at HI position
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI4HI-PSI4LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCD3 (X2)
         IF (((Y1) .LT. ZERO) .AND. ((Y2) .GT. ZERO)) GOTO 20  ! (Y1*Y2.LT.ZERO) (slc.1.2012)
C         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION FOUND 
C
      YHI= Y1                      ! Save Y-value at Hi position
      IF (ABS(Y2) .LT. EPS) THEN   ! X2 IS A SOLUTION 
         RETURN
C
C *** { YLO, YHI } < 0.0 THE SOLUTION IS ALWAYS UNDERSATURATED WITH NH3
C Physically I dont know when this might happen, but I have put this
C branch in for completeness. I assume there is no solution; all NO3 goes to the
C gas phase.
C
      ELSE IF (YLO.LT.ZERO .AND. YHI.LT.ZERO) THEN
         P4 = TINY ! PSI4LO ! CHI4
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         YY = FUNCD3(P4)
         GOTO 50
C
C *** { YLO, YHI } > 0.0 THE SOLUTION IS ALWAYS SUPERSATURATED WITH NH3
C This happens when Sul.Rat. = 2.0, so some NH4+ from sulfate evaporates
C and goes to the gas phase ; so I redefine the LO and HI limits of PSI4
C and proceed again with root tracking.
C
      ELSE IF (YLO.GT.ZERO .AND. YHI.GT.ZERO) THEN
         PSI4HI = PSI4LO
         PSI4LO = PSI4LO - 0.1*(PSI1+PSI2) ! No solution; some NH3 evaporates
         IF (PSI4LO.LT.-(PSI1+PSI2)) THEN
            CALL PUSHERR (0001, 'CALCD3')  ! WARNING ERROR: NO SOLUTION
            RETURN
         ELSE
            MOLAL(5) = PSI2              ! Include sulfate in initial water calculation
            MOLAL(6) = ZERO
            MOLAL(3) = PSI1
            MOLAL(7) = PSI1
            CALL CALCMR                  ! Initial water
            GOTO 60                        ! Redo root tracking
         ENDIF
      ENDIF
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCD3 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*ABS(X1)) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCD3')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCD3 (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
      RETURN
C
C *** END OF SUBROUTINE CALCD3 ******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** FUNCTION FUNCD3
C *** CASE D3 
C     FUNCTION THAT SOLVES THE SYSTEM OF EQUATIONS FOR CASE D3 ; 
C     AND RETURNS THE VALUE OF THE ZEROED FUNCTION IN FUNCD3.
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCD3 (P4)
      INCLUDE 'isrpia.inc'
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
      FRST   = .TRUE.
      CALAIN = .TRUE.
      PSI4   = P4
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
         A2   = XK7*(WATER/GAMA(4))**3.0
         A3   = XK4*R*TEMP*(WATER/GAMA(10))**2.0
         A4   = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
         A7   = XKW *RH*WATER*WATER
C
         PSI3 = A3*A4*CHI3*(CHI4-PSI4) - PSI1*(2.D0*PSI2+PSI1+PSI4)
         PSI3 = PSI3/(A3*A4*(CHI4-PSI4) + 2.D0*PSI2+PSI1+PSI4) 
         PSI3 = MIN(MAX(PSI3, ZERO), CHI3)
C
         BB   = PSI4 - PSI3
CCCOLD         AHI  = 0.5*(-BB + SQRT(BB*BB + 4.d0*A7)) ! This is correct also
CCC         AHI  =2.0*A7/(BB+SQRT(BB*BB + 4.d0*A7)) ! Avoid overflow when HI->0
         DENM = BB+SQRT(BB*BB + 4.d0*A7)
         IF (DENM.LE.TINY) THEN       ! Avoid overflow when HI->0
            ABB  = ABS(BB)
            DENM = (BB+ABB) + 2.0*A7/ABB ! Taylor expansion of SQRT
         ENDIF
         AHI = 2.0*A7/DENM
C
C *** SPECIATION & WATER CONTENT ***************************************
C
         MOLAL (1) = AHI                             ! HI
         MOLAL (3) = PSI1 + PSI4 + 2.D0*PSI2         ! NH4I
         MOLAL (5) = PSI2                            ! SO4I
         MOLAL (6) = ZERO                            ! HSO4I
         MOLAL (7) = PSI3 + PSI1                     ! NO3I
         CNH42S4   = CHI2 - PSI2                     ! Solid (NH4)2SO4
         CNH4NO3   = ZERO                            ! Solid NH4NO3
         GHNO3     = CHI3 - PSI3                     ! Gas HNO3
         GNH3      = CHI4 - PSI4                     ! Gas NH3
         CALL CALCMR                                 ! Water content
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
C *** CALCULATE OBJECTIVE FUNCTION ************************************
C
20    CONTINUE
CCC      FUNCD3= NH4I/HI/MAX(GNH3,TINY)/A4 - ONE 
      FUNCD3= MOLAL(3)/MOLAL(1)/MAX(GNH3,TINY)/A4 - ONE 
      RETURN
C
C *** END OF FUNCTION FUNCD3 ********************************************
C
      END
C========================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCD1A
C *** CASE D1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4NO3
C
C     THE SOLID (NH4)2SO4 IS CALCULATED FROM THE SULFATES, WHILE NH4NO3
C     IS CALCULATED FROM NH3-HNO3 EQUILIBRIUM. 'ZE' IS THE AMOUNT OF
C     NH4NO3 THAT VOLATIZES WHEN ALL POSSILBE NH4NO3 IS INITIALLY IN
C     THE SOLID PHASE.
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCD1A
      INCLUDE 'isrpia.inc'
C
C *** SETUP PARAMETERS ************************************************
C
      PARM    = XK10/(R*TEMP)/(R*TEMP)
C
C *** CALCULATE NH4NO3 THAT VOLATIZES *********************************
C
      CNH42S4 = W(2)                                    
      X       = MAX(ZERO, MIN(W(3)-2.0*CNH42S4, W(4)))  ! MAX NH4NO3
      PS      = MAX(W(3) - X - 2.0*CNH42S4, ZERO)
      OM      = MAX(W(4) - X, ZERO)
C
      OMPS    = OM+PS
      DIAK    = SQRT(OMPS*OMPS + 4.0*PARM)              ! DIAKRINOUSA
      ZE      = MIN(X, 0.5*(-OMPS + DIAK))              ! THETIKI RIZA
C
C *** SPECIATION *******************************************************
C
      CNH4NO3 = X  - ZE    ! Solid NH4NO3
      GNH3    = PS + ZE    ! Gas NH3
      GHNO3   = OM + ZE    ! Gas HNO3
C
      RETURN
C
C *** END OF SUBROUTINE CALCD1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCG5
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCG5
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.   
      CHI1   = 0.5*W(1)
      CHI2   = MAX (W(2)-CHI1, ZERO)
      CHI3   = ZERO
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
      CHI5   = W(4)
      CHI6   = W(5)
C 
      PSI1   = CHI1
      PSI2   = CHI2
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = CHI2/M0(4) + CHI1/M0(2)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCG5A (X1)
      IF (CHI6.LE.TINY) GOTO 50  
ccc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
ccc      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCG5A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) THEN
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCG5A (PSI6LO)
      ENDIF
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCG5A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCG5')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCG5A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
         MOLAL(6) = DELTA                               ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCG5 *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCG5A
C *** CASE G5
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM POOR (SODRAT < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCG5A (X)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /CASEG/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA,
     &               PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7,
     &               A1,   A2,   A3,   A4,   A5,   A6,   A7
C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A2  = XK7 *(WATER/GAMA(4))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      AKK = A4*A6
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = 2.0D0*PSI1                          ! NAI
      MOLAL (3) = 2.0*PSI2 + PSI4                     ! NH4I
      MOLAL (4) = PSI6                                ! CLI
      MOLAL (5) = PSI2 + PSI1                         ! SO4I
      MOLAL (6) = ZERO
      MOLAL (7) = PSI5                                ! NO3I
C
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)              ! Gas NH3
      GHNO3     = MAX(CHI5 - PSI5, TINY)              ! Gas HNO3
      GHCL      = MAX(CHI6 - PSI6, TINY)              ! Gas HCl
C
      CNH42S4   = ZERO                                ! Solid (NH4)2SO4
      CNH4NO3   = ZERO                                ! Solid NH4NO3
      CNH4CL    = ZERO                                ! Solid NH4Cl
C
      CALL CALCMR                                     ! Water content
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    FUNCG5A = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCG5A *******************************************
C
      END

C========================================================================
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCH6
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCH6
      INCLUDE 'isrpia.inc'
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
      CALAOU = .TRUE.   
      CHI1   = W(2)                                ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)       
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY                  
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCH6A (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50  
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX 
         Y2 = FUNCH6A (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) THEN
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCH6A (PSI6LO)
      ENDIF
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCH6A (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCH6')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCH6A (X3)
C 
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCH6 ******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE FUNCH6A
C *** CASE H6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; SODIUM RICH (SODRAT >= 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NH4CL, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCH6A (X)
      INCLUDE 'isrpia.inc'
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
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8 
      FRST   = .TRUE.
      CALAIN = .TRUE. 
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
      A7  = XK8 *(WATER/GAMA(1))**2.0
      A8  = XK9 *(WATER/GAMA(3))**2.0
      A9  = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MAX(PSI5, TINY)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = BB*BB-4.d0*CC
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(PSI4,CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
      MOLAL (3) = PSI4                                  ! NH4I
      MOLAL (4) = PSI6 + PSI7                           ! CLI
      MOLAL (5) = PSI2 + PSI1                           ! SO4I
      MOLAL (6) = ZERO                                  ! HSO4I
      MOLAL (7) = PSI5 + PSI8                           ! NO3I
C
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C 
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = MAX(CHI7 - PSI7, ZERO)
      CNANO3    = MAX(CHI8 - PSI8, ZERO)
      CNA2SO4   = MAX(CHI1 - PSI1, ZERO) 
C
      CALL CALCMR                                    ! Water content
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    FUNCH6A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCH6A *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI6
C *** CASE I6
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : (NH4)2SO4, NA2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCI6
      INCLUDE 'isrpia.inc'
C
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL CALCI1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCI1 run
      CHI2 = CLC    
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC   
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A6 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      BB   = PSI2 + PSI4 + PSI5 + A6                    ! PSI6
      CC   =-A6*(PSI2 + PSI3 + PSI1)
      DD   = BB*BB - 4.D0*CC
      PSI6 = 0.5D0*(-BB + SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (1) = PSI6                                    ! HI
      MOLAL (2) = 2.D0*PSI4 + PSI3                        ! NAI
      MOLAL (3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1            ! NH4I
      MOLAL (5) = PSI2 + PSI4 + PSI5 + PSI6               ! SO4I
      MOLAL (6) = PSI2 + PSI3 + PSI1 - PSI6               ! HSO4I
      CLC       = ZERO
      CNAHSO4   = ZERO
      CNA2SO4   = CHI4 - PSI4
      CNH42S4   = ZERO
      CNH4HS4   = ZERO
      CALL CALCMR                                         ! Water content
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
C *** END OF SUBROUTINE CALCI6 *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCI1A
C *** CASE I1 ; SUBCASE 1
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCI1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CNA2SO4 = 0.5D0*W(1)
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      FRSO4   = MAX(W(2)-CNA2SO4, ZERO)
C
      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
C
      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
      ENDIF
C
C *** CALCULATE GAS SPECIES *********************************************
C
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCI1A *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCJ3
C *** CASE J3
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, FREE ACID (SULRAT < 1.0)
C     2. THERE IS ONLY A LIQUID PHASE
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE CALCJ3
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      LAMDA  = MAX(W(2) - W(3) - W(1), TINY)  ! FREE H2SO4
      CHI1   = W(1)                           ! NA TOTAL as NaHSO4
      CHI2   = W(3)                           ! NH4 TOTAL as NH4HSO4
      PSI1   = CHI1
      PSI2   = CHI2                           ! ALL NH4HSO4 DELIQUESCED
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A3 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      BB   = A3+LAMDA                        ! KAPA
      CC   =-A3*(LAMDA + PSI1 + PSI2)
      DD   = BB*BB-4.D0*CC
      KAPA = 0.5D0*(-BB+SQRT(DD))
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (1) = LAMDA + KAPA                 ! HI
      MOLAL (2) = PSI1                         ! NAI
      MOLAL (3) = PSI2                         ! NH4I
      MOLAL (4) = ZERO                         ! CLI
      MOLAL (5) = KAPA                         ! SO4I
      MOLAL (6) = LAMDA + PSI1 + PSI2 - KAPA   ! HSO4I
      MOLAL (7) = ZERO                         ! NO3I
C
      CNAHSO4   = ZERO
      CNH4HS4   = ZERO
C
      CALL CALCMR                              ! Water content
C
C *** CALCULATE ACTIVITIES OR TERMINATE INTERNAL LOOP *****************
C
      IF (FRST.AND.CALAOU .OR. .NOT.FRST.AND.CALAIN) THEN
         CALL CALCACT     
      ELSE
         GOTO 50
      ENDIF
10    CONTINUE
C 
50    RETURN
C
C *** END OF SUBROUTINE CALCJ3 ******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCO7
C *** CASE O7
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCO7
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,
     &               PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,
     &               A5, A6, A7, A8, A9
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU = .TRUE.
      CHI9   = MIN (W(6), W(2))                     ! CCASO4
      SO4FR  = MAX (W(2)-CHI9, ZERO)
      CAFR   = MAX (W(6)-CHI9, ZERO)
      CHI7   = MIN (0.5D0*W(7), SO4FR)              ! CK2SO4
      FRK    = MAX (W(7) - 2.D0*CHI7, ZERO)
      SO4FR  = MAX (SO4FR - CHI7, ZERO)
      CHI1   = MIN (0.5D0*W(1), SO4FR)              ! NA2SO4
      NAFR   = MAX (W(1) - 2.D0*CHI1, ZERO)
      SO4FR  = MAX (SO4FR - CHI1, ZERO)
      CHI8   = MIN (W(8), SO4FR)                    ! CMGSO4
      FRMG    = MAX(W(8) - CHI8, ZERO)
      SO4FR   = MAX(SO4FR - CHI8, ZERO)
      CHI3   = ZERO
      CHI5   = W(4)
      CHI6   = W(5)
      CHI2   = MAX (SO4FR, ZERO)
      CHI4   = MAX (W(3)-2.D0*CHI2, ZERO)
C
      PSI1   = CHI1
      PSI2   = CHI2
      PSI3   = ZERO
      PSI4   = ZERO
      PSI5   = ZERO
      PSI6   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
      WATER  = CHI2/M0(4) + CHI1/M0(2) + CHI7/M0(17) + CHI8/M0(21)
      WATER  = MAX (WATER , TINY)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCO7 (X1)
      IF (CHI6.LE.TINY) GOTO 50
ccc      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
ccc      IF (WATER .LE. TINY) RETURN                    ! No water
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCO7 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) THEN
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCO7 (PSI6LO)
      ENDIF 
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCO7 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCO7')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCO7 (X3)
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN  ! If quadrat.called
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                    ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                    ! SO4  EFFECT
         MOLAL(6) = DELTA                               ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCO7 *******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE FUNCO7
C *** CASE O7
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. (Rsulfate > 2.0 ; R(Cr+Na) < 2.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: NH4NO3, NH4CL, (NH4)2SO4, MgSO4, NA2SO4, K2SO4
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCO7 (X)
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA
      COMMON /CASEO/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5,
     &               PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4,
     &               A5, A6, A7, A8, A9
C
C *** SETUP PARAMETERS ************************************************
C
      PSI6   = X
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
C
C
      IF (CHI5.GE.TINY) THEN
         PSI5 = PSI6*CHI5/(A6/A5*(CHI6-PSI6) + PSI6)
         PSI5 = MIN (PSI5,CHI5)
      ELSE
         PSI5 = TINY
      ENDIF
C
CCC      IF(CHI4.GT.TINY) THEN
      IF(W(2).GT.TINY) THEN       ! Accounts for NH3 evaporation
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6) - 2.d0*PSI2/A4
         DD   = MAX(BB*BB-4.d0*CC,ZERO)           ! Patch proposed by Uma Shankar, 19/11/01
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MAX (MIN (PSI4,CHI4), ZERO)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL (2) = 2.0D0*PSI1                       ! Na+
      MOLAL (3) = 2.0D0*PSI2 + PSI4                ! NH4I
      MOLAL (4) = PSI6                             ! CLI
      MOLAL (5) = PSI1+PSI2+PSI7+PSI8              ! SO4I
      MOLAL (6) = ZERO                             ! HSO4
      MOLAL (7) = PSI5                             ! NO3I
      MOLAL (8) = ZERO                             ! CaI
      MOLAL (9) = 2.0D0*PSI7                       ! KI
      MOLAL (10)= PSI8                             ! Mg
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
CCC      MOLAL (1) = MAX(CHI5 - PSI5, TINY)*A5/PSI5   ! HI
       SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &             -MOLAL(9)-2.D0*MOLAL(10)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C
C *** CALCULATE GAS / SOLID SPECIES (LIQUID IN MOLAL ALREADY) *********
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNA2SO4  = ZERO
      CNH42S4  = ZERO
      CNH4NO3  = ZERO
      CNH4Cl   = ZERO
      CK2SO4   = ZERO
      CMGSO4   = ZERO
      CCASO4   = CHI9
C
C *** CALCULATE MOLALR ARRAY, WATER AND ACTIVITIES **********************
C
      CALL CALCMR
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
20    FUNCO7 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
CCC         FUNCG5A = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCO7 *******************************************
C
      END
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCM8
C *** CASE M8
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCM8
      INCLUDE 'isrpia.inc'
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
      CALAOU = .TRUE.
      CHI11  = MIN (W(6), W(2))                    ! CCASO4
      SO4FR  = MAX(W(2)-CHI11, ZERO)
      CAFR   = MAX(W(6)-CHI11, ZERO)
      CHI9   = MIN (0.5D0*W(7), SO4FR)             ! CK2S04
      FRK    = MAX(W(7)-2.D0*CHI9, ZERO)
      SO4FR  = MAX(SO4FR-CHI9, ZERO)
      CHI10  = MIN (W(8), SO4FR)                  ! CMGSO4
      FRMG   = MAX(W(8)-CHI10, ZERO)
      SO4FR  = MAX(SO4FR-CHI10, ZERO)
      CHI1   = MAX (SO4FR,ZERO)                    ! CNA2SO4
      CHI2   = ZERO                                ! CNH42S4
      CHI3   = ZERO                                ! CNH4CL
      FRNA   = MAX (W(1)-2.D0*CHI1, ZERO)
      CHI8   = MIN (FRNA, W(4))                    ! CNANO3
      CHI4   = W(3)                                ! NH3(g)
      CHI5   = MAX (W(4)-CHI8, ZERO)               ! HNO3(g)
      CHI7   = MIN (MAX(FRNA-CHI8, ZERO), W(5))    ! CNACL
      CHI6   = MAX (W(5)-CHI7, ZERO)               ! HCL(g)
C
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCM8 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCM8 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) THEN
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCM8 (PSI6LO)
      ENDIF
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCM8 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCM8')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCM8 (X3)
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCM8 ******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE FUNCM8
C *** CASE M8
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr < 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: NH4NO3, NH4CL, NANO3, NACL, MgSO4, NA2SO4, K2SO4
C
C *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** WRITTEN BY  CHRISTOS FOUNTOUKIS & ATHANASIOS NENES
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCM8 (X)
      INCLUDE 'isrpia.inc'
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
      PSI6   = X
      PSI1   = CHI1
      PSI2   = ZERO
      PSI3   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = CHI9
      PSI10  = CHI10
      PSI11  = ZERO
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
C      A1  = XK5 *(WATER/GAMA(2))**3.0
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
C      A7  = XK8 *(WATER/GAMA(1))**2.0
C      A8  = XK9 *(WATER/GAMA(3))**2.0
C      A11 = XK1*WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6+PSI7) - A6/A5*PSI8*(CHI6-PSI6-PSI3)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6-PSI3) + PSI6 + PSI7)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL (2) = PSI8 + PSI7 + 2.D0*PSI1               ! NAI
      MOLAL (3) = PSI4                                  ! NH4I
      MOLAL (4) = PSI6 + PSI7                           ! CLI
      MOLAL (5) = PSI2 + PSI1 + PSI9 + PSI10            ! SO4I
      MOLAL (6) = ZERO                                  ! HSO4I
      MOLAL (7) = PSI5 + PSI8                           ! NO3I
      MOLAL (8) = PSI11                                 ! CAI
      MOLAL (9) = 2.D0*PSI9                             ! KI
      MOLAL (10)= PSI10                                 ! MGI
C
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &            - MOLAL(9) - 2.D0*MOLAL(10)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH42S4   = ZERO
      CNH4NO3   = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CNA2SO4   = ZERO
      CK2SO4    = ZERO
      CMGSO4    = ZERO
      CCASO4    = CHI11
C
      CALL CALCMR                                    ! Water content
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
C20    FUNCM8 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
20    FUNCM8 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCM8 *******************************************
C
      END

C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCP13
C *** CASE P13
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
C                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCP13
      INCLUDE 'isrpia.inc'
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
      CALAOU  = .TRUE.
      CHI11   = MIN (W(2), W(6))                    ! CCASO4
      FRCA    = MAX (W(6) - CHI11, ZERO)
      FRSO4   = MAX (W(2) - CHI11, ZERO)
      CHI9    = MIN (FRSO4, 0.5D0*W(7))             ! CK2SO4
      FRK     = MAX (W(7) - 2.D0*CHI9, ZERO)
      FRSO4   = MAX (FRSO4 - CHI9, ZERO)
      CHI10   = FRSO4                               ! CMGSO4
      FRMG    = MAX (W(8) - CHI10, ZERO)
      CHI7    = MIN (W(1), W(5))                    ! CNACL
      FRNA    = MAX (W(1) - CHI7, ZERO)
      FRCL    = MAX (W(5) - CHI7, ZERO)
      CHI12   = MIN (FRCA, 0.5D0*W(4))              ! CCANO32
      FRCA    = MAX (FRCA - CHI12, ZERO)
      FRNO3   = MAX (W(4) - 2.D0*CHI12, ZERO)
      CHI17   = MIN (FRCA, 0.5D0*FRCL)              ! CCACL2
      FRCA    = MAX (FRCA - CHI17, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI17, ZERO)
      CHI15   = MIN (FRMG, 0.5D0*FRNO3)             ! CMGNO32
      FRMG    = MAX (FRMG - CHI15, ZERO)
      FRNO3   = MAX (FRNO3 - 2.D0*CHI15, ZERO)
      CHI16   = MIN (FRMG, 0.5D0*FRCL)              ! CMGCL2
      FRMG    = MAX (FRMG - CHI16, ZERO)
      FRCL    = MAX (FRCL - 2.D0*CHI16, ZERO)
      CHI8    = MIN (FRNA, FRNO3)                   ! CNANO3
      FRNA    = MAX (FRNA - CHI8, ZERO)
      FRNO3   = MAX (FRNO3 - CHI8, ZERO)
      CHI14   = MIN (FRK, FRCL)                     ! CKCL
      FRK     = MAX (FRK - CHI14, ZERO)
      FRCL    = MAX (FRCL - CHI14, ZERO)
      CHI13   = MIN (FRK, FRNO3)                    ! CKNO3
      FRK     = MAX (FRK - CHI13, ZERO)
      FRNO3   = MAX (FRNO3 - CHI13, ZERO)
C
      CHI5    = FRNO3                               ! HNO3(g)
      CHI6    = FRCL                                ! HCL(g)
      CHI4    = W(3)                                ! NH3(g)
C
      CHI3    = ZERO                                ! CNH4CL
      CHI1    = ZERO
      CHI2    = ZERO
C
      PSI6LO = TINY
      PSI6HI = CHI6-TINY    ! MIN(CHI6-TINY, CHI4)
C
C *** INITIAL VALUES FOR BISECTION ************************************
C
      X1 = PSI6LO
      Y1 = FUNCP13 (X1)
      IF (ABS(Y1).LE.EPS .OR. CHI6.LE.TINY) GOTO 50
C
C *** ROOT TRACKING ; FOR THE RANGE OF HI AND LO **********************
C
      DX = (PSI6HI-PSI6LO)/FLOAT(NDIV)
      DO 10 I=1,NDIV
         X2 = X1+DX
         Y2 = FUNCP13 (X2)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y2).LT.ZERO) GOTO 20  ! (Y1*Y2.LT.ZERO)
         X1 = X2
         Y1 = Y2
10    CONTINUE
C
C *** NO SUBDIVISION WITH SOLUTION; IF ABS(Y2)<EPS SOLUTION IS ASSUMED
C
      IF (ABS(Y2) .GT. EPS) THEN
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y2 = FUNCP13 (PSI6LO)
      ENDIF
      GOTO 50
C
C *** PERFORM BISECTION ***********************************************
C
20    DO 30 I=1,MAXIT
         X3 = 0.5*(X1+X2)
         CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
         Y3 = FUNCP13 (X3)
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
30    CONTINUE
      CALL PUSHERR (0002, 'CALCP13')    ! WARNING ERROR: NO CONVERGENCE
C
C *** CONVERGED ; RETURN **********************************************
C
40    X3 = 0.5*(X1+X2)
      CALL RSTGAMP            ! reinitialize activity coefficients (slc.1.2012)
      Y3 = FUNCP13 (X3)
C
C *** CALCULATE HSO4 SPECIATION AND RETURN *******************************
C
50    CONTINUE
      IF (MOLAL(1).GT.TINY .AND. MOLAL(5).GT.TINY) THEN
         CALL CALCHS4 (MOLAL(1), MOLAL(5), ZERO, DELTA)
         MOLAL(1) = MOLAL(1) - DELTA                     ! H+   EFFECT
         MOLAL(5) = MOLAL(5) - DELTA                     ! SO4  EFFECT
         MOLAL(6) = DELTA                                ! HSO4 EFFECT
      ENDIF
C
      RETURN
C
C *** END OF SUBROUTINE CALCP13 ******************************************
C
      END

C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE FUNCP13
C *** CASE P13
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE POOR (SULRAT > 2.0) ; Rcr+Na >= 2.0 ; Rcr > 2)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CaSO4
C     4. Completely dissolved: CA(NO3)2, CACL2, K2SO4, KNO3, KCL, MGSO4,
C                              MG(NO3)2, MGCL2, NANO3, NACL, NH4NO3, NH4CL
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      DOUBLE PRECISION FUNCTION FUNCP13 (X)
      INCLUDE 'isrpia.inc'
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
      PSI6   = X
      PSI1   = ZERO
      PSI2   = ZERO
      PSI3   = ZERO
      PSI4   = ZERO
      PSI7   = CHI7
      PSI8   = CHI8
      PSI9   = CHI9
      PSI10  = CHI10
      PSI11  = ZERO
      PSI12  = CHI12
      PSI13  = CHI13
      PSI14  = CHI14
      PSI15  = CHI15
      PSI16  = CHI16
      PSI17  = CHI17
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A4  = (XK2/XKW)*R*TEMP*(GAMA(10)/GAMA(5))**2.0
      A5  = XK4 *R*TEMP*(WATER/GAMA(10))**2.0
      A6  = XK3 *R*TEMP*(WATER/GAMA(11))**2.0
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      PSI5 = CHI5*(PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17) -
     &       A6/A5*(PSI8+2.D0*PSI12+PSI13+2.D0*PSI15)*(CHI6-PSI6)
      PSI5 = PSI5/(A6/A5*(CHI6-PSI6) + PSI6 + PSI7 + PSI14 +
     &       2.D0*PSI16 + 2.D0*PSI17)
      PSI5 = MIN(MAX(PSI5, TINY),CHI5)
C
      IF (W(3).GT.TINY .AND. WATER.GT.TINY) THEN  ! First try 3rd order soln
         BB   =-(CHI4 + PSI6 + PSI5 + 1.d0/A4)
         CC   = CHI4*(PSI5+PSI6)
         DD   = MAX(BB*BB-4.d0*CC,ZERO)
         PSI4 =0.5d0*(-BB - SQRT(DD))
         PSI4 = MIN(MAX(PSI4,ZERO),CHI4)
      ELSE
         PSI4 = TINY
      ENDIF
C
C *** CALCULATE SPECIATION *********************************************
C
      MOLAL (2) = PSI8 + PSI7                                     ! NAI
      MOLAL (3) = PSI4                                            ! NH4I
      MOLAL (4) = PSI6 + PSI7 + PSI14 + 2.D0*PSI16 + 2.D0*PSI17   ! CLI
      MOLAL (5) = PSI9 + PSI10                                    ! SO4I
      MOLAL (6) = ZERO                                            ! HSO4I
      MOLAL (7) = PSI5 + PSI8 + 2.D0*PSI12 + PSI13 + 2.D0*PSI15   ! NO3I
      MOLAL (8) = PSI11 + PSI12 + PSI17                           ! CAI
      MOLAL (9) = 2.D0*PSI9 + PSI13 + PSI14                       ! KI
      MOLAL (10)= PSI10 + PSI15 + PSI16                           ! MGI
C
C *** CALCULATE H+ *****************************************************
C
C      REST  = 2.D0*W(2) + W(4) + W(5)
CC
C      DELT1 = 0.0d0
C      DELT2 = 0.0d0
C      IF (W(1)+W(6)+W(7)+W(8).GT.REST) THEN
CC
CC *** CALCULATE EQUILIBRIUM CONSTANTS **********************************
CC
C      ALFA1 = XK26*RH*(WATER/1.0)                   ! CO2(aq) + H2O
C      ALFA2 = XK27*(WATER/1.0)                      ! HCO3-
CC
C      X     = W(1)+W(6)+W(7)+W(8) - REST            ! EXCESS OF CRUSTALS EQUALS CO2(aq)
CC
C      DIAK  = SQRT( (ALFA1)**2.0 + 4.0D0*ALFA1*X)
C      DELT1 = 0.5*(-ALFA1 + DIAK)
C      DELT1 = MIN ( MAX (DELT1, ZERO), X)
C      DELT2 = ALFA2
C      DELT2 = MIN ( DELT2, DELT1)
C      MOLAL(1) = DELT1 + DELT2                      ! H+
C      ELSE
C
C *** NO EXCESS OF CRUSTALS CALCULATE H+ *******************************
C
      SMIN      = 2.d0*MOLAL(5)+MOLAL(7)+MOLAL(4)-MOLAL(2)-MOLAL(3)
     &            - MOLAL(9) - 2.D0*MOLAL(10) - 2.D0*MOLAL(8)
      CALL CALCPH (SMIN, HI, OHI)
      MOLAL (1) = HI
C      ENDIF
C
      GNH3      = MAX(CHI4 - PSI4, TINY)
      GHNO3     = MAX(CHI5 - PSI5, TINY)
      GHCL      = MAX(CHI6 - PSI6, TINY)
C
      CNH4NO3   = ZERO
      CNH4CL    = ZERO
      CNACL     = ZERO
      CNANO3    = ZERO
      CK2SO4    = ZERO
      CMGSO4    = ZERO
      CCASO4    = CHI11
      CCANO32   = ZERO
      CKNO3     = ZERO
      CKCL      = ZERO
      CMGNO32   = ZERO
      CMGCL2    = ZERO
      CCACL2    = ZERO
C
      CALL CALCMR                                    ! Water content
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
C *** CALCULATE FUNCTION VALUE FOR OUTER LOOP ***************************
C
C20    FUNCP13 = MOLAL(3)*MOLAL(4)/GHCL/GNH3/A6/A4 - ONE
20    FUNCP13 = MOLAL(1)*MOLAL(4)/GHCL/A6 - ONE
C
      RETURN
C
C *** END OF FUNCTION FUNCP13 *******************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE CALCL9
C *** CASE L9
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SULRAT < 2.0)
C     2. SOLID & LIQUID AEROSOL POSSIBLE
C     3. SOLIDS POSSIBLE : CASO4
C     4. COMPLETELY DISSOLVED: NH4HSO4, NAHSO4, LC, (NH4)2SO4, KHSO4, MGSO4, NA2SO4, K2SO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCL9
      INCLUDE 'isrpia.inc'
      DOUBLE PRECISION LAMDA
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
C
C *** FIND DRY COMPOSITION **********************************************
C
      CALL CALCL1A
C
C *** SETUP PARAMETERS ************************************************
C
      CHI1 = CNH4HS4               ! Save from CALCL1 run
      CHI2 = CLC
      CHI3 = CNAHSO4
      CHI4 = CNA2SO4
      CHI5 = CNH42S4
      CHI6 = CK2SO4
      CHI7 = CMGSO4
      CHI8 = CKHSO4
C
      PSI1 = CNH4HS4               ! ASSIGN INITIAL PSI's
      PSI2 = CLC
      PSI3 = CNAHSO4
      PSI4 = CNA2SO4
      PSI5 = CNH42S4
      PSI6 = CK2SO4
      PSI7 = CMGSO4
      PSI8 = CKHSO4
C
      CALAOU = .TRUE.              ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A9 = XK1 *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.
C
C  CALCULATE DISSOCIATION QUANTITIES
C
      BB   = PSI7 + PSI6 + PSI5 + PSI4 + PSI2 + A9              ! LAMDA
      CC   = -A9*(PSI8 + PSI1 + PSI2 + PSI3)
      DD   = MAX(BB*BB - 4.D0*CC, ZERO)
      LAMDA= 0.5D0*(-BB + SQRT(DD))
      LAMDA= MIN(MAX (LAMDA, TINY), PSI8+PSI3+PSI2+PSI1)
C
C *** CALCULATE SPECIATION ********************************************
C
      MOLAL(1) = LAMDA                                            ! HI
      MOLAL(2) = 2.D0*PSI4 + PSI3                                 ! NAI
      MOLAL(3) = 3.D0*PSI2 + 2.D0*PSI5 + PSI1                     ! NH4I
      MOLAL(5) = PSI2 + PSI4 + PSI5 + PSI6 + PSI7 + LAMDA         ! SO4I
      MOLAL(6) = PSI2 + PSI3 + PSI1 + PSI8 - LAMDA                ! HSO4I
      MOLAL(9) = PSI8 + 2.0D0*PSI6                                ! KI
      MOLAL(10)= PSI7                                             ! MGI
C
      CLC      = ZERO
      CNAHSO4  = ZERO
      CNA2SO4  = ZERO
      CNH42S4  = ZERO
      CNH4HS4  = ZERO
      CK2SO4   = ZERO
      CMGSO4   = ZERO
      CKHSO4   = ZERO
C
      CALL CALCMR                                         ! Water content

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
C *** END OF SUBROUTINE CALCL9 *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCL1A
C *** CASE L1A
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE RICH, NO FREE ACID (1.0 <= SO4RAT < 2.0)
C     2. SOLID AEROSOL ONLY
C     3. SOLIDS POSSIBLE : K2SO4, CASO4, MGSO4, KHSO4, NH4HSO4, NAHSO4, (NH4)2SO4, NA2SO4, LC
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCL1A
      INCLUDE 'isrpia.inc'
C
C *** CALCULATE NON VOLATILE SOLIDS ***********************************
C
      CCASO4  = MIN (W(6), W(2))                    ! CCASO4
      FRSO4   = MAX(W(2) - CCASO4, ZERO)
      CAFR    = MAX(W(6) - CCASO4, ZERO)
      CK2SO4  = MIN (0.5D0*W(7), FRSO4)             ! CK2SO4
      FRK     = MAX(W(7) - 2.D0*CK2SO4, ZERO)
      FRSO4   = MAX(FRSO4 - CK2SO4, ZERO)
      CNA2SO4 = MIN (0.5D0*W(1), FRSO4)             ! CNA2SO4
      FRNA    = MAX(W(1) - 2.D0*CNA2SO4, ZERO)
      FRSO4   = MAX(FRSO4 - CNA2SO4, ZERO)
      CMGSO4  = MIN (W(8), FRSO4)                   ! CMGSO4
      FRMG    = MAX(W(8) - CMGSO4, ZERO)
      FRSO4   = MAX(FRSO4 - CMGSO4, ZERO)
C
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CNH42S4 = ZERO
      CKHSO4  = ZERO
C
      CLC     = MIN(W(3)/3.D0, FRSO4/2.D0)
      FRSO4   = MAX(FRSO4-2.D0*CLC, ZERO)
      FRNH4   = MAX(W(3)-3.D0*CLC,  ZERO)
C
      IF (FRSO4.LE.TINY) THEN
         CLC     = MAX(CLC - FRNH4, ZERO)
         CNH42S4 = 2.D0*FRNH4

      ELSEIF (FRNH4.LE.TINY) THEN
         CNH4HS4 = 3.D0*MIN(FRSO4, CLC)
         CLC     = MAX(CLC-FRSO4, ZERO)
C         IF (CK2SO4.GT.TINY) THEN
C            FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
C           CKHSO4 = 2.D0*FRSO4
C            CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
C         ENDIF
C         IF (CNA2SO4.GT.TINY) THEN
C            FRSO4   = MAX(FRSO4-CKHSO4/2.D0, ZERO)
C            CNAHSO4 = 2.D0*FRSO4
C            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
C         ENDIF
C
         IF (CNA2SO4.GT.TINY) THEN
            FRSO4  = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CNAHSO4 = 2.D0*FRSO4
            CNA2SO4 = MAX(CNA2SO4-FRSO4, ZERO)
         ENDIF
         IF (CK2SO4.GT.TINY) THEN
            FRSO4   = MAX(FRSO4-CNH4HS4/3.D0, ZERO)
            CKHSO4 = 2.D0*FRSO4
            CK2SO4 = MAX(CK2SO4-FRSO4, ZERO)
       ENDIF
      ENDIF
C
C *** CALCULATE GAS SPECIES ********************************************
C
      GHNO3 = W(4)
      GHCL  = W(5)
      GNH3  = ZERO
C
      RETURN
C
C *** END OF SUBROUTINE CALCL1A *****************************************
C
      END
C
C
C=======================================================================
C
C *** ISORROPIA CODE II
C *** SUBROUTINE CALCK4
C *** CASE K4
C
C     THE MAIN CHARACTERISTICS OF THIS REGIME ARE:
C     1. SULFATE SUPER RICH, FREE ACID (SO4RAT < 1.0)
C     2. THERE IS BOTH A LIQUID & SOLID PHASE
C     3. SOLIDS POSSIBLE : CASO4
C
C *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY CHRISTOS FOUNTOUKIS AND ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CALCK4
      INCLUDE 'isrpia.inc'
C
      DOUBLE PRECISION LAMDA, KAPA
      COMMON /CASEK/ CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3,
     &               A1,   A2,   A3,   A4
C
C *** SETUP PARAMETERS ************************************************
C
      CALAOU =.TRUE.               ! Outer loop activity calculation flag
      FRST   = .TRUE.
      CALAIN = .TRUE.
C
      CHI1   = W(3)                !  Total NH4 initially as NH4HSO4
      CHI2   = W(1)                !  Total NA initially as NaHSO4
      CHI3   = W(7)                !  Total K initially as KHSO4
      CHI4   = W(8)                !  Total Mg initially as MgSO4
C
      LAMDA  = MAX(W(2) - W(3) - W(1) - W(6) - W(7) - W(8), TINY)  ! FREE H2SO4
      PSI1   = CHI1                            ! ALL NH4HSO4 DELIQUESCED
      PSI2   = CHI2                            ! ALL NaHSO4 DELIQUESCED
      PSI3   = CHI3                            ! ALL KHSO4 DELIQUESCED
      PSI4   = CHI4                            ! ALL MgSO4 DELIQUESCED
C
C *** SOLVE EQUATIONS ; WITH ITERATIONS FOR ACTIVITY COEF. ************
C
      DO 10 I=1,NSWEEP
C
      A4 = XK1  *WATER/GAMA(7)*(GAMA(8)/GAMA(7))**2.0
C
      BB   = A4+LAMDA+PSI4                               ! KAPA
      CC   =-A4*(LAMDA + PSI3 + PSI2 + PSI1) + LAMDA*PSI4
      DD   = MAX(BB*BB-4.D0*CC, ZERO)
      KAPA = 0.5D0*(-BB+SQRT(DD))
C
C *** SAVE CONCENTRATIONS IN MOLAL ARRAY ******************************
C
      MOLAL (1) = MAX(LAMDA + KAPA, TINY)                         ! HI
      MOLAL (2) = PSI2                                            ! NAI
      MOLAL (3) = PSI1                                            ! NH4I
      MOLAL (5) = MAX(KAPA + PSI4, ZERO)                          ! SO4I
      MOLAL (6) = MAX(LAMDA + PSI1 + PSI2 + PSI3 - KAPA, ZERO)    ! HSO4I
      MOLAL (9) = PSI3                                            ! KI
      MOLAL (10)= PSI4                                            ! MGI
C
      CNH4HS4 = ZERO
      CNAHSO4 = ZERO
      CKHSO4  = ZERO
      CCASO4  = W(6)
      CMGSO4  = ZERO
C
      CALL CALCMR                                      ! Water content
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
C *** END OF SUBROUTINE CALCK4
C
      END

