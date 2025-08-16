C=======================================================================
C  DISEASE, Subroutine, Gustavo de Angelo Luca e Izael Martins Fattori Jr
C  Universidade de São Paulo - ESALQ USP
C  
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  07/17/2023 Written.
C  11/15/2024 Revised.
C  05/20/2025 Fungicide logic.      
C  08/10/2025 Logic/robustness fixes (cohorts, IR, LAF, LWD, cum. LAI)
C-----------------------------------------------------------------------
      SUBROUTINE DISEASE_LEAF (DYNAMIC,
     &    CONTROL, ISWITCH, Tmin, Tmax, RH, LAI_TOTAL,    ! Input
     &    ESP_LAT_HIST, SUP_INF_LIST, LAI_INF_LIST,       ! Input/State
     &    YRDOY, YREMRG, NVEG0, YREND,                    ! Input
     &    DISEASE_LAI)                                    ! Output
C-----------------------------------------------------------------------
          USE ModuleDefs     
          IMPLICIT NONE
          EXTERNAL F_RH, F_LWD, T_DEV, F_IR, F_DS, F_CANSPO, F_DEW,
     &             F_TAVG, F_LR, F_LS, F_PS, F_LAF, F_LAR, F_PPSR,
     &             F_PPSR_POP, APPLY_FUNGICIDE, CALC_DVIP,
     &             READ_DISEASE_PARAMETERS
          SAVE
C-----------------------------------------------------------------------
C  Switches / constants
          LOGICAL, PARAMETER :: USE_FUNGICIDE  = .FALSE.
          LOGICAL, PARAMETER :: USE_WTH_RH     = .TRUE.
          REAL,    PARAMETER :: LAI_MIN_START  = 0.5
          INTEGER, PARAMETER :: MAXDAYS        = 200
          REAL,    PARAMETER :: EPS            = 1.0E-6
C-----------------------------------------------------------------------
          
          INTEGER DAS, DYNAMIC
          
          REAL RH, LWD, Tdew, Es, E 
          
          REAL FT_D, FT, T, Tmin, Tmax
          REAL FT_G
          REAL TMIN_G, TOT_G, TMAX_G
          REAL TMIN_D, TOT_D, TMAX_D
      
          REAL IR
          REAL FSS, LAI
          REAL DS
          REAL LR
          REAL LS
          REAL IS
          REAL LA, LAF, LESIONAGEOPT
          REAL PREV_IS
          REAL LAI_TOTAL, DISEASE_LAI
          
          INTEGER k, DAE, DAE_START
          INTEGER KMAX, DAE_IDX, PREV_IDX
          
          INTEGER YRDOY, YREMRG, PLANT_LIVE, NVEG0, YREND
          INTEGER DISEASE_LIVE
      
          REAL NDS, LESION_S, KVERHULST, RVERHULST
          REAL YMAX, COF_A, COF_B
          REAL LDMIN, LESLIFEMAX
          REAL Lesion_Rate
          REAL ESP_INOC_SEC
          REAL HEALTH_LAI, HEALTH_LAI_AVAIL, NEW_LOSS_TODAY
          REAL INF_AREA_K
          
C--------- Population (individuals) & potential rate per area -----------
          REAL SPORES_YEST, PREV_LATENTS, INF_COUNT_PREV, NPREV_POP
          REAL POT_SPO_PER_AREA, PS_K
          
          REAL, DIMENSION(MAXDAYS,5) :: ESP_LAT_HIST 
          REAL, DIMENSION(MAXDAYS)   :: SUP_INF_LIST
          REAL, DIMENSION(MAXDAYS)   :: LAI_INF_LIST
          
          REAL,    DIMENSION(MAXDAYS) :: ADMITTED_AREA
          
          INTEGER DVIP_pts(7)  
          INTEGER idx, SUM7, BufferDays, ResidualDays, NSprays
          LOGICAL FungActive
          INTEGER DVIP_today
          REAL FUNG_EFFICIENCY
          
!-----------------------------------------------------------------------
!         Constructed types
          
          TYPE (ControlType) CONTROL
          TYPE (SwitchType)  ISWITCH
          
          DYNAMIC = CONTROL % DYNAMIC
          DAS     = CONTROL % DAS

!***********************************************************************
!  RUNINIT — called once (season start)
!***********************************************************************
          
      IF (DYNAMIC .EQ. RUNINIT) THEN

          IS = 0.0
          PLANT_LIVE = 0
          DAE = 0

          ESP_LAT_HIST = 0.0
          SUP_INF_LIST = 0.0
          LAI_INF_LIST = 0.0
          ADMITTED_AREA= 0.0
          
          DISEASE_LAI  = 0.0
          
          idx = 1
          DVIP_pts(:) = 0
          SUM7 = 0
          BufferDays = 0
          ResidualDays = 0
          FungActive = .FALSE.
          NSprays = 0
          FUNG_EFFICIENCY = 0.723

          SPORES_YEST = 0.0

          OPEN(23, FILE='C:\DSSAT48\Soybean\DISEASE_DEVELOPMENT.OUT',
     &        STATUS='replace')
          WRITE(23, 24)
   24     FORMAT('     RUN   YYDOY     DAS  LAI_HEALTH    LA_DISEASE',
     &     '  LA_INFECT   NEW_LOSS       LWD         RH      FAT_TEMP',
     &     '     LAI_TOTAL   SUM7  NSPRAYS FUNG_ACT')

!***********************************************************************
!  RATE — called every day
!***********************************************************************
          
      ELSEIF (DYNAMIC .EQ. RATE) THEN

!----- Emergence gate ---------------------------------------------------
          
          IF ((CONTROL%DAS .GT. NVEG0) .AND. (PLANT_LIVE .EQ. 0)) THEN
               DAE = 1
               PLANT_LIVE   = 1
               DISEASE_LIVE = 0
               
               IS = 0.0
               ESP_LAT_HIST = 0.0
               SUP_INF_LIST = 0.0
               LAI_INF_LIST = 0.0
               ADMITTED_AREA= 0.0
               SPORES_YEST  = 0.0
               
               CALL READ_DISEASE_PARAMETERS(NDS, LESION_S, KVERHULST,
     &            YMAX, COF_A, COF_B, RVERHULST, TMIN_G, TOT_G, TMAX_G,
     &            TMIN_D, TOT_D, TMAX_D, LDmin, LESIONAGEOPT,LESLIFEMAX,
     &            DAE_START)
          END IF
              
!----- Disease activation (starts after DAE_START) ----------------------
          
          IF ((PLANT_LIVE .EQ. 1) .AND. (DAE .GE. DAE_START)) THEN 
              DISEASE_LIVE = 1
          END IF
          
          IF (DISEASE_LIVE .EQ. 1) THEN
              
              CALL F_TAVG(Tmax, Tmin, T)
              CALL F_DEW (T, Tmin, Tmax, Tdew)
              IF (.NOT. USE_WTH_RH) THEN
                 CALL F_RH(RH, Tdew, Es, E, T)
              END IF
              CALL F_LWD(RH, LWD)
              
!----- Healthy LAI from yesterday’s cumulative loss ---------------------
              
              IF (DAE .GT. 1) THEN
                 PREV_IDX   = MIN(DAE-1,MAXDAYS)
                 HEALTH_LAI = LAI_TOTAL - LAI_INF_LIST(PREV_IDX)
              ELSE
                 HEALTH_LAI = LAI_TOTAL
                 PREV_IDX   = 1
              END IF
              HEALTH_LAI = MAX(0.0, HEALTH_LAI)

!----- Fungicide decision module (optional) -----------------------------
              
              CALL CALC_DVIP(LWD, T, DVIP_today)
              CALL APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7,
     &            BufferDays, FungActive, ResidualDays, NSprays,
     &            USE_FUNGICIDE, HEALTH_LAI)
              
!----- Thermal response and inoculum deposition -------------------------
              
              CALL T_DEV(T, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D,TMAX_D,
     &                  FT, FT_D, FT_G)
              CALL F_CANSPO(HEALTH_LAI, NDS, LESION_S, FSS)
              CALL F_DS(FSS, NDS, DS)
              
              CALL F_IR(FT, LWD, YMAX, COF_A, COF_B, IR)
              IF ((RH .LT. 85.0) .AND. (LWD .LT. 10.0)) IR = 0.1 * IR
              IF (FungActive) IR = IR * (1.0 - FUNG_EFFICIENCY)

!----- No healthy leaf → no new deposits today --------------------------
              
              IF (HEALTH_LAI .LE. 0.0) THEN
                 IR = 0.0
                 DS = 0.0
              END IF

              CALL F_LS(IR, DS, LS)

!----- Register today's latent spores (cohort at DAE) -------------------
              
              DAE_IDX = MIN(DAE,MAXDAYS)
              ESP_LAT_HIST(DAE_IDX, 1) = ESP_LAT_HIST(DAE_IDX, 1) + LS 

!----- Day accumulators --------------------------------------------------
              
              IS             = 0.0      ! infectious surface (m2 m-2) today
              ESP_INOC_SEC   = 0.0      ! total secondary spores today
              NEW_LOSS_TODAY = 0.0      ! newly activated surface today

!----- Yesterday carry-over for population/logistic ---------------------
              
              PREV_IS  = 0.0
              IF (DAE .GT. 1) PREV_IS = SUP_INF_LIST(PREV_IDX)

!----- Lesion aging rate (constant within the day) ----------------------
              
              CALL F_LAR (FT_D, LESLIFEMAX, Lesion_Rate)

!----- Build yesterday's population (individuals) -----------------------
              
              PREV_LATENTS = 0.0
              IF (DAE .GT. 1) THEN
                 KMAX = MIN(DAE-1, MAXDAYS)
                 DO k = 1, KMAX
                    IF (ESP_LAT_HIST(k,3) .EQ. 0.0) THEN
                       PREV_LATENTS = PREV_LATENTS + ESP_LAT_HIST(k,1)
                    END IF
                 END DO
              END IF

              INF_COUNT_PREV = 0.0
              IF (DAE .GT. 1) THEN
                 INF_COUNT_PREV = PREV_IS / MAX(LESION_S, EPS)
              END IF

              NPREV_POP = MAX(SPORES_YEST,0.0) + MAX(PREV_LATENTS,0.0)
     &                   + MAX(INF_COUNT_PREV,0.0)

!----- Potential spore rate per unit sporulating area (yesterday) -------
              
              CALL F_PPSR_POP(KVERHULST, RVERHULST, NPREV_POP,
     &                        PREV_IS, HEALTH_LAI, POT_SPO_PER_AREA)

!----- Cohort loop -------------------------------------------------------
              
              KMAX = MIN(DAE,MAXDAYS)
              DO k = 1, KMAX

                  CALL F_LR (FT_D, LDmin, LR)
                  ESP_LAT_HIST(k, 2) = ESP_LAT_HIST(k, 2) + LR
                  
!--------- Latent -> infectious (happens once) --------------------------
                  
                  IF ((ESP_LAT_HIST(k, 2) .GE. 1.0) .AND.
     &                (ESP_LAT_HIST(k, 3) .EQ. 0.0)) THEN
                       ESP_LAT_HIST(k, 3) = 1.0
                       ! allocate effective area respecting remaining healthy leaf
                       INF_AREA_K = ESP_LAT_HIST(k,1) * LESION_S
                       INF_AREA_K = MAX(INF_AREA_K, 0.0)
                       HEALTH_LAI_AVAIL = MAX(HEALTH_LAI - 
     &                                         NEW_LOSS_TODAY,0.0)
                       ADMITTED_AREA(k)=MIN(INF_AREA_K,HEALTH_LAI_AVAIL)
                       NEW_LOSS_TODAY=NEW_LOSS_TODAY + ADMITTED_AREA(k)
                  END IF
                  
!--------- After infectious ------------------------------------------------
     
                  IF (ESP_LAT_HIST(k, 3) .EQ. 1.0) THEN 
                      ESP_LAT_HIST(k,4) = ESP_LAT_HIST(k,4) +Lesion_Rate
                      LA = ESP_LAT_HIST(k,4)

                      CALL F_LAF(LA, LAF, LESIONAGEOPT)

!--------- Infectious area of the cohort today (only while LA<=1) -------
                      
                      INF_AREA_K = 0.0
                      IF (LA .LE. 1.0) INF_AREA_K = ADMITTED_AREA(k)
                      IF (LAI_TOTAL .LT. LAI_MIN_START) INF_AREA_K = 0.0

                      IS = IS + INF_AREA_K
                      
!--------- Secondary production while LA<=1 and substrate remains -------
                      
                      IF (LA .LE. 1.0) THEN
                          LAI  = MAX(HEALTH_LAI - NEW_LOSS_TODAY, 0.0)
                          IF(LAI.GT.0.0.AND.POT_SPO_PER_AREA.GT.0.0)THEN
                              PS_K = POT_SPO_PER_AREA * INF_AREA_K * 
     &                        FT_D * LAF
                              ESP_INOC_SEC =ESP_INOC_SEC + MAX(PS_K,0.0)
                          END IF
                      END IF
                  END IF
              END DO

!----- Store infectious surface today & add secondary inoculum ----------
              
              SUP_INF_LIST(DAE_IDX) = IS

              IF (ESP_INOC_SEC .GT. 0.0) THEN
                  ESP_LAT_HIST(DAE_IDX, 1) = ESP_LAT_HIST(DAE_IDX, 1) + 
     &                                     ESP_INOC_SEC
                  ESP_LAT_HIST(DAE_IDX, 5) = ESP_LAT_HIST(DAE_IDX, 5) + 
     &                                     ESP_INOC_SEC
              END IF

!----- Pool of spores of the day (life=1 day) ---------------------------
              
              SPORES_YEST = ESP_INOC_SEC
              
!----- Cumulative removed LAI (m2 m-2) ----------------------------------
              
              IF (DAE .GT. 1) THEN
                 LAI_INF_LIST(DAE_IDX) = LAI_INF_LIST(PREV_IDX) + 
     &                                   NEW_LOSS_TODAY
              ELSE
                 LAI_INF_LIST(DAE_IDX) = NEW_LOSS_TODAY
              END IF

!----- Internal var uses cm2 m-2 (legacy), output will divide by 10000 --
              
              DISEASE_LAI = MIN(LAI_TOTAL,LAI_INF_LIST(DAE_IDX))*10000.0
          
          ELSE
              HEALTH_LAI = LAI_TOTAL
          END IF
          
!----- Advance time ------------------------------------------------------
          
          IF (PLANT_LIVE .EQ. 1) DAE = DAE + 1
          IF (YREND .EQ. YRDOY)  PLANT_LIVE = 0
          
!***********************************************************************
!  OUTPUT — called on output events
!***********************************************************************
          
      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN

            WRITE(23, '(3I8, 9F12.3, I8, L2)') CONTROL%RUN, YRDOY,
     &  CONTROL%DAS, HEALTH_LAI, DISEASE_LAI/10000.0, IS,NEW_LOSS_TODAY,
     &  LWD, RH, FT, LAI_TOTAL, REAL(SUM7), NSprays, FungActive

!***********************************************************************
!  SEASEND — called once (season end)
!***********************************************************************
            
      ELSEIF (DYNAMIC .EQ. SEASEND) THEN
            PLANT_LIVE   = 0
            DISEASE_LIVE = 0
            ESP_LAT_HIST = 0.0
            SUP_INF_LIST = 0.0
            LAI_INF_LIST = 0.0
            ADMITTED_AREA= 0.0
            DISEASE_LAI  = 0.0
            DAE = 1
            SPORES_YEST  = 0.0
            CLOSE(23)
      ENDIF

      END SUBROUTINE DISEASE_LEAF

!-----------------------------------------------------------------------
!  SUBROUTINES
!-----------------------------------------------------------------------

! ------ READ PARAMETERS
!  keep parameter file consistent with columns order below.
     
      SUBROUTINE READ_DISEASE_PARAMETERS(NDS, LESION_S, KVERHULST,
     &                                  YMAX, COF_A, COF_B, 
     &                                  RVERHULST, TMIN_G, TOT_G,TMAX_G,
     &                                  TMIN_D, TOT_D, TMAX_D,    
     &                                  LDmin, LESIONAGEOPT, LESLIFEMAX,
     &                                  DAE_START)

          IMPLICIT NONE
          
          REAL    NDS, LESION_S, KVERHULST, RVERHULST
          REAL    YMAX, COF_A, COF_B
          REAL    TMIN_G, TOT_G, TMAX_G
          REAL    TMIN_D, TOT_D, TMAX_D
          REAL    LDmin, LESIONAGEOPT, LESLIFEMAX
          INTEGER DAE_START

          CHARACTER(LEN=6)  ID       
          CHARACTER(LEN=20) VRNAME   
          INTEGER UNIT, IOS

          UNIT = 10

          OPEN(UNIT, FILE='C:\DSSAT48\Soybean\disease_parameters.txt', 
     &        STATUS='OLD', ACTION='READ', IOSTAT=IOS)

          READ(UNIT, *)
          READ(UNIT, *)
          READ(UNIT, *)
          READ(UNIT, *)
          READ(UNIT, *)

          READ(UNIT,*) ID, VRNAME, 
     &             NDS, LESION_S, KVERHULST, RVERHULST, 
     &             YMAX, COF_A, COF_B,                  
     &             TMIN_G, TOT_G, TMAX_G,               
     &             TMIN_D, TOT_D, TMAX_D,               
     &             LDmin, LESIONAGEOPT, LESLIFEMAX,
     &             DAE_START

          CLOSE(UNIT)

          RETURN
      END SUBROUTINE READ_DISEASE_PARAMETERS

! ------ AVERAGE TEMPERATURE

      SUBROUTINE F_TAVG(Tmax, Tmin, T)
          USE ModuleDefs
          IMPLICIT NONE
          REAL Tmax, Tmin, T
          T = (Tmax + Tmin) / 2.0
      END SUBROUTINE F_TAVG
      
! ------ DEW POINT TEMPERATURE (empirical)
      SUBROUTINE F_DEW(T, Tmin, Tmax, Tdew)
          USE ModuleDefs
          IMPLICIT NONE
          REAL T, Tmin, Tmax, Tdew
          Tdew = (-0.036*T) + (0.9679*Tmin)+(0.0072*(Tmax-Tmin)) +1.0111
      END SUBROUTINE F_DEW
      
! ------ RELATIVE HUMIDITY (Tetens form) -> RH in %
      SUBROUTINE F_RH(RH, Tdew, Es, E, T)
          USE ModuleDefs
          IMPLICIT NONE
          REAL RH, Tdew, Es, E, T
          Es = EXP(17.625*T    /(243.04+T))
          E  = EXP(17.625*Tdew /(243.04+Tdew))
          RH = MAX((E / Es) * 100.0, 0.0)
      END SUBROUTINE F_RH
      
! ------ LEAF WETNESS DURATION (no side-effect on RH)
      SUBROUTINE F_LWD(RH, LWD)
          USE ModuleDefs
          IMPLICIT NONE
          REAL RH, LWD
          REAL RH_LOC
          
          RH_LOC = RH
          IF (RH_LOC .LE. 1.0) RH_LOC = RH_LOC * 100.0
          RH_LOC = MIN(MAX(RH_LOC, 0.0), 100.0)
          
          LWD = MAX(31.31 / (1.0 + EXP(-((RH_LOC - 85.17) / 9.13))),0.0)
          LWD = MIN(LWD, 24.0)
      END SUBROUTINE F_LWD

! ------ TEMPERATURE DEVELOPMENT FUNCTION (beta response)
      SUBROUTINE T_DEV (T, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &                  FT, FT_D, FT_G)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FT_D, FT_G, FT, T
          REAL TMIN_D, TOT_D, TMAX_D, TMIN_G, TOT_G, TMAX_G
          
          FT_G = ((TMAX_G-T)/(TMAX_G-TOT_G)) * ((T-TMIN_G)/
     &           (TOT_G-TMIN_G))**((TOT_G-TMIN_G)/(TMAX_G-TOT_G)) 
          
          IF (T .GT. TMIN_D) THEN
              FT_D = ((TMAX_D-T)/(TMAX_D-TOT_D)) * ((T-TMIN_D)/
     &               (TOT_D-TMIN_D))**((TOT_D-TMIN_D)/(TMAX_D-TOT_D))
          ELSE
              FT_D = 0.0
          END IF
      
          IF ((T .LT. 12.0) .OR. (T .GT. 32.0)) THEN
              FT = 0.0
          ELSE
              FT = MAX(MIN(FT_G * FT_D, 1.0), 0.0)
          END IF
      END SUBROUTINE T_DEV

! ------ INFECTION RATE
      SUBROUTINE F_IR(FT, LWD, Ymax, A, B, IR)
          USE ModuleDefs
          IMPLICIT NONE
          REAL Ymax, A, B, IR, LWD, FT
          IR = Ymax * FT * (1.0 - EXP(-(A * LWD)**B))
      END SUBROUTINE F_IR
    
! ------ CANOPY SUPPORTED SPORES (capacity fraction)
      SUBROUTINE F_CANSPO(LAI, NDS, LESION_S, FSS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FSS, LAI, NDS, LESION_S
          IF (NDS .LE. 0.0) THEN
              FSS = 0.0
          ELSE
              FSS = LAI / (NDS * LESION_S)
              FSS = MIN(MAX(FSS, 0.0), 1.0)
          END IF
      END SUBROUTINE F_CANSPO

! ------ DEPOSITED SPORES
      SUBROUTINE F_DS(FSS, NDS, DS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FSS, NDS, DS
          DS = FSS * NDS
      END SUBROUTINE F_DS

! ------ LATENCY RATE (uses FT_D; no re-computation)
      SUBROUTINE F_LR(FT_D, LDmin, LR)
          USE ModuleDefs
          IMPLICIT NONE
          REAL LR, FT_D, LDmin
          LR = MAX(FT_D, 0.0) / MAX(LDmin, 1.0E-6)
      END SUBROUTINE F_LR
    
! ------ LATENT SPORES
      SUBROUTINE F_LS(IR, DS, LS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL IR, DS, LS
          LS = MAX(IR * DS, 0.0)
      END SUBROUTINE F_LS

! ------ SPORES PRODUCTION (secondary inoculum) — kept for compatibility
      SUBROUTINE F_PS(PPSR, LESION_S, FT_D, LAF, PS)
          USE ModuleDefs
          IMPLICIT NONE
          REAL PS, PPSR, LESION_S, FT_D, LAF
          PS = PPSR * LESION_S * FT_D * LAF
      END SUBROUTINE F_PS 
    
! ------ LESION AGE RATE
      SUBROUTINE F_LAR (FT_D, LESLIFEMAX, Lesion_Rate)
          USE ModuleDefs
          IMPLICIT NONE
          REAL FT_D, LESLIFEMAX, Lesion_Rate
          IF (FT_D .LE. 0.0) THEN
              Lesion_Rate = 0.0
          ELSE
              Lesion_Rate = FT_D / MAX(LESLIFEMAX, 1.0E-6)
          END IF
      END SUBROUTINE F_LAR
          
! ------ LESION AGE FACTOR (0..1)
      SUBROUTINE F_LAF(LA, LAF, LESIONAGEOPT)
          USE ModuleDefs
          IMPLICIT NONE
          REAL LAF, LA, LESIONAGEOPT, TMP
          IF (LA .LT. LESIONAGEOPT) THEN
              TMP = LA / MAX(LESIONAGEOPT, 1.0E-6)
          ELSEIF (LA .EQ. LESIONAGEOPT) THEN
              TMP = 1.0
          ELSE 
              TMP = 1.0 - (LA-LESIONAGEOPT)/MAX(1.0-LESIONAGEOPT,1.0E-6)
          END IF
          LAF = MAX(MIN(TMP, 1.0), 0.0)
      END SUBROUTINE F_LAF

! ------ POTENTIAL SPORE RATE (legacy logistic, area-based)
!  (not used in core; kept for compatibility)
      SUBROUTINE F_PPSR(KVERHULST, RVERHULST, INFECTIOUS_S, PREV_IS,LAI,
     &                  PPSR)
          USE ModuleDefs
          IMPLICIT NONE
          REAL PPSR, KVERHULST, RVERHULST, INFECTIOUS_S, PREV_IS, LAI
          REAL CARRY
          IF (INFECTIOUS_S .EQ. 0.0) THEN
              PPSR = 0.0   
          ELSE
              CARRY = KVERHULST * LAI
              PPSR = (PREV_IS + (RVERHULST * PREV_IS *(CARRY-PREV_IS)) /
     &               MAX(KVERHULST, 1.0E-6)) / INFECTIOUS_S
              PPSR = MAX(PPSR, 0.0)
          END IF
      END SUBROUTINE F_PPSR

! ------ POTENTIAL SPORE RATE PER UNIT SPORULATING AREA 

      SUBROUTINE F_PPSR_POP(KVERHULST, RVERHULST, NPREV, INF_SURF_PREV,
     &                       LAI_SUSC, POT_SPO_PER_AREA)
          USE ModuleDefs
          IMPLICIT NONE
          REAL KVERHULST, RVERHULST, NPREV, INF_SURF_PREV, LAI_SUSC
          REAL POT_SPO_PER_AREA
          REAL KTOTAL, DN, NEW_SPORES
          REAL EPSL
          EPSL = 1.0E-6

          POT_SPO_PER_AREA = 0.0

          IF (LAI_SUSC .LE. 0.0) RETURN
          KTOTAL = MAX(KVERHULST * LAI_SUSC, EPSL)

          IF (NPREV .LE. 0.0) RETURN

          DN         = RVERHULST * NPREV * (KTOTAL - NPREV) / KTOTAL
          NEW_SPORES = MAX(DN, 0.0)

          IF (INF_SURF_PREV .GT. 0.0) THEN
             POT_SPO_PER_AREA = NEW_SPORES / INF_SURF_PREV
          ELSE
             POT_SPO_PER_AREA = 0.0
          END IF
      END SUBROUTINE F_PPSR_POP
      
! ------ DVIP (Daily infection-probability index)
      SUBROUTINE CALC_DVIP(LWD, T, DVIP)
          USE ModuleDefs
          IMPLICIT NONE
          REAL LWD, T
          INTEGER DVIP
          DVIP = 0
          IF (T .LT. 15.0) THEN  
            IF (LWD .GT. 14.1) THEN        
              DVIP = 2
            ELSEIF (LWD .GE. 11.1 .AND. LWD .LE. 14.0) THEN 
              DVIP = 1
            ELSE
              DVIP = 0
            ENDIF
          ELSEIF (T .GE. 15.0 .AND. T .LT. 20.0) THEN 
            IF (LWD .GT. 17.1) THEN     
              DVIP = 3
            ELSEIF (LWD .GE. 13.1 .AND. LWD .LE. 17.0) THEN
              DVIP = 2
            ELSEIF (LWD .GE. 7.1  .AND. LWD .LE. 13.0) THEN 
              DVIP = 1
            ELSE
              DVIP = 0
            ENDIF
          ELSEIF (T .GE. 20.0 .AND. T .LT. 25.0) THEN 
            IF (LWD .GT. 17.1) THEN        
              DVIP = 3
            ELSEIF (LWD .GE. 10.1 .AND. LWD .LE. 17.0) THEN 
              DVIP = 2
            ELSEIF (LWD .GE. 7.1  .AND. LWD .LE. 10.0) THEN  
              DVIP = 1
            ELSE
              DVIP = 0
            ENDIF
          ELSE
            IF (LWD .GT. 17.1) THEN        
              DVIP = 3
            ELSEIF (LWD .GE. 11.1 .AND. LWD .LE. 17.0) THEN 
              DVIP = 2
            ELSEIF (LWD .GE. 7.1  .AND. LWD .LE. 11.0) THEN  
              DVIP = 1
            ELSE
              DVIP = 0
            ENDIF
          ENDIF 
      END SUBROUTINE CALC_DVIP

! ------ FUNGICIDE DECISION/APPLICATIONS
      SUBROUTINE APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7, 
     &             BufferDays, FungActive, ResidualDays, NSprays,
     &             USE_FUNGICIDE, HEALTH_LAI)
          USE ModuleDefs
          IMPLICIT NONE
          INTEGER DVIP_today, DVIP_pts(7), idx, SUM7
          INTEGER BufferDays, ResidualDays, NSprays
          LOGICAL FungActive, USE_FUNGICIDE
          REAL HEALTH_LAI
          LOGICAL CAN_SPRAY
          
          SUM7 = SUM7 - DVIP_pts(idx) + DVIP_today
          DVIP_pts(idx) = DVIP_today
          idx = MOD(idx,7) + 1
          
          IF (BufferDays .GT. 0) BufferDays = BufferDays - 1
          
          IF (FungActive) THEN
              ResidualDays = ResidualDays - 1
              IF (ResidualDays .LE. 0) FungActive = .FALSE.
          END IF
          
          CAN_SPRAY = (HEALTH_LAI .GE. 1.5)
          
          IF (SUM7 .GE. 6 .AND. BufferDays .EQ. 0 .AND. CAN_SPRAY) THEN
              NSprays = NSprays + 1
              BufferDays = 16
              
              IF (USE_FUNGICIDE) THEN
                  FungActive   = .TRUE.
                  ResidualDays = 14
              ELSE
                  FungActive   = .FALSE.
                  ResidualDays = 0 
              END IF
          END IF
      END SUBROUTINE APPLY_FUNGICIDE

!=======================================================================

!***********************************************************************
!  Variable listing (fully updated — 16 Aug 2025)
!***********************************************************************
! --------------------------- Arguments --------------------------------
! DYNAMIC           : DSSAT phase flag (RUNINIT/RATE/OUTPUT/SEASEND)
! CONTROL           : ControlType with %DAS (days after sowing), %RUN, etc.
! ISWITCH           : SwitchType (not used here; kept for interface parity)
! Tmin, Tmax (°C)   : Daily min/max air temperature
! RH (%)            : Relative humidity (either from weather or computed)
! LAI_TOTAL (m2 m-2): Canopy leaf area index (total)
! ESP_LAT_HIST(:,:) : State array (MAXDAYS x 5) for cohorts; columns:
!                     (1)=latent amount (spores m-2),
!                     (2)=latent progress (0..1),
!                     (3)=infectious flag (0/1),
!                     (4)=lesion age (d),
!                     (5)=secondary inoculum credited today (spores m-2)
! SUP_INF_LIST(:)   : Infectious surface by day (m2 m-2) → LA_INFECT
! LAI_INF_LIST(:)   : Cumulative removed LAI by day (m2 m-2) → LA_DISEASE
! YRDOY             : Current date (YYDOY)
! YREMRG            : Date of emergence (YYDOY) [not used internally]
! NVEG0             : DAS threshold for emergence gate
! YREND             : Harvest/end date (YYDOY)
! DISEASE_LAI       : OUTPUT — cumulative diseased area (cm2 m-2; printed as m2 m-2)

! --------------------------- Parameters -------------------------------
! USE_FUNGICIDE     : Toggle for applying fungicide (logic retained)
! USE_WTH_RH        : Use RH from weather (.TRUE.) or compute via dew point
! LAI_MIN_START     : Minimal LAI to allow activity (substrate safeguard)
! MAXDAYS           : Max internal history length (days)
! EPS               : Small epsilon for safe divisions

! --------------------------- Locals (scalars) -------------------------
! DAS               : Days after sowing (from CONTROL)
! T (°C)            : Daily mean air temperature
! FT                : Combined temperature response (0..1)
! FT_D, FT_G        : Post-infection and germination temperature responses
! IR (0..1)         : Infection rate for the day
! FSS (0..1)        : Canopy fraction supporting spores (capacity fraction)
! DS (spores m-2)   : Deposited spores today
! LR (d-1)          : Latency progress rate
! LS (spores m-2)   : Latent spores added to today's cohort
! IS (m2 m-2)       : Infectious surface today
! LA (d)            : Lesion age of a cohort
! LAF (0..1)        : Lesion-age factor
! PREV_IS (m2 m-2)  : Infectious surface yesterday
! LWD (h)           : Leaf wetness duration (capped at 24 h)
! Tdew (°C)         : Dew point temperature (for optional RH calc)
! Es, E             : Saturation and actual vapor pressure (Tetens; RH calc)
! HEALTH_LAI (m2 m-2)      : Healthy/susceptible LAI at day start
! HEALTH_LAI_AVAIL (m2 m-2): Remaining healthy LAI available for new lesions
! NEW_LOSS_TODAY (m2 m-2)  : Newly activated (removed) LAI today
! INF_AREA_K (m2 m-2)      : Infectious area allocated to cohort k today
! Lesion_Rate (d-1)        : Lesion aging rate
! ESP_INOC_SEC (spores m-2): Total secondary inoculum produced today
! SPORES_YEST (spores m-2) : Secondary spores produced yesterday (1-day pool)
! PREV_LATENTS (spores m-2): Sum of latent individuals (yesterday)
! INF_COUNT_PREV (#)       : Count of infectious lesions yesterday
! NPREV_POP (#)            : Population yesterday (spores + latents + lesions)
! POT_SPO_PER_AREA (sp m-2 d-1 per m2 m-2):
!                     Potential secondary spores per unit sporulating area
! PS_K (spores m-2)  : Secondary spores contributed by cohort k today
! NDS (spores m-2)   : Primary dispersed spores (from parameter file)
! LESION_S (m2)      : Average lesion surface area
! KVERHULST (#/LAI)  : Carrying capacity per unit LAI (population space)
! RVERHULST (d-1)    : Intrinsic logistic rate (population space)
! YMAX, COF_A, COF_B : Infection response parameters
! TMIN_G/TOT_G/TMAX_G (°C): Temperature response for germination
! TMIN_D/TOT_D/TMAX_D (°C): Temperature response for development
! LDMin (d)          : Minimum latency (used to scale LR)
! LESLIFEMAX (d)     : Max lesion lifespan (caps aging)
! DAE                : Days after emergence (internal counter)
! DAE_START          : DAE at which disease becomes active
! KMAX               : Upper bound for cohort loop (≤ MIN(DAE,MAXDAYS))
! DAE_IDX, PREV_IDX  : Indices for today and yesterday in history arrays
! PLANT_LIVE         : 1 while crop is alive (after emergence gate)
! DISEASE_LIVE       : 1 after disease activation gate (DAE ≥ DAE_START)
! DVIP_today         : Daily infection-probability class (0..3)
! idx (1..7)         : Circular index for 7-day DVIP buffer
! SUM7               : Sum of last seven DVIP classes
! BufferDays (d)     : Spray buffer to avoid back-to-back applications
! ResidualDays (d)   : Fungicide residual protection counter
! NSprays (#)        : Number of sprays applied
! FungActive (L)     : Whether fungicide residual is currently active
! FUNG_EFFICIENCY    : Proportional reduction in IR while active

! --------------------------- Locals (arrays) --------------------------
! ESP_LAT_HIST(MAXDAYS,5):
!   (1) latent amount (spores m-2) created on day k
!   (2) latent progress (0..1) accumulated for cohort k
!   (3) infectious flag for cohort k (0 or 1)
!   (4) lesion age (d) for cohort k (advances only after infectious)
!   (5) secondary inoculum tallied for cohort k on its creation day
!
! SUP_INF_LIST(MAXDAYS):
!   Infectious surface per day (m2 m-2); yesterday’s value used to derive
!   infectious counts and as denominator for population-based logistic.
!
! LAI_INF_LIST(MAXDAYS):
!   Cumulative removed/sick LAI per day (m2 m-2).
!
! ADMITTED_AREA(MAXDAYS):
!   Infectious area actually admitted to cohort k (m2 m-2), capped by
!   remaining healthy LAI at activation time.

! --------------------------- Subroutines ------------------------------
! F_TAVG      : Mean temperature
! F_DEW       : Empirical dew point (for optional RH path)
! F_RH        : RH from Tetens (uses T and Tdew)
! F_LWD       : Leaf wetness duration from RH (0..24 h)
! T_DEV       : Beta-type temperature responses (FT_G, FT_D, FT)
! F_IR        : Infection rate vs FT and LWD (Ymax, A, B)
! F_CANSPO    : Canopy spore support fraction (capacity)
! F_DS        : Deposited spores = FSS * NDS
! F_LR        : Latency progress rate (FT_D / LDMin)
! F_LS        : Latent spores = IR * DS (≥ 0)
! F_LAF       : Lesion-age factor (triangular shape peaking at LESIONAGEOPT)
! F_LAR       : Lesion aging rate (FT_D / LESLIFEMAX; 0 if FT_D≤0)
! F_PPSR      : Legacy area-based logistic (kept for compatibility)
! F_PS        : Legacy secondary inoculum per lesion (compatibility)
! F_PPSR_POP  : Population-space logistic → spores per unit sporulating area
! CALC_DVIP   : Daily infection-probability index class (0..3)
! APPLY_FUNGICIDE:
!   7-day risk sum trigger; optional spray; residual IR reduction window.

! --------------------------- Files/IO ---------------------------------
! Unit 23      : Writes daily diagnostics to
!                C:\DSSAT48\Soybean\DISEASE_DEVELOPMENT.OUT
! Output cols  : RUN, YYDOY, DAS, LAI_HEALTH, LA_DISEASE, LAI_INFECT,
!                NEW_LOSS, LWD, RH, FAT_TEMP(=FT), LAI_TOTAL, SUM7, 
!                NSPRAYS, FUNG_ACT (logical)
