C=======================================================================
C  DISEASE, Subroutine, Gustavo de Angelo Luca e Izael Martins Fattori Jr
C  Universidade de São Paulo - ESALQ USP
C  
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  07/17/2023 Written.
C  11/15/2024 Revised.      
C-----------------------------------------------------------------------
C  Called : 
C  Calls  : 
C========================================================================
      SUBROUTINE DISEASE_LEAF (DYNAMIC,
     &    CONTROL, ISWITCH, Tmin, Tmax,RH, LAI_TOTAL,    !Input
     &    ESP_LAT_HIST, SUP_INF_LIST, LAI_INF_LIST,      !Input
     &    YRDOY, YREMRG, NVEG0, YREND,                   !Input
     &    DISEASE_LAI)                                   !Output
                                                         !adicionar RH como parametro para usar rhum media do wth
C-----------------------------------------------------------------------
          !Definitions of constructed variable types
          USE ModuleDefs     
          IMPLICIT NONE
          EXTERNAL F_RH, F_LWD, T_DEV, F_IR, F_DS,F_CANSPO,F_DEW,F_TAVG,
     &    F_LR, F_LS, F_PS, F_LAF, F_LAR,F_PPSR,READ_DISEASE_PARAMETERS,
     &    APPLY_FUNGICIDE, CALC_DVIP 
          SAVE
C-----------------------------------------------------------------------
          LOGICAL, PARAMETER :: USE_FUNGICIDE = .FALSE.  !if true, protective effect of fungicide is active;
C-----------------------------------------------------------------------
          
          INTEGER DAS, DYNAMIC, RUN
          
          !Variables for Leaf Wetness Duration and Relative Humidity
          REAL RH, LWD, Tdew, Es, E 
          
          !Variables for Temperature Function
          REAL FT_D, FT, T, Tmin, Tmax
          REAL FT_G
          REAL TMIN_G
          REAL TOT_G
          REAL TMAX_G
          REAL TMIN_D
          REAL TOT_D
          REAL TMAX_D
      
          !Variable for infection rate
          REAL IR
      
          ! Variables for canopy suported spores 
          REAL FSS, LAI
      
          ! Variables for deposited spores
          REAL DS
      
          ! Variable for latency rate
          REAL LR
      
          ! Variable for latent spores
          REAL LS
              
          ! Variable for spores production
          REAL PPSR, IS, PS
      
          ! Variables for lesion age factor
          REAL LA, LAF, LESIONAGEOPT
      
          ! Variables for potencial spore production
          REAL PREV_IS, INFECTIOUS_S 
      
          ! Variables for secondary spores
          REAL ESP_INOC_SEC_LESION_DAY
      
          ! LAI calculations
          REAL LAI_TOTAL, DISEASE_LAI
          
          INTEGER k, DAE, DAE_START
          
          INTEGER YRDOY, YREMRG, PLANT_LIVE, NVEG0, YREND
          
          INTEGER DISEASE_LIVE
      
          ! Disease Variables 
          REAL NDS
          REAL LESION_S
          REAL KVERHULST
          REAL RVERHULST
          REAL YMAX
          REAL COF_A
          REAL COF_B
          REAL LDMIN
          REAL LESLIFEMAX
          REAL Lesion_Rate
          REAL ESP_INOC_SEC
          REAL HEALTH_LAI
          
          REAL, DIMENSION(200,5) :: ESP_LAT_HIST 
          REAL, DIMENSION(200) :: SUP_INF_LIST
          REAL, DIMENSION(200) :: LAI_INF_LIST
          
          ! Variables for Fungicide Control
          INTEGER DVIP_pts(7)  
          INTEGER idx          
          INTEGER SUM7         
          INTEGER BufferDays   
          INTEGER ResidualDays 
          INTEGER NSprays    
          LOGICAL FungActive
          INTEGER DVIP_today
          REAL FUNG_EFFICIENCY
          
!-----------------------------------------------------------------------
!         Define constructed variable types based on definitions in
!         ModuleDefs.for.
          
          TYPE (ControlType) CONTROL
          TYPE (SwitchType) ISWITCH
          
!         Transfer values from constructed data types into local variables.
          DYNAMIC = CONTROL % DYNAMIC
          DAS = CONTROL % DAS
          RUN = CONTROL % RUN

!***********************************************************************
!     Run Initialization - Called once per simulation
!***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
!-----------------------------------------------------------------------             

          IS = 0.
          INFECTIOUS_S = 0.
          PLANT_LIVE = 0

          ESP_LAT_HIST = 0.
          SUP_INF_LIST = 0.
          LAI_INF_LIST = 0.
          
          DISEASE_LAI = 0.
          
          idx = 1
          DVIP_pts(:) = 0
          SUM7 = 0
          BufferDays = 0
          ResidualDays = 0
          FungActive = .FALSE.
          NSprays = 0
          FUNG_EFFICIENCY = 0.723 !Circular tecnica embrapa 2023-2024 (Godoy et al.,2024)
          
          OPEN(23, FILE='C:\DSSAT48\Soybean\DISEASE_DEVELOPMENT.OUT',
     &        STATUS='replace')

          WRITE(23, 24)
   24     FORMAT(          
     x     '     RUN   YYDOY     DAS  LAI_HEALTH  ',
     x     'LA_DISEASE        LWD          RH    FAT_TEMP   LAI_TOTAL',
     x     '    SUM7  NSPRAYS FUNG_ACT')

 
      ELSEIF (DYNAMIC .EQ. RATE) THEN

!***********************************************************************
!     EMERGENCE CALCULATIONS - 
!***********************************************************************
          IF ((CONTROL%DAS .GT. NVEG0) .AND. 
     &    (PLANT_LIVE .EQ. 0)) THEN

               DAE = 1
               PLANT_LIVE = 1
               
               DISEASE_LIVE = 0.
               
               IS = 0.
               INFECTIOUS_S = 0.
      
               ESP_LAT_HIST = 0.
               SUP_INF_LIST = 0.
               LAI_INF_LIST = 0.
               
               CALL READ_DISEASE_PARAMETERS(NDS, LESION_S, KVERHULST,
     &        YMAX, COF_A, COF_B, RVERHULST, TMIN_G, TOT_G, TMAX_G,
     &        TMIN_D, TOT_D, TMAX_D, LDmin, LESIONAGEOPT, LESLIFEMAX,
     &        DAE_START)

          END IF
!-----------------------------------------------------------------------     
              
          IF ((PLANT_LIVE == 1) .AND. (DAE .GE. DAE_START)) THEN
              DISEASE_LIVE = 1
          END IF
          
          IF ((DISEASE_LIVE .EQ. 1)) THEN
              
              CALL F_TAVG(Tmax, Tmin, T)
              
              CALL F_DEW (T, Tmin, Tmax, Tdew)
              
              !comentar a call F_RH para usar a umidade relativa do wth
              !CALL F_RH(RH, Tdew, Es, E, T)
              
              CALL F_LWD(RH, LWD)
              
              IF (DAE > 1) THEN
                  HEALTH_LAI = LAI_TOTAL - LAI_INF_LIST(DAE-1)
              ELSE
                  HEALTH_LAI = LAI_TOTAL
              END IF
              
              HEALTH_LAI  = MAX(0.,HEALTH_LAI)
              
             !---------------- Fungicide Control --------------  
              CALL CALC_DVIP(LWD, T, DVIP_today)
              
              CALL APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7,
     &            BufferDays, FungActive, ResidualDays, NSprays,
     &            USE_FUNGICIDE, HEALTH_LAI)
             !-----------------------------------------------------
              
              CALL T_DEV(T,TMIN_G,TOT_G,TMAX_G,TMIN_D,TOT_D,TMAX_D,FT, 
     &        FT_D, FT_G)

              CALL F_CANSPO(HEALTH_LAI, NDS, LESION_S, FSS)

              CALL F_DS(FSS, NDS, DS)
              
              IF ((RH < 80.0) .AND. (LWD < 8.0)) THEN
                  IR = 0.1 * IR
              ELSE
                  CALL F_IR(FT, LWD, YMAX, COF_A, COF_B, IR)
              END IF
              
              !Block new infections if fungicide is active -> protective effect
              IF (FungActive) THEN
                  IR = IR * (1.0 - FUNG_EFFICIENCY)
              END IF

              CALL F_LS(IR, DS, LS)

              ESP_LAT_HIST(DAE, 1) = ESP_LAT_HIST(DAE, 1) + LS 

              IS = 0.0
              INFECTIOUS_S = 0.0
              ESP_INOC_SEC = 0.0

              DO k=1, DAE
                  CALL F_LR (FT_D, LDmin, LR, T, TMAX_D, TMIN_D, TOT_D)
                  ESP_LAT_HIST(DAE, 2) = ESP_LAT_HIST(DAE, 2) + LR
                  
                  CALL F_LAR (FT_D, LESLIFEMAX, Lesion_Rate)
                  ESP_LAT_HIST(DAE,4)= ESP_LAT_HIST(DAE,4) +Lesion_Rate 
                  LA = ESP_LAT_HIST (DAE,4) 
                  
                  IF ((ESP_LAT_HIST(DAE, 2) >= 1) .AND.
     &                (ESP_LAT_HIST(DAE, 3) == 0)) THEN
                  
                       ESP_LAT_HIST(DAE, 3) = 1
                  END IF
                    
                  IF (ESP_LAT_HIST(DAE, 3) == 1) THEN 
                      INFECTIOUS_S = ESP_LAT_HIST(DAE, 1) * LESION_S 
                      IS = IS + INFECTIOUS_S
                      
                      IF (LAI_TOTAL .LT. 0.5) THEN
                          IS = 0 !evita progressao elevada da doença quando iaf é mto baixo
                      END IF
                  END IF
                  
                  CALL F_LAF(LA, LAF, LESIONAGEOPT)
                  
                  IF (PREV_IS == 0.) THEN
                      IF (DAE > 1) THEN
                          PREV_IS = ESP_LAT_HIST (DAE - 1., 1)
                      ELSE
                          PREV_IS = 0.
                      END IF     
                  END IF
                  
                  CALL F_PPSR(KVERHULST, RVERHULST,INFECTIOUS_S,
     &            PREV_IS, LAI, PPSR) 
                 
            
                  IF (ESP_LAT_HIST(DAE,4) <= 1.0) THEN
                      CALL F_PS (PPSR,LESION_S, FT_D, LAF,PS) 
                      ESP_INOC_SEC_LESION_DAY = PS
            
                  ELSE
                      ESP_INOC_SEC_LESION_DAY = 0.0   
                  END IF
                  
                  ESP_INOC_SEC = ESP_INOC_SEC + ESP_INOC_SEC_LESION_DAY
                        
              END DO
              SUP_INF_LIST(DAE) = IS

              ! Add secondary spores to the daily spores 
              ESP_LAT_HIST(DAE, 1) = ESP_LAT_HIST(DAE, 1) + ESP_INOC_SEC
              ESP_LAT_HIST(DAE,5) = ESP_LAT_HIST(DAE,5) + ESP_INOC_SEC
            
              LAI_INF_LIST(DAE) = LAI_INF_LIST(DAE) + IS
        
              DISEASE_LAI = LAI_INF_LIST(DAE)
              DISEASE_LAI = min(LAI_TOTAL, DISEASE_LAI) 
              DISEASE_LAI = DISEASE_LAI * 10000. 
          
          ELSE
              HEALTH_LAI = LAI_TOTAL
          
          END IF
          
          IF (PLANT_LIVE .EQ. 1) THEN
              DAE = DAE + 1
          END IF
          
          IF (YREND .EQ. YRDOY) THEN
              PLANT_LIVE = 0
          END IF
          
!***********************************************************************
!     OUTPUT
!***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN

            ! Write output
            WRITE(23, '(3I8, 7F12.3, I8, L2)') CONTROL%RUN, YRDOY,
     &           CONTROL%DAS, HEALTH_LAI, DISEASE_LAI/10000., LWD, RH, 
     &           FT, LAI_TOTAL, REAL(SUM7), NSprays, FungActive

      ELSEIF (DYNAMIC .EQ. SEASEND) THEN

            PLANT_LIVE = 0
            DISEASE_LIVE = 0
            
            ESP_LAT_HIST = 0.
            SUP_INF_LIST = 0.
            LAI_INF_LIST = 0.
            
            DISEASE_LAI = 0.
            DAE = 1
            CLOSE(23)
        
!***********************************************************************
!     END OF DYNAMIC IF CONSTRUCT
!***********************************************************************
      ENDIF
!***********************************************************************
      END SUBROUTINE DISEASE_LEAF

!-----------------------------------------------------------------
! SUBROUTINES 
!-----------------------------------------------------------------

! ------ READ PARAMETERS
      SUBROUTINE READ_DISEASE_PARAMETERS(NDS, LESION_S,KVERHULST,
     &                                  YMAX, COF_A, COF_B, 
     &                                  RVERHULST,TMIN_G, TOT_G,TMAX_G,
     &                                  TMIN_D, TOT_D, TMAX_D,    
     &                                  LDmin, LESIONAGEOPT,LESLIFEMAX,
     &                                  DAE_START)

          IMPLICIT NONE
          
          REAL NDS, LESION_S, KVERHULST, RVERHULST
          REAL YMAX, COF_A, COF_B
          REAL TMIN_G, TOT_G, TMAX_G
          REAL TMIN_D, TOT_D, TMAX_D
          REAL LDmin, LESIONAGEOPT, LESLIFEMAX
          INTEGER DAE_START

          CHARACTER(LEN=6) ID       
          CHARACTER(LEN=20) VRNAME   
          INTEGER UNIT, IOS

          UNIT = 10

          OPEN(UNIT, FILE="C:\DSSAT48\Soybean\disease_parameters.txt", 
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
          
          T = (Tmax + Tmin) / 2.
          
      END SUBROUTINE F_TAVG
      
! ------ DEW POINT TEMPERATURE
      SUBROUTINE F_DEW(T, Tmin, Tmax, Tdew)
          USE ModuleDefs
          IMPLICIT NONE
          
          REAL T, Tmin, Tmax, Tdew
          
          Tdew = (-0.036*T)+(0.9679*Tmin)+(0.0072*(Tmax-Tmin))+1.0111
          
      END SUBROUTINE F_DEW
      
! ------ RELATIVE HUMIDITY
      SUBROUTINE F_RH(RH, Tdew, Es, E, T)
          USE ModuleDefs
          IMPLICIT NONE
          
          REAL RH, Tdew, Es, E, T
          
          Es = EXP(17.625*T/(243.04+T))
          E = EXP((17.625 * Tdew) / (243.04 + Tdew))
          RH = max((E / Es) * 100.0, 0.0)
          
      END SUBROUTINE F_RH
      
! ------ LEAF WETNESS DURATION
      SUBROUTINE F_LWD(RH, LWD)
          USE ModuleDefs
          IMPLICIT NONE

          REAL RH, LWD
          
          IF (RH <= 1.0) THEN
              RH = RH * 100.0
              RH = MIN(MAX(RH, 0.0), 100.0)
          END IF
          
          LWD = max(31.31 / (1.0 + exp(-((RH - 85.17) / 9.13))), 0.0)
           
      END SUBROUTINE F_LWD

! ------ TEMPERATURE DEVELOPMENT FUNCTION
      SUBROUTINE T_DEV 
     &   (T,TMIN_G,TOT_G,TMAX_G,TMIN_D,TOT_D,TMAX_D,FT, FT_D, FT_G)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL FT_D, FT_G, FT, T
          REAL TMIN_D, TOT_D, TMAX_D, TMIN_G, TOT_G, TMAX_G
          
          FT_G = ((TMAX_G-T)/(TMAX_G-TOT_G))*((T-TMIN_G)/
     &    (TOT_G-TMIN_G))**((TOT_G-TMIN_G)/(TMAX_G-TOT_G)) 
          
          IF (T > TMIN_D) THEN
              FT_D = ((TMAX_D-T)/(TMAX_D-TOT_D)) * ((T-TMIN_D)/
     &    (TOT_D-TMIN_D))**((TOT_D-TMIN_D)/(TMAX_D-TOT_D))
          ELSE
              FT_D = 0.
          END IF
      
          IF ((T < 12.0) .OR. (T > 32.0)) then
              FT = 0.0
          ELSE
              FT = max(min(FT_G * FT_D, 1.0), 0.0)
          END IF
            
      END SUBROUTINE T_DEV

! ------ INFECTION RATE
      SUBROUTINE F_IR(FT, LWD, Ymax, A, B, IR)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL Ymax, A, B, IR, LWD, FT
          
          IR = Ymax * FT * (1.0 - exp(-(A * LWD)**B))
            
      END SUBROUTINE F_IR
    
! ------ CANOPY SUPORTED SPORES
      SUBROUTINE F_CANSPO(LAI, NDS, LESION_S, FSS)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL FSS, LAI, NDS, LESION_S
          
          IF (NDS == 0.) THEN
              FSS = 0.
          ELSE
              FSS = LAI / (NDS * LESION_S)
              FSS = min(max(FSS, 0.0), 1.0)
          END IF
        
      END SUBROUTINE F_CANSPO

! ------ DEPOSITED SPORES
      SUBROUTINE F_DS(FSS, NDS,DS)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL FSS, NDS, DS
    
          DS = FSS * NDS
        
      END SUBROUTINE F_DS

! ------ LATENCY RATE
      SUBROUTINE F_LR(FT_D, LDmin, LR, T, TMAX_D, TMIN_D, TOT_D)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL LR, FT_D, LDmin, T, TMAX_D, TMIN_D, TOT_D
          
          
          IF (T > TMIN_D) THEN
              FT_D = ((TMAX_D-T)/(TMAX_D-TOT_D)) * ((T-TMIN_D)/
     &    (TOT_D-TMIN_D))**((TOT_D-TMIN_D)/(TMAX_D-TOT_D))
          ELSE
              FT_D = 0.
          END IF

          LR = FT_D / LDmin
        
      END SUBROUTINE F_LR
    
! ------ LATENT SPORES    
      SUBROUTINE F_LS(IR, DS, LS)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL IR, DS, LS
    
          LS = max((IR * DS), 0.0)
        
      END SUBROUTINE F_LS

! ------ SPORES PRODUCTION (2º INOCULUM PROD)
      SUBROUTINE F_PS(PPSR, LESION_S, FT_D, LAF,PS)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL PS, PPSR, LESION_S, FT_D, LAF

          PS = PPSR * LESION_S * FT_D * LAF

      END SUBROUTINE F_PS 
    
! ------- LESION AGE RATE
      SUBROUTINE F_LAR (FT_D, LESLIFEMAX, Lesion_Rate)
          USE ModuleDefs
          IMPLICIT NONE
          
          REAL FT_D, LESLIFEMAX, Lesion_Rate
          
          IF (FT_D .LE. 0.) THEN
              Lesion_Rate = 0.
          ELSE
              Lesion_Rate = 1.0 / (FT_D * LESLIFEMAX)
          END IF

      END SUBROUTINE F_LAR
          
! ------ LESION AGE FACTOR  
      SUBROUTINE F_LAF(LA, LAF, LESIONAGEOPT)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL LAF, LA, LESIONAGEOPT
            
          IF (LA < LESIONAGEOPT) THEN
              LAF = LA / LESIONAGEOPT
                
          ELSEIF (LA == LESIONAGEOPT) THEN
              LAF = 1.0
              
          ELSE 
              LAF = 1.0 - (LA - LESIONAGEOPT) / (1.0 - LESIONAGEOPT)
              
          END IF
        
      END SUBROUTINE F_LAF

! ------ POTENCIAL SPORES PRODUCTION RATE
      SUBROUTINE F_PPSR(KVERHULST,RVERHULST,INFECTIOUS_S, PREV_IS, LAI,
     &    PPSR)
          USE ModuleDefs
          IMPLICIT NONE
        
          REAL PPSR, KVERHULST, RVERHULST, INFECTIOUS_S, PREV_IS, 
     &    LAI
          
          IF (INFECTIOUS_S == 0.0) THEN
              PPSR = 0.0   
                    
          ELSE
              PPSR = (PREV_IS + (RVERHULST*PREV_IS*
     &    (KVERHULST*LAI-PREV_IS))/KVERHULST)/INFECTIOUS_S
                      
              PPSR = max(PPSR, 0.0)
          END IF

      END SUBROUTINE F_PPSR
      
! ------ DAILY VALUE OF INFECTION PROBABILITY - DVIP (Fungicide)
      SUBROUTINE CALC_DVIP(LWD, T, DVIP)
          USE ModuleDefs
          IMPLICIT NONE
      
          REAL LWD, T
          INTEGER DVIP
      
          DVIP = 0
          
          IF (T < 15.0) THEN  
          IF (LWD > 14.1) THEN        
            DVIP = 2
          ELSEIF (LWD >= 11.1 .AND. LWD <= 14.0) THEN 
            DVIP = 1
          ELSEIF (LWD < 11.0) THEN    
            DVIP = 0
          ENDIF

        ELSEIF (T >= 15.0 .AND. T < 20.0) THEN 
          IF (LWD > 17.1) THEN     
            DVIP = 3
          ELSEIF (LWD >= 13.1 .AND. LWD <= 17.0) THEN
            DVIP = 2
          ELSEIF (LWD >= 7.1 .AND. LWD <= 13.0) THEN 
            DVIP = 1
          ELSEIF (LWD < 7.0) THEN   
            DVIP = 0
          ENDIF

        ELSEIF (T >= 20.0 .AND. T < 25.0) THEN 
          IF (LWD > 17.1) THEN        
            DVIP = 3
          ELSEIF (LWD >= 10.1 .AND. LWD <= 17.0) THEN 
            DVIP = 2
          ELSEIF (LWD >= 7.1 .AND. LWD <= 10.0) THEN  
            DVIP = 1
          ELSEIF (LWD < 7.0) THEN     
            DVIP = 0
          ENDIF

        ELSEIF (T >= 25.0) THEN            
          IF (LWD > 17.1) THEN        
            DVIP = 3
          ELSEIF (LWD >= 11.1 .AND. LWD <= 17.0) THEN 
            DVIP = 2
          ELSEIF (LWD >= 7.1 .AND. LWD <= 11.0) THEN  
            DVIP = 1
          ELSEIF (LWD < 7.0) THEN    
            DVIP = 0
          ENDIF
          
        ENDIF 
      
      END SUBROUTINE CALC_DVIP
      
! ------ FUNGICIDE DECICION APPLICATIONS
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
          
          IF (BufferDays > 0) BufferDays = BufferDays - 1
          
          IF (FungActive) THEN
              ResidualDays = ResidualDays - 1
              IF (ResidualDays <= 0) FungActive = .FALSE.
          END IF
          
          CAN_SPRAY = (HEALTH_LAI >= 1.5 )
          
          IF (SUM7 >= 6 .AND. BufferDays == 0 .AND. CAN_SPRAY) THEN
              NSprays = NSprays + 1
              BufferDays = 16 !periodo de carência para nova aplicação (Godoy et al. 2024)
              
              IF (USE_FUNGICIDE) THEN
                  FungActive = .TRUE.
                  ResidualDays = 14 !periodo residual de efeito protetivo do produto (Beruski et al. 2020)
              ELSE
                  FungActive = .FALSE.
                  ResidualDays = 0 
              END IF
              
          END IF
          
      END SUBROUTINE APPLY_FUNGICIDE 
      
!=======================================================================

!***********************************************************************
!     Variable listing (updated 21 May 2025)
!***********************************************************************
! RH = Relative humidity (%)
! LWD = Leaf wetness duration (h)
! T_dew = Dew point temperature (oC)
! Vpsat_o = Saturation Vapor Pressure at Dew Point Temperature (kPa)
! Vpsat_a = Saturation Vapor Pressure at Air Temperature (kPa)
!disease_parameters.txt = input parameters values
!
! T = Daily air temperature (oC)
! TMIN_G = Minimal temperature for spore germination (oC)
! TOT_G = Optial temperature for spore germination (oC)
! TMAX_G = Maximum temperature for spore germination (oC)
! TMIN_D = Minimal temperature for fungal development after germination (oC)
! TOT_D = Optimal temperature for fungal development after germination (oC)
! TMAX_D = Maximum temperature for fungal development after germination (oC)
! FT = Temperature function (0-1)
! FT_D = Temperature function for fungal development after infection (0-1)
! FT_G = Temperature function for substract spores that doesn´t germinate (0-1)
!
! Ymax = Infection efficiency (0-1)
! A = Intrinsic increase rate (0-1)
! B = Intrinsic response to LWD (0-1)
! IR = Infection rate (0-1)
!     
! NDS = Number of dispersed spores (spores/m2)
! LESION_S = Average lesion surface (m2)
! FSS = Fraction of suported spores (0-1)
! DS = Deposited spores (0-1)
!
! LR = Latency rate (0-1)
! LDMin = Minimal days for latency to end (days)
! LS = Latent spores (spores/m2)
! 
! LAF = Lesion age factor (0-1)
! PS = Producted spores (spores/m2)
! PPSR = Rate of potencial spores production (0-1)
!
! LESLIFEMAX = Maximum lesions lifespan (days)
! LA = Lesion age (days)
!
! RVERHULST = Number of individuals per day
! KVERHULST = Maximum number of individuals per LAI unit 
! PREV_IS = Previous infected surface
! INFECTIOUS_S = Infectious surface
! ESP_INOC_SEC_LESION_DAY = Secondary spores producted per day
!
! DAE = Days after emergence
! DISEASE_LAI = Infected and removed LAI (cm2/m2)
!
! DVIP_today = Daily infection-probability score (0–3)
! DVIP_pts = Circular buffer holding the last seven DVIP scores
! idx = Current position inside the circular buffer (1–7)
! SUM7 = Moving sum of DVIP scores over the previous seven days
! BufferDays = Mandatory waiting period (days) between successive sprays
! ResidualDays = Remaining days of preventive protection conferred by the last spray
! NSprays = Cumulative number of fungicide applications in the season
! FungActive = Logical flag; TRUE while preventive protection is active
! USE_FUNGICIDE = Global switch enabling/disabling the chemical control module
! CAN_SPRAY = Logical guard that forbids new sprays when healthy LAI < 1 m² m⁻²
! or when ≤ 15 days remain before the end of the season