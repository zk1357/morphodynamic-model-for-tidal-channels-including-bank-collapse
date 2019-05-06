#include 'define.h'
    
SUBROUTINE HYDRODYNAMIC
    USE TIME_DIN
    USE INTERFACE_DEFINITIONS
#ifdef PRE_HYD        
    USE HYDRO
#endif

#ifdef PRE_AD
    USE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
    USE BANK_EROSION
    USE DEF_EXCHANGE
    USE SED_EXCHANGE
#endif

#ifdef PRE_BANK_COLLAPSE 
    USE BANK_COLLAPSE
#endif

    IMPLICIT NONE
    INTEGER   , PARAMETER                    :: NN0 = SELECTED_REAL_KIND(8)
    REAL(NN0) , DIMENSION(:,:) , POINTER     :: PP1 , PP2
    REAL(NN0) , DIMENSION(:,:) , POINTER     :: QQ1 , QQ2 
    REAL(NN0) , DIMENSION(:,:) , POINTER     :: HH1 , HH2
    REAL(NN0) , DIMENSION(:,:) , POINTER     :: CC1 , CC2
    REAL(NN0) , DIMENSION(:,:) , POINTER     :: DXX , DYY
    REAL(NN0) , DIMENSION(:,:) , POINTER     :: TEMP
    REAL(NN0) , DIMENSION(:)   , ALLOCATABLE :: OUT_PUT
    LOGICAL                                  :: UPSWEEP , DOWNSWEEP
    REAL(NN0)                                :: T0 , TEMP1 , TEMP2 , U , V , UV
    INTEGER                                  :: I1 , I2 , IP , J1 , J2 , JP , I , J , TEMP_N , II , JJ , TEMP_STEP
    INTEGER                                  :: STATUS , CALCU_DAY , BEFORE_DAY , DAY_START , DAY_END , &
                                             &  DAY , HOUR , HOUR_STEP , DAY_STEP , SEC_STEP
    CHARACTER(LEN=80)                        :: ERR_MSG
    CHARACTER*10                             :: FILENAME , XF
    INTEGER                                  :: STEP ,TSTEP
    REAL                                     :: TT1 , TT2 , TTC , BEFORE_STEP
    
    NULLIFY(PP1,PP2,QQ1,QQ2,HH1,HH2,CC1,CC2,TEMP)
    
    CALL CPU_TIME(TT1)

    I1=M_START+1
    I2=M_END  -1
    J1=N_START+1
    J2=N_END  -1
    T0=0.D0
    UPSWEEP = .TRUE.
    DOWNSWEEP = .FALSE.
    IP   = 1
    JP   = 1
    STEP = 0

!POINTER
    PP1=>P1
    PP2=>P2
    QQ1=>Q1
    QQ2=>Q2
    HH1=>H1
    HH2=>H2   
    DXX=>DX1
    DYY=>DY1

    CALCU_DAY  = 2700  !计算天数 模型中一年为360天
    BEFORE_DAY = 1   !提前计算天数，保持流场稳定
    DAY_START  = 0    !水流输出起始天
    DAY_END    = 0    !水流输出起始天
    !DAY        = 1
    !HOUR       = 1
    !HOUR_STEP  = 0
    !DAY_STEP   = 0
    !SEC_STEP   = 0
    TSTEP=INT((CALCU_DAY+BEFORE_DAY)*86400/DT)
    BEFORE_STEP = INT(BEFORE_DAY*86400/DT)
    GLOBTIME%GLOBSTEP%BEFSTEP=BEFORE_STEP
    
LOOP:   DO
            ! time
            CALL TIME_COUNTING
            !write(*,*) GLOBTIME%GLOBSTEP%SECSTEP
            ! flow stability
            !call flow_stability
                
            IF(TSTEP<GLOBTIME%GLOBSTEP%SECSTEP) EXIT LOOP
           
            !X-SWEEP
            
            T0=T0+0.5D0*DT
            
            !UPDATE BOUNDRY CONDITION
            CALL BND_WATERLEVEL(T0,HH2)
            
            !DRYING AND FLOODING PROCESS
            CALL DRYING_AND_FLOODING_X(PP1,HH1,QQ2)
            !call FLOW_TEST_M;PAUSE
            !CONTINUITY AND MOMENTUM EQUATION SOLVING
            CALL CONTINUITY_MOMENTUM_X(PP1,PP2,HH1,HH2,QQ1,QQ2,T0,DXX,DYY,I1,I2,IP,UPSWEEP,DOWNSWEEP)
            
            !call FLOW_TEST_M;PAUSE
                        
            !open(101, file='pp1.txt');open(102, file='pp2.txt');open(103, file='p1.txt');open(104, file='p2.txt')
            !WRITE(101,'(<N_END>e10.2)') ((pp1(I,J),J=N_START, N_END),I=M_START, M_END)
            !WRITE(102,'(<N_END>e10.2)') ((pp2(I,J),J=N_START, N_END),I=M_START, M_END)
            !WRITE(103,'(<N_END>e10.2)') ((p1(I,J),J=N_START, N_END),I=M_START, M_END)
            !WRITE(104,'(<N_END>e10.2)') ((p2(I,J),J=N_START, N_END),I=M_START, M_END)            
            !close(101);close(102);close(103);close(104)
            !pause  
            
            !EXCHANGING,水深每个时间步长交换两次，即h1表示水深；p、q交换一次，双数步长时p1表示流速，否则为p2.
            TEMP=>HH1
            HH1=>HH2
            HH2=>TEMP
            TEMP=>QQ1
            QQ1 =>QQ2
            QQ2 =>TEMP             
          
#ifdef PRE_AD    
            !AD EQUATION SOLVING
            IF ( GLOBTIME%GLOBSTEP%SECSTEP > INT( 0*BEFORE_STEP ) ) CALL TWO_D_SEDIMENTX
#endif 
            
            !Y-SWEEP
            
            T0=T0+0.5D0*DT
            
            !UPDATE BOUNDRY CONDITION
            CALL BND_WATERLEVEL(T0,HH2)  
            
            !DRYING AND FLOODING PROCESS
            CALL DRYING_AND_FLOODING_Y(QQ1,HH1,PP2)
            
            !CONTINUITY AND MOMENTUM EQUATION SOLVING
            CALL CONTINUITY_MOMENTUM_Y(QQ1,QQ2,HH1,HH2,PP1,PP2,T0,DYY,DXX,J1,J2,JP,UPSWEEP,DOWNSWEEP)
            
            !call FLOW_TEST_M;PAUSE
            
            !EXCHANGING
            TEMP=>HH1
            HH1=>HH2
            HH2=>TEMP
            TEMP=>PP1
            PP1 =>PP2
            PP2 =>TEMP            
            
            !open(101, file='pp1.txt');open(102, file='pp2.txt');open(103, file='p1.txt');open(104, file='p2.txt')
            !WRITE(101,'(<N_END>e10.2)') ((pp1(I,J),J=N_START, N_END),I=M_START, M_END)
            !WRITE(102,'(<N_END>e10.2)') ((pp2(I,J),J=N_START, N_END),I=M_START, M_END)
            !WRITE(103,'(<N_END>e10.2)') ((p1(I,J),J=N_START, N_END),I=M_START, M_END)
            !WRITE(104,'(<N_END>e10.2)') ((p2(I,J),J=N_START, N_END),I=M_START, M_END)            
            !close(101);close(102);close(103);close(104)
            !pause 
            
#ifdef PRE_AD     
            !AD EQUATION SOLVING
            IF ( GLOBTIME%GLOBSTEP%SECSTEP > INT( 0*BEFORE_STEP ) ) CALL TWO_D_SEDIMENTY
#endif 
            
#ifdef PRE_MOR 
            !BED ELEVATION UPDATING
            IF ( GLOBTIME%GLOBSTEP%SECSTEP > INT( BEFORE_STEP ) ) CALL BED_ELEVATION_UPDATING
#endif 

            !CYCLING EXCHANGE
            UPSWEEP = .NOT. UPSWEEP
            DOWNSWEEP = .NOT. DOWNSWEEP
            I1 = M_START+M_END-I1
            I2 = M_START+M_END-I2
            IP = - IP
            J1 = N_START+N_END-J1
            J2 = N_START+N_END-J2
            JP = - JP

            ! COUNTING OUTPUT AND COUPLING OTHER PROCESSES 
            ! EVERY TIME STEP
            IF ( GLOBTIME%GLOBSTEP%SECSTEP >  BEFORE_STEP  ) THEN  
                
#ifdef PRE_BANK_EROSION  
                ! CALCULATE LATERAL FLOW EROSION
                CALL FLOW_INCISION
#endif      

                !EVERY HOUR (MOR TIME)               
                IF ( GLOBTIME%GLOBMIN/60==GLOBTIME%GLOBMIN/60.0 .AND. GLOBTIME%GLOBSEC==0 ) THEN                                      
                    
#ifdef PRE_BANK_COLLAPSE      
                    ! SIMULATE BANK COLLAPSE
                    CALL BANK_COLLAPSE_SIMULATION
                    ! NOTE THE EFFECT OF MOR_FACTOR, CHOOSE THE SUITABLE TIME INTERVAL                    
#endif                                                  

                END IF
                
                ! OUTPUT
                CALL OUTPUT(INT(BEFORE_STEP))
                
                ! change MOR_FACTOR
                if ( GLOBTIME%GLOBDAY==2000 ) then
                    MOR_FACTOR=10
                    MOR_FACTOR_BANK=10
                end if
                

            END IF   
            

        END DO LOOP
        
                          
END SUBROUTINE HYDRODYNAMIC