#include 'define.h'
    
SUBROUTINE INITIALIZATION
        USE TIME_DIN 
#ifdef PRE_HYD        
        USE HYDRO
#endif

#ifdef PRE_AD
        USE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
        USE BANK_EROSION
        USE EXCHANGE
        USE DEF_EXCHANGE
        USE SED_EXCHANGE
#endif

#ifdef PRE_BANK_COLLAPSE
        USE BANK_COLLAPSE
#endif

        IMPLICIT NONE    
        INTEGER :: I , J , K , P 
        REAL*8  :: X1 , X2 , X3 , Y1 , Y2 , Y3

        !time count         
        GLOBTIME%GLOBSEC =0
        GLOBTIME%GLOBMIN =0 
        GLOBTIME%GLOBHOUR=0
        GLOBTIME%GLOBHALF=0
        GLOBTIME%GLOBDAY =0
        GLOBTIME%GLOBMON =0
        GLOBTIME%GLOBYEAR=0
        GLOBTIME%GLOBSTEP%SECSTEP=0
        GLOBTIME%GLOBSTEP%HOURSTEP=0
        GLOBTIME%GLOBSTEP%DAYSTEP=0
        GLOBTIME%GLOBSTEP%BEFSTEP=0
        
        ! MOR_FACTOR
        MOR_FACTOR=1
        MOR_FACTOR_BANK=1
        
#ifdef PRE_HYD        
        !MESH SIZE
        DX1(:,:)=20.D0
        DY1(:,:)=5.D0

        !CREEK LOCATION
        I=INT(WIDTH_CREEK_INI/2/DY1(1,1))+1
        TCI_START=-I+INT(M_END+M_START)/2
        TCI_END  = I+INT(M_END+M_START)/2
        TCJ_START= N_END-1
        TCJ_END  = TCJ_START-5*(NUMBER_BV-10)/2+1
        BANK_H   = 1.D0 
#endif
        
#ifdef PRE_BANK_EROSION 
        !ZOOMIN   = BANK_H/2.D0

        !INITIAL BVP
        J=5
        K=(J-1)/2
        DO I=1 , NUMBER_BV/2   !-2详见附图2
            BVP(I,1)=TCI_START+1
            BVP(I,2)=TCJ_START-J*(I-1)-K
            BVP(I+NUMBER_BV/2,1)=TCI_END-1
            BVP(I+NUMBER_BV/2,2)=BVP(I,2)
        END DO
        BVP_IN=BVP
#endif
        
        !GIVING WATER DEPTH
         CALL INI_TOPOGRAPHY(3) !数字含义见子程序  
         TCJ_END  = TCJ_START-5*NUMBER_BV/2+1  !需要计算过渡段的边壁侵蚀，重新赋值
         
        !断电重跑
         !CALL RE_TOPOGRAPHY
         
#ifdef PRE_BANK_EROSION          
         ! 岸壁高度赋值
         DO I=1, NUMBER_BV/2
             BANK_H(I)=ABS( DEP( BVP(I,1)-1, BVP(I,2) )-DEP( BVP(I,1), BVP(I,2) ) )
             BANK_H(I+NUMBER_BV/2)=ABS( DEP( BVP(I+NUMBER_BV/2,1)+1, BVP(I+NUMBER_BV/2,2) )-DEP( BVP(I+NUMBER_BV/2,1), BVP(I+NUMBER_BV/2,2) ) )
         END DO
#endif

        !输出模型参数,matlab作图使用
        CALL PARAMETER_OUT      
        
#ifdef PRE_HYD
        !初始地形作图
        CALL OUTPUT_MOR_M
!PAUSE
        !初始水深流速
        H1(:,:) = DEP(:,:)
        DO I= M_START , M_END                              
            DO J= N_START , N_END
                IF (H1(I,J)<=0.2D0 .AND. H1(I,J)>-90.D0) H1(I,J)=0.2D0
            END DO
        END DO 
        H2(:,:)=H1(:,:)
        P1(:,:) = 0.D0
        P2(:,:) = 0.D0
        Q1(:,:) = 0.D0
        Q2(:,:) = 0.D0
        QA_Q(:,:) = 0.D0
        QA_A(:,:) = 0.D0
        QA_Q_P = 0.D0
        QA_A_P = 0.D0       
#endif

#ifdef PRE_HYD
        !NDD READY
        DO I= M_START , M_END                              
            DO J= N_START , N_END                          
               NDD(I,J) = 0
               IF(I==M_END)        NDD(I,J) = -1
               IF(I==M_START)      NDD(I,J) = -1  
               IF(J==N_START)      NDD(I,J) =  2               
               IF(J==N_END)        NDD(I,J) = -1               
               IF(DEP(I,J)<-90.D0) NDD(I,J) = -1             
            END DO                                         
        END DO  
 !       open(1,file='ndd.txt')
 !write(1,'(101i4)')((ndd(i,j),i=1,101),j=1,401)
 !close(1)
 !pause
        !VDD READY
        VDDX=-1 ; VDDY=-1 
     
#endif

#ifdef PRE_AD
        C1(:,:) = 0.D0
        C2(:,:) = 0.D0
#endif

#ifdef PRE_BANK_EROSION
        !JUDGE_BANK READY
        JUDGE_BANK=1
        JUDGE_TIME=-1
        BF=0
        
        !INITIAL LOW WATER LEVEL
        LHL =0.D0
        LHLB=0.D0
        STATION_BANK=0.D0

        !SED_EXCHANGE READY
        S_BANK_EROSION=0.D0
        S_BANK_COLLAPSE=0.D0
        
        !INITIAL BANK  READY  
        JUDGE_TEMP      =0
        NUMBER_LOAD     =1
        NUMBER_BANKCYCLE=1
        LENGTH_COLLAPSE =0.D0
        LENGTH_FLOW     =0.D0
        LENGTH_RETREATS =0.D0
        LENGTH_EROSION_DAY=0.D0
        LENGTH_COLLAPSE_DAY=0.D0
        RETREATS_TOTAL  =0.D0
        UPDATE_C        =0.D0
        BVP_VELOCITY    =0.D0
        QE_OUT          =0.D0
        LENGTH_COLLAPSE_BEFORE=0.D0
        BANK_START_VELOCITY  =SQRT( TE_BANK/FLOWDENSITY_BANK/CD_BANK )
        
        !EXCHANGE READY
        EXCHANGE1%GRID_BANK_H = 2.D0
        EXCHANGE1%GRID_BANK_W = 10.D0
        EXCHANGE1%ZOOMIN = BANK_H/EXCHANGE1%GRID_BANK_H
        EXCHANGE1%DIG_NUMBER = 0
        POSITION=0
        DO I =1, NUMBER_BV
            POSITION(I,:)=EXCHANGE1%GRID_BANK_W*EXCHANGE1%ZOOMIN(I)
        END DO
        EXCHANGE1%JUDGE_BANK = 1
        EXCHANGE1%FLOWEROSION_CHECK = 0
        EXCHANGE1%NUMBER_LOAD=0
        EXCHANGE1%NUMBER_BANKCYCLE=1
        EXCHANGE1%FLOW_EXCESS_JUDGE=0
        EXCHANGE1%COLLAPSE_TIME=0
        EXCHANGE1%COLLAPSE_DEF=0
        EXCHANGE1%COLLAPSE_VOLUME=0
        EXCHANGE1%FLOWEROSION_VOLUME=0
        EXCHANGE1%BANKEROSION_VOLUME=0
        EXCHANGE1%SHIELD_VOLUME=0
#endif

#ifdef PRE_BANK_COLLAPSE

        ! READ GRID AND MESH INFORMATION
        OPEN (UNIT=41, FILE='初始网格节点坐标.TXT',  STATUS='OLD')
        OPEN (UNIT=42, FILE='初始网格单元节点号.TXT',STATUS='OLD')   
        READ (42,*) ((EXCHANGE1%MESH(K,J),J=1,3),K=1,EXCHANGE1%ELEMENT)    
        READ (41,*) ((EXCHANGE1%GRID(I,J),J=1,2),I=1,EXCHANGE1%NODE)
        CLOSE(41)
        CLOSE(42)
        
        ! READ GETECHNICAL PARAMETER
        OPEN (UNIT=40, FILE='土体参数.TXT',STATUS='OLD')
        READ (40,*) BANKCOLLAPSE%MU , BANKCOLLAPSE%COH , BANKCOLLAPSE%FRI
        CLOSE(40)
        
        ! SET INITIAL VALUE
        BANKCOLLAPSE%T_LE            = 1.D0
        BANKCOLLAPSE%MODULUS_INITIAL = 5.D6
        BANKCOLLAPSE%MODULUS         = BANKCOLLAPSE%MODULUS_INITIAL
        BANKCOLLAPSE%TSTRESS         = 0.D0
        BANKCOLLAPSE%STATION         = 0
        BANKCOLLAPSE%F               = 0.D0
        BANKCOLLAPSE%WY              = 0.D0
        BANKCOLLAPSE%RESTRICTION     = 1
        BANKCOLLAPSE%MATCH           = 0
        BANKCOLLAPSE%COLL_TO_BED     = 0.5
        
        ! FIND MESH AROUND NODE
        DO I = 1 , EXCHANGE1%NODE
            K=1
            DO J = 1 , EXCHANGE1%ELEMENT
                IF (I==EXCHANGE1%MESH(J,1) .OR. &
                &   I==EXCHANGE1%MESH(J,2) .OR. &
                &   I==EXCHANGE1%MESH(J,3)) THEN
                    BANKCOLLAPSE%MATCH(I,K)=J    
                    K=K+1
                END IF
            END DO
        END DO
        
        ! SET BOUNDARY CONDITION FOR BANK
        DO I = 1 , EXCHANGE1%NODE
            IF (EXCHANGE1%GRID(I,1)==0 ) THEN
                BANKCOLLAPSE%RESTRICTION(2*I-1)=0
!                RESTRICTION(2*I  )=0
            END IF
            IF (EXCHANGE1%GRID(I,2)==0) THEN
                BANKCOLLAPSE%RESTRICTION(2*I-1)=0
                BANKCOLLAPSE%RESTRICTION(2*I  )=0
            END IF
        END DO 
        
        ! SET INITIAL BANK SIZE
        CALL MESH_DIVISION_INITIAL
#endif

    END
    
    
    SUBROUTINE INI_TOPOGRAPHY(CHA)

#ifdef PRE_HYD        
        USE HYDRO
#endif
        
        INTEGER :: CHA ! 1:均一水深; 2:均一坡度; 3:均一坡度+初始中间潮沟; 4:Roelvink2006算例地形; 5~7:Maanen2013算例地形; 8:Roelvink2006算例地形修改;
        REAL*8  :: X , Y
        
        DEP(:      ,    :)= DEPOUTSEA                      
        DEP(:      ,N_END)=-99.99D0                       
        DEP(M_START,    :)=-99.99D0                       
        DEP(M_END  ,    :)=-99.99D0   
        
        IF ( CHA==2 ) THEN
        !Giving slope,closed boundary
            DO I= M_START+1 , M_END-1                             
                DO J= N_START+1 , N_END-1  
                    call random_number(X)
                    call random_number(Y)
                    IF (Y>0.5) X=-X                     
                    IF (J<41) THEN
                        DEP(I,J)= DEPOUTSEA
                    ELSE
                        DEP(I,J)= DEPOUTSEA-SLOPE*DX1(I,J)*(J-1)+X*0.2      
                    END IF
                END DO                                            
            END DO  
            
        ELSE IF ( CHA==3 ) THEN            
!#ifdef PRE_BANK_EROSION 
            !在滩面中间加BANK_H深的潮沟,详见附图2
            DO I= M_START , M_END                                     
                DO J= N_START , N_END
                    IF (J<=75) THEN
                        DEP(I,J)= DEPOUTSEA-SLOPE*DX1(I,J)*(J-1)
                    ELSE
                        DEP(I,J)=DEPOUTSEA-SLOPE*DX1(I,J)*75-SLOPE1*DX1(I,J)*(J-75)
                    END IF
                    IF (I>TCI_START .AND. I<TCI_END .AND. &
                    &   J<=TCJ_START .AND. J> TCJ_END) DEP(I,J)=DEP(I,J)+BANK_H(1) 
                     IF (I==TCI_START .AND. &
                    &   J<=TCJ_START .AND. J> TCJ_END) DEP(I,J)=DEP(I,J)                    
                     IF (I==TCI_END .AND. &
                    &   J<=TCJ_START .AND. J> TCJ_END) DEP(I,J)=DEP(I,J)                     
                END DO                                            
            END DO   
        
            !加坡度潮沟
            DO I= TCI_START+1 , TCI_END-1                                     
                DO J= TCJ_END , TCJ_END-25 , -1                              
                    DEP(I,J)=DEP(I,J)+0.002*DX1(1,1)*(J-TCJ_END+25)    
                END DO 
            END DO   
            
            !坡度潮沟扩大
            !P=1
            !DO J=TCJ_END-1 , TCJ_END-50 , -5
            !    DO K=0, 4
            !        DO I=TCI_START+1-P , TCI_END-1+P 
            !            DEP(I,J-K)=1
            !        END DO
            !    END DO
            !    P=P+2
            !END DO
            
            !!
            !!改初始潮沟深度
            !DO I = 51 , 60                                       
            !    BANK_H(I)=BANK_H(1)-(I-51)*0.1
            !    BANK_H(I+NUMBER_BV/2)=BANK_H(1)-(I-51)*0.1
            !END DO
!#endif
        ELSE IF ( CHA==31 ) THEN            
!#ifdef PRE_BANK_EROSION 
            !在滩面中间加BANK_H深的潮沟,详见附图2
            DO I= M_START+1 , M_END-1                                     
                DO J= N_START+1 , N_END-1
                    IF (J<=100) THEN
                        DEP(I,J)= DEPOUTSEA
                    ELSE
                        DEP(I,J)=-1.D0
                    END IF
                    IF (I>TCI_START .AND. I<TCI_END .AND. &
                    &   J<=TCJ_START .AND. J> 100) DEP(I,J)=DEP(I,J)+BANK_H(1) 
                END DO                                            
            END DO   
        
!#endif            
        ELSE IF ( CHA==4 ) THEN
            call random_seed()
            DO I= M_START+1 , M_END-1                             
                DO J= N_START+1 , N_END-1  

                    call random_number(X)
                    call random_number(Y)
                    IF (Y>0.5) X=-X                   
                    IF (J<=51) THEN
                        DEP(I,J) = DEPOUTSEA-SLOPE*DX1(I,J)*(J-1) 
                    ELSE
                        DEP(I,J) = 2.D0+X*0.1
                    END IF
                    IF ( I>=(M_START+M_END)/2+20 .OR. I<=(M_START+M_END)/2-20 ) THEN
                        IF ( J==52 ) THEN
                            DEP(I,J) = -99.9D0
                        END IF
                    END IF
                END DO                                            
            END DO 
 
        ELSE IF ( CHA==8 ) THEN
            call random_seed()
            DO I= M_START+1 , M_END-1                             
                DO J= N_START+1 , N_END-1  

                    call random_number(X)
                    call random_number(Y)
                    IF (Y>0.5) X=-X                   
                    IF (J<=51) THEN
                        DEP(I,J) = DEPOUTSEA-SLOPE*DX1(I,J)*(J-1) 
                    ELSE
                        DEP(I,J) = 2.D0+X*0.1
                    END IF
                    IF ( I>=(M_START+M_END)/2+20 .OR. I<=(M_START+M_END)/2-20 ) THEN
                        IF ( J<=71 .AND. J>=51 ) THEN
                            DEP(I,J) = -99.9D0
                        END IF
                    END IF
                END DO                                            
            END DO            
            
        ELSE IF ( CHA==5 ) THEN
            call random_seed()
            DO I= M_START+1 , M_END-1                             
                DO J= N_START+1 , N_END-1  
                    call random_number(X)
                    call random_number(Y)
                    IF (Y>0.5)    X=-X                   
                    IF (J<=40)  THEN
                        DEP(I,J) = 8.D0+X*0.1
                    ELSE IF (J<=61) THEN 
                        DEP(I,J) = 8-0.003*DX1(I,J)*(J-41)+X*0.1 
                    ELSE IF (J<=91) THEN
                        DEP(I,J) = 2.D0+X*0.1
                    ELSE
                        DEP(I,J) = 2-0.0005*DX1(I,J)*(J-91)+X*0.1
                    END IF
                    
                    IF ( I>=(M_START+M_END)/2+25 .OR. I<=(M_START+M_END)/2-25 ) THEN
                        IF ( J<=81 .AND. J>=71 ) THEN
                            DEP(I,J) = -99.9D0
                        END IF
                    END IF
                    
                END DO                                            
            END DO   
            
        END IF
        
    END SUBROUTINE 
    
    SUBROUTINE RE_TOPOGRAPHY
    
        USE TIME_DIN 
#ifdef PRE_HYD        
        USE HYDRO
#endif

#ifdef PRE_AD
        USE Morphodynamics
#endif

#ifdef PRE_BANK_EROSION
        USE BANK_EROSION
        USE EXCHANGE
        USE SED_EXCHANGE
#endif
    
        IMPLICIT NONE
    
        INTEGER            :: I , J , K , STEP
        CHARACTER*40       :: FILENAME
    
        FILENAME='001000MORM.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME
    
        OPEN(600, FILE= FILENAME, FORM='UNFORMATTED', ACCESS='DIRECT' , RECL=2)
        K=1
        DO I = M_START+1 , M_END-1
            DO J = N_START+1 , N_END-1
                READ(600 , REC = K ) DEP(I,J)
                K=K+1
            END DO
        END DO
    
        CLOSE(600)
    
    
    END SUBROUTINE