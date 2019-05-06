#include 'define.h'    
!输出模型参数,matlab作图使用,文件号从600开始
SUBROUTINE PARAMETER_OUT

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
    REAL*8 , DIMENSION(:,:) :: A(12)
    INTEGER                 :: I
    
    OPEN(600 , FILE= 'MATLAB_PLOT\HYDRO_PARAMETERS_REAL.RES'    , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2) !输出水流模块参数
    OPEN(601 , FILE= 'MATLAB_PLOT\HYDRO_PARAMETERS_INTEGER.RES' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2) !输出水流模块参数
    
    !需要输出的模型参数包括：振幅、边界水深、底床启动应力、沉积应力、沉速、底床冲刷系数、边界含沙量、岸壁启动应力、岸壁冲刷系数、BVP、I J节点数   
    A=0.D0

#ifdef PRE_HYD        
    A(1)=AMP
    A(2)=DEPOUTSEA
#endif
    
#ifdef PRE_AD
    A(3)=TE
    A(4)=TD
    A(5)=WS
    A(6)=ME
    A(7)=C_BOUNDARY    
#endif

#ifdef PRE_BANK_EROSION
    A(8)=TE_BANK    
    A(9)=ME_BANK    
#endif

    WRITE(600 , REC= 1 )  A(1)    
    WRITE(600 , REC= 2 )  A(2) 
    WRITE(600 , REC= 3 )  A(3) 
    WRITE(600 , REC= 4 )  A(4)
    WRITE(600 , REC= 5 )  A(5) 
    WRITE(600 , REC= 6 )  A(6)
    WRITE(600 , REC= 7 )  A(7)
    WRITE(600 , REC= 8 )  A(8)
    WRITE(600 , REC= 9 )  A(9)
#ifdef PRE_BANK_EROSION
    WRITE(601 , REC= 3 )  NUMBER_BV/2
#endif    
    WRITE(601 , REC= 1 )  M_END
    WRITE(601 , REC= 2 )  N_END
    
    CLOSE(600)
    CLOSE(601)

#ifdef PRE_BANK_EROSION
!输出水流模块参数
    OPEN(600 , FILE= 'MATLAB_PLOT\BVP_PARAMETERS.RES'    , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2) 
    DO I = 1 , NUMBER_BV/2
        WRITE(600 , REC=I) BVP(I,2)
    END DO
    CLOSE(600)
#endif    
    
END SUBROUTINE

!计算Q~A关系    
SUBROUTINE CACULATE_Q_A(STEP)

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
    
    INTEGER            :: I , J , STEP
    REAL*8             :: TEMPQ  
    
    DO J = N_START+1 , N_END-1
        TEMPQ=0.D0
        DO I = M_START+1 , M_END-1
        !DO I = TCI_START+1 , TCI_END-1    
            !IF (H1(I,J)<=0.05) CYCLE
            TEMPQ = TEMPQ + P1(I,J)*DY1(1,1)
            QA_A(STEP,J) = QA_A(STEP,J) + H1(I,J)*DT*DY1(1,1)/43200.
        END DO
        IF ( QA_Q(STEP,J)<TEMPQ ) QA_Q(STEP,J) = TEMPQ
    END DO
    
    
    END SUBROUTINE
 
SUBROUTINE CACULATE_Q_A_POINT(STEP,I)

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
    
    INTEGER            :: I , J , STEP , K
    REAL*8             :: TEMPQ  
    
    IF (I==51) K=1
    IF (I==25) K=2
    TEMPQ=0.D0
    DO J = N_START , N_END
        TEMPQ = P1(I,J)*DY1(1,1)
        IF ( QA_Q_P(STEP,K,J)<TEMPQ ) QA_Q_P(STEP,K,J) = TEMPQ
        QA_A_P(STEP,K,J) = QA_A_P(STEP,K,J) + H1(I,J)*DT*DY1(1,1)/43200.
    END DO

    
    
END SUBROUTINE
       
SUBROUTINE OUTPUT_Q_A

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
    
    INTEGER            :: I , J , K , P
    CHARACTER*40       :: FILENAME , FILENAME1

    FILENAME='QA.RES'
    FILENAME='MATLAB_PLOT\'//FILENAME
    FILENAME1='QA_P.RES'
    FILENAME1='MATLAB_PLOT\'//FILENAME1
    
    OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    OPEN(601 , FILE= FILENAME1 , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2) 
    
    K=1
    DO I = 1 , 30
        DO J = N_START , N_END-1
            WRITE( 600 , REC = K ) QA_Q(I,J)
            K=K+1
        END DO
    END DO
    DO I = 1 , 30
        DO J = N_START , N_END-1
            WRITE( 600 , REC = K ) QA_A(I,J)
            K=K+1
        END DO
    END DO   
    
    K=1
    DO P = 1 , 2
        DO I = 1 , 30
            DO J = N_START , N_END-1
                WRITE( 601 , REC = K ) QA_Q_P(I,P,J)
                K=K+1
            END DO
        END DO
    END DO
    DO P = 1 , 2
        DO I = 1 , 30
            DO J = N_START , N_END-1
                WRITE( 601 , REC = K ) QA_A_P(I,P,J)
                K=K+1
            END DO
        END DO    
    END DO
    
    
    CLOSE(600)
    CLOSE(601)
    
END SUBROUTINE   
           
!MATLAB单点、剖面水动力信息作图
SUBROUTINE OUTPUT_HYD_M(NODE_I,NODE_J,STEP,POINT_OR_PROFILE)
    
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
    
    INTEGER               :: I , J , K , NODE_I , NODE_J , STEP , POINT_OR_PROFILE
    REAL*8                :: U , V , UV
    REAL*8 , DIMENSION(:) :: U_F(N_END) , V_F(N_END) , UV_F(N_END)
    CHARACTER*40          :: FILENAME
    
    ! POINT_OR_PROFILE, 0单点、 1 NODE_I剖面、 2 NODE_J剖面
    IF ( POINT_OR_PROFILE==0 ) THEN
        
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') NODE_I
        WRITE(FILENAME(5:8) , '(I4.4)') NODE_J
        FILENAME(9:16)='HYD0.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME
        IF ( STEP==1 ) THEN
            OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*5)
        ELSE
            OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*5)
        END IF  
        
        CALL UV_CACULATION( NODE_I , NODE_J , U , V , UV )
        
        WRITE( 600 , REC = STEP ) H1(NODE_I,NODE_J) , P1(NODE_I,NODE_J) , Q1(NODE_I,NODE_J) , U , V
        
    ELSE IF ( POINT_OR_PROFILE==1 ) THEN
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') NODE_I
        WRITE(FILENAME(5:8) , '(I4.4)') NODE_J
        FILENAME(9:16)='HYD1.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME       
        IF ( STEP==1 ) THEN
            OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=10*N_END)
        ELSE
            OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=10*N_END)
        END IF 
        
        DO J = N_START , N_END
            CALL UV_CACULATION( NODE_I , J , U_F(J) , V_F(J) , UV_F(J) )
        END DO
        
        WRITE( 600 , REC = STEP ) (H1(NODE_I,J),J=N_START,N_END) , (P1(NODE_I,J),J=N_START,N_END) , (Q1(NODE_I,J),J=N_START,N_END) , &
            &                     (U_F(J),J=N_START,N_END) , (V_F(J),J=N_START,N_END)
              
    ELSE IF ( POINT_OR_PROFILE==2 ) THEN
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') NODE_I
        WRITE(FILENAME(5:8) , '(I4.4)') NODE_J
        FILENAME(9:16)='HYD2.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME       
        IF ( STEP==1 ) THEN
            OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=10*N_END)
        ELSE
            OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=10*N_END)
        END IF 
        
        DO I = N_START , N_END
            CALL UV_CACULATION( I , NODE_J , U_F(I) , V_F(I) , UV_F(I) )
        END DO
        
        WRITE( 600 , REC = STEP ) (H1(I,NODE_J),I=M_START,M_END) , (P1(I,NODE_J),I=M_START,M_END) , (Q1(I,NODE_J),I=M_START,M_END) , &
            &                     (U_F(I),I=M_START,M_END) , (V_F(I),I=M_START,M_END)        
        
    END IF
    
    CLOSE(600)
    
    END SUBROUTINE OUTPUT_HYD_M

#ifdef PRE_BANK_EROSION
!MATLAB BVP点水流信息作图
SUBROUTINE OUTPUT_HYD_BVP_M(STEP)
    
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
    
    INTEGER            :: I , J , K , NODE_I , NODE_J , STEP
    REAL*8             :: U , V , UV
    CHARACTER*40       :: FILENAME
    
    FILENAME='BVP_VELOCITY.RES'
    FILENAME='MATLAB_PLOT\'//FILENAME
    
    IF ( STEP==1 ) THEN
        OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    ELSE
        OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    END IF    
    
    K=1
    
    DO I = 1 , NUMBER_BV/2
        U=2.D0*P1( BVP(I,1),BVP(I,2) )/( H1( BVP(I,1),BVP(I,2)+1 )+H1( BVP(I,1),BVP(I,2) ))
        WRITE( 600 , REC = (STEP-1)*NUMBER_BV/2+K ) U
        K=K+1
    END DO
    
    CLOSE(600)
    
    END SUBROUTINE OUTPUT_HYD_BVP_M

!输出BVP处的侧向侵蚀、底床侵蚀潮周期累积量，岸壁高度   
 SUBROUTINE OUTPUT_QE
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
    
    INTEGER            :: I , J , K , P, STEP
    CHARACTER*40       :: FILENAME , FILENAME1

    FILENAME='QE.RES'
    FILENAME='MATLAB_PLOT\'//FILENAME
    FILENAME1='BANK_H.RES'
    FILENAME1='MATLAB_PLOT\'//FILENAME1
    
    STEP=(GLOBTIME%GLOBSTEP%DAYSTEP-1)/10+1
    
    IF ( STEP==1 ) THEN
        OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        OPEN(601 , FILE= FILENAME1 , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    ELSE
        OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        OPEN(601 , FILE= FILENAME1 , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    END IF  
    
    K=3*NUMBER_BV/2*(STEP-1)+1
    DO I = 1 , 3*NUMBER_BV/2
        WRITE( 600 , REC = K ) QE_OUT(I,STEP)
        K=K+1
    END DO
    
    K=NUMBER_BV/2*(STEP-1)+1
    DO I = 1 , NUMBER_BV/2
        WRITE( 601 , REC = K ) BANK_H(I)
        K=K+1
    END DO    
    
    
    CLOSE(600)
    CLOSE(601)
    
END SUBROUTINE     
#endif   
    
#ifdef PRE_AD
!MATLAB单点、剖面含沙量信息作图
SUBROUTINE OUTPUT_SED_M(NODE_I,NODE_J,STEP,POINT_OR_PROFILE)
    
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
    
    INTEGER               :: I , J , K , NODE_I , NODE_J , STEP , &
        &                    POINT_OR_PROFILE, BVPN
    REAL*8                :: CHEZY, U, V, UV, TB
    CHARACTER*40          :: FILENAME
    
    ! POINT_OR_PROFILE, 0单点、 1 NODE_I剖面、 2 NODE_J剖面
    IF ( POINT_OR_PROFILE==0 ) THEN
        
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') NODE_I
        WRITE(FILENAME(5:8) , '(I4.4)') NODE_J
        FILENAME(9:16)='SED0.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME
        IF ( STEP==1 ) THEN
            OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        ELSE
            OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        END IF 
        
        WRITE( 600 , REC = STEP ) C1(NODE_I,NODE_J)
                 
        
    ELSE IF ( POINT_OR_PROFILE==1 ) THEN
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') NODE_I
        WRITE(FILENAME(5:8) , '(I4.4)') NODE_J
        FILENAME(9:16)='SED1.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME       
        IF ( STEP==1 ) THEN
            OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*N_END)
        ELSE
            OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2*N_END)
        END IF 
        
        WRITE( 600 , REC = STEP ) (C1(NODE_I,J),J=N_START,N_END)
              
    ELSE IF ( POINT_OR_PROFILE==2 ) THEN
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') NODE_I
        WRITE(FILENAME(5:8) , '(I4.4)') NODE_J
        FILENAME(9:16)='SED2.RES'
        FILENAME='MATLAB_PLOT\'//FILENAME       
        IF ( STEP==1 ) THEN
            OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=10*N_END)
        ELSE
            OPEN(600 , FILE= FILENAME , STATUS='OLD'     , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=10*N_END)
        END IF 
        
        WRITE( 600 , REC = STEP ) (C1(I,NODE_J),I=M_START,M_END)      
        
    END IF
    
    CLOSE(600)
    
    END SUBROUTINE OUTPUT_SED_M

#endif    

!MATLAB地貌演变作图
SUBROUTINE OUTPUT_MOR_M
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
    
    FILENAME=''
    WRITE(FILENAME(1:6) , '(I6.6)') GLOBTIME%GLOBSTEP%DAYSTEP
    FILENAME(7:14)='MORM.RES'
    FILENAME='MATLAB_PLOT\'//FILENAME
    
    OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    K=1
    DO I = M_START+1 , M_END-1
        DO J = N_START+1 , N_END-1 ,1
            WRITE(600 , REC = K ) DEP(I,J)
            K=K+1
        END DO
    END DO
    
    CLOSE(600)


    END SUBROUTINE OUTPUT_MOR_M  

#ifdef PRE_BANK_EROSION

!!潮沟拓宽演变过程作图
!SUBROUTINE OUTPUT_TC_WIDEN_M
!    
!#ifdef PRE_HYD        
!    USE HYDRO
!#endif
!
!#ifdef PRE_AD
!    USE Morphodynamics
!#endif
!
!#ifdef PRE_BANK_EROSION
!    USE BANK_EROSION
!    USE EXCHANGE
!    USE SED_EXCHANGE
!#endif
!    
!    IMPLICIT NONE
!    INTEGER            :: I , J , K , STEP
!    CHARACTER*40       :: FILENAME , FILENAME1
!    
!    FILENAME=''
!    FILENAME(1:8)='TCWR.RES'
!    FILENAME='MATLAB_PLOT\'//FILENAME    
!    
!    FILENAME1=''
!    FILENAME1(1:8)='TCWL.RES'
!    FILENAME1='MATLAB_PLOT\'//FILENAME1    
!    
!    IF (STEP==1) THEN
!        OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!        OPEN(601 , FILE= FILENAME1 , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!    ELSE
!        OPEN(600 , FILE= FILENAME , STATUS='OLD' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!        OPEN(601 , FILE= FILENAME1 , STATUS='OLD' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!    END IF
!    
!#ifdef PRE_BANK_COLLAPSE
!    DO I = 1 , NUMBER_BV/2
!        J=I+NUMBER_BV/2
!        WRITE(601 , REC= I+NUMBER_BV/2*(STEP-1)) -0.5D0*WIDTH_CREEK_INI+DY1(1,1)*(M_END-M_START)*0.5D0-RETREATS_TOTAL( I )
!        WRITE(600 , REC= I+NUMBER_BV/2*(STEP-1))  0.5D0*WIDTH_CREEK_INI+DY1(1,1)*(M_END-M_START)*0.5D0+RETREATS_TOTAL( I)
!    END DO
!    GOTO 100
!#endif 
!
!    DO I = 1 , NUMBER_BV/2
!        J=I+NUMBER_BV/2
!        WRITE(601 , REC= I+NUMBER_BV/2*(STEP-1)) -0.5D0*WIDTH_CREEK_INI+DY1(1,1)*(M_END-M_START)*0.5D0-LENGTH_FLOW( I , NUMBER_BANKCYCLE(I) )
!        WRITE(600 , REC= I+NUMBER_BV/2*(STEP-1))  0.5D0*WIDTH_CREEK_INI+DY1(1,1)*(M_END-M_START)*0.5D0+LENGTH_FLOW( J , NUMBER_BANKCYCLE(J) )
!    END DO
!    
!100 CONTINUE
!    
!    CLOSE(600)
!    CLOSE(601)
!    
!    END SUBROUTINE OUTPUT_TC_WIDEN_M

!潮沟断面演变过程作图
!SUBROUTINE OUTPUT_CS_M(II,JJ,STEP)
!#ifdef PRE_HYD        
!    USE HYDRO
!#endif
!
!#ifdef PRE_AD
!    USE Morphodynamics
!#endif
!
!#ifdef PRE_BANK_EROSION
!    USE BANK_EROSION
!    USE EXCHANGE
!    USE SED_EXCHANGE
!#endif
!    
!    IMPLICIT NONE
!    INTEGER            :: I , J , K , II , JJ , STEP
!    CHARACTER*40       :: FILENAME    
!
!    FILENAME=''
!    WRITE(FILENAME(1:4) , '(I4.4)') II
!    WRITE(FILENAME(5:8) , '(I4.4)') JJ
!    FILENAME(9:14)='CS.RES'
!    FILENAME='MATLAB_PLOT\'//FILENAME
!    
!    IF (STEP==1) THEN
!        OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!    ELSE
!        OPEN(600 , FILE= FILENAME , STATUS='OLD' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!    END IF
!    
!    K=1
!    
!    IF (II==0) THEN
!        DO I = M_START+1 , M_END-1
!            WRITE(600 , REC= K+(M_END-1-M_START)*(STEP-1) ) DEP(I,JJ)
!            k=k+1
!        END DO
!    ELSE
!        DO I = N_START+1 , N_END-1
!            WRITE(600 , REC= K+(N_END-1-N_START)*(STEP-1) ) DEP(II,J)
!            k=k+1
!        END DO
!    END IF
!         
!    CLOSE(600)
!    
!    END SUBROUTINE OUTPUT_CS_M


!输出每个断面每天COLLAPSE EROSION    
SUBROUTINE OUTPUT_CONTRIBUTION_M
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
    CHARACTER*40       :: FILENAME , FILENAME1
    
    FILENAME=''
    FILENAME(1:8)='EROS.RES'
    FILENAME='MATLAB_PLOT\'//FILENAME  
    
    FILENAME1=''
    FILENAME1(1:8)='COLL.RES'
    FILENAME1='MATLAB_PLOT\'//FILENAME1         
    
    STEP=GLOBTIME%GLOBYEAR*360+GLOBTIME%GLOBMON*30+GLOBTIME%GLOBDAY
    IF (STEP==1) THEN
        OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        OPEN(601 , FILE= FILENAME1 , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    ELSE
        OPEN(600 , FILE= FILENAME , STATUS='OLD' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
        OPEN(601 , FILE= FILENAME1 , STATUS='OLD' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
    END IF
    
    DO I = 1 , NUMBER_BV
        
        J=I+NUMBER_BV*(STEP-1)
        
        WRITE( 600 , REC= J ) LENGTH_EROSION_DAY(I) 
        WRITE( 601 , REC= J ) LENGTH_COLLAPSE_DAY(I) 
        LENGTH_EROSION_DAY(I) =0.D0
        LENGTH_COLLAPSE_DAY(I)=0.D0
        
    END DO

    CLOSE(600)
    CLOSE(601)
    
END SUBROUTINE OUTPUT_CONTRIBUTION_M
    
#endif    

#ifdef PRE_BANK_COLLAPSE

!!每块岸壁COLLAPSE EROSION所占比例
!SUBROUTINE OUTPUT_COLLAPSE_M
!    
!#ifdef PRE_HYD        
!    USE HYDRO
!#endif
!
!#ifdef PRE_AD
!    USE Morphodynamics
!#endif
!
!#ifdef PRE_BANK_EROSION
!    USE BANK_EROSION
!    USE EXCHANGE
!    USE SED_EXCHANGE
!#endif
!    
!    IMPLICIT NONE
!    INTEGER            :: I , J , K , STEP
!    CHARACTER*40       :: FILENAME
!    
!    FILENAME=''
!    FILENAME(1:12)='COLLAPSE.RES'
!    FILENAME='MATLAB_PLOT\'//FILENAME     
!    
!    !读取顺序 几块岸壁--总侵蚀距离---水流、坍塌侵蚀距离
!    
!    OPEN(600 , FILE= FILENAME , STATUS='REPLACE' , FORM='UNFORMATTED' , ACCESS='DIRECT' , RECL=2)
!    K=1
!    DO I = 1 , NUMBER_BV/2
!        WRITE(600 , REC= K) NUMBER_BANKCYCLE(I)
!        K=K+1       
!    END DO
!
!    DO I = 1 , NUMBER_BV/2
!        WRITE(600 , REC= K) RETREATS_TOTAL(I)
!        K=K+1        
!    END DO
!    
!    DO I = 1 , NUMBER_BV/2
!        DO J = 1 , NUMBER_BANKCYCLE(I)
!            WRITE(600 , REC= K) LENGTH_FLOW(I,J)
!            K=K+1
!            WRITE(600 , REC= K) LENGTH_COLLAPSE(I,J)
!            K=K+1           
!        END DO
!    END DO
!    
!    CLOSE(600)
!    
!    END SUBROUTINE OUTPUT_COLLAPSE_M

#endif