#include 'define.h'  
#ifdef PRE_BANK_EROSION  

SUBROUTINE MESH_DIVISION
            
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
    
    IMPLICIT NONE


    INTEGER   , PARAMETER                    :: NFI=SELECTED_REAL_KIND(8)
    REAL(NFI) , ALLOCATABLE , DIMENSION(:,:) :: GRID , MESH !记录网格节点信息
    REAL(NFI) , ALLOCATABLE , DIMENSION(:,:) :: MEA
    REAL(NFI) , ALLOCATABLE , DIMENSION(:)   :: A_X , A_Y 
    INTEGER   , ALLOCATABLE , DIMENSION(:,:) :: LINE
    INTEGER   , ALLOCATABLE , DIMENSION(:)   :: DIG_NUMBER
    INTEGER                                  :: NODE , ELEMENT ,TEMP
    INTEGER                                  :: I ,J, K , NUMBER1 , NUMBER2 , NUMBER3
    REAL(NFI)                                :: X1 , X2 , X3 , X4 , Y1 , Y2 , Y3 , Y4 ,&
    &                                           TEMP_X
    CHARACTER*40                             :: FILENAME , FILENAME_PATH
    
!----------------------------------------------------------------------
    ALLOCATE(      A_X  (EXCHANGE1%NODE) )  !记录侵蚀线外节点坐标
    ALLOCATE(      A_Y  (EXCHANGE1%NODE) )
    ALLOCATE(     LINE(100,2) )  !记录侵蚀线各个线段
    
LOOP0:  DO I = 1 , NUMBER_BV
    
            !防止潮间带下部发生岸壁不够用,超出给定岸壁范围
            IF ( EXCHANGE1%FLOWEROSION_VOLUME( I, EXCHANGE1%NUMBER_BANKCYCLE(I) )/DX1(1,1)/BANK_H(I) >= BANK_H(I) ) THEN
                EXCHANGE1%FLOW_EXCESS_JUDGE(I) = 1
                CYCLE
            END IF
            
            IF ( EXCHANGE1%FLOWEROSION_CHECK(I)==1 ) CYCLE        
            
            IF (BANK_H(I)<=BANK_EROSION_LIMIT) CYCLE
            
            !IF (I==59) THEN
            !    WRITE(*,*) BANK_H(I), EXCHANGE1%FLOWEROSION_VOLUME( I, EXCHANGE1%NUMBER_BANKCYCLE(I) )/DX1(1,1)/BANK_H(I)
            !END IF
            
            !IF ( JUDGE_BANK(I)==1 ) THEN !初次数据自己写
            !
            !    EXCHANGE1%EXCHANGE_BANK_SIZE(I,:,:)=0.D0  
            !    
            !    DO J = 1 , EXCHANGE1%ELEMENT
            !        X1=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,1),1) 
            !        X2=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,2),1) 
            !        X3=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,3),1) 
            !        Y1=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,1),2) 
            !        Y2=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,2),2) 
            !        Y3=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(EXCHANGE1%MESH(J,3),2)
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,10)=0.5D0*ABS(Y3*(X2-X1)+Y1*(X3-X2)+Y2*(X1-X3))
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,1)=EXCHANGE1%MESH(J,1) 
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,2)=EXCHANGE1%MESH(J,2) 
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,3)=EXCHANGE1%MESH(J,3)
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,4)=X1 
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,5)=Y1
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,6)=X2
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,7)=Y2
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,8)=X3 
            !        EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,9)=Y3
            !    END DO            
            !END IF
        
            A_X       =0.D0
            A_Y       =0.D0
            LINE      =0
            EXCHANGE1%DIG_NUMBER(I,:)=0
            NUMBER1=0 !记录侵蚀线节点个数
            NUMBER2=0 !记录开挖节点个数
LOOP4:      DO J = 1 , NUMBER_OB  !侵蚀点编号从底往上依次增加
                A_X(J) = POSITION(I,J)
                A_Y(J) = (J-1)*EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID_BANK_H/(NUMBER_OB-1)  !岸壁网格为(10m*2m)*ZOOM(I)
                NUMBER1=NUMBER1+1
                IF ( (EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID_BANK_W-A_X(J))<0.01 ) THEN
                    A_X(J) = EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID_BANK_W
                    EXIT LOOP4
                END IF
            END DO LOOP4

            IF (NUMBER1==1) GOTO 100
        
            DO J = 1 , NUMBER1-1
                LINE(J,1)=J
                LINE(J,2)=J+1
            END DO

            Y1=A_Y(NUMBER1)
            X2=A_X(1)
            DO J = 1 , EXCHANGE1%NODE     !此过程设计点与多边形位置的判断，详见附图1
                IF (EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,1)>=X2 .AND. EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,2)<=Y1) THEN
                    DO K = 1 , NUMBER1-1
                        IF ( EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,2)>=A_Y(LINE(K,1)) .AND. EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,2)<=A_Y(LINE(K,2)) ) THEN !记录交点的X坐标
                            X3=A_X(LINE(K,1)) ; Y3=A_Y(LINE(K,1))
                            X4=A_X(LINE(K,2)) ; Y4=A_Y(LINE(K,2))
                            TEMP_X=(X4-X3)/(Y4-Y3)*(EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,2)-Y3)+X3        !记录交点的X坐标
                            IF ( TEMP_X<=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,1) ) THEN
                                NUMBER2=NUMBER2+1
                                A_X(NUMBER1+NUMBER2)=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,1)
                                A_Y(NUMBER1+NUMBER2)=EXCHANGE1%ZOOMIN(I)*EXCHANGE1%GRID(J,2)
                                EXCHANGE1%DIG_NUMBER(I,J)       =1
                            END IF
                        END IF
                    END DO
                END IF
            END DO

100         CONTINUE

            NUMBER3=0
            !找出开挖单元
LOOP1:      DO J = 1 , EXCHANGE1%ELEMENT
                IF ( EXCHANGE1%DIG_NUMBER(I, INT(EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,1)) ) ==1 .AND. &
                &    EXCHANGE1%DIG_NUMBER(I, INT(EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,2)) ) ==1 .AND. &
                &    EXCHANGE1%DIG_NUMBER(I, INT(EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,3)) ) ==1       ) THEN
                    NUMBER3=NUMBER3+1
                    IF ( EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,11)==0.D0 ) EXCHANGE1%EXCHANGE_BANK_SIZE(I,J,11)=1.D0  !后续需把1==》2
                END IF
            END DO LOOP1
        END DO LOOP0

        DEALLOCATE( A_X  )
        DEALLOCATE( A_Y  )
        DEALLOCATE( LINE )  

        RETURN
    END SUBROUTINE MESH_DIVISION  
    
#endif

