#include 'define.h'
  
#ifdef PRE_BANK_EROSION
SUBROUTINE FLOW_MESH_UPDATING
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

#ifdef PRE_BANK_COLLAPSE
    USE LINEAR_ELASTIC
#endif

    IMPLICIT NONE
    INTEGER   , PARAMETER                :: FMU = SELECTED_REAL_KIND(8)        
    REAL(FMU)                            :: TEMP1 , TEMP2
    INTEGER                              :: I , J , K , PP , II
    INTEGER                              :: TEMP_N

!------------------------------------------更新地貌,主要为岸壁位置

    IF (GLOBTIME%GLOBHOUR==9) THEN !低水位时更换网格
    PP=3
    LOOP0:  DO I = 1 , NUMBER_BV/2
                TEMP1=0.D0
                TEMP2=0.D0
                TEMP_N=0
                TEMP_N=INT((0.5*DY1(1,1)+LENGTH_FLOW( I , NUMBER_BANKCYCLE(I) ))/DY1(1,1))
                II=I+NUMBER_BV/2
#ifdef PRE_BANK_COLLAPSE                        
                TEMP_N=INT( ( 0.5*DY1(1,1)+RETREATS_TOTAL(I) ) /DY1(1,1) )
#endif
                !IF (GLOBTIME%GLOBSEC ==0  .AND. &
                !&    GLOBTIME%GLOBMIN ==0 .AND. &
                !&    GLOBTIME%GLOBHOUR==9 .AND. &
                !&    GLOBTIME%GLOBDAY ==2 .AND. &
                !&    I>=11 .AND. I<=20)  TEMP_N=6  
                
                IF ( TEMP_N > ABS (BVP_IN(I,1)-BVP(I,1)) ) THEN  
                
                    DO J = 1 , 5
                        !更新地形
                        DEP( BVP(I,1)-1 , BVP(I,2)-PP+J ) = DEP( BVP(I,1) , BVP(I,2)-PP+J ) 
                        DEP( BVP(II,1)+1 , BVP(II,2)-PP+J ) = DEP( BVP(II,1) , BVP(II,2)-PP+J )
                        !更新水深
                        H1 ( BVP(I,1)-1 , BVP(I,2)-PP+J ) = H1 ( BVP(I,1)   , BVP(I,2)-PP+J ) - DEP( BVP(I,1) , BVP(I,2)-PP+J ) + DEP( BVP(I,1)-1 , BVP(I,2)-PP+J )
                        H1 ( BVP(II,1)+1 , BVP(II,2)-PP+J ) = H1 ( BVP(II,1)   , BVP(II,2)-PP+J ) - DEP( BVP(II,1) , BVP(II,2)-PP+J ) + DEP( BVP(II,1)+1 , BVP(II,2)-PP+J )
                        H2 ( BVP(I,1)-1 , BVP(I,2)-PP+J ) = H2 ( BVP(I,1)   , BVP(I,2)-PP+J ) - DEP( BVP(I,1) , BVP(I,2)-PP+J ) + DEP( BVP(I,1)-1 , BVP(I,2)-PP+J )
                        H2 ( BVP(II,1)+1 , BVP(II,2)-PP+J ) = H2 ( BVP(II,1)   , BVP(II,2)-PP+J ) - DEP( BVP(II,1) , BVP(II,2)-PP+J ) + DEP( BVP(II,1)+1 , BVP(II,2)-PP+J )                    
                    
                        !!更新地形
                        !DEP( BVP(I,1)-1 , BVP(I,2)-PP+J ) = DEP( BVP(I,1)-1 , BVP(I,2)-PP+J ) + BANK_H(I) 
                        !DEP( BVP(II,1)+1 , BVP(II,2)-PP+J ) = DEP( BVP(II,1)+1 , BVP(II,2)-PP+J ) + BANK_H(II)
                        !!更新水深
                        !H1 ( BVP(I,1)-1 , BVP(I,2)-PP+J ) = H1 ( BVP(I,1)   , BVP(I,2)-PP+J ) - DEP( BVP(I,1) , BVP(I,2)-PP+J ) + DEP( BVP(I,1)-1 , BVP(I,2)-PP+J )
                        !H1 ( BVP(II,1)+1 , BVP(II,2)-PP+J ) = H1 ( BVP(II,1)   , BVP(II,2)-PP+J ) - DEP( BVP(II,1) , BVP(II,2)-PP+J ) + DEP( BVP(II,1)+1 , BVP(II,2)-PP+J )
                    
                        !坍塌土块重分布
                        !DO K = 1 , 3
                        !    !1/4
                        !    DEP( BVP(I,1)-2+K , BVP(I,2)-PP+J ) = DEP( BVP(I,1)-1+K , BVP(I,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/4.D0
                        !    DEP( BVP(II,1)+2-K , BVP(II,2)-PP+J ) = DEP( BVP(II,1)+1-K , BVP(II,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/4.D0
                        !END DO
                        !    !1/6
                        !DEP( BVP(I+1,1) , BVP(I+1,2)-PP+J ) = DEP( BVP(I+1,1) , BVP(I+1,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/6.D0
                        !DEP( BVP(I-1,1) , BVP(I-1,2)-PP+J ) = DEP( BVP(I-1,1) , BVP(I-1,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/6.D0
                        !DEP( BVP(II+1,1) , BVP(II+1,2)-PP+J ) = DEP( BVP(II+1,1) , BVP(II+1,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/6.D0
                        !DEP( BVP(II-1,1) , BVP(II-1,2)-PP+J ) = DEP( BVP(II-1,1) , BVP(II-1,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/6.D0                    
                        !    !1/12
                        !DEP( BVP(I+2,1) , BVP(I+2,2)-PP+J ) = DEP( BVP(I+2,1) , BVP(I+2,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/12.D0
                        !DEP( BVP(I-2,1) , BVP(I-2,2)-PP+J ) = DEP( BVP(I-2,1) , BVP(I-2,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/12.D0
                        !DEP( BVP(II+2,1) , BVP(II+2,2)-PP+J ) = DEP( BVP(II+2,1) , BVP(II+2,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/12.D0
                        !DEP( BVP(II-2,1) , BVP(II-2,2)-PP+J ) = DEP( BVP(II-2,1) , BVP(II-2,2)-PP+J ) - UPDATE_C(I)*BANK_H(I)/5.D0/12.D0
                    
                    END DO
                    
                    !含沙量给为0
                    !C1(:,:) = 0.D0
                    !C2(:,:) = 0.D0
                    
                    !平滑
                    !CALL SMOOTH
                
                    BVP(I,1)=BVP(I,1)-1
                    BVP(II,1)=BVP(II,1)+1
                
                    !WRITE(*,*) UPDATE_C(I)
                    !PAUSE
                
                    UPDATE_C(I)=0.D0
            
                BANK_H(I)=MAX( ABS( DEP(BVP(I,1)-1 , BVP(I,2)) - DEP(BVP(I,1) , BVP(I,2)) ) , ABS( DEP(BVP(I,1)+1 , BVP(I,2)) - DEP(BVP(I,1) , BVP(I,2)) ) )           
                BANK_H(II)=MAX( ABS( DEP(BVP(II,1)-1 , BVP(II,2)) - DEP(BVP(II,1) , BVP(II,2)) ) , ABS( DEP(BVP(II,1)+1 , BVP(II,2)) - DEP(BVP(II,1) , BVP(II,2)) ) )
            
                END IF
            
            END DO LOOP0
        END IF
END SUBROUTINE
#endif
    
#ifdef PRE_BANK_COLLAPSE    
SUBROUTINE BANK_UPDATING(STEP,DSTEP)
    
    USE BANK_EROSION
    USE EXCHANGE
    USE HYDRO
    
    IMPLICIT NONE
    INTEGER   , PARAMETER                :: BU = SELECTED_REAL_KIND(8)
    !REAL(BU) , DIMENSION(:,:)            :: TIME_DIF(NUMBER_BV/2,300) !记录每次破坏在涨潮还是落潮
    INTEGER                              :: I , J , K , STEP , DSTEP
    REAL(BU)                             :: TEMP1      
    CHARACTER*40                         :: FILENAME 
    
    DO I = 1 , NUMBER_BV/2
        IF ( JUDGE_BANK(I)==1 ) THEN !根据应力状态判断是否更换为初始岸壁网格
            POSITION(I,:)=10.D0
        END IF
    END DO

    IF ( DSTEP/10==DSTEP/10. ) THEN
        FILENAME=''
        WRITE(FILENAME(1:4) , '(I4.4)') DSTEP
        FILENAME(5:15)='FAILURE.TXT'
        
        OPEN(UNIT=600 , FILE=FILENAME)
LOOP1:  DO I = 1 , NUMBER_BV/2
            WRITE(600,*)  BF(I,1) , BF(I,2) , BF(I,3)
        END DO LOOP1
        CLOSE(600)
    END IF

END SUBROUTINE
#endif