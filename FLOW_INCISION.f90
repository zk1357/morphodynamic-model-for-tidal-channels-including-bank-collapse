#include 'define.h'
    
#ifdef PRE_BANK_EROSION    
SUBROUTINE FLOW_INCISION !文件号从40起
 
    USE TIME_DIN
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
    
    IMPLICIT NONE
    INTEGER   , PARAMETER                    :: NFI=SELECTED_REAL_KIND(8)
    INTEGER                                  :: I , J , K , STEP , II
    INTEGER                                  :: NODE , ELEMENT
    REAL(NFI)                                :: U , V , UV , TB_BANK , TEMP , TEMP1 , TEMP_SH, QE_BANK
    REAL(NFI)                                :: Y , CHEZY, ERO_VOLUME
    CHARACTER*30                             :: FILENAME

    ! 计算侧侵蚀距离     
    DO I= 1 , NUMBER_BV
        
        !II=I+NUMBER_BV/2
        
        ! WHETHER BANK EROSION HAPPENS
        IF ( EXCHANGE1%FLOWEROSION_CHECK(I)==1 ) CYCLE        
        IF (BANK_H(I)<=BANK_EROSION_LIMIT) CYCLE
        
        CALL UV_CACULATION(BVP(I,1),BVP(I,2),U,V,UV)
#ifdef PRE_CHEZY         
        !谢才公式
        CHEZY=CHE
#endif        
#ifdef PRE_MANNING  
        !MANNING公式
        CHEZY=MANNING*H1( BVP(I,1) , BVP(I,2) )**(1.D0/6.D0)
#endif
        TB_BANK= 1000*9.8D0/CHEZY**2*UV**2
        QE_BANK= ME_BANK*MAX( (TB_BANK/TE_BANK-1.D0) , 0.D0 )
        TEMP   = MOR_FACTOR_BANK*QE_BANK*DT/SANDDENSITY_BANK
        S_BANK_EROSION(I)=0.D0
        Y    = 0.D0
        
        IF (TEMP<=0.D0) CYCLE
            
#ifdef PRE_BANK_COLLAPSE
        ! LATERAL FLOW EROSION VOLUME
        IF ( H1( BVP(I,1) , BVP(I,2) ) >= BANK_H(I) ) THEN
            ERO_VOLUME= TEMP*BANK_H(I)*DX1(1,1)
        ELSE
            ERO_VOLUME= TEMP*H1( BVP(I,1) , BVP(I,2) )*DX1(1,1)
        END IF

        ! EFFECT OF TOE SHIELD
        IF ( S_BANK_COLLAPSE(I) > 0 ) THEN
            IF ( S_BANK_COLLAPSE(I) >= ERO_VOLUME ) THEN
                S_BANK_COLLAPSE(I)= S_BANK_COLLAPSE(I)-ERO_VOLUME
                ! 记录每个岸壁BANK SHIELD
                EXCHANGE1%SHIELD_VOLUME(I, EXCHANGE1%NUMBER_BANKCYCLE(I))=EXCHANGE1%SHIELD_VOLUME(I, EXCHANGE1%NUMBER_BANKCYCLE(I))+ERO_VOLUME
                ! 记录每天BANK SHIELD
                EXCHANGE1%BANKEROSION_VOLUME(I,3)=EXCHANGE1%BANKEROSION_VOLUME(I,3)+ERO_VOLUME
                ERO_VOLUME= 0.D0
            ELSE
                ERO_VOLUME= TEMP-S_BANK_COLLAPSE(I)
                ! 记录每个岸壁BANK SHIELD
                EXCHANGE1%SHIELD_VOLUME(I, EXCHANGE1%NUMBER_BANKCYCLE(I))=EXCHANGE1%SHIELD_VOLUME(I, EXCHANGE1%NUMBER_BANKCYCLE(I))+S_BANK_COLLAPSE(I)
                ! 记录每天BANK SHIELD
                EXCHANGE1%BANKEROSION_VOLUME(I,3)=EXCHANGE1%BANKEROSION_VOLUME(I,3)+ S_BANK_COLLAPSE(I)
                S_BANK_COLLAPSE(I)=0.D0
            END IF
        END IF
        
        ! REST LATERAL FLOW EROSION VOLUME
        IF ( H1( BVP(I,1) , BVP(I,2) ) >= BANK_H(I) ) THEN
            TEMP= ERO_VOLUME/BANK_H(I)/DX1(1,1)
        ELSE
            TEMP= ERO_VOLUME/H1( BVP(I,1) , BVP(I,2) )/DX1(1,1)
        END IF
#endif   
        
        IF ( H1( BVP(I,1) , BVP(I,2) ) >= BANK_H(I) ) THEN
            Y = BANK_H(I)
            POSITION(I,:)=POSITION(I,:)-TEMP
        ELSE
            DO J= 1 , NUMBER_OB 
                IF ( H1( BVP(I,1) , BVP(I,2) ) > (J-1)*BANK_H(I)/(NUMBER_OB-1) ) THEN
                    POSITION(I,J)=POSITION(I,J)-TEMP
                    Y = (J-1)*BANK_H(I)/(NUMBER_OB-1)
                END IF
            END DO
        END IF
        
        ! 记录水流侵蚀引起的SSC增加
        S_BANK_EROSION(I) = QE_BANK*DT/DY1(1,1)*Y/H1( BVP(I,1) , BVP(I,2) ) !此公式正确，推导过几次了
        !S_BANK_EROSION(II) = S_BANK_EROSION(I)
        
        ! 记录水流侵蚀体积（每个断面的每块岸壁）
        EXCHANGE1%FLOWEROSION_VOLUME( I, EXCHANGE1%NUMBER_BANKCYCLE(I) ) = TEMP*Y*DX1(1,1) + &
     &  EXCHANGE1%FLOWEROSION_VOLUME( I, EXCHANGE1%NUMBER_BANKCYCLE(I) )      
        
     !   EXCHANGE1%FLOWEROSION_VOLUME( II, EXCHANGE1%NUMBER_BANKCYCLE(II) ) =  &
     !&  EXCHANGE1%FLOWEROSION_VOLUME( I, EXCHANGE1%NUMBER_BANKCYCLE(I) )
        
        !记录每天FLOW EROSION 
        EXCHANGE1%BANKEROSION_VOLUME(I, 1) = EXCHANGE1%BANKEROSION_VOLUME(I, 1) + TEMP*Y*DX1(1,1)
        
    END DO
    
    RETURN
    
END SUBROUTINE FLOW_INCISION
#endif