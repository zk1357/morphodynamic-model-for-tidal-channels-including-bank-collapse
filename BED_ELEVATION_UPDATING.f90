#include 'define.h'
#ifdef PRE_MOR

SUBROUTINE BED_ELEVATION_UPDATING
    
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
    REAL*8                              :: DTT , TEMP , QD , TB , TB_BANK , QE , QE_BANK , &
                                        &  A , B , CX , CY , U , V  , UV , &
                                        &  DXX , DYY , HL , HR, CHEZY
    INTEGER                             :: I , J , K , P , JC
    
    ! TRANSFER TO SSC坡脚掩护泥沙体积需计算 
        ! BED LEVEL CHANGES INDUCED BY BED EROSION
LOOP0:  DO I = M_START+1 , M_END-1
            DO J = N_START+1 , N_END-1
                !IF ( I>TCI_START .AND. I<TCI_END ) CYCLE
                IF ( NDD(I,J)==-1 ) CYCLE
                
                !水深小于临界值，不算
                IF (H1(I,J)<= FLOODINGDEPTHX) CYCLE
                    
                CALL UV_CACULATION(I,J,U,V,UV)
                
                TEMP = 0.D0

#ifdef PRE_CHEZY 
                !CHEZY公式
                CHEZY=CHE
#endif      

#ifdef PRE_MANNING 
                !MANNING公式
                CHEZY=MANNING*H1(I,J)**(1.D0/6.D0)
#endif     

#ifdef PRE_AD
                TB= FLOWDENSITY*9.8D0/CHEZY**2*U**2                
                QE= MAX( 0.D0 , ME*(TB/TE-1.D0) )
                QD= C1(I,J)*WS*MAX(0.D0,(1.D0-TB/TD))
                !右边除以(1-p),p为孔隙率，取值0.4
                TEMP= MOR_FACTOR*(QE-QD)*DT/SANDDENSITY/(1-PORO)
#endif

#ifdef PRE_BANK_COLLAPSE
                ! EFFECT OF TOE SHIELD
                CALL COLLAPSE_TOE_SHIELD(I, J, TEMP)
#endif
                !!!!!!!!!!!!!!!!!!!!!!!!
                !TEMP=0.001819/20.D0*((2*P1(I,J)/(H1(I,J)+H1(I,J+1)))**5-(2*P1(I,J-1)/(H1(I,J)+H1(I,J-1)))**5)


                !!!!!!!!!!!!!!!!!!!!!!!
                
                DEP(I,J)=DEP(I,J)+TEMP
                
            END DO
        END DO LOOP0      
        
        
#ifdef PRE_BANK_EROSION
        !根据DRY CELL EROSION METHOD,更新岸壁DEPTH
        DO I = 1, NUMBER_BV
            ! WHETHER BANK EROSION HAPPENS
            IF ( EXCHANGE1%FLOWEROSION_CHECK(I)==1 ) CYCLE
            
            TEMP=MOR_FACTOR_BANK*S_BANK_EROSION(I)*H1( BVP(I,1), BVP(I,2) )*DX1(1,1)*DY1(1,1)/SANDDENSITY_BANK
            
#ifdef PRE_BANK_COLLAPSE
            ! EFFECT OF TOE SHIELD
            IF ( S_BANK_COLLAPSE(I) > 0 ) THEN
                IF ( S_BANK_COLLAPSE(I) >= TEMP ) THEN
                    S_BANK_COLLAPSE(I)= S_BANK_COLLAPSE(I)-TEMP
                    TEMP= 0.D0
                ELSE
                    TEMP= TEMP-S_BANK_COLLAPSE(I)
                    S_BANK_COLLAPSE(I)=0.D0
                END IF
            END IF
#endif      

            ! BED LEVEL CHANGES INDUCED BY LATERAL FLOW EROSION
            DO K = 1 , 5
                IF ( I<=NUMBER_BV/2 ) THEN
                    DEP( BVP(I,1)-1, BVP(I,2)-3+K )=DEP( BVP(I,1)-1, BVP(I,2)-3+K )+TEMP/DX1(1,1)/DY1(1,1)     
                ELSE
                    DEP( BVP(I,1)+1, BVP(I,2)-3+K )=DEP( BVP(I,1)+1, BVP(I,2)-3+K )+TEMP/DX1(1,1)/DY1(1,1)
                END IF
            END DO
            
            ! CHANGE BANK HEIGHT
            IF ( I<=NUMBER_BV/2 ) THEN
                BANK_H(I)=-DEP( BVP(I,1)-1, BVP(I,2) )+DEP( BVP(I,1), BVP(I,2) )
            ELSE
                BANK_H(I)=-DEP( BVP(I,1)+1, BVP(I,2) )+DEP( BVP(I,1), BVP(I,2) )    
            END IF
            
            ! CHANGE POSITION OF BVP
            IF ( BANK_H(I)<=BANK_EROSION_LIMIT ) THEN 
                IF ( I<=NUMBER_BV/2 ) THEN
                    BVP(I,1)=BVP(I,1)-1
                    BANK_H(I)=-DEP( BVP(I,1)-1, BVP(I,2) )+DEP( BVP(I,1), BVP(I,2) )
                ELSE
                    BVP(I,1)=BVP(I,1)+1
                    BANK_H(I)=-DEP( BVP(I,1)+1, BVP(I,2) )+DEP( BVP(I,1), BVP(I,2) )
                END IF
                
                IF ( BANK_H(I) <= BANK_EROSION_LIMIT ) EXCHANGE1%FLOWEROSION_CHECK(I)=1
                
#ifdef PRE_BANK_COLLAPSE                
                ! BANK UPDATA
                CALL BANK_UPDATA(I)
#endif    

            END IF         
            
            !边壁侵蚀含沙量归零    
            S_BANK_EROSION(I)=0.d0
            
        END DO     
#endif
    
        RETURN
END SUBROUTINE BED_ELEVATION_UPDATING

#endif