#include 'define.h' 
#ifdef PRE_BANK_COLLAPSE    
! 计算坡脚掩护作用    
SUBROUTINE COLLAPSE_TOE_SHIELD(NUMBER_I, NUMBER_J, SHIELD)
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

#ifdef PRE_BANK_COLLAPSE
        USE BANK_COLLAPSE
#endif
    
        IMPLICIT NONE
        REAL*8               :: TEMP, SHIELD, X
        INTEGER              :: I, J, K, NUMBER_I, NUMBER_J
        CHARACTER*40         :: FILENAME, FILENAME1
        
        DO I = 1, NUMBER_BV
            
            IF ( S_BANK_COLLAPSE(I) == 0 ) CYCLE
            
            IF ( BVP(I, 1) .NE. NUMBER_I ) CYCLE
            
            IF ( BVP(I, 2)-2 <= NUMBER_J .AND. BVP(I, 2)+2 >= NUMBER_J  ) THEN
                
                TEMP= SHIELD*DX1(1,1)*DY1(1,1)
                
                IF ( S_BANK_COLLAPSE(I) > TEMP ) THEN
                    S_BANK_COLLAPSE(I)= S_BANK_COLLAPSE(I)-TEMP
                    TEMP = 0.D0
                ELSE
                    TEMP= TEMP-S_BANK_COLLAPSE(I)
                    S_BANK_COLLAPSE(I) = 0.D0
                END IF
                
            END IF
            
        END DO
        

        RETURN
    END SUBROUTINE COLLAPSE_TOE_SHIELD
#endif    