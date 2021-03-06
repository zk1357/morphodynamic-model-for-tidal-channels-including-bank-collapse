#include 'define.h'
#ifdef PRE_BANK_COLLAPSE    
SUBROUTINE FLOWEROSION_FAILURE(NUMBER_I)
    USE TIME_DIN
#ifdef PRE_BANK_EROSION
    USE BANK_EROSION
    USE DEF_EXCHANGE
#endif

#ifdef PRE_BANK_COLLAPSE
    USE BANK_COLLAPSE
#endif
    
    IMPLICIT NONE
    INTEGER  :: NUMBER_I
    
    EXCHANGE1%COLLAPSE_TIME(NUMBER_I, EXCHANGE1%NUMBER_BANKCYCLE(NUMBER_I))= GLOBTIME%GLOBSEC 
        
    EXCHANGE1%COLLAPSE_DEF(NUMBER_I, EXCHANGE1%NUMBER_BANKCYCLE(NUMBER_I)) = 2
        
    EXCHANGE1%COLLAPSE_VOLUME(NUMBER_I, EXCHANGE1%NUMBER_BANKCYCLE(NUMBER_I)) = 0
        
    EXCHANGE1%FLOW_EXCESS_JUDGE(NUMBER_I) = 0.D0
        
    EXCHANGE1%JUDGE_BANK(NUMBER_I) = 1

    END SUBROUTINE
#endif    