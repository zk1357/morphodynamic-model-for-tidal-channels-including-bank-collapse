#include 'define.h'

SUBROUTINE BED_SY_SMOOTH
    
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
    REAL*8                              :: TEMP
    INTEGER                             :: I, J, K, JC
    
    JC=(M_START+M_END)/2
LOOP0:  DO J = N_START+1 , N_END-1
            DO I = 1 , M_END-JC-1
                TEMP=(DEP(JC+I,J)+DEP(JC-I,J))/2.D0
                DEP(JC+I,J)=TEMP
                DEP(JC-I,J)=TEMP
            END DO
        END DO LOOP0
       


    
        RETURN
END SUBROUTINE BED_SY_SMOOTH
