#include 'define.h'   
#ifdef PRE_BANK_COLLAPSE
SUBROUTINE BANK_COLLAPSE_SIMULATION
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
#endif

#ifdef PRE_BANK_COLLAPSE
    USE BANK_COLLAPSE
#endif 

    IMPLICIT NONE
    
    INTEGER :: I, J, NUM_J
    
    
    ! SIMULATE BANK PROFILE
    CALL MESH_DIVISION
    
    DO I = 1, NUMBER_BV
    !DO I = 40, 100, 50
        !WRITE(*,*) I
        ! WHETHER BANK EROSION HAPPENS
        IF ( EXCHANGE1%FLOWEROSION_CHECK(I)==1 ) CYCLE
        
        ! WHETHER EXCESS FLOW EROSION HAPPENS
        IF ( EXCHANGE1%FLOW_EXCESS_JUDGE(I) .EQ. 1 ) THEN
            CALL FLOWEROSION_FAILURE(I)
            !WRITE(*,*) 1
            CYCLE
        END IF
                
        !WRITE(*,*) BANK_H(I)
        !CALL OUTPUT_POSITION
        
        ! WHETHER BANK EROSION HAPPENS
        IF (BANK_H(I)<=BANK_COLLAPSE_LIMIT) CYCLE
        
        !WRITE(*,*) 2
        
        ! SIMULATE EXTERNAL FORCE
        CALL BANK_EXTERNAL_FORCE(I, NUM_J)

        !WRITE(*,*) 3
        
        ! NO MESH IS REMOVED
        IF (NUM_J <= 0 ) CYCLE
        
LOOP2:  DO J = 1, NUM_J
            
            ! CALCULATE STRESS CHANGES
            CALL LINEAR_ELASTIC_STIFFNESS(I, NUM_J)
            
            ! FAILURE JUDGEMENT
            IF ( EXCHANGE1%JUDGE_BANK(I) == 1 ) EXIT LOOP2
            
        END DO LOOP2
        
    END DO
    
    DO I =1, NUMBER_BV
    !DO I = 40, 100, 50
        ! VARIATE EXCHANGE AND RESET
        CALL COLLAPSE_EXCHANGES(I)
    END DO
    
    !CALL OUTPUT_MOR_M
    
    !CALL OUTPUT_BANKEROSION_VOLUME
    
    END SUBROUTINE BANK_COLLAPSE_SIMULATION  
#endif    
    


