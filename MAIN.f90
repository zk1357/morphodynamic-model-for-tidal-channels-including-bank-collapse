#include 'define.h'
    
PROGRAM MAIN !�ļ��űൽ900
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

    INTEGER:: I
    
    !�������
    CALL ALLOCATION
    
    !����ʼ����
    CALL INITIALIZATION
        
#ifdef PRE_BANK_COLLAPSE
    
    !DO I = 59, 59, 10
    DO I = 1, NUMBER_BV
        ! SIMULATE INITIAL BANK STRESS FIELD
        CALL INITIAL_BANK_STRESS(I)    
        
        ! OUTPUT ZOOMIN
        CALL OUTPUT_ZOOMIN(I)
    END DO
    
    !DO I = 1, 20
    !    GLOBTIME%GLOBSTEP%DAYSTEP=GLOBTIME%GLOBSTEP%DAYSTEP+1
    !    IF (I>1) THEN
    !        POSITION(1,:)=POSITION(1,1:10)-0.02
    !    END IF
    !    
    !    !̮������
    !    CALL BANK_COLLAPSE_SIMULATION
    !    
    !END DO
#endif    
    
    ! HYDRODYNAMIC MODEL
    CALL HYDRODYNAMIC
    
END PROGRAM MAIN